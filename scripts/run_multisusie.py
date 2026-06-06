"""
Multi-ancestry fine-mapping using MultiSuSiE.

For each clumped locus (union of lead SNPs across ancestries), extracts GWAS
summary statistics and LD matrices per ancestry, then runs MultiSuSiE's
multisusie_rss to produce cross-ancestry credible sets.
"""

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

import MultiSuSiE


def available_memory_mb():
    """Slurm allocated memory (MB) when under Slurm, else total system memory."""
    slurm = os.environ.get("SLURM_MEM_PER_NODE", "").strip()
    if slurm:
        try:
            value = float(slurm)
            if value > 0:
                return value
        except ValueError:
            pass
    if sys.platform == "linux":
        try:
            with open("/proc/meminfo", encoding="utf-8") as meminfo:
                line = meminfo.readline()
            kb = int("".join(c for c in line if c.isdigit()))
            return kb / 1024
        except (OSError, ValueError):
            pass
    elif sys.platform == "darwin":
        try:
            raw = subprocess.check_output(
                ["sysctl", "-n", "hw.memsize"], stderr=subprocess.DEVNULL
            )
            return int(raw.strip()) / (1024 ** 2)
        except (subprocess.SubprocessError, ValueError):
            pass
    return None


def available_cpus():
    """Logical CPUs on this host."""
    cpus = os.cpu_count()
    return max(1, cpus if cpus is not None else 1)


def calculate_parallelism(max_workers=10, memory_per_worker_mb=8000):
    slurm = os.environ.get("SLURM_CPUS_ON_NODE", "").strip()
    if slurm:
        try:
            return max(1, int(slurm) - 1)
        except ValueError:
            pass
    cpus = available_cpus()
    mem = available_memory_mb()
    if mem is not None and mem > 0:
        by_mem = int(mem // memory_per_worker_mb)
        return max(1, min(max_workers, by_mem, cpus))
    return max(1, min(max_workers, cpus))


def parse_args():
    parser = argparse.ArgumentParser(description="Run MultiSuSiE multi-ancestry fine-mapping")
    parser.add_argument("--gwas_files", nargs="+", required=True,
                        help="Standardised GWAS files (one per ancestry)")
    parser.add_argument("--clumped_files", nargs="+", required=True,
                        help="Plink clumped files (one per ancestry)")
    parser.add_argument("--ancestries", nargs="+", required=True,
                        help="Ancestry codes matching 1000 Genomes panels (e.g. EUR EAS)")
    parser.add_argument("--sample_sizes", nargs="+", type=int, required=True,
                        help="GWAS sample sizes (one per ancestry)")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for per-locus results")
    parser.add_argument("--window_kb", type=int, default=1000,
                        help="Half-width of fine-mapping window in kb")
    parser.add_argument("--max_causal", type=int, default=10,
                        help="Maximum number of causal signals (L)")
    parser.add_argument("--coverage", type=float, default=0.95,
                        help="Credible set coverage")
    parser.add_argument("--min_abs_corr", type=float, default=0.5,
                        help="Minimum absolute correlation for credible sets")
    parser.add_argument("--thousand_genomes_dir", required=True,
                        help="Path to 1000 Genomes reference panel directory")
    parser.add_argument("--completion_file", required=True,
                        help="Sentinel file written on success")
    return parser.parse_args()


def read_gwas(gwas_file):
    """Read a standardised GWAS TSV (possibly gzipped)."""
    return pd.read_csv(gwas_file, sep="\t", dtype={"CHR": str})


def read_clumped(clumped_file):
    """Read plink --clump output (whitespace-delimited)."""
    try:
        df = pd.read_csv(clumped_file, sep=r"\s+", engine="python")
        if "SNP" not in df.columns or len(df) == 0:
            return pd.DataFrame(columns=["SNP", "CHR", "BP"])
        return df[["SNP", "CHR", "BP"]]
    except Exception:
        return pd.DataFrame(columns=["SNP", "CHR", "BP"])


def get_union_loci(clumped_dfs, window_bp):
    """Merge lead SNPs across ancestries into non-overlapping loci."""
    all_leads = pd.concat(clumped_dfs, ignore_index=True)
    if len(all_leads) == 0:
        return []

    all_leads["CHR"] = all_leads["CHR"].astype(str)
    all_leads["BP"] = all_leads["BP"].astype(int)
    all_leads = all_leads.sort_values(["CHR", "BP"]).reset_index(drop=True)

    loci = []
    current_chr = None
    current_start = None
    current_end = None
    current_leads = []

    for _, row in all_leads.iterrows():
        bp = int(row["BP"])
        chrom = str(row["CHR"])
        if current_chr is None or chrom != current_chr or bp - current_end > window_bp:
            if current_chr is not None:
                loci.append({
                    "chr": current_chr,
                    "start": current_start,
                    "end": current_end,
                    "center": int((current_start + current_end) / 2),
                    "leads": current_leads
                })
            current_chr = chrom
            current_start = bp - window_bp
            current_end = bp + window_bp
            current_leads = [row["SNP"]]
        else:
            current_end = max(current_end, bp + window_bp)
            current_leads.append(row["SNP"])

    if current_chr is not None:
        loci.append({
            "chr": current_chr,
            "start": current_start,
            "end": current_end,
            "center": int((current_start + current_end) / 2),
            "leads": current_leads
        })

    return loci


def _read_plink_ld_matrix(ld_file, n_snps):
    """Parse plink --r square output into an n_snps x n_snps matrix."""
    try:
        R = np.loadtxt(ld_file)
    except (OSError, ValueError):
        return None

    if n_snps < 2:
        return None
    if R.ndim == 0:
        return None
    if R.ndim == 1:
        if len(R) == n_snps * n_snps:
            R = R.reshape(n_snps, n_snps)
        else:
            return None
    if R.shape[0] != n_snps or R.shape[1] != n_snps:
        return None
    return R


def compute_ld_matrix(rsids, chrom, ancestry, thousand_genomes_dir):
    """Compute LD correlation matrix from 1000 Genomes via plink."""
    if len(rsids) < 2:
        return None, None

    bfile = os.path.join(thousand_genomes_dir, ancestry)
    tmpdir = tempfile.mkdtemp(prefix="ld_")
    snp_file = os.path.join(tmpdir, "snps.txt")
    out_prefix = os.path.join(tmpdir, "ld")

    try:
        with open(snp_file, "w", encoding="utf-8") as f:
            f.write("\n".join(rsids) + "\n")

        cmd = [
            "plink1.9",
            "--bfile", bfile,
            "--chr", str(chrom),
            "--extract", snp_file,
            "--r", "square",
            "--out", out_prefix,
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        ld_file = out_prefix + ".ld"

        if result.returncode != 0 or not os.path.exists(ld_file):
            return None, None

        bim_file = out_prefix + ".bim"
        if os.path.exists(bim_file):
            bim = pd.read_csv(bim_file, sep="\t", header=None, usecols=[1], names=["SNP"])
            snp_order = bim["SNP"].tolist()
        else:
            bim_ref = bfile + ".bim"
            bim = pd.read_csv(
                bim_ref, sep="\t", header=None, usecols=[0, 1, 3],
                names=["CHR", "SNP", "BP"],
            )
            bim = bim[(bim["CHR"].astype(str) == str(chrom)) & (bim["SNP"].isin(rsids))]
            bim = bim.sort_values("BP")
            snp_order = bim["SNP"].tolist()

        R = _read_plink_ld_matrix(ld_file, len(snp_order))
        if R is None:
            return None, None

        return R, snp_order
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _validate_multisusie_inputs(b_list, s_list, R_list, ancestries, locus_label):
    """Ensure every population has the same variant count before MultiSuSiE."""
    p = int(R_list[0].shape[0])
    for i, ancestry in enumerate(ancestries):
        b_arr = np.asarray(b_list[i])
        s_arr = np.asarray(s_list[i])
        r_arr = np.asarray(R_list[i])
        if (
            b_arr.ndim != 1
            or s_arr.ndim != 1
            or r_arr.ndim != 2
            or b_arr.shape[0] != p
            or s_arr.shape[0] != p
            or r_arr.shape != (p, p)
        ):
            print(
                f"  MultiSuSiE input mismatch for {locus_label} "
                f"({ancestry}): b={b_arr.shape}, s={s_arr.shape}, "
                f"R={r_arr.shape}, expected P={p}",
                file=sys.stderr,
            )
            return False
    return True


def estimate_var_y(se, eaf, n):
    """Estimate outcome variance from GWAS SE, EAF, and sample size."""
    se = np.asarray(se, dtype=float)
    if eaf is None:
        return 1.0
    eaf = np.asarray(eaf, dtype=float)
    valid = np.isfinite(se) & (se > 0) & np.isfinite(eaf) & (eaf > 0) & (eaf < 1)
    if valid.sum() < 2:
        return 1.0
    oneover = 1.0 / (se[valid] ** 2)
    nvx = 2 * n * eaf[valid] * (1 - eaf[valid])
    denom = np.dot(oneover, oneover)
    if denom <= 0:
        return 1.0
    var_y = np.dot(nvx, oneover) / denom
    return float(var_y) if var_y > 0 else 1.0


def assign_credible_set_columns(sets, n_variants):
    """Map MultiSuSiE credible sets to IN_CS and CS_ID per variant."""
    in_cs = np.zeros(n_variants, dtype=int)
    cs_id = np.full(n_variants, np.nan)
    if sets is None or not isinstance(sets, (list, tuple)) or len(sets) == 0:
        return in_cs, cs_id

    cs_list = sets[0]
    include_mask = sets[3] if len(sets) > 3 else [True] * len(cs_list)
    cs_counter = 0
    for l_idx, indices in enumerate(cs_list):
        if l_idx >= len(include_mask) or not include_mask[l_idx] or len(indices) == 0:
            continue
        cs_counter += 1
        for idx in np.asarray(indices, dtype=int):
            if 0 <= idx < n_variants:
                in_cs[idx] = 1
                if np.isnan(cs_id[idx]):
                    cs_id[idx] = cs_counter
    return in_cs, cs_id


def build_multisusie_result_df(fit, final_rsids, ref_gwas, ancestries):
    """Build per-variant output table from a MultiSuSiE fit."""
    n = len(final_rsids)
    result_df = pd.DataFrame({
        "SNP": ref_gwas["SNP"].values if "SNP" in ref_gwas.columns else final_rsids,
        "CHR": ref_gwas["CHR"].values,
        "BP": ref_gwas["BP"].values,
        "RSID": final_rsids,
    })

    if hasattr(fit, "pip") and fit.pip is not None and len(fit.pip) == n:
        result_df["JOINT_PIP"] = fit.pip
    else:
        result_df["JOINT_PIP"] = np.nan

    in_cs, cs_id = assign_credible_set_columns(
        getattr(fit, "sets", None), n
    )
    result_df["IN_CS"] = in_cs
    result_df["CS_ID"] = pd.array(cs_id, dtype="Int64")

    if hasattr(fit, "coef") and fit.coef is not None:
        for k, ancestry in enumerate(ancestries):
            if fit.coef.shape[0] > k and len(fit.coef[k]) == n:
                result_df[f"BETA_{ancestry}"] = fit.coef[k]

    return result_df


def build_credible_sets_df(fit, final_rsids):
    """Summarise each MultiSuSiE credible set for locus_credible_sets.tsv."""
    sets = getattr(fit, "sets", None)
    columns = ["CS_ID", "CS_SIZE", "PURITY", "MAX_PIP_RSID", "MAX_PIP", "LBF"]
    if sets is None or not isinstance(sets, (list, tuple)) or len(sets) == 0:
        return pd.DataFrame(columns=columns)

    cs_list = sets[0]
    purity_list = sets[1] if len(sets) > 1 else [np.nan] * len(cs_list)
    include_mask = sets[3] if len(sets) > 3 else [True] * len(cs_list)
    pip = (
        np.asarray(fit.pip, dtype=float)
        if hasattr(fit, "pip") and fit.pip is not None
        else np.full(len(final_rsids), np.nan)
    )
    lbf = (
        np.asarray(fit.lbf, dtype=float)
        if hasattr(fit, "lbf") and fit.lbf is not None
        else np.full(len(cs_list), np.nan)
    )

    rows = []
    cs_counter = 0
    for l_idx, indices in enumerate(cs_list):
        if l_idx >= len(include_mask) or not include_mask[l_idx] or len(indices) == 0:
            continue
        cs_counter += 1
        indices = np.asarray(indices, dtype=int)
        cs_pip = pip[indices]
        best_local = int(np.nanargmax(cs_pip))
        best_idx = int(indices[best_local])
        rows.append({
            "CS_ID": cs_counter,
            "CS_SIZE": len(indices),
            "PURITY": purity_list[l_idx] if l_idx < len(purity_list) else np.nan,
            "MAX_PIP_RSID": final_rsids[best_idx],
            "MAX_PIP": float(cs_pip[best_local]),
            "LBF": float(lbf[l_idx]) if l_idx < len(lbf) else np.nan,
        })

    return pd.DataFrame(rows, columns=columns)


def extract_locus_data(gwas_df, chrom, start, end):
    """Extract GWAS data for a genomic window."""
    mask = (
        (gwas_df["CHR"].astype(str) == str(chrom)) &
        (gwas_df["BP"].astype(float) >= start) &
        (gwas_df["BP"].astype(float) <= end) &
        gwas_df["BETA"].notna() &
        gwas_df["SE"].notna() &
        (gwas_df["SE"].astype(float) > 0)
    )
    return gwas_df[mask].copy()


def run_multisusie_for_locus(locus, gwas_dfs, ancestries, sample_sizes,
                             thousand_genomes_dir, max_causal, coverage,
                             min_abs_corr):
    """Run MultiSuSiE on a single locus across ancestries."""
    chrom = locus["chr"]
    start = locus["start"]
    end = locus["end"]

    ancestry_data = []
    for i, (gwas_df, ancestry, n) in enumerate(zip(gwas_dfs, ancestries, sample_sizes)):
        locus_gwas = extract_locus_data(gwas_df, chrom, start, end)
        if len(locus_gwas) < 2 or "RSID" not in locus_gwas.columns:
            return None
        locus_gwas = locus_gwas[locus_gwas["RSID"].notna() & (locus_gwas["RSID"].str.len() > 0)]
        if len(locus_gwas) < 2:
            return None
        ancestry_data.append({
            "gwas": locus_gwas,
            "ancestry": ancestry,
            "n": n
        })

    shared_rsids = set(ancestry_data[0]["gwas"]["RSID"])
    for ad in ancestry_data[1:]:
        shared_rsids &= set(ad["gwas"]["RSID"])

    if len(shared_rsids) < 2:
        return None

    shared_rsids_list = sorted(shared_rsids)

    ancestry_ld = []
    for ad in ancestry_data:
        gwas_sub = ad["gwas"][ad["gwas"]["RSID"].isin(shared_rsids_list)].copy()
        gwas_sub = gwas_sub.drop_duplicates(subset="RSID")
        gwas_sub = gwas_sub.set_index("RSID").loc[shared_rsids_list].reset_index()

        R, snp_order = compute_ld_matrix(
            shared_rsids_list, chrom, ad["ancestry"], thousand_genomes_dir
        )
        if R is None:
            return None

        ancestry_ld.append({
            "gwas_sub": gwas_sub,
            "R": R,
            "snp_order": snp_order,
            "ancestry": ad["ancestry"],
            "n": ad["n"],
        })

    final_rsids = set(shared_rsids_list)
    for entry in ancestry_ld:
        final_rsids &= set(entry["snp_order"])
    final_rsids = sorted(final_rsids)
    if len(final_rsids) < 2:
        return None

    b_list = []
    s_list = []
    R_list = []
    n_list = []
    varY_list = []

    locus_label = f"chr{chrom}:{start}-{end}"
    for entry in ancestry_ld:
        snp_order = entry["snp_order"]
        snp_to_idx = {rsid: idx for idx, rsid in enumerate(snp_order)}
        ld_idx = [snp_to_idx[r] for r in final_rsids]
        R_sub = np.ascontiguousarray(entry["R"][np.ix_(ld_idx, ld_idx)], dtype=np.float64)

        gwas_sub = entry["gwas_sub"].set_index("RSID").reindex(final_rsids)
        if gwas_sub["BETA"].isna().any() or gwas_sub["SE"].isna().any():
            return None
        beta = np.ascontiguousarray(
            gwas_sub["BETA"].astype(float).to_numpy().ravel(), dtype=np.float64
        )
        se = np.ascontiguousarray(
            gwas_sub["SE"].astype(float).to_numpy().ravel(), dtype=np.float64
        )
        eaf = gwas_sub["EAF"].astype(float).values if "EAF" in gwas_sub.columns else None

        b_list.append(beta)
        s_list.append(se)
        R_list.append(R_sub)
        n_list.append(entry["n"])
        varY_list.append(estimate_var_y(se, eaf, entry["n"]))

    if not _validate_multisusie_inputs(b_list, s_list, R_list, ancestries, locus_label):
        return None

    try:
        fit = MultiSuSiE.multisusie_rss(
            b_list=b_list,
            s_list=s_list,
            R_list=R_list,
            varY_list=varY_list,
            population_sizes=n_list,
            L=max_causal,
            coverage=coverage,
            min_abs_corr=min_abs_corr,
            variant_ids=final_rsids,
            low_memory_mode=False,
            float_type=np.float64,
        )
    except Exception as e:
        print(f"  MultiSuSiE failed for locus chr{chrom}:{start}-{end}: {e}", file=sys.stderr)
        return None

    ref_gwas = ancestry_data[0]["gwas"]
    ref_gwas = ref_gwas[ref_gwas["RSID"].isin(final_rsids)]
    ref_gwas = ref_gwas.drop_duplicates(subset="RSID")
    ref_gwas = ref_gwas.set_index("RSID").loc[final_rsids].reset_index()

    variant_df = build_multisusie_result_df(fit, final_rsids, ref_gwas, ancestries)
    cs_df = build_credible_sets_df(fit, final_rsids)
    return variant_df, cs_df


_WORKER_CTX = {}


def _init_worker(ctx):
    global _WORKER_CTX
    _WORKER_CTX = ctx


def _run_locus_job(job):
    """Process one locus in a worker process; writes output files."""
    i, n_loci, locus = job
    ctx = _WORKER_CTX
    locus_label = f"chr{locus['chr']}:{locus['start']}-{locus['end']}"

    result = run_multisusie_for_locus(
        locus=locus,
        gwas_dfs=ctx["gwas_dfs"],
        ancestries=ctx["ancestries"],
        sample_sizes=ctx["sample_sizes"],
        thousand_genomes_dir=ctx["thousand_genomes_dir"],
        max_causal=ctx["max_causal"],
        coverage=ctx["coverage"],
        min_abs_corr=ctx["min_abs_corr"],
    )

    if result is None or len(result[0]) == 0:
        return i, n_loci, locus_label, False

    variant_df, cs_df = result
    safe_name = f"{locus['chr']}_{locus['center']}"
    output_dir = ctx["output_dir"]
    out_file = os.path.join(output_dir, f"{safe_name}_finemap.tsv.gz")
    variant_df.to_csv(out_file, sep="\t", index=False, compression="gzip")
    cs_file = os.path.join(output_dir, f"{safe_name}_locus_credible_sets.tsv")
    cs_df.to_csv(cs_file, sep="\t", index=False)
    return i, n_loci, locus_label, True


def main():
    args = parse_args()

    n_ancestries = len(args.ancestries)
    if n_ancestries < 2:
        print("ERROR: MultiSuSiE requires at least 2 ancestries. "
              "Use the standard finemap pipeline for single-ancestry.",
              file=sys.stderr)
        sys.exit(1)

    if len(args.gwas_files) != n_ancestries:
        print("ERROR: Number of GWAS files must match number of ancestries.", file=sys.stderr)
        sys.exit(1)
    if len(args.clumped_files) != n_ancestries:
        print("ERROR: Number of clumped files must match number of ancestries.", file=sys.stderr)
        sys.exit(1)
    if len(args.sample_sizes) != n_ancestries:
        print("ERROR: Number of sample sizes must match number of ancestries.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    window_bp = args.window_kb * 1000

    print(f"Loading {n_ancestries} GWAS datasets...")
    gwas_dfs = [read_gwas(f) for f in args.gwas_files]
    clumped_dfs = [read_clumped(f) for f in args.clumped_files]

    print("Merging lead SNPs across ancestries into loci...")
    loci = get_union_loci(clumped_dfs, window_bp)
    print(f"  {len(loci)} loci to fine-map")

    if len(loci) == 0:
        print("No loci to fine-map.")
        _write_completion(args.completion_file, 0)
        return

    worker_ctx = {
        "gwas_dfs": gwas_dfs,
        "ancestries": args.ancestries,
        "sample_sizes": args.sample_sizes,
        "thousand_genomes_dir": args.thousand_genomes_dir,
        "max_causal": args.max_causal,
        "coverage": args.coverage,
        "min_abs_corr": args.min_abs_corr,
        "output_dir": args.output_dir,
    }
    jobs = [(i + 1, len(loci), locus) for i, locus in enumerate(loci)]
    max_workers = min(
        calculate_parallelism(memory_per_worker_mb=6000),
        len(loci),
    )
    print(f"Running up to {max_workers} loci in parallel...")

    n_success = 0
    n_fail = 0

    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_init_worker,
        initargs=(worker_ctx,),
    ) as executor:
        futures = [executor.submit(_run_locus_job, job) for job in jobs]
        for future in as_completed(futures):
            i, n_loci, locus_label, ok = future.result()
            if ok:
                n_success += 1
                print(f"  [{i}/{n_loci}] Completed {locus_label}")
            else:
                n_fail += 1
                print(f"  [{i}/{n_loci}] Skipped {locus_label} "
                      f"(insufficient shared data or MultiSuSiE failure)")

    print(f"\nMultiSuSiE complete: {n_success} loci succeeded, {n_fail} skipped/failed.")
    _write_completion(args.completion_file, n_success)


def _write_completion(path, n_loci):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(str(n_loci) + "\n")


if __name__ == "__main__":
    main()
