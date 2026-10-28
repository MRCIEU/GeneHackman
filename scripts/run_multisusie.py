"""
Multi-ancestry fine-mapping using MultiSuSiE.

For each clumped locus (union of lead SNPs across ancestries), extracts GWAS
summary statistics and LD matrices per ancestry, then runs MultiSuSiE's
multisusie_rss to produce cross-ancestry credible sets.
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd

import MultiSuSiE

from ld_reference import LDPack, make_positive_semidefinite


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
                        help="Ancestry codes (e.g. EUR EAS); pan_ukb also allows CSA (or SAS) and MID")
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
    parser.add_argument("--thousand_genomes_dir", default=None,
                        help="Path to 1000 Genomes reference panel directory (ld_reference 1000genomes)")
    parser.add_argument("--ld_reference", choices=["1000genomes", "pan_ukb"], default="1000genomes",
                        help="LD source: 1000 Genomes genotypes via PLINK, or the precomputed Pan-UKB LD pack")
    parser.add_argument("--ld_pack_dir", default=None,
                        help="Pan-UKB LD pack directory (ld_reference pan_ukb)")
    parser.add_argument("--build", choices=["GRCh37", "GRCh38"], default="GRCh37",
                        help="Reference build of the standardised GWAS (selects the pack's SNP keys)")
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
            "--keep-allele-order",
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


def _iter_credible_sets(sets, n_variants):
    """Yield non-empty credible sets from a MultiSuSiE fit.sets tuple."""
    if sets is None or not isinstance(sets, (list, tuple)) or len(sets) == 0:
        return

    cs_list = sets[0]
    purity_list = sets[1] if len(sets) > 1 else [np.nan] * len(cs_list)
    coverage_list = sets[2] if len(sets) > 2 else [np.nan] * len(cs_list)
    pass_filter = sets[3] if len(sets) > 3 else [True] * len(cs_list)

    cs_counter = 0
    for ser_idx, indices in enumerate(cs_list):
        if len(indices) == 0:
            continue
        indices = np.asarray(indices, dtype=int)
        indices = indices[(indices >= 0) & (indices < n_variants)]
        if len(indices) == 0:
            continue
        cs_counter += 1
        yield {
            "cs_id": cs_counter,
            "ser_idx": ser_idx,
            "indices": indices,
            "purity": purity_list[ser_idx] if ser_idx < len(purity_list) else np.nan,
            "coverage": coverage_list[ser_idx] if ser_idx < len(coverage_list) else np.nan,
            "pass_filter": bool(pass_filter[ser_idx]) if ser_idx < len(pass_filter) else True,
        }


def _fit_pip(fit, n_variants):
    if hasattr(fit, "pip") and fit.pip is not None and len(fit.pip) == n_variants:
        return np.asarray(fit.pip, dtype=float)
    return np.full(n_variants, np.nan)


def _fit_lbf(fit, n_sers):
    if hasattr(fit, "lbf") and fit.lbf is not None:
        return np.asarray(fit.lbf, dtype=float)
    return np.full(n_sers, np.nan)


def _compute_cs_purity(indices, r_matrix):
    """Minimum absolute pairwise LD correlation within a credible set."""
    indices = np.asarray(indices, dtype=int)
    if len(indices) < 2:
        return 1.0
    sub = r_matrix[np.ix_(indices, indices)]
    n = len(indices)
    min_corr = 1.0
    for i in range(n):
        for j in range(i + 1, n):
            min_corr = min(min_corr, abs(sub[i, j]))
    return float(min_corr)


def _append_variant_ancestry_columns(result_df, fit, ancestry_gwas, ancestries, variant_indices):
    """Add per-population GWAS stats and MultiSuSiE posterior effect estimates."""
    for k, ancestry in enumerate(ancestries):
        gwas_sub = ancestry_gwas[k]
        result_df[f"BETA_{ancestry}"] = gwas_sub["BETA"].iloc[variant_indices].astype(float).values
        result_df[f"SE_{ancestry}"] = gwas_sub["SE"].iloc[variant_indices].astype(float).values
        if "P" in gwas_sub.columns:
            result_df[f"P_{ancestry}"] = gwas_sub["P"].iloc[variant_indices].astype(float).values
        else:
            result_df[f"P_{ancestry}"] = np.nan

        if hasattr(fit, "coef") and fit.coef is not None and fit.coef.shape[0] > k:
            result_df[f"COEF_{ancestry}"] = fit.coef[k][variant_indices]
        else:
            result_df[f"COEF_{ancestry}"] = np.nan

        if hasattr(fit, "coef_sd") and fit.coef_sd is not None and fit.coef_sd.shape[0] > k:
            result_df[f"COEF_SD_{ancestry}"] = fit.coef_sd[k][variant_indices]
        else:
            result_df[f"COEF_SD_{ancestry}"] = np.nan


def build_locus_credible_sets_df(fit, final_rsids, r_list, ancestries, snp_keys=None):
    """One row per credible set for locus_credible_sets.tsv.

    final_rsids are the variant labels reported as TOP_SNP_RSID (RSIDs, or SNP keys when the GWAS has
    none); snp_keys, when given, are the GeneHackman SNP keys reported as TOP_SNP.
    """
    sets = getattr(fit, "sets", None)
    n = len(final_rsids)
    pip = _fit_pip(fit, n)
    n_sers = len(sets[0]) if sets is not None and len(sets) > 0 else 0
    lbf = _fit_lbf(fit, n_sers)

    rows = []
    for cs in _iter_credible_sets(sets, n):
        indices = cs["indices"]
        cs_pip = pip[indices]
        top_local = int(np.nanargmax(cs_pip))
        top_idx = int(indices[top_local])
        row = {
            "CS_ID": cs["cs_id"],
            "CS_SIZE": len(indices),
            "TOP_SNP": snp_keys[top_idx] if snp_keys is not None else final_rsids[top_idx],
            "TOP_SNP_RSID": final_rsids[top_idx],
            "TOP_SNP_PIP": float(pip[top_idx]),
            "TOP_SNP_LBF": (
                float(lbf[cs["ser_idx"]]) if cs["ser_idx"] < len(lbf) else np.nan
            ),
            "COVERAGE": cs["coverage"],
            "PURITY": cs["purity"],
            "PASS_FILTER": cs["pass_filter"],
        }
        for k, ancestry in enumerate(ancestries):
            row[f"PURITY_{ancestry}"] = _compute_cs_purity(indices, r_list[k])
        rows.append(row)

    if not rows:
        return pd.DataFrame()

    cs_df = pd.DataFrame(rows)
    purity_cols = [f"PURITY_{ancestry}" for ancestry in ancestries]
    return cs_df[
        [
            "CS_ID",
            "CS_SIZE",
            "TOP_SNP",
            "TOP_SNP_RSID",
            "TOP_SNP_PIP",
            "TOP_SNP_LBF",
            "COVERAGE",
            "PURITY",
            "PASS_FILTER",
            *purity_cols,
        ]
    ]


def build_locus_credible_set_variants_df(
    fit, final_rsids, ref_gwas, ancestry_gwas, ancestries
):
    """One row per variant in any credible set for locus_credible_set_variants.tsv."""
    sets = getattr(fit, "sets", None)
    n = len(final_rsids)
    pip = _fit_pip(fit, n)

    rows = []
    for cs in _iter_credible_sets(sets, n):
        for idx in cs["indices"]:
            rows.append({
                "CS_ID": cs["cs_id"],
                "SNP": ref_gwas["SNP"].iloc[idx] if "SNP" in ref_gwas.columns else final_rsids[idx],
                "RSID": final_rsids[idx],
                "CHROM": ref_gwas["CHR"].iloc[idx],
                "POS": ref_gwas["BP"].iloc[idx],
                "ALLELE_A1": ref_gwas["EA"].iloc[idx],
                "ALLELE_A2": ref_gwas["OA"].iloc[idx],
                "GLOBAL_PIP": float(pip[idx]),
                "_variant_idx": idx,
            })

    if not rows:
        return pd.DataFrame()

    variants_df = pd.DataFrame(rows)
    _append_variant_ancestry_columns(
        variants_df,
        fit,
        ancestry_gwas,
        ancestries,
        variants_df["_variant_idx"].to_numpy(dtype=int),
    )
    variants_df = variants_df.drop(columns="_variant_idx")
    coef_cols = [f"COEF_{ancestry}" for ancestry in ancestries]
    coef_sd_cols = [f"COEF_SD_{ancestry}" for ancestry in ancestries]
    beta_cols = [f"BETA_{ancestry}" for ancestry in ancestries]
    se_cols = [f"SE_{ancestry}" for ancestry in ancestries]
    p_cols = [f"P_{ancestry}" for ancestry in ancestries]
    return variants_df[
        [
            "CS_ID",
            "SNP",
            "RSID",
            "CHROM",
            "POS",
            "ALLELE_A1",
            "ALLELE_A2",
            "GLOBAL_PIP",
            *coef_cols,
            *coef_sd_cols,
            *beta_cols,
            *se_cols,
            *p_cols,
        ]
    ]


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


def compute_ld_from_reference(ids, chrom, ancestry, ld_options):
    """LD matrix and the ids it covers, from the configured LD reference."""
    if ld_options["ld_reference"] == "pan_ukb":
        pack = ld_options.get("pack")
        if pack is None:
            pack = ld_options["pack"] = LDPack(ld_options["ld_pack_dir"])
        r_matrix, found = pack.get_ld_matrix(ids, chrom, ancestry, ld_options["build"])
        if len(found) < 2:
            return None, None
        return r_matrix, found
    return compute_ld_matrix(ids, chrom, ancestry, ld_options["thousand_genomes_dir"])


def run_multisusie_for_locus(locus, gwas_dfs, ancestries, sample_sizes,
                             ld_options, max_causal, coverage,
                             min_abs_corr):
    """Run MultiSuSiE on a single locus across ancestries.

    Variants are matched across ancestries and to the LD reference on RSID (1000genomes) or on the
    GeneHackman SNP key (pan_ukb, which needs no RSIDs).
    """
    chrom = locus["chr"]
    start = locus["start"]
    end = locus["end"]
    id_col = "SNP" if ld_options["ld_reference"] == "pan_ukb" else "RSID"

    ancestry_data = []
    for i, (gwas_df, ancestry, n) in enumerate(zip(gwas_dfs, ancestries, sample_sizes)):
        locus_gwas = extract_locus_data(gwas_df, chrom, start, end)
        if len(locus_gwas) < 2 or id_col not in locus_gwas.columns:
            return None
        locus_gwas = locus_gwas[locus_gwas[id_col].notna() & (locus_gwas[id_col].astype(str).str.len() > 0)]
        if id_col == "SNP":
            locus_gwas = locus_gwas.assign(SNP=locus_gwas["SNP"].astype(str).str.upper())
        if len(locus_gwas) < 2:
            return None
        ancestry_data.append({
            "gwas": locus_gwas,
            "ancestry": ancestry,
            "n": n
        })

    shared_ids = set(ancestry_data[0]["gwas"][id_col])
    for ad in ancestry_data[1:]:
        shared_ids &= set(ad["gwas"][id_col])

    if len(shared_ids) < 2:
        return None

    shared_ids_list = sorted(shared_ids)

    ancestry_ld = []
    for ad in ancestry_data:
        gwas_sub = ad["gwas"][ad["gwas"][id_col].isin(shared_ids_list)].copy()
        gwas_sub = gwas_sub.drop_duplicates(subset=id_col)
        gwas_sub = gwas_sub.set_index(id_col).loc[shared_ids_list].reset_index()

        R, snp_order = compute_ld_from_reference(shared_ids_list, chrom, ad["ancestry"], ld_options)
        if R is None:
            return None

        ancestry_ld.append({
            "gwas_sub": gwas_sub,
            "R": R,
            "snp_order": snp_order,
            "ancestry": ad["ancestry"],
            "n": ad["n"],
        })

    final_ids = set(shared_ids_list)
    for entry in ancestry_ld:
        final_ids &= set(entry["snp_order"])
    final_ids = sorted(final_ids)
    if len(final_ids) < 2:
        return None

    b_list = []
    s_list = []
    R_list = []
    n_list = []
    varY_list = []

    locus_label = f"chr{chrom}:{start}-{end}"
    for entry in ancestry_ld:
        snp_order = entry["snp_order"]
        snp_to_idx = {variant: idx for idx, variant in enumerate(snp_order)}
        ld_idx = [snp_to_idx[r] for r in final_ids]
        R_sub = np.ascontiguousarray(entry["R"][np.ix_(ld_idx, ld_idx)], dtype=np.float64)
        if ld_options["ld_reference"] == "pan_ukb":
            # pairs beyond the pack radius are r = 0 and r is quantised; keep R positive semi-definite
            R_sub, _ = make_positive_semidefinite(R_sub)

        gwas_sub = entry["gwas_sub"].set_index(id_col).reindex(final_ids)
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
            variant_ids=final_ids,
            low_memory_mode=False,
            float_type=np.float64,
        )
    except Exception as e:
        print(f"  MultiSuSiE failed for locus chr{chrom}:{start}-{end}: {e}", file=sys.stderr)
        return None

    ref_gwas = ancestry_data[0]["gwas"]
    ref_gwas = ref_gwas[ref_gwas[id_col].isin(final_ids)]
    ref_gwas = ref_gwas.drop_duplicates(subset=id_col)
    ref_gwas = ref_gwas.set_index(id_col).loc[final_ids].reset_index()

    ancestry_gwas = []
    for entry in ancestry_ld:
        gwas_sub = entry["gwas_sub"].set_index(id_col).reindex(final_ids)
        ancestry_gwas.append(gwas_sub.reset_index())

    # report RSIDs where the GWAS has them, whichever id the LD was matched on
    labels = list(final_ids)
    if id_col == "SNP" and "RSID" in ref_gwas.columns:
        labels = [r if isinstance(r, str) and r else s for r, s in zip(ref_gwas["RSID"], ref_gwas["SNP"])]
    snp_keys = list(ref_gwas["SNP"]) if "SNP" in ref_gwas.columns else None
    cs_sets_df = build_locus_credible_sets_df(fit, labels, R_list, ancestries, snp_keys=snp_keys)
    cs_variants_df = build_locus_credible_set_variants_df(
        fit, labels, ref_gwas, ancestry_gwas, ancestries
    )
    return cs_sets_df, cs_variants_df


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
        ld_options=ctx["ld_options"],
        max_causal=ctx["max_causal"],
        coverage=ctx["coverage"],
        min_abs_corr=ctx["min_abs_corr"],
    )

    if result is None:
        return i, n_loci, locus_label, False

    cs_sets_df, cs_variants_df = result
    if len(cs_sets_df) == 0:
        return i, n_loci, locus_label, False
    safe_name = f"{locus['chr']}_{locus['center']}"
    output_dir = ctx["output_dir"]
    cs_file = os.path.join(output_dir, f"{safe_name}_locus_credible_sets.tsv")
    cs_variants_file = os.path.join(
        output_dir, f"{safe_name}_locus_credible_set_variants.tsv"
    )
    cs_sets_df.to_csv(cs_file, sep="\t", index=False)
    cs_variants_df.to_csv(cs_variants_file, sep="\t", index=False)
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

    if args.ld_reference == "pan_ukb" and not args.ld_pack_dir:
        print("ERROR: --ld_pack_dir is required with --ld_reference pan_ukb", file=sys.stderr)
        sys.exit(1)
    if args.ld_reference == "1000genomes" and not args.thousand_genomes_dir:
        print("ERROR: --thousand_genomes_dir is required with --ld_reference 1000genomes", file=sys.stderr)
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
        "ld_options": {"ld_reference": args.ld_reference, "thousand_genomes_dir": args.thousand_genomes_dir,
                       "ld_pack_dir": args.ld_pack_dir, "build": args.build},
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
