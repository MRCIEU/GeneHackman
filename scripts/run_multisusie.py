"""
Multi-ancestry fine-mapping using MultiSuSiE.

For each clumped locus (union of lead SNPs across ancestries), extracts GWAS
summary statistics and LD matrices per ancestry, then runs MultiSuSiE's
multisusie_rss to produce cross-ancestry credible sets.
"""

import argparse
import gzip
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

import MultiSuSiE


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


def compute_ld_matrix(rsids, chrom, ancestry, thousand_genomes_dir):
    """Compute LD correlation matrix from 1000 Genomes via plink."""
    if len(rsids) < 2:
        return None, None

    bfile = os.path.join(thousand_genomes_dir, ancestry)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".snps", delete=False) as f:
        snp_file = f.name
        f.write("\n".join(rsids) + "\n")

    out_prefix = tempfile.mktemp(prefix="ld_")
    cmd = [
        "plink1.9",
        "--bfile", bfile,
        "--chr", str(chrom),
        "--extract", snp_file,
        "--r", "square",
        "--out", out_prefix
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    ld_file = out_prefix + ".ld"

    if result.returncode != 0 or not os.path.exists(ld_file):
        _cleanup_plink(snp_file, out_prefix)
        return None, None

    bim_file = out_prefix + ".bim"
    if os.path.exists(bim_file):
        bim = pd.read_csv(bim_file, sep="\t", header=None, usecols=[1], names=["SNP"])
        snp_order = bim["SNP"].tolist()
    else:
        bim_ref = bfile + ".bim"
        bim = pd.read_csv(bim_ref, sep="\t", header=None, usecols=[0, 1, 3],
                          names=["CHR", "SNP", "BP"])
        bim = bim[(bim["CHR"].astype(str) == str(chrom)) & (bim["SNP"].isin(rsids))]
        bim = bim.sort_values("BP")
        snp_order = bim["SNP"].tolist()

    R = np.loadtxt(ld_file)
    _cleanup_plink(snp_file, out_prefix)

    if R.shape[0] != len(snp_order):
        return None, None

    return R, snp_order


def _cleanup_plink(snp_file, out_prefix):
    for f in [snp_file] + [out_prefix + ext for ext in
              (".ld", ".nosex", ".log", ".bim", ".bed", ".fam")]:
        try:
            os.unlink(f)
        except OSError:
            pass


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

    b_list = []
    s_list = []
    R_list = []
    n_list = []

    for ad in ancestry_data:
        gwas_sub = ad["gwas"][ad["gwas"]["RSID"].isin(shared_rsids_list)].copy()
        gwas_sub = gwas_sub.drop_duplicates(subset="RSID")
        gwas_sub = gwas_sub.set_index("RSID").loc[shared_rsids_list].reset_index()

        R, snp_order = compute_ld_matrix(
            shared_rsids_list, chrom, ad["ancestry"], thousand_genomes_dir
        )
        if R is None:
            return None

        ld_rsids = set(snp_order)
        common = [s for s in shared_rsids_list if s in ld_rsids]
        if len(common) < 2:
            return None

        ld_idx = [snp_order.index(s) for s in common]
        R_sub = R[np.ix_(ld_idx, ld_idx)]

        gwas_sub = gwas_sub[gwas_sub["RSID"].isin(common)]
        gwas_sub = gwas_sub.set_index("RSID").loc[common].reset_index()

        beta = gwas_sub["BETA"].astype(float).values
        se = gwas_sub["SE"].astype(float).values

        b_list.append(beta)
        s_list.append(se)
        R_list.append(R_sub)
        n_list.append(ad["n"])

    final_rsids = common
    if len(final_rsids) < 2:
        return None

    try:
        fit = MultiSuSiE.multisusie_rss(
            b_list=b_list,
            s_list=s_list,
            R_list=R_list,
            population_sizes=n_list,
            L=max_causal,
            coverage=coverage,
            min_abs_corr=min_abs_corr,
            low_memory_mode=True,
        )
    except Exception as e:
        print(f"  MultiSuSiE failed for locus chr{chrom}:{start}-{end}: {e}", file=sys.stderr)
        return None

    ref_gwas = ancestry_data[0]["gwas"]
    ref_gwas = ref_gwas[ref_gwas["RSID"].isin(final_rsids)]
    ref_gwas = ref_gwas.drop_duplicates(subset="RSID")
    ref_gwas = ref_gwas.set_index("RSID").loc[final_rsids].reset_index()

    result_df = pd.DataFrame({
        "SNP": ref_gwas["SNP"].values if "SNP" in ref_gwas.columns else final_rsids,
        "CHR": ref_gwas["CHR"].values,
        "BP": ref_gwas["BP"].values,
        "RSID": final_rsids,
    })

    if hasattr(fit, "sets") and fit.sets is not None and hasattr(fit.sets, "cs"):
        cs_membership = np.full(len(final_rsids), np.nan)
        for cs_i, cs_snps in enumerate(fit.sets.cs, start=1):
            for idx in cs_snps:
                if idx < len(cs_membership):
                    cs_membership[idx] = cs_i
        result_df["CS"] = cs_membership.astype("Int64")
    else:
        result_df["CS"] = pd.NA

    if hasattr(fit, "pip"):
        result_df["PIP"] = fit.pip

    if hasattr(fit, "lbf_variable") and fit.lbf_variable is not None:
        lbf = fit.lbf_variable
        if lbf.ndim == 2:
            for j in range(lbf.shape[0]):
                result_df[f"LBF_{j+1}"] = lbf[j, :]

    return result_df


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

    n_success = 0
    n_fail = 0

    for i, locus in enumerate(loci):
        locus_label = f"chr{locus['chr']}:{locus['start']}-{locus['end']}"
        print(f"  [{i+1}/{len(loci)}] Processing {locus_label}...")

        result = run_multisusie_for_locus(
            locus=locus,
            gwas_dfs=gwas_dfs,
            ancestries=args.ancestries,
            sample_sizes=args.sample_sizes,
            thousand_genomes_dir=args.thousand_genomes_dir,
            max_causal=args.max_causal,
            coverage=args.coverage,
            min_abs_corr=args.min_abs_corr,
        )

        if result is not None and len(result) > 0:
            safe_name = f"{locus['chr']}_{locus['center']}_finemap.tsv.gz"
            out_file = os.path.join(args.output_dir, safe_name)
            result.to_csv(out_file, sep="\t", index=False, compression="gzip")
            n_success += 1
        else:
            print(f"    Skipped {locus_label} (insufficient shared data or MultiSuSiE failure)")
            n_fail += 1

    print(f"\nMultiSuSiE complete: {n_success} loci succeeded, {n_fail} skipped/failed.")
    _write_completion(args.completion_file, n_success)


def _write_completion(path, n_loci):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(str(n_loci) + "\n")


if __name__ == "__main__":
    main()
