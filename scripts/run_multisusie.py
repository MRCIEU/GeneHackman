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
    parser.add_argument("--sample_sizes", nargs="+", type=float, required=True,
                        help="GWAS sample sizes (one per ancestry)")
    parser.add_argument("--output_dir", required=True,
                        help="Output directory for per-locus results")
    parser.add_argument("--window_kb", type=int, default=1000,
                        help="Half-width of fine-mapping window in kb")
    parser.add_argument("--max_n_ratio", type=float, default=10,
                        help="For --max_causal auto: largest GWAS N / LD reference N that still uses L = 10")
    parser.add_argument("--max_region_kb", type=float, default=3000,
                        help="Maximum merged region width, or tile width, in kb")
    parser.add_argument("--region_strategy", choices=REGION_STRATEGIES, default="merge",
                        help="Fine-mapping regions: merge, none or tile (as define_finemap_regions in R/finemap.R)")
    parser.add_argument("--ld_bfiles", nargs="+", default=None,
                        help="Per ancestry: PLINK prefix of in-sample/user genotypes for LD, or 'none'")
    parser.add_argument("--remove_ld_outliers", action="store_true",
                        help="Remove variants flagged by the LD-mismatch check (any ancestry) before fitting")
    parser.add_argument("--max_causal", default="auto",
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
    parser.add_argument("--build", choices=["GRCh37", "GRCh38"], default="GRCh38",
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
        # drop rows without a lead SNP or position (e.g. a row of NAs from trailing blank lines); such a row makes
        # pandas read CHR and BP as floats, so restore the types the GWAS use (CHR as a string, "2" not "2.0")
        df = df[["SNP", "CHR", "BP"]].dropna().reset_index(drop=True)
        df["CHR"] = df["CHR"].astype(str).str.replace(r"\.0$", "", regex=True)
        df["BP"] = df["BP"].astype("int64")
        return df
    except Exception:
        return pd.DataFrame(columns=["SNP", "CHR", "BP"])


REGION_STRATEGIES = ["merge", "none", "tile"]
LD_MISMATCH_S_THRESHOLD = 0.1  # as ld_mismatch_s_threshold in R/finemap.R
MHC_BOUNDS = {"GRCh38": (28510120, 33480577), "GRCh37": (28477797, 33448354)}


def _pve_adjusted_z(z, n):
    """z rescaled for the variance it explains, as susieR's estimate_s_rss() and kriging_rss() do."""
    z = np.nan_to_num(np.asarray(z, dtype=float))
    return np.sqrt((n - 1) / (z ** 2 + n - 2)) * z


def estimate_s_rss(z, r_matrix, n, r_tol=1e-8):
    """Port of susieR::estimate_s_rss(method = "null-mle"): the s that maximises the likelihood of
    z ~ N(0, (1 - s) R + s I).  0 when the z-scores and the LD are consistent."""
    from scipy.optimize import minimize_scalar
    z = _pve_adjusted_z(z, n)
    values, vectors = np.linalg.eigh(r_matrix)
    values = np.where(values < r_tol, 0.0, values)
    ztv = vectors.T @ z

    def negloglik(s):
        d = (1 - s) * values + s
        return 0.5 * np.sum(np.log(d)) + 0.5 * np.sum(ztv ** 2 / d)
    return float(minimize_scalar(negloglik, bounds=(0, 1), method="bounded").x)


def _mixture_weights(log_lik, iterations=500, tol=1e-8):
    """Maximum-likelihood mixture proportions by EM (susieR uses mixsqp for the same convex problem)."""
    lik = np.exp(log_lik)
    w = np.full(lik.shape[1], 1.0 / lik.shape[1])
    for _ in range(iterations):
        post = lik * w
        post /= post.sum(axis=1, keepdims=True)
        new = post.mean(axis=0)
        if np.max(np.abs(new - w)) < tol:
            return new
        w = new
    return w


def kriging_rss(z, r_matrix, n, s=None, r_tol=1e-8):
    """Port of susieR::kriging_rss(): each z-score's expectation given the others (condmean, condvar) and the
    log likelihood ratio for an allele flip (logLR).  Returns (adjusted z, condmean, condvar, logLR)."""
    from scipy.stats import norm
    if s is None:
        s = estimate_s_rss(z, r_matrix, n, r_tol)
    s = min(max(s, 0.0), 0.8 if s > 1 else 1.0)
    z = _pve_adjusted_z(z, n)
    values, vectors = np.linalg.eigh(r_matrix)
    values = np.where(values < r_tol, 0.0, values)
    with np.errstate(divide="ignore"):
        dinv = 1 / ((1 - s) * values + s)
    dinv[~np.isfinite(dinv)] = 0
    precision = (vectors * dinv) @ vectors.T
    diag = np.diag(precision)
    condmean = -(precision @ z - diag * z) / diag
    condvar = 1 / diag
    z_std_diff = (z - condmean) / np.sqrt(condvar)
    a_min = 0.8
    a_max = 2.0 if np.max(z_std_diff ** 2) < 1 else 2 * np.sqrt(np.max(z_std_diff ** 2))
    npoint = int(np.ceil(np.log2(a_max / a_min) / np.log2(1.05)))
    a_grid = 1.05 ** np.arange(-npoint, 1) * a_max
    sd = np.outer(np.sqrt(condvar), a_grid)

    def log_mixture(x, w=None):
        llik = norm.logpdf(x[:, None], scale=sd)
        factors = llik.max(axis=1, keepdims=True)
        llik = llik - factors
        if w is None:
            return llik, factors
        return np.log(np.exp(llik) @ (w + 1e-15)) + factors[:, 0]
    llik0, _ = log_mixture(z - condmean)
    w = _mixture_weights(llik0)
    log_lr = log_mixture(z + condmean, w) - log_mixture(z - condmean, w)
    return z, condmean, condvar, log_lr


def ld_mismatch_qc(z, r_matrix, n):
    """As ld_mismatch_qc() in R/ld_reference.R: (s, outlier mask with logLR > 2 and |z| > 2)."""
    z = np.asarray(z, dtype=float)
    if len(z) < 3:
        return float("nan"), np.zeros(len(z), dtype=bool)
    try:
        s = estimate_s_rss(z, r_matrix, n)
        z_adj, _, _, log_lr = kriging_rss(z, r_matrix, n, s)
        return s, (log_lr > 2) & (np.abs(z_adj) > 2)
    except Exception as e:  # diagnostics must never stop fine-mapping
        print(f"  LD mismatch check failed: {e}", file=sys.stderr)
        return float("nan"), np.zeros(len(z), dtype=bool)


def ld_mismatch_flag(s, n_outliers_region, n_outliers_cs=0, mhc=False):
    """As ld_mismatch_flag() in R/finemap.R: OUTLIER_IN_CS, LD_MISMATCH, OUTLIERS_IN_REGION, MHC, or OK."""
    flags = []
    in_cs = n_outliers_cs is not None and n_outliers_cs > 0
    if in_cs:
        flags.append("OUTLIER_IN_CS")
    if s is not None and s == s and s > LD_MISMATCH_S_THRESHOLD:
        flags.append("LD_MISMATCH")
    if not in_cs and n_outliers_region:
        flags.append("OUTLIERS_IN_REGION")
    if mhc:
        flags.append("MHC")
    return ";".join(flags) if flags else "OK"


def overlaps_mhc(chrom, start, end, build="GRCh38"):
    low, high = MHC_BOUNDS.get(build, MHC_BOUNDS["GRCh38"])
    return str(chrom) == "6" and start <= high and end >= low


def define_regions(leads, window_bp, max_region_bp=3_000_000, strategy="merge"):
    """Fine-mapping regions from lead SNPs: the same rules as define_finemap_regions() in R/finemap.R.

    merge: overlapping windows are merged, up to max_region_bp; none: one window per lead; tile: each group of
    overlapping windows is covered by tiles max_region_bp wide every max_region_bp / 3 (Weissbrod et al. 2020),
    keeping tiles that contain a lead.  Returns dicts with chr, start, end, center, leads and cluster."""
    if strategy not in REGION_STRATEGIES:
        raise ValueError(f"region strategy must be one of {REGION_STRATEGIES}, not {strategy}")
    if len(leads) == 0:
        return []
    max_region_bp = max(max_region_bp, 2 * window_bp)
    leads = leads.assign(CHR=leads["CHR"].astype(str), BP=leads["BP"].astype(float))
    leads = leads.sort_values(["CHR", "BP"], key=lambda c: c.map(_chrom_order) if c.name == "CHR" else c)
    leads = leads.reset_index(drop=True)
    regions = []

    def add(chrom, start, end, rows, cluster=None):
        regions.append({"chr": chrom, "start": int(start), "end": int(end), "center": int((start + end) / 2),
                        "leads": list(leads.loc[rows, "SNP"]), "cluster": cluster})

    if strategy == "none":
        for i, row in leads.iterrows():
            add(row["CHR"], max(0, row["BP"] - window_bp), row["BP"] + window_bp, [i])
        return regions

    cap = max_region_bp if strategy == "merge" else float("inf")
    groups = []  # [chr, start, end, rows]
    for i, row in leads.iterrows():
        chrom, bp = row["CHR"], row["BP"]
        if not groups or chrom != groups[-1][0] or bp - window_bp > groups[-1][2]:
            groups.append([chrom, max(0, bp - window_bp), bp + window_bp, [i]])
            continue
        g = groups[-1]
        new_end = max(g[2], bp + window_bp)
        if new_end - g[1] <= cap:
            g[2] = new_end
            g[3].append(i)
        elif bp <= g[2]:
            g[2] = max(g[2], g[1] + cap)
            g[3].append(i)
        else:
            groups.append([chrom, g[2] + 1, bp + window_bp, [i]])
    for k, (chrom, start, end, rows) in enumerate(groups, start=1):
        if strategy == "merge" or end - start <= max_region_bp:
            add(chrom, start, end, rows, k if strategy == "tile" else None)
            continue
        step = max_region_bp / 3
        tile_start = start
        while True:
            tile_end = min(tile_start + max_region_bp, end)
            inside = [r for r in rows if tile_start <= leads.loc[r, "BP"] <= tile_end]
            if inside:
                add(chrom, tile_start, tile_end, inside, k)
            if tile_start + max_region_bp >= end:
                break
            tile_start += step
    return regions


def _chrom_order(chrom):
    chrom = str(chrom)
    return int(chrom) if chrom.isdigit() else {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}.get(chrom.upper(), 100)


def tile_core_mask(positions, regions, index):
    """True for positions whose nearest tile centre (within the region's cluster) is region `index`."""
    positions = np.asarray(positions, dtype=float)
    cluster = regions[index].get("cluster")
    if cluster is None:
        return np.ones(len(positions), dtype=bool)
    same = [j for j, r in enumerate(regions) if r.get("cluster") == cluster and r["chr"] == regions[index]["chr"]]
    if len(same) <= 1:
        return np.ones(len(positions), dtype=bool)
    centres = np.array([(regions[j]["start"] + regions[j]["end"]) / 2 for j in same])
    nearest = np.array(same)[np.abs(positions[:, None] - centres[None, :]).argmin(axis=1)]
    return nearest == index


def get_union_loci(clumped_dfs, window_bp, max_region_bp=3_000_000, strategy="merge"):
    """Merge lead SNPs across ancestries into fine-mapping regions (see define_regions)."""
    all_leads = pd.concat(clumped_dfs, ignore_index=True)
    if len(all_leads) == 0:
        return []
    return define_regions(all_leads[["SNP", "CHR", "BP"]], window_bp, max_region_bp, strategy)


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

        # PLINK writes the matrix in .bim order; keep the panel's alleles to orient r to the GWAS effect allele
        bim = pd.read_csv(bfile + ".bim", sep=r"\s+", header=None, usecols=[0, 1, 4, 5],
                          names=["CHR", "SNP", "A1", "A2"], dtype=str)
        bim = bim[(bim["CHR"] == str(chrom)) & (bim["SNP"].isin(rsids))]
        snp_order = bim["SNP"].tolist()

        R = _read_plink_ld_matrix(ld_file, len(snp_order))
        if R is None:
            return None, None

        return R, snp_order, list(bim["A1"].str.upper()), list(bim["A2"].str.upper())
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


def build_locus_credible_sets_df(fit, final_rsids, r_list, ancestries, snp_keys=None, outlier=None):
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
            "N_LD_OUTLIERS_IN_CS": int(np.sum(outlier[indices])) if outlier is not None else 0,
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
            "N_LD_OUTLIERS_IN_CS",
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


def compute_ld_from_bfile(snp_keys, chrom, bfile):
    """LD from user-supplied PLINK genotypes, matched on position and alleles (as compute_ld_matrix_from_bfile
    in R/ld_reference.R).  r is returned for the GWAS effect allele (the first allele of each SNP key)."""
    bim = pd.read_csv(bfile + ".bim", sep=r"\s+", header=None, names=["CHR", "ID", "CM", "BP", "A1", "A2"], dtype=str)
    bim["ROW"] = np.arange(len(bim))
    bim = bim[bim["CHR"].str.replace("^chr", "", regex=True, case=False) == str(chrom)]
    key_a1 = (str(chrom) + ":" + bim["BP"] + "_" + bim["A1"].str.upper() + "_" + bim["A2"].str.upper()).tolist()
    key_a2 = (str(chrom) + ":" + bim["BP"] + "_" + bim["A2"].str.upper() + "_" + bim["A1"].str.upper()).tolist()
    count_a1 = pd.Series(key_a1).value_counts()
    count_a2 = pd.Series(key_a2).value_counts()
    index_a1 = {k: i for i, k in enumerate(key_a1)}
    index_a2 = {k: i for i, k in enumerate(key_a2)}
    rows, signs, keys = [], [], []
    for key in dict.fromkeys(str(k).upper() for k in snp_keys):
        in_a1, in_a2 = key in index_a1, key in index_a2
        if in_a1 == in_a2 or count_a1.get(key, 0) > 1 or count_a2.get(key, 0) > 1:
            continue  # no match, or ambiguous
        rows.append(index_a1[key] if in_a1 else index_a2[key])
        signs.append(1.0 if in_a1 else -1.0)
        keys.append(key)
    if len(keys) < 2:
        return None, None
    order = np.argsort(bim["ROW"].values[rows])
    ids = bim["ID"].values[rows][order]
    if len(set(bim["ID"])) != len(bim) or "." in set(ids):
        raise ValueError(f"Variant IDs in {bfile}.bim must be unique (and not '.') to extract LD")
    keys = [keys[i] for i in order]
    signs = np.array(signs)[order]
    tmpdir = tempfile.mkdtemp(prefix="user_ld_")
    try:
        snp_file = os.path.join(tmpdir, "snps.txt")
        with open(snp_file, "w", encoding="utf-8") as f:
            f.write("\n".join(ids) + "\n")
        out_prefix = os.path.join(tmpdir, "ld")
        result = subprocess.run(["plink1.9", "--bfile", bfile, "--extract", snp_file, "--r", "square",
                                 "--keep-allele-order", "--out", out_prefix], capture_output=True, text=True)
        if result.returncode != 0 or not os.path.exists(out_prefix + ".ld"):
            return None, None
        r_matrix = _read_plink_ld_matrix(out_prefix + ".ld", len(keys))
        if r_matrix is None:
            return None, None
        r_matrix = np.nan_to_num(r_matrix * np.outer(signs, signs))
        return r_matrix, keys
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def align_to_gwas_alleles(r_matrix, gwas_sub, a1, a2):
    """Flip r where the panel's A1 is the GWAS other allele; mask variants whose alleles don't match."""
    ea = gwas_sub["EA"].astype(str).str.upper().to_numpy()
    oa = gwas_sub["OA"].astype(str).str.upper().to_numpy()
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)
    same = (ea == a1) & (oa == a2)
    swapped = (ea == a2) & (oa == a1)
    signs = np.where(swapped, -1.0, 1.0)
    return r_matrix * np.outer(signs, signs), same | swapped


def ld_reference_n_samples(ancestry, ld_options, bfile=None):
    """Number of individuals behind the LD: the .fam of a user file or 1000G panel, or the pack's n_samples."""
    def count_fam(prefix):
        try:
            with open(prefix + ".fam", encoding="utf-8") as f:
                return float(sum(1 for _ in f))
        except OSError:
            return float("nan")
    if bfile:
        return count_fam(bfile)
    if ld_options["ld_reference"] == "pan_ukb":
        try:
            return float(LDPack(ld_options["ld_pack_dir"]).info(ancestry)["n_samples"])
        except Exception:
            return float("nan")
    return count_fam(os.path.join(ld_options["thousand_genomes_dir"] or "", ancestry))


def choose_max_causal(max_causal, sample_sizes, reference_sizes, max_n_ratio=10, default_l=10):
    """SuSiE L, as choose_max_causal() in R/ld_reference.R: with 'auto', L = 1 when any ancestry's GWAS N is more
    than max_n_ratio times its LD reference N, else default_l.  Returns (L, rule)."""
    if str(max_causal).lower() != "auto":
        return int(max_causal), "set by user"
    ratios = [n / ref for n, ref in zip(sample_sizes, reference_sizes) if ref and ref == ref and ref > 0]
    if len(ratios) < len(sample_sizes):
        return default_l, "auto: LD reference size unknown"
    if max(ratios) > max_n_ratio:
        return 1, f"auto: max N / N_REF = {max(ratios):.1f} > {max_n_ratio:g}"
    return default_l, f"auto: max N / N_REF = {max(ratios):.1f} <= {max_n_ratio:g}"


def locus_sample_size(locus_gwas, default_n):
    """Median per-variant N_EFF (binary traits) or N in the locus, else the input sample size."""
    for column in ("N_EFF", "N"):
        if column in locus_gwas.columns:
            n = pd.to_numeric(locus_gwas[column], errors="coerce")
            n = n[n > 0]
            if len(n):
                return float(n.median())
    return float(default_n)


def compute_ld_from_reference(ids, chrom, ancestry, ld_options, bfile=None):
    """LD matrix and the ids it covers, from the configured LD reference (or a user PLINK file).

    For 1000 Genomes the panel's alleles are returned as well: (matrix, ids, a1, a2)."""
    if bfile:
        return compute_ld_from_bfile(ids, chrom, bfile)
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
                             min_abs_corr, remove_ld_outliers=False):
    """Run MultiSuSiE on a single locus across ancestries.

    Variants are matched across ancestries and to the LD reference on RSID (1000genomes) or on the
    GeneHackman SNP key (pan_ukb, which needs no RSIDs).
    """
    chrom = locus["chr"]
    start = locus["start"]
    end = locus["end"]
    bfiles = ld_options.get("ld_bfiles") or [None] * len(ancestries)
    uses_keys = ld_options["ld_reference"] == "pan_ukb" or any(bfiles)
    id_col = "SNP" if uses_keys else "RSID"

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
            "n": locus_sample_size(locus_gwas, n),
            "bfile": bfiles[i],
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

        ld = compute_ld_from_reference(shared_ids_list, chrom, ad["ancestry"], ld_options, ad["bfile"])
        R, snp_order = ld[0], ld[1]
        if R is None:
            return None
        if len(ld) == 4:
            # matched on RSID: orient r to the GWAS effect allele and drop variants whose alleles differ
            by_id = gwas_sub.set_index(id_col).reindex(snp_order)
            R, matched = align_to_gwas_alleles(R, by_id, ld[2], ld[3])
            snp_order = [s for s, keep in zip(snp_order, matched) if keep]
            R = R[np.ix_(matched, matched)]

        ancestry_ld.append({
            "gwas_sub": gwas_sub,
            "R": R,
            "snp_order": snp_order,
            "ancestry": ad["ancestry"],
            "n": ad["n"],
            "bfile": ad["bfile"],
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
    eaf_list = []

    locus_label = f"chr{chrom}:{start}-{end}"
    for entry in ancestry_ld:
        snp_order = entry["snp_order"]
        snp_to_idx = {variant: idx for idx, variant in enumerate(snp_order)}
        ld_idx = [snp_to_idx[r] for r in final_ids]
        R_sub = np.ascontiguousarray(entry["R"][np.ix_(ld_idx, ld_idx)], dtype=np.float64)
        if ld_options["ld_reference"] == "pan_ukb" and not entry.get("bfile"):
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
        eaf_list.append(eaf)

    if not _validate_multisusie_inputs(b_list, s_list, R_list, ancestries, locus_label):
        return None

    # LD mismatch check per ancestry (each against its own LD reference)
    qc = [ld_mismatch_qc(b / s_, r, n) for b, s_, r, n in zip(b_list, s_list, R_list, n_list)]
    outlier = np.zeros(len(final_ids), dtype=bool)
    for _, mask in qc:
        outlier |= mask
    region = {"N_VARIANTS": len(final_ids), "N_LD_OUTLIERS": int(outlier.sum()),
              "LD_OUTLIERS_REMOVED": bool(remove_ld_outliers)}
    for ancestry, (s_value, mask) in zip(ancestries, qc):
        region[f"LD_MISMATCH_S_{ancestry}"] = s_value
        region[f"N_LD_OUTLIERS_{ancestry}"] = int(mask.sum())
    if remove_ld_outliers and outlier.any():
        keep = ~outlier
        if keep.sum() < 2:
            return None
        print(f"  {locus_label}: removing {int(outlier.sum())} LD outlier(s) before fitting", file=sys.stderr)
        b_list = [b[keep] for b in b_list]
        s_list = [x[keep] for x in s_list]
        R_list = [np.ascontiguousarray(r[np.ix_(keep, keep)]) for r in R_list]
        varY_list = [estimate_var_y(x, None if e is None else e[keep], n) for x, e, n in zip(s_list, eaf_list, n_list)]
        final_ids = [f for f, k in zip(final_ids, keep) if k]
        outlier = outlier[keep]

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
    cs_sets_df = build_locus_credible_sets_df(fit, labels, R_list, ancestries, snp_keys=snp_keys, outlier=outlier)
    max_s = max((q[0] for q in qc if q[0] == q[0]), default=float("nan"))
    mhc = bool(locus.get("mhc", False))
    region["LD_FLAG"] = ld_mismatch_flag(max_s, region["N_LD_OUTLIERS"], 0, mhc)
    if len(cs_sets_df):
        cs_sets_df.insert(0, "LOCUS", locus.get("name", f"{chrom}_{locus['center']}"))
        cs_sets_df["MAX_LD_MISMATCH_S"] = max_s
        cs_sets_df["LD_FLAG"] = [ld_mismatch_flag(max_s, region["N_LD_OUTLIERS"], k, mhc)
                                 for k in cs_sets_df["N_LD_OUTLIERS_IN_CS"]]
    cs_sets_df.attrs["region"] = region
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
        remove_ld_outliers=ctx.get("remove_ld_outliers", False),
    )
    region = {"LOCUS": locus.get("name", f"{locus['chr']}_{locus['center']}"), "CHR": locus["chr"],
              "START": locus["start"], "END": locus["end"], "LEADS": ",".join(map(str, locus.get("leads", []))),
              "MHC": bool(locus.get("mhc", False))}

    if result is None:
        region["STATUS"] = "skipped: too few shared variants, no LD, or MultiSuSiE failed"
        return i, n_loci, locus_label, False, region

    cs_sets_df, cs_variants_df = result
    region.update(cs_sets_df.attrs.get("region", {}))
    region["N_CS"] = len(cs_sets_df)
    region["STATUS"] = "fine-mapped"
    if len(cs_sets_df) == 0:
        return i, n_loci, locus_label, False, region
    if locus.get("cluster") is not None and "POS" in cs_variants_df.columns:
        # overlapping tiles: each variant's result is taken from the tile whose centre is nearest
        cs_variants_df["TILE_CORE"] = tile_core_mask(cs_variants_df["POS"].astype(float), ctx["loci"], i - 1)
    safe_name = locus.get("name", f"{locus['chr']}_{locus['center']}")
    output_dir = ctx["output_dir"]
    cs_file = os.path.join(output_dir, f"{safe_name}_locus_credible_sets.tsv")
    cs_variants_file = os.path.join(
        output_dir, f"{safe_name}_locus_credible_set_variants.tsv"
    )
    cs_sets_df.to_csv(cs_file, sep="\t", index=False)
    cs_variants_df.to_csv(cs_variants_file, sep="\t", index=False)
    return i, n_loci, locus_label, True, region


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
    ld_bfiles = [None if b in (None, "", "none", "NA") else b for b in (args.ld_bfiles or [None] * n_ancestries)]
    if len(ld_bfiles) != n_ancestries:
        print("ERROR: give one --ld_bfiles entry per ancestry ('none' where there is no user LD).", file=sys.stderr)
        sys.exit(1)

    print(f"Loading {n_ancestries} GWAS datasets...")
    gwas_dfs = [read_gwas(f) for f in args.gwas_files]
    clumped_dfs = [read_clumped(f) for f in args.clumped_files]

    print("Merging lead SNPs across ancestries into loci...")
    loci = get_union_loci(clumped_dfs, window_bp, args.max_region_kb * 1000, args.region_strategy)
    print(f"  {len(loci)} loci to fine-map ({args.region_strategy} regions)")
    names = [f"{locus['chr']}_{locus['center']}" for locus in loci]
    for locus, name in zip(loci, names):
        locus["name"] = name if names.count(name) == 1 else f"{name}_{locus['start']}"
        locus["mhc"] = overlaps_mhc(locus["chr"], locus["start"], locus["end"], args.build)

    ld_options = {"ld_reference": args.ld_reference, "thousand_genomes_dir": args.thousand_genomes_dir,
                  "ld_pack_dir": args.ld_pack_dir, "build": args.build, "ld_bfiles": ld_bfiles}
    gwas_sizes = [locus_sample_size(df, n) for df, n in zip(gwas_dfs, args.sample_sizes)]
    reference_sizes = [ld_reference_n_samples(a, ld_options, b) for a, b in zip(args.ancestries, ld_bfiles)]
    max_causal, l_rule = choose_max_causal(args.max_causal, gwas_sizes, reference_sizes, args.max_n_ratio)
    print(f"MultiSuSiE L = {max_causal} ({l_rule})")
    if max_causal == 1 and str(args.max_causal).lower() == "auto":
        print("WARNING: the GWAS are much larger than their LD references, so only one signal per locus is fitted. "
              "Use in-sample LD (ld_bfile), or set finemap.max_causal to fit more.", file=sys.stderr)
    settings = pd.DataFrame({
        "ANCESTRY": args.ancestries, "N_GWAS": gwas_sizes, "N_REF": reference_sizes,
        "LD_REFERENCE": [f"user:{os.path.basename(b)}" if b else args.ld_reference for b in ld_bfiles],
        "L": max_causal, "L_RULE": l_rule, "REGION_STRATEGY": args.region_strategy,
        "WINDOW_KB": args.window_kb, "MAX_REGION_KB": args.max_region_kb, "N_LOCI": len(loci)})
    settings.to_csv(os.path.join(args.output_dir, "multisusie_settings.tsv"), sep="\t", index=False)

    if len(loci) == 0:
        print("No loci to fine-map.")
        pd.DataFrame(columns=["LOCUS", "CHR", "START", "END", "STATUS"]).to_csv(
            os.path.join(args.output_dir, "multisusie_regions.tsv"), sep="\t", index=False)
        _write_completion(args.completion_file, 0)
        return

    worker_ctx = {
        "gwas_dfs": gwas_dfs,
        "ancestries": args.ancestries,
        "sample_sizes": args.sample_sizes,
        "ld_options": ld_options,
        "max_causal": max_causal,
        "coverage": args.coverage,
        "min_abs_corr": args.min_abs_corr,
        "output_dir": args.output_dir,
        "loci": loci,
        "remove_ld_outliers": args.remove_ld_outliers,
    }
    jobs = [(i + 1, len(loci), locus) for i, locus in enumerate(loci)]
    max_workers = min(
        calculate_parallelism(memory_per_worker_mb=6000),
        len(loci),
    )
    print(f"Running up to {max_workers} loci in parallel...")

    n_success = 0
    n_fail = 0
    region_rows = []

    with ProcessPoolExecutor(
        max_workers=max_workers,
        initializer=_init_worker,
        initargs=(worker_ctx,),
    ) as executor:
        futures = [executor.submit(_run_locus_job, job) for job in jobs]
        for future in as_completed(futures):
            i, n_loci, locus_label, ok, region = future.result()
            region_rows.append(region)
            if ok:
                n_success += 1
                print(f"  [{i}/{n_loci}] Completed {locus_label}")
            else:
                n_fail += 1
                print(f"  [{i}/{n_loci}] Skipped {locus_label} "
                      f"(insufficient shared data or MultiSuSiE failure)")

    regions_df = pd.DataFrame(region_rows)
    if len(regions_df):
        regions_df = regions_df.sort_values(
            ["CHR", "START"], key=lambda c: c.map(_chrom_order) if c.name == "CHR" else c)
    regions_df.to_csv(os.path.join(args.output_dir, "multisusie_regions.tsv"), sep="\t", index=False)
    print(f"\nMultiSuSiE complete: {n_success} loci succeeded, {n_fail} skipped/failed.")
    _write_completion(args.completion_file, n_success)


def _write_completion(path, n_loci):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(str(n_loci) + "\n")


if __name__ == "__main__":
    main()
