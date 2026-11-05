"""
Put the Pan-UKB LD scores into the layout GeneHackman's LDSC step reads (scripts/run_ldsc.sh):

    python 20_import_pan_ukb_ldscores.py            (settings from config.sh via the environment)

Downloads UKBB.ALL.ldscore.tar.gz (LDSC-format LD scores of HapMap3 variants for EUR, AFR, AMR, CSA, EAS, MID;
https://pan.ukbb.broadinstitute.org/docs/ld), and writes for each population
    $BUNDLE_DIR/LDSCORE/pan_ukb/<POP>/ld.l2.ldscore.gz   CHR SNP BP L2, SNP = rsID
    $BUNDLE_DIR/LDSCORE/pan_ukb/<POP>/ld.l2.M, ld.l2.M_5_50
and $BUNDLE_DIR/LDSCORE/pan_ukb/w_hm3.snplist (copied from the GRCh37 LD score set; it is keyed by rsID).

LDSC matches variants on rsID, so the scores work with GRCh37 and GRCh38 GWAS alike.  The archive's layout and
SNP column aren't documented, so the script picks each population's single-annotation scores (UKBB.<POP>.rsid.l2...,
else UKBB.<POP>.l2...; the stratified 8LDMS/25LDMS sets are skipped) and, if the SNP column is not rsIDs, maps
CHR:BP:REF:ALT (GRCh37) to rsIDs with the Pan-UKB LD pack's RSID_UKB column.  It prints what
it found; check that output (see RUNBOOK.md).
"""

import glob
import os
import shutil
import subprocess
import sys
import tarfile

import pandas as pd

POPULATIONS = ["EUR", "AFR", "AMR", "CSA", "EAS", "MID"]


def env(name, default=None):
    value = os.environ.get(name, "")
    return value if value else default


def find_population_files(extract_dir, pop):
    """The single-annotation LD scores of a population: (ldscore file, prefix of its .l2.M files).

    The archive (UKBB.ALL.ldscore/) has, per population, UKBB.<POP>.l2.ldscore.gz (SNP = CHR:BP:REF:ALT),
    UKBB.<POP>.rsid.l2.ldscore.gz (the same scores with SNP = rsID) and UKBB.<POP>.l2.M / .l2.M_5_50, plus MAF/LD
    stratified sets (<POP>.8LDMS, <POP>.25LDMS) for S-LDSC, which are not used.  The rsID file is preferred."""
    def matches(pattern):
        files = glob.glob(os.path.join(extract_dir, "**", pattern), recursive=True)
        return [f for f in files if "LDMS" not in os.path.basename(f)]
    for pattern in (f"*.{pop}.rsid.l2.ldscore.gz", f"*.{pop}.l2.ldscore.gz", f"*{pop}*.l2.ldscore.gz"):
        files = matches(pattern)
        if len(files) == 1:
            ldscore = files[0]
            return ldscore, ldscore[: -len(".l2.ldscore.gz")].replace(".rsid", "")
        if len(files) > 1:
            raise SystemExit(f"{pop}: expected one {pattern} in the archive, found {files}")
    raise SystemExit(f"{pop}: no *{pop}*.l2.ldscore.gz in the archive")


def rsid_lookup(pack_dir, pop):
    """GRCh37 CHR:BP:REF:ALT -> rsID from the Pan-UKB LD pack's variant tables."""
    frames = []
    for path in sorted(glob.glob(os.path.join(pack_dir, pop, "variants_chr*.tsv.gz"))):
        frames.append(pd.read_csv(path, sep="\t", dtype=str, usecols=["CHR", "BP37", "REF37", "ALT37", "RSID_UKB"]))
    if not frames:
        raise SystemExit(f"No LD pack variant tables for {pop} under {pack_dir}: needed to map variant IDs to rsIDs")
    table = pd.concat(frames)
    keys = table["CHR"] + ":" + table["BP37"] + ":" + table["REF37"] + ":" + table["ALT37"]
    return dict(zip(keys, table["RSID_UKB"]))


def convert_population(ldscore_file, prefix, pop, out_dir, pack_dir):
    scores = pd.read_csv(ldscore_file, sep=r"\s+")
    print(f"{pop}: {os.path.basename(ldscore_file)}, {len(scores):,} rows, columns {list(scores.columns)}")
    print(scores.head(3).to_string(index=False))
    snp_column = "SNP" if "SNP" in scores.columns else next(c for c in scores.columns if c.lower() in ("rsid", "snp"))
    l2_columns = [c for c in scores.columns if c.endswith("L2")]
    if len(l2_columns) != 1:
        raise SystemExit(f"{pop}: expected one L2 column, found {l2_columns}")
    ids = scores[snp_column].astype(str)
    share_rs = ids.str.match(r"^rs\d+$").mean()
    if share_rs > 0.95:
        rsids = ids
    else:
        lookup = rsid_lookup(pack_dir, pop)
        normalised = ids.str.replace(r"^chr", "", regex=True).str.replace(r"[_/]", ":", regex=True)
        rsids = normalised.map(lookup)
        print(f"{pop}: SNP column is not rsIDs ({100 * share_rs:.1f}% are); mapped "
              f"{rsids.notna().sum():,} of {len(ids):,} variant IDs to rsIDs with the LD pack")
    out = pd.DataFrame({"CHR": scores["CHR"], "SNP": rsids, "BP": scores["BP"], "L2": scores[l2_columns[0]]})
    out = out[out["SNP"].notna() & out["SNP"].str.startswith("rs")]
    out = out[~out["SNP"].duplicated(keep=False)]
    pop_dir = os.path.join(out_dir, pop)
    os.makedirs(pop_dir, exist_ok=True)
    out.to_csv(os.path.join(pop_dir, "ld.l2.ldscore.gz"), sep="\t", index=False, compression="gzip")
    for suffix in (".l2.M", ".l2.M_5_50"):
        if not os.path.exists(prefix + suffix):
            raise SystemExit(f"{pop}: {prefix + suffix} not found in the archive")
        shutil.copy(prefix + suffix, os.path.join(pop_dir, "ld" + suffix))
    print(f"{pop}: wrote {len(out):,} LD scores to {pop_dir}")
    return out


def main():
    work_dir = env("WORK_DIR")
    bundle_dir = env("BUNDLE_DIR")
    old_data = env("OLD_PIPELINE_DATA_DIR")
    pack_dir = env("LD_PACK_DIR")
    url = env("PAN_UKB_LDSCORE_URL")
    if not (work_dir and bundle_dir and url):
        raise SystemExit("source config.sh first (WORK_DIR, BUNDLE_DIR, PAN_UKB_LDSCORE_URL)")
    download_dir = os.path.join(work_dir, "pan_ukb_ldscores")
    os.makedirs(download_dir, exist_ok=True)
    archive = os.path.join(download_dir, os.path.basename(url))
    if not os.path.exists(archive + ".done"):
        subprocess.run(["curl", "-fSL", "-C", "-", url, "-o", archive], check=True)
        open(archive + ".done", "w").close()
    extract_dir = os.path.join(download_dir, "extracted")
    if not os.path.isdir(extract_dir):
        with tarfile.open(archive) as tar:
            members = [m for m in tar.getmembers() if not os.path.isabs(m.name) and ".." not in m.name.split("/")]
            tar.extractall(extract_dir, members=members)
    print("Archive contents:", sorted(os.path.relpath(p, extract_dir)
                                      for p in glob.glob(os.path.join(extract_dir, "**", "*"), recursive=True))[:40])

    out_dir = os.path.join(bundle_dir, "LDSCORE", "pan_ukb")
    os.makedirs(out_dir, exist_ok=True)
    hm3 = set()
    snplist_source = os.path.join(old_data or "", "LDSCORE", "b37_dbsnp156", "w_hm3.snplist")
    if old_data and os.path.exists(snplist_source):
        shutil.copy(snplist_source, os.path.join(out_dir, "w_hm3.snplist"))
        hm3 = set(pd.read_csv(snplist_source, sep=r"\s+", usecols=[0]).iloc[:, 0])
    else:
        print(f"WARNING: {snplist_source} not found: copy w_hm3.snplist into {out_dir} by hand", file=sys.stderr)
    for pop in POPULATIONS:
        out = convert_population(*find_population_files(extract_dir, pop), pop, out_dir, pack_dir)
        if hm3:
            print(f"{pop}: {out['SNP'].isin(hm3).mean():.1%} of the LD score SNPs are in w_hm3.snplist")
    return 0


if __name__ == "__main__":
    sys.exit(main())
