"""
Stage 10: list the variants to keep for one population.

Reads the downloaded Pan-UKB variant index (locus, alleles, idx) and LD score table (AF, rsid),
keeps biallelic autosomal variants with MAF >= MAF_MIN in that population, and writes
raw/<POP>/variants_raw.tsv.gz  (CHR, BP37, REF, ALT, AF, IDX_GLOBAL, RSID_UKB)
raw/<POP>/summary.json         (n_samples, counts by chromosome and MAF bin, AF coverage, schemas)

AF is the frequency of alleles[1] (ALT) as stored in the LD score table; phase 0 checks this
against 1000 Genomes before the full build relies on it.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

from build_common import AUTOSOMES, Paths, env, get_logger, init_hail, write_json

MAF_BINS = [0.0, 0.001, 0.005, 0.01, 0.05, 0.5000001]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pop", required=True)
    parser.add_argument("--maf-min", type=float, default=float(env("MAF_MIN", "0.01")))
    args = parser.parse_args()

    paths = Paths().make()
    logger = get_logger(f"variants.{args.pop}", os.path.join(paths.logs, f"variants_{args.pop}.log"))
    hl = init_hail(paths, f"variants_{args.pop}")

    index = hl.read_table(paths.variant_ht(args.pop))
    ldscore = hl.read_table(paths.ldscore_ht(args.pop))
    schemas = {"variant_index": str(index.row.dtype), "variant_index_globals": str(index.globals.dtype),
               "ldscore": str(ldscore.row.dtype)}
    logger.info(f"variant index schema: {schemas['variant_index']}; globals {schemas['variant_index_globals']}")
    logger.info(f"LD score table schema: {schemas['ldscore']}")
    for field in ("AF",):
        if field not in ldscore.row:
            raise SystemExit(f"LD score table has no '{field}' field; schema is {schemas['ldscore']}")
    n_samples = hl.eval(index.globals.n_samples) if "n_samples" in index.globals else None

    joined = index.annotate(AF=ldscore[index.key].AF,
                            RSID_UKB=(ldscore[index.key].rsid if "rsid" in ldscore.row else hl.missing(hl.tstr)))
    joined = joined.select(CHR=joined.locus.contig, BP37=joined.locus.position,
                           REF=joined.alleles[0],
                           ALT=hl.if_else(hl.len(joined.alleles) == 2, joined.alleles[1], hl.missing(hl.tstr)),
                           AF=joined.AF, IDX_GLOBAL=joined.idx, RSID_UKB=joined.RSID_UKB)
    all_path = os.path.join(paths.raw_pop(args.pop), "variants_all.tsv.bgz")
    joined.key_by().select("CHR", "BP37", "REF", "ALT", "AF", "IDX_GLOBAL", "RSID_UKB").export(all_path)

    logger.info("summarising and filtering")
    frame = pd.read_csv(all_path, sep="\t", compression="gzip", dtype={"CHR": str, "RSID_UKB": str},
                        na_values=["NA"], keep_default_na=False)
    n_total = len(frame)
    has_af = frame["AF"].notna()
    biallelic = frame["ALT"].notna()
    autosomal = frame["CHR"].isin(AUTOSOMES)
    frame["MAF"] = np.minimum(frame["AF"], 1 - frame["AF"])
    keep = has_af & biallelic & autosomal & (frame["MAF"] >= args.maf_min)

    maf_bins = pd.cut(frame.loc[has_af & autosomal, "MAF"], MAF_BINS, right=False)
    summary = {
        "pop": args.pop, "n_samples": n_samples, "maf_min": args.maf_min, "schemas": schemas,
        "n_index_variants": int(n_total), "n_with_af": int(has_af.sum()),
        "af_coverage": float(has_af.mean()) if n_total else 0.0,
        "n_not_biallelic": int((~biallelic).sum()), "n_non_autosomal": int((~autosomal).sum()),
        "n_kept": int(keep.sum()),
        "kept_by_chrom": {c: int(v) for c, v in frame.loc[keep, "CHR"].value_counts().items()},
        "maf_bins": {str(k): int(v) for k, v in maf_bins.value_counts(sort=False).items()},
    }
    if summary["af_coverage"] < 0.99:
        logger.warning(f"only {summary['af_coverage']:.1%} of index variants have AF in the LD score table; "
                       "variants without AF are dropped")

    kept = frame.loc[keep, ["CHR", "BP37", "REF", "ALT", "AF", "IDX_GLOBAL", "RSID_UKB"]].copy()
    kept["CHR_NUM"] = kept["CHR"].astype(int)
    kept = kept.sort_values(["CHR_NUM", "BP37", "IDX_GLOBAL"], kind="stable").drop(columns="CHR_NUM")
    monotonic = bool(all(np.all(np.diff(group["IDX_GLOBAL"].values) > 0)
                         for _, group in kept.groupby("CHR", sort=False)))
    summary["idx_monotonic_in_position"] = monotonic
    if not monotonic:
        raise SystemExit("Source matrix index is not increasing with genomic position; the pack layout "
                         "assumes it is. Inspect raw/<POP>/variants_all.tsv.bgz.")
    out_path = os.path.join(paths.raw_pop(args.pop), "variants_raw.tsv.gz")
    kept.to_csv(out_path + ".tmp", sep="\t", index=False, compression="gzip")
    os.replace(out_path + ".tmp", out_path)
    write_json(os.path.join(paths.raw_pop(args.pop), "summary.json"), summary)
    os.remove(all_path)
    logger.info(f"{args.pop}: kept {summary['n_kept']:,} of {n_total:,} variants (MAF >= {args.maf_min}); "
                f"n_samples = {n_samples}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
