"""
Stage 45: refresh the GRCh38 keys of a finished pack without touching its LD.

    45_update_pack_keys.py --pop EUR

After 20_annotate_variants.R has been re-run (e.g. with a newer GeneHackman liftover), this copies
BP38, SNP38 and LIFT_STATUS from raw/<POP>/variants_annot_chr<C>.tsv.gz into the pack's
variants_chr<C>.tsv.gz.  The pack rows are the annotated variants in (BP37, source index) order, as
the writer and 50_validate_pack.py assume; the script checks that BP37/REF37/ALT37 line up row by
row before replacing anything, and leaves every other column (and ld_chr<C>.bin) as it is.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

from build_common import AUTOSOMES, Paths, get_logger

UPDATED_COLUMNS = ["BP38", "SNP38", "LIFT_STATUS"]


def update_chromosome(paths, pop, chrom, logger):
    pack_file = os.path.join(paths.pack, pop, f"variants_chr{chrom}.tsv.gz")
    annot_file = os.path.join(paths.raw_pop(pop), f"variants_annot_chr{chrom}.tsv.gz")
    if not os.path.exists(pack_file):
        return None
    if not os.path.exists(annot_file):
        raise SystemExit(f"{annot_file} not found: re-run 20_annotate_variants.R --pop {pop} first")
    pack = pd.read_csv(pack_file, sep="\t", dtype=str, keep_default_na=False, na_values=[""])
    text_columns = {"CHR": str, "SNP37": str, "SNP38": str, "RSID_UKB": str, "REF37": str, "ALT37": str}
    annot = pd.read_csv(annot_file, sep="\t", dtype=text_columns, na_values=["NA"], keep_default_na=False)
    annot = annot.sort_values(["BP37", "IDX_GLOBAL"], kind="stable").reset_index(drop=True)
    if len(annot) != len(pack):
        raise SystemExit(f"{pop} chr{chrom}: {len(pack)} pack rows but {len(annot)} annotated variants")
    for column in ("BP37", "REF37", "ALT37"):
        left = pack[column].astype(str).values
        right = annot[column].astype(np.int64 if column == "BP37" else str).astype(str).values
        if not np.array_equal(left, right):
            raise SystemExit(f"{pop} chr{chrom}: pack rows and annotated variants differ in {column}")
    before = pack["LIFT_STATUS"].value_counts().to_dict()
    pack["BP38"] = ["" if pd.isna(bp) else str(int(round(bp))) for bp in annot["BP38"]]
    pack["SNP38"] = annot["SNP38"].fillna("").values
    pack["LIFT_STATUS"] = annot["LIFT_STATUS"].values
    tmp = pack_file + ".tmp.gz"
    pack.to_csv(tmp, sep="\t", index=False, compression="gzip")
    os.replace(tmp, pack_file)
    after = pack["LIFT_STATUS"].value_counts().to_dict()
    logger.info(f"{pop} chr{chrom}: LIFT_STATUS {before} -> {after}")
    return after


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pop", required=True)
    args = parser.parse_args()
    paths = Paths().make()
    logger = get_logger(f"update_keys.{args.pop}", os.path.join(paths.logs, f"update_keys_{args.pop}.log"))
    updated = [c for c in AUTOSOMES if update_chromosome(paths, args.pop, c, logger) is not None]
    logger.info(f"{args.pop}: GRCh38 keys refreshed for {len(updated)} chromosomes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
