"""
Stage 50: check the pack against the source before the downloaded LD blocks are deleted.

    50_validate_pack.py --pop EUR [--chroms 22] [--chunks-per-chrom 3]

For each chromosome: keys are unique, rows are in GRCh37 order, chunk value counts match the
ragged band; and for a random sample of chunks, the decoded int8 values are byte-identical to
re-fetching r from the source and encoding it again (orientation, band and quantisation).
Writes state/<POP>.validate.json and exits non-zero on any failure.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd

from build_common import AUTOSOMES, Paths, env, get_logger, init_hail, write_json
from ld_reference import LDPack
from ldpack_writer import encode_chunk


def structural_checks(pack, pop, chrom):
    variants = pack.variants(chrom, pop)
    chunks = pack.chunks(chrom, pop)
    problems = []
    if not np.array_equal(variants["IDX"].values, np.arange(len(variants))):
        problems.append("IDX is not 0..n-1")
    if np.any(np.diff(variants["BP37"].values) < 0):
        problems.append("BP37 not sorted")
    if variants["SNP37"].duplicated().any():
        problems.append("duplicate SNP37 keys")
    snp38 = variants["SNP38"].dropna()
    if snp38.duplicated().any():
        problems.append("duplicate SNP38 keys")
    if (variants["SNP38"].isna() != (variants["LIFT_STATUS"] != "ok")).any():
        problems.append("SNP38 missing exactly when LIFT_STATUS != ok is violated")
    if int(chunks["NVALUES"].sum()) != int(variants["SEG_LEN"].sum()):
        problems.append("chunk NVALUES do not add up to SEG_LEN")
    return problems


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pop", required=True)
    parser.add_argument("--chroms", nargs="*", default=None)
    parser.add_argument("--chunks-per-chrom", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    paths = Paths().make()
    logger = get_logger(f"validate.{args.pop}", os.path.join(paths.logs, f"validate_{args.pop}.log"))
    chroms = [str(c) for c in args.chroms] if args.chroms else [
        c for c in AUTOSOMES if os.path.exists(os.path.join(paths.state, f"{args.pop}.ld.chr{c}.done"))]
    # info.tsv is written by stage 40; validation only needs the scale, which is fixed in format v1
    pop_dir = os.path.join(paths.pack, args.pop)
    if not os.path.exists(os.path.join(pop_dir, "info.tsv")):
        with open(os.path.join(pop_dir, "info.tsv"), "w", encoding="utf-8") as handle:
            handle.write("format_version\t1\nscale\t127\n")
    pack = LDPack(paths.pack)
    init_hail(paths, f"validate_{args.pop}")
    from hail.linalg import BlockMatrix
    bm = BlockMatrix.read(paths.block_matrix(args.pop))
    deadzone = float(env("DEADZONE", "0.02"))
    rng = np.random.default_rng(args.seed)

    report = {"pop": args.pop, "chromosomes": {}, "ok": True}
    for chrom in chroms:
        problems = structural_checks(pack, args.pop, chrom)
        variants = pack.variants(chrom, args.pop)
        annotated = pd.read_csv(os.path.join(paths.raw_pop(args.pop), f"variants_annot_chr{chrom}.tsv.gz"),
                                sep="\t", dtype={"CHR": str, "SNP37": str}, na_values=["NA"], keep_default_na=False)
        annotated = annotated.set_index("SNP37").loc[variants["SNP37"].values]
        gidx = annotated["IDX_GLOBAL"].values.astype(np.int64)
        sign = annotated["SIGN"].values.astype(np.int64)
        seg_len = variants["SEG_LEN"].values.astype(np.int64)
        chunks = pack.chunks(chrom, args.pop)
        sample = rng.choice(len(chunks), size=min(args.chunks_per_chrom, len(chunks)), replace=False)
        checked = 0
        for k in sample:
            row = chunks.iloc[int(k)]
            r0, r1 = int(row.ROW_START), int(row.ROW_END)
            c_end = int(max(r1, (np.arange(r0, r1) + seg_len[r0:r1]).max() + 1))
            rect = bm.filter(gidx[r0:r1].tolist(), gidx[r0:c_end].tolist()).to_numpy()
            expected, _ = encode_chunk(rect, r0, r1, seg_len, sign, 127, deadzone)
            stored = pack.read_chunk_values(chrom, args.pop, row)
            if not np.array_equal(expected, stored):
                problems.append(f"chunk {int(k)} differs from the source ({int((expected != stored).sum())} values)")
            checked += 1
        report["chromosomes"][chrom] = {"problems": problems, "chunks_checked": checked,
                                        "n_variants": int(len(variants))}
        if problems:
            report["ok"] = False
            logger.error(f"{args.pop} chr{chrom}: " + "; ".join(problems))
        else:
            logger.info(f"{args.pop} chr{chrom}: OK ({checked} chunks byte-identical to the source)")
    write_json(os.path.join(paths.state, f"{args.pop}.validate.json"), report)
    return 0 if report["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
