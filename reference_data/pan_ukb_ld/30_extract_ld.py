"""
Stage 30: read r from the downloaded Pan-UKB LD blocks and write the pack, chromosome by chromosome.

    30_extract_ld.py band  --pop EUR [--chroms 22] [--measure]
    30_extract_ld.py dense --pop EUR --region 6:35600000-37600000 --out phase0/EUR_region.npz

band:  for each chromosome, rows are the kept variants (stage 20 output) in GRCh37 order; r is
       fetched from the source in windows of FETCH_ROWS rows, cut into the ragged +/-RADIUS_KB
       band, oriented to the canonical EA allele, quantised and xz-compressed (see ldpack_writer).
       Every chunk is checkpointed, so a killed run resumes where it stopped; a finished chromosome
       leaves state/<POP>.ld.chr<C>.done.  --measure (phase 0) also records the compressed size of
       each chunk under alternative dead-zones and gzip, for the size report.
dense: full float r (EA-oriented) for all kept variants in a region, for the phase 0 check of
       how much LD lies beyond the stored radius.
"""

import argparse
import gzip
import json
import lzma
import os
import sys
import time

import numpy as np
import pandas as pd

from build_common import AUTOSOMES, Paths, env, get_logger, init_hail, write_json
from ldpack_writer import ChromosomeWriter, encode_chunk

MEASURE_DEADZONES = (0.0, 0.02, 0.03)


def read_annotated(paths, pop, chrom):
    path = os.path.join(paths.raw_pop(pop), f"variants_annot_chr{chrom}.tsv.gz")
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found: run 20_annotate_variants.R --pop {pop} first")
    return pd.read_csv(path, sep="\t", dtype={"CHR": str, "SNP37": str, "SNP38": str, "RSID_UKB": str},
                       na_values=["NA"], keep_default_na=False)


def windows(pending, bounds, fetch_rows):
    """Group consecutive pending chunks into fetch windows of at most fetch_rows rows."""
    groups, current = [], []
    for k in pending:
        if current and (k != current[-1] + 1 or bounds[k][1] - bounds[current[0]][0] > fetch_rows):
            groups.append(current)
            current = []
        current.append(k)
    if current:
        groups.append(current)
    return groups


def measure_chunk(rect, r0, r1, writer):
    sizes = {}
    for deadzone in MEASURE_DEADZONES:
        values, _ = encode_chunk(rect, r0, r1, writer.seg_len, writer.sign, writer.scale, deadzone)
        raw = values.tobytes()
        sizes[f"xz_dz{deadzone}"] = len(lzma.compress(raw, preset=6))
        sizes[f"gzip_dz{deadzone}"] = len(gzip.compress(raw, compresslevel=6))
        sizes["n_values"] = int(values.size)
        sizes[f"zero_fraction_dz{deadzone}"] = float((values == 0).mean()) if values.size else 0.0
    return sizes


def run_band(args, paths, logger, hl):
    from hail.linalg import BlockMatrix
    pop = args.pop
    radius_bp = int(float(env("RADIUS_KB", "500")) * 1000)
    chunk_rows = int(env("CHUNK_ROWS", "2048"))
    fetch_rows = int(env("FETCH_ROWS", "8192"))
    deadzone = float(env("DEADZONE", "0.02"))
    bm = BlockMatrix.read(paths.block_matrix(pop))
    pop_pack = os.path.join(paths.pack, pop)
    pop_work = os.path.join(paths.chunks, pop)
    chunks_written = 0

    for chrom in args.chroms:
        done_marker = os.path.join(paths.state, f"{pop}.ld.chr{chrom}.done")
        if os.path.exists(done_marker):
            logger.info(f"{pop} chr{chrom}: already done")
            continue
        variants = read_annotated(paths, pop, chrom)
        writer = ChromosomeWriter(pop_pack, pop_work, chrom, variants, radius_bp, chunk_rows=chunk_rows,
                                  deadzone=deadzone)
        gidx = writer.variants["IDX_GLOBAL"].values.astype(np.int64)
        pending = writer.pending_chunks()
        n_chunks = len(writer.bounds)
        logger.info(f"{pop} chr{chrom}: {len(writer.variants):,} variants, {n_chunks} chunks, "
                    f"{len(pending)} to do, {int(writer.seg_len.sum()):,} stored pairs")
        measure_file = os.path.join(paths.phase0, f"measure_{pop}_chr{chrom}.jsonl") if args.measure else None
        started = time.time()
        done_rows = 0
        todo_rows = sum(writer.bounds[k][1] - writer.bounds[k][0] for k in pending)
        for group in windows(pending, writer.bounds, fetch_rows):
            r0, r1 = writer.bounds[group[0]][0], writer.bounds[group[-1]][1]
            c_end = max(writer.column_end(k) for k in group)
            fetch_start = time.time()
            rect = bm.filter(gidx[r0:r1].tolist(), gidx[r0:c_end].tolist()).to_numpy()
            fetch_seconds = time.time() - fetch_start
            for k in group:
                k0, k1 = writer.bounds[k]
                sub = rect[k0 - r0:k1 - r0, k0 - r0:writer.column_end(k) - r0]
                meta = writer.write_chunk(k, sub)
                chunks_written += 1
                if args.stop_after_chunks is not None and chunks_written >= args.stop_after_chunks:
                    logger.warning(f"stopping after {chunks_written} chunks (--stop-after-chunks)")
                    sys.exit(3)
                if measure_file:
                    record = {"pop": pop, "chrom": chrom, "chunk": k, "fetch_seconds": fetch_seconds / len(group),
                              "nbytes": meta["nbytes"], "n_nan": meta["n_nan"]}
                    record.update(measure_chunk(sub, k0, k1, writer))
                    with open(measure_file, "a", encoding="utf-8") as handle:
                        handle.write(json.dumps(record) + "\n")
                if meta["n_nan"]:
                    logger.warning(f"{pop} chr{chrom} chunk {k}: {meta['n_nan']} NaN values stored as 0")
            done_rows += r1 - r0
            elapsed = time.time() - started
            eta_hours = elapsed / done_rows * (todo_rows - done_rows) / 3600 if done_rows else float("nan")
            logger.info(f"{pop} chr{chrom}: {done_rows:,}/{todo_rows:,} rows "
                        f"({fetch_seconds:.0f}s for this window; ETA {eta_hours:.1f} h for this chromosome)")
        totals = writer.finalise()
        totals["seconds"] = time.time() - started
        write_json(os.path.join(paths.state, f"{pop}.ld.chr{chrom}.json"), totals)
        open(done_marker, "w").close()
        logger.info(f"{pop} chr{chrom}: finished, {totals['bytes'] / 1e6:.1f} MB")


def run_dense(args, paths, logger, hl):
    from hail.linalg import BlockMatrix
    chrom, span = args.region.split(":")
    start, end = (int(x) for x in span.split("-"))
    variants = read_annotated(paths, args.pop, chrom)
    region = variants[(variants["BP37"] >= start) & (variants["BP37"] <= end)]
    region = region.sort_values(["BP37", "IDX_GLOBAL"], kind="stable")
    gidx = region["IDX_GLOBAL"].values.astype(np.int64).tolist()
    logger.info(f"{args.pop} {args.region}: dense extraction of {len(gidx)} variants")
    bm = BlockMatrix.read(paths.block_matrix(args.pop))
    upper = bm.filter(gidx, gidx).to_numpy()
    upper = np.triu(upper, 1)
    r_alt = upper + upper.T + np.eye(len(gidx))
    sign = region["SIGN"].values.astype(np.float64)
    r_ea = r_alt * np.outer(sign, sign)
    out = args.out if os.path.isabs(args.out) else os.path.join(paths.work, args.out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez_compressed(out, r=r_ea.astype(np.float32), snp37=region["SNP37"].values.astype(str),
                        snp38=region["SNP38"].fillna("").values.astype(str), bp37=region["BP37"].values,
                        eaf=region["EAF"].values)
    # the same matrix for R (62_truncation_finemap.R): raw little-endian float32 plus the variant list
    stem = out[:-4] if out.endswith(".npz") else out
    r_ea.astype("<f4").tofile(stem + ".f32")
    region[["SNP37", "SNP38", "BP37", "EAF"]].to_csv(stem + ".variants.tsv", sep="\t", index=False)
    logger.info(f"wrote {out} and {stem}.f32")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["band", "dense"])
    parser.add_argument("--pop", required=True)
    parser.add_argument("--chroms", nargs="*", default=None)
    parser.add_argument("--measure", action="store_true")
    parser.add_argument("--region", default=None)
    parser.add_argument("--out", default=None)
    parser.add_argument("--stop-after-chunks", type=int, default=None,
                        help="testing only: exit after writing this many chunks, as if the run were killed")
    args = parser.parse_args()

    paths = Paths().make()
    if args.chroms:
        args.chroms = [str(c) for c in args.chroms]
    else:  # every autosome with kept variants (stage 20 writes one file per chromosome)
        args.chroms = [c for c in AUTOSOMES
                       if os.path.exists(os.path.join(paths.raw_pop(args.pop), f"variants_annot_chr{c}.tsv.gz"))]
    logger = get_logger(f"ld.{args.pop}", os.path.join(paths.logs, f"ld_{args.pop}.log"))
    hl = init_hail(paths, f"ld_{args.pop}")
    if args.command == "band":
        run_band(args, paths, logger, hl)
    else:
        if not args.region or not args.out:
            raise SystemExit("dense needs --region and --out")
        run_dense(args, paths, logger, hl)
    return 0


if __name__ == "__main__":
    sys.exit(main())
