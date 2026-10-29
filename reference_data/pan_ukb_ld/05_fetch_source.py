"""
Download the parts of the Pan-UKB LD release needed for the pack, one population at a time.

    05_fetch_source.py tables --pop EUR                  variant index + LD score tables (small)
    05_fetch_source.py plan   --pop EUR                  bytes of LD blocks the full pack needs (no download)
    05_fetch_source.py blocks --pop EUR [--chroms 22]    download only the LD blocks the pack needs
    05_fetch_source.py blocks --pop EUR --region 6:35600000-37600000 --radius-kb 2000
    05_fetch_source.py clean-blocks --pop EUR            delete downloaded LD blocks (keeps tables)

The Pan-UKB matrices are Hail BlockMatrix directories (one file per stored block).  "blocks" works
out which blocks hold r for kept variants (MAF >= MAF_MIN) within RADIUS_KB of each other, downloads
only those files, and writes a local metadata.json listing exactly the blocks present, so Hail reads
the local copy as a sparse block matrix.  Re-running resumes: complete files are skipped.
"""

import argparse
import json
import math
import os
import shutil
import sys

import numpy as np
import pandas as pd

from build_common import (AUTOSOMES, Paths, download_files, env, free_bytes, get_logger, open_store,
                          read_json, source_names, write_json)
from ld_reference import compute_seg_len


def fetch_tables(store, paths, pop, logger, threads):
    names = source_names(pop)
    for table in ("variant_ht", "ldscore_ht"):
        prefix = names[table]
        items = store.list(prefix + "/")
        if not items:
            raise SystemExit(f"Nothing found at {prefix}/ in {env('SOURCE_URI')}. Check SOURCE_URI and the "
                             f"Pan-UKB release layout (see README.md).")
        download_files(store, items, os.path.join(paths.source_pop(pop), prefix), prefix, logger, threads)
    logger.info(f"{pop}: tables present under {paths.source_pop(pop)}")


def read_block_metadata(store, paths, pop):
    prefix = source_names(pop)["block_matrix"]
    cache = os.path.join(paths.source_pop(pop), f"{prefix}.source_metadata.json")
    meta = read_json(cache)
    if meta is None:
        os.makedirs(paths.source_pop(pop), exist_ok=True)
        meta = json.loads(store.read_bytes(f"{prefix}/metadata.json"))
        write_json(cache, meta)
    return meta


def list_part_sizes(store, paths, pop, logger):
    prefix = source_names(pop)["block_matrix"]
    cache = os.path.join(paths.source_pop(pop), f"{prefix}.part_sizes.json")
    sizes = read_json(cache)
    if sizes is None:
        logger.info(f"{pop}: listing LD block files (one-off, may take a few minutes)")
        sizes = {key.rsplit("/", 1)[-1]: size for key, size in store.list(f"{prefix}/parts/")}
        write_json(cache, sizes)
    return sizes


def kept_variants(paths, pop, chroms):
    path = os.path.join(paths.raw_pop(pop), "variants_raw.tsv.gz")
    if not os.path.exists(path):
        raise SystemExit(f"{path} not found: run 10_extract_variants.py --pop {pop} first")
    variants = pd.read_csv(path, sep="\t", dtype={"CHR": str})
    if chroms:
        variants = variants[variants["CHR"].isin(chroms)]
    return variants


def needed_blocks_for_band(variants, radius_bp, block_size):
    """(block_row, block_col) pairs holding r for kept variant pairs within radius_bp, upper triangle."""
    max_col = {}
    for _, chrom_variants in variants.groupby("CHR", sort=False):
        chrom_variants = chrom_variants.sort_values(["BP37", "IDX_GLOBAL"], kind="stable")
        gidx = chrom_variants["IDX_GLOBAL"].values.astype(np.int64)
        seg = compute_seg_len(chrom_variants["BP37"].values, radius_bp)
        gstop = gidx[np.arange(gidx.size) + seg]
        block_row = gidx // block_size
        block_col = gstop // block_size
        frame = pd.DataFrame({"br": block_row, "bc": block_col}).groupby("br")["bc"].max()
        for br, bc in frame.items():
            max_col[int(br)] = max(max_col.get(int(br), int(br)), int(bc))
    return {(br, bc) for br, top in max_col.items() for bc in range(br, top + 1)}


def needed_blocks_for_region(variants, chrom, start, end, block_size):
    region = variants[(variants["CHR"] == str(chrom)) & (variants["BP37"] >= start) & (variants["BP37"] <= end)]
    if region.empty:
        raise SystemExit(f"No kept variants in region {chrom}:{start}-{end}")
    lo = int(region["IDX_GLOBAL"].min()) // block_size
    hi = int(region["IDX_GLOBAL"].max()) // block_size
    return {(br, bc) for br in range(lo, hi + 1) for bc in range(br, hi + 1)}


def map_blocks_to_parts(blocks, meta):
    """Translate (block_row, block_col) to source part files, using Hail's column-major block index."""
    n_block_rows = math.ceil(meta["nRows"] / meta["blockSize"])
    wanted = sorted(br + bc * n_block_rows for br, bc in blocks)
    filtered = meta.get("maybeFiltered")
    parts = meta["partFiles"]
    if filtered is None:
        if len(parts) != n_block_rows * math.ceil(meta["nCols"] / meta["blockSize"]):
            raise SystemExit("Unexpected BlockMatrix layout: partFiles does not match the block grid")
        return {b: parts[b] for b in wanted if b < len(parts)}
    if len(filtered) != len(parts):
        raise SystemExit("Unexpected BlockMatrix layout: maybeFiltered and partFiles differ in length")
    position = {block: i for i, block in enumerate(filtered)}
    return {b: parts[position[b]] for b in wanted if b in position}


def write_local_metadata(paths, pop, meta, present_blocks):
    """metadata.json for the local copy, listing only the blocks downloaded so far."""
    local_dir = paths.block_matrix(pop)
    present = sorted(present_blocks)
    position = ({b: i for i, b in enumerate(meta["maybeFiltered"])} if meta.get("maybeFiltered") is not None
                else {b: b for b in present})
    local = dict(meta)
    local["maybeFiltered"] = present
    local["partFiles"] = [meta["partFiles"][position[b]] for b in present]
    write_json(os.path.join(local_dir, "metadata.json"), local)
    open(os.path.join(local_dir, "_SUCCESS"), "w").close()
    write_json(os.path.join(local_dir, "blocks_present.json"), present)


def fetch_blocks(store, paths, pop, logger, threads, chroms=None, region=None, radius_kb=None, plan_only=False):
    meta = read_block_metadata(store, paths, pop)
    block_size = int(meta["blockSize"])
    variants = kept_variants(paths, pop, chroms if not region else None)
    if region:
        chrom, span = region.split(":")
        start, end = (int(x) for x in span.split("-"))
        blocks = needed_blocks_for_region(variants, chrom, start, end, block_size)
    else:
        radius_bp = int(float(radius_kb or env("RADIUS_KB", "1000")) * 1000)
        blocks = needed_blocks_for_band(variants, radius_bp, block_size)
    parts = map_blocks_to_parts(blocks, meta)
    sizes = list_part_sizes(store, paths, pop, logger)
    prefix = source_names(pop)["block_matrix"]
    items = [(f"{prefix}/parts/{name}", sizes[name]) for name in parts.values() if name in sizes]
    missing = [name for name in parts.values() if name not in sizes]
    if missing:
        raise SystemExit(f"{len(missing)} block files listed in metadata are missing from the source")
    total = sum(size for _, size in items)
    logger.info(f"{pop}: {len(blocks)} blocks needed ({len(parts)} stored in the source), "
                f"{total / 1e12:.2f} TB")
    if plan_only:
        return {"pop": pop, "blocks": len(parts), "bytes": total}

    local_dir = paths.block_matrix(pop)
    already = {os.path.basename(k) for k, s in items
               if os.path.exists(os.path.join(local_dir, "parts", os.path.basename(k)))
               and os.path.getsize(os.path.join(local_dir, "parts", os.path.basename(k))) == s}
    remaining = sum(s for k, s in items if os.path.basename(k) not in already)
    margin = float(env("MIN_FREE_GB", "500")) * 1e9
    free = free_bytes(paths.source)
    if remaining + margin > free:
        raise SystemExit(f"Not enough disk: need {remaining / 1e12:.2f} TB + {margin / 1e12:.2f} TB margin, "
                         f"{free / 1e12:.2f} TB free under {paths.source}. Free space, or lower MIN_FREE_GB.")
    download_files(store, items, local_dir, prefix, logger, threads)
    present = set(read_json(os.path.join(local_dir, "blocks_present.json"), []))
    present.update(parts.keys())
    write_local_metadata(paths, pop, meta, present)
    logger.info(f"{pop}: {len(present)} blocks present locally in {local_dir}")
    return {"pop": pop, "blocks": len(parts), "bytes": total}


def clean_blocks(paths, pop, logger):
    local_dir = paths.block_matrix(pop)
    if os.path.isdir(local_dir):
        shutil.rmtree(local_dir)
        logger.info(f"{pop}: removed downloaded LD blocks {local_dir}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["tables", "plan", "blocks", "clean-blocks"])
    parser.add_argument("--pop", required=True)
    parser.add_argument("--chroms", nargs="*", default=None, help="chromosomes (default: all autosomes)")
    parser.add_argument("--region", default=None, help="CHR:START-END for a dense extraction (phase 0)")
    parser.add_argument("--radius-kb", type=float, default=None)
    parser.add_argument("--threads", type=int, default=int(env("DOWNLOAD_THREADS", "16")))
    args = parser.parse_args()

    paths = Paths().make()
    logger = get_logger(f"fetch.{args.pop}", os.path.join(paths.logs, f"fetch_{args.pop}.log"))
    store = open_store()
    chroms = [str(c) for c in args.chroms] if args.chroms else AUTOSOMES
    if args.command == "tables":
        fetch_tables(store, paths, args.pop, logger, args.threads)
    elif args.command == "plan":
        result = fetch_blocks(store, paths, args.pop, logger, args.threads, chroms=chroms,
                              radius_kb=args.radius_kb, plan_only=True)
        write_json(os.path.join(paths.state, f"{args.pop}.download_plan.json"), result)
        print(f"{args.pop}: {result['blocks']} block files, {result['bytes'] / 1e12:.2f} TB to download")
    elif args.command == "blocks":
        fetch_blocks(store, paths, args.pop, logger, args.threads, chroms=chroms, region=args.region,
                     radius_kb=args.radius_kb)
    elif args.command == "clean-blocks":
        clean_blocks(paths, args.pop, logger)
    return 0


if __name__ == "__main__":
    sys.exit(main())
