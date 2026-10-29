"""
Stage 40: write info.tsv for each population in the pack, the pack MANIFEST.json (settings,
per-population counts, sha256 of every file) and the pack README.

    40_finalise_pack.py                # every population with at least one finished chromosome
    40_finalise_pack.py --require-all  # fail unless every population has all 22 autosomes
"""

import argparse
import datetime
import hashlib
import os
import shutil
import subprocess
import sys

from build_common import AUTOSOMES, HERE, REPO, Paths, env, get_logger, read_json, write_json
from ldpack_writer import base_info, write_info
from ld_reference import FORMAT_NAME, FORMAT_VERSION

SOURCE_DESCRIPTION = ("Pan-UKB LD matrices (UKBB.<POP>.ldadj.bm), https://pan.ukbb.broadinstitute.org/docs/ld, "
                      "CC BY 4.0; usage may be restricted by UK Biobank")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit():
    try:
        return subprocess.check_output(["git", "-C", REPO, "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--require-all", action="store_true")
    parser.add_argument("--skip-checksums", action="store_true")
    args = parser.parse_args()

    paths = Paths().make()
    logger = get_logger("finalise", os.path.join(paths.logs, "finalise.log"))
    pops = env("POPS", "EUR AFR CSA EAS MID AMR").split()
    manifest = {"format": FORMAT_NAME, "format_version": FORMAT_VERSION,
                "created": datetime.datetime.now().isoformat(timespec="seconds"),
                "genehackman_commit": git_commit(), "source": SOURCE_DESCRIPTION,
                "settings": {"maf_min": float(env("MAF_MIN", "0.01")), "radius_kb": float(env("RADIUS_KB", "1000")),
                             "deadzone": float(env("DEADZONE", "0.02")), "chunk_rows": int(env("CHUNK_ROWS", "2048")),
                             "scale": 127, "compression": "xz", "index_build": "GRCh37",
                             "key_builds": ["GRCh37", "GRCh38"]},
                "populations": {}}
    for pop in pops:
        done = [c for c in AUTOSOMES if os.path.exists(os.path.join(paths.state, f"{pop}.ld.chr{c}.done"))]
        if not done:
            continue
        expected = [c for c in AUTOSOMES
                    if os.path.exists(os.path.join(paths.raw_pop(pop), f"variants_annot_chr{c}.tsv.gz"))]
        if args.require_all and sorted(done, key=int) != sorted(expected, key=int):
            raise SystemExit(f"{pop}: only {len(done)} of {len(expected)} chromosomes finished")
        summary = read_json(os.path.join(paths.raw_pop(pop), "summary.json"), {})
        n_samples = summary.get("n_samples")
        pop_dir = os.path.join(paths.pack, pop)
        info = base_info(n_samples, int(float(env("RADIUS_KB", "1000")) * 1000), float(env("MAF_MIN", "0.01")),
                         float(env("DEADZONE", "0.02")), f"Pan-UKB {pop} ldadj", int(env("CHUNK_ROWS", "2048")))
        info["chromosomes"] = ",".join(done)
        write_info(pop_dir, info)
        totals = {c: read_json(os.path.join(paths.state, f"{pop}.ld.chr{c}.json"), {}) for c in done}
        entry = {"n_samples": n_samples, "chromosomes": done,
                 "n_variants": sum(t.get("n_variants", 0) for t in totals.values()),
                 "ld_bytes": sum(t.get("bytes", 0) for t in totals.values()),
                 "n_stored_pairs": sum(t.get("n_values", 0) for t in totals.values()),
                 "source_variants_kept": summary.get("n_kept"), "af_coverage": summary.get("af_coverage")}
        if not args.skip_checksums:
            entry["sha256"] = {name: sha256(os.path.join(pop_dir, name)) for name in sorted(os.listdir(pop_dir))
                               if not name.startswith(".")}
        manifest["populations"][pop] = entry
        logger.info(f"{pop}: {len(done)} chromosomes, {entry['n_variants']:,} variants, "
                    f"{entry['ld_bytes'] / 1e9:.2f} GB")
    if not manifest["populations"]:
        raise SystemExit("No finished chromosomes to finalise")
    write_json(os.path.join(paths.pack, "MANIFEST.json"), manifest)
    shutil.copyfile(os.path.join(HERE, "PACK_README.md"), os.path.join(paths.pack, "README.md"))
    logger.info(f"pack finalised in {paths.pack}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
