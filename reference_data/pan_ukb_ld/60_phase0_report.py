"""
Stage 60: phase 0 report — measured counts, sizes and run times, extrapolated to the full pack.

    60_phase0_report.py            writes phase0/report.md and phase0/report.json

Uses, per population: raw/<POP>/summary.json (variant counts, AF coverage, n_samples),
raw/<POP>/annotate_summary.tsv (liftover outcomes), phase0/measure_<POP>_chr<C>.jsonl (compressed
bytes per stored value under each encoding, Hail fetch time per chunk), state/<POP>.download_plan.json
(LD blocks to download for the full pack) and, if present, phase0/<POP>_region.npz (LD beyond the
stored radius).  Stored-pair counts for every chromosome come from the annotated variant tables, so
the genome-wide size estimate is (bytes per value on the phase 0 chromosome) x (genome-wide pairs).
"""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

from build_common import AUTOSOMES, Paths, env, get_logger, read_json, write_json
from ld_reference import compute_seg_len


def genome_pairs(paths, pop, radius_bp):
    total, per_chrom = 0, {}
    for chrom in AUTOSOMES:
        path = os.path.join(paths.raw_pop(pop), f"variants_annot_chr{chrom}.tsv.gz")
        if not os.path.exists(path):
            continue
        bp = pd.read_csv(path, sep="\t", usecols=["BP37"])["BP37"].sort_values().values
        per_chrom[chrom] = int(compute_seg_len(bp, radius_bp).sum())
        total += per_chrom[chrom]
    return total, per_chrom


def region_beyond_radius(npz_path, radius_bp):
    data = np.load(npz_path, allow_pickle=False)
    r, bp = data["r"].astype(np.float64), data["bp37"].astype(np.int64)
    distance = np.abs(bp[:, None] - bp[None, :])
    upper = np.triu(np.ones_like(r, dtype=bool), 1)
    beyond = upper & (distance > radius_bp)
    return {"n_variants": int(bp.size), "pairs_beyond_radius": int(beyond.sum()),
            "max_abs_r_beyond_radius": float(np.abs(r[beyond]).max()) if beyond.any() else 0.0,
            "pairs_beyond_radius_abs_r_gt_0.1": int((np.abs(r[beyond]) > 0.1).sum()),
            "pairs_beyond_radius_abs_r_gt_0.3": int((np.abs(r[beyond]) > 0.3).sum()),
            "share_of_r2_mass_beyond_radius": float((r[beyond] ** 2).sum() / max(1e-12, (r[upper] ** 2).sum()))}


def main():
    paths = Paths().make()
    logger = get_logger("phase0_report", os.path.join(paths.logs, "phase0_report.log"))
    radius_bp = int(float(env("RADIUS_KB", "500")) * 1000)
    deadzone = float(env("DEADZONE", "0.02"))
    pops = env("POPS", "EUR AFR CSA EAS MID AMR").split()
    report = {"radius_kb": radius_bp / 1000, "deadzone": deadzone, "maf_min": float(env("MAF_MIN", "0.01")),
              "populations": {}}
    for pop in pops:
        summary = read_json(os.path.join(paths.raw_pop(pop), "summary.json"))
        if summary is None:
            continue
        entry = {"n_samples": summary.get("n_samples"), "n_index_variants": summary.get("n_index_variants"),
                 "af_coverage": summary.get("af_coverage"), "n_kept_maf": summary.get("n_kept"),
                 "maf_bins": summary.get("maf_bins")}
        annotate = os.path.join(paths.raw_pop(pop), "annotate_summary.tsv")
        if os.path.exists(annotate):
            entry["liftover"] = dict(pd.read_csv(annotate, sep="\t").values.tolist())
        pairs, _ = genome_pairs(paths, pop, radius_bp)
        entry["genome_stored_pairs"] = pairs
        measures = [json.loads(line) for f in glob.glob(os.path.join(paths.phase0, f"measure_{pop}_chr*.jsonl"))
                    for line in open(f, encoding="utf-8")]
        if measures:
            frame = pd.DataFrame(measures).drop_duplicates(["chrom", "chunk"], keep="last")
            n_values = frame["n_values"].sum()
            encodings = [c for c in frame.columns if c.startswith(("xz_dz", "gzip_dz"))]
            entry["phase0_bytes_per_value"] = {e: float(frame[e].sum() / n_values) for e in encodings}
            entry["estimated_pack_gb"] = {e: pairs * v / 1e9 for e, v in entry["phase0_bytes_per_value"].items()}
            entry["zero_fraction"] = {c: float((frame[c] * frame["n_values"]).sum() / n_values)
                                      for c in frame.columns if c.startswith("zero_fraction")}
            seconds_per_value = frame["fetch_seconds"].sum() / n_values
            entry["estimated_extract_hours"] = pairs * seconds_per_value / 3600
            entry["n_nan_values"] = int(frame["n_nan"].sum())
        plan = read_json(os.path.join(paths.state, f"{pop}.download_plan.json"))
        if plan:
            entry["download_tb"] = plan["bytes"] / 1e12
        validate = read_json(os.path.join(paths.state, f"{pop}.validate.json"))
        if validate:
            entry["validation_ok"] = validate.get("ok")
        npz = os.path.join(paths.phase0, f"{pop}_region.npz")
        if os.path.exists(npz):
            entry["region_check"] = region_beyond_radius(npz, radius_bp)
        report["populations"][pop] = entry
    write_json(os.path.join(paths.phase0, "report.json"), report)

    lines = ["# Pan-UKB LD pack: phase 0 report", "",
             f"Radius {radius_bp // 1000} kb, MAF >= {report['maf_min']}, dead-zone {deadzone}.", "",
             "| Pop | LD n | kept variants (MAF) | AF coverage | stored pairs (genome) | est. pack GB (xz, dead-zone) "
             "| est. download TB | est. extract h | validated |",
             "|---|---|---|---|---|---|---|---|---|"]
    key = f"xz_dz{deadzone}"

    def fmt(value, spec):
        return "" if value is None else format(value, spec)

    for pop, e in report["populations"].items():
        cells = [pop, str(e.get("n_samples")), fmt(e.get("n_kept_maf"), ","), fmt(e.get("af_coverage"), ".3f"),
                 fmt(e.get("genome_stored_pairs"), ","), fmt(e.get("estimated_pack_gb", {}).get(key), ".1f"),
                 fmt(e.get("download_tb"), ".2f"), fmt(e.get("estimated_extract_hours"), ".0f"),
                 str(e.get("validation_ok", ""))]
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "Size by encoding (GB, genome-wide estimate):", ""]
    for pop, e in report["populations"].items():
        if "estimated_pack_gb" in e:
            lines.append(f"- {pop}: " + ", ".join(f"{k} {v:.1f}" for k, v in sorted(e["estimated_pack_gb"].items())))
    for pop, e in report["populations"].items():
        if "region_check" in e:
            lines += ["", f"LD beyond {radius_bp // 1000} kb in the {pop} test region "
                          f"{env('PHASE0_REGION', '')}: " + json.dumps(e["region_check"])]
    lines += ["", "Liftover outcomes and full numbers: report.json"]
    with open(os.path.join(paths.phase0, "report.md"), "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    logger.info(f"wrote {os.path.join(paths.phase0, 'report.md')}")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
