"""
Write a small synthetic copy of the Pan-UKB LD release, in the same Hail formats and layout, to test
the pack build end to end without network access:

    make_synthetic_source.py <out_dir>        (run with the build environment's python: needs Hail)

<out_dir>/UKBB.<POP>.ldadj.variant.ht   locus<GRCh37>, alleles, idx; globals n_samples, pop
<out_dir>/UKBB.<POP>.ldscore.ht         locus, alleles, rsid, varid, AF, ld_score
<out_dir>/UKBB.<POP>.ldadj.bm           r on ALT dosages, upper triangle, banded (sparse blocks)
<out_dir>/truth_<POP>.npz               the full r matrix and variant list, for verify_synthetic_pack.py

Includes the awkward cases the build must handle: rare variants (filtered by MAF), a multi-allelic
row, two biallelic variants at one position, variants with no AF, chrX, and alleles where ALT sorts
before REF (so the EA orientation flips the sign of r).
"""

import os
import sys

import numpy as np
import pandas as pd

BASES = np.array(["A", "C", "G", "T"])
# GRCh37 stretches that lift over cleanly to GRCh38
LAYOUT = [("21", 33_000_000, 260), ("22", 20_000_000, 340), ("X", 40_000_000, 20)]
BLOCK_SIZE = 64
BAND_INDICES = 160


def simulate(rng, n_samples):
    rows = []
    for contig, start, n in LAYOUT:
        bp = start + np.cumsum(rng.integers(300, 6000, size=n))
        for i in range(n):
            ref = rng.choice(BASES)
            alt = rng.choice(BASES[BASES != ref])
            rows.append({"contig": contig, "position": int(bp[i]), "alleles": [str(ref), str(alt)]})
    frame = pd.DataFrame(rows)
    # two biallelic variants at one position, and one multi-allelic row
    k = 40
    frame.at[k + 1, "position"] = frame.at[k, "position"]
    ref = frame.at[k, "alleles"][0]
    other = [b for b in BASES if b not in frame.at[k, "alleles"]][0]
    frame.at[k + 1, "alleles"] = [ref, str(other)]
    frame.at[k + 5, "alleles"] = frame.at[k + 5, "alleles"] + ["T" if "T" not in frame.at[k + 5, "alleles"] else "A"]

    n_variants = len(frame)
    latent = rng.standard_normal((n_samples, n_variants))
    for j in range(1, n_variants):
        latent[:, j] = 0.9 * latent[:, j - 1] + np.sqrt(1 - 0.81) * latent[:, j]
    freq = rng.uniform(0.02, 0.5, size=n_variants)
    rare = rng.random(n_variants) < 0.12
    freq[rare] = rng.uniform(0.002, 0.008, size=rare.sum())
    hap1 = np.zeros((n_samples, n_variants))
    hap2 = np.zeros((n_samples, n_variants))
    for j in range(n_variants):
        threshold = np.quantile(latent[:, j], 1 - freq[j])
        hap1[:, j] = latent[:, j] > threshold
        hap2[:, j] = np.roll(latent[:, j], 7) > threshold
    dosage = hap1 + hap2
    dosage[0, :] = np.where(dosage[:, :].sum(axis=0) == 0, 1, dosage[0, :])  # no monomorphic columns
    r = np.corrcoef(dosage.T)
    frame["AF"] = dosage.mean(axis=0) / 2
    frame["idx"] = np.arange(n_variants)
    return frame, r


def main(out_dir):
    import hail as hl
    from hail.linalg import BlockMatrix
    hl.init(master="local[2]", quiet=True, tmp_dir=os.path.join(out_dir, "tmp"), default_reference="GRCh37",
            log=os.path.join(out_dir, "hail.log"))
    rng = np.random.default_rng(7)
    os.makedirs(out_dir, exist_ok=True)
    for pop, n_samples in (("EUR", 900), ("CSA", 300)):
        frame, r = simulate(rng, n_samples)
        records = [{"contig": row.contig, "position": int(row.position), "alleles": list(row.alleles),
                    "idx": int(row.idx), "AF": float(row.AF)} for row in frame.itertuples()]
        table = hl.Table.parallelize(records, hl.tstruct(contig=hl.tstr, position=hl.tint32,
                                                         alleles=hl.tarray(hl.tstr), idx=hl.tint64, AF=hl.tfloat64))
        table = table.annotate(locus=hl.locus(table.contig, table.position))
        index = table.select("locus", "alleles", "idx").key_by("locus", "alleles")
        index = index.annotate_globals(n_samples=n_samples, pop=pop)
        index.write(os.path.join(out_dir, f"UKBB.{pop}.ldadj.variant.ht"), overwrite=True)

        no_af = set(rng.choice(len(frame), size=4, replace=False).tolist())
        ldscore = table.annotate(rsid=hl.str("rs") + hl.str(table.idx + 1000), varid=hl.str(table.idx),
                                 AF=hl.if_else(hl.literal(no_af).contains(hl.int32(table.idx)),
                                               hl.missing(hl.tfloat64), table.AF),
                                 ld_score=hl.float64(1.0))
        ldscore = ldscore.select("locus", "alleles", "rsid", "varid", "AF", "ld_score").key_by("locus", "alleles")
        ldscore.write(os.path.join(out_dir, f"UKBB.{pop}.ldscore.ht"), overwrite=True)

        bm = BlockMatrix.from_numpy(r, block_size=BLOCK_SIZE).sparsify_band(0, BAND_INDICES)
        bm.write(os.path.join(out_dir, f"UKBB.{pop}.ldadj.bm"), overwrite=True)
        np.savez(os.path.join(out_dir, f"truth_{pop}.npz"), r=r, contig=frame["contig"].values.astype(str),
                 position=frame["position"].values, alleles=np.array([",".join(a) for a in frame["alleles"]]),
                 af=frame["AF"].values, no_af=np.array(sorted(no_af)))
        print(f"{pop}: {len(frame)} variants written")


if __name__ == "__main__":
    main(sys.argv[1])
