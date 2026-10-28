"""
Check a pack built from make_synthetic_source.py against the truth it saved.

    verify_synthetic_pack.py <source_dir> <work_dir>     (RADIUS_KB, DEADZONE, MAF_MIN from the environment)

For every population and chromosome in the pack: the kept variants are exactly the biallelic,
autosomal variants with AF present and MAF >= MAF_MIN; EA/OA/EAF/SNP37 follow GeneHackman's
canonical allele order; and every decoded r equals the true r, oriented to EA, restricted to the
radius, with the dead-zone applied and quantised to int8 (exact equality).
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
sys.path.insert(0, os.path.join(HERE, ".."))
from ld_reference import LDPack  # noqa: E402
from ldpack_writer import quantise  # noqa: E402


def main(source_dir, work_dir):
    radius = int(float(os.environ.get("RADIUS_KB", "500")) * 1000)
    deadzone = float(os.environ.get("DEADZONE", "0.02"))
    maf_min = float(os.environ.get("MAF_MIN", "0.01"))
    pack = LDPack(os.path.join(work_dir, "pack"))
    failures = 0
    for pop in ("EUR", "CSA"):
        truth = np.load(os.path.join(source_dir, f"truth_{pop}.npz"))
        r_true = truth["r"]
        alleles = [a.split(",") for a in truth["alleles"]]
        no_af = set(truth["no_af"].tolist())
        for chrom in ("21", "22"):
            try:
                variants = pack.variants(chrom, pop)
            except FileNotFoundError:
                continue
            on_chrom = [i for i in range(len(alleles)) if truth["contig"][i] == chrom]
            expected_kept = [i for i in on_chrom if len(alleles[i]) == 2 and i not in no_af
                             and min(truth["af"][i], 1 - truth["af"][i]) >= maf_min]
            kept = []
            for row in variants.itertuples():
                ref, alt = row.REF37, row.ALT37
                match = [i for i in on_chrom if truth["position"][i] == row.BP37 and alleles[i] == [ref, alt]]
                kept.append(match[0])
                ea, oa = sorted([ref, alt])
                assert (row.EA, row.OA) == (ea, oa), (row.EA, row.OA, ea, oa)
                assert row.SNP37 == f"{chrom}:{row.BP37}_{ea}_{oa}"
                eaf = truth["af"][match[0]] if ea == alt else 1 - truth["af"][match[0]]
                assert abs(row.EAF - eaf) < 1e-5, (row.EAF, eaf)
            if sorted(kept) != expected_kept:
                print(f"FAIL {pop} chr{chrom}: kept variants differ from the MAF/biallelic/AF filter")
                failures += 1
            sign = np.array([1 if variants.EA.values[k] == variants.ALT37.values[k] else -1
                             for k in range(len(variants))])
            expected = r_true[np.ix_(kept, kept)] * np.outer(sign, sign)
            bp = variants["BP37"].values
            expected = np.where(np.abs(bp[:, None] - bp[None, :]) <= radius, expected, 0.0)
            expected = quantise(expected, 127, deadzone).astype(np.float64) / 127
            np.fill_diagonal(expected, 1.0)
            decoded = pack.matrix_for_rows(chrom, pop, np.arange(len(variants)))
            if not np.allclose(decoded, expected, atol=1e-12, rtol=0):
                print(f"FAIL {pop} chr{chrom}: decoded LD differs from truth "
                      f"(max {np.abs(decoded - expected).max():.4f})")
                failures += 1
            else:
                flipped = int((sign < 0).sum())
                print(f"OK   {pop} chr{chrom}: {len(variants)} variants ({flipped} with EA = REF), "
                      f"{int(variants['SEG_LEN'].sum())} stored pairs, LD identical to truth; "
                      f"liftover {variants['LIFT_STATUS'].value_counts().to_dict()}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
