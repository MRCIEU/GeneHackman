"""
Build the tiny LD pack fixture used by the R and Python unit tests.

    python reference_data/pan_ukb_ld/make_test_pack.py tests/testthat/data/ld_pack_tiny

It simulates genotypes for two populations (EUR and CSA) on chr21 and chr22, computes r on the
ALT allele, and writes a pack with a 40 kb radius and 16-row chunks.  Alongside the pack it writes
truth_<POP>_chr<C>.tsv.gz: the full float r matrix oriented to the canonical EA allele (rows and
columns in pack order), for the readers' tests.  Deterministic for a fixed seed.
"""

import json
import os
import shutil
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ldpack_writer import ChromosomeWriter, base_info, write_info  # noqa: E402

RADIUS_BP = 40_000
CHUNK_ROWS = 16
BASES = np.array(["A", "C", "G", "T"])


def canonical(chrom, bp, ref, alt):
    """Mirror GeneHackman's standardise_alleles: EA is the alphabetically smaller allele."""
    ea, oa = (alt, ref) if alt < ref else (ref, alt)
    sign = 1 if ea == alt else -1
    return ea, oa, sign, f"{chrom}:{bp}_{ea}_{oa}".upper()


def simulate_chromosome(rng, chrom, n_variants, n_samples, start_bp):
    # positions and alleles are the same in every population (seeded by chromosome), genotypes are not
    layout = np.random.default_rng(int(chrom))
    bp = start_bp + np.cumsum(layout.integers(500, 4000, size=n_variants))
    if n_variants > 10:
        bp[10] = bp[9]  # a multi-allelic site: two variants at the same position
    latent = rng.standard_normal((n_samples, n_variants))
    for j in range(1, n_variants):  # local LD that decays along the chromosome
        latent[:, j] = 0.85 * latent[:, j - 1] + np.sqrt(1 - 0.85 ** 2) * latent[:, j]
    freq = rng.uniform(0.05, 0.5, size=n_variants)
    threshold = np.array([np.quantile(latent[:, j], 1 - freq[j]) for j in range(n_variants)])
    hap1 = (latent > threshold).astype(float)
    hap2 = (np.roll(latent, 1, axis=0) > threshold).astype(float)
    dosage_alt = hap1 + hap2
    r_alt = np.corrcoef(dosage_alt.T)

    ref = layout.choice(BASES, size=n_variants)
    alt = np.array([layout.choice(BASES[BASES != r]) for r in ref])
    if n_variants > 10:  # second ALT allele at the multi-allelic site
        ref[10] = ref[9]
        alt[10] = [b for b in BASES if b not in (ref[9], alt[9])][0]

    rows = []
    for i in range(n_variants):
        ea, oa, sign, snp37 = canonical(chrom, bp[i], ref[i], alt[i])
        bp38 = int(bp[i] + 1000)
        snp38 = f"{chrom}:{bp38}_{ea}_{oa}"
        status = "ok"
        if i in (5, 6):
            bp38, snp38, status = None, None, "unmapped"
        eaf = dosage_alt[:, i].mean() / 2
        eaf = eaf if sign == 1 else 1 - eaf
        rows.append({"CHR": chrom, "BP37": int(bp[i]), "BP38": bp38, "EA": ea, "OA": oa,
                     "EAF": round(float(eaf), 4), "SNP37": snp37, "SNP38": snp38,
                     "LIFT_STATUS": status, "REF37": ref[i], "ALT37": alt[i],
                     "RSID_UKB": f"rs{chrom}{i:04d}", "SIGN": sign, "IDX_GLOBAL": i})
    return pd.DataFrame(rows), r_alt


def main(out_dir):
    rng = np.random.default_rng(20261028)
    for pop, n_samples in (("EUR", 400), ("CSA", 120)):
        pop_dir = os.path.join(out_dir, pop)
        work_dir = os.path.join(out_dir, ".work", pop)
        for chrom, n_variants, start in (("21", 40, 9_400_000), ("22", 60, 16_050_000)):
            variants, r_alt = simulate_chromosome(rng, chrom, n_variants, n_samples, start)
            writer = ChromosomeWriter(pop_dir, work_dir, chrom, variants, RADIUS_BP,
                                      chunk_rows=CHUNK_ROWS, deadzone=0.02)
            order = writer.variants["IDX_GLOBAL"].values
            r_sorted = r_alt[np.ix_(order, order)]
            for k in writer.pending_chunks():
                r0, _ = writer.bounds[k]
                c_end = writer.column_end(k)
                r1 = writer.bounds[k][1]
                writer.write_chunk(k, r_sorted[r0:r1, r0:c_end])
            writer.finalise()
            sign = writer.variants["SIGN"].values
            truth = r_sorted * np.outer(sign, sign)
            np.fill_diagonal(truth, 1.0)
            pd.DataFrame(truth).to_csv(os.path.join(out_dir, f"truth_{pop}_chr{chrom}.tsv.gz"),
                                       sep="\t", index=False, header=False, float_format="%.6f",
                                       compression="gzip")
        write_info(pop_dir, base_info(n_samples, RADIUS_BP, 0.01, 0.02, "synthetic test fixture",
                                      CHUNK_ROWS))
    with open(os.path.join(out_dir, "MANIFEST.json"), "w", encoding="utf-8") as handle:
        json.dump({"format": "genehackman-ldpack", "format_version": 1, "source": "synthetic test fixture",
                   "populations": {"EUR": {}, "CSA": {}}}, handle)
        handle.write("\n")
    # the work directory only holds resumable checkpoints; the fixture does not need it
    shutil.rmtree(os.path.join(out_dir, ".work"), ignore_errors=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "tests/testthat/data/ld_pack_tiny")
