"""Tests for scripts/ld_reference.py (run: python -m unittest discover tests/python)."""

import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
sys.path.insert(0, os.path.join(REPO, "reference_data", "pan_ukb_ld"))

from ld_reference import (LDPack, beyond_band_fraction, compute_seg_len,  # noqa: E402
                          make_positive_semidefinite, resolve_ancestry)
from ldpack_writer import ChromosomeWriter, quantise  # noqa: E402

PACK = os.path.join(REPO, "tests", "testthat", "data", "ld_pack_tiny")
RADIUS = 40_000
TOLERANCE = 0.5 / 127 + 1e-6


def expected_matrix(pack, pop, chrom):
    variants = pack.variants(chrom, pop)
    truth = pd.read_csv(os.path.join(PACK, f"truth_{pop}_chr{chrom}.tsv.gz"), sep="\t", header=None).values
    bp = variants["BP37"].values
    expected = np.where(np.abs(bp[:, None] - bp[None, :]) <= RADIUS, truth, 0.0)
    expected[np.abs(expected) < 0.02] = 0
    np.fill_diagonal(expected, 1)
    return expected


class TestLDPackReader(unittest.TestCase):
    def setUp(self):
        self.pack = LDPack(PACK)

    def test_matches_source_within_quantisation(self):
        for pop in ("EUR", "CSA"):
            for chrom in ("21", "22"):
                variants = self.pack.variants(chrom, pop)
                snps = list(variants["SNP37"].values[::-1])
                r_matrix, found = self.pack.get_ld_matrix(snps, chrom, pop, "GRCh37")
                self.assertEqual(found, snps)
                self.assertTrue(np.allclose(r_matrix, r_matrix.T))
                idx = [int(np.flatnonzero(variants["SNP37"].values == s)[0]) for s in found]
                expected = expected_matrix(self.pack, pop, chrom)[np.ix_(idx, idx)]
                # the truth file is rounded to 6 decimals, so skip values sitting on the dead-zone edge
                truth = pd.read_csv(os.path.join(PACK, f"truth_{pop}_chr{chrom}.tsv.gz"), sep="\t",
                                    header=None).values[np.ix_(idx, idx)]
                clear = np.abs(np.abs(truth) - 0.02) > 1e-5
                self.assertLessEqual(np.abs(r_matrix - expected)[clear].max(), TOLERANCE)

    def test_grch38_keys_skip_unmapped(self):
        variants = self.pack.variants("22", "EUR")
        mapped = variants[variants["LIFT_STATUS"] == "ok"]
        unmapped = variants[variants["LIFT_STATUS"] != "ok"]
        wanted = list(mapped["SNP38"].values[:5]) + ["22:1_A_C", unmapped["SNP37"].values[0]]
        r38, found = self.pack.get_ld_matrix(wanted, "22", "EUR", "GRCh38")
        self.assertEqual(found, list(mapped["SNP38"].values[:5]))
        r37, _ = self.pack.get_ld_matrix(list(mapped["SNP37"].values[:5]), "22", "EUR", "GRCh37")
        np.testing.assert_array_equal(r38, r37)

    def test_sas_alias_and_band_fraction(self):
        self.assertEqual(resolve_ancestry("sas"), "CSA")
        variants = self.pack.variants("21", "SAS")
        self.assertEqual(beyond_band_fraction(self.pack, variants["SNP37"].values[:3], "21", "SAS", "GRCh37"), 0.0)
        far = [variants["SNP37"].values[0], variants["SNP37"].values[-1]]
        self.assertEqual(beyond_band_fraction(self.pack, far, "21", "SAS", "GRCh37"), 1.0)

    def test_positive_semidefinite_fix(self):
        fixed, ridge = make_positive_semidefinite(np.array([[1, .9, 0], [.9, 1, .9], [0, .9, 1]]))
        self.assertGreater(ridge, 0)
        self.assertGreaterEqual(np.linalg.eigvalsh(fixed).min(), -1e-9)
        _, no_ridge = make_positive_semidefinite(np.eye(3))
        self.assertEqual(no_ridge, 0.0)

    @unittest.skipUnless(shutil.which(os.environ.get("GH_RSCRIPT", "Rscript")),
                         "Rscript not available for the R/Python parity check")
    def test_identical_to_r_reader(self):
        variants = self.pack.variants("22", "EUR")
        snps = list(variants["SNP38"].dropna().values[::2])
        r_py, found = self.pack.get_ld_matrix(snps, "22", "EUR", "GRCh38")
        out = tempfile.mktemp(suffix=".tsv")
        code = (
            "suppressMessages(pkgload::load_all('{repo}', quiet = TRUE));"
            "res <- get_ld_matrix_from_pack(readLines('{snps}'), '22', 'EUR', 'GRCh38', '{pack}');"
            "data.table::fwrite(as.data.frame(res$matrix), '{out}', sep = '\\t', col.names = FALSE)"
        )
        snp_file = tempfile.mktemp(suffix=".txt")
        with open(snp_file, "w", encoding="utf-8") as handle:
            handle.write("\n".join(snps) + "\n")
        subprocess.run([os.environ.get("GH_RSCRIPT", "Rscript"), "-e",
                        code.format(repo=REPO, snps=snp_file, pack=PACK, out=out)], check=True)
        r_r = pd.read_csv(out, sep="\t", header=None).values
        np.testing.assert_allclose(r_py, r_r, rtol=0, atol=1e-12)


class TestLDPackWriter(unittest.TestCase):
    def test_seg_len_and_quantise(self):
        np.testing.assert_array_equal(compute_seg_len([100, 200, 300, 1_000], 200), [2, 1, 0, 0])
        np.testing.assert_array_equal(quantise([0.5, -1.2, 0.01, np.nan], deadzone=0.02), [64, -127, 0, 0])

    def test_resume_after_interruption(self):
        work = tempfile.mkdtemp()
        rng = np.random.default_rng(1)
        n = 30
        variants = pd.DataFrame({
            "CHR": "1", "BP37": np.arange(n) * 1000 + 10_000, "EA": "A", "OA": "G", "EAF": 0.3,
            "SNP37": [f"1:{i * 1000 + 10_000}_A_G" for i in range(n)], "SNP38": None, "BP38": None,
            "LIFT_STATUS": "unmapped", "REF37": "A", "ALT37": "G", "SIGN": -1})
        r_full = np.corrcoef(rng.standard_normal((50, n)).T)
        pop_dir = os.path.join(work, "pack", "EUR")

        writer = ChromosomeWriter(pop_dir, os.path.join(work, "w"), "1", variants, 5_000, chunk_rows=8)
        first = writer.pending_chunks()[0]
        r0, r1 = writer.bounds[first]
        writer.write_chunk(first, r_full[r0:r1, r0:writer.column_end(first)])

        resumed = ChromosomeWriter(pop_dir, os.path.join(work, "w"), "1", variants, 5_000, chunk_rows=8)
        self.assertEqual(len(resumed.pending_chunks()), len(resumed.bounds) - 1)
        for k in resumed.pending_chunks():
            r0, r1 = resumed.bounds[k]
            resumed.write_chunk(k, r_full[r0:r1, r0:resumed.column_end(k)])
        totals = resumed.finalise()
        self.assertEqual(totals["n_variants"], n)

        with self.assertRaises(ValueError):
            ChromosomeWriter(pop_dir, os.path.join(work, "w"), "1", variants, 10_000, chunk_rows=8)
        shutil.rmtree(work)


if __name__ == "__main__":
    unittest.main()
