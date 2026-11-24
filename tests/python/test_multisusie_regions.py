"""Fine-mapping regions, L rule and user-supplied LD in run_multisusie.py (python -m unittest discover tests/python)."""

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

try:
    import run_multisusie  # noqa: E402
except ImportError:  # MultiSuSiE not installed
    run_multisusie = None


def leads(bps, chrom="1"):
    return pd.DataFrame({"SNP": [f"s{i}" for i in range(len(bps))], "CHR": chrom, "BP": bps})


@unittest.skipIf(run_multisusie is None, "MultiSuSiE is not installed")
class TestReadClumped(unittest.TestCase):
    def test_na_row_from_trailing_blank_lines_is_dropped_and_types_kept(self):
        tmp = tempfile.mkdtemp()
        try:
            path = os.path.join(tmp, "x.clumped")
            with open(path, "w") as f:
                f.write("CHR F SNP BP P TOTAL\n2 1 rs1 826282 1e-10 5\n16 1 rs2 868869 2e-9 3\nNA\n")
            df = run_multisusie.read_clumped(path)
            self.assertEqual(list(df["SNP"]), ["rs1", "rs2"])
            self.assertEqual(list(df["CHR"]), ["2", "16"])
            self.assertEqual(list(df["BP"]), [826282, 868869])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


@unittest.skipIf(run_multisusie is None, "MultiSuSiE is not installed")
class TestRegions(unittest.TestCase):
    def test_merge_is_capped_like_the_r_module(self):
        # the same case as test_finemap.R: leads every 400 kb, 500 kb windows, 1.5 Mb cap
        regions = run_multisusie.define_regions(leads([1e6, 1.4e6, 1.8e6, 2.2e6, 2.6e6]), 5e5, 1.5e6)
        self.assertTrue(all(r["end"] - r["start"] <= 1.5e6 for r in regions))
        self.assertEqual(sorted(s for r in regions for s in r["leads"]), [f"s{i}" for i in range(5)])
        for a, b in zip(regions, regions[1:]):
            self.assertLess(a["end"], b["start"])

    def test_none_keeps_one_window_per_lead(self):
        regions = run_multisusie.define_regions(leads([1e6, 1.5e6]), 5e5, 3e6, "none")
        self.assertEqual([(r["start"], r["end"]) for r in regions], [(500000, 1500000), (1000000, 2000000)])

    def test_tiles_overlap_and_assign_each_position_once(self):
        regions = run_multisusie.define_regions(leads(np.arange(2e6, 7.1e6, 1e6)), 1e6, 3e6, "tile")
        self.assertGreater(len(regions), 1)
        self.assertTrue(all(r["end"] - r["start"] <= 3e6 for r in regions))
        positions = np.arange(1.5e6, 7.5e6, 1e5)
        owners = np.zeros(len(positions), dtype=int)
        for i, r in enumerate(regions):
            inside = (positions >= r["start"]) & (positions <= r["end"])
            owners += run_multisusie.tile_core_mask(positions, regions, i) & inside
        self.assertTrue(np.all(owners == 1))

    def test_chromosomes_are_kept_apart(self):
        both = pd.concat([leads([1e6], "2"), leads([1e6], "10")])
        regions = run_multisusie.define_regions(both, 5e5)
        self.assertEqual([r["chr"] for r in regions], ["2", "10"])

    def test_max_causal_is_used_as_given(self):
        self.assertEqual(run_multisusie.resolve_max_causal("10"), 10)
        self.assertEqual(run_multisusie.resolve_max_causal(3), 3)
        for bad in ("auto", "0", "2.5"):
            with self.assertRaises(ValueError):
                run_multisusie.resolve_max_causal(bad)

    def test_small_ld_reference_gives_a_warning(self):
        warning = run_multisusie.ld_reference_size_warning("EUR", 443172, 503)
        self.assertIn("EUR: GWAS sample size (443,172) is 881 times that of the LD reference (503 individuals)", warning)
        self.assertEqual(run_multisusie.ld_reference_size_warning("EUR", 4000, 503), "")
        self.assertEqual(run_multisusie.ld_reference_size_warning("EAS", 4000, float("nan")), "")

    def test_locus_sample_size_prefers_n_eff(self):
        df = pd.DataFrame({"N": [10, 20, 30], "N_EFF": [4, 5, None]})
        self.assertEqual(run_multisusie.locus_sample_size(df, 99), 4.5)
        self.assertEqual(run_multisusie.locus_sample_size(df[["N"]], 99), 20)
        self.assertEqual(run_multisusie.locus_sample_size(pd.DataFrame({"B": [1]}), 99), 99)

    def test_allele_alignment(self):
        r = np.array([[1, 0.5], [0.5, 1]])
        gwas = pd.DataFrame({"EA": ["A", "C"], "OA": ["G", "T"]})
        aligned, matched = run_multisusie.align_to_gwas_alleles(r, gwas, ["A", "T"], ["G", "C"])
        self.assertTrue(matched.all())
        self.assertEqual(aligned[0, 1], -0.5)


KRIGING = os.path.join(REPO, "tests", "python", "data", "kriging")


@unittest.skipIf(run_multisusie is None, "MultiSuSiE is not installed")
class TestLdMismatch(unittest.TestCase):
    """The Python port of susieR's estimate_s_rss / kriging_rss against susieR 0.14.2's output (fixture made in R:
    40 simulated variants, one causal, variant 5's allele flipped, N = 50,000)."""

    def setUp(self):
        self.r = np.loadtxt(os.path.join(KRIGING, "R.tsv"))
        self.z = np.loadtxt(os.path.join(KRIGING, "z.txt"))
        self.expected = pd.read_csv(os.path.join(KRIGING, "k.tsv"), sep="\t")

    def test_matches_susier(self):
        s = run_multisusie.estimate_s_rss(self.z, self.r, 50000)
        self.assertAlmostEqual(s, 0.782257, places=5)
        z_adj, condmean, condvar, log_lr = run_multisusie.kriging_rss(self.z, self.r, 50000, s)
        np.testing.assert_allclose(condmean, self.expected["condmean"], atol=1e-5)
        np.testing.assert_allclose(condvar, self.expected["condvar"], atol=1e-5)
        # mixture weights by EM rather than mixsqp: logLR agrees closely, not exactly
        np.testing.assert_allclose(log_lr, self.expected["logLR"], atol=0.05)

    def test_flipped_allele_is_the_only_outlier(self):
        s, outlier = run_multisusie.ld_mismatch_qc(self.z, self.r, 50000)
        self.assertEqual(list(np.where(outlier)[0]), [4])
        self.assertGreater(s, 0.5)

    def test_flags(self):
        self.assertEqual(run_multisusie.ld_mismatch_flag(0.01, 0, 0), "OK")
        self.assertEqual(run_multisusie.ld_mismatch_flag(0.3, 2, 1, True), "OUTLIER_IN_CS;LD_MISMATCH;MHC")
        self.assertEqual(run_multisusie.ld_mismatch_flag(float("nan"), 1, 0), "OUTLIERS_IN_REGION")
        self.assertTrue(run_multisusie.overlaps_mhc("6", 30e6, 31e6))
        self.assertFalse(run_multisusie.overlaps_mhc("6", 40e6, 41e6))


@unittest.skipIf(run_multisusie is None or shutil.which("plink1.9") is None, "needs MultiSuSiE and plink1.9")
class TestUserLd(unittest.TestCase):
    def test_bfile_ld_is_oriented_to_the_gwas_effect_allele(self):
        rng = np.random.default_rng(4)
        n, p = 300, 5
        hap = lambda: np.cumsum(rng.random((n, p)) < 0.2, axis=1) % 2  # noqa: E731
        genotypes = hap() + hap()
        tmp = tempfile.mkdtemp()
        try:
            prefix = os.path.join(tmp, "user")
            a1, a2 = ["A", "C", "G", "A", "T"], ["G", "T", "T", "C", "C"]
            with open(prefix + ".map", "w") as f:
                for j in range(p):
                    f.write(f"1 v{j} 0 {(j + 1) * 1000}\n")
            with open(prefix + ".ped", "w") as f:
                for i in range(n):
                    alleles = []
                    for j in range(p):
                        g = genotypes[i, j]
                        alleles += [a1[j] if g >= 1 else a2[j], a1[j] if g == 2 else a2[j]]
                    f.write(f"f{i} i{i} 0 0 0 -9 " + " ".join(alleles) + "\n")
            subprocess.run(["plink1.9", "--file", prefix, "--make-bed", "--out", prefix],
                           capture_output=True, check=True)
            bim = pd.read_csv(prefix + ".bim", sep=r"\s+", header=None, dtype=str)
            ea = [min(x, y) for x, y in zip(bim[4], bim[5])]
            oa = [max(x, y) for x, y in zip(bim[4], bim[5])]
            keys = [f"1:{bp}_{e}_{o}" for bp, e, o in zip(bim[3], ea, oa)]
            r_matrix, found = run_multisusie.compute_ld_from_bfile(keys, "1", prefix)
            self.assertEqual(found, keys)
            dosage = np.column_stack([genotypes[:, j] if a1[j] == ea[j] else 2 - genotypes[:, j] for j in range(p)])
            np.testing.assert_allclose(r_matrix, np.corrcoef(dosage.T), atol=1e-4)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
