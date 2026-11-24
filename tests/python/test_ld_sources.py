"""Tests for snakemake/util/ld_sources.py, where each pipeline step gets LD (run: python -m pytest tests/python)."""

import os
import sys
import tempfile
import unittest

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "snakemake", "util"))

import ld_sources  # noqa: E402
from ld_sources import LDSourceError, resolve  # noqa: E402


class TestResolve(unittest.TestCase):
    def test_gwas_then_step_then_top_level(self):
        self.assertEqual(resolve("clumping", top="pan_ukb"), "pan_ukb")
        self.assertEqual(resolve("clumping", top="pan_ukb", step_value="1000genomes"), "1000genomes")
        self.assertEqual(resolve("clumping", top="pan_ukb", step_value="1000genomes", gwas_value="/g/x"), "/g/x")
        self.assertEqual(resolve("clumping", top=" pan_ukb ", step_value="", gwas_value=None), "pan_ukb")

    def test_nothing_set_is_empty(self):
        self.assertEqual(resolve("finemap"), "")
        self.assertIn("No LD reference is set for fine-mapping (GWAS a)", ld_sources.not_set_message("finemap", " (GWAS a)"))
        self.assertIn("finemap: {ld_reference: ...}", ld_sources.not_set_message("finemap"))

    def test_paths_fill_ancestry_and_environment(self):
        os.environ["LD_SOURCES_TEST_DIR"] = "/data/geno"
        self.assertEqual(resolve("finemap", top="$LD_SOURCES_TEST_DIR/cohort_{ancestry}", ancestry="EAS"),
                         "/data/geno/cohort_EAS")
        self.assertEqual(resolve("finemap", top="relative/prefix"), os.path.abspath("relative/prefix"))
        with self.assertRaises(LDSourceError):
            resolve("finemap", top="/data/cohort_{ancestry}", ancestry="")

    def test_keywords_are_not_paths(self):
        self.assertFalse(ld_sources.is_genotypes("pan_ukb"))
        self.assertFalse(ld_sources.is_genotypes("1000genomes"))
        self.assertFalse(ld_sources.is_genotypes(""))
        self.assertTrue(ld_sources.is_genotypes("panukb"))  # a typo is read as a path, and then not found


class TestGenotypeFiles(unittest.TestCase):
    def test_missing_files_are_listed(self):
        with tempfile.TemporaryDirectory() as d:
            prefix = os.path.join(d, "geno")
            with self.assertRaises(LDSourceError) as error:
                ld_sources.check_genotypes(prefix, "clumping", " (GWAS a)")
            self.assertIn(prefix + ".bed", str(error.exception))
            self.assertIn("pan_ukb", str(error.exception))
            for ext in (".bed", ".bim", ".fam"):
                open(prefix + ext, "w").close()
            ld_sources.check_genotypes(prefix, "clumping")
            ld_sources.check_genotypes("pan_ukb", "clumping")  # keywords aren't checked here

    def test_ld_bfile_message_points_to_ld_reference(self):
        self.assertIn("ld_reference on the GWAS", ld_sources.ld_bfile_message("a.tsv.gz"))


class TestDescribe(unittest.TestCase):
    def test_one_line_per_step_in_order(self):
        lines = ld_sources.describe({"finemap": {"a": "/g/x", "b": "pan_ukb"}, "clumping": {"a": "1000genomes",
                                                                                           "b": "1000genomes"}},
                                    ldsc_scores="pan_ukb")
        self.assertEqual(lines[0], "clumping: 1000 Genomes panel (1000genomes)")
        self.assertEqual(lines[1], "fine-mapping: a: your genotypes (/g/x); b: Pan-UK Biobank LD (pan_ukb)")
        self.assertTrue(lines[2].startswith("LDSC: pan_ukb LD scores"))


if __name__ == "__main__":
    unittest.main()
