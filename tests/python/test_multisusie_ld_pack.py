"""MultiSuSiE with the Pan-UKB LD pack (run: python -m unittest discover tests/python; needs MultiSuSiE)."""

import os
import shutil
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "scripts"))
PACK = os.path.join(REPO, "tests", "testthat", "data", "ld_pack_tiny")

try:
    import run_multisusie  # noqa: E402
except ImportError:  # MultiSuSiE not installed
    run_multisusie = None

from ld_reference import LDPack  # noqa: E402


def simulated_gwas(pack, pop, causal, build="GRCh37"):
    variants = pack.variants("22", pop)
    variants = variants[variants["LIFT_STATUS"] == "ok"].reset_index(drop=True)
    key = "SNP37" if build == "GRCh37" else "SNP38"
    r_matrix, found = pack.get_ld_matrix(list(variants[key]), "22", pop, build)
    z = r_matrix[:, causal] * 10
    return pd.DataFrame({"SNP": found, "CHR": "22", "BP": variants["BP37" if build == "GRCh37" else "BP38"].values,
                         "EA": variants["EA"], "OA": variants["OA"], "EAF": variants["EAF"],
                         "BETA": z * 0.01, "SE": 0.01, "P": 1e-3})


def write_bfile(prefix, genotypes, chrom, bp, a1, a2):
    """PLINK 1 binary files from an individuals x variants matrix of A1 counts (0, 1, 2)."""
    n, m = genotypes.shape
    codes = {2: 0b00, 1: 0b10, 0: 0b11}
    with open(prefix + ".bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        for j in range(m):
            packed = bytearray((n + 3) // 4)
            for i in range(n):
                packed[i // 4] |= codes[int(genotypes[i, j])] << (2 * (i % 4))
            f.write(bytes(packed))
    with open(prefix + ".bim", "w") as f:
        for j in range(m):
            f.write(f"{chrom}\tvariant{j}\t0\t{bp[j]}\t{a1[j]}\t{a2[j]}\n")
    with open(prefix + ".fam", "w") as f:
        for i in range(n):
            f.write(f"f{i} i{i} 0 0 0 -9\n")


def simulated_genotypes(n, m, seed=1):
    """A1 counts with LD between neighbouring variants (two haplotypes per person from an AR(1) latent)."""
    rng = np.random.default_rng(seed)
    haplotypes = []
    for _ in range(2):
        latent = np.zeros((n, m))
        latent[:, 0] = rng.standard_normal(n)
        for j in range(1, m):
            latent[:, j] = 0.9 * latent[:, j - 1] + np.sqrt(1 - 0.81) * rng.standard_normal(n)
        haplotypes.append((latent > rng.uniform(-0.8, 0.8, m)).astype(int))
    return haplotypes[0] + haplotypes[1]


@unittest.skipIf(run_multisusie is None, "MultiSuSiE is not installed")
class TestMultiSusieWithPack(unittest.TestCase):
    def test_locus_fine_maps_on_snp_keys_without_rsids(self):
        pack = LDPack(PACK)
        for build in ("GRCh37", "GRCh38"):
            gwas_eur = simulated_gwas(pack, "EUR", 25, build)
            gwas_sas = simulated_gwas(pack, "SAS", 25, build)  # SAS reads the pack's CSA data
            bp = gwas_eur["BP"].values
            locus = {"chr": "22", "start": int(bp[10]), "end": int(bp[40]), "center": int(bp[25])}
            options = {"ld_reference": "pan_ukb", "ld_pack_dir": PACK, "build": build,
                       "thousand_genomes_dir": None}
            result = run_multisusie.run_multisusie_for_locus(
                locus, [gwas_eur, gwas_sas], ["EUR", "SAS"], [50000, 20000], options,
                max_causal=3, coverage=0.95, min_abs_corr=0.5)
            self.assertIsNotNone(result)
            credible_sets, variants = result
            self.assertGreaterEqual(len(credible_sets), 1)
            self.assertIn(gwas_eur["SNP"].iloc[25], set(variants["SNP"]))
            self.assertEqual(credible_sets["TOP_SNP"].iloc[0], gwas_eur["SNP"].iloc[25])

    def test_rsids_are_reported_when_present(self):
        pack = LDPack(PACK)
        gwas_eur = simulated_gwas(pack, "EUR", 25)
        gwas_csa = simulated_gwas(pack, "CSA", 25)
        gwas_eur["RSID"] = [f"rs{i}" for i in range(len(gwas_eur))]
        bp = gwas_eur["BP"].values
        locus = {"chr": "22", "start": int(bp[10]), "end": int(bp[40]), "center": int(bp[25])}
        options = {"ld_reference": "pan_ukb", "ld_pack_dir": PACK, "build": "GRCh37", "thousand_genomes_dir": None}
        credible_sets, variants = run_multisusie.run_multisusie_for_locus(
            locus, [gwas_eur, gwas_csa], ["EUR", "CSA"], [50000, 20000], options, 3, 0.95, 0.5)
        self.assertEqual(credible_sets["TOP_SNP_RSID"].iloc[0], "rs25")
        self.assertEqual(credible_sets["LD_FLAG"].iloc[0], "OK")  # simulated from the pack's own LD
        self.assertIn("LD_MISMATCH_S_EUR", credible_sets.attrs["region"])
        self.assertEqual(credible_sets["TOP_SNP"].iloc[0], gwas_eur["SNP"].iloc[25])


@unittest.skipIf(run_multisusie is None or shutil.which("plink1.9") is None, "MultiSuSiE or plink1.9 not installed")
class TestMultiSusieMixedLdReferences(unittest.TestCase):
    """Each ancestry has its own LD reference; with the pack for one ancestry, variants are matched on SNP keys and
    PLINK genotypes (the user's own, or the 1000 Genomes panel) are read on position and alleles."""

    def setUp(self):
        self.dir = tempfile.mkdtemp()
        pack = LDPack(PACK)
        self.gwas_eur = simulated_gwas(pack, "EUR", 25)
        variants = pack.variants("22", "EUR")
        variants = variants[variants["LIFT_STATUS"] == "ok"].reset_index(drop=True)
        variants = variants[variants["SNP37"].isin(self.gwas_eur["SNP"])].reset_index(drop=True)
        genotypes = simulated_genotypes(400, len(variants))
        # .bim alleles in either order: r must still come out for the GWAS effect allele
        flip = np.arange(len(variants)) % 3 == 0
        a1 = np.where(flip, variants["OA"], variants["EA"])
        a2 = np.where(flip, variants["EA"], variants["OA"])
        self.prefix = os.path.join(self.dir, "EAS")
        write_bfile(self.prefix, genotypes, 22, variants["BP37"].astype(int).values, a1, a2)
        dosage_ea = np.where(flip, 2 - genotypes, genotypes).astype(float)
        r_eas = np.corrcoef(dosage_ea, rowvar=False)
        z = r_eas[:, 25] * 10
        self.gwas_eas = pd.DataFrame({"SNP": variants["SNP37"], "CHR": "22", "BP": variants["BP37"].values,
                                      "EA": variants["EA"], "OA": variants["OA"], "EAF": dosage_ea.mean(0) / 2,
                                      "BETA": z * 0.01, "SE": 0.01, "P": 1e-3})
        bp = self.gwas_eur["BP"].values
        self.locus = {"chr": "22", "start": int(bp[10]), "end": int(bp[40]), "center": int(bp[25])}

    def tearDown(self):
        shutil.rmtree(self.dir)

    def run_locus(self, sources):
        options = {"ld_references": sources, "ld_pack_dir": PACK, "build": "GRCh37", "thousand_genomes_dir": self.dir}
        return run_multisusie.run_multisusie_for_locus(
            self.locus, [self.gwas_eur, self.gwas_eas], ["EUR", "EAS"], [50000, 20000], options, 3, 0.95, 0.5)

    def test_pack_with_users_genotypes(self):
        result = self.run_locus(["pan_ukb", self.prefix])
        self.assertIsNotNone(result)
        credible_sets, variants = result
        self.assertEqual(credible_sets["TOP_SNP"].iloc[0], self.gwas_eur["SNP"].iloc[25])

    def test_pack_with_1000_genomes_read_on_position(self):
        ld = run_multisusie.compute_ld_from_reference(list(self.gwas_eas["SNP"]), "22", "EAS",
                                                      {"thousand_genomes_dir": self.dir}, "1000genomes",
                                                      ids_are_keys=True)
        self.assertEqual(len(ld), 2)  # matched on SNP keys, so no panel alleles to align
        self.assertIsNotNone(self.run_locus(["pan_ukb", "1000genomes"]))

    def test_labels_and_sizes(self):
        self.assertEqual(run_multisusie.ld_reference_label(self.prefix), "user:EAS")
        self.assertEqual(run_multisusie.ld_reference_label("pan_ukb"), "pan_ukb")
        self.assertEqual(run_multisusie.ld_reference_n_samples("EAS", {}, self.prefix), 400)


if __name__ == "__main__":
    unittest.main()
