"""MultiSuSiE with the Pan-UKB LD pack (run: python -m unittest discover tests/python; needs MultiSuSiE)."""

import os
import sys
import unittest

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
        self.assertEqual(credible_sets["TOP_SNP"].iloc[0], gwas_eur["SNP"].iloc[25])


if __name__ == "__main__":
    unittest.main()
