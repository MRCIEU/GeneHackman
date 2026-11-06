# GRCh38 reference data (bundle 1.2.0)

Scripts that build the GRCh38 reference data GeneHackman uses from release 1.2.0, when outputs default to GRCh38:

| Script | Builds |
|---|---|
| `10_build_1000g_panel.sh` (+ `11_assign_rsids.R`) | `genomic_data/1000genomes/b38_dbsnp156/`: the 1000G 30x (native GRCh38) panel, same people and rules as the GRCh37 phase 3 panel |
| `20_import_pan_ukb_ldscores.py` | `LDSCORE/pan_ukb/`: Pan-UKB LD scores (HapMap3, six ancestries) for LDSC |
| `31_qtl_from_gpm.R` | `<qtl>/<dataset>/GRCh38/`: cis QTLs of eQTLGen, BrainMeta and UKB-PPP (MR instruments and coloc regions), from the genotype-phenotype map 1.0.0 study data |
| `30_lift_qtl_data.R` | `<qtl>/<dataset>/GRCh37/`: the same QTL data lifted to GRCh37, with GeneHackman's liftover QC (optional; not published) |
| `40_validate_grch38.R` | `validation_report.md`: pack keys vs native GRCh38, liftover vs native positions, panel continuity, bundle contents |
| `50_assemble_bundle.sh` | the 1.2.0 bundle (hard links), and prints the upload commands. Users download it with [`download_reference_data.sh`](../../download_reference_data.sh) |

Settings: `config.sh`, overridden in `local_config.sh` (not tracked). Step-by-step instructions, with what to expect at each step: [`RUNBOOK.md`](RUNBOOK.md). The Pan-UKB LD pack itself is built by [`../pan_ukb_ld`](../pan_ukb_ld/README.md).
