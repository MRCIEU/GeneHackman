# Pipeline input and workflows

GeneHackman runs as Snakemake workflows (`.smk` files under `snakemake/`). Each workflow consumes a **YAML** input file and writes under **`PROJECT_DIR/data/`** and **`PROJECT_DIR/results/`** (set in `.env`).

**Run:**

```bash
./run_pipeline.sh snakemake/<workflow>.smk [path/to/input.yaml]
```

The second argument defaults to `input.yaml`. With bare `snakemake`, pass `--config genehackman_input=path/to/input.yaml`.

Starter YAML files live in [`snakemake/input_templates/`](snakemake/input_templates/). Each pipeline doc below inlines the template and explains **only the fields used by that workflow**.

**Profiles:** Snakemake `--profile` directories live under `snakemake/profiles/` (see [PLATFORM_SETUP.md](PLATFORM_SETUP.md), [README.md](README.md)).

**Harmonisation background:** [README — How it works](README.md#how-it-works).

---

## Pipelines

| Workflow | Doc | Role |
| -------- | --- | ---- |
| [`standardise_gwas.smk`](snakemake/standardise_gwas.smk) | [standardise_gwas.md](snakemake/standardise_gwas.md) | Harmonise columns and reference build only. |
| [`compare_gwases.smk`](snakemake/compare_gwases.smk) | [compare_gwases.md](snakemake/compare_gwases.md) | Standardise → clump → heterogeneity, LDSC, expected vs observed. |
| [`disease_progression.smk`](snakemake/disease_progression.smk) | [disease_progression.md](snakemake/disease_progression.md) | Incident vs subsequent GWAS; collider bias correction. |
| [`finemap.smk`](snakemake/finemap.smk) | [finemap.md](snakemake/finemap.md) | Standardise → clump → SuSiE or MultiSuSiE fine-mapping. |
| [`coloc.smk`](snakemake/coloc.smk) | [coloc.md](snakemake/coloc.md) | Finemap ≥2 GWAS → pairwise coloc. |
| [`qtl_mr.smk`](snakemake/qtl_mr.smk) | [qtl_mr.md](snakemake/qtl_mr.md) | One outcome GWAS → MR vs QTL panel → coloc on MR hits. |

---

## Common conventions

These apply to every pipeline unless a workflow doc says otherwise.

**File formats:** `gwases[].file` supports VCF, CSV, TSV, and TXT (`.gz` / `.zip`). Absolute paths are recommended.

**Column mapping:** set `columns` to a header map or a preset string (`metal`, `gwama`, `opengwas_phewas`, …). See [`inst/extdata/predefined_column_maps.csv`](inst/extdata/predefined_column_maps.csv).

**Minimum columns after mapping:**

- Required: `CHR`, `BP`, `EA`, `OA`
- Effect: `BETA`+`SE`, or `OR`+`OR_LB`+`OR_UB`, or `Z`
- P-value: `P` or `LOG_P`

Default logical names: `SNP`, `CHR`, `BP`, `EA`, `OA`, `EAF`, `BETA`, `SE`, `OR`, `P`, `LOG_P`, `Z`, `OR_LB`, `OR_UB`, `RSID`, `N`, `N_CASES`, `ENSEMBL_ID`, `GENE_NAME`.

**Allele harmonisation:** root `flip_alleles` defaults to `true`. Only [`standardise_gwas.smk`](snakemake/standardise_gwas.smk) accepts `flip_alleles: false`.

**Clumping:** pipelines that clump require `plink_clump_arguments`. If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus).

**LD reference:** `ld_reference` chooses where clumping, fine-mapping (SuSiE and MultiSuSiE), LocusZoom and EAF filling get LD:

| Value | Source | Variants are matched on | Ancestries |
|---|---|---|---|
| `1000genomes` (default) | 1000 Genomes Phase 3 genotypes, via PLINK | RSID | EUR, EAS, AFR, AMR, SAS |
| `pan_ukb` | Precomputed Pan-UKB LD: 6,000–420,000 people per ancestry, MAF ≥ 1%, ±500 kb, GRCh37 and GRCh38 | SNP key (`CHR:BP_EA_OA`) | EUR, AFR, AMR, CSA, EAS, MID (`SAS` maps to `CSA`) |

The Pan-UKB pack must be at `PIPELINE_DATA_DIR/genomic_data/pan_ukb_ld/`; see [`reference_data/pan_ukb_ld/README.md`](reference_data/pan_ukb_ld/README.md) for how it is built.

With `pan_ukb`:

- **Clumping** follows PLINK 1.9's `--clump` rules on the stored LD, and `--clump-kb` is capped at the pack's 500 kb radius.
- **Fine-mapping** variant pairs further apart than 500 kb are set to r = 0. `finemap_regions.tsv` reports the share of such pairs and any ridge added to keep the LD matrix positive semi-definite.
- **`compare_gwases`** still needs 1000G ancestry codes, because LDSC uses the 1000 Genomes LD scores.

**RSIDs and genome build:** RSIDs are looked up in the reference for the output build (`output.build`). Partial RSID population (the default for pipelines that clump) uses the GRCh37 1000 Genomes marker map, so for GRCh38 output it only works when the GWAS was lifted over from GRCh37; otherwise set `populate_rsid: true` to use the full dbSNP lookup for GRCh38.

---

## Standard intermediate paths

Most pipelines write harmonised GWAS and clumps to fixed locations:

| Path | Description |
| ---- | ----------- |
| `data/gwas/<prefix>_std.tsv.gz` | Standardised GWAS (`<prefix>` = input file basename). |
| `data/clumped_snps/<prefix>.clumped` | PLINK clump output (when the pipeline clumps). |

Fine-mapping completion markers:

| Path | Description |
| ---- | ----------- |
| `results/finemap/<prefix>/finemap_complete_<prefix>.txt` | Single-ancestry SuSiE done (one line = locus count). |
| `results/finemap/multi_ancestry/finemap_complete_<gwas_prefixes>.txt` | Multi-ancestry MultiSuSiE done; label is sorted GWAS prefixes joined with `_`. |

---

Pipeline DOI: [10.5281/zenodo.10624713](https://doi.org/10.5281/zenodo.10624713)
