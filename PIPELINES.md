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
