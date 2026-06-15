# Pipeline input and workflows

GeneHackman runs as Snakemake workflows (`.smk` files under `snakemake/`). Each workflow consumes a **YAML** input file and writes under **`PROJECT_DIR/data/`** and **`PROJECT_DIR/results/`** (set in `.env`).

**Run:**

```bash
./run_pipeline.sh snakemake/<workflow>.smk [path/to/input.yaml]
```

The second argument defaults to `input.yaml`. With bare `snakemake`, pass `--config genehackman_input=path/to/input.yaml`.

Copy starter YAML from [`snakemake/input_templates/`](snakemake/input_templates/). Validate syntax with [yamllint](https://www.yamllint.com/) if needed.

**Profiles:** Snakemake `--profile` directories live under `snakemake/profiles/` (see [PLATFORM_SETUP.md](PLATFORM_SETUP.md), [README.md](README.md)).

**Harmonisation background:** [README — How it works](README.md#how-it-works).

---

## Pipelines

Each workflow has a dedicated doc next to its `.smk` file with **input** and **output** details.

| Workflow | Doc | Role |
| -------- | --- | ---- |
| [`standardise_gwas.smk`](snakemake/standardise_gwas.smk) | [standardise_gwas.md](snakemake/standardise_gwas.md) | Harmonise columns and reference build only. |
| [`compare_gwases.smk`](snakemake/compare_gwases.smk) | [compare_gwases.md](snakemake/compare_gwases.md) | Standardise → clump → heterogeneity, LDSC, expected vs observed. |
| [`disease_progression.smk`](snakemake/disease_progression.smk) | [disease_progression.md](snakemake/disease_progression.md) | Incident vs subsequent GWAS; collider bias correction. |
| [`finemap.smk`](snakemake/finemap.smk) | [finemap.md](snakemake/finemap.md) | Standardise → clump → SuSiE or MultiSuSiE fine-mapping. |
| [`coloc.smk`](snakemake/coloc.smk) | [coloc.md](snakemake/coloc.md) | Finemap ≥2 GWAS → pairwise coloc. |
| [`qtl_mr.smk`](snakemake/qtl_mr.smk) | [qtl_mr.md](snakemake/qtl_mr.md) | One outcome GWAS → MR vs QTL panel → coloc on MR hits. |

---

## Shared YAML schema

These fields apply across workflows unless a pipeline doc says otherwise.

### `gwases` array

Each element describes one input GWAS file.

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Absolute path recommended. Supports VCF, CSV, TSV, TXT (`.gz` / `.zip`). |
| `columns` | No | Maps logical names (`CHR`, `BETA`, …) to file headers, or a preset string (`metal`, `gwama`, `opengwas_phewas`, …). See [`inst/extdata/predefined_column_maps.csv`](inst/extdata/predefined_column_maps.csv). |
| `N` | Varies | Sample size. **Required** for `compare_gwases` (LDSC) and `qtl_mr`. Used by fine-mapping where applicable. |
| `ancestry` | Varies | `AFR`, `AMR`, `EAS`, `EUR`, `SAS`. **Required** when the pipeline clumps or fine-maps, and whenever `populate_eaf: true`. |
| `build` | No | Input build: `GRCh36`, `GRCh37`, `GRCh38`. Default `GRCh37`. |
| `populate_rsid` | No | Per-GWAS override of root `populate_rsid`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |
| `flip_alleles` | No | Inherits root default `true`. **`false` allowed only for [`standardise_gwas.smk`](snakemake/standardise_gwas.smk).** |
| `study_type` | No | `continuous` (default) or `categorical`. Used by `qtl_mr` coloc. |
| `remove_extra_columns` | No | Default `false`. |

**Column content** (after mapping):

- Required: `CHR`, `BP`, `EA`, `OA`
- Effect: `BETA`+`SE`, or `OR`+`OR_LB`+`OR_UB`, or `Z`
- P-value: `P` or `LOG_P`

Default logical column names: `SNP`, `CHR`, `BP`, `EA`, `OA`, `EAF`, `BETA`, `SE`, `OR`, `P`, `LOG_P`, `Z`, `OR_LB`, `OR_UB`, `RSID`, `N`, `N_CASES`, `ENSEMBL_ID`, `GENE_NAME`.

### Root-level fields

| Field | Default | Notes |
| ----- | ------- | ----- |
| `is_test` | `false` | Tests / logging. |
| `populate_rsid` | `false` | Pipeline-wide RSID default. With clumping and `false` everywhere, RSID population is **partial** unless overridden. |
| `populate_eaf` | `false` | Pipeline-wide EAF fill default. |
| `flip_alleles` | `true` | Pipeline-wide allele harmonisation. `false` only on `standardise_gwas.smk`. |
| `output.build` | `GRCh37` | Target build after liftover (standardisation step). |
| `output.columns` | default map | Rename standard columns on output. |
| `plink_clump_arguments` | — | **Required** on pipelines that clump. Passed to `plink1.9 --clump`. [PLINK clump options](https://zzz.bwh.harvard.edu/plink/clump.shtml). |

### Optional blocks

#### `finemap` (finemap, coloc, qtl_mr)

| Key | Meaning | Default |
| --- | ------- | ------- |
| `window_kb` | Half-width of fine-mapping window around each lead (kb) | `500` |
| `max_causal` | SuSiE / MultiSuSiE max causal signals (`L`) | `10` |
| `coverage` | Credible set coverage | `0.95` |
| `min_abs_corr` | Minimum \|r\| for credible-set purity | `0.5` |

#### `coloc` (coloc only)

| Key | Meaning | Default |
| --- | ------- | ------- |
| `overlap_kb` | Max lead–lead distance (kb) for overlapping signals | `1000` |
| `p1` | Prior: SNP associated with trait 1 | `1e-4` |
| `p2` | Prior: SNP associated with trait 2 | `1e-4` |
| `p12` | Prior: SNP associated with both traits | `5e-6` |

#### `qtl` (qtl_mr only)

| Key | Meaning |
| --- | ------- |
| `dataset` | QTL source (`metabrain`, `eqtlgen`, …) |
| `subcategory` | Panel within dataset (`cortex`, `cis`, …) |
| `exposures` | Optional list of exposures; empty = all in panel |

#### `output.adjusted_gwas` (disease_progression only)

| Key | Meaning |
| --- | ------- |
| `type` | `slopehunter`, `cwls`, or `mr_ivw` (default `slopehunter`) |
| `p_val` | Threshold: `0.1`, `0.01`, `0.001`, `1e-5` (default `0.001`) |

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

## References

- Expected vs observed: [doi:10.1038/nature17671](https://doi.org/10.1038/nature17671)
- LDSC: [github.com/bulik/ldsc](https://github.com/bulik/ldsc)
- coloc: [chr1swallace.github.io/coloc](https://chr1swallace.github.io/coloc/)
- Pipeline DOI: [10.5281/zenodo.10624713](https://doi.org/10.5281/zenodo.10624713)
