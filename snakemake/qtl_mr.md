# qtl_mr

Standardise, clump, and fine-map **one** outcome GWAS, run Mendelian randomization against a QTL panel, plot MR results, and run  colocalisation on significant MR hits.

**Workflow file:** [`qtl_mr.smk`](qtl_mr.smk)  
**Example input:** [`input_templates/qtl_mr.yaml`](input_templates/qtl_mr.yaml)

## Run

```bash
./run_pipeline.sh snakemake/qtl_mr.smk path/to/input.yaml
```

Requires **`QTL_DATA_DIR`** in `.env` (or the default layout under `PIPELINE_DATA_DIR/qtl_datasets/`). See [`.env_example`](../.env_example) and [PLATFORM_SETUP.md](../PLATFORM_SETUP.md).

## Input

### GWAS count

- **Exactly one** outcome GWAS under `gwases:`.

### Per-GWAS fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | Yes | Outcome GWAS summary statistics. |
| `ancestry` | **Yes** | LD reference for clumping and fine-mapping. |
| `N` | **Yes** | Sample size (used in MR coloc LBF conversion). |
| `study_type` | No | `continuous` (default) or `categorical` — controls QTL LBF conversion in the coloc step. |
| `columns` | No | Column map or preset. |
| `populate_rsid` / `populate_eaf` | No | See [PIPELINES.md](../PIPELINES.md#shared-yaml-schema). |

### Root-level fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | PLINK clump settings. |
| `finemap.*` | No | Same keys as [finemap.md](finemap.md#root-level-fields). |
| `qtl.dataset` | **Yes** | QTL source, e.g. `metabrain`, `eqtlgen`. |
| `qtl.subcategory` | **Yes** | Sub-panel, e.g. `cortex` (MetaBrain) or `cis` (eQTLGen). Combined with `dataset` in output names. |
| `qtl.exposures` | No | List of exposure names to test; if empty, all exposures in the panel are run. |

### QTL dataset examples

| `dataset` | Example `subcategory` values |
| --------- | ------------------------------ |
| `metabrain` | `basalganglia`, `cerebellum`, `cortex`, `hippocampus`, `spinalcord` |
| `eqtlgen` | `cis` |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise outcome GWAS
2. PLINK clump
3. SuSiE RSS fine-mapping (`run_finemap.R`)
4. MR vs QTL panel (`run_mr_against_qtl_data.R`)
5. Volcano plot of MR results
6.  coloc on significant MR results (GWAS LBF vs QTL LBF from `convert_z_to_lbf`)
7. HTML report

## Output

Let `<prefix>` = outcome GWAS file prefix, `<qtl>` = `{dataset}_{subcategory}`.

### Data (`data/`)

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | Standardised outcome GWAS. |
| `clumped_snps/<prefix>.clumped` | Clumped leads. |

### Results (`results/`)

| Path | Description |
| ---- | ----------- |
| `finemap/<prefix>/finemap_complete_<prefix>.txt` | Fine-mapping completion marker. |
| `finemap/<prefix>/<CHR>_<BP>_finemap.tsv.gz` | Per-locus SuSiE output. |
| `mr/<prefix>_<qtl>.tsv.gz` | MR results vs QTL panel. |
| `mr/coloc_<prefix>_<qtl>.tsv` |  coloc for significant MR hits. |
| `mr/result_<qtl>_<prefix>.html` | HTML report. |
| `plots/volcano_plot_<prefix>_<qtl>.png` | Volcano plot of MR results. |

## Notes

- E2E tests for this pipeline run only when `QTL_DATA_DIR` is set (see [CONTRIBUTING.md](../CONTRIBUTING.md)).
- Fine-mapping uses the same single-ancestry SuSiE path as [`finemap.smk`](finemap.smk) and [`coloc.smk`](coloc.smk).
