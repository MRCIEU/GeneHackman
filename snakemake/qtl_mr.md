# qtl_mr

Standardise, clump, and fine-map **one** outcome GWAS, run Mendelian randomization against a QTL panel, plot MR results, and run colocalisation on significant MR hits.

**Workflow file:** [`qtl_mr.smk`](qtl_mr.smk)

## Run

```bash
./run_pipeline.sh snakemake/qtl_mr.smk path/to/input.yaml
```

Requires **`QTL_DATA_DIR`** in `.env` (or the default layout under `PIPELINE_DATA_DIR/qtl_datasets/`). See [`.env_example`](../.env_example) and [PLATFORM_SETUP.md](../PLATFORM_SETUP.md).

## Input

**GWAS count:** exactly one outcome GWAS under `gwases:`.

Copy from [`input_templates/qtl_mr.yaml`](input_templates/qtl_mr.yaml).

### Example YAML

```yaml
gwases:
  - file: /path/to/outcome_gwas.tsv.gz
    columns:
      CHR: CHR
      BP: BP
      EA: EA
      OA: OA
      EAF: EAF
      P: P
      BETA: BETA
      SE: SE
    ancestry: EUR
    N: 100000
    study_type: continuous
plink_clump_arguments: "--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"
finemap:
  window_kb: 1000
  max_causal: 10
  coverage: 0.95
  min_abs_corr: 0.5
qtl:
  dataset: metabrain
  subcategory: cortex
  exposures: []
populate_rsid: false
```

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Outcome GWAS summary statistics. |
| `ancestry` | **Yes** | LD reference for clumping and fine-mapping. |
| `N` | **Yes** | Sample size (used in MR coloc LBF conversion). |
| `study_type` | No | `continuous` (default) or `categorical` — controls QTL LBF conversion in the coloc step. |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump`. |

### `finemap`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `finemap.window_kb` | No | Half-width of fine-mapping window around each lead (kb). Default `500`. |
| `finemap.max_causal` | No | Max causal signals per locus. Default `10`. |
| `finemap.coverage` | No | Credible set coverage. Default `0.95`. |
| `finemap.min_abs_corr` | No | Minimum \|r\| for credible-set purity. Default `0.5`. |

### `qtl`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `qtl.dataset` | **Yes** | QTL source, e.g. `metabrain`, `eqtlgen`. |
| `qtl.subcategory` | **Yes** | Sub-panel, e.g. `cortex` (MetaBrain) or `cis` (eQTLGen). Combined with `dataset` in output names. |
| `qtl.exposures` | No | List of exposure names to test; if empty (`[]`), all exposures in the panel are run. |

**Example subcategories:**

| `dataset` | Example `subcategory` values |
| --------- | ------------------------------ |
| `metabrain` | `basalganglia`, `cerebellum`, `cortex`, `hippocampus`, `spinalcord` |
| `eqtlgen` | `cis` |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise outcome GWAS
2. PLINK clump
3. SuSiE RSS fine-mapping (`run_finemap.R`)
4. MR vs QTL panel (`run_mr_against_qtl_data.R`)
5. Volcano plot of MR results
6. Coloc on significant MR results (GWAS LBF vs QTL LBF from `convert_z_to_lbf`)
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
| `mr/coloc_<prefix>_<qtl>.tsv` | Coloc for significant MR hits. |
| `mr/result_<qtl>_<prefix>.html` | HTML report. |
| `plots/volcano_plot_<prefix>_<qtl>.png` | Volcano plot of MR results. |

## Notes

- E2E tests for this pipeline run only when `QTL_DATA_DIR` is set (see [CONTRIBUTING.md](../CONTRIBUTING.md)).
- Fine-mapping uses the same single-ancestry SuSiE path as [`finemap.smk`](finemap.smk) and [`coloc.smk`](coloc.smk).
