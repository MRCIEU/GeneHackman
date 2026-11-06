# qtl_mr

Standardise, clump, and fine-map **one** outcome GWAS, run Mendelian randomization against the **cis** QTLs of a QTL dataset, plot MR results, and run colocalisation on significant MR hits.

**Workflow file:** [`qtl_mr.smk`](qtl_mr.smk)

## Run

```bash
./run_pipeline.sh snakemake/qtl_mr.smk path/to/input.yaml
```

Requires the QTL dataset and **`QTL_DATA_DIR`** in `.env`: download it with `./download_reference_data.sh --qtl <dataset>` (see the [README](../README.md#3-download-reference-data)), which prints the `QTL_DATA_DIR` to set.

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
plink_clump_arguments: "--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"
finemap:
  window_kb: 1000
  max_causal: 10
  coverage: 0.95
  min_abs_corr: 0.5
qtl:
  dataset: eqtlgen
  exposures: []
  study_type: continuous
populate_rsid: false
```

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Outcome GWAS summary statistics. |
| `ancestry` | **Yes** | LD reference for clumping and fine-mapping. |
| `N` | **Yes** | GWAS sample size (used for fine-mapping). |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`; lifted to `output.build`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `output`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `output.build` | No | Build of all outputs (GWAS, clumps, fine-mapping, plots, MR and coloc). Must be `GRCh38` (the default), because the QTL data are GRCh38 only. GWASes in GRCh37 are lifted to GRCh38. |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump`. If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus). |

### `finemap`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `finemap.window_kb` | No | Half-width of fine-mapping window around each lead (kb). Default `500`. |
| `finemap.max_causal` | No | Max causal signals per locus. Default `10`. |
| `finemap.coverage` | No | Credible set coverage. Default `0.95`. |
| `finemap.min_abs_corr` | No | Minimum \|r\| for credible-set purity. Default `0.5`. |
| `finemap.max_region_kb` | No | Maximum width (kb) of a region made by merging overlapping lead windows. Default `3000`. See [`finemap.md`](finemap.md). |

### `qtl`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `qtl.dataset` | **Yes** | QTL dataset: `eqtlgen`, `brainmeta` or `ukb_pqtl` (see [QTL datasets](#qtl-datasets)). |
| `qtl.subcategory` | No | Only `cis` is supported (the default). Trans QTLs are not offered: they are too often pleiotropic to be MR instruments. |
| `qtl.exposures` | No | List of exposure names to test; if empty (`[]`), all exposures in the panel are run. |
| `qtl.study_type` | No | Trait type of the **QTL** data, used to convert QTL summary statistics to log Bayes factors in the coloc step: `continuous` (default; expression and protein levels) or `categorical`. |
| `qtl.N` | No | QTL sample size, used in the coloc step only when the QTL files have no `N` column. If neither is available, coloc is skipped for that exposure. |

### `coloc`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `coloc.p1` | No | Prior: SNP associated with the GWAS trait. Default `1e-4`. |
| `coloc.p2` | No | Prior: SNP associated with the QTL trait. Default `1e-4`. |
| `coloc.p12` | No | Prior: SNP associated with both traits. Default `5e-6`. |

### QTL datasets

All three are cis QTLs in European-ancestry samples. They are GRCh38 only (`QTL_DATA_DIR/<dataset>/GRCh38/`), so the pipeline runs with `output.build: GRCh38`. They were built from the standardised (pre-imputation) study data of the genotype-phenotype map 1.0.0 by [`reference_data/grch38/31_qtl_from_gpm.R`](../reference_data/grch38/31_qtl_from_gpm.R).

| `dataset` | Exposures | Tissue | N used | Source |
| --------- | --------- | ------ | - | ------ |
| `eqtlgen` | gene expression (cis-eQTLs) | blood | 27,750 | eQTLGen (Võsa et al., *Nat Genet* 2021) |
| `brainmeta` | gene expression (cis-eQTLs) | brain (meta-analysis) | 2,443 | BrainMeta v2 (Qi et al., *Nat Genet* 2022) |
| `ukb_pqtl` | plasma protein levels (cis-pQTLs, ±1 Mb of the gene) | plasma | 34,557 | UK Biobank Pharma Proteomics Project (Sun et al., *Nature* 2023) |

Instruments are each exposure's cis lead variants with P < 5e-8, one per LD block. An exposure with several instruments gets an IVW estimate, otherwise a Wald ratio. For coloc, the QTL summary statistics of the lead's LD block are used. N is one value per dataset, as recorded by the genotype-phenotype map (eQTLGen's per-variant N, up to 31,684, is not available there).

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
4. MR vs the cis QTLs of the dataset (`run_mr_against_qtl_data.R`)
5. Volcano plot of MR results
6. Coloc on significant MR results (GWAS LBF vs QTL LBF from `convert_z_to_lbf`)
7. HTML report

## Output

Let `<prefix>` = outcome GWAS file prefix, `<qtl>` = `{dataset}_cis`.

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
