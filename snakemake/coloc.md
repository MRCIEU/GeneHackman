# coloc

Standardise, clump, and fine-map **≥2** GWAS (single ancestry only), then run **pairwise** colocalisation (`coloc::coloc.bf_bf`) on overlapping finemapped signals and generate locus zoom plots.

**Workflow file:** [`coloc.smk`](coloc.smk)

## Run

```bash
./run_pipeline.sh snakemake/coloc.smk path/to/input.yaml
```

## Input

**GWAS count:** at least two under `gwases:`.

**Ancestry:** all GWAS must share the same `ancestry`. Multi-ancestry coloc is **not** supported (use [`finemap.smk`](finemap.smk) with MultiSuSiE for cross-population fine-mapping).

Copy from [`input_templates/coloc.yaml`](input_templates/coloc.yaml).

### Example YAML

```yaml
gwases:
  - file: /path/to/trait1.tsv.gz
    columns: {}
    N: 100000
    ancestry: EUR
  - file: /path/to/trait2.tsv.gz
    columns: {}
    N: 80000
    ancestry: EUR
plink_clump_arguments: "--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"
finemap:
  window_kb: 1000
  max_causal: 10
  coverage: 0.95
  min_abs_corr: 0.5
coloc:
  overlap_kb: 1000
  p1: 1e-4
  p2: 1e-4
  p12: 5e-6
populate_rsid: false
```

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Input summary-statistics file. |
| `ancestry` | **Yes** | Must be **identical** across all GWAS. |
| `N` | Recommended | Sample size for fine-mapping. |
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

### `coloc`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `coloc.overlap_kb` | No | Two signals “overlap” if lead positions are within this distance (kb). Default `1000`. |
| `coloc.p1` | No | Prior: SNP associated with trait 1. Default `1e-4`. |
| `coloc.p2` | No | Prior: SNP associated with trait 2. Default `1e-4`. |
| `coloc.p12` | No | Prior: SNP associated with both traits. Default `5e-6`. |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS
2. PLINK clump per GWAS
3. SuSiE RSS fine-mapping per GWAS (`run_finemap.R`)
4. Pairwise coloc on overlapping loci (`run_coloc.R`)
5. Locus zoom plots for coloc hits (`run_locus_zoom.R`)
6. HTML report

## Output

### Data (`data/`)

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | Standardised GWAS. |
| `clumped_snps/<prefix>.clumped` | Clumped leads per GWAS. |

### Results (`results/`)

| Path | Description |
| ---- | ----------- |
| `finemap/<prefix>/finemap_complete_<prefix>.txt` | Fine-mapping completion marker per GWAS. |
| `finemap/<prefix>/<CHR>_<BP>_finemap.tsv.gz` | Per-locus SuSiE output with LBF columns. |
| `coloc/coloc_results.tsv` | Pairwise coloc results (sorted by PP.H4 in the HTML report). |
| `coloc/result_coloc.html` | Summary HTML report. |
| `coloc/locus_zoom/` | Locus zoom plot outputs. |
| `coloc/locus_zoom/locus_zoom_complete_<gwas_prefixes>.txt` | Snakemake completion marker for locus zoom step. |

## Notes

- Coloc requires LBF columns in finemap outputs; multi-ancestry MultiSuSiE outputs are not compatible with this pipeline.
- See [coloc documentation](https://chr1swallace.github.io/coloc/) for interpretation of posterior probabilities.
