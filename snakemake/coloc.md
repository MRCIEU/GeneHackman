# coloc

Standardise, clump, and fine-map **≥2** GWAS (single ancestry only), then run **pairwise**  colocalisation (`coloc::coloc.bf_bf`) on overlapping finemapped signals and generate locus zoom plots.

**Workflow file:** [`coloc.smk`](coloc.smk)  
**Example input:** [`input_templates/coloc.yaml`](input_templates/coloc.yaml)

## Run

```bash
./run_pipeline.sh snakemake/coloc.smk path/to/input.yaml
```

## Input

### GWAS count

- **At least two** GWAS under `gwases:`.

### Ancestry rules

- **All GWAS must share the same `ancestry`.** Multi-ancestry coloc is **not** supported (use [`finemap.smk`](finemap.smk) with MultiSuSiE for cross-population fine-mapping).

### Per-GWAS fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | Yes | Input summary-statistics file. |
| `ancestry` | **Yes** | Must be identical across all GWAS. |
| `N` | Recommended | Sample size for fine-mapping. |
| `columns` | No | Column map or preset. |
| `populate_rsid` / `populate_eaf` | No | See [PIPELINES.md](../PIPELINES.md#shared-yaml-schema). |

### Root-level fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | PLINK clump settings. |
| `finemap.*` | No | Same keys as [finemap.md](finemap.md#root-level-fields). |
| `coloc.overlap_kb` | No | Two signals “overlap” if lead positions are within this distance (kb). Default `1000`. |
| `coloc.p1` | No | Prior: SNP associated with trait 1. Default `1e-4`. |
| `coloc.p2` | No | Prior: SNP associated with trait 2. Default `1e-4`. |
| `coloc.p12` | No | Prior: SNP associated with both traits. Default `5e-6`. |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS
2. PLINK clump per GWAS
3. SuSiE RSS fine-mapping per GWAS (`run_finemap.R`)
4. Pairwise  coloc on overlapping loci (`run_coloc.R`)
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
| `finemap/<prefix>/finemap_complete.txt` | Fine-mapping completion marker per GWAS. |
| `finemap/<prefix>/<CHR>_<BP>_finemap.tsv.gz` | Per-locus SuSiE output with LBF columns. |
| `coloc/coloc_results.tsv` | Pairwise coloc results (sorted by PP.H4 in the HTML report). |
| `coloc/result_coloc.html` | Summary HTML report. |
| `coloc/locus_zoom/` | Locus zoom plot outputs. |
| `coloc/locus_zoom/locus_zoom_complete.txt` | Snakemake completion marker for locus zoom step. |

## Notes

- Coloc requires LBF columns in finemap outputs; multi-ancestry MultiSuSiE outputs are not compatible with this pipeline.
- See [coloc documentation](https://chr1swallace.github.io/coloc/) for interpretation of posterior probabilities.
