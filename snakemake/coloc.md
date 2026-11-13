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
ld_reference: pan_ukb
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
| `ld_reference` | No | This GWAS's LD source for every step, overriding the step and top-level settings: `pan_ukb`, `1000genomes`, or the PLINK prefix (`.bed/.bim/.fam`, output build) of its own genotypes, e.g. those it was computed from. |
| `N` | Recommended | Sample size for fine-mapping. |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`; lifted to `output.build`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `output`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `output.build` | No | Build of all outputs (GWAS, clumps, fine-mapping, plots, MR and coloc). Default `GRCh38` (`GRCh37` needs `download_reference_data.sh --grch37`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump`. If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus). |
| `ld_reference` | **Yes** | Where LD comes from, with no default: `pan_ukb` (Pan-UK Biobank LD, about 420,000 Europeans and 6,000–8,000 or fewer for other ancestries), `1000genomes` (about 500 people per super-population) or the PLINK prefix of your own genotypes. Used by every LD step here (strand check and EAF filling, clumping, fine-mapping and the regional plots); `harmonise: {ld_reference: ...}`, `clumping:`, `finemap:`, `locus_zoom:` or a GWAS's own `ld_reference` sets a different source for one step or one GWAS. See [PIPELINES.md](../PIPELINES.md#common-conventions). |

### `finemap`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `finemap.window_kb` | No | Half-width of fine-mapping window around each lead (kb). Default `500`. |
| `finemap.max_causal` | No | Max causal signals per locus. Default `10`. |
| `finemap.coverage` | No | Credible set coverage. Default `0.95`. |
| `finemap.min_abs_corr` | No | Minimum \|r\| for credible-set purity. Default `0.5`. |
| `finemap.max_region_kb` | No | Maximum width (kb) of a region made by merging overlapping lead windows. Default `3000`. See [`finemap.md`](finemap.md). |

### `coloc`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `coloc.overlap_kb` | No | Two fine-mapped loci are tested if their lead positions are within this distance (kb), or if their fine-mapping regions overlap. Default `1000`. |
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
| `finemap/<prefix>/<CHR>_<BP>_finemap.tsv.gz` | Per-region SuSiE output with PIP and LBF columns. |
| `finemap/<prefix>/finemap_regions.tsv`, `finemap_credible_sets.tsv` | Region and credible-set summaries (see [`finemap.md`](finemap.md)). |
| `coloc/coloc_results.tsv` | Pairwise coloc results (sorted by PP.H4 in the HTML report). |
| `coloc/result_coloc.html` | Summary HTML report. |
| `coloc/locus_zoom/` | Locus zoom plot outputs. |
| `coloc/locus_zoom/locus_zoom_complete_<gwas_prefixes>.txt` | Snakemake completion marker for locus zoom step. |

## Notes

- Coloc requires LBF columns in finemap outputs; multi-ancestry MultiSuSiE outputs are not compatible with this pipeline.
- Only variants with log Bayes factors in both loci (i.e. present in the LD panel for both fine-mapping runs) are used.
- See [coloc documentation](https://chr1swallace.github.io/coloc/) for interpretation of posterior probabilities.
- **Sample overlap.** coloc assumes the two traits are measured in independent samples (Giambartolomei et al. 2014). Traits from the same or overlapping cohorts (e.g. nested case definitions with shared controls) violate this; same-biobank colocalisations are common, but read the posterior probabilities with caution, and note that a spurious signal from shared samples (or from the shared LD reference) can appear in both traits. The LDSC cross-trait intercept from `compare_gwases.smk` (`GCOV_INTERCEPT`, near 0 without overlap) shows how much overlap there is.
- `label:` per GWAS names the traits in the results and plots.
- **LD flags.** `coloc_results.tsv` has `ld_flag1`/`ld_flag2`, the fine-mapping `LD_FLAG` of each trait's region, and `ld_check` (TRUE when either isn't `OK`). A colocalization from a region with LD mismatch may be spurious; the report counts them and lists the flagged credible sets.
