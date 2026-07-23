# compare_gwases

Standardise and clump multiple GWAS, then compare them with heterogeneity metrics, LDSC heritability / genetic correlation, and expected-vs-observed replication statistics.

**Workflow file:** [`compare_gwases.smk`](compare_gwases.smk)

## Run

```bash
./run_pipeline.sh snakemake/compare_gwases.smk path/to/input.yaml
```

## Input

**GWAS count:** two or more under `gwases:`.

Copy from [`input_templates/compare_gwases.yaml`](input_templates/compare_gwases.yaml).

### Example YAML

```yaml
gwases:
  - file: /path/to/gwas1.tsv.gz
    columns:
      SNP: rsid
      CHR: chr
      BP: pos
      EA: a1
      OA: a2
      EAF: freq
      P: p
      BETA: beta
      SE: se
    N: 10000
    ancestry: EUR
  - file: /path/to/gwas2.tsv.gz
    columns: {}
    N: 10000
    ancestry: AFR
  - file: /path/to/gwas3.tsv.gz
    columns: {}
    N: 10000
    ancestry: EAS
plink_clump_arguments: "--clump-p1 0.00000005"
populate_rsid: false
```

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Input summary-statistics file. |
| `N` | **Yes** | Sample size; required for LDSC. |
| `ancestry` | **Yes** | `AFR`, `AMR`, `EAS`, `EUR`, or `SAS` — LD reference for clumping and LDSC grouping. |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`; lifted to `output.build` if set. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. With clumping and `false`, RSID population is **partial** unless overridden. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump` (e.g. `--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000`). |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |
| `output.build` | No | Target build after standardisation. Default `GRCh37`. |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS → `data/gwas/<prefix>_std.tsv.gz`
2. PLINK clump per GWAS (ancestry-matched 1000 Genomes panel)
3. Expected vs observed replication (`compare_observed_vs_expected_gwas.R`)
4. Cross-GWAS heterogeneity scores and plots (`calculate_heterogeneity.R`)
5. LDSC per **unique ancestry** present in the input (`run_ldsc.sh`)
6. HTML report (`result_compare_gwases.html`)

## Output

Paths under **`PROJECT_DIR`** unless noted.

### Data (`data/`)

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | Standardised GWAS per input. |
| `clumped_snps/<prefix>.clumped` | PLINK clump output per GWAS. |

### Results (`results/`)

| Path | Description |
| ---- | ----------- |
| `gwas_comparison/expected_vs_observed_outcomes.tsv` | Outcome-level replication metrics. |
| `gwas_comparison/expected_vs_observed_variants.tsv` | Variant-level replication metrics. |
| `gwas_comparison/heterogeneity_scores.tsv` | Heterogeneity statistics across GWAS pairs. |
| `gwas_comparison/result_compare_gwases.html` | Summary HTML report. |
| `plots/snp_heterogeneity_plot.png` | Heterogeneity overview plot. |
| `plots/snps_with_heterogeneity_forest_plot.png` | Per-SNP forest plot for heterogeneous SNPs. |
| `ldsc/results_<ancestry>.log` | LDSC log per ancestry (also other LDSC artefacts under `ldsc/`). |

## Notes

- LDSC runs once per **distinct** `ancestry` value; all GWAS sharing an ancestry are analysed together.
- References: expected vs observed ([Nature 2016](https://doi.org/10.1038/nature17671)); LDSC ([bulik/ldsc](https://github.com/bulik/ldsc)).
