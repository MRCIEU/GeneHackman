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
| `N` | **Yes** | Sample size, for LDSC when the file has no per-variant `N_EFF` or `N` column (those take precedence). |
| `study_type` | No | `continuous` (default) or `binary`. |
| `N_case`, `N_control` | No | Binary traits: case and control counts, for the effective sample size. |
| `prevalence` | No | Binary traits: population prevalence, for liability-scale heritability. |
| `label` | No | Name of the GWAS in plots and tables (default: the file name). |
| `ancestry` | **Yes** | `AFR`, `AMR`, `EAS`, `EUR`, or `SAS` — LD reference for clumping and LDSC grouping. These 1000 Genomes codes are needed even with `ld_reference: pan_ukb`, because LDSC uses the 1000 Genomes LD scores (`SAS` clumps against Pan-UKB `CSA`). |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`; lifted to `output.build`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. With clumping and `false`, RSID population is **partial** unless overridden. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump` (e.g. `--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000`). If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus). |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |
| `output.build` | No | Target build after standardisation. Default `GRCh38` (`GRCh37` needs `download_reference_data.sh --grch37`). |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS → `data/gwas/<prefix>_std.tsv.gz`
2. Clump per GWAS: PLINK on the ancestry-matched 1000 Genomes panel, or the Pan-UKB LD pack with `ld_reference: pan_ukb`
3. Expected vs observed replication (`compare_observed_vs_expected_gwas.R`)
4. Cross-GWAS heterogeneity scores and plots (`calculate_heterogeneity.R`)
5. LDSC per **unique ancestry** present in the input (`run_ldsc.sh`): heritability of each GWAS and genetic correlation of **every pair** within the ancestry. LDSC's own `--rg a,b,c` estimates only the first trait against the rest, so it is run once per GWAS. Binary traits use the effective sample size (`N_EFF`) and, with `prevalence`, report liability-scale h² with a sample prevalence of 0.5 (Grotzinger et al. 2023). Genetic correlations are flagged as unreliable when either trait's h² z-score is below 4 or |rg| > 1. LDSC is not run between ancestries: cross-ancestry genetic correlation needs dedicated methods (e.g. Popcorn) and is not part of the pipeline. Don't give LDSC genomic-control-corrected summary statistics.
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
| `ldsc/results_<ancestry>_h2.tsv` | Heritability per GWAS: observed and (binary traits with prevalence) liability-scale h², h² z-score, intercept, ratio, mean χ², and flags. |
| `ldsc/results_<ancestry>_rg.tsv` | Genetic correlation for every pair: rg, SE, P, cross-trait intercept (an indication of sample overlap), both traits' h² z-scores and `RG_RELIABLE`. |
| `ldsc/results_<ancestry>.log` | All LDSC logs for the ancestry, concatenated (the individual logs are beside it). |

## Notes

- LDSC runs once per **distinct** `ancestry` value; all GWAS sharing an ancestry are analysed together.
- References: expected vs observed ([Nature 2016](https://doi.org/10.1038/nature17671)); LDSC ([bulik/ldsc](https://github.com/bulik/ldsc)).
