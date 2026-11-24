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
ld_reference: pan_ukb
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
| `ancestry` | **Yes** | LD population for clumping, and the LDSC grouping: `AFR`, `AMR`, `EAS`, `EUR`, `SAS` (also `CSA` and `MID` with Pan-UKB; `SAS` uses Pan-UKB's `CSA`). |
| `ld_reference` | No | This GWAS's LD source for every step, overriding the step and top-level settings: `pan_ukb`, `1000genomes`, or the PLINK prefix (`.bed/.bim/.fam`, output build) of its own genotypes, e.g. those it was computed from. |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`; lifted to `output.build`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. With clumping and `false`, RSID population is **partial** unless overridden. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump` (e.g. `--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000`). If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus). |
| `ld_reference` | **Yes** | Where LD comes from, with no default: `pan_ukb` (Pan-UK Biobank LD, about 420,000 Europeans and 6,000–8,000 or fewer for other ancestries), `1000genomes` (about 500 people per super-population) or the PLINK prefix of your own genotypes. Used by every LD step here (strand check and EAF filling, and clumping; not LDSC, which uses `ldsc.ld_scores`); `harmonise: {ld_reference: ...}`, `clumping:` or a GWAS's own `ld_reference` sets a different source for one step or one GWAS. See [PIPELINES.md](../PIPELINES.md#common-conventions). |
| `ldsc.ld_scores` | No | LD scores for LDSC, the exception to `ld_reference`: `pan_ukb` (default, Pan-UKB HapMap3 LD scores) whichever `ld_reference` is, or `1000genomes_b37` only to reproduce release 1.1.0 results (those scores were computed with a 1-kb window instead of 1 cM). |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |
| `output.build` | No | Target build after standardisation. Default `GRCh38` (`GRCh37` needs `download_reference_data.sh --grch37`). |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS → `data/gwas/<prefix>_std.tsv.gz`
2. Clump per GWAS with its LD source (`ld_reference`): the Pan-UKB LD pack, PLINK on the ancestry-matched 1000 Genomes panel, or PLINK on your own genotypes
3. Expected vs observed replication (`compare_observed_vs_expected_gwas.R`)
4. Cross-GWAS heterogeneity scores and plots (`calculate_heterogeneity.R`)
5. LDSC per **unique ancestry**, with the Pan-UKB LD scores (`ldsc.ld_scores`). This is the one step that doesn't follow `ld_reference`: LD-score regression needs precomputed LD scores, and the Pan-UKB scores are the only correct ones in the reference data, so they are used whichever LD reference the other steps use. It runs for each ancestry present in the input (`run_ldsc.sh`): heritability of each GWAS and genetic correlation of **every pair** within the ancestry. LDSC's own `--rg a,b,c` estimates only the first trait against the rest, so it is run once per GWAS. Binary traits use the effective sample size (`N_EFF`) and, with `prevalence`, report liability-scale h² with a sample prevalence of 0.5 (Grotzinger et al. 2023). Genetic correlations are flagged as unreliable when either trait's h² z-score is below 4 or |rg| > 1. LDSC is not run between ancestries: cross-ancestry genetic correlation needs dedicated methods (e.g. Popcorn) and is not part of the pipeline. Don't give LDSC genomic-control-corrected summary statistics.
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
