# finemap

Standardise, clump, and fine-map GWAS loci with **SuSiE** (single ancestry) or **MultiSuSiE** (multi-ancestry). No trait-trait colocalisation or MR.

**Workflow file:** [`finemap.smk`](finemap.smk)

## Run

```bash
./run_pipeline.sh snakemake/finemap.smk path/to/input.yaml
```

## Input

**GWAS count:** one or more under `gwases:`.

**Ancestry rules:**

| Mode | Condition |
| ---- | --------- |
| **Single-ancestry SuSiE** | All GWAS share the same `ancestry`, **or** only one GWAS. |
| **Multi-ancestry MultiSuSiE** | Each GWAS has a **distinct** `ancestry` (no duplicates). Requires **≥2** ancestries. |

Mixed duplicates (e.g. two EUR + one EAS) **fail at startup**.

Copy from [`input_templates/finemap.yaml`](input_templates/finemap.yaml). Multi-ancestry example: [`tests/testthat/data/snakemake_inputs/finemap_multi_ancestry.yaml`](../tests/testthat/data/snakemake_inputs/finemap_multi_ancestry.yaml).

### Example YAML

```yaml
gwases:
  - file: /path/to/gwas1.tsv.gz
    columns: {}
    N: 100000
    ancestry: EUR
ld_reference: pan_ukb
plink_clump_arguments: "--clump-p1 5e-8 --clump-r2 0.01 --clump-kb 1000"
finemap:
  window_kb: 1000
  max_causal: 10
  coverage: 0.95
  min_abs_corr: 0.5
populate_rsid: false
```

Add more entries under `gwases:` for additional traits or ancestries.

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Input summary-statistics file. |
| `ancestry` | **Yes** | LD reference population for clumping and fine-mapping: `AFR`, `AMR`, `EAS`, `EUR`, `SAS` with `ld_reference: 1000genomes`; also `CSA` and `MID` with `ld_reference: pan_ukb` (`SAS` maps to `CSA`). |
| `ld_reference` | No | This GWAS's LD source for every step, overriding the step and top-level settings: `pan_ukb`, `1000genomes`, or the PLINK prefix (`.bed/.bim/.fam`, output build) of its own genotypes, e.g. those it was computed from. |
| `N` | Recommended | Sample size, used where the file has no per-variant `N` (a per-variant `N` column takes precedence). Fine-mapping uses the median per-variant `N_EFF` (binary traits) or `N` of each region. |
| `study_type` | No | `continuous` (default) or `binary` (case-control). |
| `N_case`, `N_control` | No | Binary traits: numbers of cases and controls, for the effective sample size N_eff = 4 / (1/N_case + 1/N_control) when the file has no `N_EFF` or `N_CASES`/`N_CONTROLS` columns (a file with `N_CASES` and a per-variant total `N` gets controls = `N` − `N_CASES`). For a meta-analysis, give a per-variant `N_EFF` column (the sum of the cohorts' effective sample sizes) if you have one. |
| `label` | No | Name of the GWAS in plots and reports (default: the file name). |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`; lifted to `output.build`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump`. [PLINK clump options](https://zzz.bwh.harvard.edu/plink/clump.shtml). If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus). |

### `finemap`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `finemap.window_kb` | No | Half-width of fine-mapping window around each lead (kb). Default **1000**. |
| `finemap.max_causal` | No | SuSiE `L`, the maximum number of signals per region. Default `10`. It is never changed by the pipeline; see the LD-reference warning below. |
| `finemap.max_n_ratio` | No | Warn when the GWAS sample size is more than this many times the LD reference's. Default `10`. |
| `finemap.region_strategy` | No | `merge` (default), `none` or `tile`; see [Fine-mapping regions](#fine-mapping-regions). The same rules are used for SuSiE and MultiSuSiE. |
| `finemap.remove_ld_outliers` | No | Remove variants flagged by the LD-mismatch check and refit. Default `false` (flagged variants are only reported). |
| `finemap.coverage` | No | Credible set coverage. Default `0.95`. |
| `finemap.min_abs_corr` | No | Minimum \|r\| for credible-set purity. Default `0.5`. |
| `finemap.max_region_kb` | No | Maximum width (kb) of a merged region, or the tile width with `region_strategy: tile`. Default `3000`. |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. With clumping and `false`, RSID population is **partial** unless overridden per GWAS. |
| `populate_eaf` | No | Default `false`. |
| `output.build` | No | Target build after standardisation. Default `GRCh38` (`GRCh37` needs `download_reference_data.sh --grch37`). |
| `ld_reference` | **Yes** | Where LD comes from, with no default: `pan_ukb` (Pan-UK Biobank LD, about 420,000 Europeans and 6,000–8,000 or fewer for other ancestries), `1000genomes` (about 500 people per super-population) or the PLINK prefix of your own genotypes. Used by every LD step here (strand check and EAF filling, clumping and fine-mapping); `harmonise: {ld_reference: ...}`, `clumping:`, `finemap:` or a GWAS's own `ld_reference` sets a different source for one step or one GWAS. See [PIPELINES.md](../PIPELINES.md#common-conventions). |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS
2. PLINK clump per GWAS
3. Fine-map:
   - **Single ancestry:** `run_finemap.R` → SuSiE RSS per GWAS
   - **Multi ancestry:** `run_multisusie.py` → MultiSuSiE across all GWAS jointly

### LD reference and its limits

LD for fine-mapping comes from the source chosen with `ld_reference` (required, no default; `finemap: {ld_reference: ...}` or a GWAS's own `ld_reference` overrides it): the Pan-UK Biobank LD (`pan_ukb`; Karczewski et al. 2025), the 1000 Genomes panel of the GWAS's ancestry (`1000genomes`, about 500 people per super-population), or your own genotypes. Reference-panel LD can miscalibrate fine-mapping, increasingly so as the GWAS gets larger than the panel: credible sets lose coverage and spurious secondary signals appear (Benner et al. 2017; Zou et al. 2022; Kanai et al. 2022). The pipeline therefore:

- warns, in the log, in `finemap_regions.tsv` (`LD_SIZE_WARNING`) and in the HTML report, when the GWAS sample size is more than `finemap.max_n_ratio` (default 10) times that of the LD reference. Nothing is changed automatically: the number of signals is `finemap.max_causal` (default 10). With a much larger GWAS than reference, consider in-sample LD (`ld_reference` on the GWAS), a larger reference, or `max_causal: 1`, for which the PIPs don't depend on the LD matrix (Weissbrod et al. 2020);
- aligns LD to the GWAS effect alleles: variants are matched to the panel on position and alleles (Pan-UKB, user LD) or rsID plus alleles (1000 Genomes), r is flipped where the panel counts the other allele, and variants whose alleles don't match are removed (`N_ALLELE_MISMATCH`);
- checks each region for inconsistency between z-scores and LD with SuSiE-RSS's diagnostics (`susieR::estimate_s_rss` and `kriging_rss`; Zou et al. 2022). `LD_MISMATCH_S` is the estimated inconsistency (0 = consistent; above about 0.1–0.2 suggests mismatch) and `N_LD_OUTLIERS` the variants whose z-score is unlikely given the others (`logLR > 2` and `|z| > 2`), flagged per variant in `LD_OUTLIER`. A flipped allele is the typical outlier.

Each credible set gets an `LD_FLAG`: `OUTLIER_IN_CS` (an LD outlier is in the set: treat the set as suspect), `LD_MISMATCH` (the region's `LD_MISMATCH_S` is above 0.1: credible sets beyond the lead one may be spurious), `OUTLIERS_IN_REGION`, `MHC`, or `OK`. By default nothing is removed; `finemap.remove_ld_outliers: true` removes the flagged variants and fits once more. The same checks run in multi-ancestry fine-mapping, against each ancestry's LD reference (a Python port of the susieR functions, tested against susieR's output). The HTML report lists the flagged credible sets first.

Treat multi-signal results from reference LD as exploratory.

### In-sample and user-supplied LD

Give a GWAS its own `ld_reference:`, the PLINK prefix of the genotypes it was computed from, in the output build. Every LD step then uses them for that GWAS: the strand check and EAF filling, clumping and fine-mapping (SuSiE and MultiSuSiE), as well as the regional plots in the colocalization pipeline. To use them for fine-mapping only, set them as `finemap: {ld_reference: ...}` instead (with `{ancestry}` in the path when GWAS of several ancestries each have their own genotypes). In multi-ancestry fine-mapping each ancestry can have a different source; when any of them is the Pan-UKB LD or your own genotypes, variants are matched on position and alleles, including in the 1000 Genomes panel. Inside a biobank's trusted research environment (where individual-level data can't leave), run the pipeline there from the container with the genotypes mounted; only summary-level outputs are written. To convert a `pgen` or VCF: `plink2 --pfile <prefix> --make-bed --out <prefix>`. Variants whose `.bim` ID is duplicated (or `.`) are left out; to keep them, set unique IDs, e.g. `plink2 --bfile <prefix> --set-all-var-ids @:#:\$r:\$a --make-bed --out <new prefix>`. The genotypes must be in a folder the container can see (your home folder, `PIPELINE_DATA_DIR`, or the project's `data` or `results` folders).

### Fine-mapping regions

Each clumped lead SNP gets a window of `window_kb` either side. With `region_strategy`:

- **`merge`** (default): windows that overlap on the same chromosome are merged into one region, so a signal is fine-mapped once rather than once per nearby lead SNP. A merged region never grows beyond `max_region_kb`: if adding a lead would exceed it, the region is closed and the next region starts immediately after it, so regions never overlap and every lead belongs to exactly one region.
- **`none`**: one window per lead, not merged. Windows of nearby leads overlap, so a signal can be reported more than once.
- **`tile`**: overlapping windows are grouped as for `merge` but without a size limit, and each group wider than `max_region_kb` is covered by tiles `max_region_kb` wide, one every `max_region_kb / 3` (3 Mb tiles every 1 Mb by default, as in Weissbrod et al. 2020). Tiles without a lead are dropped. Each variant's result should be taken from the tile whose centre is nearest: `TILE_CORE` is TRUE there, and `LEAD_IN_TILE_CORE` marks the credible sets to keep.

Each region is named `<CHR>_<BP>` after its most significant lead SNP. Regions overlapping the MHC (chr6:28.5–33.5 Mb) are flagged (`MHC`): LD there is long-range and population-specific, and fine-mapping results are unreliable. The regions, their lead SNPs and the outcome of fine-mapping each one are listed in `finemap_regions.tsv`, with the strategy and the `L` used.

## Output

Snakemake `rule all` tracks completion markers only; additional per-locus files are written alongside them.

Note that the pipeline does not filter out any credible sets based on criteria.  It merely runs finemapping, and lets the user decide
what should be kept.

### Data (`data/`)

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | Standardised GWAS. |
| `clumped_snps/<prefix>.clumped` | Clumped leads per GWAS. |

### Results — single-ancestry SuSiE (`results/finemap/<prefix>/`)

| Path | Description |
| ---- | ----------- |
| `finemap/<prefix>/finemap_complete_<prefix>.txt` | **Snakemake target.** One line = number of fine-mapping regions processed. |
| `finemap/<prefix>/<CHR>_<BP>_finemap.tsv.gz` | Per-region fine-mapping table: every GWAS variant in the region, with SuSiE `Z`, posterior inclusion probability `PIP`, credible-set membership `CS` and per-credible-set log Bayes factors `LBF_1`, `LBF_2`, … (empty for variants not in the LD panel). |
| `finemap/<prefix>/finemap_regions.tsv` | One row per region: bounds, lead SNP(s), `STATUS` (fine-mapped, or why it was skipped), number of variants with LD, number of credible sets, SuSiE convergence, `LD_REFERENCE`, `N_REF` (people the LD comes from), `N_GWAS` (sample size used), `L`, `LD_SIZE_WARNING` (set when the GWAS is more than `max_n_ratio` times larger than the LD reference), `LD_MISMATCH_S`, `N_LD_OUTLIERS`, `N_ALLELE_MISMATCH`, `MHC`, `REGION_STRATEGY`, `FRAC_BEYOND_LD_RADIUS` (share of variant pairs further apart than the pack's 1 Mb radius, set to r = 0) and `LD_MIN_EIGENVALUE` (smallest eigenvalue of the pack LD matrix; below 0 means its negative eigenvalues were set to 0 to make it positive semi-definite). |
| `finemap/<prefix>/finemap_credible_sets.tsv` | One row per credible set: size, coverage, purity (min/mean/median \|r\|), lead variant and its PIP, log Bayes factor, `LD_MISMATCH_S`, `N_LD_OUTLIERS_IN_CS` and `LD_FLAG`. |
| `finemap/result_finemap_<gwas_prefixes>.html` | **HTML report**: settings, a summary per GWAS (L used, credible sets, flagged credible sets, regions with LD mismatch or outliers, MHC regions), credible sets with flagged ones first, and the region table. |

### Results — multi-ancestry MultiSuSiE (`results/finemap/multi_ancestry/`)

Full documentation for interpreting MultiSuSiE results can be [found here](https://deepwiki.com/jordanero/MultiSuSiE/4.3-interpreting-results)

| Path | Description |
| ---- | ----------- |
| `finemap/multi_ancestry/finemap_complete_<gwas_prefixes>.txt` | **Snakemake target.** One line = number of successfully fine-mapped loci. `<gwas_prefixes>` is sorted input GWAS file prefixes joined with `_`. |
| `finemap/multi_ancestry/multisusie_regions.tsv` | One row per region: bounds, status, number of credible sets, `LD_MISMATCH_S_<ancestry>` and `N_LD_OUTLIERS_<ancestry>`, `MHC`, `LD_FLAG`. Credible-set files also get `N_LD_OUTLIERS_IN_CS`, `MAX_LD_MISMATCH_S` and `LD_FLAG`. |
| `finemap/multi_ancestry/multisusie_settings.tsv` | Per ancestry: GWAS N, LD reference and its N, and the `L`, region strategy and window used. |
| `finemap/multi_ancestry/<CHR>_<center>_locus_credible_sets.tsv` | One row per credible set: size, top variant PIP/LBF, `COVERAGE`, MultiSuSiE `PURITY`, `PASS_FILTER`, and per-ancestry `PURITY_*` diagnostics. |
| `finemap/multi_ancestry/<CHR>_<center>_locus_credible_set_variants.tsv` | One row per variant in any credible set: coordinates, alleles, global PIP, per-ancestry `COEF`/`COEF_SD`, and input `BETA`/`SE`/`P`. |

Locus names use the clump window centre (`<CHR>_<center>`).

## Notes

- For pairwise colocalisation between traits, use [`coloc.smk`](coloc.smk) instead.
- Multi-ancestry mode does **not** produce per-trait LBF columns suitable for `coloc::coloc.bf_bf` on individual GWAS directories.
- LD matrices come from each GWAS's fine-mapping LD source: the ancestry-matched 1000 Genomes panel of the output build (`PIPELINE_DATA_DIR/genomic_data/1000genomes/b38_dbsnp156/` for GRCh38), the Pan-UKB LD pack, or your own genotypes.
