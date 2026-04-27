# Pipeline Input

All piplines will require a GWAS object or list of objects, they have the following properties

### GWAS Object:
* `file`: Full path to file.  Mandatory
* `ancestry`: one of `AFR, AMR, EAS, EUR, SAS`.  Mandatory for some pipelines
* `build`: one of `GRCh36, GRCh37, GRCh38`.  Default: `GRCh37`
* `columns`: object of column name maps (or string of predefined map).  Explained below
* `populate_rsid`: boolean value.  Populates RSID column if it doesn't exist. Default `false`
* `populate_eaf`: boolean value.  When `true`, missing `EAF` is filled from the 1000 Genomes LD reference panel (`<THOUSAND_GENOMES_DIR>/<ancestry>.bim` + `.frq`; generate `.frq` with `plink --bfile <prefix> --freq --out <prefix>`). Only variants with `is.na(EAF)` are updated. Default `false`. **Requires `ancestry` on that GWAS** (same field as `compare_gwases`: `AFR`, `AMR`, `EAS`, `EUR`, `SAS`); there is no default ancestry.

Pipeline-level `populate_eaf` (under the root JSON object, next to `populate_rsid`) applies to every GWAS unless a GWAS sets its own `populate_eaf`. Each GWAS that ends up with `populate_eaf` enabled must still include its own `"ancestry"` in the JSON.

**GWAS Columns:**

With each GWAS file, you can specify column names ex. `{"P":"pval", ...}`, if you do not specify header names it will assume your GWAS has default names

Default names : `SNP, CHR, BP, EA, OA, EAF, BETA, SE, OR, P, LOG_P, Z, OR_LB, OR_UB, RSID, N, N_CASES, ENSEMBL_ID, GENE_NAME`

* Mandatory Columns: `CHR, BP, EA, OA`
Effect Column Options.  One of these sets are mandatory:
  * `BETA, SE`
  * `OR, OR_LB, OR_UB`
  * `Z`
* P-value Column: `P or LOG_P`

Alternatively, `columns` accepts a string of some of the more common output formats (ex: `metal`, `gwama`).  There are also a list of [predefined common column maps](../inst/extdata/predefined_column_maps.csv) here.

### standardise_gwas 

All pipelines will standardise each GWAS before running the subsequent steps.  The `SNP` field will be recalculated as `CHR:POS_EA_OA`, where EA and OA are ordered alphabetically, and the subsequent BETA and EAF will be adjusted accordingly. Optional `populate_eaf` (per GWAS or pipeline) runs **after** this harmonisation and joins on `CHR`/`BP` to the reference panel.

* `n GWAS objects`: See above for GWAS Object explanation
* `output`:
  * `build`: one of the supported reference builds.  Default: `GRCh37`
  * `columns`: same format as input columns.  Either a object {} or predefined map string

## disease_progression 

* `2 GWAS objects`: Incident and Subsequent GWAS.  See above for GWAS Object explanation
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)
* `output`: 2 fields in `adjusted_gwas`, where `type` can be `slopehunter`, `cwls`, or `mr_ivw`, and `p_val` can be any of `0.1, 0.01, 0.001, 1e-5`

## compare_gwases 

* `n GWAS objects`: See above for GWAS Object explanation
  * Please also explicitly include `N` in the GWAS object, for use in the LDSC tool
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

## finemap

* `n GWAS objects`: See above for GWAS Object explanation. Each GWAS requires `ancestry`.
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call. [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)
* `finemap` (optional):
  * `window_kb`: half-width of fine-mapping window in kb (default 1000 = ±1 Mb)
  * `max_causal`: maximum number of causal signals per locus (SuSiE L, default 10)
  * `coverage`: credible set coverage (default 0.95)
  * `min_abs_corr`: minimum absolute correlation for credible set purity (default 0.5)

## coloc

Runs standardisation, clumping, and fine-mapping on all input GWASes, then performs BF-BF colocalization (`coloc::coloc.bf_bf`) on any overlapping finemapped signals (within ±`overlap_kb` kb) across all trait pairs.

* `n GWAS objects` (at least 2): See above for GWAS Object explanation. Each GWAS requires `ancestry`.
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call
* `finemap` (optional): same options as the finemap pipeline above
* `coloc` (optional):
  * `overlap_kb`: distance in kb to define overlapping signals (default 1000 = ±1 Mb)
  * `p1`: prior probability a SNP is associated with trait 1 (default 1e-4)
  * `p2`: prior probability a SNP is associated with trait 2 (default 1e-4)
  * `p12`: prior probability a SNP is associated with both traits (default 5e-6)

## qtl_mr

QTL summary statistics are read from **`QTL_DATA_DIR`** (see `.env` / `.env_example`): same tree as under `PIPELINE_DATA_DIR/qtl_datasets` when **`QTL_DATA_DIR`** is unset (`pqtl`, `metabrain`, `eqtlgen`, …). Point **`QTL_DATA_DIR`** at a separate mount or bucket path to download only the QTL data you need.

This pipeline now includes clumping and fine-mapping of the GWAS before MR. Colocalization uses finemapped LBF values from the GWAS and converts QTL Z-scores to LBF via `convert_z_to_lbf`, then runs `coloc::coloc.bf_bf`.

* `1 GWAS object`: See above for GWAS Object explanation. Requires `ancestry`.
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call
* `finemap` (optional): same options as the finemap pipeline above
* `qtl`:
  * `dataset`: see below for options
  * `subcategories`: list of subcategories to run, see below for options
    * If left empty, all subcategories will be run
  * `exposures`: list of specific exposures to be run. If left empty all exposures will be run in the dataset
    * example: `["IL6", "CCL2"]`
  * `study_type`: study type for LBF conversion, either `continuous` (default) or `categorical`

Available Datasets and Subcategories:
* dataset: `metabrain`
  * subcategories: `basalganglia, cerebellum, cortex, hippocampus, spinalcord`
* dataset: `eqtlgen`
  * subcategories: `cis` (future: `trans`)
