# Pipeline Input

All piplines will require a GWAS object or list of objects, they have the following properties

### GWAS Object:
* `file`: Full path to file.  Mandatory
* `ancestry`: one of `AFR, AMR, EAS, EUR, SAS`.  Mandatory for some pipelines
* `build`: one of `GRCh36, GRCh37, GRCh38`.  Default: `GRCh37`
* `columns`: object of column name maps (or string of predefined map).  Explained below
* `populate_rsid`: boolean value.  Populates RSID column if it doesn't exist. Default `false`

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

All pipelines will standardise each GWAS before running the subsequent steps.  The `SNP` field will be recalculated as `CHR:POS_EA_OA`, where EA and OA are ordered alphabetically, and the subsequent BETA and EAF will be adjusted accordingly

* `n GWAS objects`: See above for GWAS Object explanation
* `output`:
  * `build`: one of the supported reference builds.  Default: `GRCh37`
  * `columns`: same format as input columns.  Either a object {} or predefined map string

## collider_bias 

* `2 GWAS objects`: Incident and Subsequent GWAS.  See above for GWAS Object explanation
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

## compare_gwases 

* `n GWAS objects`: See above for GWAS Object explanation
  * Please also explicitly include `N` in the GWAS object, for use in the LDSC tool
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

## qtl_mr

* `1 GWAS object`: See above for GWAS Object explanation
* `qtl`:
  * `dataset`: see below for options
  * `subcategories`: list of subcategories to run, see below for options
    * If left empty, all subcategories will be run
  * `exposures`: list of specific exposures to be run. If left empty all exposures will be run in the dataset
    * example: `["IL6", "CCL2"]`

Available Datasets and Subcategories:
* dataset: `metabrain`
  * subcategories: `basalganglia, cerebellum, cortex, hippocampus, spinalcord`
* dataset: `eqtlgen`
  * subcategories: `cis` (future: `trans`)
