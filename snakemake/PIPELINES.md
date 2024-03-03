# Pipeline Input

All piplines will require a GWAS or list of GWASes, they have the following keys

GWAS Object:
* `file`: Full path to file.  Mandatory
* `ancestry`: one of `AFR, AMR, EAS, EUR, SAS`.  Mandatory for some pipelines
* `build`: one of `GRCh36, GRCh37, GRCh38`.  Default: `GRCh37`
* `columns`: object of column name maps (or string of predefined map).  Explained below
* `populate_rsid`: boolean value.  Populates RSID column if it doesn't exist. Default `false`

### GWAS Columns:

With each GWAS file, you can specify column names ex. `{"P":"pval", ...}`, if you do not specify header names it will assume your GWAS has default names

* Mandatory Columns: `CHR, BP, EA, OA`
Effect Column Options.  One of these sets are mandatory:
  * `BETA, SE`
  * `OR, OR_LB, OR_UB`
  * `Z`
* P-value Column: `P or LOG_P`

Alternatively, `columns` accepts a string of some of the more common output formats (ex: `metal`, `gwama`).  There are also a list of [predefined common column maps](../inst/extdata/predefined_column_maps.csv) here.

### standardise_gwas 

All pipelines will standardise each GWAS before running the subsequent steps.  The `SNP` field will be recalculated as `CHR:POS_EA_OA`, where EA and OA are ordered alphabetically, and the subsequent BETA and EAF will be adjusted accordingly

* `n GWAS objects`: GWAS summary statistics files and optional column name map
* `output`: 

## collider_bias 

### input.json

* `2 GWAS objects`: Incident and Subsequent GWAS summary statistics file and optional column name map
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

## compare_gwases 

### input.json
* `n GWAS objects`: GWAS summary statistics files and optional column name map
* `plink_clump_arguments`: arguments that are fed into the `plink --clump` call.  [Options here](https://zzz.bwh.harvard.edu/plink/clump.shtml)

## MR for QTLs 

### input.json
* `1 GWAS object`: GWAS summary statistics file, ancestry, and optional column name map
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
  * subcategories: `cis, trans`
