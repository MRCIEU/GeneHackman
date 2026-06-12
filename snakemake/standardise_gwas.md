# standardise_gwas

Harmonise one or more GWAS summary-statistic files into the GeneHackman standard format. No PLINK clumping or downstream analysis.

**Workflow file:** [`standardise_gwas.smk`](standardise_gwas.smk)  
**Example input:** [`input_templates/standardise_gwases.yaml`](input_templates/standardise_gwases.yaml)

## Run

```bash
./run_pipeline.sh snakemake/standardise_gwas.smk path/to/input.yaml
```

## Input

### GWAS count

- **One or more** GWAS files under `gwases:`.

### Per-GWAS fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | Yes | Path to VCF, CSV, TSV, or TXT (`.gz` / `.zip` supported). |
| `columns` | No | Column rename map or preset name (`metal`, `gwama`, `opengwas_phewas`, …). See [`inst/extdata/predefined_column_maps.csv`](../inst/extdata/predefined_column_maps.csv). |
| `N` | No | Sample size; written when the file has no `N` column. |
| `build` | No | Input reference build: `GRCh36`, `GRCh37`, `GRCh38`. Default `GRCh37`. |
| `ancestry` | If `populate_eaf: true` | One of `AFR`, `AMR`, `EAS`, `EUR`, `SAS`. |
| `populate_rsid` | No | Per-GWAS override of root `populate_rsid`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes `.frq` for `ancestry`. |
| `flip_alleles` | No | Default inherits root (`true`). **`false` is allowed only in this pipeline.** |
| `remove_extra_columns` | No | Default `false`. |

### Root-level fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `output.build` | No | Target reference build after liftover. Default `GRCh37`. |
| `output.columns` | No | Rename standard columns on output. |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |
| `flip_alleles` | No | Default `true`. |
| `is_test` | No | Used by tests / logging. |

### Source data requirements

After column mapping, each row needs:

- **Required:** `CHR`, `BP`, `EA`, `OA`
- **Effect:** `BETA`+`SE`, or `OR`+`OR_LB`+`OR_UB`, or `Z`
- **P-value:** `P` or `LOG_P`

Standardisation converts OR→BETA, builds harmonised SNP id `CHR:BP_EA_OA` (alphabetical EA/OA when `flip_alleles: true`), optional liftover, optional RSID / gene annotation / EAF fill.

## Output

All paths are under **`PROJECT_DIR/data/`** (see `.env`).

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | One harmonised file per input GWAS (`<prefix>` = basename of input file). |

Intermediate VCF inputs are converted to TSV under `data/gwas/` before standardisation.

## Notes

- This is the **only** pipeline that accepts `flip_alleles: false`.
- Use `output.build` for liftover (e.g. input `build: GRCh38`, `output.build: GRCh37`).
- Downstream pipelines (`finemap`, `coloc`, …) include the same standardisation step automatically; run this pipeline alone when you only need harmonised files.
