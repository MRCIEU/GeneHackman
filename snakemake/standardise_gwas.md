# standardise_gwas

Harmonise one or more GWAS summary-statistic files into the GeneHackman standard format. No PLINK clumping or downstream analysis.

**Workflow file:** [`standardise_gwas.smk`](standardise_gwas.smk)

## Run

```bash
./run_pipeline.sh snakemake/standardise_gwas.smk path/to/input.yaml
```

## Input

**GWAS count:** one or more under `gwases:`.

Copy from [`input_templates/standardise_gwases.yaml`](input_templates/standardise_gwases.yaml) and edit paths.

### Example YAML

```yaml
gwases:
  - file: /path/to/gwas.tsv.gz
    columns:
      CHR: CHR
      BP: BP
      EA: EA
      OA: OA
      EAF: EAF
      P: P
      BETA: BETA
      SE: SE
    build: GRCh38
    ancestry: EUR
    populate_eaf: false
output:
  build: GRCh37
populate_rsid: false
```

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Path to VCF, CSV, TSV, or TXT (`.gz` / `.zip` supported). |
| `columns` | No | Maps logical names to file headers, or a preset (`metal`, `gwama`, `opengwas_phewas`, …). See [`inst/extdata/predefined_column_maps.csv`](../inst/extdata/predefined_column_maps.csv). |
| `build` | No | Input reference build: `GRCh36`, `GRCh37`, `GRCh38`. Default `GRCh37`. |
| `ancestry` | If `populate_eaf: true` | One of `AFR`, `AMR`, `EAS`, `EUR`, `SAS`. |
| `populate_rsid` | No | Per-GWAS override of root `populate_rsid`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes `.frq` for `ancestry`. |
| `flip_alleles` | No | Inherits root default `true`. **`false` is allowed only in this pipeline.** |
| `N` | No | Sample size; written when the file has no `N` column. |
| `remove_extra_columns` | No | Default `false`. |

### `output`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `output.build` | No | Target reference build after liftover. Default `GRCh37`. |
| `output.columns` | No | Rename standard columns on output. |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |
| `flip_alleles` | No | Default `true`. |
| `is_test` | No | Used by tests / logging. |

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
- Downstream pipelines include the same standardisation step automatically; run this pipeline alone when you only need harmonised files.
