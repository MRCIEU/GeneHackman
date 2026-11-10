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
  build: GRCh38
populate_rsid: false
```

### `gwases[]`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `file` | **Yes** | Path to VCF, CSV, TSV, or TXT (`.gz` / `.zip` supported). |
| `columns` | No | Maps logical names to file headers, or a preset (`metal`, `gwama`, `opengwas`, `opengwas_phewas`, `gwas_catalog` for the GWAS Catalog's GWAS-SSF format, …). See [`inst/extdata/predefined_column_maps.csv`](../inst/extdata/predefined_column_maps.csv). **Column names are not detected automatically**: without `columns`, the file must already use GeneHackman's names (`CHR`, `BP`, `EA`, `OA`, `BETA`, `SE`, `EAF`, `P`, `N`, …). |
| `study_type` | No | `continuous` (default) or `binary`. Binary traits get an effective sample size column `N_EFF`: kept from the file if present, else from `N_CASES`/`N_CONTROLS` columns (or `N_CASES` and a per-variant total `N`, with controls = `N` − `N_CASES`), else from `N_case`/`N_control` below. |
| `N_case`, `N_control` | No | Binary traits: numbers of cases and controls. |
| `build` | No | Input reference build: `GRCh36`, `GRCh37`, `GRCh38`. Default `GRCh37`, because most published GWAS still are; standardisation stops if the positions match the other build much better (checked against the Pan-UKB LD pack). An already standardised `*_std.tsv.gz` input must declare `build` equal to `output.build`. |
| `ancestry` | If `populate_eaf: true` | One of `AFR`, `AMR`, `EAS`, `EUR`, `SAS`; with `ld_reference: pan_ukb` EAF comes from the Pan-UKB LD pack (matched on SNP keys) and `CSA`, `MID` are also allowed. |
| `populate_rsid` | No | Per-GWAS override of root `populate_rsid`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes `.frq` for `ancestry`. |
| `flip_alleles` | No | Inherits root default `true`. **`false` is allowed only in this pipeline.** |
| `N` | No | Sample size; written when the file has no `N` column. |
| `remove_extra_columns` | No | Default `false`. |

### `output`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `output.build` | No | Target reference build after liftover. Default `GRCh38`; `GRCh37` needs `download_reference_data.sh --grch37` (see README). |
| `output.columns` | No | Rename standard columns on output. |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. |
| `populate_eaf` | No | Default `false`. |
| `flip_alleles` | No | Default `true`. |
| `is_test` | No | Used by tests / logging. |
| `liftover.complement_minus_strand` | No | Default `false`: variants that map to the minus strand are dropped. `true` keeps non-palindromic ones with complemented alleles. |
| `liftover.max_drop_fraction` | No | Default `0.05`: stop if liftover drops more than this share of variants (usually a wrong `build`). |
| `liftover.check_build` | No | Default `true`: check the declared `build` against the Pan-UKB LD pack's GRCh37 and GRCh38 keys. |

Standardisation converts OR→BETA, builds harmonised SNP id `CHR:BP_EA_OA` (alphabetical EA/OA when `flip_alleles: true`), optional liftover, optional RSID / gene annotation / EAF fill.

## Output

All paths are under **`PROJECT_DIR/data/`** (see `.env`).

| Path | Description |
| ---- | ----------- |
| `gwas/<prefix>_std.tsv.gz` | One harmonised file per input GWAS (`<prefix>` = basename of input file). |

Intermediate VCF inputs are converted to TSV under `data/gwas/` before standardisation.

## Harmonisation rules

- **Alleles** are put in alphabetical order (`EA` < `OA`, also for indels and I/D coding), flipping `BETA`, `Z` and `EAF`; the `SNP` key is `CHR:BP_EA_OA`. `flip_alleles: false` restores the input orientation on output.
- **Effects:** `OR` with `OR_LB`/`OR_UB` gives `BETA = log(OR)` and its SE; `Z` alone (no `BETA`/`SE`) gives standardised effects `SE = 1/sqrt(2p(1-p)(N + Z²))`, `BETA = Z × SE` (needs `EAF` and `N`); `LOG_P` gives `P`; `P = 0` is set to the smallest positive double, with a warning.
- **Sample size:** a per-variant `N` in the file takes precedence; the input `N` only fills missing values.
- **Strand** (when `ancestry` is given and its LD reference is installed; `harmonise.strand_check`, default true): alleles are compared with the reference panel (`<output>.allele_qc.tsv` gives the counts). Non-palindromic SNVs that match only on the other strand are complemented; A/T and C/G SNVs are flipped when their EAF is on the other side of 0.5 from the reference's, and are ambiguous when either minor allele frequency is 0.42 or more (kept, or removed with `harmonise.drop_ambiguous_palindromic: true`); variants whose EAF differs from the reference's by more than 0.2 are reported. Without `EAF` in the file, palindromic SNVs can't be checked.
- **Build:** liftover with QC, see the liftover section of [PIPELINES.md](../PIPELINES.md).

The rules are tested against known answers in `tests/testthat/test_harmonisation_known_answers.R`.

## Notes

- This is the **only** pipeline that accepts `flip_alleles: false`.
- Liftover runs whenever `build` differs from `output.build` (by default GRCh37 → GRCh38). Variants whose new position can't be trusted are dropped, each with a reason: not mapped, mapped to an alt/unplaced contig or another chromosome, mapped to the minus strand, not returning to the start when lifted back (conversion-unstable positions), or two input positions landing on one variant. Counts per reason are in `gwas/<prefix>_std.liftover_qc.tsv` and the dropped rows in `gwas/<prefix>_std.liftover_dropped.tsv.gz`.
- Downstream pipelines include the same standardisation step automatically; run this pipeline alone when you only need harmonised files.
