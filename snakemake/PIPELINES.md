# Pipeline input (`input.yaml`)

**Execution profiles:** Snakemake **`--profile`** directories live under **`snakemake/profiles/`** from the repo root (e.g. **`snakemake/profiles/local/`**, **`snakemake/profiles/slurm/`**), selected via **`SNAKEMAKE_PROFILE`** when using **`./run_pipeline.sh`**. Details: [PLATFORM_SETUP.md](../PLATFORM_SETUP.md), [README.md](../README.md).

---

All pipelines consume a **YAML** file. Use **`./run_pipeline.sh <workflow>.smk [path/to/input.yaml]`** — the second argument is optional and defaults to **`input.yaml`** in the working directory. If you invoke **`snakemake` directly**, pass **`--config genehackman_input=path/to/input.yaml`** (or omit for default **`input.yaml`**). Copy an example from [`input_templates/`](input_templates/).

Indented nesting defines hierarchy (`gwases:` lists `- item` blocks). Use `true` / `false` for booleans. Quote a string if it could be parsed as a YAML 1.1 keyword (e.g. column names `yes`, `no`). Validate syntax with [yamllint](https://www.yamllint.com/) if needed.

Set **`PROJECT_DIR`** in **`.env`**: the pipeline uses **`PROJECT_DIR/data/`** (GWAS, clumps, …) and **`PROJECT_DIR/results/`** (finemap, coloc, plots, …).

See also: [GWAS harmonisation README](../README.md#how-it-works).

---

## GWAS objects (`gwases` array)

Each element describes one GWAS file.


| Field                  | Required | Notes                                                                                                                                                                                                                           |
| ---------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `file`                 | **Yes**  | Absolute path (recommended) or path relative to the working directory from which Snakemake runs.                                                                                                                                |
| `columns`              | No       | Maps logical names (`CHR`, `BETA`, …) to header names in your file, or a preset string (e.g. `metal`, `gwama`, `opengwas`). Omit to use defaults. See [predefined_column_maps.csv](../inst/extdata/predefined_column_maps.csv). |
| `N`                    | Varies   | Sample size; **required** for compare_gwases (LDSC). Used where finemapping/coloc need `N`.                                                                                                                                     |
| `ancestry`             | Varies   | One of `AFR`, `AMR`, `EAS`, `EUR`, `SAS`. **Required** for pipelines that use LD (clumping, finemap, coloc, qtl_mr) and whenever `populate_eaf` is enabled for that GWAS.                                                       |
| `build`                | No       | `GRCh36`, `GRCh37`, `GRCh38`. Default `GRCh37`.                                                                                                                                                                                 |
| `populate_rsid`        | No       | Per-GWAS override: `true` / `false`. With clumping and `false` at both levels, RSID population is **partial** (see below).                                                                                                      |
| `populate_eaf`         | No       | Per-GWAS: if `true`, fills missing `EAF` from the 1000 Genomes reference (`<PIPELINE_DATA_DIR>/genomic_data/1000genomes/b37_dbsnp156/` … `<ancestry>.bim` + `.frq`). Requires `ancestry`.                                                                               |
| `flip_alleles`         | No       | If omitted, inherited from root (default **`true`**). **`false`** is **allowed only for [`standardise_gwas.smk`](standardise_gwas.smk)** (other pipelines fail at YAML load). With `false`, EA/OA order is preserved; **partial** RSID + `flip_alleles: false` remains invalid (use **`none`** or **`full`** RSID). |
| `remove_extra_columns` | No       | Default `false`.                                                                                                                                                                                                                |


**Effect / P-value columns** (see also top-level `output`):

- Required: `CHR`, `BP`, `EA`, `OA`
- One of: `BETA`+`SE`, or `OR`+`OR_LB`+`OR_UB`, or `Z`
- P-value: `P` or `LOG_P`

Default column names: `SNP, CHR, BP, EA, OA, EAF, BETA, SE, OR, P, LOG_P, Z, OR_LB, OR_UB, RSID, N, N_CASES, ENSEMBL_ID, GENE_NAME`.

---

## Root-level fields (shared)

These may appear at the **root** of the YAML document (alongside `gwases`).


| Field                          | Default                              | Notes                                                                                                                                                                                        |
| ------------------------------ | ------------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `is_test`                      | `false`                              | Used by tests / logging.                                                                                                                                                                     |
| `populate_rsid`                | `false`                              | Pipeline-wide default; each GWAS can override. **Clumping + false** here and on each GWAS ⇒ RSID population is **partial** (chrpos-based) unless you set `populate_rsid: true` where needed. |
| `populate_eaf`                 | `false`                              | Pipeline-wide default; each GWAS can override. Any GWAS with `populate_eaf: true` must include `ancestry`.                                                                                   |
| `flip_alleles`                 | `true`                               | Pipeline-wide default. Root or per-GWAS **`false`** only for **`standardise_gwas.smk`**; downstream pipelines require default harmonised alleles (`true`).                                                                                                      |
| `output`                       | `build: GRCh37` + default column map | Only `standardise_gwas` / harmonisation: target `build` and/or output column names.                                                                                                          |
| `finemap`                      | See below                            | Present in **finemap**, **coloc**, and **qtl_mr** inputs. Omitted keys get defaults from the parser (e.g. `window_kb` defaults to **500** unless you set it under `finemap`).                |
| `coloc`                        | See below                            | Only **coloc** pipeline.                                                                                                                                                                     |
| `plink_clump_arguments`        | —                                    | String passed to `plink --clump` for pipelines that clump. [PLINK clump options](https://zzz.bwh.harvard.edu/plink/clump.shtml).                                                             |
| `qtl`                          | —                                    | Only **qtl_mr** (see below).                                                                                                                                                                 |
| `output` (disease_progression) | —                                    | `adjusted_gwas`: `type` (`slopehunter`, `cwls`, `mr_ivw`) and `p_val`.                                                                                                                       |


### `finemap` block (optional)

Used by **finemap**, **coloc**, and **qtl_mr**. If the block is missing, defaults are filled by `parse_pipeline_input` (e.g. `window_kb` **500** unless set).


| Key            | Meaning                                                              | Default (if missing) |
| -------------- | -------------------------------------------------------------------- | -------------------- |
| `window_kb`    | Half-width of the fine-mapping window around each clumped lead (kb). | `500`                |
| `max_causal`   | SuSiE `L` (max causal signals per locus).                            | `10`                 |
| `coverage`     | Credible set coverage.                                               | `0.95`               |
| `min_abs_corr` | Minimum |r| for credible-set purity.                                 | `0.5`                |


### `coloc` block (optional, **coloc** pipeline only)


| Key          | Meaning                                                                                                          | Default |
| ------------ | ---------------------------------------------------------------------------------------------------------------- | ------- |
| `overlap_kb` | Two finemapped signals on the same chromosome are “overlapping” if lead positions are within this distance (kb). | `1000`  |
| `p1`         | Coloc prior: SNP associated with trait 1.                                                                        | `1e-4`  |
| `p2`         | Coloc prior: SNP associated with trait 2.                                                                        | `1e-4`  |
| `p12`        | Coloc prior: SNP associated with both traits.                                                                    | `5e-6`  |


---

## Pipeline overview


| Snakemake file                                          | Example template                                                                       | One-line role                                                                       |
| ------------------------------------------------------- | -------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| [`standardise_gwas.smk`](standardise_gwas.smk) | [`input_templates/standardise_gwases.yaml`](input_templates/standardise_gwases.yaml) | Harmonise GWAS columns and builds only (no PLINK). |
| [`compare_gwases.smk`](compare_gwases.smk) | [`input_templates/compare_gwases.yaml`](input_templates/compare_gwases.yaml) | Standardise → clump → heterogeneity / LDSC / expected–observed plots. |
| [`disease_progression.smk`](disease_progression.smk) | [`input_templates/disease_progression.yaml`](input_templates/disease_progression.yaml) | Incident vs subsequent GWAS: collider bias, plots, comparisons. |
| [`finemap.smk`](finemap.smk) | [`input_templates/finemap.yaml`](input_templates/finemap.yaml) | Standardise → clump → SuSiE RSS fine-mapping per GWAS (no pairwise coloc). |
| [`coloc.smk`](coloc.smk) | [`input_templates/coloc.yaml`](input_templates/coloc.yaml) | Same as finemap prep for ≥2 GWAS, then pairwise `coloc::coloc.bf_bf` + HTML report. |
| [`qtl_mr.smk`](qtl_mr.smk) | [`input_templates/qtl_mr.yaml`](input_templates/qtl_mr.yaml) | One GWAS: MR vs a QTL panel, volcano plot, BF–BF coloc for MR hits. |


---

## `standardise_gwas`

Harmonisation only (`snakemake/standardise_gwas.smk`).

- **GWAS**: one or more objects (see above).
- **Root**: optional `populate_rsid`, `populate_eaf`, `flip_alleles`, `output.build` / `output.columns`, `is_test`. This is the **only** pipeline that allows **`flip_alleles: false`**; when combined with **`populate_eaf`**, EAF merge runs **before** the SNP id is rebuilt without allele flipping.
- **Outputs**: `{DATA_DIR}/gwas/<prefix>_std.tsv.gz` per input file.

---

## `compare_gwases`

- **GWAS**: two or more; include `**N`** on each for LDSC.
- **Root**: `plink_clump_arguments`, optional `populate_rsid` / `populate_eaf`.
- **Behaviour**: same harmonisation as other pipelines, then PLINK clumping, heterogeneity, LDSC *h*² / *r*g, expected vs observed, HTML report.

---

## `disease_progression`

- **GWAS**: exactly **two** (incident and subsequent).
- **Root**: `plink_clump_arguments`, optional `populate_rsid`, and `**output.adjusted_gwas`**: `type` ∈ `slopehunter`, `cwls`, `mr_ivw`; `p_val` ∈ `0.1`, `0.01`, `0.001`, `1e-5`.

---

## `finemap`

SuSiE fine-mapping per GWAS — no MR, no trait–trait coloc.

- **GWAS**: one or more; each needs `**ancestry`** (LD reference for clumping and LD matrix for `susieR::susie_rss`).
- **Root**: `plink_clump_arguments`, optional `populate_rsid` / `populate_eaf`, optional `**finemap`** block (see above). Example template sets `window_kb` to `1000`.
- **Outputs** (under `RESULTS_DIR`): `finemap/<prefix>/` per GWAS (per-locus `*_finemap.tsv.gz`, plus `finemap_complete.txt` with one line = expected lead count `n`). Snakemake registers only that completion file as the finemap rule output (`FINEMAP_COMPLETE_TXT_PATTERN` in `snakemake/util/constants.smk`).

---

## `coloc`

Standardise → clump → finemap for **all** input GWAS, then **pairwise** BF–BF colocalisation between overlapping finemapped signals (`coloc::coloc.bf_bf`).

- **GWAS**: **at least two**; each needs `**ancestry`**.
- **Root**: `plink_clump_arguments`, optional `populate_rsid` / `populate_eaf`, optional `**finemap`** and `**coloc**` blocks (see tables above).
- **Outputs** (under `RESULTS_DIR`): `coloc/coloc_results.tsv`, `coloc/result_coloc.html` (summary table sorted by **PP.H4.abf**; includes a **cross-ancestry disclaimer** when more than one ancestry code is used across GWAS).

---

## `qtl_mr`

- **GWAS**: exactly **one** outcome GWAS; `**ancestry`** required.
- **Environment**: QTL summary data from `**QTL_DATA_DIR`** (see `.env` / `.env_example`); if unset, defaults follow `PIPELINE_DATA_DIR/qtl_datasets` layout (`pqtl`, `metabrain`, `eqtlgen`, …).
- **Root**: `plink_clump_arguments`, optional `**finemap`** (same keys as finemap pipeline), and `**qtl**`:
  - `dataset` — e.g. `metabrain`, `eqtlgen`
  - `subcategory` — one sub-run, e.g. `cortex` (MetaBrain) or `cis` (eQTLGen); combined with `dataset` in output names
  - `exposures` — optional list of exposure names; if empty, all exposures in that dataset run
  - `study_type` — `continuous` (default) or `categorical` (affects LBF conversion from QTL *Z*-scores)

**Dataset / subcategory examples** (see also repository docs):

- **metabrain**: `basalganglia`, `cerebellum`, `cortex`, `hippocampus`, `spinalcord`
- **eqtlgen**: `cis` (trans may be added later)

Pipeline flow: standardise → clump → finemap → MR → volcano → BF–BF coloc on significant MR results (GWAS LBF vs QTL LBF from `convert_z_to_lbf`).

---

## References

- Expected vs observed: doi:10.1038/nature17671  
- LDSC: [github.com/bulik/ldsc](https://github.com/bulik/ldsc)  
- `coloc`: [chr1swallace.github.io/coloc](https://chr1swallace.github.io/coloc/)  
- Pipeline DOI: [10.5281/zenodo.10624713](https://doi.org/10.5281/zenodo.10624713)

