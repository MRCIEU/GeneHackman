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
| `ancestry` | **Yes** | LD reference for clumping and fine-mapping (`AFR`, `AMR`, `EAS`, `EUR`, `SAS`). |
| `N` | Recommended | Sample size passed to SuSiE / MultiSuSiE. |
| `columns` | No | Column map or preset. |
| `build` | No | Input build. Default `GRCh37`. |
| `populate_rsid` | No | Per-GWAS override; default inherits root `false`. |
| `populate_eaf` | No | Fill missing `EAF` from 1000 Genomes (requires `ancestry`). |

### `plink_clump_arguments`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `plink_clump_arguments` | **Yes** | Passed to `plink1.9 --clump`. [PLINK clump options](https://zzz.bwh.harvard.edu/plink/clump.shtml). |

### `finemap`

| Field | Required | Notes |
| ----- | -------- | ----- |
| `finemap.window_kb` | No | Half-width of fine-mapping window around each lead (kb). Default **500** in parser; templates often use **1000**. |
| `finemap.max_causal` | No | Max causal signals per locus (`L`). Default `10`. |
| `finemap.coverage` | No | Credible set coverage. Default `0.95`. |
| `finemap.min_abs_corr` | No | Minimum \|r\| for credible-set purity. Default `0.5`. |

### Root fields

| Field | Required | Notes |
| ----- | -------- | ----- |
| `populate_rsid` | No | Default `false`. With clumping and `false`, RSID population is **partial** unless overridden per GWAS. |
| `populate_eaf` | No | Default `false`. |
| `output.build` | No | Target build after standardisation. Default `GRCh37`. |

`flip_alleles: false` is **not** allowed.

## Workflow

1. Standardise each GWAS
2. PLINK clump per GWAS
3. Fine-map:
   - **Single ancestry:** `run_finemap.R` → SuSiE RSS per GWAS
   - **Multi ancestry:** `run_multisusie.py` → MultiSuSiE across all GWAS jointly

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
| `finemap/<prefix>/finemap_complete_<prefix>.txt` | **Snakemake target.** One line = number of clumped loci processed. |
| `finemap/<prefix>/<CHR>_<BP>_finemap.tsv.gz` | Per-locus fine-mapping table (LBF columns, credible-set membership, posterior effects). |

### Results — multi-ancestry MultiSuSiE (`results/finemap/multi_ancestry/`)

Full documentation for interpreting MultiSuSiE results can be [found here](https://deepwiki.com/jordanero/MultiSuSiE/4.3-interpreting-results)

| Path | Description |
| ---- | ----------- |
| `finemap/multi_ancestry/finemap_complete_<gwas_prefixes>.txt` | **Snakemake target.** One line = number of successfully fine-mapped loci. `<gwas_prefixes>` is sorted input GWAS file prefixes joined with `_`. |
| `finemap/multi_ancestry/<CHR>_<center>_locus_credible_sets.tsv` | One row per credible set: size, top variant PIP/LBF, per-ancestry purity, coverage. |
| `finemap/multi_ancestry/<CHR>_<center>_locus_credible_set_variants.tsv` | One row per variant in any credible set: coordinates, alleles, global PIP, per-ancestry `COEF`/`COEF_SD`, and input `BETA`/`SE`/`P`. |

Locus names use the clump window centre (`<CHR>_<center>`).

## Notes

- For pairwise colocalisation between traits, use [`coloc.smk`](coloc.smk) instead.
- Multi-ancestry mode does **not** produce per-trait LBF columns suitable for `coloc::coloc.bf_bf` on individual GWAS directories.
- LD matrices use ancestry-matched 1000 Genomes panels under `PIPELINE_DATA_DIR/genomic_data/1000genomes/b37_dbsnp156/`.
