# Pipeline input and workflows

GeneHackman runs as Snakemake workflows (`.smk` files under `snakemake/`). Each workflow consumes a **YAML** input file and writes under **`PROJECT_DIR/data/`** and **`PROJECT_DIR/results/`** (set in `.env`).

**Run:**

```bash
./run_pipeline.sh snakemake/<workflow>.smk [path/to/input.yaml]
```

The second argument defaults to `input.yaml`. With bare `snakemake`, pass `--config genehackman_input=path/to/input.yaml`.

Starter YAML files live in [`snakemake/input_templates/`](snakemake/input_templates/). Each pipeline doc below inlines the template and explains **only the fields used by that workflow**.

**Profiles:** Snakemake `--profile` directories live under `snakemake/profiles/` (see [PLATFORM_SETUP.md](PLATFORM_SETUP.md), [README.md](README.md)).

**Harmonisation background:** [README — How it works](README.md#how-it-works).

---

## Pipelines

| Workflow | Doc | Role |
| -------- | --- | ---- |
| [`standardise_gwas.smk`](snakemake/standardise_gwas.smk) | [standardise_gwas.md](snakemake/standardise_gwas.md) | Harmonise columns and reference build only. |
| [`compare_gwases.smk`](snakemake/compare_gwases.smk) | [compare_gwases.md](snakemake/compare_gwases.md) | Standardise → clump → heterogeneity, LDSC, expected vs observed. |
| [`disease_progression.smk`](snakemake/disease_progression.smk) | [disease_progression.md](snakemake/disease_progression.md) | Incident vs subsequent GWAS; collider bias correction. |
| [`finemap.smk`](snakemake/finemap.smk) | [finemap.md](snakemake/finemap.md) | Standardise → clump → SuSiE or MultiSuSiE fine-mapping. |
| [`coloc.smk`](snakemake/coloc.smk) | [coloc.md](snakemake/coloc.md) | Finemap ≥2 GWAS → pairwise coloc. |
| [`qtl_mr.smk`](snakemake/qtl_mr.smk) | [qtl_mr.md](snakemake/qtl_mr.md) | One outcome GWAS → MR vs QTL panel → coloc on MR hits. |

---

## Common conventions

These apply to every pipeline unless a workflow doc says otherwise.

**File formats:** `gwases[].file` supports VCF, CSV, TSV, and TXT (`.gz` / `.zip`). Absolute paths are recommended.

**Column mapping:** set `columns` to a header map or a preset string (`metal`, `gwama`, `opengwas_phewas`, …). See [`inst/extdata/predefined_column_maps.csv`](inst/extdata/predefined_column_maps.csv).

**Minimum columns after mapping:**

- Required: `CHR`, `BP`, `EA`, `OA`
- Effect: `BETA`+`SE`, or `OR`+`OR_LB`+`OR_UB`, or `Z`
- P-value: `P` or `LOG_P`

Default logical names: `SNP`, `CHR`, `BP`, `EA`, `OA`, `EAF`, `BETA`, `SE`, `OR`, `P`, `LOG_P`, `Z`, `OR_LB`, `OR_UB`, `RSID`, `N`, `N_CASES`, `N_CONTROLS`, `N_EFF`, `ENSEMBL_ID`, `GENE_NAME`. For binary traits, per-variant `N_CASES` with `N_CONTROLS`, or `N_CASES` with a total `N`, give a per-variant effective sample size.

**Allele harmonisation:** root `flip_alleles` defaults to `true`. Only [`standardise_gwas.smk`](snakemake/standardise_gwas.smk) accepts `flip_alleles: false`.

**Clumping:** pipelines that clump require `plink_clump_arguments`. If `--clump-r2` or `--clump-kb` are not given, GeneHackman adds `--clump-r2 0.01 --clump-kb 1000` (PLINK's own defaults, 0.5 and 250 kb, leave several correlated lead SNPs per locus).

**LD reference:** `ld_reference` is **required** (there is no default) in every pipeline that clumps, fine-maps or checks alleles against a reference panel. It chooses where clumping, fine-mapping (SuSiE and MultiSuSiE), LocusZoom, strand checks and EAF filling get LD:

| Value | Source | Variants are matched on | Ancestries |
|---|---|---|---|
| `1000genomes` | 1000 Genomes genotypes via PLINK, in the output build: the 30x high-coverage release for GRCh38, phase 3 for GRCh37 | RSID | EUR, EAS, AFR, AMR, SAS |
| `pan_ukb` | Precomputed Pan-UKB LD: 6,000–420,000 people per ancestry, MAF ≥ 1%, ±1 Mb, GRCh37 and GRCh38 | SNP key (`CHR:BP_EA_OA`) | EUR, AFR, AMR, CSA, EAS, MID (`SAS` maps to `CSA`) |

The sample YAMLs use `pan_ukb`, the larger reference. Reference LD can miscalibrate fine-mapping when the GWAS is much larger than the panel, so when a GWAS is more than `finemap.max_n_ratio` (default 10) times larger than its LD reference, the fine-mapping, colocalisation and QTL pipelines warn in the log and in their HTML reports; they don't change any setting.

The Pan-UKB pack must be at `PIPELINE_DATA_DIR/genomic_data/pan_ukb_ld/`; see [`reference_data/pan_ukb_ld/README.md`](reference_data/pan_ukb_ld/README.md) for how it is built.

With `pan_ukb`:

- **Clumping** follows PLINK 1.9's `--clump` rules on the stored LD, and `--clump-kb` is capped at the pack's radius (1,000 kb, so the default `--clump-kb 1000` is not cut).
- **Fine-mapping:** variant pairs further apart than the radius are set to r = 0. Because r is stored rounded, the LD matrix is made positive semi-definite by setting its negative eigenvalues to 0. `finemap_regions.tsv` reports the share of pairs beyond the radius (`FRAC_BEYOND_LD_RADIUS`) and the smallest eigenvalue before that fix (`LD_MIN_EIGENVALUE`).
**Genome build:** outputs are in `output.build`, **GRCh38 by default**. Each GWAS declares its own `build` (default GRCh37) and is lifted to the output build during standardisation, with QC: variants whose new position can't be trusted are dropped and counted (see [standardise_gwas.md](snakemake/standardise_gwas.md)). The declared build is checked against the Pan-UKB LD pack before lifting. Every reference resource is then used in the output build: the 1000 Genomes panel, the RSID maps, the Pan-UKB LD keys, the LocusZoom gene tracks (Ensembl 86 for GRCh38, 75 for GRCh37) and the QTL data. `output.build: GRCh37` is still supported, but needs the GRCh37 files (`download_reference_data.sh --grch37`); `qtl_mr` is GRCh38 only, because the QTL data are. An already standardised `*_std.tsv.gz` input must declare `build` equal to `output.build`.

**RSIDs:** RSIDs are looked up in the reference for the output build. Partial RSID population (the default for pipelines that clump) uses the 1000 Genomes marker map of that build; `populate_rsid: true` uses the full dbSNP lookup.

**LDSC** (`compare_gwases`): `ldsc.ld_scores` chooses the LD scores, `pan_ukb` (default: Pan-UKB HapMap3 LD scores for EUR, AFR, AMR, CSA/SAS, EAS, MID) or `1000genomes_b37` (the earlier 1000 Genomes scores, downloaded with `--grch37`). LDSC matches variants on RSID, so either works with either output build.

---

## Standard intermediate paths

Most pipelines write harmonised GWAS and clumps to fixed locations:

| Path | Description |
| ---- | ----------- |
| `data/gwas/<prefix>_std.tsv.gz` | Standardised GWAS (`<prefix>` = input file basename). |
| `data/clumped_snps/<prefix>.clumped` | PLINK clump output (when the pipeline clumps). |

Fine-mapping completion markers:

| Path | Description |
| ---- | ----------- |
| `results/finemap/<prefix>/finemap_complete_<prefix>.txt` | Single-ancestry SuSiE done (one line = locus count). |
| `results/finemap/multi_ancestry/finemap_complete_<gwas_prefixes>.txt` | Multi-ancestry MultiSuSiE done; label is sorted GWAS prefixes joined with `_`. |

---

Pipeline DOI: [10.5281/zenodo.10624713](https://doi.org/10.5281/zenodo.10624713)
