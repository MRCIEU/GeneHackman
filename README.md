# GeneHackman

![CI Tests](https://github.com/MRCIEU/GeneHackman/actions/workflows/main.yml/badge.svg)

## Introduction

GeneHackman is a reproducible workflow suite for genetic epidemiology from GWAS summary statistics. It harmonises heterogeneous summary-statistic files, runs common downstream analyses (clumping, fine-mapping, colocalisation, Mendelian randomisation, and cross-cohort comparison), and writes structured outputs with HTML reports.

The workflows are implemented in [Snakemake](https://snakemake.readthedocs.io/) and executed inside a pinned [Docker / Apptainer image](https://hub.docker.com/repository/docker/mrcieu/genehackman/general) so that the same steps can run on a Linux server or a batch scheduler (ex. Slurm). Analysis logic lives in an R package and companion scripts; reference panels (1000 Genomes LD, liftover chains, LDSC weights) and optional QTL panels are fetched separately.

GeneHackman accompanies an application note: each pipeline corresponds to a distinct analysis scenario, with YAML configuration, documented inputs and outputs, and end-to-end tests on small public example data.

**Why use it**

- **Data Harmonisation** — map arbitrary GWAS column layouts to a shared schema, optional liftover, RSID population, EAF population, and alphabetically ordered effect alleles.
- **Reproducible pipelines** — standardise only, or chain standardisation with clumping, SuSiE / MultiSuSiE fine-mapping, BF-BF colocalisation, LDSC, or QTL MR.
- **Reproducible by default** — containerised runtime, versioned reference data, and explicit completion markers per analysis run.

If you use GeneHackman in published work, please cite the [Zenodo record](https://doi.org/10.5281/zenodo.10624713) and the application note (when available).

## Available Pipelines

There are **six** Snakemake pipelines (grouped as two tables of three). Each pipeline is a `.smk` file under `snakemake/`; see [PIPELINES.md](PIPELINES.md) for YAML inputs and parameters.

### Pipelines — table 1

| standardise_gwas.smk                                  | compare_gwases.smk                                                        | disease_progression.smk                                               |
|-------------------------------------------------------|---------------------------------------------------------------------------|----------------------------------------------------------------------|
| Takes in any of: VCF, CSV, TSV, TXT (also zip/gz)     | Runs **standardise_gwas** for each GWAS, then pairwise comparison tooling | Runs **standardise_gwas** for incident and subsequent GWASes        |
| Liftover with QC to GRCh38 (default) or GRCh37        | PLINK clumping                                                           | Runs collider-bias-aware analyses (e.g. SlopeHunter, MR-IVW) and compares results |
| Optional RSID and EAF fill from reference panels      | Heterogeneity across GWASes; LDSC *h*² and *r*<sub>g</sub>                | Miami plots of unadjusted vs adjusted GWAS                          |
| Converts z-scores and odds ratios to BETA/SE          | Expected vs observed replication metrics                                 | Expected vs observed (before and after collider adjustment)          |
| Harmonised SNP ID = `CHR:BP_EA_OA` (EA/OA sorted)     | HTML report of LDSC, plots, and tables                                    | HTML report; optional instructions to refit subsequent GWAS from collider output    |
| Optional gene ↔ Ensembl mapping                       |                                                                           |                                                                      |

### Pipelines — table 2

| qtl_mr.smk                                                                 | finemap.smk                                                                 | coloc.smk                                                                 |
|----------------------------------------------------------------------------|----------------------------------------------------------------------------|---------------------------------------------------------------------------|
| Runs **standardise_gwas**, clumping, and **SuSiE** fine-mapping on one outcome GWAS | Same **standardise** + **clump** + **SuSiE** path for **each** input GWAS (one or more) | Same **standardise** + **clump** + **SuSiE** for **≥2** GWASes (required for colocalization) |
| Runs Mendelian randomization vs the cis QTLs of eQTLGen, BrainMeta or UKB-PPP | Per-locus fine-mapping using summary stats and **ancestry-matched LD** (PLINK reference); outputs credible sets and LBF columns per locus | **Pairwise** `coloc::coloc.bf_bf` on overlapping finemapped signals (same chr, leads within ±`overlap_kb` kb) across all trait pairs |
| Volcano plot of MR results; BF-BF coloc for exposures that pass MR FDR | Finemap-only: no MR or coloc between traits (use when you only need SuSiE outputs) | Full coloc table + **HTML report** (`result_coloc.html`) |
| Requires **QTL_DATA_DIR** (see `.env_example`) for QTL files              | Each GWAS must declare **ancestry** (for LD)                                | Configurable `finemap` and `coloc` priors/overlap; see [PIPELINES.md](PIPELINES.md) |

## Quick start

### 1. Clone the repository

```bash
git clone https://github.com/MRCIEU/GeneHackman.git
cd GeneHackman
```

### 2. Create the conda environment

Ensure you [have conda installed and initialised](https://www.anaconda.com/docs/getting-started/miniconda/install/overview), then:

```bash
conda env create -f environment.yml
conda activate genehackman
```

### 3. Download reference data

The reference data (release 1.2.0, GRCh38) are public, and [`download_reference_data.sh`](download_reference_data.sh) downloads only what you need. It needs `curl` only, no Google Cloud tools.

```bash
# European ancestry, no QTL data: about 20 GB
./download_reference_data.sh --dir /path/to/my_pipeline_data

# More ancestries, and QTL data for the qtl_mr pipeline
./download_reference_data.sh --dir /path/to/my_pipeline_data --ancestries EUR,EAS --qtl eqtlgen,ukb_pqtl

# See what a command would download, and how big it is, without downloading anything
./download_reference_data.sh --dir /path/to/my_pipeline_data --ancestries all --qtl all --dry-run
```

The script prints the `.env` lines to set at the end (`PIPELINE_DATA_DIR`, and `QTL_DATA_DIR` if you downloaded QTL data). To add ancestries or datasets later, run it again with the new options: files already there are kept, and an interrupted download resumes. Every file is checked against its md5. Run `./download_reference_data.sh --help` for all options.

| Option | Default | What it downloads | Size |
|---|---|---|---|
| (always) | | liftover chains, gene name map, 1000 Genomes RSID marker map, HapMap3 SNP list | 0.2 GB |
| `--ancestries` | `EUR` | per ancestry: the 1000 Genomes 30x panel (GRCh38), the Pan-UKB LD pack (see [PIPELINES.md](PIPELINES.md)) and the Pan-UKB LD scores. Choose from `EUR EAS AFR AMR SAS MID`, or `all` | EUR 7.8 GB, EAS 7.3, SAS 8.7, AMR 16.9, MID 17.5, AFR 27.1 |
| `--no-dbsnp` | dbSNP included | leaves out dbSNP 156 (GRCh38), which only `populate_rsid: true` uses. The default RSID fill uses the 1000 Genomes marker map | saves 11.7 GB |
| `--qtl` | `none` | QTL datasets for `qtl_mr` (GRCh38): `eqtlgen`, `brainmeta`, `ukb_pqtl`, or `all`. They go in `--qtl-dir` (default: `QTL_DATA_DIR` from `.env`, else `<dir>/qtl_data`) | eqtlgen 1.8 GB, brainmeta 2.3, ukb_pqtl 0.4 |
| `--grch37` | not included | the GRCh37 files for `output.build: GRCh37`: the 1000 Genomes phase 3 panel and the earlier LD scores for your ancestries, plus dbSNP GRCh37 unless `--no-dbsnp` | about 1.9 GB per ancestry, plus 11.4 GB dbSNP |

Notes:
- **Ancestries.** `MID` exists only in Pan-UKB, so it works with `ld_reference: pan_ukb` but not `1000genomes`. South Asian is `SAS` in 1000 Genomes and `CSA` in Pan-UKB; `SAS` downloads both.
- **GRCh37.** GWASes in GRCh37 don't need `--grch37`: they are lifted to GRCh38 during standardisation. It is only for GRCh37 *output*. The QTL data are GRCh38 only, so `qtl_mr` always runs with `output.build: GRCh38`.
- **Missing data.** If a pipeline needs data you haven't downloaded (another ancestry, dbSNP, GRCh37 files or a QTL dataset), it stops before running and prints the `download_reference_data.sh` command that adds it.
- **Pan-UKB data.** The Pan-UKB LD pack and LD scores are derived from the [Pan-UKB](https://pan.ukbb.broadinstitute.org/downloads) LD matrices and LD scores, which Pan-UKB releases openly. Their use is bound by the UK Biobank terms of use; don't attempt to re-identify participants. If you use them, acknowledge Pan-UKB and UK Biobank and cite Karczewski et al., [Nat Genet 57, 2408–2417 (2025)](https://doi.org/10.1038/s41588-025-02335-7).
- **Earlier releases.** Release 1.1.0's data (GRCh37) stay at `gs://genehackman/1.1.0/` and `gs://genehackman-qtl/1.1.0/` for reproducing earlier results; see that release's README.

### 4. Configure `.env` and input YAML

Copy the template and edit paths for your machine:

```bash
cp .env_example .env
```

Variables match [.env_example](.env_example):

**Required**

| Variable | Purpose |
|----------|---------|
| **`PROJECT_DIR`** | Root folder for this analysis. The pipeline writes inputs under **`PROJECT_DIR/data/`** and results under **`PROJECT_DIR/results/`**. |
| **`PIPELINE_DATA_DIR`** | Directory containing **`genomic_data/`** and **`LDSCORE/`** (from step 3). Also used for the Apptainer image cache: **`PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<version>.sif`**, where `<version>` comes from **`DOCKER_VERSION`** or defaults to **`Version:`** in [`DESCRIPTION`](DESCRIPTION). The reference data come from `download_reference_data.sh` (step 3). If the SIF is missing and the directory is writable, `run_pipeline.sh` builds it from `docker://mrcieu/genehackman:<version>`. |

**Optional**

| Variable | Purpose |
|----------|---------|
| **`DOCKER_VERSION`** | Docker/Apptainer image tag (e.g. `1.2.0` or `develop`). Defaults to **`Version:`** in [`DESCRIPTION`](DESCRIPTION) when unset. |
| **`QTL_DATA_DIR`** | QTL data for **`qtl_mr`**, from `download_reference_data.sh --qtl` (step 3). Leave empty for other pipelines. |
| **`SNAKEMAKE_PROFILE`** | Snakemake profile directory. Default: **`snakemake/profiles/local/`** (local Apptainer). On a cluster, use e.g. **`snakemake/profiles/slurm/`**. |
| **`APPTAINER_MODULE`** | Environment module to load Apptainer/Singularity on HPC (when not using the `local` profile). |
| **`SLURM_PARTITION`** | Slurm partition. If unset, `run_pipeline.sh` tries `sinfo`, else **`compute`**. |
| **`SLURM_ACCOUNT`** | Slurm account. If unset, `run_pipeline.sh` tries to infer from `sacctmgr`. |

Example `.env`:

```bash
PROJECT_DIR=/path/to/my_project
PIPELINE_DATA_DIR=/path/to/my_pipeline_data/
# DOCKER_VERSION=1.2.0  # optional; defaults to Version: in DESCRIPTION
QTL_DATA_DIR=/path/to/my_qtl_data/
SNAKEMAKE_PROFILE=snakemake/profiles/local/
```

**Input YAML**

* Copy a template: `cp snakemake/input_templates/compare_gwases.yaml input.yaml`
* Examples for every pipeline: [`snakemake/input_templates/`](snakemake/input_templates/)
* Field reference: [PIPELINES.md](PIPELINES.md)
* Pass the YAML as the second argument to `run_pipeline.sh`, or use **`input.yaml`** in the working directory. When calling **`snakemake`** directly, pass **`--config genehackman_input=/path/to/file.yaml`** and ensure **`PROJECT_DIR`** is set in **`.env`**.

### 5. Run a pipeline

Smoke test with bundled example input:

```bash
./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas.yaml
```

Run your analysis:

```bash
./run_pipeline.sh snakemake/<pipeline>.smk [path/to/input.yaml]
```

Snakemake profiles live under **`snakemake/profiles/`** — see [PLATFORM_SETUP.md](PLATFORM_SETUP.md) for **Linux**, **HPC (Slurm)**, and **macOS** (not officially supported).

- Default: **`SNAKEMAKE_PROFILE=snakemake/profiles/local/`** (Apptainer on the current machine).
- On a cluster: e.g. **`SNAKEMAKE_PROFILE=snakemake/profiles/slurm/`**, or add a site-specific profile beside the bundled ones.
- `run_pipeline.sh` is a convenience wrapper around **`snakemake`**; see the [Snakemake documentation](https://snakemake.readthedocs.io/) for advanced options.
- Long jobs: consider a terminal multiplexer (e.g. `tmux`) so SSH disconnects do not kill the run.

### Trusted research environments

Biobank trusted research environments (TREs) keep individual-level data inside a controlled computing environment, and many let users bring their own software containers. GeneHackman runs from one container image and a YAML file, and the pipelines don't download anything while they run, so it can be set up once and then used without internet access:

1. **Container.** Outside the TRE, build the image (`apptainer build genehackman_1.2.0.sif docker://mrcieu/genehackman:1.2.0`) and copy it into the TRE as `PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_1.2.0.sif`. `run_pipeline.sh` uses the image if it is there, and only tries to download it if it is missing.
2. **Reference data.** Download what you need with `./download_reference_data.sh` outside the TRE (for example `--ancestries EUR`), and copy the directory in as `PIPELINE_DATA_DIR`.
3. **Snakemake.** Snakemake runs on the host from the `genehackman` conda environment (`environment.yml`): install it from the TRE's package mirror, or move it in with `conda-pack`. Use the `local` profile, or add a profile for the TRE's job scheduler under `snakemake/profiles/`.
4. **In-sample LD.** Give each GWAS its own `ld_reference:`, the PLINK prefix of the genotypes the GWAS was computed from, so that every LD step (strand check, clumping, fine-mapping, colocalisation and regional plots) uses in-sample LD rather than a reference panel (see [PIPELINES.md](PIPELINES.md#common-conventions) and [snakemake/finemap.md](snakemake/finemap.md)). LDSC still uses the Pan-UKB LD scores.
5. **What leaves the TRE.** The outputs are summary-level (harmonised summary statistics, credible sets, colocalisation, MR and LDSC results, and the HTML reports); exporting them follows the TRE's own rules. Usage logging is off unless you turn it on (below).

### Usage logging (opt-in)

GeneHackman doesn't send anything anywhere by default. If you set `GENEHACKMAN_USAGE_LOGGING=true` in `.env` and the developers' credentials file is installed, each run appends one row to the developers' usage sheet: user name, host name, pipeline, start time, run time and outcome, and for failed runs the input YAML and the error log, so the developers can help. Leave it unset to log nothing.

### Platform setup

See **[PLATFORM_SETUP.md](PLATFORM_SETUP.md)** for Linux, HPC (Slurm), and macOS setup (macOS is not officially supported).

### Contributing

See **[CONTRIBUTING.md](CONTRIBUTING.md)** for development setup, Docker rebuilds, tests, and adding new pipelines. Bug reports and feature requests: [GitHub issues](https://github.com/MRCIEU/GeneHackman/issues).
