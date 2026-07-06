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
| Optional liftover (e.g. GRCh38 → GRCh37)              | PLINK clumping                                                           | Runs collider-bias-aware analyses (e.g. SlopeHunter, MR-IVW) and compares results |
| Optional RSID and EAF fill from reference panels      | Heterogeneity across GWASes; LDSC *h*² and *r*<sub>g</sub>                | Miami plots of unadjusted vs adjusted GWAS                          |
| Converts z-scores and odds ratios to BETA/SE          | Expected vs observed replication metrics                                 | Expected vs observed (before and after collider adjustment)          |
| Harmonised SNP ID = `CHR:BP_EA_OA` (EA/OA sorted)     | HTML report of LDSC, plots, and tables                                    | HTML report; optional instructions to refit subsequent GWAS from collider output    |
| Optional gene ↔ Ensembl mapping                       |                                                                           |                                                                      |

### Pipelines — table 2

| qtl_mr.smk                                                                 | finemap.smk                                                                 | coloc.smk                                                                 |
|----------------------------------------------------------------------------|----------------------------------------------------------------------------|---------------------------------------------------------------------------|
| Runs **standardise_gwas**, clumping, and **SuSiE** fine-mapping on one outcome GWAS | Same **standardise** + **clump** + **SuSiE** path for **each** input GWAS (one or more) | Same **standardise** + **clump** + **SuSiE** for **≥2** GWASes (required for colocalization) |
| Runs Mendelian randomization vs a chosen QTL panel (e.g. eQTLGen, MetaBrain) | Per-locus fine-mapping using summary stats and **ancestry-matched LD** (PLINK reference); outputs credible sets and LBF columns per locus | **Pairwise** `coloc::coloc.bf_bf` on overlapping finemapped signals (same chr, leads within ±`overlap_kb` kb) across all trait pairs |
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

Reference files (1000 Genomes LD panels, liftover chains, LDSC weights) are hosted on Google Cloud Storage under a **version prefix** that matches the pipeline release, e.g. `gs://genehackman/1.1.0. [Install gsutil](https://docs.cloud.google.com/storage/docs/gsutil_install) if needed.

**Required for most pipelines**

```bash
gsutil -m rsync -r "gs://genehackman/1.1.0/" /path/to/my_pipeline_data/
```

Set `PIPELINE_DATA_DIR=/path/to/my_pipeline_data/` in `.env`.

**Optional — QTL MR pipeline only**

```bash
gsutil -m rsync -r "gs://genehackman-qtl/1.1.0/" /path/to/my_qtl_data/
```

Set `QTL_DATA_DIR=/path/to/my_qtl_data/`, or place QTL data under `PIPELINE_DATA_DIR/qtl_datasets/` (see [PLATFORM_SETUP.md](PLATFORM_SETUP.md)).

Example `.env` paths:

```bash
PIPELINE_DATA_DIR=/path/to/my_pipeline_data/
QTL_DATA_DIR=/path/to/my_qtl_data/
```

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
| **`PIPELINE_DATA_DIR`** | Directory containing **`genomic_data/`** and **`LDSCORE/`** (from step 3). Also used for the Apptainer image cache: **`PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<version>.sif`**, where `<version>` comes from **`DOCKER_VERSION`** or defaults to **`Version:`** in [`DESCRIPTION`](DESCRIPTION). Download reference data from **`gs://genehackman/1.1.0/`** with the same version string. If the SIF is missing and the directory is writable, `run_pipeline.sh` builds it from `docker://mrcieu/genehackman:<version>`. |

**Optional**

| Variable | Purpose |
|----------|---------|
| **`DOCKER_VERSION`** | Docker/Apptainer image tag (e.g. `1.1.0` or `develop`). Defaults to **`Version:`** in [`DESCRIPTION`](DESCRIPTION) when unset. |
| **`QTL_DATA_DIR`** | Path to **`genehackman-qtl`** data for **`qtl_mr`**. Leave empty for other pipelines. |
| **`SNAKEMAKE_PROFILE`** | Snakemake profile directory. Default: **`snakemake/profiles/local/`** (local Apptainer). On a cluster, use e.g. **`snakemake/profiles/slurm/`**. |
| **`APPTAINER_MODULE`** | Environment module to load Apptainer/Singularity on HPC (when not using the `local` profile). |
| **`SLURM_PARTITION`** | Slurm partition. If unset, `run_pipeline.sh` tries `sinfo`, else **`compute`**. |
| **`SLURM_ACCOUNT`** | Slurm account. If unset, `run_pipeline.sh` tries to infer from `sacctmgr`. |

Example `.env`:

```bash
PROJECT_DIR=/path/to/my_project
PIPELINE_DATA_DIR=/path/to/my_pipeline_data/
# DOCKER_VERSION=1.1.0  # optional; defaults to Version: in DESCRIPTION
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

### Platform setup

See **[PLATFORM_SETUP.md](PLATFORM_SETUP.md)** for Linux, HPC (Slurm), and macOS setup (macOS is not officially supported).

### Contributing

See **[CONTRIBUTING.md](CONTRIBUTING.md)** for development setup, Docker rebuilds, tests, and adding new pipelines. Bug reports and feature requests: [GitHub issues](https://github.com/MRCIEU/GeneHackman/issues).
