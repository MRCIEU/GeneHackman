# GeneHackman 

![CI Tests](https://github.com/MRCIEU/GeneHackman/actions/workflows/main.yml/badge.svg)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10624713.svg)](https://doi.org/10.5281/zenodo.10624713)

A pipeline for performing common genetic epidemiology tasks at the University of Bristol.

Goals:
* Implement comment steps in GWAS investigations to create reproducible, more efficient research
* Reusable pipelines that can be utilised on different projects
* Shared code and steps that can be updated according to the latest knowledge and practices

## Available Pipelines

There are **six** Snakemake pipelines (grouped as two tables of three). Each pipeline is a `.smk` file under `snakemake/`; see [PIPELINES.md](snakemake/PIPELINES.md) for YAML inputs and parameters.

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
| Volcano plot of MR results; **BF-BF coloc** for exposures that pass MR FDR | Finemap-only: no MR or coloc between traits (use when you only need SuSiE outputs) | Full coloc table + **HTML report** (`result_coloc.html`), including a disclaimer when ancestries differ between GWASes |
| Requires **QTL_DATA_DIR** (see `.env_example`) for QTL files              | Each GWAS must declare **ancestry** (for LD)                                | Configurable `finemap` and `coloc` priors/overlap; see [PIPELINES.md](snakemake/PIPELINES.md) |

## Onboarding

### 1. Clone the repository into your personal space on BlueCrystal 4
`git clone git@github.com:MRCIEU/GeneHackman.git && cd GeneHackman`

### 2. Ensure you [have conda installed and initialised before activating](https://www.anaconda.com/docs/getting-started/miniconda/install/overview)

```conda env create --file environment.yml```

or if you have already created the environment

```conda activate genehackman```

### 3. Get reference data (Google Cloud Storage)

The data to run the pipelines have been split into two buckets, the mandatory bucket, and QTL bucket (only needed if you want to run MR-QTL pipeline).  To download, [install gsutil](https://docs.cloud.google.com/storage/docs/gsutil_install)

* Mandatory: `gs://genehackman`
  * `gsutil -m rsync -r gs://genehackman/ /path/to/my_pipeline_data/`
  * Update `PIPELINE_DATA_DIR` to `/path/to/my_pipeline_data/` in the .env file
* Optional: `gs://genehackman-qtl` 
  * `gsutil -m rsync -r gs://genehackman-qtl/ /path/to/my_qtl_data/`
  * Update `QTL_DATA_DIR` to `/path/to/my_qtl_data/` in the .env file

To copy only selected prefixes instead of the QTL bucket, use `gsutil -m cp -r gs://genehackman-qtl/SOME_PREFIX/ ...` as needed.  For example, you may only be interested in cis, not trans data.

Then point your **`.env`** at those directories (trailing slashes are fine):

```bash
PIPELINE_DATA_DIR=/path/to/my_pipeline_data/
QTL_DATA_DIR=/path/to/my_qtl_data/
```

Alternatives if you don’t set **`QTL_DATA_DIR`** directly: keep the same directory layout under **`PIPELINE_DATA_DIR/qtl_datasets/`**, or follow [PLATFORM_SETUP.md](PLATFORM_SETUP.md) (bind mounts and defaulting **`QTL_DATA_DIR`** to **`PIPELINE_DATA_DIR/qtl_datasets`**).

### 4. Populate `.env` and `input.yaml` files

Copy the template and edit paths for your machine:

```bash
cp .env_example .env
```

Variables match [.env_example](.env_example):

**Mandatory**

| Variable | Purpose |
|----------|---------|
| **`PROJECT_DIR`** | Root folder for this analysis. The pipeline uses **`PROJECT_DIR/data/`** for inputs (GWAS, clumps, …) and **`PROJECT_DIR/results/`** for outputs (finemap, coloc, plots, …). |
| **`PIPELINE_DATA_DIR`** | Path where you unpacked the shared **`genehackman`** reference data (see §3). Also used for the Apptainer image: **`PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<DOCKER_VERSION>.sif`**. If that SIF is missing and the directory is writable, `run_pipeline.sh` builds it from `docker://mrcieu/genehackman:<DOCKER_VERSION>`. |
| **`DOCKER_VERSION`** | Docker/Apptainer image tag (e.g. `1.1.0` or `develop`). Must match the SIF filename and the image you intend to run. |

**Optional**

| Variable | Purpose |
|----------|---------|
| **`QTL_DATA_DIR`** | Path to **`genehackman-qtl`** data if you run **`qtl_mr`**. Leave empty for other pipelines. You can instead place QTL data under **`PIPELINE_DATA_DIR/qtl_datasets/`** (see [PLATFORM_SETUP.md](PLATFORM_SETUP.md)). |
| **`SNAKEMAKE_PROFILE`** | Snakemake profile directory. Default in `.env_example`: **`snakemake/profiles/local/`** (local Apptainer). On HPC use e.g. **`snakemake/profiles/slurm/`**. |
| **`APPTAINER_MODULE`** | Environment module name to load Apptainer/Singularity before running on HPC (only used when the profile is not `local`). |
| **`SLURM_PARTITION`** | Slurm partition for cluster jobs. If unset, `run_pipeline.sh` tries to detect one from `sinfo`, else uses **`compute`**. |
| **`SLURM_ACCOUNT`** | Slurm account for cluster jobs. If unset, `run_pipeline.sh` tries to infer it from `sacctmgr`. |

Example `.env`:

```bash
PROJECT_DIR=/path/to/my_project
PIPELINE_DATA_DIR=/path/to/my_pipeline_data/
DOCKER_VERSION=1.1.0
QTL_DATA_DIR=/path/to/my_qtl_data/
SNAKEMAKE_PROFILE=snakemake/profiles/local/
```

**`input.yaml`**

* Example: `cp snakemake/input_templates/compare_gwases.yaml input.yaml`
* Each pipeline has its own shape; examples live under [`snakemake/input_templates/`](snakemake/input_templates/).
* See [PIPELINES.md](snakemake/PIPELINES.md) for all fields.
* Pass the input YAML as the **second** argument to `run_pipeline.sh`, or rely on **`input.yaml`** in the working directory. To call **`snakemake`** yourself without the wrapper, set **`PROJECT_DIR`** in **`.env`** and export **`DATA_DIR="${PROJECT_DIR%/}/data"`** and **`RESULTS_DIR="${PROJECT_DIR%/}/results"`** (or load the same paths Snakemake uses) so shell rules and profile bind mounts resolve; also pass **`--config genehackman_input=/path/to/file.yaml`**.

### 5. Run the pipeline

To ensure you have configured everything correctly, you can run a test pipeline from `run_test_pipelines.sh`

`./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas.yaml`

To run your pipeline:

`./run_pipeline.sh snakemake/<specific_pipeline>.smk <optional_input_file.yaml>`

Snakemake execution profiles (**`--profile`**) live under **`snakemake/profiles/`** (see also [PLATFORM_SETUP.md](PLATFORM_SETUP.md)).

* By default `run_pipeline.sh` uses **`SNAKEMAKE_PROFILE=snakemake/profiles/local/`** (local Apptainer). On HPC, set e.g. **`SNAKEMAKE_PROFILE=snakemake/profiles/slurm/`** (generic Slurm) or add a site-specific directory beside **`snakemake/profiles/local/`** and **`snakemake/profiles/slurm/`**.
* `run_pipeline.sh` is just a convience wrapper around the `snakemake` command, if you want to do anything out of the ordinary, [please read up on snakemake](https://snakemake.readthedocs.io/en/v7.26.0/)
* If there are errors while running the pipeline, you can find error messages either directly on the screen, or in slurm log file that is outputted on error
* It is recommended that you run the your pipeline [inside a tmux session](https://github.com/MRCIEU/GeneHackman/wiki/Common-Errors#ssh-disconnection-while-pipeline-is-running).

## How it works:

The standard column naming for GWASes are:

| CHR | BP  | EA  | OA  | BETA | SE  | P   | EAF | SNP | RSID |
|-----|-----|-----|-----|------|-----|-----|-----|-----|:-----|

A full list of names and default values [can be found here](inst/extdata/predefined_column_maps.csv)

There are 2 main components to the pipeline
1. Snakemake to define the steps to complete for each pipeline.
2. Docker / Singularity container with installed languages (R and python), packages, os libraries, and code

The pipeline can be run either on its own, or via your institutions HPC.  Each snakemake step spins up a singularity container inside an HPC job (ex. slurm).  Each step can specify different cpu/memory requirements.

### Repository Organisation

* `R` directory holds R package code that can also be called and reused by any step in the pipeline (accessed by a cli script)
* `scripts` directory holds the scripts that can be easily called by snakemake (`Rscript example.R --input_ex example_input`)
* `snakemake` directory: workflow `.smk` files, **`snakemake/profiles/`** (Snakemake `--profile`: local vs Slurm/HPC defaults), **`input_templates/`**, and shared **`util/`** code between pipelines
* `docker` directory holds the information for creating the docker image that the pipeline runs
* `tests` directory holds all R tests, and end to end pipeline tests

### Platform Setup

**Running on macOS, Linux, Slurm, or PBS?** See **[PLATFORM_SETUP.md](PLATFORM_SETUP.md)** for platform-specific setup (Apptainer/Lima, SIF cache, Snakemake profiles under **`snakemake/profiles/`**, `qsub` template).

### Making changes

See **[CONTRIBUTING.md](CONTRIBUTING.md)** for development setup, Docker rebuilds, unit tests, and end-to-end tests. You can also contact open an issue in this repo.
