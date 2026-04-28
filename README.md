# GeneHackman 

![CI Tests](https://github.com/MRCIEU/GeneHackman/actions/workflows/main.yml/badge.svg)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10624713.svg)](https://doi.org/10.5281/zenodo.10624713)

A pipeline for performing common genetic epidemiology tasks at the University of Bristol.

Goals:
* Implement comment steps in GWAS investigations to create reproducible, more efficient research
* Reusable pipelines that can be utilised on different projects
* Shared code and steps that can be updated according to the latest knowledge and practices

## Available Pipelines

There are **six** Snakemake pipelines (grouped as two tables of three). Each pipeline is a `.smk` file under `snakemake/`; see [PIPELINES.md](snakemake/PIPELINES.md) for JSON inputs and parameters.

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
| Volcano plot of MR results; **BF–BF coloc** for exposures that pass MR FDR | Finemap-only: no MR or coloc between traits (use when you only need SuSiE outputs) | Full coloc table + **HTML report** (`result_coloc.html`), including a disclaimer when ancestries differ between GWASes |
| Requires **QTL_DATA_DIR** (see `.env_example`) for QTL files              | Each GWAS must declare **ancestry** (for LD)                                | Configurable `finemap` and `coloc` priors/overlap; see [PIPELINES.md](snakemake/PIPELINES.md) |

## Onboarding

**Running on macOS, Linux, Slurm, or PBS?** See **[PLATFORM_SETUP.md](PLATFORM_SETUP.md)** for platform-specific setup (Apptainer/Lima, SIF cache, HPC profiles, `qsub` template).

### 1. Clone the repository into your personal space on BlueCrystal 4
`git clone git@github.com:MRCIEU/GeneHackman.git && cd GeneHackman`

### 2. Ensure you [have conda installed and initialised before activating](https://www.acrc.bris.ac.uk/protected/hpc-docs/software/python_conda.html)

```conda env create --file environment.yml```

or if you have already created the environment

```conda activate genehackman```

### 3. Populate .env and input.json files

`cp .env_example .env`
* populate the DATA_DIR, RESULTS_DIR and RDFS_DIR environment variables in .env file
These should probably be in your *work* or *scratch* space (`/user/work/{userid}/...`)
* PIPELINE_DATA_DIR: data that can be found here
* QTL_DATA_DIR (optional): if you want to run the QTL MR pipieline
* RDFS_DIR is optional.  All generated files can be copied automatically.  Please ensure the path
ends in `working/`
* **Container cache:** If there is no pre-built `genehackman_<version>.sif` under `PIPELINE_DATA_DIR`, Snakemake pulls the `docker://` image and caches the SIF under `.snakemake/singularity` by default. Set **`GENEHACKMAN_SINGULARITY_PREFIX`** in `.env` to use another directory (e.g. scratch). Running `snakemake` without `./run_pipeline.sh`? Pass `--singularity-prefix /path` or add `singularity-prefix:` to your profile `config.yaml`.

### Fill out input.json file
* Ex: `cp snakemake/input_templates/compare_gwases.json input.json`
* Each pipeline (as defined in `snakemake` directory) has its own input format.
  * [Here are example pipelines here, copy to input.json](snakemake/input_templates/)
  * [Documentation per pipeline](snakemake/PIPELINES.md)
* You can either copy into input.json, or supply the file into the script from another location

### 4. Run the pipeline

`./run_pipeline.sh snakemake/<specific_pipeline>.smk <optional_input_file.json>`

* By default `run_pipeline.sh` uses **local Docker** (`snakemake/local/`). On HPC, set e.g. `GENEHACKMAN_PROFILE=snakemake/bp1/` (or `bc4/`, `slurm_singularity/`).
* `run_pipeline.sh` is just a convience wrapper around the `snakemake` command, if you want to do anything out of the ordinary, [please read up on snakemake](https://snakemake.readthedocs.io/en/v7.26.0/)
* If there are errors while running the pipeline, you can find error messages either directly on the screen, or in slurm log file that is outputted on error
* It is recommended that you run the your pipeline [inside a tmux session](https://github.com/MRCIEU/GeneHackman/wiki/Common-Errors#ssh-disconnection-while-pipeline-is-running).

## How it works:

The standard column naming for GWASes are:

| CHR | BP  | EA  | OA  | BETA | SE  | P   | EAF | SNP | RSID |
|-----|-----|-----|-----|------|-----|-----|-----|-----|:-----|

A full list of names and default values [can be found here](inst/extdata/predefined_column_maps.csv)

There are 3 main components to the pipeline
1. Snakemake to define the steps to complete for each pipeline.
2. Docker / Singularity container with installed languages (R and python), packages, os libraries, and code
3. Slurm: each snakemake step spins up a singularity container inside a slurm job.  Each step can specify different slurm requirements.

### Repository Organisation

* `R` directory holds R package code that can also be called and reused by any step in the pipeline (accessed by a cli script)
* `scripts` directory holds the scripts that can be easily called by snakemake (`Rscript example.R --input_ex example_input`)
* `snakemake` directory, which defines the pipeline steps and configuration, and shared code between pipelines
* `docker` directory holds the information for creating the docker image that the pipeline runs
* `tests` directory holds all R tests, and a end to end pipeline test script 

### Docker image

Pipelines use the image **`mrcieu/genehackman`** (tag usually matches the **`Version:`** field in [`DESCRIPTION`](DESCRIPTION) and **`DOCKER_VERSION`** in `.env`; see [`snakemake/util/common.smk`](snakemake/util/common.smk)).

**Build** (run from the **repository root** so paths in the Dockerfile resolve correctly):

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t mrcieu/genehackman:$(grep '^Version:' DESCRIPTION | awk '{print $2}') .
```

To tag a specific version explicitly (for example `1.0.0`):

```bash
docker build --platform linux/amd64 -f docker/Dockerfile -t mrcieu/genehackman:1.0.0 .
```

The image is **x86_64/amd64-only** (Posit R `.deb`, PLINK, Miniconda, etc.). Use **`--platform linux/amd64`** when building so it works on **ARM** laptops (e.g. Apple Silicon) as well as on amd64 Linux; on native amd64 machines the flag is optional but harmless.

The image installs OS packages and tools in [`docker/Dockerfile`](docker/Dockerfile), R packages via [`docker/requirements.R`](docker/requirements.R), and Python packages via [`docker/requirements.txt`](docker/requirements.txt). After changing those files, rebuild the image before running pipelines that depend on the new dependencies.

**Publish** (maintainers with access to the `mrcieu` organisation on Docker Hub):

```bash
docker push mrcieu/genehackman:<tag>
```

Use the same `<tag>` you built locally. HPC users without Docker can still obtain the image via Apptainer/Singularity (for example `singularity build genehackman_<tag>.sif docker://mrcieu/genehackman:<tag>`), as in [`run_pipeline.sh`](run_pipeline.sh).

### Making changes

If you want to make any additions / changes please contact andrew.elmore@bristol.ac.uk, or open an issue in this repo.
