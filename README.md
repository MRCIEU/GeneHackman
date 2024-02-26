# gwaspipeline

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10624713.svg)](https://doi.org/10.5281/zenodo.10624713)

A pipeline for creating reusable steps for genetic epidemiology projects at the University of Bristol.

Goals:
* Reproducible more efficient research 
* Reusable pipelines that can be utilised on different projects
* Shared code and steps that can be updated according to the latest knowledge and practices

## Onboarding

Steps to start using the pipeline
### 1. Clone the repository into your personal space on BlueCrystal 4
`git clone git@github.com:MRCIEU/gwaspipeline.git && cd gwaspipeline`

### 2. Create and activate conda environment
Ensure you [have conda installed and initialised](https://www.acrc.bris.ac.uk/protected/hpc-docs/software/python_conda.html) before running these commands:
```
conda env create --file environment.yml #only needed the first time
conda activate gwaspipeline
```
### 3. Populate .env and input.json files

### `cp .env_example .env`
* populate the DATA_DIR, RESULTS_DIR and RDFS_DIR environment variables in .env file
These should probably be in your *work* or *scratch* space (`/user/work/userid/...`)
* RDFS_DIR is optional.  All generated files can be copied automatically.  Please ensure the path
ends in `working/`

### Fill out input.json file
* Each pipeline (as defined in `snakemake` directory) has its own input format.  [There is documentation per pipeline here](snakemake/PIPELINES.md)
* With each GWAS, you can specify header names ex. `{"P":"your_gwas_pval_col", ...}`, if you do not specify header names it will assume your GWAS has the headers below.
* You can either copy into input.json, or supply the file into the script from another location

### 4. Run the pipeline
### `./run_pipeline.sh snakemake/<specific_pipeline>.smk <optional_input_file.json>`

If there are errors while running the pipeline, you can find error messages either directly on the screen, or there may be a slurm log file that is outputted on error, which you can look at for clues.

## How it works:

The standard column naming for GWASes are:

|           | CHR | BP  | EA  | OA  | BETA | SE  | P   | EAF | SNP | RSID |
|-----------|-----|-----|-----|-----|------|-----|-----|-----|-----|:-----|
| Mandatory | x   | x   | x   | x   |      |     |     |     |     |      |

* x denotes mandatory columns

There are 3 main components to the pipeline
1. Snakemake to define the steps to complete for each pipeline.
2. Docker / Singularity container with installed languages (R and python), packages, os libraries, and code
3. Slurm: each snakemake step spins up a singularity container inside a slurm job.  Each step can specify different slurm requirements.

## Repository Organisation

* `docker` directory holds the information for creating the docker image that the pipeline runs
* `scripts` directory holds the scripts that can be easily called by snakemake (`Rscript example.R --input_file example.txt`)
* `R` directory holds R package code that can also be called and reused by any step in the pipeline (accessed by a cli script)
* `snakemake` directory, which defines the pipeline steps and configuration, and shared code between pipelines

## Making changes

To make changes to the docker images, you will have to contact andrew.r.elmore@gmail.com, since the storage of the docker image is tied to his docker account at the moment.
