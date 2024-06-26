# GeneHackman 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10624713.svg)](https://doi.org/10.5281/zenodo.10624713)

A pipeline for performing common genetic epidemiology tasks at the University of Bristol.

Goals:
* Implement comment steps in GWAS investigations to create reproducible, more efficient research
* Reusable pipelines that can be utilised on different projects
* Shared code and steps that can be updated according to the latest knowledge and practices

## Available Pipelines

Here is a list of available pipelines, and the steps they run

| standardise_gwas.smk                                    | compare_gwases.smk                     | disease_progression.smk                               | qtl_mr.smk                                      |
|---------------------------------------------------------|----------------------------------------|-------------------------------------------------------|-------------------------------------------------|
| Takes in any of: vcf, csv, tsv, txt (and zip, gz)       | All steps in 'standardise_gwas.smk'    | All steps in 'standardise_gwas.smk'                   | All steps in 'standardise_gwas.smk'             |
| Convert Reference Build<br/>(GRCh38 -> GRCh37)          | PLINK clumping                         | Runs some collider bias corrections, compares results | Run MR against top hits of specific QTL dataset |
| Populate RSID from CHR, BP, EA, and OA                  | Calculate heterogeneity between GWASes | Miami Plot of Collider Bias Results                   | Volcano Plot of Results                         |
| Converts z-scores and OR to BETA                        | LDSC h2 and rg                         | Expected vs. Observed Comparison                      | Run coloc of significant top hit MR results     |
| Auto-populate GENE ID <-> ENSEMBL ID                    | Expected vs. Observed Comparison       |                                                       |                                                 |
| Unique SNP = CHR:BP_EA_OA<br/> (EA < OA alphabetically) |                                        |                                                       |                                                 |

## Onboarding

### 1. Clone the repository into your personal space on BlueCrystal 4
`git clone git@github.com:MRCIEU/GeneHackman.git && cd GeneHackman`

### 2. Ensure you [have conda installed and initialised before activating](https://www.acrc.bris.ac.uk/protected/hpc-docs/software/python_conda.html)

`conda activate /mnt/storage/private/mrcieu/data/genomic_data/pipeline/genehackman`

(you can alternatively create your own conda environment if you like: `conda env create --file environment.yml`)

### 3. Populate .env and input.json files

### `cp .env_example .env`
* populate the DATA_DIR, RESULTS_DIR and RDFS_DIR environment variables in .env file
These should probably be in your *work* or *scratch* space (`/user/work/userid/...`)
* RDFS_DIR is optional.  All generated files can be copied automatically.  Please ensure the path
ends in `working/`

### Fill out input.json file
* Ex: `cp snakemake/input_templates/compare_gwases.json input.json`
* Each pipeline (as defined in `snakemake` directory) has its own input format.
  * [Here are example pipelines here, copy to input.json](snakemake/input_templates/)
  * [Documentation per pipeline](snakemake/PIPELINES.md)
* You can either copy into input.json, or supply the file into the script from another location

### 4. Run the pipeline

`./run_pipeline.sh snakemake/<specific_pipeline>.smk <optional_input_file.json>`

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

### Making changes

If you want to make any additions / changes please contact andrew.elmore@bristol.ac.uk, or open an issue in this repo.
