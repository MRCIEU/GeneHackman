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

### 3. Get reference data (Google Cloud Storage)

Much of what the pipelines expect (LD references, LDSC masks, MetaBrain/QTL summaries, …) ships in two buckets. Copy them somewhere on your machine or HPC scratch space (**do not commit data to git**).

| Bucket | Typical use | Env var to set after download |
|--------|--------------|-------------------------------|
| Mandatory: `gs://genehackman` | Core pipeline bundle (e.g. 1000 Genomes LD panels, genomic assets, LDSC helpers under the layout expected under `PIPELINE_DATA_DIR`) | **`PIPELINE_DATA_DIR`** → absolute path of the folder that mirrors this bucket root (same internal directory names as on GCS). |
| Optional: `gs://genehackman-qtl` | QTL summary statistic trees used by **`qtl_mr`** | **`QTL_DATA_DIR`** → absolute path of the folder containing that hierarchy. Needed for **`qtl_mr`**; ignored by other pipelines. |

**Authenticate** (private buckets — use whichever applies): [install Google Cloud SDK](https://cloud.google.com/sdk/docs/install), then `gcloud auth login` so `gsutil` can read the buckets you’ve been granted.

**Download with `gsutil`** (multipart copy; adjust local paths):

```bash
mkdir -p /path/to/my_pipeline_data /path/to/my_qtl_data
gsutil -m rsync -r gs://genehackman/ /path/to/my_pipeline_data/
gsutil -m rsync -r gs://genehackman-qtl/ /path/to/my_qtl_data/
```

To copy only selected prefixes instead of the QTL bucket, use `gsutil -m cp -r gs://genehackman-qtl/SOME_PREFIX/ ...` as needed.  For example, you may only be interested in cis, not trans data.

**Download via web UI:** Open [Google Cloud Console → Cloud Storage](https://console.cloud.google.com/storage/browser), select the project/storage account you were given access to, open bucket **`genehackman`** or **`genehackman-qtl`**, and download files or folders in the browser (or use **“Activate Cloud Shell”** and run `gsutil` there, then drag files out if convenient).

Then point your **`.env`** at those directories (trailing slashes are fine):

```bash
PIPELINE_DATA_DIR=/path/to/my_pipeline_data/
QTL_DATA_DIR=/path/to/my_qtl_data/
```

Alternatives if you don’t set **`QTL_DATA_DIR`** directly: keep the same directory layout under **`PIPELINE_DATA_DIR/qtl_datasets/`**, or follow [PLATFORM_SETUP.md](PLATFORM_SETUP.md) (bind mounts and defaulting **`QTL_DATA_DIR`** to **`PIPELINE_DATA_DIR/qtl_datasets`**).

### 4. Populate `.env` and `input.yaml` files

`cp .env_example .env`
* Populate **`DATA_DIR`** and **`RESULTS_DIR`** — usually under *work* or *scratch* (e.g. `/user/work/{userid}/...`).
* Set **`PIPELINE_DATA_DIR`** to the path where you unpacked **`genehackman`** (see §3).
* Set **`QTL_DATA_DIR`** to the path where you unpacked **`genehackman-qtl`** if you run **`qtl_mr`** (can be left empty otherwise; see [.env_example](.env_example)).
* **Container cache:** If there is no pre-built `genehackman_<version>.sif` under `PIPELINE_DATA_DIR`, Snakemake pulls the `docker://` image and caches the SIF under `.snakemake/singularity` by default. Set **`GENEHACKMAN_SINGULARITY_PREFIX`** in `.env` to use another directory (e.g. scratch). Running `snakemake` without `./run_pipeline.sh`? Pass `--singularity-prefix /path` or add `singularity-prefix:` to your profile `config.yaml`.

**`input.yaml`**

* Example: `cp snakemake/input_templates/compare_gwases.yaml input.yaml`
* Each pipeline has its own shape; examples live under [`snakemake/input_templates/`](snakemake/input_templates/).
* See [PIPELINES.md](snakemake/PIPELINES.md) for all fields.
* Pass the input YAML as the **second** argument to `run_pipeline.sh`, or rely on **`input.yaml`** in the working directory. To call **`snakemake`** yourself without the wrapper, use **`--config genehackman_input=/path/to/file.yaml`**.

### 5. Run the pipeline

`./run_pipeline.sh snakemake/<specific_pipeline>.smk <optional_input_file.yaml>`

* By default `run_pipeline.sh` uses **local Docker** (`snakemake/local/`). On HPC, set e.g. `GENEHACKMAN_PROFILE=snakemake/bp1/` (or `bc4/`, `slurm_singularity/`).
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
