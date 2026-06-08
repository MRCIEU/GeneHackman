# Contributing to GeneHackman

This guide covers how to change code, update the Docker image, and run the test suite before opening a pull request.

For running pipelines in production, see [README.md](README.md) and [PLATFORM_SETUP.md](PLATFORM_SETUP.md). For pipeline inputs and parameters, see [snakemake/PIPELINES.md](snakemake/PIPELINES.md).

## Development setup

1. Clone the repository and create the conda environment:

   ```bash
   git clone git@github.com:MRCIEU/GeneHackman.git
   cd GeneHackman
   conda env create -f environment.yml
   conda activate genehackman
   ```

2. Copy and edit `.env` (see [.env_example](.env_example)):

   ```bash
   cp .env_example .env
   ```

   For development you need at least:

   - **`PROJECT_DIR`** — absolute path where test outputs go (e.g. a scratch folder; the pipeline writes to `PROJECT_DIR/data/` and `PROJECT_DIR/results/`).
   - **`PIPELINE_DATA_DIR`** — absolute path to the reference data bundle from `gs://genehackman` (1000 Genomes LD panels, LDSC assets, etc.).
   - **`DOCKER_VERSION`** — image tag to run (e.g. `1.1.0` or `develop`).

   Use **absolute paths** in `.env`. Relative paths (e.g. `QTL_DATA_DIR=hi`) break Apptainer bind mounts with errors like `destination must be an absolute path`.

3. Install the R package locally (for unit tests outside Docker):

   ```bash
   Rscript -e "devtools::install()"
   ```

## Repository layout

| Path | Role |
|------|------|
| `R/` | R package functions used by pipeline steps |
| `scripts/` | CLI entry points called from Snakemake (`Rscript …`, `python …`) |
| `snakemake/` | Workflow `.smk` files, `profiles/`, `input_templates/`, shared `util/` |
| `docker/` | `Dockerfile`, `requirements.R`, `requirements.txt` |
| `tests/testthat/` | Unit tests and small test GWAS files |
| `tests/e2e_tests/` | End-to-end Snakemake test runner |
| `inst/` | Package data (e.g. column maps) |

Snakemake profiles bind-mount `R/`, `scripts/`, and `inst/` from the repo into the container, so **changes to R and script code take effect without rebuilding the image** on the next pipeline run. New R or Python **dependencies** still require a Docker rebuild.

## Making code changes

### R package (`R/`)

- Follow existing style: **data.table** / **dplyr** patterns already in the file, **roxygen2** docs for exported functions.
- Regenerate documentation when you change exports:

  ```bash
  Rscript -e "devtools::document()"
  ```

- Snakemake rules typically call thin wrappers in `scripts/` that load the package and parse CLI args.

### Scripts (`scripts/`)

- Keep scripts as CLI wrappers; put reusable logic in `R/`.
- Python scripts (e.g. `run_multisusie.py`) should stay compatible with packages in `docker/requirements.txt`.

### Snakemake (`snakemake/`)

- Shared helpers live in `snakemake/util/` (`common.smk`, `constants.smk`, rules under `snakemake/rules/`).
- Add or update an example input under `snakemake/input_templates/` when you change required YAML fields.
- Site-specific cluster settings belong in new profiles under `snakemake/profiles/` (copy `local/` or `slurm/` as a template).

### Conventions

- **Ancestry** codes must be one of: `EUR`, `EAS`, `AFR`, `AMR`, `SAS`.
- **Finemap** (`finemap.smk`): ancestries must be either all the same (single-ancestry SuSiE) or all distinct (multi-ancestry MultiSuSiE). Mixed duplicates fail at startup.
- **Coloc** (`coloc.smk`): all GWAS inputs must share the same ancestry.

## Docker changes

The pipeline runs inside **`mrcieu/genehackman`** (Apptainer/Singularity on HPC, Docker locally).

| File | Purpose |
|------|---------|
| `docker/Dockerfile` | Base OS, R, PLINK, liftOver, LDSC, PHESANT, bcftools |
| `docker/requirements.R` | CRAN/Bioconductor R dependencies |
| `docker/requirements.txt` | Python dependencies (Snakemake, MultiSuSiE, …) |

The Dockerfile copies only `DESCRIPTION`, `docker/requirements.R`, and `docker/requirements.txt` before installing dependencies, so **edits to those files invalidate the dependency layer** and rebuild quickly without copying the whole repo first.

### Build locally

From the **repository root**:

```bash
docker build --platform linux/amd64 -f docker/Dockerfile \
  -t mrcieu/genehackman:$(grep '^Version:' DESCRIPTION | awk '{print $2}') .
```

The image is **linux/amd64 only**. Use `--platform linux/amd64` on Apple Silicon.

After changing `DESCRIPTION` (new R package in `Imports:`), update `docker/requirements.R` or rely on `remotes::install_deps("docker", …)` picking up new imports.

After changing Python deps, edit `docker/requirements.txt` and rebuild.

### Publish (maintainers)

```bash
docker push mrcieu/genehackman:<tag>
```

Bump **`Version:`** in `DESCRIPTION` and **`DOCKER_VERSION`** in `.env` together so the SIF name (`genehackman_<version>.sif`) matches the image tag.

HPC users without Docker pull the same image via `run_pipeline.sh`, which builds or uses `$PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<DOCKER_VERSION>.sif`.

## Unit tests

Unit tests use **testthat** and live in `tests/testthat/`.

### Run all package tests

```bash
# Inside the conda env, with the package installed:
Rscript -e "devtools::test()"

# Or:
Rscript tests/testthat.R
```

### Full package check (what CI runs)

Runs `R CMD check`–style validation (examples, vignettes, namespace, etc.):

```bash
Rscript -e "devtools::check()"
```

CI runs this inside `mrcieu/genehackman:develop` (see [.github/workflows/main.yml](.github/workflows/main.yml)).

### Writing tests

- Add new test files as `tests/testthat/test_<topic>.R`.
- Use `testthat::local_mocked_bindings()` to mock external tools (PLINK, SuSiE, liftover) where the existing tests do.
- Small GWAS fixtures are under `tests/testthat/data/`.

## End-to-end pipeline tests

End-to-end tests run real Snakemake workflows against tiny test GWAS files via Apptainer.

### Prerequisites

- `.env` configured with valid **`PROJECT_DIR`** and **`PIPELINE_DATA_DIR`** (reference data required for LD, liftover, etc.).
- Apptainer/Singularity available (see [PLATFORM_SETUP.md](PLATFORM_SETUP.md)).
- Conda env activated.

### Run all e2e tests

```bash
./tests/e2e_tests/run_test_pipelines.sh
```

This script runs `run_pipeline.sh` with `-F` (force rerun) for:

| Pipeline | Test input |
|----------|------------|
| `standardise_gwas.smk` | `tests/testthat/data/snakemake_inputs/standardise_gwas.yaml` |
| `disease_progression.smk` | `tests/testthat/data/snakemake_inputs/disease_progression.yaml` |
| `compare_gwases.smk` | `tests/testthat/data/snakemake_inputs/compare_gwases.yaml` |
| `finemap.smk` | `finemap.yaml` and `finemap_multi_ancestry.yaml` |
| `coloc.smk` | `tests/testthat/data/snakemake_inputs/coloc.yaml` |
| `qtl_mr.smk` | `qtl_mr_eqtlgen.yaml` (only if **`QTL_DATA_DIR`** is set in `.env`) |

On success it writes **`tests/testing_complete.txt`** with a line like:

```text
SUCCESS: All tests passed on branch: your-branch-name
```

### CI requirement

Pull requests must include an updated **`tests/testing_complete.txt`** from a successful run on **your branch**. GitHub Actions checks that:

1. The file exists.
2. On non-`main` branches, the file contains the branch name.

Run the e2e script on your feature branch, then commit `tests/testing_complete.txt` with your other changes.

### Run a single pipeline test

```bash
./run_pipeline.sh snakemake/finemap.smk \
  tests/testthat/data/snakemake_inputs/finemap.yaml -F
```

Useful flags: `--dry-run`, `--unlock`, `-n` (dry run), `-R <rule>` (rerun specific rule).

## Suggested workflow before a PR

1. Make changes on a feature branch.
2. Run unit tests: `Rscript -e "devtools::test()"` (or `devtools::check()` for a fuller pass).
3. Run e2e tests: `./tests/e2e_tests/run_test_pipelines.sh`.
4. Commit code changes **and** `tests/testing_complete.txt`.
5. Open a pull request against `main`.

If you change Docker dependencies, note the new image tag in the PR description and confirm you have rebuilt (or that maintainers will publish) the matching `mrcieu/genehackman` image.

## Cutting a release (maintainers)

Releases tie together three versioned artefacts:

| Artefact | Where | Format |
|----------|--------|--------|
| R package | `DESCRIPTION` → `Version:` | `1.2.0` (no `v` prefix) |
| Docker / Apptainer image | Docker Hub `mrcieu/genehackman` | tag `1.2.0` (matches `Version:`) |
| Git tag | GitHub | `v1.2.0` (`v` + same semver) |

Users set **`DOCKER_VERSION=1.2.0`** in `.env`; Snakemake looks for `genehackman_1.2.0.sif` under `PIPELINE_DATA_DIR/genomic_data/pipeline/`.

### Before you release

1. Merge all intended changes to **`main`**.
2. Confirm CI is green on `main` ([Actions](https://github.com/MRCIEU/GeneHackman/actions)).
3. Run the full test suite on `main`:

   ```bash
   git checkout main && git pull
   Rscript -e "devtools::check()"
   ./tests/e2e_tests/run_test_pipelines.sh
   ```

4. Commit `tests/testing_complete.txt` on `main` if the e2e run updated it.

### 1. Bump the version

Edit **`Version:`** in [`DESCRIPTION`](DESCRIPTION) to the new semver (e.g. `1.2.0`).

Update the example default in [`.env_example`](.env_example):

```bash
DOCKER_VERSION=1.2.0
```

Regenerate R docs if exports changed:

```bash
Rscript -e "devtools::document()"
```

Commit on `main` (or via PR):

```bash
git add DESCRIPTION .env_example
git commit -m "Bump version to 1.2.0"
git push origin main
```

### 2. Build and publish the Docker image

From the repository root, on a machine with Docker Hub access to **`mrcieu`**:

```bash
VERSION=1.2.0

docker build --platform linux/amd64 -f docker/Dockerfile \
  -t mrcieu/genehackman:${VERSION} .

docker push mrcieu/genehackman:${VERSION}
```

Optional: refresh the rolling **`develop`** tag used by CI (`mrcieu/genehackman:develop` in [`.github/workflows/main.yml`](.github/workflows/main.yml)):

```bash
docker tag mrcieu/genehackman:${VERSION} mrcieu/genehackman:develop
docker push mrcieu/genehackman:develop
```

### 3. Tag the release in Git

Create an annotated tag on `main` pointing at the version bump commit:

```bash
git checkout main && git pull
git tag -a v1.2.0 -m "Release 1.2.0"
git push origin v1.2.0
```

Tags use a **`v` prefix** (e.g. `v1.0.0`); Docker tags do **not** (`1.2.0`).

### 4. Create the GitHub release

Using the GitHub CLI:

```bash
gh release create v1.2.0 \
  --title "1.2.0" \
  --notes "$(cat <<'EOF'
## Summary
- …

## Docker
`docker pull mrcieu/genehackman:1.2.0`

## Citation
https://doi.org/10.5281/zenodo.10624713
EOF
)"
```

Or in the browser: **GitHub → Releases → Draft a new release** → choose tag `v1.2.0`, title `1.2.0`, and add release notes (changes since the previous tag, Docker pull command, any breaking changes).

### 5. Zenodo archive

The project is archived on Zenodo ([10.5281/zenodo.10624713](https://doi.org/10.5281/zenodo.10624713)). If the [Zenodo–GitHub integration](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content) is enabled for this repository, publishing the GitHub release should trigger a new Zenodo version automatically. Otherwise, upload the release manually on Zenodo and note the new version DOI in the GitHub release.

### After release

Tell users to:

1. Set **`DOCKER_VERSION`** in `.env` to the new version.
2. Pull or build the SIF, e.g. delete an old `genehackman_*.sif` and re-run `run_pipeline.sh` (it builds from `docker://mrcieu/genehackman:<version>` if the file is missing), or on HPC:

   ```bash
   singularity build genehackman_1.2.0.sif docker://mrcieu/genehackman:1.2.0
   ```

## Getting help

- Open a [GitHub issue](https://github.com/MRCIEU/GeneHackman/issues) for bugs or feature requests.
- Contact **andrew.elmore at bristol dot ac uk** for Bristol-internal coordination.
