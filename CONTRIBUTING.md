# Contributing to GeneHackman

This guide covers how to change code, update the Docker image, and run the test suite before opening a pull request.

For running pipelines in production, see [README.md](README.md) and [PLATFORM_SETUP.md](PLATFORM_SETUP.md). For pipeline inputs and parameters, see [PIPELINES.md](PIPELINES.md).

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
   - **`PIPELINE_DATA_DIR`** — absolute path to the reference data bundle from `gs://genehackman:1.1.0` (1000 Genomes LD panels, LDSC assets, etc.).

   **`DOCKER_VERSION`** is optional; it defaults to **`Version:`** in [`DESCRIPTION`](DESCRIPTION). Set it in `.env` only when you need a different image tag (e.g. `develop`).

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

## Adding a new pipeline

Use this checklist when you add a new top-level workflow under `snakemake/` (a new `.smk` that runs a distinct analysis end-to-end).

### 1. Create the workflow file

Add `snakemake/<pipeline_name>.smk` and `snakemake/<pipeline_name.md>` documentation at the repo root of the Snakemake tree (not under `rules/`). Follow the pattern used by existing workflows:

```python
include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "my_pipeline"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)  # if the workflow clumps

onstart:
    print("##### My Pipeline #####")

rule all:
    input: ...  # every final output Snakemake must build

include: "rules/standardise_rule.smk"   # reuse where appropriate
# include: "rules/clumping_rule.smk"
# include: "rules/finemap_rule.smk"

onsuccess:
    onsuccess(pipeline_name, files_created, results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
```

**Shared building blocks**

| Include | When to use |
|---------|-------------|
| `rules/standardise_rule.smk` | Almost always — harmonises each GWAS to `data/gwas/<prefix>_std.tsv.gz`. |
| `rules/clumping_rule.smk` | When PLINK clumping is required (`pipeline_includes_clumping=True` in `parse_pipeline_input`). |
| `rules/finemap_rule.smk` / `rules/finemap_multi_ancestry_rule.smk` | When SuSiE or MultiSuSiE fine-mapping is part of the workflow. |

Put logic that might be reused across workflows in `snakemake/rules/`. Keep pipeline-specific rules in the main `.smk` or a dedicated `rules/<pipeline>_rule.smk` included from there.

Call **`parse_pipeline_input()`** early. It loads the YAML path from `GENEHACKMAN_INPUT` / `--config genehackman_input=…`, validates `.env`, and attaches per-GWAS fields (`prefix`, `standardised_gwas`, `clumped_file`, column maps, etc.) on `pipeline.gwases`.

### 2. Add CLI entry points

Snakemake rules should call thin wrappers, not inline R/Python:

- **R:** add `scripts/my_step.R` that `source("load.R")`, parses args with **argparser**, and calls a function in `R/`.
- **Python:** add `scripts/my_step.py` and list any new packages in `docker/requirements.txt`.

Export new R functions from the package (`NAMESPACE`) and run `devtools::document()` when you add roxygen.

### 3. Define inputs and outputs

**YAML input**

- Add `snakemake/input_templates/<pipeline>.yaml` with sensible defaults and comments.
- Add a tiny fixture under `tests/testthat/data/snakemake_inputs/` for e2e runs (`is_test: true` is fine).
- If the pipeline needs new root-level YAML keys, extend `parse_pipeline_input()` in `snakemake/util/common.smk` (defaults, validation, and error messages belong there).

**Outputs**

- Write under **`PROJECT_DIR/results/`** (`RESULTS_DIR`) or **`PROJECT_DIR/data/`** (`DATA_DIR`) — do not hard-code user-specific paths.
- Register every deliverable in `rule all` so Snakemake knows when the run is complete.
- For **completion sentinel files** (`*_complete*.txt`), name them after the GWAS run (see `gwas_run_label()`, `FINEMAP_COMPLETE_TXT_PATTERN`, and `multi_finemap_complete_file()` in `snakemake/util/common.smk`). A generic `finemap_complete.txt` in a shared folder will block reruns when users reuse the same results directory for different inputs.

**Wildcards**

- Per-GWAS outputs usually key off `wildcards.prefix`, set from `file_prefix(g.file)` during YAML parsing.
- Use helpers in `common.smk` (`standardised_gwas_name()`, etc.) rather than duplicating path logic.

### 4. Document the pipeline

1. Add **`snakemake/<pipeline_name>.md`** next to the `.smk` with **Input** and **Output** sections (see existing files such as [`snakemake/finemap.md`](snakemake/finemap.md)).
2. Add a row to the pipeline table in [`PIPELINES.md`](PIPELINES.md) linking to the new doc.
3. Optionally add a one-line summary to the pipeline tables in [`README.md`](README.md).

### 5. Test

**Unit tests** — mock external tools (PLINK, SuSiE, liftover) and test R/Python logic in `tests/testthat/`.

**Dry run**

```bash
./run_pipeline.sh snakemake/my_pipeline.smk tests/testthat/data/snakemake_inputs/my_pipeline.yaml -n
```

**End-to-end** — append a line to [`tests/e2e_tests/run_test_pipelines.sh`](tests/e2e_tests/run_test_pipelines.sh):

```bash
./run_pipeline.sh snakemake/my_pipeline.smk tests/testthat/data/snakemake_inputs/my_pipeline.yaml -F
```

Run the full e2e script before opening a PR and commit the updated **`tests/testing_complete.txt`**.

### 6. Review checklist before opening a PR

- [ ] `rule all` lists every required output; no orphan rules.
- [ ] Example YAML and [`PIPELINES.md`](PIPELINES.md) / `snakemake/<pipeline>.md` updated.
- [ ] New R exports documented; `devtools::test()` passes.
- [ ] E2e entry added (unless the pipeline needs data you cannot ship in the repo — document why).
- [ ] New Python deps added to `docker/requirements.txt`; note in the PR if a new Docker image is required.
- [ ] Completion markers and other Snakemake targets are **run-specific** when outputs share a directory across analyses.

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

Bump **`Version:`** in `DESCRIPTION` when releasing; the pipeline defaults to that tag for the SIF name (`genehackman_<version>.sif`) and Docker pull. Set **`DOCKER_VERSION`** in `.env` only to override (e.g. `develop`).

HPC users without Docker pull the same image via `run_pipeline.sh`, which builds or uses `$PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<version>.sif`.

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

Users on release **`1.2.0`** get image tag **`1.2.0`** from **`Version:`** in `DESCRIPTION` by default; Snakemake looks for `genehackman_1.2.0.sif` under `PIPELINE_DATA_DIR/genomic_data/pipeline/`. Override with **`DOCKER_VERSION=1.2.0`** (or another tag) in `.env` if needed.

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

Edit **`Version:`** in [`DESCRIPTION`](DESCRIPTION) to the new semver (e.g. `1.2.0`). The pipeline and `run_pipeline.sh` use that value for the Docker/Apptainer image tag unless **`DOCKER_VERSION`** is set in `.env`.

Optionally document an override example in [`.env_example`](.env_example):

```bash
# DOCKER_VERSION=1.2.0
```

Regenerate R docs if exports changed:

```bash
Rscript -e "devtools::document()"
```

Commit on `main` (or via PR):

```bash
git add DESCRIPTION
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

1. Pull the new release (or check out a tag whose **`DESCRIPTION`** `Version:` matches the image you want). Set **`DOCKER_VERSION`** in `.env` only if you need a tag other than that default.
2. Pull or build the SIF, e.g. delete an old `genehackman_*.sif` and re-run `run_pipeline.sh` (it builds from `docker://mrcieu/genehackman:<version>` if the file is missing), or on HPC:

   ```bash
   singularity build genehackman_1.2.0.sif docker://mrcieu/genehackman:1.2.0
   ```

## Getting help

- Open a [GitHub issue](https://github.com/MRCIEU/GeneHackman/issues) for bugs or feature requests.
- Contact **andrew.elmore at bristol dot ac uk** for Bristol-internal coordination.
