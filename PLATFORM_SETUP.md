# Running GeneHackman by platform

Snakemake orchestrates steps; container execution is selected by **`SNAKEMAKE_PROFILE`**. Local runs can use either **Docker** (`snakemake/profiles/docker/`) or **Apptainer/Singularity** (`snakemake/profiles/apptainer/`). Cluster runs use Apptainer/Singularity profiles such as `snakemake/profiles/slurm/`. Profiles live under **`snakemake/profiles/`** and define how jobs run plus which host paths are bind-mounted into the container.

General prerequisites:

1. **Conda env:** `conda env create -f environment.yml` then `conda activate genehackman`. The conda env provides Snakemake for Apptainer profiles and the wrapper utilities for all profiles.
2. **Reference data:** download the current release data from **`gs://genehackman/1.1.0/`** into **`PIPELINE_DATA_DIR`** so `genomic_data/` and `LDSCORE/` sit directly under that directory:

   ```bash
   gsutil -m rsync -r gs://genehackman/1.1.0/ /path/to/my_pipeline_data/
   ```

   If a future release updates code and reference data together, use the matching bucket prefix (for example `gs://genehackman/1.2.0/`) and update the docs accordingly.

3. **`.env`:** copy from `.env_example` and set **`PROJECT_DIR`** (the pipeline uses **`PROJECT_DIR/data/`** and **`PROJECT_DIR/results/`**) and **`PIPELINE_DATA_DIR`** (path from step 2). Optionally set **`QTL_DATA_DIR`** when large QTL mirrors live on another volume or object-store mount; if omitted, `run_pipeline.sh` defaults it to `PIPELINE_DATA_DIR/qtl_datasets`.
4. **Input YAML:** see `snakemake/input_templates/` and [PIPELINES.md](PIPELINES.md).

**`./run_pipeline.sh`** loads `.env`, picks a **profile** (**`SNAKEMAKE_PROFILE`**, default `snakemake/profiles/apptainer/`), then runs Snakemake. Pass the **`.smk` workflow first**, then optional input YAML, for example:

```bash
./run_pipeline.sh snakemake/standardise_gwas.smk path/to/input.yaml
```

---

## Choosing a Runtime

| Profile | Runtime | Use case |
|---------|---------|----------|
| `snakemake/profiles/docker/` | Local Docker | Workstations where Docker is available and you do not want to build/use a SIF. |
| `snakemake/profiles/apptainer/` | Local Apptainer/Singularity | Workstations or servers with Apptainer/Singularity. This remains the default. |
| `snakemake/profiles/slurm/` | Slurm + Apptainer/Singularity | Cluster execution; no Docker-on-Slurm support is provided. |
| `snakemake/profiles/uob-bp1/` | Site-specific Slurm + Apptainer/Singularity | Example site profile; edit for your own cluster. |

Set the profile in `.env`:

```bash
SNAKEMAKE_PROFILE=snakemake/profiles/docker/
```

`run_pipeline.sh` only builds or uses `genehackman_<version>.sif` for Apptainer/Singularity profiles. With the Docker profile it runs Snakemake inside `mrcieu/genehackman:<version>` directly.

---

## Local Docker

Use this when Docker is installed and you want a standalone local run without Apptainer/Singularity:

```bash
SNAKEMAKE_PROFILE=snakemake/profiles/docker/
./run_pipeline.sh snakemake/compare_gwases.smk path/to/input.yaml
```

The wrapper mounts the repository, `PROJECT_DIR` outputs, `PIPELINE_DATA_DIR`, optional `QTL_DATA_DIR`, and then runs Snakemake inside `mrcieu/genehackman:<version>`. It uses **`DOCKER_VERSION`** if set, otherwise **`Version:`** in `DESCRIPTION`.

On Apple Silicon, the wrapper defaults to `--platform linux/amd64` because the published image is amd64. Override with **`GENEHACKMAN_DOCKER_PLATFORM`** only if a multi-architecture image is published later.

This is local-only Docker support. The Slurm profiles continue to use Apptainer/Singularity.

---

## macOS

### Apptainer and Snakemake

Snakemake looks for an executable named **`singularity`** on `PATH`. Many installs only provide **`apptainer`**. Shell **aliases** in `~/.zshrc` (e.g. `alias singularity='limactl shell apptainer'`) are **not** visible to `run_pipeline.sh` because it uses **bash** and does not load your zsh aliases.

Typical setups:

- **Lima VM** with Apptainer (e.g. `limactl shell apptainer`) — run the pipeline from an environment where **`apptainer`** or **`singularity`** is a real binary on `PATH`, or use a small **`~/bin/singularity`** wrapper that calls `limactl shell <vm> -- apptainer "$@"`.
- **Homebrew Apptainer** (if available for your OS version).

Customising Lima VM for macOS
```
mountType: virtiofs
mounts:
- location: ~
  mountPoint: ~
-location /tmp/lima
  mountPOint: /tmp/lima
```


### SIF cache and `/Users` mounts

If the workflow pulls a `docker://` image, Snakemake caches the `.sif` under **`--singularity-prefix`** (default: `.snakemake/singularity` in the repo). On **macOS + Lima**, paths under **`/Users/...`** are often VM mounts where **creating large SIF files fails** with errors like **read-only file system**.

- Set **`SINGULARITY_DIR`** in `.env` to a **writable** location **inside the Linux VM** (commonly a path under that VM’s **`/tmp`**, e.g. `/tmp/genehackman_snakemake_singularity`), or to native VM disk — not only a bind-mounted folder that cannot handle the build.
- Alternatively, place a pre-built **`$PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<version>.sif`** (underscore; `<version>` defaults to **`Version:`** in `DESCRIPTION`, overridable via **`DOCKER_VERSION`**) so Snakemake does not need to build on a problematic mount.
- Avoid **colons** in host-side SIF filenames (`genehackman:1.1.0.sif`); macOS can reject them. The repo’s Snakemake helper uses **`genehackman_<version>.sif`**.

### Apple Silicon (arm64)

The published image may be **linux/amd64**. Use Lima with **Rosetta / x86 Linux** (or run on an amd64 Linux machine / HPC) if pulls or runs fail with architecture warnings.

---

## Linux (workstation or single server)

For Docker, use `snakemake/profiles/docker/` as described above.

For Apptainer/Singularity:

1. Install **Apptainer** or **SingularityCE** (distribution packages or upstream instructions).
2. Ensure **`singularity`** or **`apptainer`** is on `PATH` (some sites symlink `singularity` → `apptainer`).
3. Use **`snakemake/profiles/apptainer/`** for local execution (same as default in `run_pipeline.sh`):

   ```bash
   export SNAKEMAKE_PROFILE=snakemake/profiles/apptainer/
   ./run_pipeline.sh snakemake/compare_gwases.smk
   ```

4. Adjust **`snakemake/profiles/apptainer/config.yaml`** `singularity-args` if your data live outside **`PROJECT_DIR`** (or **`PIPELINE_DATA_DIR`** / **`QTL_DATA_DIR`**); profiles bind **`DATA_DIR`** and **`RESULTS_DIR`** (derived from **`PROJECT_DIR`**). Add more `-B host:host` pairs if needed.

---

## Slurm (HPC)

The repo ships Slurm-oriented profiles (paths and partitions are **Bristol / MRC IEU-oriented**; **edit** them for your site):

| Profile directory | Notes |
|-------------------|--------|
| `snakemake/profiles/slurm/` | Slurm + Apptainer: make sure Slurm `account` and `partition` are set correctly for your HPC |
| `snakemake/profiles/apptainer/` | Runs Apptainer locally, no job invocation |
| `snakemake/profiles/docker/` | Runs Docker locally, no job invocation and no SIF build |
| `snakemake/profiles/uob-bp1/` | Specific example for UoB BP1 HPC |

**Typical usage:**

```bash
export SNAKEMAKE_PROFILE=snakemake/profiles/slurm/
./run_pipeline.sh snakemake/compare_gwases.smk
```

`run_pipeline.sh` runs **`module load ${APPTAINER_MODULE}`** for non-local Apptainer profiles — set **`APPTAINER_MODULE`** in `.env` to match your cluster. The Docker profile is local only and does not load Apptainer modules.

**Customize** each profile under **`snakemake/profiles/`** (each directory has a **`config.yaml`**) as needed:

- **`Slurm account / partition`:** In **`.env`**, optional **`SLURM_ACCOUNT`** overrides **`sbatch --account`**; if unset, the profile uses **`sacctmgr … | grep … | head -n1`**. **`SLURM_PARTITION`**, if unset, is inferred by **`run_pipeline.sh`** from **`sinfo -h -o '%P'`** (the partition marked with **`*`**, e.g. **`compute*`** → **`compute`**); if **`sinfo`** is unavailable or returns nothing, the fallback is **`compute`**. Passed to Snakemake as **`--default-resources partition=…`**.
- **`cluster:`** block: `sbatch` options (`partition`, `account`, walltime, memory, `--output` log directory).
- **`singularity-args`:** The generic profile binds repo code, **`$HOME`**, **`DATA_DIR`** and **`RESULTS_DIR`** (paths under **`PROJECT_DIR`**, set by `run_pipeline.sh`), **`PIPELINE_DATA_DIR`**, **`/tmp`**, and adds **`QTL_DATA_DIR`** only when non-empty (see **`GENEHACKMAN_EXTRA_SINGULARITY_BINDS`** in `run_pipeline.sh`). Add more `-B host:host` pairs for shared reference data on your site.

**Snakemake version:** `environment.yml` pins **Snakemake 7.x**. Migrating to **Snakemake 8+** changes cluster syntax (`cluster-generic` executor + plugins); this repository’s profiles target the v7 style unless you have already updated them.

---

## PBS Pro / OpenPBS / Torque

There is **no** checked-in PBS profile. You can add one by mirroring the Slurm profiles with a **`cluster:`** line that calls your scheduler, for example (Snakemake 7 — adjust directives to your site):

```yaml
cluster: >-
  qsub -N {rule}-{wildcards}
  -l walltime={resources.time}
  -l select=1:ncpus={threads}:mem={resources.mem}
  -o /path/to/your/pbs_logs/{rule}_%I.out
  -j oe
```

Requirements:

- **`qsub`** must return the job ID on stdout in a form Snakemake can track.
- For large workflows, set **`--cluster-status`** (or Snakemake’s DRMAA integration) if your site documents it; otherwise use Snakemake’s polling defaults where supported.

Copy `snakemake/slurm_singularity/config.yaml` as a template, replace the **`cluster:`** section with **`qsub`**, and add **`singularity-args`** appropriate for PBS compute nodes (shared filesystem mounts).

## Environment variables (quick reference)

| Variable | Purpose |
|----------|---------|
| **`SNAKEMAKE_PROFILE`** | Snakemake profile path (default `snakemake/profiles/apptainer/`). Local Docker: `snakemake/profiles/docker/`. Example HPC: `snakemake/profiles/slurm/` |
| **`SLURM_ACCOUNT`** | Optional (**`profiles/slurm`**). Sets **`sbatch --account`**; omitted ⇒ existing **`sacctmgr`** derivation in `config.yaml` |
| **`SLURM_PARTITION`** | Optional (**`profiles/slurm`**). Overrides **`sinfo`** default-partition detection in **`run_pipeline.sh`** (see **`sinfo`** `*` suffix); ultimate fallback **`compute`** |

The **pipeline YAML path** is not configured via `.env`. Use **`./run_pipeline.sh <workflow>.smk [path/to/input.yaml]`** (defaults to **`input.yaml`**) or run **`snakemake`** with **`--config genehackman_input=path/to/input.yaml`** (see [PIPELINES.md](PIPELINES.md)).

Use **`.env`** for values you want loaded every time `./run_pipeline.sh` runs (`export $(cat .env | xargs)`).

---

## Further reading

- Pipeline input schema: [PIPELINES.md](PIPELINES.md)
- Main readme: [`README.md`](README.md)
- Snakemake profiles: [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
- Snakemake cluster execution: [docs](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
