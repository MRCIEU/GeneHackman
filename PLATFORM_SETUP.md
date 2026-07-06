# Running GeneHackman by platform

GeneHackman is developed and tested on **Linux** and **HPC (Slurm)**. Those are the supported platforms.

Before following the platform setup, please follow the general prerequisites and setup process in the [README.md](README.md).
---

## Linux

Use this for a workstation or single Linux server.

1. Install **Apptainer** or **SingularityCE** (distribution packages or upstream instructions).
2. Ensure **`singularity`** or **`apptainer`** is on `PATH` (some sites symlink `singularity` → `apptainer`).
3. Set the local profile in `.env` (this is the default):

   ```bash
   SNAKEMAKE_PROFILE=snakemake/profiles/local/
   ./run_pipeline.sh snakemake/compare_gwases.smk path/to/input.yaml
   ```

4. Adjust **`snakemake/profiles/local/config.yaml`** `singularity-args` if data live outside **`PROJECT_DIR`**, **`PIPELINE_DATA_DIR`**, or **`QTL_DATA_DIR`**. Profiles bind  **`PROJECT_DIR`**. Add more `-B host:host` pairs if needed.

On first run, `run_pipeline.sh` builds or uses **`$PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<version>.sif`** (or caches under `.snakemake/singularity/` if that directory is not writable).

---

## HPC (Slurm)

Use this on a shared cluster with Slurm and Apptainer/Singularity.

The repo ships Slurm-oriented profiles. Paths and partitions in the bundled examples are **site-specific** — copy and edit them for your cluster.

| Profile directory | Notes |
| ----------------- | ----- |
| `snakemake/profiles/slurm/` | Generic Slurm + Apptainer |
| `snakemake/profiles/local/` | Local Apptainer only (no scheduler) |

**Typical usage:**

```bash
SNAKEMAKE_PROFILE=snakemake/profiles/slurm/
./run_pipeline.sh snakemake/compare_gwases.smk path/to/input.yaml
```

For non-local profiles, `run_pipeline.sh` runs **`module load ${APPTAINER_MODULE}`** — set **`APPTAINER_MODULE`** in `.env` to match your cluster.

**Customize** each profile’s **`config.yaml`** as needed:

- **`SLURM_ACCOUNT`** in `.env` sets **`sbatch --account`**. If unset, the Slurm profile falls back to **`sacctmgr`** output.
- **`SLURM_PARTITION`**, if unset, is inferred by `run_pipeline.sh` from **`sinfo`** (partition marked with `*`); fallback **`compute`**.
- **`cluster:`** block: `sbatch` options (partition, account, walltime, memory, log paths).
- **`singularity-args`:** binds repo code, **`$HOME`**, **`DATA_DIR`**, **`RESULTS_DIR`**, **`PIPELINE_DATA_DIR`**, **`/tmp`**, and **`QTL_DATA_DIR`** when set (via **`GENEHACKMAN_EXTRA_SINGULARITY_BINDS`** in `run_pipeline.sh`). Add `-B host:host` pairs for shared reference data on your site.

**Snakemake version:** `environment.yml` pins **Snakemake 7.x**. Migrating to **Snakemake 8+** changes cluster syntax; bundled profiles target v7 unless you have updated them.

### PBS / Torque / OpenPBS

There is **no PBS profile** in this repository — we do not have access to a PBS cluster to develop or test one. If you use PBS, please **contribute** a profile under `snakemake/profiles/` (copy `snakemake/profiles/slurm/config.yaml` as a starting point and replace the **`cluster:`** block with your `qsub` invocation). See [CONTRIBUTING.md](CONTRIBUTING.md).

---

## macOS

**macOS is not officially supported.** We do not test GeneHackman on Mac, and Apptainer/Singularity behaviour on macOS is fragile compared with Linux.

Running on a Mac **is possible** if you accept extra setup and troubleshooting. Common approaches:

### Apptainer via Lima (or similar)

Snakemake looks for an executable named **`singularity`** on `PATH`. Shell **aliases** in `~/.zshrc` (e.g. `alias singularity='limactl shell apptainer'`) are **not** visible to `run_pipeline.sh` because it uses **bash** and does not load zsh aliases.

Typical setups:

- **Lima VM** with Apptainer — run from an environment where **`apptainer`** or **`singularity`** is a real binary on `PATH`, or add a **`~/bin/singularity`** wrapper that calls `limactl shell <vm> -- apptainer "$@"`.
- **Homebrew Apptainer** (if available for your OS version).

Example Lima mount config (adjust for your VM):

```yaml
mountType: virtiofs
mounts:
  - location: ~
    mountPoint: ~
  - location: /tmp/lima
    mountPoint: /tmp/lima
```

Use **`SNAKEMAKE_PROFILE=snakemake/profiles/local/`** as on Linux.

### SIF cache and `/Users` mounts

Snakemake caches pulled images under **`.snakemake/singularity/`** (or **`--singularity-prefix`**). On **macOS + Lima**, paths under **`/Users/...`** are often VM mounts where **creating large `.sif` files fails** (e.g. read-only file system).

- Set **`SINGULARITY_DIR`** in `.env` to a **writable path inside the Linux VM** (e.g. `/tmp/genehackman_snakemake_singularity`), not only a bind mount that cannot hold the build.
- Or place a pre-built **`$PIPELINE_DATA_DIR/genomic_data/pipeline/genehackman_<version>.sif`** so Snakemake does not build on a problematic mount.
- Avoid **colons** in SIF filenames on macOS; the repo uses **`genehackman_<version>.sif`**.

### Apple Silicon (arm64)

The published image is **linux/amd64**. Use Lima with **Rosetta / x86 Linux**, or run on amd64 Linux / HPC, if pulls or runs fail with architecture errors.

For a supported experience, use **Linux** or **HPC** instead of macOS.

---

## Environment variables (quick reference)

| Variable | Purpose |
| -------- | ------- |
| **`SNAKEMAKE_PROFILE`** | Profile path (default `snakemake/profiles/local/`). HPC example: `snakemake/profiles/slurm/` |
| **`APPTAINER_MODULE`** | Module to load on HPC when not using the local profile |
| **`SLURM_ACCOUNT`** | Optional; sets `sbatch --account` for Slurm profiles |
| **`SLURM_PARTITION`** | Optional; overrides `sinfo` partition detection in `run_pipeline.sh` |

The pipeline input YAML path is **not** set in `.env`. Pass it as the second argument to `run_pipeline.sh`, or use **`--config genehackman_input=...`** with bare `snakemake` (see [PIPELINES.md](PIPELINES.md)).

---

## Further reading

- Pipeline inputs: [PIPELINES.md](PIPELINES.md)
- Main readme: [README.md](README.md)
- Snakemake profiles: [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
- Snakemake cluster execution: [documentation](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
