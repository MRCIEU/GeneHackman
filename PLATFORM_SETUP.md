# Running GeneHackman by platform

Snakemake orchestrates steps; each step usually runs inside an **Apptainer/Singularity** image (`docker://mrcieu/genehackman` or a local `.sif`). **Profiles** under `snakemake/` choose how jobs run (local cores vs cluster scheduler) and which host paths are **bind-mounted** into the container.

General prerequisites:

1. **Conda env:** `conda env create -f environment.yml` then `conda activate genehackman`
2. **`.env`:** copy from `.env_example` and set at least `DATA_DIR`, `RESULTS_DIR`, `PIPELINE_DATA_DIR` (and paths your pipeline needs for genomic / 1000G data). Optionally set **`QTL_DATA_DIR`** when large QTL mirrors live on another volume or object-store mount; if omitted, `run_pipeline.sh` defaults it to `PIPELINE_DATA_DIR/qtl_datasets`.
3. **Input JSON:** see `snakemake/input_templates/` and `snakemake/PIPELINES.md`.

**`./run_pipeline.sh`** loads `.env`, picks a **profile** (`GENEHACKMAN_PROFILE`, default `snakemake/local/`), then runs Snakemake. Pass the **`.smk` workflow first**, then optional input JSON, for example:

```bash
./run_pipeline.sh snakemake/standardise_gwas.smk path/to/input.json
```

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
- Alternatively, place a pre-built **`$PIPELINE_DATA_DIR/genehackman_<version>.sif`** (underscore; see `DESCRIPTION` / `DOCKER_VERSION`) so Snakemake does not need to build on a problematic mount.
- Avoid **colons** in host-side SIF filenames (`genehackman:1.0.0.sif`); macOS can reject them. The repo’s Snakemake helper uses **`genehackman_<version>.sif`**.

### Apple Silicon (arm64)

The published image may be **linux/amd64**. Use Lima with **Rosetta / x86 Linux** (or run on an amd64 Linux machine / HPC) if pulls or runs fail with architecture warnings.

---

## Linux (workstation or single server)

1. Install **Apptainer** or **SingularityCE** (distribution packages or upstream instructions).
2. Ensure **`singularity`** or **`apptainer`** is on `PATH` (some sites symlink `singularity` → `apptainer`).
3. Use **`snakemake/local/`** for local execution (same as default in `run_pipeline.sh`):

   ```bash
   export GENEHACKMAN_PROFILE=snakemake/local/
   ./run_pipeline.sh snakemake/compare_gwases.smk
   ```

4. Adjust **`snakemake/local/config.yaml`** `singularity-args` if your data live outside `DATA_DIR` / `RESULTS_DIR` / `PIPELINE_DATA_DIR` / `QTL_DATA_DIR` (profiles bind these four; add more `-B host:host` pairs if needed).
5. Large pulls: set **`GENEHACKMAN_SINGULARITY_PREFIX`** to fast scratch (e.g. `/scratch/$USER/snakemake_singularity`).

---

## Slurm (HPC)

The repo ships Slurm-oriented profiles (paths and partitions are **Bristol / MRC IEU–oriented**; **edit** them for your site):

| Profile directory | Notes |
|-------------------|--------|
| `snakemake/bp1/` | Example: binds `/bp1/mrcieu1/data/`, uses `sbatch`, partition `compute` in `config.yaml` |
| `snakemake/bc4/` | Similar; binds `/mnt/storage/private/mrcieu/data/` in the version checked into this repo |
| `snakemake/slurm_singularity/` | Slurm + Singularity; site-specific binds |

**Typical usage:**

```bash
export GENEHACKMAN_PROFILE=snakemake/bp1/
./run_pipeline.sh snakemake/compare_gwases.smk
```

`run_pipeline.sh` runs **`module load apptainer/1.3.1-ksax`** when the profile is **not** `snakemake/local/*` — **change the module name** in `run_pipeline.sh` to match your cluster.

**You must customize** each profile’s `config.yaml`:

- **`cluster:`** block: `sbatch` options (`partition`, `account`, walltime, memory, `--output` log directory).
- **`singularity-args`:** `-B` mounts for your home, scratch, and shared reference data (replacing Bristol-specific paths). Profiles include **`${QTL_DATA_DIR}`**; use **`./run_pipeline.sh`** (or export **`QTL_DATA_DIR`** yourself, defaulting to **`PIPELINE_DATA_DIR/qtl_datasets`**) so the bind is never empty.

Log files: profiles often write Slurm stdout under **`/user/work/$USER/slurm_logs/`**; create that tree or change `--output` to your scratch.

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

Copy `snakemake/slurm_singularity/config.yaml` (or `bp1`) as a template, replace the **`cluster:`** section with **`qsub`**, and add **`singularity-args`** appropriate for PBS compute nodes (shared filesystem mounts).

---

## Environment variables (quick reference)

| Variable | Purpose |
|----------|---------|
| **`GENEHACKMAN_PROFILE`** | Snakemake profile path (default `snakemake/local/`). Example HPC: `snakemake/bp1/` |
| **`INPUT_FILE`** | Path to pipeline JSON (also set by `run_pipeline.sh` from the second argument) |
| **`GENEHACKMAN_SINGULARITY_PREFIX`** | Snakemake `--singularity-prefix`: where pulled `.sif` images are stored |
| **`GENEHACKMAN_LIMA_INSTANCE`** | Lima VM name when using a `limactl shell <name> -- apptainer` wrapper (default often `apptainer`) |
| **`PIPELINE_LOG_DIR`** | Optional override for log-directory hints on errors (see `snakemake/util/constants.smk`) |

Use **`.env`** for values you want loaded every time `./run_pipeline.sh` runs (`export $(cat .env | xargs)`).

---

## Further reading

- Pipeline input schema: [`snakemake/PIPELINES.md`](snakemake/PIPELINES.md)
- Main readme: [`README.md`](README.md)
- Snakemake profiles: [documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles)
- Snakemake cluster execution: [docs](https://snakemake.readthedocs.io/en/stable/executing/cluster.html)
