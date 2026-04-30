#!/bin/bash
set -e

export GRPC_VERBOSITY=NONE
export GRPC_TRACE=

if [[ $# -lt 1 ]] ; then
  echo """
  Error: You have to provide at least 1 argument:
    PIPELINE_FILE (ex. snakemake/standardise_gwas.smk)
    INPUT_FILE YAML (optional, defaults to input.yaml; omit when passing only Snakemake flags)

  Examples:
    ./run_pipeline.sh snakemake/coloc.smk
    ./run_pipeline.sh snakemake/coloc.smk my_input.yaml
    ./run_pipeline.sh snakemake/coloc.smk --unlock
    ./run_pipeline.sh snakemake/coloc.smk my_input.yaml --dry-run

  Default profile is local Apptainer (snakemake/profiles/local/). For Slurm, set
  SNAKEMAKE_PROFILE=snakemake/profiles/slurm/ or create a new profile under snakemake/profiles/
  """
  exit 1
fi

if [ -f .env ]
then
  export $(cat .env | xargs)
else
  echo "Error: .env file missing"
  exit 1
fi

if [[ -z "${DATA_DIR:-}" || -z "${RESULTS_DIR:-}" || -z "${PIPELINE_DATA_DIR:-}" ]]; then
  echo "Error: DATA_DIR, RESULTS_DIR, and PIPELINE_DATA_DIR must be set in .env"
  exit 1
fi
export PIPELINE_LOG_DIR="${DATA_DIR%/}/snakemake_logs"

if [[ -z "${SLURM_PARTITION:-}" ]]; then
  if command -v sinfo >/dev/null 2>&1; then
    SLURM_PARTITION="$(sinfo -h -o '%P' 2>/dev/null | grep '\*' | head -n1 | tr -d '*[:space:]')"
  fi
  SLURM_PARTITION="${SLURM_PARTITION:-compute}"
fi
export SLURM_PARTITION

export GENEHACKMAN_EXTRA_SINGULARITY_BINDS=""
_trim_qtl="$(echo "${QTL_DATA_DIR:-}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
if [[ -n "${_trim_qtl}" ]]; then
  export GENEHACKMAN_EXTRA_SINGULARITY_BINDS="-B ${_trim_qtl%/}:${_trim_qtl%/}"
fi
unset _trim_qtl

# Default: local Apptainer. HPC e.g.: SNAKEMAKE_PROFILE=snakemake/profiles/slurm/ or snakemake/slurm_singularity/
PROFILE="${SNAKEMAKE_PROFILE:-snakemake/profiles/local/}"
if [[ -z "${SNAKEMAKE_PROFILE:-}" && ( "${GENEHACKMAN_LOCAL:-}" == "1" || "${GENEHACKMAN_LOCAL:-}" == "true" ) ]]; then
  PROFILE="snakemake/profiles/local/"
fi
if [[ "${PROFILE}" == snakemake/profiles/local/* ]] || [[ "${PROFILE}" == snakemake/local/* ]]; then
  export GENEHACKMAN_LOCAL=1
else
  unset GENEHACKMAN_LOCAL
fi

SMK_FILE=$1
shift
# Second arg is the input YAML only if it is present and not a Snakemake/cli flag (e.g. --unlock).
if [[ $# -ge 1 && "${1}" != -* ]]; then
  PIPELINE_INPUT="$1"
  shift
else
  PIPELINE_INPUT="input.yaml"
fi

if [[ "${PROFILE}" != snakemake/profiles/local/* ]] && [[ "${PROFILE}" != snakemake/local/* ]]; then
  module load ${APPTAINER_MODULE}
fi

# Image tag: default GENEHACKMAN_VERSION from DOCKER_VERSION so pull and SIF basename stay aligned.
GENEHACKMAN_VERSION="${GENEHACKMAN_VERSION:-${DOCKER_VERSION:-}}"
SIF_VERSION="${DOCKER_VERSION:-${GENEHACKMAN_VERSION:-}}"
if [[ -z "${SIF_VERSION}" ]]; then
  echo "Error: Set DOCKER_VERSION in .env (e.g. 1.0.0) so the SIF name matches the Docker image tag."
  exit 1
fi

SIF_NAME="genehackman_${SIF_VERSION}.sif"
PIPELINE_GENOMIC_DIR="${PIPELINE_DATA_DIR%/}/pipeline"
PIPELINE_GENOMIC_DIR="${PIPELINE_GENOMIC_DIR%/}"
SIF_PATH="${PIPELINE_GENOMIC_DIR}/${SIF_NAME}"

echo "SIF_PATH: ${SIF_PATH}"
if [[ -f "${SIF_PATH}" ]] || [[ -w "${PIPELINE_GENOMIC_DIR}" ]] ; then
  mkdir -p "${PIPELINE_GENOMIC_DIR}"
else
  SIF_PATH=".snakemake/singularity/${SIF_NAME}"
fi  

SINGULARITY_DOCKER_REFERENCE="docker://mrcieu/genehackman:${GENEHACKMAN_VERSION}"

# mrcieu/genehackman is linux/amd64 only. On ARM (Apple Silicon), Apptainer defaults to arm64 and
# fails with: no child with platform linux/arm64 in index. Force amd64 for the OCI/docker pull.
SINGULARITY_BUILD_ARCH_ARGS=()
if [[ "${GENEHACKMAN_SINGULARITY_NO_ARCH:-}" != "1" && "${GENEHACKMAN_SINGULARITY_NO_ARCH:-}" != "true" ]]; then
  case "$(uname -m)" in
    aarch64|arm64)
      SINGULARITY_BUILD_ARCH_ARGS=(--arch "${GENEHACKMAN_SINGULARITY_BUILD_ARCH:-amd64}")
      ;;
  esac
fi

if [[ ! -f "${SIF_PATH}" ]]; then
  echo "SIF file not found: ${SIF_PATH}"
  mkdir -p "${PIPELINE_GENOMIC_DIR}"
  echo "Building container with singularity from: ${SINGULARITY_DOCKER_REFERENCE}"
  if [[ "${#SINGULARITY_BUILD_ARCH_ARGS[@]}" -gt 0 ]]; then
    echo "(host is ARM: using singularity build ${SINGULARITY_BUILD_ARCH_ARGS[*]} for amd64 image)"
  fi
  # Build next to the final path (same bind-mount as SIF_DIR). Do not use host /tmp: if singularity
  # runs inside Lima/VM, VM /tmp is not the Mac's /tmp, so a follow-up cp from /tmp would fail.
  SIF_BUILD_TMP="${SIF_PATH}.tmp.$$"
  if singularity build "${SINGULARITY_BUILD_ARCH_ARGS[@]}" "${SIF_BUILD_TMP}" "${SINGULARITY_DOCKER_REFERENCE}"; then
    mv -f "${SIF_BUILD_TMP}" "${SIF_PATH}"
  else
    rm -f "${SIF_BUILD_TMP}"
    exit 1
  fi
else
  echo "Using pre-built SIF file: ${SIF_PATH}"
fi


echo "Running pipeline with profile: ${PROFILE}"

snakemake --snakefile "${SMK_FILE}" --profile "${PROFILE}" --config "genehackman_input=${PIPELINE_INPUT}" "$@"
