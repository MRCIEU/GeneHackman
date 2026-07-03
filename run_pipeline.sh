#!/bin/bash
set -e

export GRPC_VERBOSITY=NONE
export GRPC_TRACE=
export APPTAINER_SILENT=true

if [[ $# -lt 1 || "$1" =~ "help" ]] ; then
  echo """
  Error: You have to provide at least 1 argument:
    PIPELINE_FILE (ex. snakemake/standardise_gwas.smk)
    INPUT_FILE YAML (optional, defaults to input.yaml; omit when passing only Snakemake flags)

  Examples:
    ./run_pipeline.sh snakemake/coloc.smk
    ./run_pipeline.sh snakemake/coloc.smk my_input.yaml
    ./run_pipeline.sh snakemake/coloc.smk --unlock
    ./run_pipeline.sh snakemake/coloc.smk my_input.yaml --dry-run

  Default profile is local Apptainer (snakemake/profiles/apptainer/). For local Docker, set
  SNAKEMAKE_PROFILE=snakemake/profiles/docker/. For Slurm, set
  SNAKEMAKE_PROFILE=snakemake/profiles/slurm/ or create a new profile under snakemake/profiles/
  """
  exit 1
fi

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENV_FILE="$(pwd)/.env"

if [ -f "${ENV_FILE}" ]
then
  export $(grep -vE '^[[:space:]]*#' "${ENV_FILE}" | xargs)
else
  echo "Error: .env file missing"
  exit 1
fi

if [[ -z "${DOCKER_VERSION:-}" ]]; then
  DOCKER_VERSION="$(grep -E '^Version:' "${REPO_ROOT}/DESCRIPTION" 2>/dev/null | awk '{print $2}' | head -n1)"
  export DOCKER_VERSION
fi

# GWAS working dirs: always PROJECT_DIR/data and PROJECT_DIR/results (exported as DATA_DIR / RESULTS_DIR for shell rules and Singularity binds).
if [[ -z "${PROJECT_DIR:-}" || -z "${PIPELINE_DATA_DIR:-}" ]]; then
  echo "Error: PROJECT_DIR and PIPELINE_DATA_DIR must be set in .env."
  exit 1
fi
_pd="$(echo "${PROJECT_DIR}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s|/*$||')"
export DATA_DIR="${_pd}/data"
export RESULTS_DIR="${_pd}/results"
unset _pd
export PIPELINE_LOG_DIR="${DATA_DIR%/}/snakemake_logs"
mkdir -p $PIPELINE_LOG_DIR $DATA_DIR $RESULTS_DIR

if [[ -z "${SLURM_PARTITION:-}" ]]; then
  if command -v sinfo >/dev/null 2>&1; then
    SLURM_PARTITION="$(sinfo -h -o '%P' 2>/dev/null | grep '\*' | head -n1 | tr -d '*[:space:]')"
  fi
  SLURM_PARTITION="${SLURM_PARTITION:-compute}"
fi
export SLURM_PARTITION

if [[ -z "${SLURM_ACCOUNT:-}" ]] && command -v sacctmgr >/dev/null 2>&1; then
  export USER=$(whoami)
  export ACCOUNT_ID=$(sacctmgr show user withassoc format=account where user="$USER" | grep '[0-9]' | head -n1)
  export SLURM_ACCOUNT="${ACCOUNT_ID:-}"
fi

export GENEHACKMAN_EXTRA_SINGULARITY_BINDS=""
_trim_qtl="$(echo "${QTL_DATA_DIR:-}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
if [[ -n "${_trim_qtl}" ]]; then
  export GENEHACKMAN_EXTRA_SINGULARITY_BINDS="-B ${_trim_qtl%/}:${_trim_qtl%/}"
fi
unset _trim_qtl

# Default: local Apptainer. Local Docker: SNAKEMAKE_PROFILE=snakemake/profiles/docker/.
# HPC e.g.: SNAKEMAKE_PROFILE=snakemake/profiles/slurm/.
PROFILE="${SNAKEMAKE_PROFILE:-snakemake/profiles/apptainer/}"
if [[ -z "${SNAKEMAKE_PROFILE:-}" && ( "${GENEHACKMAN_LOCAL:-}" == "1" || "${GENEHACKMAN_LOCAL:-}" == "true" ) ]]; then
  PROFILE="snakemake/profiles/apptainer/"
fi
PROFILE_TRIMMED="${PROFILE%/}"
if [[ "${PROFILE_TRIMMED}" == "snakemake/profiles/docker" ]]; then
  CONTAINER_RUNTIME="docker"
  export GENEHACKMAN_LOCAL=1
elif [[ "${PROFILE}" == snakemake/profiles/apptainer/* ]] || [[ "${PROFILE}" == snakemake/apptainer/* ]] || [[ "${PROFILE}" == snakemake/profiles/local/* ]] || [[ "${PROFILE}" == snakemake/local/* ]]; then
  CONTAINER_RUNTIME="apptainer"
  export GENEHACKMAN_LOCAL=1
else
  CONTAINER_RUNTIME="apptainer"
  unset GENEHACKMAN_LOCAL
fi
export GENEHACKMAN_CONTAINER_RUNTIME="${CONTAINER_RUNTIME}"

SMK_FILE=$1
shift
# Second arg is the input YAML only if it is present and not a Snakemake/cli flag (e.g. --unlock).
if [[ $# -ge 1 && "${1}" != -* ]]; then
  PIPELINE_INPUT="$1"
  shift
else
  PIPELINE_INPUT="input.yaml"
fi

if [[ "${CONTAINER_RUNTIME}" == "apptainer" && "${PROFILE}" != snakemake/profiles/apptainer/* ]] && [[ "${PROFILE}" != snakemake/apptainer/* ]] && [[ "${PROFILE}" != snakemake/profiles/local/* ]] && [[ "${PROFILE}" != snakemake/local/* ]]; then
  module load ${APPTAINER_MODULE}
fi

# Image tag: DOCKER_VERSION defaults to DESCRIPTION Version:; override in .env if needed.
GENEHACKMAN_VERSION="${GENEHACKMAN_VERSION:-${DOCKER_VERSION:-}}"
SIF_VERSION="${DOCKER_VERSION:-${GENEHACKMAN_VERSION:-}}"
if [[ -z "${SIF_VERSION}" ]]; then
  echo "Error: could not determine Docker image version. Set DOCKER_VERSION in .env or ensure DESCRIPTION contains Version:."
  exit 1
fi

CONTAINER_IMAGE="mrcieu/genehackman:${GENEHACKMAN_VERSION}"

if [[ "${CONTAINER_RUNTIME}" == "apptainer" ]]; then
  SIF_NAME="genehackman_${SIF_VERSION}.sif"
  PIPELINE_GENOMIC_DIR="${PIPELINE_DATA_DIR%/}/genomic_data/pipeline"
  PIPELINE_GENOMIC_DIR="${PIPELINE_GENOMIC_DIR%/}"
  SIF_PATH="${PIPELINE_GENOMIC_DIR}/${SIF_NAME}"

  echo "SIF_PATH: ${SIF_PATH}"
  if [[ ! -f "${SIF_PATH}" ]]; then
    if ! mkdir -p "${PIPELINE_GENOMIC_DIR}" 2>/dev/null || [[ ! -w "${PIPELINE_GENOMIC_DIR}" ]]; then
      SIF_PATH=".snakemake/singularity/${SIF_NAME}"
      echo "Note: caching SIF under ${SIF_PATH} (${PIPELINE_GENOMIC_DIR} missing or not writable)"
    fi
  fi
  mkdir -p "$(dirname "${SIF_PATH}")"

  SINGULARITY_DOCKER_REFERENCE="docker://${CONTAINER_IMAGE}"

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
else
  if ! command -v docker >/dev/null 2>&1; then
    echo "Error: Docker profile selected but docker is not on PATH."
    exit 1
  fi
  echo "Using Docker image: ${CONTAINER_IMAGE}"
fi


echo "Running pipeline with profile: ${PROFILE}"

export GENEHACKMAN_INPUT="${PIPELINE_INPUT}"
if [[ "${CONTAINER_RUNTIME}" == "docker" ]]; then
  DOCKER_PLATFORM="${GENEHACKMAN_DOCKER_PLATFORM:-linux/amd64}"
  DOCKER_MOUNTS=(
    -v "${REPO_ROOT}:/workspace"
    -v "${REPO_ROOT}/R:/home/R"
    -v "${REPO_ROOT}/scripts:/home/scripts"
    -v "${REPO_ROOT}/inst:/home/inst"
    -v "${HOME}:${HOME}"
    -v "${DATA_DIR}:${DATA_DIR}"
    -v "${RESULTS_DIR}:${RESULTS_DIR}"
    -v "${PIPELINE_DATA_DIR}:${PIPELINE_DATA_DIR}"
    -v "/tmp:/tmp"
  )
  if [[ -n "${WORK:-}" ]]; then
    DOCKER_MOUNTS+=(-v "${WORK}:${WORK}")
  fi
  if [[ -n "${QTL_DATA_DIR:-}" ]]; then
    DOCKER_MOUNTS+=(-v "${QTL_DATA_DIR%/}:${QTL_DATA_DIR%/}")
  fi
  if [[ "${PIPELINE_INPUT}" = /* ]]; then
    DOCKER_MOUNTS+=(-v "$(dirname "${PIPELINE_INPUT}")":"$(dirname "${PIPELINE_INPUT}")")
  fi

  docker run --rm \
    --platform "${DOCKER_PLATFORM}" \
    --user "$(id -u):$(id -g)" \
    --env-file "${ENV_FILE}" \
    -e "HOME=${HOME}" \
    -e "PROJECT_DIR=${PROJECT_DIR}" \
    -e "DATA_DIR=${DATA_DIR}" \
    -e "RESULTS_DIR=${RESULTS_DIR}" \
    -e "PIPELINE_DATA_DIR=${PIPELINE_DATA_DIR}" \
    -e "PIPELINE_LOG_DIR=${PIPELINE_LOG_DIR}" \
    -e "QTL_DATA_DIR=${QTL_DATA_DIR:-}" \
    -e "DOCKER_VERSION=${DOCKER_VERSION}" \
    -e "GENEHACKMAN_VERSION=${GENEHACKMAN_VERSION}" \
    -e "GENEHACKMAN_INPUT=${PIPELINE_INPUT}" \
    -e "GENEHACKMAN_CONTAINER_RUNTIME=${GENEHACKMAN_CONTAINER_RUNTIME}" \
    "${DOCKER_MOUNTS[@]}" \
    -w /workspace \
    "${CONTAINER_IMAGE}" \
    snakemake --snakefile "/workspace/${SMK_FILE}" --profile "/workspace/${PROFILE}" "$@"
else
  snakemake --snakefile "${SMK_FILE}" --profile "${PROFILE}" "$@"
fi
