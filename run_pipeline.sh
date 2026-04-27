#!/bin/bash
set -e

export GRPC_VERBOSITY=NONE
export GRPC_TRACE=

if [[ $# -lt 1 ]] ; then
  echo """
  Error: You have to provide at least 1 argument:
    PIPELINE_FILE (ex. snakemake/standardise_gwas.smk)
    INPUT_FILE (optional, defaults to input.json)

  Default profile is local Docker (snakemake/local/). HPC (Slurm/Apptainer): set
  GENEHACKMAN_PROFILE=snakemake/bp1/ or snakemake/bc4/ (or snakemake/slurm_singularity/).
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

# Default: local Docker. HPC: GENEHACKMAN_PROFILE=snakemake/bp1/ (or bc4/, slurm_singularity/).
PROFILE="${GENEHACKMAN_PROFILE:-snakemake/local/}"
if [[ -z "${GENEHACKMAN_PROFILE:-}" && ( "${GENEHACKMAN_LOCAL:-}" == "1" || "${GENEHACKMAN_LOCAL:-}" == "true" ) ]]; then
  PROFILE="snakemake/local/"
fi
if [[ "${PROFILE}" == snakemake/local/* ]]; then
  export GENEHACKMAN_LOCAL=1
else
  unset GENEHACKMAN_LOCAL
fi

SMK_FILE=$1
export INPUT_FILE="${2:-input.json}"
ADDITIONAL_ARGS="${@:3}"

if [[ "${PROFILE}" != snakemake/local/* ]]; then
  module load apptainer/1.3.1-ksax
fi

SIF_DIR="${SINGULARITY_DIR:-.snakemake/singularity}"
SIF_DIR="${SIF_DIR%/}"
# Image tag: default GENEHACKMAN_VERSION from DOCKER_VERSION so pull and SIF basename stay aligned.
GENEHACKMAN_VERSION="${GENEHACKMAN_VERSION:-${DOCKER_VERSION:-}}"
SIF_VERSION="${DOCKER_VERSION:-${GENEHACKMAN_VERSION:-}}"
if [[ -z "${SIF_VERSION}" ]]; then
  echo "Error: Set DOCKER_VERSION in .env (e.g. 1.0.0) so the SIF name matches the Docker image tag."
  exit 1
fi
SIF_NAME="genehackman_${SIF_VERSION}.sif"
SIF_PATH="${SIF_DIR}/${SIF_NAME}"

# Source for singularity build: default pulls docker://mrcieu/genehackman:<tag> from Hub.
# Set GENEHACKMAN_SINGULARITY_DOCKER_REFERENCE=docker-daemon://mrcieu/genehackman:develop (or other tag)
# to build the SIF from your local Docker image (same image ID, no Hub pull / no rebuild of layers).
SINGULARITY_DOCKER_REFERENCE="${GENEHACKMAN_SINGULARITY_DOCKER_REFERENCE:-docker://mrcieu/genehackman:${GENEHACKMAN_VERSION}}"

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
  mkdir -p "${SIF_DIR}"
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

snakemake --snakefile "${SMK_FILE}" --profile "${PROFILE}" ${ADDITIONAL_ARGS}
