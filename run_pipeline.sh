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

SIF_DIR="${SINGULARITY_PREFIX:-.snakemake/singularity}"
SIF_PATH="${SIF_DIR}/genehackman:${DOCKER_VERSION}.sif"

if [[ ! -f "${SIF_PATH}" ]]; then
  mkdir -p "${SIF_DIR}"
  echo "SIF file not found: ${SIF_PATH}"
  echo "Building container with singularity..."
  singularity build "${SIF_PATH}" "docker://mrcieu/genehackman:${GENEHACKMAN_VERSION}"
fi


echo "Running pipeline with profile: ${PROFILE}"

snakemake --snakefile "${SMK_FILE}" --profile "${PROFILE}" ${ADDITIONAL_ARGS}
