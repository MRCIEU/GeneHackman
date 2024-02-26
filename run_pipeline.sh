#!/bin/bash
set -e

if [[ $# -lt 1 ]] ; then
  echo """
  Error: You have to provide at least 1 argument:
    PIPELINE_FILE (ex. snakemake/standardise_gwas.smk)
    INPUT_FILE (defaults to input.json)
  """
  exit 1
fi

#if [$(hostname)...] set PROFILE accordingly, for when multiple HPCs are 
PROFILE=snakemake/bc4/

SMK_FILE=$1
export INPUT_FILE=$2
ADDITIONAL_ARGS=$3

module load apps/singularity/3.8.3
snakemake --snakefile ${SMK_FILE} --profile ${PROFILE} ${ADDITIONAL_ARGS}
