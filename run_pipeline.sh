#!/bin/bash
set -e
SMK_FILE=$1
INPUT_FILE=$2
#if [$(hostname)...] set profile accordingly...

conda activate gwaspipeline || echo "Error: Please ensure conda is initialised and 'gwaspipeline' env is available" && exit 1
module load apps/singularity/3.8.3
tmux new-session -d -s snakemake setenv INPUT_FILE "${INPUT_FILE}" "snakemake --snakefile ${SMK_FILE} --profile snakemake/bc4/"