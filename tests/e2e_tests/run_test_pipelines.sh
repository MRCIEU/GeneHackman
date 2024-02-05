#!/bin/bash
set -e

./run_pipeline.sh snakemake/disease_progression.smk tests/testthat/data/snakemake_inputs/disease_progression.json -F
./run_pipeline.sh snakemake/compare_gwases.smk tests/testthat/data/snakemake_inputs/compare_gwases.json -F

