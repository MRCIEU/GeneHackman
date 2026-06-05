#!/bin/bash
set -e

./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas.yaml -F
./run_pipeline.sh snakemake/disease_progression.smk tests/testthat/data/snakemake_inputs/disease_progression.yaml -F
./run_pipeline.sh snakemake/compare_gwases.smk tests/testthat/data/snakemake_inputs/compare_gwases.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_multi_ancestry.yaml -F
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc.yaml -F

if [[ -n "${QTL_DATA_DIR}" ]]; then
  ./run_pipeline.sh snakemake/qtl_mr.smk tests/testthat/data/snakemake_inputs/qtl_mr_eqtlgen.yaml -F
fi

branch_name=$(git rev-parse --abbrev-ref HEAD)
status_message="SUCCESS: All tests passed on branch: $branch_name"
echo "$status_message" > tests/testing_complete.txt