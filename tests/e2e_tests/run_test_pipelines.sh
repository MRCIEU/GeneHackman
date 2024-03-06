#!/bin/bash
set -e

./run_pipeline.sh snakemake/qtl_mr.smk tests/testthat/data/snakemake_inputs/qtl_mr_metabrain.json -F
./run_pipeline.sh snakemake/qtl_mr.smk tests/testthat/data/snakemake_inputs/qtl_mr_eqtlgen.json
exit 0
./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas.json -F &
./run_pipeline.sh snakemake/disease_progression.smk tests/testthat/data/snakemake_inputs/disease_progression.json -F & 
./run_pipeline.sh snakemake/compare_gwases.smk tests/testthat/data/snakemake_inputs/compare_gwases.json -F
