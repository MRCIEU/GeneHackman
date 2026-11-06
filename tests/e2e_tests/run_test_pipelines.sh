#!/bin/bash
set -e

./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas.yaml -F
./run_pipeline.sh snakemake/disease_progression.smk tests/testthat/data/snakemake_inputs/disease_progression.yaml -F
./run_pipeline.sh snakemake/compare_gwases.smk tests/testthat/data/snakemake_inputs/compare_gwases.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_multi_ancestry.yaml -F
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc.yaml -F

# The same pipelines with the Pan-UKB LD pack (ld_reference: pan_ukb), downloaded with each ancestry
# (PIPELINE_DATA_DIR/genomic_data/pan_ukb_ld)
./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas_pan_ukb.yaml -F
./run_pipeline.sh snakemake/disease_progression.smk tests/testthat/data/snakemake_inputs/disease_progression_pan_ukb.yaml -F
./run_pipeline.sh snakemake/compare_gwases.smk tests/testthat/data/snakemake_inputs/compare_gwases_pan_ukb.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_pan_ukb.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_multi_ancestry_pan_ukb.yaml -F
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc_pan_ukb.yaml -F

# GRCh37 output (opt-in; needs download_reference_data.sh --grch37), with both LD references
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_grch37.yaml -F
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc_grch37.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_pan_ukb_grch37.yaml -F
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc_pan_ukb_grch37.yaml -F

if [[ -n "${QTL_DATA_DIR}" ]]; then
  for qtl_dataset in eqtlgen brainmeta ukb_pqtl; do
    ./run_pipeline.sh snakemake/qtl_mr.smk "tests/testthat/data/snakemake_inputs/qtl_mr_${qtl_dataset}.yaml" -F
  done
fi

branch_name=$(git rev-parse --abbrev-ref HEAD)
status_message="SUCCESS: All tests passed on branch: $branch_name"
echo "$status_message" > tests/testing_complete.txt