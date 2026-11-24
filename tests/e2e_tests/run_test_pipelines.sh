#!/bin/bash
set -e

./run_pipeline.sh snakemake/standardise_gwas.smk tests/testthat/data/snakemake_inputs/standardise_gwas.yaml -F
./run_pipeline.sh snakemake/disease_progression.smk tests/testthat/data/snakemake_inputs/disease_progression.yaml -F
./run_pipeline.sh snakemake/compare_gwases.smk tests/testthat/data/snakemake_inputs/compare_gwases.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_multi_ancestry.yaml -F
# coloc also covers choosing LD per step and per GWAS: one GWAS uses a PLINK prefix as its own genotypes for the
# strand check, clumping, fine-mapping and regional plot, and the other takes its regional-plot LD from Pan-UKB
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc.yaml -F

# The Pan-UKB LD pack (ld_reference: pan_ukb, downloaded with each ancestry): coloc covers standardisation with
# the strand check against the pack, clumping, SuSiE fine-mapping and colocalization (R reader); multi-ancestry
# fine-mapping covers MultiSuSiE, which reads the pack with the Python reader.
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc_pan_ukb.yaml -F
./run_pipeline.sh snakemake/finemap.smk tests/testthat/data/snakemake_inputs/finemap_multi_ancestry_pan_ukb.yaml -F

# One run with GRCh37 output (opt-in; needs download_reference_data.sh --grch37): conversion to GRCh37, the
# GRCh37 1000 Genomes panel, fine-mapping, colocalization and GRCh37 gene annotation in the regional plots.
./run_pipeline.sh snakemake/coloc.smk tests/testthat/data/snakemake_inputs/coloc_grch37.yaml -F

# QTL MR (needs QTL_DATA_DIR): one dataset, as the three use the same code
if [[ -n "${QTL_DATA_DIR}" ]]; then
  ./run_pipeline.sh snakemake/qtl_mr.smk tests/testthat/data/snakemake_inputs/qtl_mr_eqtlgen.yaml -F
fi

branch_name=$(git rev-parse --abbrev-ref HEAD)
status_message="SUCCESS: All tests passed on branch: $branch_name"
echo "$status_message" > tests/testing_complete.txt