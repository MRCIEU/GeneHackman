set -e
cp tests/testthat/data/snakemake_inputs/disease_progression.json input.json
snakemake --snakefile snakemake/disease_progression.smk --profile snakemake/bc4/ -F

cp tests/testthat/data/snakemake_inputs/compare_gwases.json input.json
snakemake --snakefile snakemake/compare_gwases.smk --profile snakemake/bc4/ -F

rm input.json
