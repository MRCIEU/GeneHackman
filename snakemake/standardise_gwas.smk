include: "util/common.smk"
singularity: docker_container

pipeline_name = Path(__file__).name
pipeline = parse_pipeline_input()

onstart:
    print("##### GWAS Standardisation Pipeline #####")

std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases])

include: "rules/standardise_rule.smk"

onsuccess:
    onsuccess(pipeline_name, test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, test=pipeline.is_test)
