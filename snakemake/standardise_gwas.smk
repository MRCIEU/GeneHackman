include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "standardise_gwas"
pipeline = parse_pipeline_input(allow_flip_false=True)

onstart:
    print("##### GWAS Standardisation Pipeline #####")

std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases])

include: "rules/standardise_rule.smk"

onsuccess:
    onsuccess(
        pipeline_name,
        is_test=pipeline.is_test,
        files_created=expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases])
    )

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
