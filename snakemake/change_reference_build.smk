include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "standardise_gwas"
pipeline = parse_pipeline_input()

onstart:
    print("##### GWAS Standardisation Pipeline #####")

std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases])

rule change_reference_build:
    params:
        input_gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).file,
        input_build = lambda wildcards: getattr(pipeline, wildcards.prefix).build,
        output_build = lambda wildcards: getattr(pipeline, wildcards.prefix).build,
    threads: 8
    resources:
        mem = lambda wildcards: f"{getattr(pipeline, wildcards.prefix).standardised_memory}G",
        time = '04:00:00'
    output: std_file_pattern
    shell:
        """
        INPUT_GWAS={params.input_gwas}
        if [[ {params.input_gwas} =~ .vcf ]]; then
            INPUT_GWAS=$DATA_DIR/gwas/$(basename "{params.input_gwas}" | sed  s/.vcf.*/\.tsv/g)
            ./vcf_to_tsv.sh {params.input_gwas} {params.vcf_columns} $INPUT_GWAS
        fi

        Rscript change_snp_identifiers.R \
            --input_gwas $INPUT_GWAS \
            --output_gwas {output} \
            --input_build {params.input_build} \
            --output_build {params.output_build} \
            --input_columns {params.input_columns} \
            --populate_rsid {params.populate_rsid}
        """



onsuccess:
    onsuccess(pipeline_name, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
