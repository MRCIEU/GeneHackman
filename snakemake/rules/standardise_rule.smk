rule standardise_gwases:
    params:
        input_gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).file,
        N = lambda wildcards: getattr(pipeline, wildcards.prefix).N,
        remove_extra_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).remove_extra_columns,
        vcf_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).vcf_columns,
        input_build = lambda wildcards: getattr(pipeline, wildcards.prefix).build,
        input_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).input_columns,
        output_columns = lambda wildcards: getattr(pipeline, wildcards.prefix).output_columns,
        populate_rsid = lambda wildcards: getattr(pipeline, wildcards.prefix).populate_rsid.value
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

        Rscript standardise_gwas.R \
            --input_gwas $INPUT_GWAS \
            --output_gwas {output} \
            --N {params.N} \
            --remove_extra_columns {params.remove_extra_columns} \
            --input_build {params.input_build} \
            --output_build {pipeline.output.build} \
            --input_columns {params.input_columns} \
            --output_columns {params.output_columns} \
            --populate_rsid {params.populate_rsid}
        """

