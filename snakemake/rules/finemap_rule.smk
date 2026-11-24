rule run_finemapping:
    threads: 8
    resources:
        mem = "48G",
        time = "24:00:00"
    input:
        gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).standardised_gwas,
        clumped_file = lambda wildcards: getattr(pipeline, wildcards.prefix).clumped_file
    params:
        ancestry = lambda wildcards: getattr(pipeline, wildcards.prefix).ancestry,
        n = lambda wildcards: getattr(pipeline, wildcards.prefix).analysis_N,
        finemap_dir = lambda wildcards: getattr(pipeline, wildcards.prefix).finemap_dir,
        window_kb = finemap_window_kb,
        max_causal = finemap_max_causal,
        coverage = finemap_coverage,
        min_abs_corr = finemap_min_abs_corr,
        max_region_kb = finemap_max_region_kb,
        ld_reference = lambda wildcards: ld_source(getattr(pipeline, wildcards.prefix), "finemap"),
        build = pipeline.output.build,
        region_strategy = pipeline.finemap.region_strategy,
        max_n_ratio = pipeline.finemap.max_n_ratio,
        remove_ld_outliers = "TRUE" if pipeline.finemap.remove_ld_outliers else "FALSE"
    output:
        finemap_complete = FINEMAP_COMPLETE_TXT_PATTERN
    shell:
        """
        Rscript run_finemap.R \
            --gwas_filename {input.gwas} \
            --clumped_filename {input.clumped_file} \
            --ancestry {params.ancestry} \
            --N {params.n} \
            --output_finemap_dir {params.finemap_dir} \
            --completion_file {output.finemap_complete} \
            --window_kb {params.window_kb} \
            --max_causal {params.max_causal} \
            --coverage {params.coverage} \
            --min_abs_corr {params.min_abs_corr} \
            --max_region_kb {params.max_region_kb} \
            --ld_reference '{params.ld_reference}' \
            --build {params.build} \
            --region_strategy {params.region_strategy} \
            --max_n_ratio {params.max_n_ratio} \
            --remove_ld_outliers {params.remove_ld_outliers}
        """
