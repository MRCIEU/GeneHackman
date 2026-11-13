rule run_multi_ancestry_finemapping:
    threads: 8
    resources:
        mem = "64G",
        time = "48:00:00"
    input:
        gwas_files = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_file for g in pipeline.gwases]
    params:
        gwas_files = " ".join([g.standardised_gwas for g in pipeline.gwases]),
        clumped_files = " ".join([g.clumped_file for g in pipeline.gwases]),
        ancestries = " ".join([g.ancestry for g in pipeline.gwases]),
        sample_sizes = " ".join([str(g.analysis_N) for g in pipeline.gwases]),
        ld_references = " ".join(["'" + ld_source(g, "finemap") + "'" for g in pipeline.gwases]),
        output_dir = RESULTS_DIR + "finemap/multi_ancestry",
        window_kb = finemap_window_kb,
        max_causal = finemap_max_causal,
        coverage = finemap_coverage,
        min_abs_corr = finemap_min_abs_corr,
        max_region_kb = finemap_max_region_kb,
        region_strategy = pipeline.finemap.region_strategy,
        max_n_ratio = pipeline.finemap.max_n_ratio,
        remove_ld_outliers = "--remove_ld_outliers" if pipeline.finemap.remove_ld_outliers else "",
        build = pipeline.output.build,
        panel_dir = thousand_genomes_dir(pipeline.output.build)
    output:
        finemap_complete = MULTI_FINEMAP_COMPLETE
    shell:
        """
        python run_multisusie.py \
            --gwas_files {params.gwas_files} \
            --clumped_files {params.clumped_files} \
            --ancestries {params.ancestries} \
            --sample_sizes {params.sample_sizes} \
            --output_dir {params.output_dir} \
            --window_kb {params.window_kb} \
            --max_causal {params.max_causal} \
            --coverage {params.coverage} \
            --min_abs_corr {params.min_abs_corr} \
            --max_region_kb {params.max_region_kb} \
            --region_strategy {params.region_strategy} \
            --max_n_ratio {params.max_n_ratio} \
            --ld_references {params.ld_references} {params.remove_ld_outliers} \
            --thousand_genomes_dir {params.panel_dir} \
            --ld_pack_dir {LD_PACK_DIR} \
            --build {params.build} \
            --completion_file {output.finemap_complete}
        """
