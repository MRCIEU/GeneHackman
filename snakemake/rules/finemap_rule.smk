rule run_finemapping:
    threads: 6
    resources:
        mem = "48G",
        time = "24:00:00"
    input:
        gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).standardised_gwas,
        clumped_file = lambda wildcards: getattr(pipeline, wildcards.prefix).clumped_file
    params:
        ancestry = lambda wildcards: getattr(pipeline, wildcards.prefix).ancestry,
        n = lambda wildcards: getattr(pipeline, wildcards.prefix).N,
        window_kb = finemap_window_kb,
        max_causal = finemap_max_causal,
        coverage = finemap_coverage,
        min_abs_corr = finemap_min_abs_corr
    output:
        finemap_complete = FINEMAP_COMPLETE_TXT_PATTERN
    shell:
        """
        OUT_FINEMAP="$(dirname {output.finemap_complete})"
        Rscript run_finemap.R \
            --gwas_filename {input.gwas} \
            --clumped_filename {input.clumped_file} \
            --ancestry {params.ancestry} \
            --N {params.n} \
            --output_finemap_dir "$OUT_FINEMAP" \
            --window_kb {params.window_kb} \
            --max_causal {params.max_causal} \
            --coverage {params.coverage} \
            --min_abs_corr {params.min_abs_corr}
        """
