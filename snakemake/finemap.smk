include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "finemap"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)

onstart:
    print("##### SuSiE Fine-Mapping Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)

finemap_opts = getattr(pipeline, "finemap", SimpleNamespace())
finemap_window_kb = getattr(finemap_opts, "window_kb", 500)
finemap_max_causal = getattr(finemap_opts, "max_causal", 10)
finemap_coverage = getattr(finemap_opts, "coverage", 0.95)
finemap_min_abs_corr = getattr(finemap_opts, "min_abs_corr", 0.5)

std_file_pattern = standardised_gwas_name("{prefix}")
lbf_pattern = RESULTS_DIR + "finemap/{prefix}_lbf.tsv"
credible_set_pattern = RESULTS_DIR + "finemap/{prefix}_credible_sets.tsv.gz"

for g in pipeline.gwases:
    g.lbf_file = RESULTS_DIR + "finemap/" + file_prefix(g.file) + "_lbf.tsv"
    g.credible_set_file = RESULTS_DIR + "finemap/" + file_prefix(g.file) + "_credible_sets.tsv.gz"

rule all:
    input:
        expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
        [g.lbf_file for g in pipeline.gwases],
        [g.credible_set_file for g in pipeline.gwases]

include: "rules/standardise_rule.smk"
include: "rules/clumping_rule.smk"

rule run_finemapping:
    resources:
        mem = "16G",
        time = "04:00:00"
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
        lbf = lbf_pattern,
        credible_sets = credible_set_pattern
    shell:
        """
        Rscript run_finemap.R \
            --gwas_filename {input.gwas} \
            --clumped_filename {input.clumped_file} \
            --ancestry {params.ancestry} \
            --N {params.n} \
            --output_lbf_file {output.lbf} \
            --output_credible_set_file {output.credible_sets} \
            --window_kb {params.window_kb} \
            --max_causal {params.max_causal} \
            --coverage {params.coverage} \
            --min_abs_corr {params.min_abs_corr}
        """

onsuccess:
    files = [g.lbf_file for g in pipeline.gwases] + [g.credible_set_file for g in pipeline.gwases]
    onsuccess(pipeline_name, files, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
