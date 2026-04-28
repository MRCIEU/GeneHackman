include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "coloc"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)

onstart:
    print("##### Colocalization Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)

finemap_opts = getattr(pipeline, "finemap", SimpleNamespace())
finemap_window_kb = getattr(finemap_opts, "window_kb", 1000)
finemap_max_causal = getattr(finemap_opts, "max_causal", 10)
finemap_coverage = getattr(finemap_opts, "coverage", 0.95)
finemap_min_abs_corr = getattr(finemap_opts, "min_abs_corr", 0.5)

coloc_opts = getattr(pipeline, "coloc", SimpleNamespace())
coloc_overlap_kb = getattr(coloc_opts, "overlap_kb", 1000)
coloc_p1 = getattr(coloc_opts, "p1", 1e-4)
coloc_p2 = getattr(coloc_opts, "p2", 1e-4)
coloc_p12 = getattr(coloc_opts, "p12", 5e-6)

std_file_pattern = standardised_gwas_name("{prefix}")
finemap_dir_pattern = RESULTS_DIR + "finemap/{prefix}"

for g in pipeline.gwases:
    g.finemap_dir = RESULTS_DIR + "finemap/" + g.prefix

trait_names = [g.prefix for g in pipeline.gwases]
finemap_dirs_str = " ".join([g.finemap_dir for g in pipeline.gwases])
trait_names_str = " ".join(trait_names)

coloc_results = RESULTS_DIR + "coloc/coloc_results.tsv"
results_file = RESULTS_DIR + "coloc/result_coloc.html"

trait_bar = "|".join(trait_names)
ancestry_bar = "|".join([g.ancestry for g in pipeline.gwases])

rule all:
    input:
        expand(std_file_pattern, prefix=trait_names),
        expand(finemap_dir_pattern, prefix=trait_names),
        coloc_results,
        results_file

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
        finemap_dir = directory(finemap_dir_pattern)
    shell:
        """
        Rscript run_finemap.R \
            --gwas_filename {input.gwas} \
            --clumped_filename {input.clumped_file} \
            --ancestry {params.ancestry} \
            --N {params.n} \
            --output_finemap_dir {output.finemap_dir} \
            --window_kb {params.window_kb} \
            --max_causal {params.max_causal} \
            --coverage {params.coverage} \
            --min_abs_corr {params.min_abs_corr}
        """

rule run_coloc:
    resources:
        mem = "16G",
        time = "02:00:00"
    input:
        finemap_dirs = expand(finemap_dir_pattern, prefix=trait_names)
    params:
        finemap_dirs = finemap_dirs_str,
        trait_names = trait_names_str,
        overlap_kb = coloc_overlap_kb,
        p1 = coloc_p1,
        p2 = coloc_p2,
        p12 = coloc_p12
    output:
        coloc_results = coloc_results
    shell:
        """
        Rscript run_coloc.R \
            --finemap_dirs {params.finemap_dirs} \
            --trait_names {params.trait_names} \
            --overlap_kb {params.overlap_kb} \
            --p1 {params.p1} \
            --p2 {params.p2} \
            --p12 {params.p12} \
            --output_file {output.coloc_results}
        """

rmd_report_params = {
    "coloc_results": coloc_results,
    "traits_ordered": trait_bar,
    "ancestries_ordered": ancestry_bar,
    "overlap_kb": str(coloc_overlap_kb),
    "p1": str(coloc_p1),
    "p2": str(coloc_p2),
    "p12": str(coloc_p12),
}
rmd_report_param_string = turn_dict_into_cli_string(rmd_report_params)

rule create_results_file:
    threads: 4
    resources:
        mem = "8G",
    input:
        coloc_results
    output: results_file
    shell:
        """
        Rscript create_results_file.R \
            --rmd_file /home/R/markdown/coloc.Rmd \
            --params '{rmd_report_param_string}' \
            --output_file "{output}"
        """

onsuccess:
    files = [g.finemap_dir for g in pipeline.gwases] + [coloc_results, results_file]
    onsuccess(pipeline_name, files, results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
