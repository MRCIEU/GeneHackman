include: "util/common.smk"
container: get_docker_container()

pipeline_name = "qtl_mr"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)
gwas = pipeline.gwases[0]

onstart:
    print("##### QTL MR Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)

finemap_opts = getattr(pipeline, "finemap", SimpleNamespace())
finemap_window_kb = getattr(finemap_opts, "window_kb", 1000)
finemap_max_causal = getattr(finemap_opts, "max_causal", 10)
finemap_coverage = getattr(finemap_opts, "coverage", 0.95)
finemap_min_abs_corr = getattr(finemap_opts, "min_abs_corr", 0.5)

qtl_opts = getattr(pipeline, "qtl", SimpleNamespace())
study_type = getattr(qtl_opts, "study_type", "continuous")

qtl_name = pipeline.qtl.dataset + "_" + pipeline.qtl.subcategory
if not hasattr(pipeline.qtl, "exposures"): pipeline.qtl.exposures = []
exposures_string = " ".join(pipeline.qtl.exposures)

gwas_prefix = file_prefix(gwas.file)
gwas.finemap_dir = RESULTS_DIR + "finemap/" + gwas_prefix
finemap_dir_pattern = RESULTS_DIR + "finemap/{prefix}"

mr_results = RESULTS_DIR + "mr/" + gwas_prefix + "_" + qtl_name + ".tsv.gz"
coloc_results = RESULTS_DIR + "mr/coloc_" + gwas_prefix + "_" + qtl_name + ".tsv"
volcano_plot = RESULTS_DIR + "plots/volcano_plot_" + gwas_prefix  + "_" + qtl_name + ".png"
results_file = RESULTS_DIR + "mr/result_" + qtl_name + "_" + gwas_prefix  + ".html"


std_file_pattern = standardised_gwas_name("{prefix}")
rule all:
    input: gwas.standardised_gwas, gwas.finemap_dir, mr_results, volcano_plot, coloc_results, results_file

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

rule run_mr_against_qtl_datasets:
    threads: 4
    resources:
        mem = "24G"
    input:
        gwas = gwas.standardised_gwas
    params:
        ancestry = gwas.ancestry,
        exposures = f"--exposures {exposures_string}" if exposures_string else ""
    output: mr_results
    shell:
        """
        Rscript run_mr_against_qtl_data.R \
            --gwas_filename {input.gwas} \
            --ancestry {params.ancestry} \
            --dataset {pipeline.qtl.dataset} \
            --subcategory {pipeline.qtl.subcategory} {params.exposures} \
            --output_file {output}
        """

rule create_volcano_plot_of_mr_results:
    threads: 2
    resources:
        mem = "8G"
    input: mr_results
    output: volcano_plot
    shell:
        """
        Rscript volcano.R \
            --results_filename {input} \
            --title "Volcano Plot of MR results against {pipeline.qtl.dataset}" \
            --output_file {output}
        """

rule run_coloc_analysis_of_significant_mr_results:
    threads: 2
    resources:
        mem = "16G"
    input:
        mr_results = mr_results,
        finemap_dir = gwas.finemap_dir
    params:
        N = gwas.N,
        study_type = study_type,
        exposures = f"--exposures {exposures_string}" if exposures_string else ""
    output: coloc_results
    shell:
        """
        Rscript coloc_of_mr_results.R \
            --mr_results_filename {input.mr_results} \
            --finemap_dir {input.finemap_dir} \
            --N {params.N} \
            --study_type {params.study_type} \
            --dataset {pipeline.qtl.dataset} {params.exposures} \
            --output_file {output}
        """

files_created = {
    "gwas": gwas.standardised_gwas,
    "finemap_dir": gwas.finemap_dir,
    "mr_results": mr_results,
    "volcano_plot": volcano_plot,
    "coloc_results": coloc_results
}
results_string = turn_dict_into_cli_string(files_created)

rule create_results_file:
    threads: 4
    resources:
        mem = "8G",
    input: list(files_created.values())
    output: results_file
    shell:
        """
        Rscript create_results_file.R \
            --rmd_file /home/R/markdown/mr_for_qtls.Rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(pipeline_name, list(files_created.values()), results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
