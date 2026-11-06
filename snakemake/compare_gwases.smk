import pathlib

include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "compare_gwases"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)

onstart:
    print("##### GWAS Comparison Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries, pipeline.ld_reference)
# LDSC LD scores: ldsc.ld_scores is pan_ukb (default) or 1000genomes_b37; matched on rsID, so any output build
if not hasattr(pipeline, "ldsc"): pipeline.ldsc = SimpleNamespace()
pipeline.ldsc.ld_scores = str(getattr(pipeline.ldsc, "ld_scores", default_ld_scores)).strip()
if pipeline.ldsc.ld_scores not in LDSC_DIRS:
    stop_with_message(f"ldsc.ld_scores must be one of {', '.join(LDSC_DIRS)}, not {pipeline.ldsc.ld_scores}.")
if not all(a in ld_scores_ancestries[pipeline.ldsc.ld_scores] for a in ancestries):
    stop_with_message(f"With ldsc.ld_scores: {pipeline.ldsc.ld_scores}, each GWAS's ancestry must be one of: "
                      f"{', '.join(ld_scores_ancestries[pipeline.ldsc.ld_scores])}.")
check_reference_data(pipeline, pipeline_includes_clumping=True, ld_scores=pipeline.ldsc.ld_scores,
                     ld_score_ancestries=ancestries)


def ld_score_population(ancestry):
    """Folder of the LD score set for an ancestry (the Pan-UKB scores call South Asian CSA)."""
    if pipeline.ldsc.ld_scores == "pan_ukb":
        return pan_ukb_ancestry_aliases.get(ancestry, ancestry)
    return ancestry


#List of output files
expected_vs_observed_results = RESULTS_DIR + "gwas_comparison/expected_vs_observed_outcomes.tsv"
expected_vs_observed_variants = RESULTS_DIR + "gwas_comparison/expected_vs_observed_variants.tsv"
heterogeneity_scores = RESULTS_DIR + "gwas_comparison/heterogeneity_scores.tsv"
heterogeneity_plot = RESULTS_DIR + "plots/snp_heterogeneity_plot.png"
heterogeneity_snp_comparison = RESULTS_DIR + "plots/snps_with_heterogeneity_forest_plot.png"
results_file = RESULTS_DIR + "gwas_comparison/result_compare_gwases.html"

std_file_pattern = standardised_gwas_name("{prefix}")
ldsc_result_pattern = RESULTS_DIR + "ldsc/results_{ancestry}.log"
ldsc_directory = RESULTS_DIR + "ldsc/"
pathlib.Path(ldsc_directory).mkdir(parents=True, exist_ok=True)

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
           expand(ldsc_result_pattern, ancestry=ancestries),
           expected_vs_observed_results,
           expected_vs_observed_variants,
           heterogeneity_scores,
           heterogeneity_plot,
           heterogeneity_snp_comparison,
           results_file

include: "rules/standardise_rule.smk"
include: "rules/clumping_rule.smk"

rule compare_observed_vs_expected_gwas:
    resources:
        mem = f"{len(pipeline.gwases)*16}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_file for g in pipeline.gwases]
    output:
        results = expected_vs_observed_results,
        variants = expected_vs_observed_variants
    shell:
        """
        Rscript compare_observed_vs_expected_gwas.R  \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --result_output {output.results} \
            --variants_output {output.variants}
        """


rule heterogeneity_between_gwases:
    resources:
        mem = f"{len(pipeline.gwases)*12}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_file for g in pipeline.gwases]
    output:
        heterogeneity_scores = heterogeneity_scores,
        heterogeneity_plot = heterogeneity_plot,
        heterogeneity_snp_comparison  = heterogeneity_snp_comparison
    shell:
        """
        Rscript calculate_heterogeneity.R \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --ancestry_list {ancestries} \
            --heterogeneity_scores_output {output.heterogeneity_scores} \
            --heterogeneity_plot_output {output.heterogeneity_plot} \
            --heterogeneity_plot_per_snp_output {output.heterogeneity_snp_comparison}
        """


rule calculate_ldsc_and_genetic_correlation:
    resources:
        mem = "16G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases]
    params:
        gwases = lambda wildcards: ','.join([g.standardised_gwas for g in pipeline.gwases if g.ancestry == wildcards.ancestry]),
        ns = lambda wildcards: ','.join([str(g.N) for g in pipeline.gwases if g.ancestry == wildcards.ancestry]),
        ancestry = lambda wildcards: ld_score_population(wildcards.ancestry),
        ld_score_dir = LDSC_DIRS[pipeline.ldsc.ld_scores]
    output: ldsc_result_pattern
    shell:
        """
        ./run_ldsc.sh {params.gwases} {params.ns} {params.ancestry} {output} {params.ld_score_dir}
        """


files_created = {
    "results": expected_vs_observed_results,
    "variants": expected_vs_observed_variants,
    "heterogeneity_scores": heterogeneity_scores,
    "heterogeneity_plot": heterogeneity_plot,
    "heterogeneity_snp_comparison": heterogeneity_snp_comparison,
    "ldsc_directory": ldsc_directory
}
results_string = turn_dict_into_cli_string(files_created)


rule create_results_file:
    threads: 4
    resources:
        mem = "8G",
    input: list(files_created.values()), expand(ldsc_result_pattern, ancestry=ancestries)
    output: results_file
    shell:
        """
        Rscript create_results_file.R \
            --rmd_file /home/R/markdown/compare_gwases.Rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(pipeline_name, list(files_created.values()), results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
