include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "disease_progression"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)
incident = pipeline.gwases[0]
subsequent = pipeline.gwases[1]

onstart:
    print("##### Pipeline to Calculate Slope and Apply Correction on Collider Bias #####")

ancestries = [incident.ancestry, subsequent.ancestry]
validate_ancestries(ancestries)

collider_bias_results = RESULTS_DIR + "collider_bias/" + subsequent.prefix + "_collider_bias_results.tsv"
harmonised_effects = RESULTS_DIR + "collider_bias/" + subsequent.prefix + "_harmonised_effects.tsv.gz"
slopehunter_results = RESULTS_DIR + "collider_bias/" + subsequent.prefix + "_slopehunter.tsv.gz"

unadjusted_miami_plot = RESULTS_DIR + "plots/" + subsequent.prefix + "_miami_plot.png"
slopehunter_adjusted_miami_plot = RESULTS_DIR + "plots/" + file_prefix(slopehunter_results) + "_miami_plot.png"
results_file = RESULTS_DIR + "collider_bias/result_" + incident.prefix + "_" + subsequent.prefix + ".html"

original_expected_vs_observed_results = RESULTS_DIR + "collider_bias/original_expected_vs_observed_outcomes.tsv"
original_expected_vs_observed_variants = RESULTS_DIR + "collider_bias/original_expected_vs_observed_variants.tsv"
adjusted_expected_vs_observed_results = RESULTS_DIR + "collider_bias/adjusted_expected_vs_observed_outcomes.tsv"
adjusted_expected_vs_observed_variants = RESULTS_DIR + "collider_bias/adjusted_expected_vs_observed_variants.tsv"

std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[incident.prefix, subsequent.prefix]),
        collider_bias_results, slopehunter_results, harmonised_effects, unadjusted_miami_plot,
        slopehunter_adjusted_miami_plot, original_expected_vs_observed_results, original_expected_vs_observed_variants,
        adjusted_expected_vs_observed_results, adjusted_expected_vs_observed_variants, results_file

include: "rules/standardise_rule.smk"
include: "rules/clumping_rule.smk"


rule collider_bias_correction:
    threads: 8
    resources:
        mem = "48G"
    input:
        incidence_gwas = incident.standardised_gwas,
        subsequent_gwas = subsequent.standardised_gwas,
        clumped_file = incident.clumped_file
    output:
        results = collider_bias_results,
        slophunter_adjusted = slopehunter_results,
        harmonised_effects_results_file = harmonised_effects
    shell:
        """
        Rscript correct_for_collider_bias.R \
            --incidence_gwas {input.incidence_gwas} \
            --subsequent_gwas {input.subsequent_gwas} \
            --clumped_file {input.clumped_file} \
            --collider_bias_results_output {output.results} \
            --harmonised_effects_output {output.harmonised_effects_results_file} \
            --collider_bias_slopehunter_output {output.slophunter_adjusted}
        """

rule unadjusted_miami_plot:
    threads: 4
    resources:
        mem = "16G",
        time = "02:00:00"
    input:
        first_gwas = incident.standardised_gwas,
        second_gwas = subsequent.standardised_gwas,
    output: unadjusted_miami_plot
    shell:
        """
        Rscript miami.R \
            --first_gwas {input.first_gwas} \
            --second_gwas {input.second_gwas} \
            --miami_filename {output} \
            --title "Comparing Incidence and Subsequent GWAS"
        """

rule slopehunter_adjusted_miami_plot:
    threads: 4
    resources:
        mem = "16G",
        time = "02:00:00"
    input:
        first_gwas = incident.standardised_gwas,
        second_gwas = slopehunter_results
    output: slopehunter_adjusted_miami_plot
    shell:
        """
        Rscript miami.R \
            --first_gwas {input.first_gwas} \
            --second_gwas {input.second_gwas} \
            --miami_filename {output} \
            --title "Comparing Incidence and SlopeHunter Adjusted Subsequent GWAS"
        """

rule compare_original_observed_vs_expected_gwas:
    resources:
        mem = "32G"
    input:
        gwases = [incident.standardised_gwas, subsequent.standardised_gwas],
        clumped_files = [incident.clumped_file]
    output:
        results = original_expected_vs_observed_results,
        variants = original_expected_vs_observed_variants
    shell:
        """
        Rscript compare_observed_vs_expected_gwas.R  \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --result_output {output.results} \
            --variants_output {output.variants}
        """

rule compare_adjusted_observed_vs_expected_gwas:
    resources:
        mem = "32G"
    input:
        gwases = [incident.standardised_gwas, slopehunter_results],
        clumped_files = [incident.clumped_file]
    output:
        results = adjusted_expected_vs_observed_results,
        variants = adjusted_expected_vs_observed_variants
    shell:
        """
        Rscript compare_observed_vs_expected_gwas.R  \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --result_output {output.results} \
            --variants_output {output.variants}
        """

files_created = {
    "incidence_gwas": incident.standardised_gwas,
    "subsequent_gwas": subsequent.standardised_gwas,
    "clumped_snps": incident.clumped_file,
    "clumped_subsequent": subsequent.clumped_file,
    "collider_bias_results": collider_bias_results,
    "harmonised_gwas": harmonised_effects,
    "slopehunter_results": slopehunter_results,
    "unadjuested_miami_plot": unadjusted_miami_plot,
    "slopehunter_adjusted_miami_plot": slopehunter_adjusted_miami_plot,
    "original_expected_vs_observed": original_expected_vs_observed_results,
    "adjusted_expected_vs_observed": adjusted_expected_vs_observed_results
}
results_string = turn_dict_into_cli_string(files_created)

rule create_results_file:
    threads: 4
    resources:
        mem = "16G",
    input: list(files_created.values())
    output: results_file
    shell:
        """
        Rscript create_results_file.R \
            --rmd_file /home/R/markdown/disease_progression.Rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(pipeline_name, list(files_created.values()), results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
