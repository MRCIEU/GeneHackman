import pathlib

include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "drug_target_mr"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)

if (pipeline.genes != None or len(pipeline.genes) == 0):
    raise Exception("No genes provided")

onstart:
    print("##### Drug Target MR Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)

#List of output files
insturmental_variable_results = RESULTS_DIR + "drug_target_mr/instrumental_variable_results.tsv"
# ...

std_file_pattern = standardised_gwas_name("{prefix}")

rule all:
    input: expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
           insturmental_variable_results

include: "rules/standardise_rule.smk"
include: "rules/clumping_rule.smk"

rule extract_instrumental_variables:
    resources:
        mem = f"{len(pipeline.gwases)*16}G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases],
        clumped_files = [g.clumped_file for g in pipeline.gwases]
    output:
        results = insturmental_variable_results 
    shell:
        """
        Rscript extract_instrumental_variables.R  \
            --gwas_filenames {input.gwases} \
            --clumped_filenames {input.clumped_files} \
            --result_output {output.results}
        """


files_created = {
    "results": insturmental_variable_results ,
}
results_string = turn_dict_into_cli_string(files_created)


onsuccess:
    onsuccess(pipeline_name, list(files_created.values()), results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
