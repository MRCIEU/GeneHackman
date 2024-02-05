include: "common.smk"
singularity: docker_container

pipeline = parse_pipeline_input("input.json")
additional_mandatory_columns = ["N"]

onstart:
    print("##### QTL MR Pipeline #####")

pipeline.gwas.columns = resolve_gwas_columns(pipeline.gwas.file, pipeline.gwas.columns, additional_mandatory_columns)
pipeline.gwas.standardised_gwas = standardised_gwas_name(pipeline.gwas.file)
qtl_name = pipeline.qtl.dataset + "_" + pipeline.qtl.subcategory
exposures_string = join(" ", pipeline.qtl.exposures)

mr_results = RESULTS_DIR + "mr/" + file_prefix(pipeline.gwas.file) + "_" + qtl_name + ".tsv.gz"
coloc_results = RESULTS_DIR + "mr/coloc_" + file_prefix(pipeline.gwas.file) + "_" + qtl_name + ".tsv"
volcano_plot = RESULTS_DIR + "plots/volcano_plot" + file_prefix(pipeline.gwas.file) + "_" + qtl_name + ".png"
results_file = RESULTS_DIR + "mr/result_" + qtl_name + "_" + file_prefix(pipeline.gwas.file) + ".html"


rule all:
    input: mr_results, volcano_plot, coloc_results, results_file

rule standardise_gwases:
    threads: 4
    resources:
        mem = "8G"
    input: pipeline.gwas.file
    output: pipeline.gwas.standardised_gwas
    shell:
        """
        Rscript standardise_gwas.R \
            --input_gwas {input} \
            --output_gwas {output} \
            --input_columns {pipeline.gwas.columns} \
            --populate_rsid NO
        """


rule run_mr_against_qtl_datasets:
    input:
        gwas = pipeline.gwas.standardised_gwas,
        ancestry = pipeline.gwas.ancestry
    output: mr_results
    shell:
        """
        Rscript run_mr_against_qtl_data.R --gwas_filename {input.gwas} \
            --ancestry {input.ancestry} \
            --dataset {pipeline.qtl.dataset} \
            --exposures {exposures_string} \
            --subcategory {pipeline.qtl.subcategory} \
            --output {output}
        """

rule create_volcano_plot_of_mr_results:
    input: mr_results
    output: volcano_plot
    shell:
        """
        Rscript volcano.R --results_filename {input} \
            --title "Volcano Plot of MR results against {pipeline.qtl.dataset}" \
            --output_file {output}
        """

rule run_coloc_analysis_of_significant_mr_results:
    input:
        gwas = pipeline.gwas.file,
        mr_results = mr_results
    output: coloc_results
    shell:
        """
        Rscript coloc_of_mr_results.R \
            --mr_results_filename {input.mr_results} \
            --gwas_filename {input.gwas} \
            --qtl_dataset {pipeline.qtl.dataset} \
            --exposures {exposures_string} \
            --output_file {output} 
        """

files_created = {
    "gwas": pipeline.gwas.standardised_gwas,
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
            --rmd_file /home/R/markdown/mr_for_qtls.rmd \
            --params {results_string} \
            --output_file {output}
        """

onsuccess:
    onsuccess(list(files_created.values()), results_file)

onerror:
    onerror_message()
