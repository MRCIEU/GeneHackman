include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "finemap"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)

onstart:
    print("##### SuSiE Fine-Mapping Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries, pipeline.ld_reference)
check_reference_data(pipeline, pipeline_includes_clumping=True)
unique_ancestries = list(set(ancestries))
if len(unique_ancestries) > 1 and len(unique_ancestries) < len(ancestries):
    raise ValueError(
        "Finemap requires either the same ancestry for all GWAS inputs, "
        "or a distinct ancestry per GWAS (no duplicates). "
        f"Found ancestries: {ancestries}"
    )
is_multi_ancestry = len(unique_ancestries) > 1

# defaults are set in parse_pipeline_input (util/common.smk)
finemap_opts = pipeline.finemap
finemap_window_kb = finemap_opts.window_kb
finemap_max_causal = finemap_opts.max_causal
finemap_coverage = finemap_opts.coverage
finemap_min_abs_corr = finemap_opts.min_abs_corr
finemap_max_region_kb = finemap_opts.max_region_kb

std_file_pattern = standardised_gwas_name("{prefix}")

for g in pipeline.gwases:
    g.finemap_dir = RESULTS_DIR + "finemap/" + g.prefix

MULTI_FINEMAP_COMPLETE = multi_finemap_complete_file(pipeline.gwases)
MULTI_FINEMAP_DIR = RESULTS_DIR + "finemap/multi_ancestry"
results_file = RESULTS_DIR + "finemap/result_finemap_" + gwas_run_label(pipeline.gwases) + ".html"

if is_multi_ancestry:
    rule all:
        input:
            expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
            MULTI_FINEMAP_COMPLETE,
            results_file
else:
    rule all:
        input:
            expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
            expand(FINEMAP_COMPLETE_TXT_PATTERN, prefix=[g.prefix for g in pipeline.gwases]),
            results_file

include: "rules/standardise_rule.smk"
include: "rules/clumping_rule.smk"

if is_multi_ancestry:
    include: "rules/finemap_multi_ancestry_rule.smk"
else:
    include: "rules/finemap_rule.smk"

report_params = {
    "finemap_dirs": "none" if is_multi_ancestry else "|".join(g.finemap_dir for g in pipeline.gwases),
    "labels": "none" if is_multi_ancestry else "|".join(re.sub(r"[|,=\s]+", "_", g.label) for g in pipeline.gwases),
    "multi_ancestry_dir": MULTI_FINEMAP_DIR if is_multi_ancestry else "none",
}
report_param_string = turn_dict_into_cli_string(report_params)

rule create_results_file:
    threads: 4
    resources:
        mem = "8G",
    input:
        MULTI_FINEMAP_COMPLETE if is_multi_ancestry else expand(FINEMAP_COMPLETE_TXT_PATTERN, prefix=[g.prefix for g in pipeline.gwases])
    output: results_file
    shell:
        """
        Rscript create_results_file.R \
            --rmd_file /home/R/markdown/finemap.Rmd \
            --params '{report_param_string}' \
            --output_file "{output}"
        """

onsuccess:
    if is_multi_ancestry:
        files = [RESULTS_DIR + "finemap/multi_ancestry"]
    else:
        files = [g.finemap_dir for g in pipeline.gwases]
    onsuccess(pipeline_name, files, results_file, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
