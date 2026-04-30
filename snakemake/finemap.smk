include: "util/common.smk"
singularity: get_docker_container()

pipeline_name = "finemap"
pipeline = parse_pipeline_input(pipeline_includes_clumping=True)

onstart:
    print("##### SuSiE Fine-Mapping Pipeline #####")

ancestries = list([g.ancestry for g in pipeline.gwases])
validate_ancestries(ancestries)

finemap_opts = getattr(pipeline, "finemap", SimpleNamespace())
finemap_window_kb = getattr(finemap_opts, "window_kb", 1000)
finemap_max_causal = getattr(finemap_opts, "max_causal", 10)
finemap_coverage = getattr(finemap_opts, "coverage", 0.95)
finemap_min_abs_corr = getattr(finemap_opts, "min_abs_corr", 0.5)

std_file_pattern = standardised_gwas_name("{prefix}")
finemap_dir_pattern = RESULTS_DIR + "finemap/{prefix}"

for g in pipeline.gwases:
    g.finemap_dir = RESULTS_DIR + "finemap/" + g.prefix

rule all:
    input:
        expand(std_file_pattern, prefix=[g.prefix for g in pipeline.gwases]),
        expand(finemap_dir_pattern, prefix=[g.prefix for g in pipeline.gwases])

include: "rules/standardise_rule.smk"
include: "rules/clumping_rule.smk"
include: "rules/finemap_rule.smk"

onsuccess:
    files = [g.finemap_dir for g in pipeline.gwases]
    onsuccess(pipeline_name, files, is_test=pipeline.is_test)

onerror:
    onerror_message(pipeline_name, is_test=pipeline.is_test)
