include: "common.smk"
singularity: docker_container

output_file = RESULTS_DIR + "/test_output.log"
output_file_prefix = RESULTS_DIR + "/test_output.log"

onstart:
    print("Test pipeline for gwaspipeline package")

rule all:
    input: expand(RESULTS_DIR + "{prefix}.log", ["test1", "test2"])

rule is_everything_installed:
    input: expand("{hi}", ["Dockerfile", "Dockerfile_app"])
    output: RESULTS_DIR + "{prefix}.log"
    shell:
        """
        echo {hi}
        echo {{print({prefix})}}
        plink1.9 --version >> {output}
        """

