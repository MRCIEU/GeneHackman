clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
    os.makedirs(clump_dir)

for g in pipeline.gwases:
    g.clumped_snp_prefix = clump_dir + file_prefix(g.file)
    g.clumped_file = g.clumped_snp_prefix + ".clumped"
    setattr(pipeline, g.prefix, g)

CLUMPED_FILE_PATTERN = clump_dir + "{prefix}.clumped"

rule clump_gwas:
    threads: 4
    resources:
        mem = "8G"
    input:
        gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).standardised_gwas
    params:
        ancestry = lambda wildcards: getattr(pipeline, wildcards.prefix).ancestry,
        out_prefix = lambda wildcards: getattr(pipeline, wildcards.prefix).clumped_snp_prefix
    output:
        clumped = CLUMPED_FILE_PATTERN
    shell:
        """
        plink1.9 --bfile {THOUSAND_GENOMES_DIR}/{params.ancestry} \
            --threads {AVAILABLE_CPUS} \
            --clump {input.gwas} \
            --clump-snp-field RSID \
            {pipeline.plink_clump_arguments} \
            --out {params.out_prefix} || echo "{default_clump_headers}" > {output.clumped}
        """
