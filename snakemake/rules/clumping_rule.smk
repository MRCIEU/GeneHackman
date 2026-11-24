clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
    os.makedirs(clump_dir)

for g in pipeline.gwases:
    g.clumped_snp_prefix = clump_dir + file_prefix(g.file)
    g.clumped_file = g.clumped_snp_prefix + ".clumped"
    setattr(pipeline, g.prefix, g)

CLUMPED_FILE_PATTERN = clump_dir + "{prefix}.clumped"


def clumping_memory(wildcards):
    """The LD pack and the 1000 Genomes panel need little memory; the user's own genotypes may be biobank-sized."""
    source = ld_source(getattr(pipeline, wildcards.prefix), "clumping")
    return {"pan_ukb": "16G", "1000genomes": "8G"}.get(source, "32G")


# Clumping with the LD reference of the clumping step (ld_reference; see snakemake/util/ld_sources.py): the Pan-UKB
# LD pack (PLINK 1.9 --clump rules, R/clumping.R), the 1000 Genomes panel (PLINK, matched on RSID) or the user's own
# genotypes (PLINK, matched on position and alleles).  The leading space keeps argparser from reading "--clump-..."
# as one of its own options.
rule clump_gwas:
    threads: 4
    resources:
        mem = clumping_memory
    input:
        gwas = lambda wildcards: getattr(pipeline, wildcards.prefix).standardised_gwas
    params:
        ancestry = lambda wildcards: getattr(pipeline, wildcards.prefix).ancestry,
        ld_reference = lambda wildcards: ld_source(getattr(pipeline, wildcards.prefix), "clumping"),
        build = pipeline.output.build,
        clump_arguments = pipeline.plink_clump_arguments
    output:
        clumped = CLUMPED_FILE_PATTERN
    shell:
        """
        Rscript clump_gwas.R \
            --gwas {input.gwas} \
            --ancestry '{params.ancestry}' \
            --ld_reference '{params.ld_reference}' \
            --build {params.build} \
            --clump_arguments " {params.clump_arguments}" \
            --threads {threads} \
            --output {output.clumped}
        """
