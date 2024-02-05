clump_dir = DATA_DIR + "clumped_snps/"
if not os.path.isdir(clump_dir):
    os.makedirs(clump_dir)

for g in pipeline.gwases:
    g.clumped_snp_prefix = clump_dir + file_prefix(g.file)
    g.clumped_file = g.clumped_snp_prefix + ".clumped"
    setattr(pipeline, g.prefix, g)

rule clump_gwases:
    resources:
        mem = "4G"
    input:
        gwases = [g.standardised_gwas for g in pipeline.gwases]
    output: [g.clumped_file for g in pipeline.gwases]
    params:
        clumped_snp_prefixes = list([g.clumped_snp_prefix for g in pipeline.gwases])
    shell:
        """
        gwases=({input.gwases})
        ancestries=({ancestries})
        clumped_snp_prefixes=({params.clumped_snp_prefixes})

        for i in "${{!gwases[@]}}"
        do
            ancestry=${{ancestries[$i]}}
            plink1.9 --bfile {THOUSAND_GENOMES_DIR}/$ancestry \
                --clump ${{gwases[$i]}} \
                --clump-snp-field RSID \
                {pipeline.plink_clump_arguments} \
                --out ${{clumped_snp_prefixes[$i]}} || echo "{default_clump_headers}" > {output}
        done
        """
