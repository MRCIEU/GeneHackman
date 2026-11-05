# Settings for building GeneHackman's GRCh38 reference data (reference bundle 1.2.0).  Sourced by every
# script here; override any value by exporting it first or in local_config.sh (not tracked by git).

# Where downloads and intermediate files go
export WORK_DIR="${WORK_DIR:-/local-scratch/data/genehackman_grch38}"
# The reference data being assembled: genomic_data/ and LDSCORE/ are written under it
export BUNDLE_DIR="${BUNDLE_DIR:-$WORK_DIR/bundle_1.2.0}"
# An existing GeneHackman 1.1.0 reference data directory (GRCh37): source of w_hm3.snplist, the b37 panel
# for comparisons, and the chains
export OLD_PIPELINE_DATA_DIR="${OLD_PIPELINE_DATA_DIR:-}"
# dbSNP 156 GRCh38 .fst files (chr1.fst ... as read by R/chrpos_to_rsid.R), i.e. a b38_dbsnp156 directory
export DBSNP_B38_DIR="${DBSNP_B38_DIR:-}"
# The finished Pan-UKB LD pack (for the build check and the GRCh38 key concordance check)
export LD_PACK_DIR="${LD_PACK_DIR:-/local-scratch/data/pan_ukb_ld/pack}"
# QTL data (cis only): built from the genotype-phenotype map 1.0.0 study data into
# NEW_QTL_DATA_DIR/<dataset>/<build>/{top_hits,gwas}
export NEW_QTL_DATA_DIR="${NEW_QTL_DATA_DIR:-$WORK_DIR/qtl_1.2.0}"
export GPM_STUDIES="${GPM_STUDIES:-/local-scratch/projects/genotype-phenotype-map/results/1.0.0/studies_processed.tsv.gz}"
export GPM_GENE_RANGES="${GPM_GENE_RANGES:-/local-scratch/projects/genotype-phenotype-map/data/gene_ranges_hg38.tsv}"

# 1000 Genomes 30x on GRCh38 (NYGC; Byrska-Bishop et al. 2022): phased SNV/INDEL/SV panel
export KG30X_URL="${KG30X_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV}"
export KG30X_FILE_PATTERN="${KG30X_FILE_PATTERN:-1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz}"
# The 2,504 unrelated phase 3 samples with their populations (the same people as GeneHackman's GRCh37 panel)
export KG_PANEL_URL="${KG_PANEL_URL:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel}"
export SUPERPOPULATIONS="${SUPERPOPULATIONS:-EUR EAS AFR AMR SAS}"
export CHROMOSOMES="${CHROMOSOMES:-$(seq -s ' ' 1 22)}"
export MAF_MIN="${MAF_MIN:-0.01}"
export THREADS="${THREADS:-8}"

# Pan-UKB LD scores (HapMap3, LDSC format)
export PAN_UKB_LDSCORE_URL="${PAN_UKB_LDSCORE_URL:-https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.ALL.ldscore.tar.gz}"

# Tools: plink2, plink 1.9 and bcftools on the PATH (see environment.yml), and an Rscript that can load
# GeneHackman, e.g. the pipeline container:
#   R_CMD="apptainer exec --bind /local-scratch /path/genehackman_develop.sif Rscript"
export R_CMD="${R_CMD:-Rscript}"
export PLINK2="${PLINK2:-plink2}"
export PLINK19="${PLINK19:-plink}"
export GENEHACKMAN_DIR="${GENEHACKMAN_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"

if [[ -f "$(dirname "${BASH_SOURCE[0]}")/local_config.sh" ]]; then
  # shellcheck disable=SC1091
  source "$(dirname "${BASH_SOURCE[0]}")/local_config.sh"
fi
