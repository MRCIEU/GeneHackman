# Settings for building the GeneHackman Pan-UKB LD pack.  Sourced by run_build.sh and the stage
# scripts; every value can be overridden by exporting it before running (or in local_config.sh,
# which is sourced last if present and is not tracked by git).

# Where everything is written: downloads, intermediate files, the finished pack, logs.
export WORK_DIR="${WORK_DIR:-/local-scratch/data/pan_ukb_ld}"

# Populations, processed one at a time (each one's downloaded LD blocks are deleted before the next).
export POPS="${POPS:-EUR AFR CSA EAS MID AMR}"

# Pack contents
export MAF_MIN="${MAF_MIN:-0.01}"          # keep variants with minor allele frequency >= this, per population
export RADIUS_KB="${RADIUS_KB:-1000}"      # store r for variant pairs at most this far apart (GRCh37)
export DEADZONE="${DEADZONE:-0.02}"        # store |r| below this as 0
export CHUNK_ROWS="${CHUNK_ROWS:-2048}"    # rows per compressed chunk in the pack
export FETCH_ROWS="${FETCH_ROWS:-8192}"    # rows fetched from the source matrix per Hail job

# Source: Pan-UKB LD release, public AWS bucket (anonymous access).  SOURCE_URI may also be a local
# directory laid out like the bucket (used by the tests).
export SOURCE_URI="${SOURCE_URI:-s3://pan-ukb-us-east-1/ld_release}"
export DOWNLOAD_THREADS="${DOWNLOAD_THREADS:-16}"
export MIN_FREE_GB="${MIN_FREE_GB:-500}"   # refuse to download if less than this would remain free
export KEEP_SOURCE_BLOCKS="${KEEP_SOURCE_BLOCKS:-0}"  # 1 = keep downloaded LD blocks after a population is done

# Phase 0 (feasibility run): chromosome to build in full, and a region extracted densely at a wide radius
export PHASE0_CHROM="${PHASE0_CHROM:-22}"
export PHASE0_REGION="${PHASE0_REGION:-6:35600000-37600000}"
export PHASE0_REGION_POP="${PHASE0_REGION_POP:-EUR}"
export PHASE0_REGION_RADIUS_KB="${PHASE0_REGION_RADIUS_KB:-2000}"
export PHASE0_GWAS="${PHASE0_GWAS:-}"      # optional standardised GWAS with a signal in PHASE0_REGION

# Software
export CONDA_ENV_NAME="${CONDA_ENV_NAME:-gh_ldpack}"
export GENEHACKMAN_DIR="${GENEHACKMAN_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
export HAIL_CORES="${HAIL_CORES:-8}"
export HAIL_DRIVER_MEMORY="${HAIL_DRIVER_MEMORY:-32g}"
# How to run R for the annotation stage: "Rscript" (needs the GeneHackman R dependencies) or, to use
# the pipeline container, e.g.  R_CMD="apptainer exec --bind /local-scratch /path/genehackman_1.1.0.sif Rscript"
export R_CMD="${R_CMD:-Rscript}"
export LIFTOVER_BIN="${LIFTOVER_BIN:-}"    # empty = liftOver on PATH, else $WORK_DIR/bin/liftOver (fetched by setup)
export CHAIN_FILE="${CHAIN_FILE:-}"        # empty = $WORK_DIR/liftover_data/genomic_data/liftover/hg19ToHg38.over.chain.gz

if [[ -f "$(dirname "${BASH_SOURCE[0]}")/local_config.sh" ]]; then
  # shellcheck disable=SC1091
  source "$(dirname "${BASH_SOURCE[0]}")/local_config.sh"
fi
