#!/usr/bin/env bash
# One-off setup for building the Pan-UKB LD pack:
#   - a Python/Java/Hail environment in $WORK_DIR/env (micromamba is fetched if no conda is available)
#   - liftOver and the hg19ToHg38 chain file (copied from $PIPELINE_DATA_DIR/genomic_data/liftover when
#     present, so the pack uses exactly the chain GeneHackman uses, else downloaded from UCSC)
#   - smoke tests of Hail and of the R side (GeneHackman's standardise_alleles / liftover)
# Safe to re-run.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"
mkdir -p "$WORK_DIR"/{bin,logs,liftover_data/genomic_data/liftover}

case "$(uname -s)-$(uname -m)" in
  Linux-x86_64) MAMBA_PLATFORM=linux-64; UCSC_PLATFORM=linux.x86_64 ;;
  Linux-aarch64) MAMBA_PLATFORM=linux-aarch64; UCSC_PLATFORM="" ;;
  Darwin-arm64) MAMBA_PLATFORM=osx-arm64; UCSC_PLATFORM=macOSX.arm64 ;;
  Darwin-x86_64) MAMBA_PLATFORM=osx-64; UCSC_PLATFORM=macOSX.x86_64 ;;
  *) echo "Unsupported platform $(uname -sm)"; exit 1 ;;
esac

# --- environment
if [[ ! -x "$WORK_DIR/env/bin/python" ]]; then
  if command -v micromamba >/dev/null; then MAMBA=micromamba
  elif command -v mamba >/dev/null; then MAMBA=mamba
  elif command -v conda >/dev/null; then MAMBA=conda
  else
    echo "No conda found; fetching micromamba into $WORK_DIR/bin"
    curl -Ls "https://micro.mamba.pm/api/micromamba/$MAMBA_PLATFORM/latest" | tar -xj -C "$WORK_DIR" bin/micromamba
    MAMBA="$WORK_DIR/bin/micromamba"
    export MAMBA_ROOT_PREFIX="$WORK_DIR/mamba"
  fi
  echo "Creating $WORK_DIR/env with $MAMBA"
  "$MAMBA" env create -y -p "$WORK_DIR/env" -f "$HERE/environment.yml" 2>/dev/null \
    || "$MAMBA" env create -p "$WORK_DIR/env" -f "$HERE/environment.yml"
fi
export JAVA_HOME="$WORK_DIR/env/lib/jvm"
export PATH="$WORK_DIR/env/bin:$PATH"

# --- liftOver (the environment provides bioconda's ucsc-liftover; download only if that failed)
if [[ -z "$LIFTOVER_BIN" ]] && ! command -v liftOver >/dev/null && [[ ! -x "$WORK_DIR/bin/liftOver" ]] \
   && [[ -n "$UCSC_PLATFORM" ]]; then
  curl -fsSL "https://hgdownload.soe.ucsc.edu/admin/exe/$UCSC_PLATFORM/liftOver" -o "$WORK_DIR/bin/liftOver"
  chmod +x "$WORK_DIR/bin/liftOver"
fi

# --- chain file (GeneHackman's copy when available)
chain="$WORK_DIR/liftover_data/genomic_data/liftover/hg19ToHg38.over.chain.gz"
if [[ ! -s "$chain" ]]; then
  if [[ -n "$CHAIN_FILE" ]]; then cp "$CHAIN_FILE" "$chain"
  elif [[ -n "${PIPELINE_DATA_DIR:-}" && -s "$PIPELINE_DATA_DIR/genomic_data/liftover/hg19ToHg38.over.chain.gz" ]]; then
    cp "$PIPELINE_DATA_DIR/genomic_data/liftover/hg19ToHg38.over.chain.gz" "$chain"
  else
    curl -fsSL https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -o "$chain"
  fi
fi

# --- smoke tests
echo "Hail:"
"$WORK_DIR/env/bin/python" - <<EOF
import os, hail as hl
hl.init(master="local[2]", quiet=True, tmp_dir="$WORK_DIR/tmp", log="$WORK_DIR/logs/hail_setup.log")
print("  hail", hl.version(), "range_table count =", hl.utils.range_table(10).count())
EOF
echo "Anonymous S3 access:"
"$WORK_DIR/env/bin/python" - <<EOF
import sys
sys.path.insert(0, "$HERE")
from build_common import open_store
store = open_store("$SOURCE_URI")
try:
    print("  metadata.json bytes:", len(store.read_bytes("UKBB.EUR.ldadj.bm/metadata.json")))
except Exception as error:
    print("  could not read the Pan-UKB bucket:", error)
EOF
echo "R side (GeneHackman load + liftOver):"
# shellcheck disable=SC2086
$R_CMD -e 'suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE)); cat("  GeneHackman loads; standardise_alleles available:", exists("standardise_alleles"), "\n")' \
  || echo "  R_CMD cannot load GeneHackman: set R_CMD to use the pipeline container (see config.sh)"
ls -l "$chain"
command -v liftOver || ls -l "$WORK_DIR/bin/liftOver" || echo "liftOver not found: set LIFTOVER_BIN"
echo "Setup done."
