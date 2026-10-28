#!/usr/bin/env bash
# End-to-end test of the pack build on a synthetic Pan-UKB-like source (no network needed):
# builds the source with Hail, runs run_build.sh phase0 then full, interrupts and resumes a run,
# and checks the result against the truth.
#
#   PYTHON=<env python with hail> R_CMD=<Rscript able to load GeneHackman> LIFTOVER_BIN=<liftOver> \
#   CHAIN_FILE=<hg19ToHg38.over.chain.gz> bash reference_data/pan_ukb_ld/tests/run_synthetic_build.sh [test_dir]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DIR="${1:-$(mktemp -d)}"
mkdir -p "$TEST_DIR"
: "${PYTHON:?set PYTHON to a python with hail}"

export WORK_DIR="$TEST_DIR/work"
export SOURCE_URI="file://$TEST_DIR/source"
export POPS="EUR CSA" RADIUS_KB=40 DEADZONE=0.02 MAF_MIN=0.01 CHUNK_ROWS=24 FETCH_ROWS=60
export MIN_FREE_GB=0 HAIL_CORES=2 HAIL_DRIVER_MEMORY=2g DOWNLOAD_THREADS=4
export PHASE0_CHROM=22 PHASE0_REGION="22:20100000-20400000" PHASE0_REGION_POP=EUR PHASE0_REGION_RADIUS_KB=300
export GENEHACKMAN_DIR="$(cd "$HERE/../../.." && pwd)"

if [[ ! -d "$TEST_DIR/source/UKBB.CSA.ldadj.bm" ]]; then
  "$PYTHON" "$HERE/make_synthetic_source.py" "$TEST_DIR/source"
fi
mkdir -p "$WORK_DIR/liftover_data/genomic_data/liftover"
cp "$CHAIN_FILE" "$WORK_DIR/liftover_data/genomic_data/liftover/hg19ToHg38.over.chain.gz"

echo "=== phase 0"
"$HERE/../run_build.sh" phase0
test -s "$WORK_DIR/phase0/report.md"
test -s "$WORK_DIR/phase0/EUR_region.npz"

echo "=== full build, killed part-way through chr21 for EUR, then resumed"
"$PYTHON" "$HERE/../05_fetch_source.py" blocks --pop EUR --chroms 21
"$PYTHON" "$HERE/../30_extract_ld.py" band --pop EUR --chroms 21 --stop-after-chunks 3 || [[ $? -eq 3 ]]
n_before=$(ls "$WORK_DIR/work/EUR/chunks_chr21/"chunk_*.json | wc -l)
echo "chunks checkpointed before the simulated kill: $n_before"
[[ "$n_before" -eq 3 && ! -f "$WORK_DIR/state/EUR.ld.chr21.done" ]]
"$HERE/../run_build.sh" full 2>&1 | tee "$TEST_DIR/full.log"
# the resumed run must redo only the chunks that were not checkpointed
grep -E "EUR chr21: [0-9,]+ variants, [0-9]+ chunks, [0-9]+ to do" "$TEST_DIR/full.log" \
  | awk '{for (i = 1; i <= NF; i++) { if ($i == "chunks,") total = $(i - 1); if ($i == "to") todo = $(i - 1) } }
         END { print "resumed chr21: " todo " of " total " chunks redone"; exit !(todo == total - 3) }'
"$HERE/../run_build.sh" status

echo "=== verify against truth"
"$PYTHON" "$HERE/verify_synthetic_pack.py" "$TEST_DIR/source" "$WORK_DIR"
"$PYTHON" -c "import json; m = json.load(open('$WORK_DIR/pack/MANIFEST.json')); \
print('MANIFEST populations:', {p: (v['chromosomes'], v['n_variants'], v['n_samples']) for p, v in m['populations'].items()})"
echo "synthetic build test passed"
