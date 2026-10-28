#!/usr/bin/env bash
# Report what this machine has for building the Pan-UKB LD pack.  Changes nothing.
#   bash reference_data/pan_ukb_ld/00_check_environment.sh
# Paste the output back so the build settings can be checked before starting.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"

say() { printf '%-28s %s\n' "$1" "$2"; }
have() { command -v "$1" >/dev/null 2>&1 && echo "yes ($(command -v "$1"))" || echo "no"; }

echo "== system"
say "host" "$(hostname)"
say "os" "$(uname -srm)"
say "cpus" "$( (nproc 2>/dev/null || sysctl -n hw.ncpu) 2>/dev/null)"
say "memory" "$( (free -g 2>/dev/null | awk '/Mem:/{print $2 " GB"}') || (sysctl -n hw.memsize | awk '{print $1/1e9 " GB"}') )"
mkdir -p "$WORK_DIR" 2>/dev/null && say "WORK_DIR" "$WORK_DIR (writable)" || say "WORK_DIR" "$WORK_DIR (NOT writable)"
say "free space under WORK_DIR" "$(df -Ph "$(dirname "$WORK_DIR")" | awk 'NR==2{print $4 " free of " $2}')"

echo "== software"
for tool in tmux python3 java conda mamba micromamba apptainer singularity Rscript liftOver plink1.9 aws curl git; do
  say "$tool" "$(have "$tool")"
done
say "java version" "$(java -version 2>&1 | head -1)"
say "python3 version" "$(python3 --version 2>&1)"
say "R GeneHackman deps" "$( (Rscript -e 'suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet=TRUE)); cat("ok")' 2>/dev/null) || echo "not usable with plain Rscript")"
say "GENEHACKMAN_DIR" "$GENEHACKMAN_DIR ($(git -C "$GENEHACKMAN_DIR" rev-parse --abbrev-ref HEAD 2>/dev/null))"
say "R_CMD" "$R_CMD"
say "PIPELINE_DATA_DIR" "${PIPELINE_DATA_DIR:-unset}"
for sif in "${PIPELINE_DATA_DIR:-/nonexistent}"/genomic_data/pipeline/*.sif; do
  [[ -f "$sif" ]] && say "GeneHackman container" "$sif"
done
[[ -x "$WORK_DIR/env/bin/python" ]] && say "build env" "$WORK_DIR/env (exists)" || say "build env" "not created yet (run 01_setup_env.sh)"

echo "== network"
url="https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.EUR.ldadj.bm/metadata.json"
say "Pan-UKB bucket (HTTPS)" "$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 "$url" || echo failed) for $url"
say "Pan-UKB listing" "$(curl -s --max-time 20 'https://pan-ukb-us-east-1.s3.amazonaws.com/?list-type=2&prefix=ld_release/UKBB.EUR&delimiter=/' \
  | grep -o '<Prefix>[^<]*</Prefix>' | sed 's/<[^>]*>//g' | tr '\n' ' ' || echo failed)"
say "conda-forge" "$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 https://conda.anaconda.org/conda-forge/channeldata.json || echo failed)"
say "PyPI" "$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 https://pypi.org/simple/hail/ || echo failed)"
say "UCSC downloads" "$(curl -s -o /dev/null -w '%{http_code}' --max-time 20 https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/ || echo failed)"
