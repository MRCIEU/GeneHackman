#!/usr/bin/env bash
# Build the GeneHackman Pan-UKB LD pack.  Every step is resumable: re-running the same command skips
# what is finished and continues the rest, so a killed or rebooted run just needs restarting.
# Run it inside tmux (run_in_tmux.sh does that) because the full build takes days.
#
#   run_build.sh phase0           chr22 for every population + the dense test region + report
#   run_build.sh full [POP ...]   whole genome, one population at a time (default: all of $POPS)
#   run_build.sh status           progress, disk use and the last log lines
#   run_build.sh finalise         (re)write info.tsv, MANIFEST.json and checksums for finished work
#
# Settings: config.sh (override with environment variables or local_config.sh).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"

PY="${PYTHON:-$WORK_DIR/env/bin/python}"
if [[ -d "$WORK_DIR/env/lib/jvm" ]]; then export JAVA_HOME="$WORK_DIR/env/lib/jvm"; fi
export PATH="$WORK_DIR/bin:$WORK_DIR/env/bin:$PATH"
mkdir -p "$WORK_DIR"/{state,logs,phase0}
LOG="$WORK_DIR/logs/run_build_$(date +%Y%m%d).log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

settings_for() {  # settings a stage's output depends on; a change forces an explicit restart
  case "$1" in
    tables) echo "source=$SOURCE_URI" ;;
    variants) echo "maf=$MAF_MIN source=$SOURCE_URI" ;;
    annotate) echo "maf=$MAF_MIN" ;;
    *) echo "maf=$MAF_MIN radius=$RADIUS_KB deadzone=$DEADZONE chunk=$CHUNK_ROWS" ;;
  esac
}

# run_stage NAME POP KIND command...: run once; a marker records completion and the settings used
run_stage() {
  local name="$1" pop="$2" kind="$3"; shift 3
  local marker="$WORK_DIR/state/$pop.$name.done" settings
  settings="$(settings_for "$kind")"
  if [[ -f "$marker" ]]; then
    if [[ "$(cat "$marker")" != "$settings" ]]; then
      log "ERROR: $pop $name was done with different settings ($(cat "$marker")) than now ($settings)."
      log "       Remove $marker (and the outputs it covers) to redo it with the new settings."
      exit 1
    fi
    log "$pop $name: already done"
    return 0
  fi
  log "$pop $name: starting: $*"
  local started=$SECONDS
  "$@" 2>&1 | tee -a "$WORK_DIR/logs/${pop}_${name}.log"
  echo "$settings" > "$marker"
  log "$pop $name: done in $(( (SECONDS - started) / 60 )) min"
}

check_environment() {
  if [[ ! -x "$PY" ]]; then
    echo "Build environment missing ($PY). Run: bash $HERE/01_setup_env.sh"; exit 1
  fi
}

prepare_population() {  # source tables, kept variants and keys: small, kept between phases
  local pop="$1"
  run_stage tables "$pop" tables "$PY" "$HERE/05_fetch_source.py" tables --pop "$pop"
  run_stage variants "$pop" variants "$PY" "$HERE/10_extract_variants.py" --pop "$pop"
  # shellcheck disable=SC2086
  run_stage annotate "$pop" annotate $R_CMD "$HERE/20_annotate_variants.R" --pop "$pop"
  if [[ ! -f "$WORK_DIR/state/$pop.download_plan.json" ]]; then
    "$PY" "$HERE/05_fetch_source.py" plan --pop "$pop" 2>&1 | tee -a "$LOG"
  fi
}

clean_blocks() {
  local pop="$1"
  if [[ "$KEEP_SOURCE_BLOCKS" == "1" ]]; then
    log "$pop: keeping downloaded LD blocks (KEEP_SOURCE_BLOCKS=1)"
  else
    "$PY" "$HERE/05_fetch_source.py" clean-blocks --pop "$pop" 2>&1 | tee -a "$LOG"
  fi
}

phase0() {
  check_environment
  local chrom="$PHASE0_CHROM"
  for pop in $POPS; do
    prepare_population "$pop"
    run_stage "blocks_chr$chrom" "$pop" ld "$PY" "$HERE/05_fetch_source.py" blocks --pop "$pop" --chroms "$chrom"
    "$PY" "$HERE/30_extract_ld.py" band --pop "$pop" --chroms "$chrom" --measure 2>&1 | tee -a "$WORK_DIR/logs/${pop}_ld.log"
    run_stage "validate_chr$chrom" "$pop" ld "$PY" "$HERE/50_validate_pack.py" --pop "$pop" --chroms "$chrom"
    if [[ "$pop" == "$PHASE0_REGION_POP" && -n "$PHASE0_REGION" ]]; then
      run_stage region_blocks "$pop" ld "$PY" "$HERE/05_fetch_source.py" blocks --pop "$pop" \
        --region "$PHASE0_REGION" --radius-kb "$PHASE0_REGION_RADIUS_KB"
      run_stage region_dense "$pop" ld "$PY" "$HERE/30_extract_ld.py" dense --pop "$pop" \
        --region "$PHASE0_REGION" --out "phase0/${pop}_region.npz"
      if [[ -n "$PHASE0_GWAS" ]]; then
        # shellcheck disable=SC2086
        $R_CMD "$HERE/62_truncation_finemap.R" --region "phase0/${pop}_region" --gwas "$PHASE0_GWAS" \
          --radius-kb "$RADIUS_KB" 2>&1 | tee -a "$LOG" || log "WARNING: truncation fine-mapping check failed"
      fi
    fi
    clean_blocks "$pop"
  done
  "$PY" "$HERE/40_finalise_pack.py" --skip-checksums 2>&1 | tee -a "$LOG"
  "$PY" "$HERE/60_phase0_report.py" 2>&1 | tee -a "$LOG"
  log "Phase 0 finished. Report: $WORK_DIR/phase0/report.md"
}

remaining_chroms() {  # chromosomes with kept variants whose LD is not finished yet for a population
  local pop="$1" c out=""
  for c in $(seq 1 22); do
    [[ -f "$WORK_DIR/raw/$pop/variants_annot_chr$c.tsv.gz" ]] || continue
    [[ -f "$WORK_DIR/state/$pop.ld.chr$c.done" ]] || out="$out $c"
  done
  echo "$out"
}

full() {
  check_environment
  local pops
  if [[ $# -eq 0 ]]; then read -r -a pops <<< "$POPS"; else pops=("$@"); fi
  for pop in "${pops[@]}"; do
    prepare_population "$pop"
    if [[ ! -f "$WORK_DIR/state/$pop.full.done" ]]; then
      # chromosomes whose LD blocks this download covers (validation needs the blocks on disk)
      local chroms_file="$WORK_DIR/state/$pop.blocks_all.chroms" todo
      if [[ ! -f "$WORK_DIR/state/$pop.blocks_all.done" ]]; then
        remaining_chroms "$pop" > "$chroms_file"
      fi
      todo="$(cat "$chroms_file")"
      if [[ -n "${todo// /}" ]]; then
        # shellcheck disable=SC2086
        run_stage blocks_all "$pop" ld "$PY" "$HERE/05_fetch_source.py" blocks --pop "$pop" --chroms $todo
        # shellcheck disable=SC2086
        "$PY" "$HERE/30_extract_ld.py" band --pop "$pop" --chroms $todo 2>&1 | tee -a "$WORK_DIR/logs/${pop}_ld.log"
        # shellcheck disable=SC2086
        run_stage validate_all "$pop" ld "$PY" "$HERE/50_validate_pack.py" --pop "$pop" --chroms $todo
      fi
      settings_for ld > "$WORK_DIR/state/$pop.full.done"
    fi
    clean_blocks "$pop"
    "$PY" "$HERE/40_finalise_pack.py" --skip-checksums 2>&1 | tee -a "$LOG"
    log "$pop: complete ($(du -sh "$WORK_DIR/pack/$pop" | cut -f1) in the pack)"
  done
  "$PY" "$HERE/40_finalise_pack.py" 2>&1 | tee -a "$LOG"
  log "Build finished for: ${pops[*]}. Pack: $WORK_DIR/pack"
  log "Copy it to the laptop with, from the laptop:"
  log "  rsync -av --progress $(hostname):$WORK_DIR/pack/ /Users/wt23152/Documents/Projects/scratch/020/data/genomic_data/pan_ukb_ld/"
}

status() {
  echo "WORK_DIR=$WORK_DIR   POPS=$POPS   radius=${RADIUS_KB}kb maf>=$MAF_MIN deadzone=$DEADZONE"
  df -Ph "$WORK_DIR" | awk 'NR==2{print "disk: " $4 " free of " $2}'
  for pop in $POPS; do
    local done_chroms stages="" in_progress=""
    done_chroms=$(ls "$WORK_DIR/state/$pop".ld.chr*.done 2>/dev/null | wc -l | tr -d ' ')
    for s in tables variants annotate "blocks_chr$PHASE0_CHROM" blocks_all validate_all full; do
      [[ -f "$WORK_DIR/state/$pop.$s.done" ]] && stages="$stages $s"
    done
    for dir in "$WORK_DIR/work/$pop"/chunks_chr*; do
      [[ -d "$dir" ]] || continue
      local chrom="${dir##*chunks_chr}"
      [[ -f "$WORK_DIR/state/$pop.ld.chr$chrom.done" ]] && continue
      in_progress="$in_progress chr$chrom:$(ls "$dir"/*.json 2>/dev/null | grep -c chunk_)chunks"
    done
    local source_size pack_size
    source_size=$(du -sh "$WORK_DIR/source/$pop" 2>/dev/null | cut -f1)
    pack_size=$(du -sh "$WORK_DIR/pack/$pop" 2>/dev/null | cut -f1)
    printf '%-4s stages:%s | LD chromosomes done: %s/22%s | source %s | pack %s\n' \
      "$pop" "${stages:- none}" "$done_chroms" "${in_progress:+ | in progress:$in_progress}" \
      "${source_size:-0}" "${pack_size:-0}"
  done
  echo "--- last log lines"
  tail -n 5 "$(ls -t "$WORK_DIR"/logs/*.log 2>/dev/null | head -1)" 2>/dev/null || true
}

case "${1:-}" in
  phase0) phase0 ;;
  full) shift; full "$@" ;;
  status) status ;;
  finalise) check_environment; "$PY" "$HERE/40_finalise_pack.py" ;;
  *) sed -n '2,13p' "$0"; exit 1 ;;
esac
