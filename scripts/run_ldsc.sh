#!/bin/bash
set -e
source /home/scripts/conda_init.sh
conda activate ldsc

if [[ $# -lt 4 || $# -gt 9 ]] ; then
  echo """
  Heritability for each GWAS and genetic correlation for every pair (not only the first against the rest,
  which is what a single ldsc --rg a,b,c does).  Arguments:
  1. Comma separated list of (standardised) GWAS
  2. Comma separated list of their sample sizes, used when a GWAS has no N_EFF or N column
  3. Ancestry: the LD score folder (ex: EUR; CSA for South Asian with the Pan-UKB scores)
  4. output log file (all LDSC logs concatenated); tables are written beside it as <prefix>_h2.tsv, <prefix>_rg.tsv
  5. LD score directory holding <ancestry>/ and w_hm3.snplist (default: \$PIPELINE_DATA_DIR/LDSCORE/pan_ukb)
  6. Comma separated study types (continuous or binary; default continuous)
  7. Comma separated sample prevalences (cases / (cases + controls)) of binary traits, NA otherwise
  8. Comma separated population prevalences of binary traits, NA otherwise (for liability-scale h2)
  9. Comma separated labels for the GWAS (default: file names)
  """
  exit 0
fi

GWASES=$1
NS=$2
ANCESTRY=$3
OUTPUT=$4
LDSC_DIR=${5:-$PIPELINE_DATA_DIR/LDSCORE/pan_ukb}
LDSC_DIR=${LDSC_DIR%/}
LDSC_DATA=$DATA_DIR/ldsc
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)

IFS=',' read -r -a GWAS_LIST <<< "$GWASES"
IFS=',' read -r -a N_LIST <<< "$NS"
N_GWAS=${#GWAS_LIST[@]}
default_list() { local value=$1; local out=(); for ((k = 0; k < N_GWAS; k++)); do out+=("$value"); done; local IFS=','; echo "${out[*]}"; }
IFS=',' read -r -a TYPE_LIST <<< "${6:-$(default_list continuous)}"
IFS=',' read -r -a SAMP_PREV_LIST <<< "${7:-$(default_list NA)}"
IFS=',' read -r -a POP_PREV_LIST <<< "${8:-$(default_list NA)}"
IFS=',' read -r -a LABEL_LIST <<< "${9:-}"
OUTPUT_PREFIX=${OUTPUT%.*}

for list_name in N_LIST TYPE_LIST SAMP_PREV_LIST POP_PREV_LIST; do
  declare -n list=$list_name
  if [[ ${#list[@]} != "$N_GWAS" ]]; then
    echo "Error: $list_name has ${#list[@]} entries, but there are $N_GWAS GWAS"
    exit 1
  fi
done

mkdir -p "$LDSC_DATA"
mkdir -p "$(dirname "$OUTPUT_PREFIX")"

# LD scores either per chromosome (<ancestry>/1.l2.ldscore.gz ...) or genome-wide in one file
# (<ancestry>/ld.l2.ldscore.gz, as reference_data/ldscores_pan_ukb writes the Pan-UKB scores)
if [[ -f "$LDSC_DIR/$ANCESTRY/ld.l2.ldscore.gz" ]]; then
  LD_SCORE_ARGS=(--ref-ld "$LDSC_DIR/$ANCESTRY/ld" --w-ld "$LDSC_DIR/$ANCESTRY/ld")
else
  LD_SCORE_ARGS=(--ref-ld-chr "$LDSC_DIR/$ANCESTRY/" --w-ld-chr "$LDSC_DIR/$ANCESTRY/")
fi

SUMSTATS_FILES=()
SAMP=()   # sample prevalence passed to LDSC for each GWAS (nan: observed scale)
POP=()
LABELS=()
for ((i = 0; i < N_GWAS; i++)); do
  gwas=${GWAS_LIST[i]}
  file_basename=$(basename -- "$gwas")
  file_prefix=${file_basename%%.*}
  LABELS+=("${LABEL_LIST[i]:-$file_prefix}")
  header=$(zcat -f "$gwas" | head -1 | tr '\t' '\n')

  # sample size: per-variant effective N (binary traits), else per-variant N, else the given N
  if grep -qx "N_EFF" <<< "$header"; then
    N_ARGS=(--N-col N_EFF)
    n_is_effective=true
  elif grep -qx "N" <<< "$header"; then
    N_ARGS=(--N-col N)
    n_is_effective=false
  else
    N_ARGS=(--N "${N_LIST[i]}")
    n_is_effective=false
  fi

  # liability scale for binary traits: with effective N the sample prevalence is 0.5 (Grotzinger et al. 2023)
  samp=nan; pop=nan
  if [[ ${TYPE_LIST[i]} == "binary" && ${POP_PREV_LIST[i]} != "NA" ]]; then
    if [[ $n_is_effective == true ]]; then
      samp=0.5
    elif [[ ${SAMP_PREV_LIST[i]} != "NA" ]]; then
      samp=${SAMP_PREV_LIST[i]}
    fi
    [[ $samp != nan ]] && pop=${POP_PREV_LIST[i]}
  fi
  SAMP+=("$samp")
  POP+=("$pop")

  python2.7 /home/ldsc/munge_sumstats.py \
      --sumstats "${gwas}" \
      "${N_ARGS[@]}" \
      --merge-alleles "${LDSC_DIR}/w_hm3.snplist" \
      --snp RSID --a1 EA --a2 OA --p P --signed-sumstats BETA,0 --frq EAF \
      --chunksize 500000 \
      --out "$LDSC_DATA/${file_prefix}"
  SUMSTATS_FILES+=("${LDSC_DATA}/${file_prefix}.sumstats.gz")
done

join_by_comma() { local IFS=','; echo "$*"; }
LOGS=()

# heritability of each GWAS, on the liability scale when its prevalences are known
for ((i = 0; i < N_GWAS; i++)); do
  PREV_ARGS=()
  [[ ${SAMP[i]} != nan ]] && PREV_ARGS=(--samp-prev "${SAMP[i]}" --pop-prev "${POP[i]}")
  python2.7 /home/ldsc/ldsc.py \
    --h2 "${SUMSTATS_FILES[i]}" \
    "${LD_SCORE_ARGS[@]}" \
    "${PREV_ARGS[@]}" \
    --out "${OUTPUT_PREFIX}_h2_$((i + 1))"
  LOGS+=("${OUTPUT_PREFIX}_h2_$((i + 1)).log")
done

# genetic correlation of every pair: GWAS i against GWAS i+1..n, for each i
for ((i = 0; i < N_GWAS - 1; i++)); do
  PREV_ARGS=()
  if [[ " ${SAMP[*]} " =~ [0-9] ]]; then
    PREV_ARGS=(--samp-prev "$(join_by_comma "${SAMP[@]:i}")" --pop-prev "$(join_by_comma "${POP[@]:i}")")
  fi
  python2.7 /home/ldsc/ldsc.py \
    --rg "$(join_by_comma "${SUMSTATS_FILES[@]:i}")" \
    "${LD_SCORE_ARGS[@]}" \
    "${PREV_ARGS[@]}" \
    --out "${OUTPUT_PREFIX}_rg_$((i + 1))"
  LOGS+=("${OUTPUT_PREFIX}_rg_$((i + 1)).log")
done

cat "${LOGS[@]}" > "$OUTPUT"

cd "$SCRIPT_DIR"   # parse_ldsc_results.R sources load.R from the scripts folder
Rscript parse_ldsc_results.R \
  --logs "$(join_by_comma "${LOGS[@]}")" \
  --sumstats "$(join_by_comma "${SUMSTATS_FILES[@]}")" \
  --labels "$(join_by_comma "${LABELS[@]}")" \
  --population_prevalences "$(join_by_comma "${POP[@]}")" \
  --ancestry "$ANCESTRY" \
  --h2_output "${OUTPUT_PREFIX}_h2.tsv" \
  --rg_output "${OUTPUT_PREFIX}_rg.tsv"
