#!/bin/bash
set -e
source /home/scripts/conda_init.sh
conda activate ldsc

if [[ $# -lt 4 || $# -gt 5 ]] ; then
  echo """
  You gotta give 4 or 5 arguments
  1. Comma separated list of (standardised) GWAS
  2. Comma separated list of Corresponding N (sample size) of GWAS
  3. Ancestry: the LD score folder (ex: EUR; CSA for South Asian with the Pan-UKB scores)
  4. output file prefix
  5. LD score directory holding <ancestry>/ and w_hm3.snplist (default: \$PIPELINE_DATA_DIR/LDSCORE/pan_ukb)
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

N_LIST=(${NS//,/ })
GWAS_LIST=(${GWASES//,/ })
OUTPUT_PREFIX=${OUTPUT/.*/}

if [[ ${#N_LIST[@]} != ${#GWAS_LIST[@]}  ]]; then
  echo "Error: GWASES and NS must be the same size"
  exit 1
fi

mkdir -p "$LDSC_DATA"
mkdir -p $(dirname $OUTPUT_PREFIX)

SUMSTATS_FILES=""
i=0
for gwas in ${GWASES//,/ } ; do
  n=${N_LIST[i]}
  file_basename=$(basename -- "$gwas")
  file_prefix=${file_basename/.*/}

  python2.7 /home/ldsc/munge_sumstats.py \
      --sumstats "${gwas}" \
      --N "$n" \
      --merge-alleles "${LDSC_DIR}/w_hm3.snplist" \
      --snp RSID --a1 EA --a2 OA --p P --signed-sumstats BETA,0 --frq EAF \
      --chunksize 500000 \
      --out "$LDSC_DATA/${file_prefix}"

  SUMSTATS="${LDSC_DATA}/${file_prefix}.sumstats.gz"
  SUMSTATS_FILES+=($SUMSTATS)
  i=$((i+1))
done

# LD scores either per chromosome (<ancestry>/1.l2.ldscore.gz ...) or genome-wide in one file
# (<ancestry>/ld.l2.ldscore.gz, as reference_data/ldscores_pan_ukb writes the Pan-UKB scores)
if [[ -f "$LDSC_DIR/$ANCESTRY/ld.l2.ldscore.gz" ]]; then
  LD_SCORE_ARGS=(--ref-ld "$LDSC_DIR/$ANCESTRY/ld" --w-ld "$LDSC_DIR/$ANCESTRY/ld")
else
  LD_SCORE_ARGS=(--ref-ld-chr "$LDSC_DIR/$ANCESTRY/" --w-ld-chr "$LDSC_DIR/$ANCESTRY/")
fi

SUMSTATS_STRING=${SUMSTATS_FILES[@]}
SUMSTATS_STRING=${SUMSTATS_STRING// /,}
SUMSTATS_STRING=$(echo $SUMSTATS_STRING| sed '0,/,/{s///}')

if [[ $GWASES =~ "," ]]; then
  python2.7 /home/ldsc/ldsc.py \
    --rg "$SUMSTATS_STRING" \
    "${LD_SCORE_ARGS[@]}" \
    --out $OUTPUT_PREFIX

    ldsc_rg_log="$OUTPUT_PREFIX".log
else
  python2.7 /home/ldsc/ldsc.py \
    --h2 "$SUMSTATS" \
    "${LD_SCORE_ARGS[@]}" \
    --out $OUTPUT_PREFIX
fi

