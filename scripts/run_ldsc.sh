#!/bin/bash
source /home/scripts/conda_init.sh
conda activate ldsc

if [[ $# -ne 4 ]] ; then
  echo """
  You gotta give 4 arguments
  1. Comma separated list of (standardised) GWAS
  2. Comma separated list of Corresponding N (sample size) of GWAS
  3. Ancestry (ex: EUR)
  4. output file prefix
  """
  exit 0
fi

GWASES=$1
NS=$2
ANCESTRY=$3
OUTPUT=$4

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
      --snp RSID --a1 EA --a2 OA --p P --signed-sumstats BETA,0 \
      --chunksize 500000 \
      --out "$LDSC_DATA/${file_prefix}"

  SUMSTATS="${LDSC_DATA}/${file_prefix}.sumstats.gz"
  SUMSTATS_FILES+=($SUMSTATS)
  i=$((i+1))
done

SUMSTATS_STRING=${SUMSTATS_FILES[@]}
SUMSTATS_STRING=${SUMSTATS_STRING// /,}
SUMSTATS_STRING=$(echo $SUMSTATS_STRING| sed '0,/,/{s///}')

if [[ $GWASES =~ "," ]]; then
  python2.7 /home/ldsc/ldsc.py \
    --rg "$SUMSTATS_STRING" \
    --ref-ld-chr "$LDSC_DIR"/"$ANCESTRY"/ \
    --w-ld-chr "$LDSC_DIR"/"$ANCESTRY"/ \
    --out $OUTPUT_PREFIX

    ldsc_rg_log="$OUTPUT_PREFIX".log
    genetic_correlation=$(awk '/p1   /' RS= "$ldsc_rg_log" | tr -s ' ' ' ')
    echo "$genetic_correlation" > "${OUTPUT_PREFIX}".tsv
else
  python2.7 /home/ldsc/ldsc.py \
    --h2 "$SUMSTATS" \
    --ref-ld-chr "$LDSC_DIR"/"$ANCESTRY"/ \
    --w-ld-chr "$LDSC_DIR"/"$ANCESTRY"/ \
    --out $OUTPUT_PREFIX
fi

