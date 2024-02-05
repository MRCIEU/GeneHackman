#!/bin/bash

if [[ $# -ne 3 ]] ; then
  echo """
  You gotta give 3 arguments
  1. vcf file name
  2. Comma separated list of header names
  3. Output filename
  """
  exit 0
fi

VCF_FILE=$1
HEADERS=$2
TSV_FILE=$3

QUERY_STRING=""
for header in ${HEADERS//,/ } ; do
  query="[%$header]\t"
  #result=$(bcftools query $VCF_FILE --format "%${header}\n" 2>&1 | head -n 1)
  #if [[ $result =~ "Error: no such tag defined" ]]; then
  #  modified_header="[%$header]\t"
  #fi
  QUERY_STRING="${QUERY_STRING}${query}"
done

QUERY_STRING="${QUERY_STRING:0:-2}\n"

echo $HEADERS | sed 's/,/\t/g' > $TSV_FILE
/home/bcftools/bcftools query $VCF_FILE --format "$QUERY_STRING" >> $TSV_FILE

