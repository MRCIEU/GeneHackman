#!/bin/bash

# Usage: ./align_alleles.sh <input_prefix> <output_prefix>

INPUT=$1
OUTPUT=$2

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: ./align_alleles.sh <input_prefix> <output_prefix>"
    exit 1
fi

echo "Step 1: Extracting alleles and determining alphabetical A1..."

# This awk command looks at Column 5 and 6 of the .bim file.
# It picks the alphabetically 'smaller' one to be the new A1.
awk '{
    if ($5 < $6) 
        print $2, $5; 
    else 
        print $2, $6; 
}' "${INPUT}.bim" > a1_list.txt

echo "Step 2: Running PLINK to update allele order..."

# --make-bed creates the new triplet
# --a1-allele forces the alleles in our list to be in the A1 position
plink --bfile "$INPUT" \
      --a1-allele a1_list.txt \
      --make-bed \
      --out "$OUTPUT"

# Cleanup temporary file
rm a1_list.txt

echo "Done! New files saved as ${OUTPUT}.bed, .bim, .fam"