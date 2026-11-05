#!/usr/bin/env bash
# Build GeneHackman's GRCh38 1000 Genomes panel from the 1000G 30x release (NYGC, native GRCh38):
#
#   bash reference_data/grch38/10_build_1000g_panel.sh
#
# Source (KG30X_SOURCE in config.sh): PLINK 2's prebuilt all_hg38 fileset of that release (default), or the
# per-chromosome phased panel VCFs.  The PLINK 2 fileset has multiallelic sites split into biallelic records, which
# --max-alleles 2 keeps; 11_assign_rsids.R drops positions whose records end up sharing an rsID.
#
# Same rules as the GRCh37 phase 3 panel it replaces (genomic_data/1000genomes/b37_dbsnp156):
#   - the 2,504 unrelated phase 3 people, split by superpopulation (EUR, EAS, AFR, AMR, SAS);
#   - autosomes, biallelic SNVs, MAF >= MAF_MIN in each superpopulation;
#   - variant IDs are dbSNP 156 rsIDs (variants without one, or sharing one, are dropped);
#   - A1 is the alphabetically smaller allele, as GeneHackman's standardised GWAS (fine-mapping signs rely on it).
# Writes to $BUNDLE_DIR/genomic_data/1000genomes/b38_dbsnp156/: <POP>.{bed,bim,fam,frq} and marker_to_rsid.tsv.gz
# (RSID, HG38).  Every step leaves a marker in $WORK_DIR/1000g_30x/state, so re-running resumes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"

TMP="$WORK_DIR/1000g_30x"
OUT="$BUNDLE_DIR/genomic_data/1000genomes/b38_dbsnp156"
STATE="$TMP/state"
mkdir -p "$TMP/vcf" "$TMP/pgen" "$TMP/bed" "$STATE" "$OUT"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
done_or_run() {  # done_or_run NAME command...: run once
  local name="$1"; shift
  if [[ -f "$STATE/$name.done" ]]; then return 0; fi
  "$@"
  touch "$STATE/$name.done"
}

for tool in "$PLINK2" "$PLINK19" curl; do
  command -v "$tool" >/dev/null || { echo "$tool not found: see environment.yml"; exit 1; }
done
[[ -n "$DBSNP_B38_DIR" && -f "$DBSNP_B38_DIR/chr1.fst" ]] || { echo "Set DBSNP_B38_DIR to the b38_dbsnp156 .fst directory"; exit 1; }

# --- samples: the 2,504 unrelated phase 3 people, by superpopulation
if [[ ! -s "$TMP/samples.panel" ]]; then
  curl -fsSL "$KG_PANEL_URL" -o "$TMP/samples.panel"
fi
awk 'NR > 1 {print $1}' "$TMP/samples.panel" > "$TMP/samples_all.txt"
for pop in $SUPERPOPULATIONS; do
  awk -v pop="$pop" 'NR > 1 && $3 == pop {print $1}' "$TMP/samples.panel" > "$TMP/samples_$pop.txt"
  log "$pop: $(wc -l < "$TMP/samples_$pop.txt") people"
done

# --- source: PLINK 2's prebuilt all_hg38 fileset (one download), or the per-chromosome VCFs
case "$KG30X_SOURCE" in
  plink2)
    P2="$TMP/plink2_hg38"
    mkdir -p "$P2"
    if [[ ! -f "$STATE/download_plink2.done" ]]; then
      log "downloading PLINK 2's all_hg38 fileset"
      curl -fsSL -C - "$PLINK2_HG38_PGEN_URL" -o "$P2/all_hg38.pgen.zst"
      curl -fsSL -C - "$PLINK2_HG38_PVAR_URL" -o "$P2/all_hg38.pvar.zst"
      curl -fsSL "$PLINK2_HG38_PSAM_URL" -o "$P2/all_hg38.psam"
      touch "$STATE/download_plink2.done"
    fi
    if [[ ! -f "$STATE/decompress_plink2.done" ]]; then
      "$PLINK2" --zst-decompress "$P2/all_hg38.pgen.zst" > "$P2/all_hg38.pgen"
      rm "$P2/all_hg38.pgen.zst"
      touch "$STATE/decompress_plink2.done"
    fi
    log "all_hg38: $(($(wc -l < "$P2/all_hg38.psam") - 1)) people;" \
      "$(grep -cwFf "$TMP/samples_all.txt" <(cut -f1 "$P2/all_hg38.psam")) of the $(wc -l < "$TMP/samples_all.txt") phase 3 people"
    ;;
  vcf) ;;
  *) echo "KG30X_SOURCE must be plink2 or vcf, not '$KG30X_SOURCE'"; exit 1 ;;
esac

# --- per chromosome: keep the 2,504 people and biallelic SNVs, IDs chr:pos:ref:alt
for chrom in $CHROMOSOMES; do
  if [[ "$KG30X_SOURCE" == plink2 ]]; then
    input=(--pfile "$TMP/plink2_hg38/all_hg38" vzs)
  else
    vcf="$TMP/vcf/${KG30X_FILE_PATTERN//\{CHR\}/$chrom}"
    if [[ ! -f "$STATE/download_chr$chrom.done" ]]; then
      log "chr$chrom: downloading"
      curl -fsSL -C - "$KG30X_URL/$(basename "$vcf")" -o "$vcf"
      curl -fsSL "$KG30X_URL/$(basename "$vcf").tbi" -o "$vcf.tbi"
      touch "$STATE/download_chr$chrom.done"
    fi
    input=(--vcf "$vcf")
  fi
  done_or_run "pgen_chr$chrom" "$PLINK2" "${input[@]}" --keep "$TMP/samples_all.txt" \
    --snps-only just-acgt --max-alleles 2 --chr "$chrom" --output-chr 26 \
    --set-all-var-ids '@:#:$r:$a' --new-id-max-allele-len 10 missing \
    --rm-dup exclude-all --threads "$THREADS" --make-pgen --out "$TMP/pgen/chr$chrom"
  for pop in $SUPERPOPULATIONS; do
    done_or_run "bed_${pop}_chr$chrom" "$PLINK2" --pfile "$TMP/pgen/chr$chrom" --keep "$TMP/samples_$pop.txt" \
      --maf "$MAF_MIN" --threads "$THREADS" --make-bed --out "$TMP/bed/${pop}_chr$chrom"
  done
done

# --- rsIDs from dbSNP 156 (GRCh38), via GeneHackman's chrpos_to_rsid(); writes the ID map, A1 file and marker map
# shellcheck disable=SC2086
done_or_run rsids $R_CMD "$HERE/11_assign_rsids.R" --bed_dir "$TMP/bed" --dbsnp_dir "$DBSNP_B38_DIR" \
  --out_dir "$TMP" --marker_map "$OUT/marker_to_rsid.tsv.gz"

# --- per superpopulation: merge chromosomes, rename to rsIDs, alphabetical A1, allele frequencies
for pop in $SUPERPOPULATIONS; do
  if [[ -f "$STATE/final_$pop.done" ]]; then continue; fi
  : > "$TMP/merge_$pop.txt"
  for chrom in $CHROMOSOMES; do echo "$TMP/bed/${pop}_chr$chrom" >> "$TMP/merge_$pop.txt"; done
  "$PLINK19" --merge-list "$TMP/merge_$pop.txt" --keep-allele-order --make-bed --out "$TMP/${pop}_merged"
  # keep only variants with an rsID, then rename them (two steps: PLINK 1.9 applies --extract before renaming)
  "$PLINK19" --bfile "$TMP/${pop}_merged" --keep-allele-order --extract <(cut -f1 "$TMP/rsid_map.txt") \
    --make-bed --out "$TMP/${pop}_kept"
  "$PLINK19" --bfile "$TMP/${pop}_kept" --keep-allele-order --update-name "$TMP/rsid_map.txt" 2 1 \
    --make-bed --out "$TMP/${pop}_rsid"
  "$PLINK19" --bfile "$TMP/${pop}_rsid" --a1-allele "$TMP/a1_alphabetical.txt" 2 1 --make-bed --out "$OUT/$pop"
  "$PLINK19" --bfile "$OUT/$pop" --freq --out "$OUT/$pop"
  rm -f "$OUT/$pop.log" "$OUT/$pop.nosex" "$TMP/${pop}_merged".* "$TMP/${pop}_kept".* "$TMP/${pop}_rsid".*
  log "$pop: $(wc -l < "$OUT/$pop.bim") variants, $(wc -l < "$OUT/$pop.fam") people"
  touch "$STATE/final_$pop.done"
done
log "GRCh38 1000 Genomes panel written to $OUT"
