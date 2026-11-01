#!/usr/bin/env bash
# Assemble reference bundle 1.2.0 (GRCh38) and the GRCh37 extras from the pieces built here and the 1.1.0 bundle:
#
#   bash reference_data/grch38/50_assemble_bundle.sh
#
#   $BUNDLE_DIR (gs://genehackman/1.2.0/):
#     genomic_data/1000genomes/b38_dbsnp156/   from 10_build_1000g_panel.sh (already written there)
#     genomic_data/dbsnp/b38_dbsnp156/         from $DBSNP_B38_DIR
#     genomic_data/pan_ukb_ld/                 from $LD_PACK_DIR
#     genomic_data/liftover/                   the 1.1.0 chains + hg18ToHg38 (UCSC)
#     genomic_data/gene_name_map.tsv and anything else in the 1.1.0 genomic_data/ that is not build-specific
#     LDSCORE/pan_ukb/                         from 20_import_pan_ukb_ldscores.py (already written there)
#   $BUNDLE_DIR-grch37 (gs://genehackman/1.2.0-grch37/):
#     genomic_data/1000genomes/b37_dbsnp156/, genomic_data/dbsnp/b37_dbsnp156/, LDSCORE/b37_dbsnp156/ from 1.1.0
#
# Files are hard-linked when source and bundle share a filesystem (no extra disk), copied otherwise.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=config.sh
source "$HERE/config.sh"
[[ -d "$OLD_PIPELINE_DATA_DIR/genomic_data" ]] || { echo "Set OLD_PIPELINE_DATA_DIR to the 1.1.0 reference data"; exit 1; }
EXTRAS_DIR="${EXTRAS_DIR:-${BUNDLE_DIR%/}-grch37}"

link_tree() {  # link_tree SOURCE DEST: hard-link (or copy) a directory tree, keeping what is already there
  local source="$1" dest="$2"
  mkdir -p "$dest"
  cp -alnT "$source" "$dest" 2>/dev/null || cp -anT "$source" "$dest"
  echo "  $dest <- $source"
}

echo "Main bundle: $BUNDLE_DIR"
old="$OLD_PIPELINE_DATA_DIR/genomic_data"
for item in "$old"/*; do
  name="$(basename "$item")"
  case "$name" in
    1000genomes|dbsnp|pipeline|pan_ukb_ld) ;;  # build-specific, container cache, or added below
    *) if [[ -d "$item" ]]; then link_tree "$item" "$BUNDLE_DIR/genomic_data/$name"
       else mkdir -p "$BUNDLE_DIR/genomic_data"; ln -f "$item" "$BUNDLE_DIR/genomic_data/$name" 2>/dev/null \
            || cp -f "$item" "$BUNDLE_DIR/genomic_data/$name"; echo "  $BUNDLE_DIR/genomic_data/$name <- $item"; fi ;;
  esac
done
chain="$BUNDLE_DIR/genomic_data/liftover/hg18ToHg38.over.chain.gz"
if [[ ! -s "$chain" ]]; then
  curl -fsSL https://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz -o "$chain"
  echo "  $chain <- UCSC"
fi
[[ -d "$DBSNP_B38_DIR" ]] && link_tree "$DBSNP_B38_DIR" "$BUNDLE_DIR/genomic_data/dbsnp/b38_dbsnp156"
[[ -f "$LD_PACK_DIR/MANIFEST.json" ]] && link_tree "$LD_PACK_DIR" "$BUNDLE_DIR/genomic_data/pan_ukb_ld"
for needed in genomic_data/1000genomes/b38_dbsnp156/EUR.bed LDSCORE/pan_ukb/EUR/ld.l2.ldscore.gz; do
  [[ -f "$BUNDLE_DIR/$needed" ]] || echo "  WARNING: $BUNDLE_DIR/$needed missing (run the earlier steps)"
done

echo "GRCh37 extras: $EXTRAS_DIR"
link_tree "$old/1000genomes/b37_dbsnp156" "$EXTRAS_DIR/genomic_data/1000genomes/b37_dbsnp156"
link_tree "$old/dbsnp/b37_dbsnp156" "$EXTRAS_DIR/genomic_data/dbsnp/b37_dbsnp156"
link_tree "$OLD_PIPELINE_DATA_DIR/LDSCORE/b37_dbsnp156" "$EXTRAS_DIR/LDSCORE/b37_dbsnp156"

du -sh "$BUNDLE_DIR" "$EXTRAS_DIR" 2>/dev/null || true
echo "Upload (from a machine with gsutil and write access):"
echo "  gsutil -m rsync -r $BUNDLE_DIR/ gs://genehackman/1.2.0/"
echo "  gsutil -m rsync -r $EXTRAS_DIR/ gs://genehackman/1.2.0-grch37/"
