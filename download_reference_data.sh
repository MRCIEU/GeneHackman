#!/usr/bin/env bash
# Download GeneHackman's reference data, and optionally its QTL data, for only the ancestries you need.
#
#   ./download_reference_data.sh --dir /path/to/pipeline_data                      # EUR, no QTL data
#   ./download_reference_data.sh --dir /path/to/pipeline_data --ancestries EUR,EAS --qtl eqtlgen
#   ./download_reference_data.sh --dry-run ...                                      # show what and how much
#
# Run it again with more ancestries or datasets to add them: files already there are skipped, and an
# interrupted download resumes.  Run ./download_reference_data.sh --help for every option.
set -euo pipefail

DATA_VERSION="1.2.0"   # the reference data release this version of GeneHackman reads
GRCH37_VERSION="1.1.0" # the release whose GRCh37 files --grch37 downloads
DATA_URL="${GENEHACKMAN_DATA_URL:-https://storage.googleapis.com/genehackman}"
QTL_URL="${GENEHACKMAN_QTL_URL:-https://storage.googleapis.com/genehackman-qtl}"
ALL_ANCESTRIES="EUR,EAS,AFR,AMR,SAS,MID"
ALL_QTL="eqtlgen,brainmeta,ukb_pqtl"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
  cat <<EOF
Usage: ./download_reference_data.sh [options]

Downloads GeneHackman's reference data (release $DATA_VERSION, GRCh38) into one directory, which becomes
PIPELINE_DATA_DIR in .env.

Options:
  --dir DIR            where to put the reference data (default: PIPELINE_DATA_DIR from .env)
  --ancestries LIST    comma-separated, from EUR EAS AFR AMR SAS MID, or "all" (default: EUR).
                       Each ancestry brings its 1000 Genomes panel, Pan-UKB LD pack and LD scores.
                       SAS is called CSA in Pan-UKB; MID has no 1000 Genomes panel (ld_reference: pan_ukb only).
  --qtl LIST           QTL datasets for the qtl_mr pipeline: none, all, or a comma-separated list
                       from eqtlgen brainmeta ukb_pqtl (default: none)
  --qtl-dir DIR        where to put the QTL data (default: QTL_DATA_DIR from .env, else DIR/qtl_data)
  --no-dbsnp           leave out dbSNP (11.7 GB). It is only used with populate_rsid: true; the default RSID
                       fill uses the 1000 Genomes marker map, which is always downloaded.
  --grch37             also download the GRCh37 files needed for output.build: GRCh37 (from release 1.1.0):
                       the 1000 Genomes phase 3 panel and earlier LD scores for the chosen ancestries, and
                       dbSNP GRCh37 unless --no-dbsnp.  QTL data are GRCh38 only.
  --jobs N             parallel downloads (default: 8)
  --verify             check the md5 of files already downloaded too (slow for large files)
  --dry-run            only show what would be downloaded, and its size
  --version V          reference data release (default: $DATA_VERSION)
  -h, --help           show this help
EOF
}

die() { echo "Error: $*" >&2; exit 1; }

env_value() {  # env_value NAME: a value from the repository's .env, if set
  [[ -f "$HERE/.env" ]] || return 0
  sed -n "s/^[[:space:]]*$1[[:space:]]*=[[:space:]]*//p" "$HERE/.env" | tail -1 | sed 's/^["'\'']//; s/["'\'']$//; s/[[:space:]]*$//'
}

DIR="" ANCESTRIES="EUR" QTL="none" QTL_DIR="" DBSNP=1 GRCH37=0 JOBS=8 VERIFY=0 DRY_RUN=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dir) DIR="${2:?--dir needs a value}"; shift 2 ;;
    --ancestries) ANCESTRIES="${2:?--ancestries needs a value}"; shift 2 ;;
    --qtl) QTL="${2:?--qtl needs a value}"; shift 2 ;;
    --qtl-dir) QTL_DIR="${2:?--qtl-dir needs a value}"; shift 2 ;;
    --no-dbsnp) DBSNP=0; shift ;;
    --grch37) GRCH37=1; shift ;;
    --jobs) JOBS="${2:?--jobs needs a value}"; shift 2 ;;
    --verify) VERIFY=1; shift ;;
    --dry-run) DRY_RUN=1; shift ;;
    --version) DATA_VERSION="${2:?--version needs a value}"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) usage >&2; die "unknown option $1" ;;
  esac
done

command -v curl >/dev/null || die "curl is needed"
[[ "$JOBS" =~ ^[1-9][0-9]*$ ]] || die "--jobs must be a positive number"

# --- what to download
DIR="${DIR:-$(env_value PIPELINE_DATA_DIR)}"
[[ -n "$DIR" ]] || die "give --dir (or set PIPELINE_DATA_DIR in .env)"
DIR="${DIR%/}"

ANCESTRIES="$(echo "$ANCESTRIES" | tr '[:lower:]' '[:upper:]' | tr ' ' ',')"
[[ "$ANCESTRIES" == "ALL" ]] && ANCESTRIES="$ALL_ANCESTRIES"
wanted_ancestries=","
for a in ${ANCESTRIES//,/ }; do
  case "$a" in
    EUR|EAS|AFR|AMR|MID) wanted_ancestries+="$a," ;;
    SAS|CSA) wanted_ancestries+="SAS,CSA," ;;  # 1000 Genomes SAS = Pan-UKB CSA
    *) die "unknown ancestry $a; choose from EUR EAS AFR AMR SAS MID" ;;
  esac
done

QTL="$(echo "$QTL" | tr '[:upper:]' '[:lower:]' | tr ' ' ',')"
[[ "$QTL" == "all" ]] && QTL="$ALL_QTL"
[[ "$QTL" == "none" ]] && QTL=""
wanted_qtl=","
for q in ${QTL//,/ }; do
  case "$q" in
    eqtlgen|brainmeta|ukb_pqtl) wanted_qtl+="$q," ;;
    *) die "unknown QTL dataset $q; choose from eqtlgen brainmeta ukb_pqtl (or all, none)" ;;
  esac
done
if [[ "$wanted_qtl" != "," ]]; then
  QTL_DIR="${QTL_DIR:-$(env_value QTL_DATA_DIR)}"
  QTL_DIR="${QTL_DIR:-$DIR/qtl_data}"
  QTL_DIR="${QTL_DIR%/}"
fi

wanted_groups=",core,1000genomes,pan_ukb_ld,ldscores,"
[[ $DBSNP == 1 ]] && wanted_groups+="dbsnp,"
if [[ $GRCH37 == 1 ]]; then
  wanted_groups+="grch37_1000genomes,grch37_ldscores,"
  [[ $DBSNP == 1 ]] && wanted_groups+="grch37_dbsnp,"
fi

# --- helpers
md5_of() {
  if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi
}
size_of() { stat -L -c %s "$1" 2>/dev/null || stat -L -f %z "$1"; }
list_sizes() {  # list_sizes ROOT: "relative path<TAB>bytes" for every file under ROOT
  [[ -d "$1" ]] || return 0
  (cd "$1" && find -L . -type f ! -name '*.part' -exec stat -L -c '%n	%s' {} + 2>/dev/null \
    || find -L . -type f ! -name '*.part' -exec stat -L -f '%N	%z' {} +) | sed 's|^\./||'
}
gb() { awk -v b="$1" 'BEGIN { printf "%.2f GB", b / 1e9 }'; }
export -f md5_of size_of

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

# list_objects BASE_URL PREFIX: "object name<TAB>bytes<TAB>md5 (hex)" for every file under PREFIX in a bucket.
# BASE_URL is https://storage.googleapis.com/<bucket> (listed with the public JSON API), or file:///<dir> for
# testing against local files (no md5, so only sizes are checked).
list_objects() {
  local base="$1" prefix="$2" token="" page api
  if [[ "$base" == file://* ]]; then
    local root="${base#file://}"
    [[ -d "$root/$prefix" ]] || return 0
    list_sizes "$root/$prefix" | awk -F'\t' -v OFS='\t' -v p="$prefix" '$1 != "" { print p $1, $2, "" }'
    return 0
  fi
  [[ "$base" =~ ^https://storage\.googleapis\.com/([^/]+)$ ]] || die "can't list $base (expected https://storage.googleapis.com/<bucket>)"
  api="https://storage.googleapis.com/storage/v1/b/${BASH_REMATCH[1]}/o?prefix=$prefix&maxResults=1000&prettyPrint=false&fields=items(name,size,md5Hash),nextPageToken"
  while :; do
    page="$(curl -fsSL --retry 3 "$api${token:+&pageToken=$token}")" \
      || die "could not list $base/$prefix. Check your internet connection."
    # one JSON object per line; names hold no quotes or braces
    printf '%s\n' "$page" | tr '{' '\n' | awk -v OFS='\t' '
      function field(name,   v) {
        if (!match($0, "\"" name "\":\"[^\"]*\"")) return ""
        v = substr($0, RSTART, RLENGTH); sub("^\"" name "\":\"", "", v); sub("\"$", "", v); return v
      }
      function b64hex(s,   alphabet, i, c, bits, nbits, byte, out) {  # GCS md5Hash (base64) -> hex
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
        for (i = 1; i <= length(s); i++) {
          c = substr(s, i, 1); if (c == "=") break
          bits = bits * 64 + index(alphabet, c) - 1; nbits += 6
          if (nbits >= 8) { nbits -= 8; byte = int(bits / 2 ^ nbits); bits -= byte * 2 ^ nbits; out = out sprintf("%02x", byte) }
        }
        return out
      }
      { name = field("name"); if (name != "" && name !~ /\/$/) print name, field("size"), b64hex(field("md5Hash")) }'
    token="$(printf '%s' "$page" | sed -n 's/.*"nextPageToken":"\([^"]*\)".*/\1/p')"
    [[ -n "$token" ]] || break
  done
}

# classify: reads list_objects output on stdin; writes "group  ancestry  path  bytes  md5  object" for the files
# GeneHackman uses (path relative to the data directory).  Groups: core (always), 1000genomes, pan_ukb_ld,
# ldscores (per ancestry), dbsnp, grch37_1000genomes, grch37_ldscores, grch37_dbsnp (from release 1.1.0), and
# the QTL datasets.
classify() {
  awk -F'\t' -v OFS='\t' -v version="$DATA_VERSION" -v grch37="$GRCH37_VERSION" '
    function out(group, ancestry) { print group, ancestry, rel, $2, ($3 == "" ? "-" : $3), $1 }  # "-": no md5
    {
      n = split($1, parts, "/"); file = parts[n]; rel = substr($1, length(parts[1]) + 2)
      panel_pop = (file ~ /^(EUR|EAS|AFR|AMR|SAS)\.(bed|bim|fam|frq)$/) ? substr(file, 1, 3) : ""
      if (parts[1] == version) {
        if (rel ~ /^genomic_data\/1000genomes\/b38_dbsnp156\//) {
          if (file == "marker_to_rsid.tsv.gz") out("core", "-"); else if (panel_pop != "") out("1000genomes", panel_pop)
        } else if (rel ~ /^genomic_data\/dbsnp\/b38_dbsnp156\//) out("dbsnp", "-")
        else if (rel ~ /^genomic_data\/pan_ukb_ld\//) out("pan_ukb_ld", n == 4 ? "-" : parts[4])
        else if (rel ~ /^LDSCORE\/pan_ukb\//) out("ldscores", n == 4 ? "-" : parts[4])
        else if (rel ~ /^genomic_data\/liftover\// || rel == "genomic_data/gene_name_map.tsv") out("core", "-")
        else if (parts[2] ~ /^(eqtlgen|brainmeta|ukb_pqtl)$/ && parts[3] == "GRCh38") out(parts[2], "-")
      } else if (parts[1] == grch37) {
        if (rel ~ /^genomic_data\/1000genomes\/b37_dbsnp156\//) {
          if (file == "marker_to_rsid.tsv.gz") out("grch37_1000genomes", "-")
          else if (panel_pop != "") out("grch37_1000genomes", panel_pop)
        } else if (rel ~ /^genomic_data\/dbsnp\/b37_dbsnp156\//) out("grch37_dbsnp", "-")
        else if (rel ~ /^LDSCORE\/b37_dbsnp156\//) {
          if (n == 4) out("grch37_ldscores", "-"); else if (parts[4] ~ /^(EUR|EAS|AFR|AMR|SAS)$/) out("grch37_ldscores", parts[4])
        }
      }
    }'
}

# select_files ROOT BASE_URL GROUPS: rows of "url  destination  bytes  md5  group  ancestry  status" from the
# classified listing on stdin
select_files() {
  local root="$1" base="$2" groups="$3"
  cat > "$WORK/listing.tsv"
  { printf '\t\n'; list_sizes "$root"; } > "$WORK/have.tsv"  # never empty, so FNR == NR only matches this file
  awk -F'\t' -v OFS='\t' -v root="$root" -v base="$base" -v groups="$groups" -v anc="$wanted_ancestries" '
    FNR == NR { have[$1] = $2; next }
    index(groups, "," $1 ",") && ($2 == "-" || index(anc, "," $2 ",")) {
      status = ($3 in have && have[$3] == $4) ? "have" : "get"
      print base "/" $6, root "/" $3, $4, $5, $1, $2, status
    }' "$WORK/have.tsv" "$WORK/listing.tsv"
}

{
  list_objects "$DATA_URL" "$DATA_VERSION/"
  if [[ $GRCH37 == 1 ]]; then list_objects "$DATA_URL" "$GRCH37_VERSION/"; fi
} | classify > "$WORK/reference.tsv"
grep -q $'^core\t' "$WORK/reference.tsv" || die "no reference data release $DATA_VERSION at $DATA_URL"
select_files "$DIR" "$DATA_URL" "$wanted_groups" < "$WORK/reference.tsv" > "$WORK/files.tsv"
for q in ${QTL//,/ }; do
  list_objects "$QTL_URL" "$DATA_VERSION/$q/GRCh38/" | classify > "$WORK/qtl.tsv"
  [[ -s "$WORK/qtl.tsv" ]] || die "no $q QTL data in release $DATA_VERSION at $QTL_URL"
  select_files "$QTL_DIR" "$QTL_URL" ",$q," < "$WORK/qtl.tsv" >> "$WORK/files.tsv"
done

if [[ $VERIFY == 1 ]]; then  # files that are there but fail their md5 are fetched again
  echo "Checking the md5 of the files already downloaded..."
  while IFS=$'\t' read -r url dest bytes md5 group ancestry status; do
    if [[ "$status" == have && "$md5" != "-" && "$(md5_of "$dest")" != "$md5" ]]; then
      echo "  md5 mismatch, will download again: $dest" >&2; rm -f "$dest"; status=get
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$url" "$dest" "$bytes" "$md5" "$group" "$ancestry" "$status"
  done < "$WORK/files.tsv" > "$WORK/verified.tsv"
  mv "$WORK/verified.tsv" "$WORK/files.tsv"
fi

# --- summary
echo "GeneHackman reference data $DATA_VERSION"
echo "  reference data: $DIR"
[[ "$wanted_qtl" != "," ]] && echo "  QTL data:       $QTL_DIR"
echo "  ancestries: ${ANCESTRIES}   QTL: ${QTL:-none}   dbSNP: $([[ $DBSNP == 1 ]] && echo yes || echo no)   GRCh37 output: $([[ $GRCH37 == 1 ]] && echo yes || echo no)"
echo
awk -F'\t' -v OFS='\t' '
  { key = $5 " " (($6 == "-") ? "(all)" : $6); total[key] += $3; files[key]++; if ($7 == "get") get[key] += $3 }
  END { for (k in total) print k, total[k], files[k], get[k] + 0 }' "$WORK/files.tsv" | sort \
  | awk -F'\t' '
    BEGIN { printf "  %-26s %10s %8s %16s\n", "", "size", "files", "to download" }
    { printf "  %-26s %7.2f GB %8d %13.2f GB\n", $1, $2 / 1e9, $3, $4 / 1e9; t += $2; n += $3; g += $4 }
    END { printf "  %-26s %7.2f GB %8d %13.2f GB\n", "total", t / 1e9, n, g / 1e9 }'
echo

need=$(awk -F'\t' '$7 == "get" { s += $3 } END { printf "%.0f", s }' "$WORK/files.tsv")
nget=$(awk -F'\t' '$7 == "get"' "$WORK/files.tsv" | wc -l | tr -d ' ')
if [[ $DRY_RUN == 1 ]]; then echo "Dry run: nothing downloaded."; exit 0; fi

if [[ "$nget" -gt 0 ]]; then
  for root in "$DIR" ${QTL_DIR:+"$QTL_DIR"}; do mkdir -p "$root"; done
  for root in "$DIR" ${QTL_DIR:+"$QTL_DIR"}; do  # enough free space where each part goes
    root_need=$(awk -F'\t' -v r="$root/" '$7 == "get" && index($2, r) == 1 { s += $3 } END { printf "%.0f", s }' "$WORK/files.tsv")
    free=$(df -Pk "$root" | awk 'NR == 2 { printf "%.0f", $4 * 1024 }')
    if awk -v n="$root_need" -v f="$free" 'BEGIN { exit !(n > f) }'; then
      die "needs $(gb "$root_need") in $root but only $(gb "$free") is free there"
    fi
  done

  fetch_one() {  # fetch_one URL DEST BYTES MD5: download to DEST.part, check size and md5, then move into place
    local url="$1" dest="$2" bytes="$3" md5="$4" attempt
    mkdir -p "$(dirname "$dest")"
    for attempt in 1 2 3; do
      if [[ ! -f "$dest.part" || "$(size_of "$dest.part")" -lt "$bytes" ]]; then
        curl -fsSL --retry 5 --retry-delay 5 -C - -o "$dest.part" "$url" || { sleep 5; continue; }
      fi
      if [[ "$(size_of "$dest.part")" == "$bytes" && ( "$md5" == "-" || "$(md5_of "$dest.part")" == "$md5" ) ]]; then
        mv -f "$dest.part" "$dest"
        if [[ "$bytes" -ge 100000000 ]]; then echo "  done $dest"; fi
        return 0
      fi
      rm -f "$dest.part"  # wrong size or checksum: start this file again
    done
    echo "  FAILED $url" >&2
    return 1
  }
  export -f fetch_one
  echo "Downloading $nget files ($(gb "$need")) with $JOBS parallel downloads..."
  if ! awk -F'\t' -v OFS='\t' '$7 == "get" { print $1, $2, $3, $4 }' "$WORK/files.tsv" \
      | tr '\t' '\n' | xargs -n 4 -P "$JOBS" bash -c 'fetch_one "$@"' _; then
    die "some files failed to download (listed above). Run the same command again to retry; finished files are kept."
  fi
fi

# --- check everything selected is now in place
missing=$(awk -F'\t' '{ print $2 "\t" $3 }' "$WORK/files.tsv" | while IFS=$'\t' read -r dest bytes; do
  [[ -f "$dest" && "$(size_of "$dest")" == "$bytes" ]] || echo "$dest"; done | wc -l | tr -d ' ')
[[ "$missing" == 0 ]] || die "$missing files are missing or the wrong size after downloading; run the command again"

abs() { (cd "$1" && pwd); }
echo "All files are in place."
echo
echo "Set these in .env:"
echo "  PIPELINE_DATA_DIR=$(abs "$DIR")"
[[ "$wanted_qtl" != "," ]] && echo "  QTL_DATA_DIR=$(abs "$QTL_DIR")"
pan_ukb=$(echo "$wanted_ancestries" | tr ',' '\n' | grep -v '^SAS$' | grep -v '^$' | tr '\n' ' ')
kg=$(echo "$wanted_ancestries" | tr ',' '\n' | grep -v -E '^(CSA|MID)$' | grep -v '^$' | tr '\n' ' ')
echo
echo "This data supports:"
echo "  ld_reference: 1000genomes for ancestries: ${kg:-none}"
echo "  ld_reference: pan_ukb for ancestries:     ${pan_ukb}(SAS is accepted for CSA)"
echo "  populate_rsid: true:                      $([[ $DBSNP == 1 ]] && echo yes || echo "no (rerun without --no-dbsnp)")"
echo "  output.build: GRCh37:                     $([[ $GRCH37 == 1 ]] && echo yes || echo "no (rerun with --grch37)")"
