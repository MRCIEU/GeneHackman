# Assign dbSNP 156 rsIDs to the GRCh38 1000 Genomes panel variants (called by 10_build_1000g_panel.sh).
#
#   Rscript 11_assign_rsids.R --bed_dir <dir with <POP>_chr<C>.bim> --dbsnp_dir <b38_dbsnp156 dir> \
#     --out_dir <dir> --marker_map <marker_to_rsid.tsv.gz>
#
# Variant IDs in the .bim files are CHR:BP:REF:ALT.  RSIDs come from GeneHackman's chrpos_to_rsid() (matching
# CHR, BP and both alleles, either orientation).  Variants without an rsID, and rsIDs that would name more than
# one variant, are dropped.  Writes:
#   <out_dir>/rsid_map.txt         old ID <TAB> rsID            (PLINK --update-name, columns 2 1)
#   <out_dir>/a1_alphabetical.txt  rsID <TAB> alphabetically smaller allele (PLINK --a1-allele, columns 2 1)
#   <marker_map>                   RSID <TAB> HG38 (GeneHackman SNP key CHR:BP_EA_OA in GRCh38)
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  i <- which(args == name)
  if (length(i) == 0) stop("missing ", name)
  args[i + 1]
}
bed_dir <- get_arg("--bed_dir")
dbsnp_dir <- normalizePath(get_arg("--dbsnp_dir"))
out_dir <- get_arg("--out_dir")
marker_map <- get_arg("--marker_map")

suppressMessages(pkgload::load_all(Sys.getenv("GENEHACKMAN_DIR"), quiet = TRUE, export_all = TRUE))
suppressPackageStartupMessages(library(data.table))

bims <- list.files(bed_dir, pattern = "_chr[0-9]+\\.bim$", full.names = TRUE)
if (length(bims) == 0) stop("No .bim files in ", bed_dir)
variants <- unique(rbindlist(lapply(bims, fread, header = FALSE, select = c(1, 2, 4, 5, 6),
                                    col.names = c("CHR", "ID", "BP", "A1", "A2"),
                                    colClasses = list(character = c(1, 2, 5, 6)))))
message(nrow(variants), " distinct variants across ", length(bims), " population/chromosome files")

# one forked process per chromosome (chrpos_to_rsid's own parallelism starts new R sessions, which can't load a
# GeneHackman that is only load_all()ed)
cores <- as.integer(Sys.getenv("THREADS", "8"))
per_chrom <- parallel::mclapply(split(variants, by = "CHR"), function(d) {
  as.data.table(chrpos_to_rsid(d, "CHR", "BP", "A1", "A2", flip = "allow", build = basename(dbsnp_dir),
                               dbsnp_dir = dirname(dbsnp_dir), parallel_cores = 1, verbose = FALSE))
}, mc.cores = cores, mc.preschedule = FALSE)
failed <- vapply(per_chrom, function(x) inherits(x, "try-error"), logical(1))
if (any(failed)) stop("rsID lookup failed for chromosome(s) ", paste(names(per_chrom)[failed], collapse = ", "), ": ",
                      per_chrom[[which(failed)[1]]])
with_rsid <- rbindlist(per_chrom)
message(sum(!is.na(with_rsid$RSID)), " of ", nrow(with_rsid), " variants have a dbSNP 156 rsID")
with_rsid <- with_rsid[!is.na(RSID)]
shared <- with_rsid[, .N, by = RSID][N > 1, RSID]
with_rsid <- with_rsid[!RSID %in% shared]
message(nrow(with_rsid), " variants with a unique rsID (", length(shared), " rsIDs naming several variants dropped)")

with_rsid[, EA := pmin(toupper(A1), toupper(A2))]
with_rsid[, OA := pmax(toupper(A1), toupper(A2))]
with_rsid[, HG38 := toupper(paste0(CHR, ":", format(BP, scientific = FALSE, trim = TRUE), "_", EA, "_", OA))]

fwrite(with_rsid[, .(ID, RSID)], file.path(out_dir, "rsid_map.txt"), sep = "\t", col.names = FALSE)
fwrite(with_rsid[, .(RSID, EA)], file.path(out_dir, "a1_alphabetical.txt"), sep = "\t", col.names = FALSE)
fwrite(with_rsid[order(as.integer(CHR), BP), .(RSID, HG38)], marker_map, sep = "\t")
message("Wrote ", marker_map)
