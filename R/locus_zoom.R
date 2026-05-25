#' Generate locus zoom plots for significant colocalization results
#'
#' For each coloc result with PP.H4.abf above the threshold, produces a stacked
#' locus zoom plot showing both traits at the shared locus, with LD computed
#' from the local 1000 Genomes reference panel.
#'
#' @param coloc_results_file path to coloc_results.tsv
#' @param gwas_files named character vector of standardised GWAS file paths,
#'   names must match trait labels used in coloc_results (trait1/trait2 columns)
#' @param ancestry ancestry code for the LD reference panel (EUR, EAS, etc.)
#' @param ens_db Ensembl database object or package name for gene annotations
#' @param output_dir directory to write per-locus PNG files
#' @param pp_h4_threshold minimum PP.H4.abf to consider significant (default 0.8)
#' @param window_kb half-width of the region to plot in kb (default 500)
#' @export
locus_zoom_coloc <- function(coloc_results_file,
                             gwas_files,
                             ancestry,
                             ens_db,
                             output_dir,
                             pp_h4_threshold = 0.8,
                             window_kb = 500) {

  message("=== Locus zoom plot generation ===")
  message("Coloc results file: ", coloc_results_file)
  message("GWAS files: ", paste(names(gwas_files), gwas_files, sep = " -> ", collapse = "\n  "))
  message("Ancestry: ", ancestry)
  message("PP.H4 threshold: ", pp_h4_threshold)
  message("Window: ±", window_kb, " kb")

  coloc_res <- vroom::vroom(coloc_results_file, show_col_types = FALSE)
  message("Total coloc results: ", nrow(coloc_res))

  significant <- coloc_res[coloc_res$PP.H4.abf >= pp_h4_threshold, ]
  message("Significant results (PP.H4.abf >= ", pp_h4_threshold, "): ", nrow(significant))

  if (nrow(significant) == 0) {
    message("No coloc results above threshold — nothing to plot")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    return(invisible(NULL))
  }

  traits_in_results <- unique(c(significant$trait1, significant$trait2))
  missing_traits <- setdiff(traits_in_results, names(gwas_files))
  if (length(missing_traits) > 0) {
    message("WARNING: Traits in coloc results but missing from gwas_files: ",
            paste(missing_traits, collapse = ", "))
    message("  Available trait names: ", paste(names(gwas_files), collapse = ", "))
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  gwas_cache <- list()
  load_gwas <- function(trait) {
    if (is.null(gwas_cache[[trait]])) {
      gwas_cache[[trait]] <<- get_file_or_dataframe(gwas_files[[trait]])
    }
    gwas_cache[[trait]]
  }

  for (i in seq_len(nrow(significant))) {
    row <- significant[i, ]
    locus_label <- paste0(row$trait1, " vs ", row$trait2, " at ", row$locus1, "/", row$locus2)
    message("Processing locus ", i, "/", nrow(significant), ": ", locus_label)

    tryCatch({
      chr <- as.integer(strsplit(row$locus1, "_")[[1]][1])
      bp1 <- as.numeric(strsplit(row$locus1, "_")[[1]][2])
      bp2 <- as.numeric(strsplit(row$locus2, "_")[[1]][2])
      center_bp <- round(mean(c(bp1, bp2)))

      window_bp <- window_kb * 1000
      xrange <- c(center_bp - window_bp, center_bp + window_bp)

      index_snp <- find_index_snp_from_ld_panel(chr, center_bp, window_bp, ancestry)
      message("  Region: chr", chr, ":", xrange[1], "-", xrange[2],
              ", index SNP: ", if (is.null(index_snp)) "NONE" else index_snp)

      plot_list <- list()
      for (trait in c(row$trait1, row$trait2)) {
        if (!trait %in% names(gwas_files)) {
          message("  WARNING: trait '", trait, "' not found in gwas_files. ",
                  "Available: ", paste(names(gwas_files), collapse = ", "))
          next
        }
        gwas <- load_gwas(trait)
        message("  Building plot for '", trait, "' (", nrow(gwas), " rows in GWAS)")
        p <- build_single_locus_plot(gwas, chr, xrange, ens_db, ancestry,
                                     index_snp, trait)
        if (!is.null(p)) {
          plot_list[[trait]] <- p
        } else {
          message("  WARNING: build_single_locus_plot returned NULL for '", trait, "'")
        }
      }

      if (length(plot_list) < 2) {
        message("  Skipping — only ", length(plot_list), " of 2 trait plots could be built")
        next
      }

      hit_info <- ""
      if (!is.na(row$hit1) && !is.na(row$hit2)) {
        hit_info <- paste0("_", row$hit1, "_vs_", row$hit2)
      }

      safe_name <- gsub("[^A-Za-z0-9._-]+", "_", paste0(
        row$trait1, "_vs_", row$trait2, "_chr", chr, "_", center_bp, hit_info
      ))
      out_file <- file.path(output_dir, paste0(safe_name, ".png"))

      stacked <- patchwork::wrap_plots(plot_list, ncol = 1)

      ggplot2::ggsave(out_file, plot = stacked,
                      width = 10, height = 5 * length(plot_list),
                      dpi = 300, limitsize = FALSE)
      message("  Saved: ", out_file)
    }, error = function(e) {
      message("  ERROR processing locus ", locus_label, ": ", conditionMessage(e))
    })
  }

  return(invisible(NULL))
}


#' Find the top SNP near a locus center from the LD reference panel
#' @keywords internal
find_index_snp_from_ld_panel <- function(chr, center_bp, window_bp, ancestry) {
  bfile <- file.path(thousand_genomes_dir, ancestry)
  bim_file <- paste0(bfile, ".bim")
  if (!file.exists(bim_file)) return(NULL)

  bim <- data.table::fread(bim_file, header = FALSE,
                           select = c(1, 2, 4),
                           col.names = c("CHR", "SNP", "BP"))
  bim <- bim[CHR == chr & BP >= (center_bp - window_bp) & BP <= (center_bp + window_bp)]
  if (nrow(bim) == 0) return(NULL)

  bim[which.min(abs(bim$BP - center_bp)), ]$SNP
}


#' Build a single ggplot2 locus zoom panel for one trait
#' @keywords internal
build_single_locus_plot <- function(gwas, chr, xrange, ens_db, ancestry,
                                    index_snp, trait_label) {
  gwas_chr <- as.numeric(gwas$CHR)
  gwas_bp <- as.numeric(gwas$BP)
  keep <- !is.na(gwas_chr) & gwas_chr == chr &
          !is.na(gwas_bp) & gwas_bp >= xrange[1] & gwas_bp <= xrange[2]
  gwas_sub <- gwas[keep, ]

  message("    Subsetting chr", chr, " ", xrange[1], "-", xrange[2],
          ": ", nrow(gwas_sub), " SNPs in window")
  if (nrow(gwas_sub) < 2) return(NULL)

  plot_df <- data.frame(
    chrom = as.numeric(gwas_sub$CHR),
    pos = as.numeric(gwas_sub$BP),
    p = as.numeric(gwas_sub$P),
    stringsAsFactors = FALSE
  )
  if ("RSID" %in% colnames(gwas_sub)) {
    plot_df$rsid <- gwas_sub$RSID
  } else if ("SNP" %in% colnames(gwas_sub)) {
    plot_df$rsid <- gwas_sub$SNP
  }
  if ("BETA" %in% colnames(gwas_sub)) plot_df$beta <- as.numeric(gwas_sub$BETA)

  plot_df <- plot_df[!is.na(plot_df$p) & plot_df$p > 0, ]
  message("    After filtering NA/zero p-values: ", nrow(plot_df), " SNPs")
  if (nrow(plot_df) < 2) return(NULL)

  ld_r2 <- compute_ld_r2_for_locus(plot_df$rsid, chr, ancestry, index_snp)
  if (!is.null(ld_r2)) {
    plot_df$r2 <- ld_r2[match(plot_df$rsid, names(ld_r2))]
    message("    LD r2 computed for ", sum(!is.na(plot_df$r2)), "/", nrow(plot_df), " SNPs")
  } else {
    message("    LD r2 computation returned NULL — plotting without LD colouring")
  }

  ld_arg <- if ("r2" %in% colnames(plot_df)) "r2" else NULL

  loc <- locuszoomr::locus(
    data = plot_df,
    seqname = as.character(chr),
    xrange = xrange,
    ens_db = ens_db,
    chrom = "chrom",
    pos = "pos",
    p = "p",
    labs = "rsid",
    index_snp = index_snp,
    LD = ld_arg
  )

  p <- locuszoomr::locus_ggplot(loc, labels = "index",
                                filter_gene_biotype = "protein_coding")
  p <- p + ggplot2::ggtitle(trait_label)
  p
}


#' Compute r2 relative to an index SNP using the local LD reference panel
#' @keywords internal
compute_ld_r2_for_locus <- function(rsids, chr, ancestry, index_snp) {
  if (is.null(index_snp) || !index_snp %in% rsids) return(NULL)

  rsids_unique <- unique(rsids[!is.na(rsids) & nchar(rsids) > 0])
  if (length(rsids_unique) < 2) return(NULL)

  ld_result <- compute_ld_matrix(rsids_unique, chr, ancestry)
  if (is.null(ld_result)) return(NULL)

  R <- ld_result$matrix
  snps <- ld_result$snps

  if (!index_snp %in% snps) return(NULL)

  idx <- which(snps == index_snp)
  r2_vec <- R[idx, ]^2
  names(r2_vec) <- snps
  r2_vec
}
