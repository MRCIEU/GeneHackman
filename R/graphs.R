#' @import ggplot2
#' @import shiny
#' @export
forest_plot <- function(table, title, output_file, y_column=NA) {
  if (!all(c("BETA", "SE") %in% names(table))) {
    stop("data frame needs to have BETA and SE columns")
  }

  first_column_name <- if(shiny::isTruthy(y_column)) y_column else colnames(table)[1]
  table$BETA <- as.numeric(table$BETA)
  table$SE <- as.numeric(table$SE)

  table$LL <- table$BETA - (1.96 * table$SE)
  table$UL <- table$BETA + (1.96 * table$SE)

  ggplot2::ggplot(table,ggplot2::aes(y = .data[[first_column_name]],
                                     x = BETA,
                                     xmin = LL,
                                     xmax = UL,
                                     color = .data[[first_column_name]])) +
    ggplot2::geom_pointrange(cex = 1) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(output_file, limitsize = F)
}

#' @import ggplot2
#' @export
grouped_forest_plot <- function(table, title, group_column, output_file, p_value_column = NA, q_stat_column = NA) {
  if (!("BETA" %in% names(table)) || !("SE" %in% names(table))) {
    stop("data frame needs to have BETA and SE named columns")
  }

  first_column_name <- colnames(table)[1]
  table$BETA <- as.numeric(table$BETA)
  table$SE <- as.numeric(table$SE)

  table$LL <- table$BETA - (1.96 * table$SE)
  table$UL <- table$BETA + (1.96 * table$SE)

  plot_thing <- ggplot2::ggplot(table,
         ggplot2::aes(y = if(!is.na(q_stat_column)) paste0(.data[[first_column_name]], "\n Q-stat=", .data[[q_stat_column]]) else .data[[first_column_name]],
             x = BETA,
             xmin = LL,
             xmax = UL,
             color=.data[[group_column]],
             fill=.data[[group_column]])) +
    ggplot2::ylab(first_column_name) +
    ggplot2::scale_colour_brewer(type="qual") +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::geom_pointrange(cex = 1, fatten = 2, position=ggplot2::position_dodge(width = 0.5)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  if (!is.na(p_value_column)) {
    plot_thing + ggplot2::geom_text(ggplot2::aes(label = paste("P=", signif(.data[[p_value_column]]), digits=3), group = .data[[group_column]]), position = ggplot2::position_dodge(width = 0.5))
  } else if (!is.na(q_stat_column)) {
    plot_thing + ggplot2::geom_text(ggplot2::aes(x = 1.5, label = .data[[q_stat_column]]))
  }

  forest_plot_height <- max(nrow(table)*50, 2000)
  ggplot2::ggsave(output_file, width = 2000, units = "px", height = forest_plot_height)
}

#' @import ggplot2
grouped_bar_chart <- function(data, title, x_column, y_column, group_column, output_file) {
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]], fill = .data[[group_column]])) +
    ggplot2::ggtitle(title) +
    ggplot2::theme_bw() +
    ggplot2::geom_bar(position = "dodge", stat = "identity") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::scale_colour_brewer(type="qual") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) #this centres the title

  ggplot2::ggsave(output_file)
}

#' manhattan_and_qq: produce manhattan and qq plot from a GWAS file
#'
#' @param gwas_filename: a file of a gwas that includes CHR, CP, P, and SNP
#' @param name: name of plots to be saved (and named as a header in graph)
#' @param save_dir: defaults to 'scratch/results'
#' @return 2 plots: one manhattan plot and one QQ plot (with lambda included)
#' @import grDevices
#' @import qqman
#' @import graphics
#' @export
manhattan_and_qq <- function(gwas_filename, manhattan_filename, qq_filename, include_qq = T) {
  manhattan_columns <- c("SNP", "CHR", "BP", "P")
  gwas <- get_file_or_dataframe(gwas_filename, columns = manhattan_columns)
  gwas <- gwas[complete.cases(gwas), ]

  grDevices::png(manhattan_filename, width = 1500, height = 500)
  qqman::manhattan(gwas, main = "Manhattan plot of GWAS")
  grDevices::dev.off()

  if (include_qq) {
    grDevices::png(qq_filename, width = 500, height = 500)

    qqman::qq(gwas$P, main = "Q-Q plot of GWAS p-values")
    lambda <- calculate_lambda_statistic(gwas)
    graphics::text(0.5, 4, paste("lambda", "=", signif(lambda, digits = 3)))

    grDevices::dev.off()
  }
}


#' miami_plot: produce miami plot of GWAS data from two gwases
#'
#' @param gwas_dataframe: a dataframe that includes CHR, CP, P, and SNP
#' @param name: name of plots to be saved (and named as a header in graph)
#' @param save_dir: defaults to 'scratch/results'
#' @import grDevices
#' @import qqman
#' @import graphics
#' @export
miami_plot <- function(first_gwas_filename,
                       second_gwas_filename,
                       miami_plot_file,
                       title = "Comparing GWASes",
                       chr = NA,
                       bp = NA,
                       range = NA) {
  show_specific_region <- !is.na(chr) & !is.na(bp) & !is.na(range)

  manhattan_columns <- c("SNP", "CHR", "BP", "P")
  first_gwas <- get_file_or_dataframe(first_gwas_filename, manhattan_columns)
  second_gwas <- get_file_or_dataframe(second_gwas_filename, manhattan_columns)
  first_gwas <- first_gwas[complete.cases(first_gwas), ]
  second_gwas <- second_gwas[complete.cases(second_gwas), ]

  png_width <- 1500
  png_height <- 800

  top_ylim <- max(-log10(second_gwas$P))
  x_range <- NULL
  x_lab <- ""
  if (show_specific_region) {
    png_width <- 900
    x_range <- c(bp - floor(range/2), bp + floor(range/2))
    x_lab <- paste("Chromosome", chr)

    first_gwas <- gwas_region(first_gwas, chr, bp, range)
    second_gwas <- gwas_region(second_gwas, chr, bp, range)
    top_ylim <-  max(-log10(second_gwas$P))
  } else {
    first_chrs <- sort(unique(first_gwas$CHR))
    second_chrs <- sort(unique(second_gwas$CHR))
    shared_chrs <- sort(intersect(first_chrs, second_chrs))

    if (!(all(shared_chrs == first_chrs) & all(shared_chrs == second_gwas))) {
      first_gwas <- subset(first_gwas, CHR %in% shared_chrs)
      second_gwas <- subset(second_gwas, CHR %in% shared_chrs)
    }
  }
  message(paste("Number of rows for each GWAS: ", nrow(first_gwas), nrow(second_gwas)))

  grDevices::png(miami_plot_file, width = png_width, height = png_height)
  graphics::par(mfrow = c(2, 1))
  graphics::par(mar = c(0, 5, 3, 3))

  if (show_specific_region) {
    qqman::manhattan(first_gwas, main = title, xlim = x_range)
  } else {
    qqman::manhattan(first_gwas, main = title)
  }

  graphics::par(mar = c(5, 5, 3, 3))
  if (show_specific_region) {
    qqman::manhattan(second_gwas, ylim = c(top_ylim, 0), xlim = x_range, xlab = x_lab, xaxt = "n")
  } else {
    qqman::manhattan(second_gwas, ylim = c(top_ylim, 0), xlab = x_lab, xaxt = "n")
  }

  grDevices::dev.off()
}

#' @import ggplot2
#' @import ggrepel
#' @export
volcano_plot <- function(results_file, title="Volcano Plot of Results", label="EXPOSURE", num_labels=30, output_file, p_val="p.adjusted")  {
  table <- get_file_or_dataframe(results_file)

  if (!all(c("BETA", p_val) %in% names(table))) {
    stop(paste("data frame needs to have BETA and", p_val, "columns"))
  }

  table <- dplyr::mutate(table, category = dplyr::case_when(
    get({{p_val}}) > 0.05 ~ "Not Significant",
    get({{p_val}}) < 0.05 & BETA < 0 ~ "Deleterious",
    get({{p_val}}) < 0.05 & BETA > 0 ~ "Protective"
  ))

  #filter label to only showing the more 'important' labels
  important_labels <- dplyr::filter(table, get({{p_val}}) < 0.05) |>
    dplyr::arrange(dplyr::desc(-log10(get({{p_val}}) * abs(BETA))))
  important_labels <- head(important_labels, num_labels)[[label]]
  table[[label]] <- ifelse(table[[label]] %in% important_labels, table[[label]], NA)

  ggplot2::ggplot(data = table, ggplot2::aes(x = BETA , y = -log10(.data[[p_val]]), col = category, label = .data[[label]])) +
    ggplot2::geom_vline(xintercept = c(-0.1, 0.1), col = "gray", linetype = 'dashed') +
    ggplot2::geom_hline(yintercept = -log10(0.05), col = "tomato2", linetype = 'dashed') +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_colour_brewer(type="qual", palette = "Dark2") +
    ggplot2::labs(color = 'Effect', x = "BETA", y = expression("-log"[10] * " P")) +
    ggplot2::scale_x_continuous(breaks = seq(-10, 10, 0.2)) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggrepel::geom_text_repel(max.overlaps = Inf, show.legend = F)

  ggplot2::ggsave(output_file)
}

#' @import ggplot2
plot_heritability_contribution_per_ancestry <- function(heterogeneity_results_qj, output_file) {
  graph_width <- max(1000, nrow(heterogeneity_results_qj) * 100)
  plot <- tidyr::gather(as.data.frame(heterogeneity_results_qj), "key", "value", -SNP)

  ggplot2::ggplot(plot, ggplot2::aes(x=SNP, y=-log10(value))) +
    ggplot2::geom_point(ggplot2::aes(colour=key)) +
    ggplot2::scale_colour_brewer(type="qual") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggplot2::ggtitle("Contribution to heterogeneity score broken down by population") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  ggplot2::ggsave(output_file, width = graph_width, units = "px")
}
