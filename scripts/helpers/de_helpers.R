# Differential expression (limma-voom) and volcano-plot helpers.
# Requires: edgeR, limma, ggplot2, ggrepel, dplyr

# Draws the voom mean-variance trend from a voom object fitted with save.plot = TRUE.
plot_meanvar <- function(v, title) {
  plot(v$voom.xy$x, v$voom.xy$y, pch = 16, cex = 0.3,
       xlab = v$voom.xy$xlab, ylab = v$voom.xy$ylab, main = title,
       cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3)
  lines(v$voom.line, col = "red", lwd = 2)
}

# Two-group limma-voom differential expression.
#   counts:    features x samples count matrix
#   group:     group label per sample (same order as columns of counts)
#   levels:    c(reference, test)
#   coef_name: name of the test coefficient
#   block:     optional patient ID per sample; if given, within-patient
#              correlation is estimated with duplicateCorrelation and used in
#              both voom and lmFit (two passes)
# Returns list(v, fit, res, consensus).
run_voom_de <- function(counts, group, levels, coef_name, meanvar_pdf, block = NULL) {
  stopifnot(length(group) == ncol(counts), all(group %in% levels))
  group <- factor(group, levels = levels)
  design <- model.matrix(~ group)
  colnames(design) <- c("Intercept", coef_name)
  rownames(design) <- colnames(counts)

  cat("Samples:", paste(levels, table(group)[levels], sep = " = ", collapse = ", "), "\n")

  dge <- edgeR::DGEList(counts = counts)
  keep <- edgeR::filterByExpr(dge, design = design)
  cat("Features passing filterByExpr:", sum(keep), "of", length(keep), "\n")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)  # TMM

  consensus <- NA_real_
  if (is.null(block)) {
    v <- limma::voom(dge, design, save.plot = TRUE)
    fit <- limma::lmFit(v, design)
  } else {
    v <- limma::voom(dge, design, save.plot = TRUE)
    corfit <- limma::duplicateCorrelation(v, design, block = block)
    v <- limma::voom(dge, design, block = block, correlation = corfit$consensus,
                     save.plot = TRUE)
    corfit <- limma::duplicateCorrelation(v, design, block = block)
    consensus <- corfit$consensus
    cat("Within-patient consensus correlation:", round(consensus, 4), "\n")
    if (consensus <= 0) warning("Consensus correlation <= 0; blocking on patient has no effect.")
    fit <- limma::lmFit(v, design, block = block, correlation = consensus)
  }
  fit <- limma::eBayes(fit)
  res <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  res$Feature_ID <- rownames(res)

  grDevices::pdf(meanvar_pdf, width = 8, height = 6)
  plot_meanvar(v, paste("voom mean-variance:", coef_name))
  grDevices::dev.off()
  plot_meanvar(v, paste("voom mean-variance:", coef_name))

  list(v = v, fit = fit, res = res, consensus = consensus)
}

# Volcano plot. Labels either the top label_n features by moderated t in each
# direction (label_by = "t") or by adjusted p-value (label_by = "p"), restricted
# to |logFC| > label_logfc.
plot_volcano <- function(res, label_col, up_name, down_name, title, n_text,
                         adj_p, abs_logfc, label_n, label_logfc, label_by = c("t", "p")) {
  label_by <- match.arg(label_by)
  plot_df <- res %>%
    dplyr::mutate(
      Significance = dplyr::case_when(
        adj.P.Val < adj_p & logFC >  abs_logfc ~ up_name,
        adj.P.Val < adj_p & logFC < -abs_logfc ~ down_name,
        TRUE ~ "NS"
      ),
      label = if (label_by == "t") {
        ifelse((rank(t) <= label_n & logFC < -label_logfc) |
               (rank(-t) <= label_n & logFC > label_logfc), .data[[label_col]], NA_character_)
      } else {
        ifelse(rank(adj.P.Val) <= label_n & abs(logFC) > label_logfc, .data[[label_col]], NA_character_)
      }
    )
  n_up <- sum(plot_df$Significance == up_name)
  n_dn <- sum(plot_df$Significance == down_name)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
    ggplot2::geom_point(alpha = 0.5, size = 1.5) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3.5, fontface = "bold",
                             max.overlaps = 25, segment.alpha = 0.4, na.rm = TRUE,
                             show.legend = FALSE) +
    ggplot2::scale_color_manual(values = stats::setNames(c("#D73027", "#4575B4", "grey70"),
                                                         c(up_name, down_name, "NS"))) +
    ggplot2::geom_hline(yintercept = -log10(adj_p), linetype = "dashed", alpha = 0.4) +
    ggplot2::geom_vline(xintercept = c(-abs_logfc, abs_logfc), linetype = "dashed", alpha = 0.4) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.title    = ggplot2::element_text(face = "bold"),
      axis.text     = ggplot2::element_text(size = 12, face = "bold"),
      legend.title  = ggplot2::element_text(size = 12, face = "italic"),
      legend.text   = ggplot2::element_text(face = "bold"),
      plot.title    = ggplot2::element_text(size = 16, face = "bold"),
      plot.subtitle = ggplot2::element_text(face = "bold"),
      plot.caption  = ggplot2::element_text(size = 9, face = "italic", color = "grey30")
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("%s | Up: %d, Down: %d", n_text, n_up, n_dn),
      caption = sprintf("Cutoffs: adj. P < %s; |log2 FC| > %s", adj_p, abs_logfc),
      x = expression(bold(log[2] ~ "Fold Change")),
      y = expression(bold(-log[10] ~ "Adjusted P-value"))
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = 16, size = 3, alpha = 1)))
}
