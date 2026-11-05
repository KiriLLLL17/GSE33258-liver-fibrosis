# =============================================================================
# GSE33258 — Differential expression across liver fibrosis stages (HCV)
#
# Pipeline: GEO download -> sample grouping by METAVIR stage -> limma
#           -> probe annotation (GPL14795) -> tables, volcano plots, Venn diagrams
#
# Output is written to ./results_GSE33258/{tables,plots}
# The STRING / Cytoscape (cytoHubba) steps were done interactively; their
# exported results are stored in ../network/
# =============================================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

pkgs <- c("GEOquery", "limma", "EnhancedVolcano", "VennDiagram", "dplyr", "tibble")
BiocManager::install(setdiff(pkgs, rownames(installed.packages())), ask = FALSE, update = FALSE)

library(GEOquery)
library(limma)
library(EnhancedVolcano)
library(VennDiagram)
library(dplyr)
library(tibble)

cat("Working directory:", getwd(), "\n")
results_dir <- file.path(getwd(), "results_GSE33258")
tables_dir  <- file.path(results_dir, "tables")
plots_dir   <- file.path(results_dir, "plots")
dir.create(results_dir, showWarnings = FALSE)
dir.create(tables_dir,  showWarnings = FALSE)
dir.create(plots_dir,   showWarnings = FALSE)
cat("Results folder:", results_dir, "\n")

# -----------------------------------------------------------------------------
# 1) Download GSE33258 (series matrix)
# -----------------------------------------------------------------------------
gse  <- getGEO("GSE33258", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
expr <- exprs(gse)
meta <- pData(gse)

# -----------------------------------------------------------------------------
# 2) Assign samples to METAVIR fibrosis stages F1–F4
# -----------------------------------------------------------------------------
source_text <- tolower(paste(meta$title, meta$source_name_ch1, sep = " | "))
meta$group <- NA_character_
meta$group[grepl("fibrosis\\s*stage\\s*1|\\bf1\\b", source_text)] <- "F1"
meta$group[grepl("fibrosis\\s*stage\\s*2|\\bf2\\b", source_text)] <- "F2"
meta$group[grepl("fibrosis\\s*stage\\s*3|\\bf3\\b", source_text)] <- "F3"
meta$group[grepl("fibrosis\\s*stage\\s*4|\\bf4\\b", source_text)] <- "F4"

if (any(is.na(meta$group))) {
  cat("Samples that could not be assigned to a stage:\n")
  print(meta[is.na(meta$group), c("geo_accession", "title", "source_name_ch1")])
}

keep <- !is.na(meta$group)
meta <- meta[keep, , drop = FALSE]
expr <- expr[, keep, drop = FALSE]
meta$group <- factor(meta$group, levels = c("F1", "F2", "F3", "F4"))

cat("Expression matrix:", nrow(expr), "probes x", ncol(expr), "samples\n")
print(table(meta$group))

# -----------------------------------------------------------------------------
# 3) limma
# NOTE on the design of GSE33258:
# this is a two-colour microarray; exprs() holds log2(Cy5/Cy3), where Cy3
# (channel 2) is a common reference — a normal human liver biopsy.
# There is therefore no separate "Normal" group of columns.
# "Normal vs F1" is tested as: is the mean log2-ratio in F1 different from 0?
# The same applies to F4. With a no-intercept design (~0 + group) each
# coefficient is exactly that group mean.
# -----------------------------------------------------------------------------
design <- model.matrix(~0 + group, data = meta)
colnames(design) <- sub("^group", "", colnames(design))
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# -----------------------------------------------------------------------------
# 4) Probe annotation (platform GPL14795)
# -----------------------------------------------------------------------------
annot <- Table(getGEO("GPL14795", destdir = "."))
annot_subset <- annot[, c("ID", "gene symbol")]
colnames(annot_subset) <- c("ProbeID", "GeneSymbol")

prepare_results <- function(fit_obj, coef_name) {
  topTable(fit_obj, coef = coef_name, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("ProbeID") %>%
    left_join(annot_subset, by = "ProbeID") %>%
    mutate(
      GeneLabel   = ifelse(is.na(GeneSymbol) | GeneSymbol == "", ProbeID, GeneSymbol),
      Significant = adj.P.Val < 0.05 & abs(logFC) > 1
    ) %>%
    arrange(adj.P.Val)
}

res_f1 <- prepare_results(fit, "F1")
res_f4 <- prepare_results(fit, "F4")

sig_f1 <- res_f1 %>% filter(Significant)
sig_f4 <- res_f4 %>% filter(Significant)

write.csv(res_f1, file.path(tables_dir, "DE_Normal_vs_F1_all.csv"), row.names = FALSE)
write.csv(res_f4, file.path(tables_dir, "DE_Normal_vs_F4_all.csv"), row.names = FALSE)
write.csv(sig_f1, file.path(tables_dir, "DE_Normal_vs_F1_significant.csv"), row.names = FALSE)
write.csv(sig_f4, file.path(tables_dir, "DE_Normal_vs_F4_significant.csv"), row.names = FALSE)

cat("Significant probes, Normal vs F1:", nrow(sig_f1), "\n")
cat("Significant probes, Normal vs F4:", nrow(sig_f4), "\n")

# -----------------------------------------------------------------------------
# 5) Volcano plots (y axis = BH-adjusted p-value)
# -----------------------------------------------------------------------------
plot_volcano <- function(df, title_text, outfile) {
  p <- EnhancedVolcano(
    df,
    lab = df$GeneLabel,
    x = "logFC",
    y = "adj.P.Val",
    ylab = bquote(~-Log[10] ~ "adjusted P"),
    pCutoff = 0.05,
    FCcutoff = 1,
    title = title_text,
    subtitle = paste0("DEGs: ", sum(df$Significant, na.rm = TRUE)),
    selectLab = character(0),
    drawConnectors = FALSE,
    pointSize = 2,
    colAlpha = 0.7
  )
  png(outfile, width = 1600, height = 1200, res = 200)
  print(p)
  dev.off()
}

plot_volcano(res_f1, "GSE33258: Normal vs F1", file.path(plots_dir, "Volcano_Normal_vs_F1.png"))
plot_volcano(res_f4, "GSE33258: Normal vs F4", file.path(plots_dir, "Volcano_Normal_vs_F4.png"))

# -----------------------------------------------------------------------------
# 6) Venn diagrams: overlap of up- and down-regulated genes between F1 and F4
# -----------------------------------------------------------------------------
up_f1   <- unique(sig_f1$GeneLabel[sig_f1$logFC >  1])
up_f4   <- unique(sig_f4$GeneLabel[sig_f4$logFC >  1])
down_f1 <- unique(sig_f1$GeneLabel[sig_f1$logFC < -1])
down_f4 <- unique(sig_f4$GeneLabel[sig_f4$logFC < -1])

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

draw_venn <- function(sets, fills, title_text, outfile) {
  png(outfile, width = 1400, height = 1100, res = 180)
  grid::grid.draw(venn.diagram(
    x = sets, filename = NULL, fill = fills, alpha = 0.6,
    cex = 1.2, cat.cex = 1.0, cat.pos = c(-20, 20), cat.dist = 0.05,
    margin = 0.1, main = title_text
  ))
  dev.off()
}

draw_venn(list("Normal vs F1" = up_f1, "Normal vs F4" = up_f4),
          c("#66c2a5", "#fc8d62"), "Up-regulated genes",
          file.path(plots_dir, "Venn_UP_F1_vs_F4.png"))
draw_venn(list("Normal vs F1" = down_f1, "Normal vs F4" = down_f4),
          c("#8da0cb", "#e78ac3"), "Down-regulated genes",
          file.path(plots_dir, "Venn_DOWN_F1_vs_F4.png"))

cat("Done: tables, volcano plots and Venn diagrams saved to results_GSE33258/\n")
