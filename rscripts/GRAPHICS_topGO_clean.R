# Graphics and GO-term overlap for host topGO enrichment

library(topGO)
library(ggplot2)
library(Rgraphviz)
library(VennDiagram)
library(grid)

# Set the directory containing the significant GO-term TSV files.
go_dir <- "D:/RNASequencingwork/annotation and enrichment/deseq2/Diff on Diff/Venn diagram"
setwd(go_dir)

# ------------------------------------------------------------
# 1. Graphics for one topGO result
# ------------------------------------------------------------
# Select the matching topGOdata object and result object from your RData session.
# Example names in the workspace include resultTopGO_PvC_BP, resultTopGO_PvC_MF,
# resultTopGO_PvC_CC, resultTopGO_AvC_BP, resultTopGO_AvC_MF and resultTopGO_AvC_CC.
# Set GOdata and resultTopGO to the matching objects before running this section.

# GOdata <- GOdata_PvC_BP
# resultTopGO <- resultTopGO_PvC_BP
contrast_label <- "Proximal vs Control"
plot_prefix <- "PvC_BP"

if (exists("GOdata") && exists("resultTopGO")) {
  allGO <- GenTable(GOdata, weightFisher = resultTopGO, topNodes = 50)
  allGO$weightFisher <- suppressWarnings(as.numeric(as.character(allGO$weightFisher)))
  allGO$FDR <- p.adjust(allGO$weightFisher, method = "BH")
  allGO$logFDR <- -log10(pmax(allGO$FDR, .Machine$double.xmin))

  significant_GO <- subset(allGO, FDR < 0.01)

  if (nrow(significant_GO) > 0) {
    bar_plot <- ggplot(significant_GO, aes(x = reorder(Term, logFDR), y = logFDR, fill = logFDR)) +
      geom_col() +
      coord_flip() +
      scale_fill_gradient(low = "blue", high = "red") +
      labs(title = paste("Top Enriched GO Terms:", contrast_label),
           x = "GO Term", y = "-log10(FDR)") +
      theme_minimal()

    dot_plot <- ggplot(significant_GO,
                       aes(x = logFDR, y = reorder(Term, logFDR),
                           size = Annotated, color = logFDR)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "blue", high = "red") +
      labs(title = paste("GO Term Enrichment:", contrast_label),
           x = "-log10(FDR)", y = "GO Term", size = "Gene Count") +
      theme_minimal()

    print(bar_plot)
    print(dot_plot)
    ggsave(paste0("GO_BarPlot_", plot_prefix, ".png"), bar_plot,
           width = 8, height = 6, dpi = 300)
    ggsave(paste0("GO_DotPlot_", plot_prefix, ".png"), dot_plot,
           width = 8, height = 6, dpi = 300)
  } else {
    message("No GO terms passed FDR < 0.01 for the selected result.")
  }

  showSigOfNodes(GOdata, score(resultTopGO), firstSigNodes = 10, useInfo = "all")
} else {
  message("Graphics section skipped: assign the matching GOdata and resultTopGO objects first.")
}

# Note: the original heatmap section used heatmap_data, but that object was not
# created in the supplied script. The heatmap code is omitted until its input
# matrix and intended values are defined.

# ------------------------------------------------------------
# 2. Merge downregulated GO terms across BP, MF and CC for a contrast
# ------------------------------------------------------------
contrast <- "PvC"  # Change to AvC or AvP as required

bp_file <- paste0(contrast, "_TopGO_BP_Significant_downregulated_BP.tsv")
mf_file <- paste0(contrast, "_TopGO_MF_Significant_downregulated_MF.tsv")
cc_file <- paste0(contrast, "_TopGO_CC_Significant_downregulated_CC.tsv")

read_go_terms <- function(file_path) {
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(character(0))
  }

  data <- read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE)
  go_column <- grep("^GO\\.ID$|GO.ID", colnames(data), ignore.case = TRUE)
  if (length(go_column) == 0) {
    message("No GO.ID column found in: ", file_path)
    return(character(0))
  }

  terms <- unique(as.character(data[[go_column[1]]]))
  terms[!is.na(terms) & nzchar(terms)]
}

bp_terms <- read_go_terms(bp_file)
mf_terms <- read_go_terms(mf_file)
cc_terms <- read_go_terms(cc_file)
all_go_terms <- unique(c(bp_terms, mf_terms, cc_terms))

print(all_go_terms)
write.table(all_go_terms, paste0(contrast, "_merged_GO_terms.txt"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------------------------------------------------
# 3. Venn diagram of GO-term overlap among the three contrasts
# ------------------------------------------------------------
contrast_names <- c("Control vs Proximal", "Control vs Contact", "Proximal vs Contact")
merged_files <- list(
  "Control vs Proximal" = "PvC_merged_GO_terms.txt",
  "Control vs Contact" = "AvC_merged_GO_terms.txt",
  "Proximal vs Contact" = "AvP_merged_GO_terms.txt"
)

read_merged_terms <- function(file_path) {
  if (!file.exists(file_path)) {
    message("File not found: ", file_path)
    return(character(0))
  }

  # These merged files are written as one GO term per line, without a header.
  terms <- readLines(file_path, warn = FALSE)
  terms <- trimws(terms)
  unique(terms[nzchar(terms)])
}

go_sets <- lapply(merged_files, read_merged_terms)
print(sapply(go_sets, length))

if (all(lengths(go_sets) > 0)) {
  venn_plot <- venn.diagram(
    x = go_sets,
    category.names = names(merged_files),
    fill = c("red", "blue", "green"),
    alpha = 0.5,
    lwd = 2,
    cex = 1.2,
    fontface = "bold",
    cat.fontface = "bold",
    main = "Overlap of GO Terms",
    filename = NULL
  )
  grid.newpage()
  grid.draw(venn_plot)
  png("GO_terms_Venn.png", width = 1800, height = 1500, res = 300)
  grid.newpage()
  grid.draw(venn_plot)
  dev.off()
} else {
  message("Venn diagram skipped: one or more merged GO-term files are missing or empty.")
}
