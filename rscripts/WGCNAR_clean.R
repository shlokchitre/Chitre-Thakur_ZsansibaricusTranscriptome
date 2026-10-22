# WGCNA, module-DEG overlay, and module-level GO analysis
# Required objects from the existing DESeq2 session:
#   dds
#   res_Contact_vs_Control
#   res_Contact_vs_Proximal
#   res_Proximal_vs_Control
#
# The workflow below keeps the original filtering, WGCNA settings, contrasts,
# selected modules, and GO analysis approach.

# -----------------------------
# 1. Paths and packages
# -----------------------------

work_dir <- "D:/RNASequencingwork/GENIE3/WGCNA"
outdir <- file.path(work_dir, "module data")
cyto_dir <- file.path(work_dir, "WGCNA_Cytoscape_modules")
deg_dir <- file.path(work_dir, "deg overlay")
go_dir <- file.path(work_dir, "module GO")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(cyto_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "DESeq2", "WGCNA", "dplyr", "tibble", "matrixStats",
  "pheatmap", "readr", "stringr", "tidyr", "topGO",
  "ComplexHeatmap", "circlize", "grid"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install these packages before running the script: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(DESeq2)
library(WGCNA)
library(dplyr)
library(tibble)
library(matrixStats)
library(readr)
library(stringr)
library(tidyr)
library(topGO)
library(ComplexHeatmap)
library(circlize)
library(grid)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

if (!exists("dds")) {
  stop("The DESeq2 object 'dds' must already be loaded in the R session.")
}

# -----------------------------
# 2. Prepare expression matrix
# -----------------------------

dds_filt <- dds[rowSums(counts(dds) >= 10) >= 3, ]
vsd <- vst(dds_filt, blind = TRUE)
expr <- assay(vsd)

# Collapse Trinity transcript IDs to gene-level IDs by averaging transcript
# expression values belonging to the same gene.
expr_gene <- expr %>%
  as.data.frame() %>%
  rownames_to_column("tx") %>%
  mutate(gene_id = sub("_i\\d+$", "", tx)) %>%
  group_by(gene_id) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("gene_id")

# Remove low-expression genes and the lowest-variance quartile.
expr_filt <- expr_gene[rowMeans(expr_gene) > 1, , drop = FALSE]
if (nrow(expr_filt) < 2) {
  stop("Too few genes remain after the mean-expression filter.")
}

gene_vars <- rowVars(as.matrix(expr_filt))
expr_filt <- expr_filt[gene_vars > quantile(gene_vars, 0.25), , drop = FALSE]

# Keep the top 10% most variable genes, as in the original workflow.
gene_vars <- rowVars(as.matrix(expr_filt))
variance_cutoff <- quantile(gene_vars, 0.90)
expr_wgcna <- expr_filt[gene_vars >= variance_cutoff, , drop = FALSE]

write.table(
  data.frame(gene_id = rownames(expr_wgcna)),
  file.path(outdir, "WGCNA_input_genes.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# WGCNA expects samples in rows and genes in columns.
datExpr <- t(as.matrix(expr_wgcna))

# Remove any samples/genes flagged by WGCNA's quality check.
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}
if (nrow(datExpr) < 3 || ncol(datExpr) < 2) {
  stop("Too few samples or genes remain after WGCNA quality checks.")
}

# -----------------------------
# 3. Choose soft-threshold power
# -----------------------------

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

png(file.path(outdir, "WGCNA_soft_threshold.png"),
    width = 1800, height = 1200, res = 180)
plot(
  sft$fitIndices[, 1], sft$fitIndices[, 2],
  xlab = "Soft threshold (power)",
  ylab = "Scale-free topology model fit",
  type = "n"
)
text(
  sft$fitIndices[, 1], sft$fitIndices[, 2],
  labels = sft$fitIndices[, 1], col = "red"
)
dev.off()

# Retained from the original analysis. Revisit if the soft-threshold plot
# indicates that another power is more appropriate.
softPower <- 8

# -----------------------------
# 4. Build network and identify modules
# -----------------------------

net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
print(table(moduleColors))

png(file.path(outdir, "WGCNA_module_dendrogram.png"),
    width = 2000, height = 1200, res = 180)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)
dev.off()

MEs <- orderMEs(net$MEs)

# -----------------------------
# 5. Module-trait relationships
# -----------------------------

traitData <- as.data.frame(colData(dds))

# Align trait rows with the samples that survived WGCNA quality checks.
if (!all(rownames(datExpr) %in% rownames(traitData))) {
  stop("Sample names in datExpr do not match row names in colData(dds).")
}
traitData <- traitData[rownames(datExpr), , drop = FALSE]

if (!"condition" %in% colnames(traitData)) {
  stop("colData(dds) must contain a 'condition' column.")
}

traitData$condition <- factor(
  traitData$condition,
  levels = c("Control", "Proximal", "Contact")
)
traitData_num <- model.matrix(~ condition - 1, data = traitData)
rownames(traitData_num) <- rownames(traitData)

moduleTraitCor <- cor(MEs, traitData_num, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

# Convert eigengene labels (ME0, ME1, ...) to their WGCNA module colors.
me_numbers <- as.integer(sub("^ME", "", rownames(moduleTraitCor)))
me_colors <- labels2colors(me_numbers)
rownames(moduleTraitCor) <- me_colors
rownames(moduleTraitPvalue) <- me_colors

write.table(
  moduleTraitCor, file.path(outdir, "module_trait_correlations.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)
write.table(
  moduleTraitPvalue, file.path(outdir, "module_trait_pvalues.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

pheatmap::pheatmap(
  moduleTraitCor,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Module–Trait Relationships",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  cellwidth = 75,
  cellheight = 20,
  filename = file.path(outdir, "module_trait_heatmap.png"),
  width = 7,
  height = max(4, nrow(moduleTraitCor) * 0.25 + 2)
)

# -----------------------------
# 6. Export module assignments and expression
# -----------------------------

genes <- colnames(datExpr)
module_assign <- data.frame(
  gene_id = genes,
  module = unname(moduleColors[genes]),
  stringsAsFactors = FALSE
)

write.table(
  module_assign, file.path(outdir, "gene_module_assignment.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

for (mod in unique(module_assign$module)) {
  mod_genes <- module_assign$gene_id[module_assign$module == mod]

  write.table(
    data.frame(gene_id = mod_genes, module = mod),
    file.path(outdir, paste0("module_", mod, "_genes.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  write.table(
    datExpr[, mod_genes, drop = FALSE],
    file.path(outdir, paste0("module_", mod, "_expression.tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
}

write.table(
  MEs, file.path(outdir, "module_eigengenes.tsv"),
  sep = "\t", quote = FALSE, col.names = NA
)

# -----------------------------
# 7. Export module networks for Cytoscape
# -----------------------------

TOM <- TOMsimilarityFromExpr(datExpr, power = softPower, TOMType = "signed")
rownames(TOM) <- colnames(datExpr)
colnames(TOM) <- colnames(datExpr)

modules <- split(module_assign$gene_id, module_assign$module)

for (mod in names(modules)) {
  if (mod == "grey") next

  mod_genes <- modules[[mod]]
  if (length(mod_genes) < 30) next

  message("Exporting Cytoscape module: ", mod, " (", length(mod_genes), " genes)")

  moduleTOM <- TOM[mod_genes, mod_genes, drop = FALSE]

  exportNetworkToCytoscape(
    moduleTOM,
    edgeFile = file.path(cyto_dir, paste0("cyto_edges_", mod, ".tsv")),
    nodeFile = file.path(cyto_dir, paste0("cyto_nodes_", mod, ".tsv")),
    weighted = TRUE,
    threshold = 0.15,
    nodeNames = mod_genes,
    nodeAttr = rep(mod, length(mod_genes))
  )
}

# -----------------------------
# 8. Overlay DESeq2 results on modules
# -----------------------------

required_results <- c(
  "res_Contact_vs_Control",
  "res_Contact_vs_Proximal",
  "res_Proximal_vs_Control"
)
missing_results <- required_results[!vapply(required_results, exists, logical(1))]
if (length(missing_results)) {
  stop(
    "Load the following DESeq2 result objects before continuing: ",
    paste(missing_results, collapse = ", ")
  )
}

deg_AvC <- as.data.frame(res_Contact_vs_Control) %>%
  rownames_to_column("gene_id")
deg_AvP <- as.data.frame(res_Contact_vs_Proximal) %>%
  rownames_to_column("gene_id")
deg_PvC <- as.data.frame(res_Proximal_vs_Control) %>%
  rownames_to_column("gene_id")

deg_list <- list(AvC = deg_AvC, PvC = deg_PvC, AvP = deg_AvP)

module_deg <- module_assign %>%
  left_join(deg_AvP, by = "gene_id")

write.table(
  module_deg, file.path(deg_dir, "WGCNA_module_DESeq2_AvP_overlay.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

module_summary <- module_deg %>%
  filter(!is.na(padj)) %>%
  group_by(module) %>%
  summarise(
    n_genes = n(),
    n_DEGs = sum(padj < 0.05, na.rm = TRUE),
    frac_DEGs = n_DEGs / n_genes,
    up = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
    down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(frac_DEGs))

write.table(
  module_summary, file.path(deg_dir, "WGCNA_module_summary_AvP.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Extract significant DEGs from selected modules for each contrast.
selected_modules <- c("blue", "green", "purple")
module_deg_sets <- setNames(vector("list", length(selected_modules)), selected_modules)

for (mod in selected_modules) {
  module_deg_sets[[mod]] <- list()
  mod_genes <- module_assign %>%
    filter(module == mod) %>%
    pull(gene_id)

  for (contrast in names(deg_list)) {
    result <- deg_list[[contrast]] %>%
      filter(gene_id %in% mod_genes, !is.na(padj), padj < 0.05)

    module_deg_sets[[mod]][[contrast]] <- result

    write.table(
      result,
      file.path(deg_dir, paste0(mod, "_DEGs_", contrast, ".tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}

# -----------------------------
# 9. Prepare gene-to-GO mapping
# -----------------------------

gene2go_file <- paste0(
  "D:/RNASequencingwork/annotation and enrichment/Host GO/",
  "deseq2/Diff on Diff_host/Gene2GO_filtered.tsv"
)

gene2go_tx <- read.delim(gene2go_file, stringsAsFactors = FALSE)

required_go_cols <- c("transcript_id", "GOs")
if (!all(required_go_cols %in% colnames(gene2go_tx))) {
  stop("Gene2GO file must contain columns named transcript_id and GOs.")
}

gene2go_gene <- gene2go_tx %>%
  mutate(gene_id = str_replace(transcript_id, "_i\\d+$", "")) %>%
  select(gene_id, GOs) %>%
  separate_rows(GOs, sep = "[,;|]") %>%
  mutate(GOs = trimws(GOs)) %>%
  filter(!is.na(GOs), GOs != "") %>%
  distinct(gene_id, GOs)

write.table(
  gene2go_gene, file.path(go_dir, "Gene2GO_geneLevel_wgcna.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

geneID2GO <- split(gene2go_gene$GOs, gene2go_gene$gene_id)
gene_universe <- colnames(datExpr)

# -----------------------------
# 10. Run TopGO for selected modules
# -----------------------------

run_topgo_module <- function(module_name) {
  message("Running TopGO for module: ", module_name)

  module_genes <- module_assign %>%
    filter(module == module_name) %>%
    pull(gene_id)

  geneList <- factor(as.integer(gene_universe %in% module_genes))
  names(geneList) <- gene_universe

  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = geneList,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO,
    nodeSize = 10
  )

  resultFisher <- runTest(
    GOdata,
    algorithm = "classic",
    statistic = "fisher"
  )

  go_table <- GenTable(
    GOdata,
    classicFisher = resultFisher,
    orderBy = "classicFisher",
    ranksOf = "classicFisher",
    topNodes = 30
  )

  write.table(
    go_table,
    file.path(go_dir, paste0("TopGO_BP_module_", module_name, ".tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  go_table
}

modules_to_run <- selected_modules
topgo_results <- lapply(modules_to_run, run_topgo_module)
names(topgo_results) <- modules_to_run

# Count module DEGs associated with each GO term.
annotate_go_with_degs <- function(go_table, module_name) {
  get_deg_count <- function(go_id, contrast) {
    genes_in_go <- gene2go_gene %>%
      filter(GOs == go_id) %>%
      pull(gene_id) %>%
      unique()

    length(intersect(genes_in_go, module_deg_sets[[module_name]][[contrast]]$gene_id))
  }

  go_table$DEGs_AvC <- vapply(go_table$GO.ID, get_deg_count, numeric(1), contrast = "AvC")
  go_table$DEGs_PvC <- vapply(go_table$GO.ID, get_deg_count, numeric(1), contrast = "PvC")
  go_table$DEGs_AvP <- vapply(go_table$GO.ID, get_deg_count, numeric(1), contrast = "AvP")
  go_table
}

topgo_annotated <- lapply(names(topgo_results), function(mod) {
  annotate_go_with_degs(topgo_results[[mod]], mod)
})
names(topgo_annotated) <- names(topgo_results)

for (mod in names(topgo_annotated)) {
  write.table(
    topgo_annotated[[mod]],
    file.path(go_dir, paste0("TopGO_BP_module_", mod, "_with_DEGs.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# -----------------------------
# 11. Add up/down DEG counts to GO terms
# -----------------------------

get_directional_genes <- function(deg_table, direction) {
  if (direction == "up") {
    deg_table %>% filter(!is.na(padj), padj < 0.05, log2FoldChange > 0) %>% pull(gene_id)
  } else {
    deg_table %>% filter(!is.na(padj), padj < 0.05, log2FoldChange < 0) %>% pull(gene_id)
  }
}

deg_dir_sets <- lapply(names(module_deg_sets), function(mod) {
  x <- module_deg_sets[[mod]]
  list(
    AvC_up = get_directional_genes(x$AvC, "up"),
    AvC_down = get_directional_genes(x$AvC, "down"),
    PvC_up = get_directional_genes(x$PvC, "up"),
    PvC_down = get_directional_genes(x$PvC, "down"),
    AvP_up = get_directional_genes(x$AvP, "up"),
    AvP_down = get_directional_genes(x$AvP, "down")
  )
})
names(deg_dir_sets) <- names(module_deg_sets)

annotate_go_with_direction <- function(go_table, module_name) {
  get_counts <- function(go_id, gene_set) {
    genes_in_go <- gene2go_gene %>%
      filter(GOs == go_id) %>%
      pull(gene_id) %>%
      unique()
    length(intersect(genes_in_go, gene_set))
  }

  for (direction_name in names(deg_dir_sets[[module_name]])) {
    go_table[[direction_name]] <- vapply(
      go_table$GO.ID,
      get_counts,
      numeric(1),
      gene_set = deg_dir_sets[[module_name]][[direction_name]]
    )
  }

  go_table
}

topgo_final <- lapply(names(topgo_results), function(mod) {
  annotate_go_with_direction(topgo_results[[mod]], mod)
})
names(topgo_final) <- names(topgo_results)

for (mod in names(topgo_final)) {
  write.table(
    topgo_final[[mod]],
    file.path(go_dir, paste0("TopGO_BP_module_", mod, "_with_DEGs_direction.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# -----------------------------
# 12. Export REVIGO-compatible up/down GO lists
# -----------------------------

make_revigo_up_down <- function(go_table, module_name) {
  p_text <- as.character(go_table$classicFisher)
  p_values <- suppressWarnings(as.numeric(p_text))
  # GenTable may display very small p-values as "< 1e-30"; use the
  # displayed upper bound for those entries in the REVIGO input.
  below_limit <- grepl("^\\s*<", p_text)
  if (any(below_limit)) {
    p_values[below_limit] <- suppressWarnings(as.numeric(
      gsub("[^0-9eE.+-]", "", p_text[below_limit])
    ))
  }

  revigo_base <- data.frame(
    GO_term = go_table$GO.ID,
    p_value = p_values,
    stringsAsFactors = FALSE
  )

  up_rows <- (go_table$AvC_up + go_table$PvC_up + go_table$AvP_up) > 0
  down_rows <- (go_table$AvC_down + go_table$PvC_down + go_table$AvP_down) > 0

  write.table(
    revigo_base[up_rows, , drop = FALSE],
    file.path(go_dir, paste0("REVIGO_BP_", module_name, "_UP.txt")),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
  write.table(
    revigo_base[down_rows, , drop = FALSE],
    file.path(go_dir, paste0("REVIGO_BP_", module_name, "_DOWN.txt")),
    sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
  )
}

for (mod in names(topgo_final)) {
  make_revigo_up_down(topgo_final[[mod]], mod)
}

# -----------------------------
# 13. Build signed GO regulation matrices
# -----------------------------

make_go_heatmap_matrix <- function(go_table) {
  go_table %>%
    select(
      GO.ID, Term,
      AvC_up, AvC_down,
      PvC_up, PvC_down,
      AvP_up, AvP_down
    ) %>%
    mutate(
      AvC = AvC_up - AvC_down,
      PvC = PvC_up - PvC_down,
      AvP = AvP_up - AvP_down
    ) %>%
    select(GO.ID, Term, AvC, PvC, AvP) %>%
    column_to_rownames("GO.ID")
}

go_matrices <- lapply(topgo_final, make_go_heatmap_matrix)

# Retain GO terms with an absolute up-minus-down score greater than 2
# in at least one contrast, as in the original workflow.
go_matrices <- lapply(go_matrices, function(mat) {
  scores <- as.matrix(mat[, c("AvC", "PvC", "AvP"), drop = FALSE])
  mat[apply(abs(scores), 1, max) > 2, , drop = FALSE]
})

common_cols <- c("AvC", "PvC", "AvP")
go_matrices <- lapply(go_matrices, function(mat) {
  mat[, common_cols, drop = FALSE]
})

for (mod in names(go_matrices)) {
  write.table(
    go_matrices[[mod]],
    file.path(go_dir, paste0("GO_heatmap_", mod, ".tsv")),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA
  )
}

# -----------------------------
# 14. Combined module GO heatmap
# -----------------------------

col_fun <- colorRamp2(c(-50, 0, 50), c("red", "white", "blue"))

make_module_heatmap <- function(mat, module_name, show_legend = FALSE) {
  Heatmap(
    mat,
    name = "GO regulation",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_side = "top",
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 0,
    width = unit(35, "mm"),
    height = unit(max(20, nrow(mat) * 3), "mm"),
    column_title = paste(toupper(substr(module_name, 1, 1)),
                         substring(module_name, 2), "module"),
    show_heatmap_legend = show_legend,
    heatmap_legend_param = list(title = "Up - Down")
  )
}

if (all(c("blue", "green", "purple") %in% names(go_matrices))) {
  ht_blue <- make_module_heatmap(go_matrices$blue, "blue", TRUE)
  ht_green <- make_module_heatmap(go_matrices$green, "green", FALSE)
  ht_purple <- make_module_heatmap(go_matrices$purple, "purple", FALSE)
  ht_all <- ht_blue + ht_green + ht_purple

  png(file.path(go_dir, "GO_regulation_heatmap_modules.png"),
      width = 3000, height = 1800, res = 200)
  draw(ht_all, heatmap_legend_side = "right", merge_legends = TRUE)
  dev.off()

  pdf(file.path(go_dir, "GO_regulation_heatmap_modules.pdf"),
      width = 14, height = 10)
  draw(ht_all, heatmap_legend_side = "right", merge_legends = TRUE)
  dev.off()
} else {
  warning("Combined heatmap skipped: blue, green, and purple matrices are not all available.")
}

message("WGCNA and module-level GO workflow completed. Outputs are in: ", work_dir)
