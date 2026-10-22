# =============================================================================
# HOST TRANSCRIPTOME: SALMON -> tximport -> DESeq2 -> plots -> GO enrichment
# =============================================================================
# This script performs gene-level DESeq2 using the host-specific Trinotate
# transcript-to-gene mapping. Edit the paths and metadata before running.
#
# Expected metadata columns: sample, condition
# Expected conditions: Control, Proximal, Contact
# Expected Trinotate columns: transcript_id, #gene_id; GOs is needed for GO.
#
# Outputs are written to output_dir. No packages are installed by this script.
# =============================================================================

# ------------------------------- 1. SETTINGS --------------------------------

base_dir <- "D:/RNASequencingwork/Salmon and ahead/salmon_quant_results"
samples_path <- file.path(base_dir, "samples_zsan_u.csv")
trinotate_path <- "D:/RNASequencingwork/annotation and enrichment/Host GO/Trinotate_host_unique.tsv"
output_dir <- file.path(base_dir, "DESeq2_host_clean_results")

quant_pattern <- "_uquant\\.sf$"  # Change to "quant\\.sf$" if files are named quant.sf
counts_from_abundance <- "lengthScaledTPM"

# Expression prefilter: retain genes with >=10 counts in at least 20% of samples,
# with a minimum of 2 samples.
min_count <- 10
min_sample_prop <- 0.20

# Differential-expression thresholds used for plots and DEG lists.
padj_cutoff <- 0.05
lfc_cutoff <- 1

# Plot settings
heatmap_top_n <- 50
plot_width <- 8
plot_height <- 6
plot_dpi <- 300

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "results"), showWarnings = FALSE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "GO_enrichment"), showWarnings = FALSE)

# ----------------------------- 2. PACKAGES ----------------------------------

required_packages <- c(
  "readr", "dplyr", "stringr", "tximport", "DESeq2",
  "ggplot2", "pheatmap", "RColorBrewer"
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

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tximport)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

# -------------------------- 3. READ METADATA --------------------------------

if (!file.exists(samples_path)) stop("Metadata file not found: ", samples_path)
if (!file.exists(trinotate_path)) stop("Host Trinotate file not found: ", trinotate_path)

samples <- readr::read_csv(samples_path, show_col_types = FALSE) |>
  as.data.frame()

required_metadata <- c("sample", "condition")
if (!all(required_metadata %in% names(samples))) {
  stop("Metadata must contain columns named: sample, condition.")
}
samples$sample <- as.character(samples$sample)
samples$condition <- as.character(samples$condition)

if (anyNA(samples$sample) || any(samples$sample == "")) {
  stop("Metadata contains missing or empty sample names.")
}
if (anyDuplicated(samples$sample)) stop("Sample names in metadata must be unique.")
if (anyNA(samples$condition) || any(samples$condition == "")) {
  stop("Metadata contains missing or empty condition labels.")
}

expected_conditions <- c("Control", "Proximal", "Contact")
unexpected_conditions <- setdiff(unique(samples$condition), expected_conditions)
if (length(unexpected_conditions)) {
  stop(
    "Unexpected condition(s): ", paste(unexpected_conditions, collapse = ", "),
    ". Expected: ", paste(expected_conditions, collapse = ", ")
  )
}
if (!all(expected_conditions %in% samples$condition)) {
  stop("All three conditions are required: Control, Proximal, Contact.")
}

samples$condition <- factor(samples$condition, levels = expected_conditions)
rownames(samples) <- samples$sample

# ---------------------- 4. DISCOVER SALMON FILES -----------------------------

all_quant_files <- list.files(
  path = base_dir,
  pattern = quant_pattern,
  recursive = TRUE,
  full.names = TRUE
)

if (!length(all_quant_files)) {
  stop("No Salmon quantification files found. Check base_dir and quant_pattern.")
}

sample_from_path <- function(path) {
  x <- basename(path)
  x <- sub("\\.sf$", "", x)
  x <- sub("_u?quant$", "", x)
  x
}

quant_samples <- vapply(all_quant_files, sample_from_path, character(1))
if (anyDuplicated(quant_samples)) {
  dup <- unique(quant_samples[duplicated(quant_samples)])
  stop("Multiple quantification files map to sample(s): ", paste(dup, collapse = ", "))
}

files <- setNames(all_quant_files, quant_samples)
missing_quant <- setdiff(samples$sample, names(files))
if (length(missing_quant)) {
  stop("No quantification file found for: ", paste(missing_quant, collapse = ", "))
}

# Use metadata order consistently throughout import, modelling, and plots.
files <- files[samples$sample]
stopifnot(identical(names(files), rownames(samples)))

message("Salmon files found for ", length(files), " metadata samples.")

# ------------------- 5. HOST TRANSCRIPT-TO-GENE MAP --------------------------

trin <- readr::read_tsv(
  trinotate_path,
  show_col_types = FALSE,
  name_repair = "minimal",
  progress = FALSE
)

if (!all(c("transcript_id", "#gene_id") %in% names(trin))) {
  stop("Host Trinotate TSV must contain columns 'transcript_id' and '#gene_id'.")
}

tx2gene <- trin |>
  transmute(
    TX = as.character(.data[["transcript_id"]]),
    GENE = as.character(.data[["#gene_id"]])
  ) |>
  filter(
    !is.na(TX), !is.na(GENE),
    TX != "", GENE != "", TX != ".", GENE != "."
  ) |>
  distinct()

if (!nrow(tx2gene)) stop("No valid transcript-to-gene mappings found.")

ambiguous_tx <- tx2gene |>
  count(TX, name = "n_genes") |>
  filter(n_genes > 1)
if (nrow(ambiguous_tx)) {
  stop(
    nrow(ambiguous_tx),
    " transcript IDs map to multiple gene IDs. Resolve the mapping before proceeding."
  )
}

readr::write_tsv(
  tx2gene,
  file.path(output_dir, "results", "tx2gene_host.tsv")
)

# ------------------------- 6. TXIMPORT + DESEQ2 ------------------------------

txi <- tximport(
  files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = counts_from_abundance
)

dds <- DESeqDataSetFromTximport(
  txi,
  colData = samples,
  design = ~ condition
)

# Retain genes with >= min_count counts in at least the specified proportion
# of samples, and never require fewer than two samples.
min_samples <- max(2L, floor(ncol(dds) * min_sample_prop))
keep <- rowSums(counts(dds) >= min_count) >= min_samples
dds <- dds[keep, ]

if (!nrow(dds)) stop("No genes remain after expression filtering.")
message("Genes retained after prefiltering: ", nrow(dds))

dds$condition <- relevel(dds$condition, ref = "Control")
dds <- DESeq(dds)

# Named list of the three requested pairwise contrasts.
contrast_specs <- list(
  Proximal_vs_Control = c("condition", "Proximal", "Control"),
  Contact_vs_Control = c("condition", "Contact", "Control"),
  Contact_vs_Proximal = c("condition", "Contact", "Proximal")
)

results_list <- lapply(
  contrast_specs,
  function(x) results(dds, contrast = x, alpha = padj_cutoff)
)
res_all <- results(dds, alpha = padj_cutoff)

# --------------------------- 7. EXPORT RESULTS -------------------------------

results_dir <- file.path(output_dir, "results")

write_result <- function(res, filename) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- df[, c("gene_id", setdiff(names(df), "gene_id"))]
  write.csv(df, file.path(results_dir, filename), row.names = FALSE, na = "")
}

write_result(res_all, "DESeq2_default_coefficient.csv")

for (nm in names(results_list)) {
  res <- results_list[[nm]]
  write_result(res, paste0("DESeq2_", nm, ".csv"))

  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  sig <- df[!is.na(df$padj) & df$padj < padj_cutoff, , drop = FALSE]
  write.csv(
    sig,
    file.path(results_dir, paste0("DESeq2_significant_", nm, ".csv")),
    row.names = FALSE,
    na = ""
  )

  # Directional DEG lists use both adjusted p-value and |log2FC| thresholds.
  sig_lfc <- df[
    !is.na(df$padj) & df$padj < padj_cutoff &
      !is.na(df$log2FoldChange) & abs(df$log2FoldChange) >= lfc_cutoff,
    , drop = FALSE
  ]
  write.csv(
    sig_lfc[sig_lfc$log2FoldChange > 0, , drop = FALSE],
    file.path(results_dir, paste0("DESeq2_upregulated_", nm, ".csv")),
    row.names = FALSE,
    na = ""
  )
  write.csv(
    sig_lfc[sig_lfc$log2FoldChange < 0, , drop = FALSE],
    file.path(results_dir, paste0("DESeq2_downregulated_", nm, ".csv")),
    row.names = FALSE,
    na = ""
  )
}

# Top up/down genes from the Control-referenced overall result table are not
# meaningful as a single contrast; export the top genes separately per contrast.
for (nm in names(results_list)) {
  df <- as.data.frame(results_list[[nm]])
  df$gene_id <- rownames(df)
  df <- df[!is.na(df$padj) & df$padj < padj_cutoff, , drop = FALSE]
  up <- df[order(df$log2FoldChange, decreasing = TRUE), , drop = FALSE]
  down <- df[order(df$log2FoldChange, decreasing = FALSE), , drop = FALSE]
  write.csv(
    head(up, 50),
    file.path(results_dir, paste0("top50_up_", nm, ".csv")),
    row.names = FALSE, na = ""
  )
  write.csv(
    head(down, 50),
    file.path(results_dir, paste0("top50_down_", nm, ".csv")),
    row.names = FALSE, na = ""
  )
}

normalized_counts <- counts(dds, normalized = TRUE)
write.table(
  data.frame(gene_id = rownames(normalized_counts), normalized_counts,
             check.names = FALSE),
  file.path(results_dir, "normalized_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

saveRDS(dds, file.path(results_dir, "dds_host_gene_level.rds"))
saveRDS(txi, file.path(results_dir, "tximport_host.rds"))
saveRDS(results_list, file.path(results_dir, "pairwise_results.rds"))
saveRDS(res_all, file.path(results_dir, "results_default_coefficient.rds"))

# ----------------------------- 8. PLOTS --------------------------------------

plots_dir <- file.path(output_dir, "plots")

# Consistent condition colours across all plots.
condition_colours <- c(
  Control = "#3B6FB6",
  Proximal = "#E69F00",
  Contact = "#009E73"
)

# 8a. PCA: variance-stabilized counts, labelled by sample.
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca_data$sample <- rownames(pca_data)

p_pca <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, colour = condition, label = sample)
) +
  geom_point(size = 3.5, alpha = 0.9) +
  geom_text(nudge_y = 0.35, size = 3, show.legend = FALSE, check_overlap = TRUE) +
  scale_colour_manual(values = condition_colours, drop = FALSE) +
  labs(
    title = "Host transcriptome PCA",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    colour = "Condition"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

ggsave(
  file.path(plots_dir, "PCA_host.png"), p_pca,
  width = plot_width, height = plot_height, dpi = plot_dpi
)

# 8b. Sample-to-sample distance heatmap.
sample_dist <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dist)
annotation_col <- data.frame(condition = samples[colnames(vsd), "condition", drop = TRUE])
rownames(annotation_col) <- colnames(vsd)

png(file.path(plots_dir, "Sample_distance_heatmap.png"),
    width = 1800, height = 1500, res = 200)
pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  clustering_distance_rows = sample_dist,
  clustering_distance_cols = sample_dist,
  border_color = NA,
  main = "Host sample-to-sample distances"
)
dev.off()

# 8c. MA and volcano plots for each pairwise contrast.
plot_volcano <- function(res, contrast_name) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$padj), , drop = FALSE]
  if (!nrow(df)) return(NULL)

  df$significance <- "Not significant"
  df$significance[df$padj < padj_cutoff & df$log2FoldChange >= lfc_cutoff] <- "Upregulated"
  df$significance[df$padj < padj_cutoff & df$log2FoldChange <= -lfc_cutoff] <- "Downregulated"
  df$neg_log10_padj <- -log10(pmax(df$padj, .Machine$double.xmin))

  ggplot(df, aes(x = log2FoldChange, y = neg_log10_padj, colour = significance)) +
    geom_point(alpha = 0.65, size = 1.3) +
    scale_colour_manual(values = c(
      "Upregulated" = "#C83E4D",
      "Downregulated" = "#2878B5",
      "Not significant" = "grey75"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", linewidth = 0.4) +
    labs(
      title = paste("Host DEGs:", gsub("_", " ", contrast_name)),
      x = "Log2 fold change",
      y = expression(-log[10]("adjusted p-value")),
      colour = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom")
}

for (nm in names(results_list)) {
  res <- results_list[[nm]]

  png(file.path(plots_dir, paste0("MA_", nm, ".png")),
      width = 1600, height = 1200, res = 200)
  plotMA(res, alpha = padj_cutoff, ylim = c(-5, 5),
         main = paste("Host:", gsub("_", " ", nm)))
  dev.off()

  p_volcano <- plot_volcano(res, nm)
  if (!is.null(p_volcano)) {
    ggsave(
      file.path(plots_dir, paste0("Volcano_", nm, ".png")),
      p_volcano, width = plot_width, height = plot_height, dpi = plot_dpi
    )
  }

  # Heatmap of the most significant genes in each contrast. Skip if none are
  # significant; never stop the entire pipeline because one contrast has none.
  df <- as.data.frame(res)
  df <- df[!is.na(df$padj) & df$padj < padj_cutoff, , drop = FALSE]
  if (nrow(df)) {
    df <- df[order(df$padj), , drop = FALSE]
    top_ids <- head(rownames(df), heatmap_top_n)
    top_ids <- intersect(top_ids, rownames(assay(vsd)))
    if (length(top_ids) >= 2) {
      heatmap_mat <- assay(vsd)[top_ids, , drop = FALSE]
      heatmap_annotation <- annotation_col
      png(
        file.path(plots_dir, paste0("Heatmap_top_DEGs_", nm, ".png")),
        width = 2000, height = 1800, res = 200
      )
      pheatmap(
        heatmap_mat,
        scale = "row",
        annotation_col = heatmap_annotation,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = length(top_ids) <= 35,
        show_colnames = TRUE,
        border_color = NA,
        color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(101),
        main = paste("Top significant host genes:", gsub("_", " ", nm))
      )
      dev.off()
    }
  } else {
    message("No significant genes for ", nm, "; skipping contrast-specific DEG heatmap.")
  }
}

# ------------------------- 9. ANNOTATION EXPORT ------------------------------

# Attach host annotation to gene-level DE results using #gene_id. This is an
# annotation join only; it does not change the DESeq2 statistics.
annotation_cols <- intersect(
  c("#gene_id", "transcript_id", "sprot_Top_BLASTX_hit",
    "sprot_Top_BLASTX_hit_evalue", "Preferred_name", "Description",
    "GOs", "KEGG_Pathway", "PFAMs", "Pfam"),
  names(trin)
)
trin_annotation <- trin[, annotation_cols, drop = FALSE]
names(trin_annotation)[names(trin_annotation) == "#gene_id"] <- "gene_id"
trin_annotation$gene_id <- as.character(trin_annotation$gene_id)

# Keep one annotation row per gene for a clean join; retain unique annotation
# strings where multiple transcript rows contribute to the same gene.
collapse_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x) & x != "" & x != "."]))
  if (!length(x)) NA_character_ else paste(x, collapse = "; ")
}
annotation_gene <- trin_annotation |>
  group_by(gene_id) |>
  summarise(
    across(everything(), collapse_unique),
    .groups = "drop"
  )

for (nm in names(results_list)) {
  df <- as.data.frame(results_list[[nm]])
  df$gene_id <- rownames(df)
  annotated <- left_join(df, annotation_gene, by = "gene_id")
  write.csv(
    annotated,
    file.path(results_dir, paste0("DESeq2_annotated_", nm, ".csv")),
    row.names = FALSE, na = ""
  )
}

# -------------------------- 10. TOPGO ENRICHMENT -----------------------------

# GO enrichment uses all tested genes as the universe, with DEGs defined by
# adjusted p-value < padj_cutoff AND |log2FC| >= lfc_cutoff. It runs separately
# for each pairwise contrast and ontology (BP, MF, CC).
#
# topGO is optional for the DESeq2 portion. If unavailable, the script reports
# this and completes without GO enrichment.

if (!requireNamespace("topGO", quietly = TRUE)) {
  message("Package 'topGO' is not installed; skipping GO enrichment.")
} else if (!("GOs" %in% names(trin))) {
  message("No 'GOs' column found in host Trinotate TSV; skipping GO enrichment.")
} else {
  suppressPackageStartupMessages(library(topGO))

  # Build gene-level GO mappings from all host transcripts belonging to a gene.
  go_rows <- trin |>
    transmute(
      gene_id = as.character(.data[["#gene_id"]]),
      GOs = as.character(.data[["GOs"]])
    ) |>
    filter(!is.na(gene_id), gene_id != "", !is.na(GOs), GOs != "", GOs != ".")

  split_go <- function(x) {
    ids <- unlist(strsplit(x, "[|,;[:space:]]+", perl = TRUE))
    unique(ids[grepl("^GO:[0-9]{7}$", ids)])
  }

  go_map <- lapply(split(go_rows$GOs, go_rows$gene_id), function(x) {
    unique(unlist(lapply(x, split_go), use.names = FALSE))
  })
  go_map <- go_map[lengths(go_map) > 0]

  for (nm in names(results_list)) {
    res <- results_list[[nm]]
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)

    # topGO universe must consist of genes tested by DESeq2 and annotated with GO.
    universe_ids <- intersect(rownames(dds), names(go_map))
    if (length(universe_ids) < 10) {
      message("Too few GO-annotated genes for ", nm, "; skipping GO enrichment.")
      next
    }

    selected <- !is.na(res_df[universe_ids, "padj"]) &
      res_df[universe_ids, "padj"] < padj_cutoff &
      !is.na(res_df[universe_ids, "log2FoldChange"]) &
      abs(res_df[universe_ids, "log2FoldChange"]) >= lfc_cutoff
    gene_list <- setNames(as.integer(selected), universe_ids)

    if (sum(gene_list) == 0) {
      message("No DEGs pass the GO selection thresholds for ", nm, "; skipping GO.")
      next
    }

    go_map_use <- go_map[universe_ids]
    for (ontology in c("BP", "MF", "CC")) {
      go_obj <- tryCatch(
        new(
          "topGOdata",
          ontology = ontology,
          allGenes = gene_list,
          geneSelectionFun = function(x) x == 1,
          annot = annFUN.gene2GO,
          gene2GO = go_map_use,
          nodeSize = 5
        ),
        error = function(e) {
          message("Could not create topGO object for ", nm, " / ", ontology,
                  ": ", conditionMessage(e))
          NULL
        }
      )
      if (is.null(go_obj)) next

      go_test <- tryCatch(
        runTest(go_obj, algorithm = "weight01", statistic = "fisher"),
        error = function(e) {
          message("topGO failed for ", nm, " / ", ontology, ": ",
                  conditionMessage(e))
          NULL
        }
      )
      if (is.null(go_test)) next

      go_table <- GenTable(
        go_obj,
        weight01Fisher = go_test,
        orderBy = "weight01Fisher",
        ranksOf = "weight01Fisher",
        topNodes = length(score(go_test))
      )

      # Add BH-adjusted p-values across terms tested by topGO in this ontology.
      raw_p <- score(go_test)
      go_table$GO_ID <- as.character(go_table$GO.ID)
      go_table$weight01Fisher <- as.numeric(go_table$weight01Fisher)
      go_table$BH_adjusted_pvalue <- p.adjust(raw_p, method = "BH")[
        match(go_table$GO_ID, names(raw_p))
      ]

      write.csv(
        go_table,
        file.path(
          output_dir, "GO_enrichment",
          paste0("topGO_", nm, "_", ontology, ".csv")
        ),
        row.names = FALSE, na = ""
      )

      # Plot only the strongest terms for a readable summary; output table retains all terms.
      plot_df <- go_table[
        !is.na(go_table$BH_adjusted_pvalue) &
          go_table$BH_adjusted_pvalue < padj_cutoff,
        , drop = FALSE
      ]
      if (nrow(plot_df)) {
        plot_df <- plot_df[order(plot_df$BH_adjusted_pvalue), , drop = FALSE]
        plot_df <- head(plot_df, 15)
        plot_df$Term <- factor(plot_df$Term, levels = rev(plot_df$Term))
        p_go <- ggplot(
          plot_df,
          aes(x = -log10(pmax(BH_adjusted_pvalue, .Machine$double.xmin)), y = Term)
        ) +
          geom_col(fill = "#3B6FB6", width = 0.75) +
          labs(
            title = paste("Host GO enrichment:", gsub("_", " ", nm), ontology),
            x = expression(-log[10]("BH-adjusted p-value")),
            y = NULL
          ) +
          theme_classic(base_size = 11)
        ggsave(
          file.path(
            output_dir, "GO_enrichment",
            paste0("GO_dotbar_", nm, "_", ontology, ".png")
          ),
          p_go, width = 9, height = max(4.5, 0.28 * nrow(plot_df) + 1.5),
          dpi = plot_dpi
        )
      }
    }
  }
}

# ----------------------------- 11. SUMMARY -----------------------------------

summary_lines <- c(
  "Host DESeq2 pipeline completed.",
  paste("Samples:", ncol(dds)),
  paste("Genes tested after filtering:", nrow(dds)),
  paste("Conditions:", paste(levels(dds$condition), collapse = ", ")),
  paste("Adjusted p-value cutoff:", padj_cutoff),
  paste("Absolute log2 fold-change cutoff for DEG lists/GO:", lfc_cutoff),
  paste("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
)
writeLines(summary_lines, file.path(output_dir, "pipeline_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
