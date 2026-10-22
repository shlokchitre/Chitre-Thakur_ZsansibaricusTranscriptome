# Venom transcript differential-expression analysis
# Inputs: venomszsan.csv and the three DESeq2 LFC-shrunken result tables.
# Run from the VENOM EXPRESSION directory, or update work_dir below.

work_dir <- "D:/RNASequencingwork/VENOM EXPRESSION"
setwd(work_dir)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

# 1. Clean transcript IDs --------------------------------------------------
venom <- read_csv("venomszsan.csv", show_col_types = FALSE) %>%
  mutate(gene_id = sub("_ORF.*", "", gene_id),
         gene_id = sub("_i[0-9]+$", "", gene_id))

write_csv(venom, "venomtranscripts_clean.csv")

# Keep one row per transcript ID and retain its venom type.
ids <- venom %>%
  select(gene_id, Type) %>%
  distinct(gene_id, .keep_all = TRUE)

# 2. Read DESeq2 LFC-shrunken results ------------------------------------
AvC <- read_csv("DESeq2_Contact_vs_Control_LFCshrunk.csv", show_col_types = FALSE)
AvP <- read_csv("DESeq2_Contact_vs_Proximal_LFCshrunk.csv", show_col_types = FALSE)
PvC <- read_csv("DESeq2_Proximal_vs_Control_LFCshrunk.csv", show_col_types = FALSE)

# Extract the matching venom transcripts from each contrast.
AvC_f <- AvC %>%
  filter(gene_id %in% ids$gene_id) %>%
  select(gene_id, log2FoldChange, padj) %>%
  rename(AvC_log2FC = log2FoldChange, AvC_padj = padj)

AvP_f <- AvP %>%
  filter(gene_id %in% ids$gene_id) %>%
  select(gene_id, log2FoldChange, padj) %>%
  rename(AvP_log2FC = log2FoldChange, AvP_padj = padj)

PvC_f <- PvC %>%
  filter(gene_id %in% ids$gene_id) %>%
  select(gene_id, log2FoldChange, padj) %>%
  rename(PvC_log2FC = log2FoldChange, PvC_padj = padj)

final <- ids %>%
  left_join(AvC_f, by = "gene_id") %>%
  left_join(AvP_f, by = "gene_id") %>%
  left_join(PvC_f, by = "gene_id") %>%
  select(gene_id, Type, everything())

write_csv(final, "filtered_DESeq_results.csv")

# 3. Remove transcripts with no available fold change in any contrast ------
fc_cols <- c("AvC_log2FC", "AvP_log2FC", "PvC_log2FC")
padj_cols <- c("AvC_padj", "AvP_padj", "PvC_padj")

results_clean <- final %>%
  filter(if_any(all_of(fc_cols), ~ !is.na(.x)))

write_csv(results_clean, "filtered_DESeq_results_noNA.csv")

# 4. Keep transcripts significant in at least one contrast -----------------
# Nonsignificant or unavailable contrasts are set to zero for visualization.
sig <- results_clean %>%
  filter(if_any(all_of(padj_cols), ~ !is.na(.x) & .x <= 0.05)) %>%
  mutate(
    AvC_log2FC = if_else(is.na(AvC_padj) | AvC_padj > 0.05, 0, AvC_log2FC),
    AvP_log2FC = if_else(is.na(AvP_padj) | AvP_padj > 0.05, 0, AvP_log2FC),
    PvC_log2FC = if_else(is.na(PvC_padj) | PvC_padj > 0.05, 0, PvC_log2FC)
  )

write_csv(sig, "filtered_DESeq_significant.csv")

# 5. Heatmap of significant venom transcripts ------------------------------
if (nrow(sig) > 0) {
  heat_df <- sig %>% arrange(Type)
  heat_matrix <- as.matrix(heat_df[, fc_cols])
  rownames(heat_matrix) <- make.unique(heat_df$gene_id)
  heat_matrix[heat_matrix > 5] <- 5
  heat_matrix[heat_matrix < -5] <- -5

  annotation_row <- data.frame(Type = heat_df$Type)
  rownames(annotation_row) <- rownames(heat_matrix)

  pheatmap(
    heat_matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = annotation_row,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "Differential expression of predicted venom transcripts",
    angle_col = 0,
    cellwidth = 36,
    filename = "venom_transcript_heatmap.png",
    width = 9,
    height = max(5, min(18, nrow(heat_matrix) * 0.18 + 2))
  )
} else {
  message("No venom transcripts met the adjusted-p-value cutoff of 0.05.")
}

# 6. Prepare long-format data for trajectory plots -------------------------
if (nrow(sig) > 0) {
  contrast_labels <- c(
    PvC_log2FC = "Proximal vs Control",
    AvC_log2FC = "Contact vs Control",
    AvP_log2FC = "Contact vs Proximal"
  )

  df_long <- sig %>%
    select(gene_id, Type, all_of(fc_cols)) %>%
    pivot_longer(
      cols = all_of(fc_cols),
      names_to = "Contrast",
      values_to = "log2FC"
    ) %>%
    mutate(
      Contrast = recode(Contrast, !!!contrast_labels),
      Contrast = factor(
        Contrast,
        levels = c("Proximal vs Control", "Contact vs Control", "Contact vs Proximal")
      )
    )

  # Individual transcript trajectories, faceted by venom type.
  trajectory_plot <- ggplot(
    df_long,
    aes(x = Contrast, y = log2FC, group = gene_id, color = Type)
  ) +
    geom_line(alpha = 0.7, linewidth = 0.7) +
    geom_point(size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    facet_wrap(~Type, scales = "free_y") +
    theme_minimal() +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Expression trajectories of venom transcripts across contrasts",
      x = "",
      y = "log2 fold change"
    )

  ggsave("venom_transcript_trajectories.png", trajectory_plot,
         width = 10, height = 7, dpi = 300)

  # Mean trajectory for each venom type.
  mean_lines <- df_long %>%
    group_by(Type, Contrast) %>%
    summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop")

  mean_plot <- ggplot(df_long, aes(x = Contrast, y = log2FC, group = gene_id)) +
    geom_line(color = "grey40", alpha = 0.25, linewidth = 0.6) +
    geom_line(
      data = mean_lines,
      aes(y = mean_log2FC, group = Type, color = Type),
      linewidth = 1.5
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Mean expression trajectories by venom type",
         x = "", y = "log2 fold change", color = "Venom type")

  ggsave("venom_mean_trajectories.png", mean_plot,
         width = 9, height = 6, dpi = 300)

  # Significance-weighted trajectory by venom type.
  fc_long <- sig %>%
    select(gene_id, Type, all_of(fc_cols)) %>%
    pivot_longer(cols = all_of(fc_cols), names_to = "Contrast", values_to = "log2FC")

  padj_long <- sig %>%
    select(gene_id, all_of(padj_cols)) %>%
    pivot_longer(cols = all_of(padj_cols), names_to = "Contrast", values_to = "padj") %>%
    mutate(Contrast = sub("_padj$", "_log2FC", Contrast))

  weight_long <- fc_long %>%
    left_join(padj_long, by = c("gene_id", "Contrast")) %>%
    mutate(
      Contrast = recode(Contrast, !!!contrast_labels),
      Contrast = factor(
        Contrast,
        levels = c("Proximal vs Control", "Contact vs Control", "Contact vs Proximal")
      ),
      weight = if_else(!is.na(padj) & padj > 0 & padj <= 1,
                       -log10(padj), NA_real_)
    )

  weighted_lines <- weight_long %>%
    group_by(Type, Contrast) %>%
    summarise(
      weighted_log2FC = if (sum(is.finite(weight) & weight > 0, na.rm = TRUE) == 0) {
        NA_real_
      } else {
        sum(log2FC * weight, na.rm = TRUE) / sum(weight[is.finite(weight) & weight > 0], na.rm = TRUE)
      },
      .groups = "drop"
    )

  weighted_plot <- ggplot(
    df_long,
    aes(x = Contrast, y = log2FC, group = gene_id, color = Type)
  ) +
    geom_line(alpha = 0.25, linewidth = 0.6) +
    geom_line(
      data = weighted_lines,
      aes(y = weighted_log2FC, group = Type, color = Type),
      linewidth = 1.8,
      na.rm = TRUE
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right"
    ) +
    labs(title = "Significance-weighted venom transcript trajectories",
         x = "", y = "Weighted log2 fold change", color = "Venom type")

  ggsave("venom_weighted_trajectories.png", weighted_plot,
         width = 9, height = 6, dpi = 300)
}

message("Venom transcript analysis complete. Output tables and plots were saved in: ", getwd())
