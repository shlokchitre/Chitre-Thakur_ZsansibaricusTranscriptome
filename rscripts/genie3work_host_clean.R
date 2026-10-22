# Host transcriptomics: GENIE3, WGCNA, DESeq2 overlays and module GO enrichment
# Run in an R session where `dds` and the three DESeq2 result objects are available,
# or set DDS_RDS_PATH and DESeq result paths below.

# -------------------------------------------------------------------------
# 0. Configuration
# -------------------------------------------------------------------------

project_dir <- "D:/RNASequencingwork/GENIE3"
setwd(project_dir)

annotation_file <- "Trinotate_host_unique.tsv"
gene2go_file <- "D:/RNASequencingwork/annotation and enrichment/Host GO/deseq2/Diff on Diff_host/Gene2GO_filtered.tsv"

# Optional RDS containing a DESeqDataSet. Leave NULL if `dds` is already loaded.
DDS_RDS_PATH <- NULL

# Optional paths to saved DESeq2 result tables. Leave NULL to use objects in R.
DESEQ_ALL_FILE <- "D:/RNASequencingwork/RESTART ENRICHMENT ANALYSIS/Deseq results host/shrinked lfc/DESeq2_shrunk_results/DESeq2_all_levels_unshrunk.csv"

# Output folders
genie3_dir <- file.path(project_dir, "host_genie3")
wgcna_dir <- file.path(project_dir, "WGCNA_host")
module_dir <- file.path(wgcna_dir, "module_data")
cyto_dir <- file.path(wgcna_dir, "Cytoscape_modules")
overlay_dir <- file.path(wgcna_dir, "deg_overlay")
go_dir <- file.path(wgcna_dir, "module_GO")

dir.create(genie3_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wgcna_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cyto_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(overlay_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(go_dir, showWarnings = FALSE, recursive = TRUE)

# Load required packages (install missing packages separately if needed).
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

if (!is.null(DDS_RDS_PATH)) {
  dds <- readRDS(DDS_RDS_PATH)
}
if (!exists("dds") || !methods::is(dds, "DESeqDataSet")) {
  stop("A DESeq2 object named `dds` must be loaded, or DDS_RDS_PATH must point to it.")
}

# -------------------------------------------------------------------------
# 1. Identify putative host transcription factors from Trinotate annotation
# -------------------------------------------------------------------------

anno <- read_tsv(annotation_file, show_col_types = FALSE, col_types = cols(.default = col_character()))
names(anno) <- sub("^#gene_id$", "gene_id", names(anno))

if (!"gene_id" %in% names(anno)) {
  stop("The host Trinotate file must contain a `#gene_id` or `gene_id` column.")
}

get_col <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit)) hit[1] else NA_character_
}

col_gene <- get_col(anno, "gene_id")
col_trans <- get_col(anno, "transcript_id")
col_pfam <- get_col(anno, "Pfam")
col_pfams <- get_col(anno, "PFAMs")
col_go_pf <- get_col(anno, "gene_ontology_Pfam")
col_gos <- get_col(anno, "GOs")
col_descr <- get_col(anno, "Description")
col_pref <- get_col(anno, "Preferred_name")
col_keggko <- get_col(anno, "KEGG_ko")
col_ko <- get_col(anno, "KO")

get_text <- function(data, col) {
  if (is.na(col)) rep("", nrow(data)) else tolower(replace_na(as.character(data[[col]]), ""))
}

anno_search <- tibble(
  gene_id = anno[[col_gene]],
  transcript_id = if (is.na(col_trans)) "" else anno[[col_trans]],
  pfam = get_text(anno, col_pfam),
  pfams = get_text(anno, col_pfams),
  go_pf = get_text(anno, col_go_pf),
  gos = get_text(anno, col_gos),
  descr = get_text(anno, col_descr),
  pref = get_text(anno, col_pref),
  keggko = get_text(anno, col_keggko),
  ko = get_text(anno, col_ko)
)

pfam_keywords <- c(
  "zf-c2h2", "zinc finger", "homeobox", "homeodomain", "bzip",
  "basic leucine zipper", "bhlh", "helix-loop-helix", "hlh", "myb",
  "myb-like", "sant", "wrky", "hsf", "gata", "t-box", "tbx", "arid",
  "forkhead", "forkhead box", "fox", "ets", "mads", "nac", "hmg_box",
  "znf", "zf-", "dnabind", "dna-binding"
)
descr_keywords <- c(
  "transcription factor", "transcriptional regulator", "dna-binding",
  "dna binding", "transcription regulator", "transcriptional activator",
  "transcriptional repressor", "sequence-specific dna binding", "tf "
)
go_keywords <- c(
  "go:0003700", "transcription factor activity", "regulation of transcription",
  "dna-binding transcription", "go:0006355", "go:0001071"
)
kegg_keywords <- c("k15602", "k09422", "k14454", "k09171", "k09180", "k09182")

pfam_regex <- regex(paste(pfam_keywords, collapse = "|"), ignore_case = TRUE)
descr_regex <- regex(paste(descr_keywords, collapse = "|"), ignore_case = TRUE)
go_regex <- regex(paste(go_keywords, collapse = "|"), ignore_case = TRUE)
kegg_regex <- regex(paste(kegg_keywords, collapse = "|"), ignore_case = TRUE)

find_tf_reasons <- function(row) {
  reasons <- character()
  if (str_detect(row$pfam, pfam_regex)) reasons <- c(reasons, "PFAM_match")
  if (str_detect(row$pfams, pfam_regex)) reasons <- c(reasons, "PFAMs_match")
  if (str_detect(row$go_pf, go_regex)) reasons <- c(reasons, "GO_PFAM_match")
  if (str_detect(row$gos, go_regex)) reasons <- c(reasons, "GOs_match")
  if (str_detect(row$descr, descr_regex)) reasons <- c(reasons, "Description_text")
  if (str_detect(row$pref, descr_regex)) reasons <- c(reasons, "Preferred_name_text")
  if (str_detect(row$keggko, kegg_regex)) reasons <- c(reasons, "KEGG_ko_match")
  if (str_detect(row$ko, kegg_regex)) reasons <- c(reasons, "KO_match")
  unique(reasons)
}

reason_list <- lapply(seq_len(nrow(anno_search)), function(i) find_tf_reasons(anno_search[i, ]))
is_tf <- lengths(reason_list) > 0
tf_ids <- unique(anno_search$gene_id[is_tf])
writeLines(tf_ids, file.path(genie3_dir, "host_TF_list.txt"))

tf_details <- anno[match(tf_ids, anno[[col_gene]]), , drop = FALSE]
tf_details$TF_match_reasons <- vapply(tf_ids, function(id) {
  idx <- which(anno_search$gene_id == id)
  paste(unique(unlist(reason_list[idx])), collapse = ";")
}, character(1))
write_csv(tf_details, file.path(genie3_dir, "host_TF_candidates_with_reasons.csv"))
message("Host genes: ", nrow(anno), "; putative TF candidates: ", length(tf_ids))

# -------------------------------------------------------------------------
# 2. Prepare expression matrix for GENIE3
# -------------------------------------------------------------------------

dds_filt <- dds[rowSums(counts(dds) >= 10) >= 3, ]
vsd <- vst(dds_filt, blind = FALSE)
expr_matrix <- assay(vsd)

# Collapse Trinity isoforms to gene-level IDs by averaging their VST values.
expr_gene <- as.data.frame(expr_matrix) %>%
  rownames_to_column("transcript_id") %>%
  mutate(gene_id = sub("_i[0-9]+$", "", transcript_id)) %>%
  group_by(gene_id) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

keep_expr <- rowMeans(expr_gene) > 1
expr_filtered <- expr_gene[keep_expr, , drop = FALSE]

gene_variance <- apply(expr_filtered, 1, var, na.rm = TRUE)
top_n <- 10000
if (sum(is.finite(gene_variance)) > top_n) {
  selected_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(top_n)]
  expr_for_genie <- expr_filtered[selected_genes, , drop = FALSE]
} else {
  expr_for_genie <- expr_filtered
}

write_tsv(
  as.data.frame(expr_for_genie) %>% rownames_to_column("gene_id"),
  file.path(genie3_dir, "host_expression_matrix_for_GENIE3.tsv")
)

# -------------------------------------------------------------------------
# 3. Run GENIE3 and export high-weight links / hub subnet
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(GENIE3)
  library(parallel)
  library(igraph)
})

tf_list <- readLines(file.path(genie3_dir, "host_TF_list.txt"), warn = FALSE)
tf_present <- intersect(tf_list, rownames(expr_for_genie))
writeLines(tf_present, file.path(genie3_dir, "host_TF_list_filtered_for_GENIE3.txt"))
message("TF candidates: ", length(tf_list), "; present in expression matrix: ", length(tf_present))

if (length(tf_present) == 0) {
  warning("No TF candidates matched the expression matrix. GENIE3 was skipped.")
} else {
  n_cores <- min(12, parallel::detectCores())
  set.seed(42)
  weight_matrix <- GENIE3(
    expr_for_genie,
    regulators = tf_present,
    nTrees = 1000,
    nCores = n_cores
  )

  links <- getLinkList(weight_matrix, threshold = 0.001)
  write_tsv(links, file.path(genie3_dir, "host_GENIE3_links.tsv"))

  top_links <- links %>%
    arrange(desc(weight)) %>%
    slice_head(prop = 0.01)
  write_tsv(top_links, file.path(genie3_dir, "host_GENIE3_top1pct_links.tsv"))

  net_graph <- graph_from_data_frame(top_links, directed = TRUE)
  hub_degree <- sort(degree(net_graph, mode = "out"), decreasing = TRUE)
  write_tsv(
    tibble(gene_id = names(hub_degree), out_degree = as.integer(hub_degree)),
    file.path(genie3_dir, "host_GENIE3_outdegree.tsv")
  )

  hub_genes <- names(hub_degree)[seq_len(min(5, length(hub_degree)))]
  hub_network <- top_links %>% filter(regulatoryGene %in% hub_genes)
  write_tsv(hub_network, file.path(genie3_dir, "host_GENIE3_top_hubs_subnet.tsv"))
}

# -------------------------------------------------------------------------
# 4. Filter and deduplicate host Trinotate annotation to genes in `dds`
# -------------------------------------------------------------------------

# This section assumes Trinotate's #gene_id values correspond to rownames(dds).
trinotate_header <- names(data.table::fread(annotation_file, nrows = 0, sep = "\t"))
if (!"#gene_id" %in% trinotate_header && !"gene_id" %in% trinotate_header) {
  stop("Could not find #gene_id/gene_id in the host Trinotate file.")
}
id_col <- if ("#gene_id" %in% trinotate_header) "#gene_id" else "gene_id"

keep_cols <- c(
  id_col, "transcript_id", "sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit",
  "Pfam", "gene_ontology_Pfam", "seed_ortholog", "eggNOG_OGs", "COG_category",
  "Description", "Preferred_name", "GOs", "KEGG_ko", "KEGG_Pathway",
  "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
  "CAZy", "BiGG_Reaction", "PFAMs", "KO"
)
keep_cols <- intersect(unique(keep_cols), trinotate_header)

trinotate <- data.table::fread(
  annotation_file, select = keep_cols, sep = "\t", header = TRUE,
  showProgress = TRUE, data.table = TRUE
)
dds_ids <- rownames(dds)
trinotate[[id_col]] <- trimws(as.character(trinotate[[id_col]]))
common_ids <- sum(trinotate[[id_col]] %in% dds_ids)

if (common_ids == 0) {
  strip_version <- function(x) sub("\\.[0-9]+$", "", x)
  dds_ids_stripped <- strip_version(dds_ids)
  trinotate_ids_stripped <- strip_version(trinotate[[id_col]])
  keep <- trinotate_ids_stripped %in% dds_ids_stripped
  if (!any(keep)) {
    stop("No Trinotate IDs match rownames(dds), even after stripping version suffixes.")
  }
  trinotate <- trinotate[keep]
  trinotate[, gene_id_match := trinotate_ids_stripped[keep]]
} else {
  trinotate <- trinotate[trinotate[[id_col]] %in% dds_ids]
  trinotate[, gene_id_match := get(id_col)]
}

# Remove exact duplicate rows, then merge remaining annotation values by matching ID.
trinotate <- unique(trinotate)
annotation_cols <- setdiff(names(trinotate), c(id_col, "gene_id_match"))
trinotate_dedup <- trinotate[, lapply(.SD, function(x) {
  values <- unique(na.omit(as.character(x)))
  values <- values[nzchar(values)]
  if (!length(values)) NA_character_ else paste(values, collapse = "|")
}), by = "gene_id_match", .SDcols = annotation_cols]
data.table::setnames(trinotate_dedup, "gene_id_match", "gene_id")
data.table::fwrite(
  trinotate_dedup,
  file.path(genie3_dir, "host_Trinotate_filtered_deduplicated.tsv"),
  sep = "\t", quote = FALSE, na = "NA"
)
message("Filtered host Trinotate rows: ", nrow(trinotate_dedup))

# -------------------------------------------------------------------------
# 5. Export log2FC and seed values for diffusion-expression analysis
# -------------------------------------------------------------------------

if (file.exists(DESEQ_ALL_FILE)) {
  deseq_all <- read_csv(DESEQ_ALL_FILE, show_col_types = FALSE)
} else if (exists("res_all")) {
  deseq_all <- as.data.frame(res_all) %>% rownames_to_column("gene_id")
} else {
  warning("No overall DESeq2 results found; diffusion-expression inputs were skipped.")
  deseq_all <- NULL
}

if (!is.null(deseq_all)) {
  if (!"gene_id" %in% names(deseq_all)) {
    deseq_all <- deseq_all %>% rownames_to_column("gene_id")
  }
  write_tsv(
    deseq_all %>% select(gene_id, log2FoldChange) %>% filter(!is.na(log2FoldChange)),
    file.path(overlay_dir, "host_All_log2fc.tsv")
  )
  seed <- deseq_all %>%
    mutate(SeedValue = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ 1,
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ -1,
      TRUE ~ 0
    )) %>%
    filter(SeedValue != 0) %>%
    select(gene_id, SeedValue)
  write_tsv(seed, file.path(overlay_dir, "host_All_seed.tsv"))
  write_tsv(
    deseq_all %>% select(gene_id, padj) %>% filter(!is.na(padj)),
    file.path(overlay_dir, "host_All_sig.tsv")
  )
}

# -------------------------------------------------------------------------
# 6. WGCNA: expression filtering, modules, trait heatmap and exports
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(WGCNA)
  library(matrixStats)
  library(pheatmap)
})
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Keep transcripts with counts >=10 in at least 3 samples; collapse isoforms above.
wgcna_counts <- counts(dds)
wgcna_counts <- wgcna_counts[rowSums(wgcna_counts >= 10) >= 3, , drop = FALSE]
wgcna_vsd <- vst(dds[rownames(wgcna_counts), ], blind = TRUE)
wgcna_expr <- assay(wgcna_vsd)

expr_gene_wgcna <- as.data.frame(wgcna_expr) %>%
  rownames_to_column("transcript_id") %>%
  mutate(gene_id = sub("_i[0-9]+$", "", transcript_id)) %>%
  group_by(gene_id) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

expr_filtered_wgcna <- expr_gene_wgcna[rowMeans(expr_gene_wgcna) > 1, , drop = FALSE]
gene_vars <- rowVars(expr_filtered_wgcna)
expr_filtered_wgcna <- expr_filtered_wgcna[
  gene_vars > quantile(gene_vars, 0.25, na.rm = TRUE), , drop = FALSE
]
gene_vars <- rowVars(expr_filtered_wgcna)
expr_wgcna <- expr_filtered_wgcna[
  gene_vars >= quantile(gene_vars, 0.90, na.rm = TRUE), , drop = FALSE
]

datExpr <- t(expr_wgcna)
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
png(file.path(wgcna_dir, "host_WGCNA_soft_threshold.png"), width = 1200, height = 800, res = 150)
plot(sft$fitIndices[, 1], sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology model fit",
     type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 2], labels = sft$fitIndices[, 1], col = "red")
dev.off()

softPower <- 8
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
write_tsv(
  tibble(gene_id = colnames(datExpr), module = moduleColors),
  file.path(module_dir, "host_gene_module_assignment.tsv")
)

png(file.path(wgcna_dir, "host_WGCNA_module_dendrogram.png"), width = 1400, height = 800, res = 150)
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

# Save per-module gene lists and expression matrices.
module_genes <- split(colnames(datExpr), moduleColors)
for (mod in names(module_genes)) {
  genes_in_module <- module_genes[[mod]]
  write_tsv(
    tibble(gene_id = genes_in_module, module = mod),
    file.path(module_dir, paste0("host_module_", mod, "_genes.tsv"))
  )
  write_tsv(
    as.data.frame(datExpr[, genes_in_module, drop = FALSE]) %>% rownames_to_column("sample"),
    file.path(module_dir, paste0("host_module_", mod, "_expression.tsv"))
  )
}
write_tsv(
  as.data.frame(net$MEs) %>% rownames_to_column("sample"),
  file.path(module_dir, "host_module_eigengenes.tsv")
)

trait_data <- as.data.frame(colData(dds))
if ("condition" %in% names(trait_data)) {
  trait_data_num <- model.matrix(~ condition - 1, data = trait_data)
  rownames(trait_data_num) <- rownames(trait_data)
  trait_data_num <- trait_data_num[rownames(datExpr), , drop = FALSE]
  module_trait_cor <- cor(net$MEs, trait_data_num, use = "pairwise.complete.obs")
  module_trait_p <- corPvalueStudent(module_trait_cor, nSamples = nrow(datExpr))

  pheatmap(
    module_trait_cor,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = "Host WGCNA module–trait relationships",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    filename = file.path(wgcna_dir, "host_module_trait_heatmap.png"),
    width = 10, height = 7
  )
  write_tsv(
    as.data.frame(module_trait_cor) %>% rownames_to_column("module_eigengene"),
    file.path(wgcna_dir, "host_module_trait_correlations.tsv")
  )
  write_tsv(
    as.data.frame(module_trait_p) %>% rownames_to_column("module_eigengene"),
    file.path(wgcna_dir, "host_module_trait_pvalues.tsv")
  )
}

# Cytoscape export: omit grey and modules with fewer than 30 genes.
TOM <- TOMsimilarityFromExpr(datExpr, power = softPower, TOMType = "signed")
dimnames(TOM) <- list(colnames(datExpr), colnames(datExpr))
for (mod in names(module_genes)) {
  genes_in_module <- module_genes[[mod]]
  if (mod == "grey" || length(genes_in_module) < 30) next
  message("Exporting host WGCNA module ", mod, " (", length(genes_in_module), " genes)")
  module_tom <- TOM[genes_in_module, genes_in_module, drop = FALSE]
  exportNetworkToCytoscape(
    module_tom,
    edgeFile = file.path(cyto_dir, paste0("host_cyto_edges_", mod, ".tsv")),
    nodeFile = file.path(cyto_dir, paste0("host_cyto_nodes_", mod, ".tsv")),
    weighted = TRUE,
    threshold = 0.15,
    nodeNames = genes_in_module,
    nodeAttr = rep(mod, length(genes_in_module))
  )
}
rm(TOM)

# -------------------------------------------------------------------------
# 7. Overlay DESeq2 contrasts on WGCNA modules
# -------------------------------------------------------------------------

# Load result objects from the current R session, or set these to data frames.
get_deg_table <- function(object_name, contrast_name) {
  if (!exists(object_name, envir = .GlobalEnv)) {
    warning("DESeq2 result object `", object_name, "` not found; skipping ", contrast_name, ".")
    return(NULL)
  }
  obj <- get(object_name, envir = .GlobalEnv)
  as.data.frame(obj) %>% rownames_to_column("gene_id")
}

deg_tables <- list(
  AvC = get_deg_table("res_Contact_vs_Control", "AvC"),
  AvP = get_deg_table("res_Contact_vs_Proximal", "AvP"),
  PvC = get_deg_table("res_Proximal_vs_Control", "PvC")
)
module_assign <- tibble(gene_id = colnames(datExpr), module = moduleColors)

for (contrast in names(deg_tables)) {
  deg <- deg_tables[[contrast]]
  if (is.null(deg)) next
  # Match DESeq2 transcript IDs to gene-level WGCNA IDs when necessary.
  deg <- deg %>% mutate(gene_id = sub("_i[0-9]+$", "", gene_id))
  module_deg <- module_assign %>% left_join(deg, by = "gene_id")
  write_tsv(module_deg, file.path(overlay_dir, paste0("host_WGCNA_module_DESeq2_", contrast, "_overlay.tsv")))

  module_summary <- module_deg %>%
    group_by(module) %>%
    summarise(
      n_genes = n(),
      n_tested = sum(!is.na(padj)),
      n_DEGs = sum(padj < 0.05, na.rm = TRUE),
      frac_DEGs = ifelse(n_tested > 0, n_DEGs / n_tested, NA_real_),
      up = sum(padj < 0.05 & log2FoldChange > 0, na.rm = TRUE),
      down = sum(padj < 0.05 & log2FoldChange < 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(frac_DEGs))
  write_tsv(module_summary, file.path(overlay_dir, paste0("host_WGCNA_module_summary_", contrast, ".tsv")))
}

# Selected modules for module-specific DEG lists and GO enrichment.
modules_to_analyse <- c("blue", "green", "purple")
deg_sets <- list()
deg_direction_sets <- list()

for (mod in modules_to_analyse) {
  module_genes_mod <- module_assign %>% filter(module == mod) %>% pull(gene_id)
  deg_sets[[mod]] <- list()
  deg_direction_sets[[mod]] <- list()

  for (contrast in names(deg_tables)) {
    deg <- deg_tables[[contrast]]
    if (is.null(deg)) next
    deg <- deg %>% mutate(gene_id = sub("_i[0-9]+$", "", gene_id))
    module_deg <- deg %>%
      filter(gene_id %in% module_genes_mod, !is.na(padj), padj < 0.05)
    write_tsv(
      module_deg,
      file.path(overlay_dir, paste0("host_", mod, "_DEGs_", contrast, ".tsv"))
    )
    deg_sets[[mod]][[contrast]] <- module_deg$gene_id
    deg_direction_sets[[mod]][[paste0(contrast, "_up")]] <- module_deg %>%
      filter(log2FoldChange > 0) %>% pull(gene_id)
    deg_direction_sets[[mod]][[paste0(contrast, "_down")]] <- module_deg %>%
      filter(log2FoldChange < 0) %>% pull(gene_id)
  }
}

# -------------------------------------------------------------------------
# 8. Module-level TopGO BP enrichment and contrast DEG counts
# -------------------------------------------------------------------------

if (file.exists(gene2go_file)) {
  suppressPackageStartupMessages(library(topGO))

  gene2go_tx <- read_tsv(gene2go_file, show_col_types = FALSE)
  if (!all(c("transcript_id", "GOs") %in% names(gene2go_tx))) {
    stop("Gene2GO input must contain columns named `transcript_id` and `GOs`.")
  }

  gene2go_gene <- gene2go_tx %>%
    mutate(gene_id = sub("_i[0-9]+$", "", transcript_id)) %>%
    separate_rows(GOs, sep = "[,;|]") %>%
    mutate(GOs = trimws(GOs)) %>%
    filter(!is.na(GOs), nzchar(GOs)) %>%
    distinct(gene_id, GOs)

  write_tsv(gene2go_gene, file.path(go_dir, "host_Gene2GO_gene_level.tsv"))
  geneID2GO <- split(gene2go_gene$GOs, gene2go_gene$gene_id)
  gene_universe <- intersect(colnames(datExpr), names(geneID2GO))

  run_topgo_module <- function(module_name) {
    module_genes <- intersect(
      module_assign %>% filter(module == module_name) %>% pull(gene_id),
      gene_universe
    )
    gene_list <- factor(as.integer(gene_universe %in% module_genes))
    names(gene_list) <- gene_universe

    go_data <- new(
      "topGOdata",
      ontology = "BP",
      allGenes = gene_list,
      annot = annFUN.gene2GO,
      gene2GO = geneID2GO,
      nodeSize = 10
    )
    fisher <- runTest(go_data, algorithm = "classic", statistic = "fisher")
    go_table <- GenTable(
      go_data,
      classicFisher = fisher,
      orderBy = "classicFisher",
      ranksOf = "classicFisher",
      topNodes = 30
    )

    # Count module DEGs belonging to each enriched GO term, by contrast and direction.
    for (contrast in names(deg_tables)) {
      for (direction in c("up", "down")) {
        set_name <- paste0(contrast, "_", direction)
        selected_ids <- deg_direction_sets[[module_name]][[set_name]]
        if (is.null(selected_ids)) selected_ids <- character()
        go_table[[set_name]] <- vapply(go_table$GO.ID, function(go_id) {
          genes_for_go <- gene2go_gene %>% filter(GOs == go_id) %>% pull(gene_id) %>% unique()
          length(intersect(genes_for_go, selected_ids))
        }, integer(1))
      }
    }

    write_tsv(go_table, file.path(go_dir, paste0("host_TopGO_BP_module_", module_name, ".tsv")))
    go_table
  }

  topgo_results <- lapply(modules_to_analyse, run_topgo_module)
  names(topgo_results) <- modules_to_analyse

  # REVIGO-compatible UP/DOWN files and signed GO count matrices.
  make_revigo_files <- function(go_table, module_name) {
    up_cols <- intersect(c("AvC_up", "PvC_up", "AvP_up"), names(go_table))
    down_cols <- intersect(c("AvC_down", "PvC_down", "AvP_down"), names(go_table))
    up_sum <- if (length(up_cols)) rowSums(go_table[, up_cols, drop = FALSE]) else rep(0, nrow(go_table))
    down_sum <- if (length(down_cols)) rowSums(go_table[, down_cols, drop = FALSE]) else rep(0, nrow(go_table))

    revigo_up <- go_table[up_sum > 0, c("GO.ID", "classicFisher"), drop = FALSE]
    revigo_down <- go_table[down_sum > 0, c("GO.ID", "classicFisher"), drop = FALSE]
    names(revigo_up) <- c("GO_term", "p_value")
    names(revigo_down) <- c("GO_term", "p_value")
    write_tsv(revigo_up, file.path(go_dir, paste0("host_REVIGO_BP_", module_name, "_UP.txt")))
    write_tsv(revigo_down, file.path(go_dir, paste0("host_REVIGO_BP_", module_name, "_DOWN.txt")))

    go_matrix <- go_table %>%
      transmute(
        GO.ID, Term,
        AvC = AvC_up - AvC_down,
        PvC = PvC_up - PvC_down,
        AvP = AvP_up - AvP_down
      ) %>%
      filter(pmax(abs(AvC), abs(PvC), abs(AvP), na.rm = TRUE) > 2)

    write_tsv(go_matrix, file.path(go_dir, paste0("host_GO_heatmap_", module_name, ".tsv")))
    go_matrix
  }

  go_matrices <- lapply(names(topgo_results), function(mod) {
    make_revigo_files(topgo_results[[mod]], mod)
  })
  names(go_matrices) <- names(topgo_results)

  # Plot one combined heatmap, keeping GO terms separate within each module.
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })

  common_cols <- c("AvC", "PvC", "AvP")
  go_matrices <- lapply(go_matrices, function(x) {
    x <- as.data.frame(x)
    x <- x[, c("GO.ID", "Term", common_cols), drop = FALSE]
    x
  })

  make_heatmap <- function(mat, module_name, show_legend = FALSE) {
    if (nrow(mat) == 0) return(NULL)
    rownames(mat) <- make.unique(paste(mat$GO.ID, mat$Term, sep = " | "))
    Heatmap(
      as.matrix(mat[, common_cols, drop = FALSE]),
      name = "Up - Down",
      col = colorRamp2(c(-50, 0, 50), c("red", "white", "blue")),
      cluster_rows = FALSE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontsize = 7),
      column_names_side = "top",
      column_names_gp = gpar(fontsize = 10, fontface = "bold"),
      column_names_rot = 0,
      width = unit(35, "mm"),
      column_title = paste(tools::toTitleCase(module_name), "module"),
      show_heatmap_legend = show_legend
    )
  }

  heatmaps <- lapply(names(go_matrices), function(mod) {
    make_heatmap(go_matrices[[mod]], mod, show_legend = mod == names(go_matrices)[1])
  })
  heatmaps <- Filter(Negate(is.null), heatmaps)
  if (length(heatmaps)) {
    png(file.path(go_dir, "host_module_GO_heatmaps.png"), width = 1800, height = 1400, res = 180)
    draw(Reduce(`+`, heatmaps), heatmap_legend_side = "right", merge_legends = TRUE)
    dev.off()
  }
} else {
  warning("Gene2GO file not found; module TopGO analysis was skipped: ", gene2go_file)
}

message("Host-only workflow completed. Outputs are under: ", project_dir)
