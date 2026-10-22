# GO term overlap and heatmap preparation
# Host transcriptome

library(dplyr)
library(readr)

# Main folders
go_dir <- "D:/RNASequencingwork/Genes and GO Result FILTERED/Host Genes and GO/GO Terms filtered"
all_dir <- file.path(go_dir, "ALL")
non_shared_dir <- file.path(go_dir, "NON SHARED")
shared_dir <- file.path(go_dir, "SHARED_EXACTLY_2_PAIRWISE")

dir.create(non_shared_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(shared_dir, showWarnings = FALSE, recursive = TRUE)


# Read the GO ID and p-value columns from the up/down files and combine them.
read_go_pair <- function(contrast, ontology) {
  up_file <- file.path(
    all_dir, paste0(contrast, "_TopGO_", ontology, "_Significant_upregulated.txt")
  )
  down_file <- file.path(
    all_dir, paste0(contrast, "_TopGO_", ontology, "_Significant_downregulated.txt")
  )

  read_one <- function(path) {
    if (!file.exists(path)) stop("File not found: ", path)

    x <- read.delim(path, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE)
    if (ncol(x) < 2) stop("Expected at least two columns in: ", path)

    x <- x[, 1:2, drop = FALSE]
    names(x) <- c("GO.ID", "PValue")
    x$GO.ID <- trimws(as.character(x$GO.ID))
    x$PValue <- suppressWarnings(as.numeric(as.character(x$PValue)))
    x <- x[!is.na(x$GO.ID) & x$GO.ID != "", , drop = FALSE]
    x
  }

  bind_rows(read_one(up_file), read_one(down_file)) |>
    distinct(GO.ID, .keep_all = TRUE)
}


# 1. GO terms found in only one contrast
contrasts <- c("AvC", "AvP", "PvC")

# The original analysis checks unique Cellular Component terms.
ontology <- "CC"
go_sets <- setNames(
  lapply(contrasts, read_go_pair, ontology = ontology),
  contrasts
)

for (contrast in contrasts) {
    other_contrasts <- setdiff(contrasts, contrast)
    other_ids <- unique(unlist(lapply(go_sets[other_contrasts], `[[`, "GO.ID")))

    unique_terms <- go_sets[[contrast]] |>
      filter(!GO.ID %in% other_ids)

    write.table(
      unique_terms,
      file.path(non_shared_dir, paste0(contrast, "_NonShared_", ontology, "_GO.txt")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
}


# 2. GO terms shared by exactly two contrasts
# The original analysis checks Biological Process terms here.
ontology <- "BP"
go_sets <- setNames(
  lapply(contrasts, read_go_pair, ontology = ontology),
  contrasts
)

pairs <- list(
    c("AvC", "AvP"),
    c("AvC", "PvC"),
    c("AvP", "PvC")
  )

  for (pair in pairs) {
    third <- setdiff(contrasts, pair)
    shared_ids <- intersect(go_sets[[pair[1]]]$GO.ID, go_sets[[pair[2]]]$GO.ID)
    shared_ids <- setdiff(shared_ids, go_sets[[third]]$GO.ID)

    # Keep the p-value from the first contrast in the pair.
    shared_terms <- go_sets[[pair[1]]] |>
      filter(GO.ID %in% shared_ids)

    out_name <- paste0(pair[1], "_", pair[2], "_", ontology, "_Shared_Only.txt")
    write.table(
      shared_terms, file.path(shared_dir, out_name),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }


# 3. Fill the heatmap table with signed p-values from the GO result files.
# Positive values are upregulated; negative values are downregulated.
heatmap_file <- file.path(go_dir, "hostGOheatmap.csv")
if (file.exists(heatmap_file)) {
  heatmap_data <- read.csv(heatmap_file, stringsAsFactors = FALSE,
                           check.names = FALSE)
  if (!"GO.Term" %in% names(heatmap_data)) {
    stop("hostGOheatmap.csv must contain a column named 'GO.Term'.")
  }

  contrast_files <- list.files(
    all_dir,
    pattern = "TopGO_(BP|MF|CC)_Significant_(upregulated|downregulated)\\.txt$",
    full.names = TRUE
  )

  # The contrast label is the first three characters of each filename.
  read_heatmap_result <- function(path) {
    x <- read.delim(path, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE)
    if (ncol(x) < 2) return(NULL)

    x <- x[, 1:2, drop = FALSE]
    names(x) <- c("GO.Term", "p_value")
    x$GO.Term <- trimws(as.character(x$GO.Term))
    x$p_value <- suppressWarnings(as.numeric(as.character(x$p_value)))

    direction <- if (grepl("_downregulated\\.txt$", path)) -1 else 1
    x$Contrast <- substr(basename(path), 1, 3)
    x$signed_p <- direction * x$p_value
    x
  }

  contrast_data <- bind_rows(lapply(contrast_files, read_heatmap_result))

  for (contrast in contrasts) {
    if (!contrast %in% names(heatmap_data)) heatmap_data[[contrast]] <- 0
    heatmap_data[[contrast]] <- 0
  }

  for (i in seq_len(nrow(heatmap_data))) {
    go_id <- trimws(as.character(heatmap_data$GO.Term[i]))

    for (contrast in contrasts) {
      match_row <- contrast_data |>
        filter(GO.Term == go_id, Contrast == contrast)

      if (nrow(match_row)) {
        heatmap_data[[contrast]][i] <- match_row$signed_p[1]
      }
    }
  }

  # Set Regulation only when a term has a non-zero value in the results.
  if (!"Regulation" %in% names(heatmap_data)) heatmap_data$Regulation <- NA_character_

  for (i in seq_len(nrow(heatmap_data))) {
    values <- unlist(heatmap_data[i, contrasts], use.names = FALSE)
    present <- values[values != 0]

    if (length(present)) {
      heatmap_data$Regulation[i] <- if (present[1] > 0) "upregulated" else "downregulated"
    } else {
      heatmap_data$Regulation[i] <- NA_character_
    }
  }

  preferred_order <- c("Go.Process", "Sub.process", "GO.Term", "Regulation", contrasts)
  heatmap_data <- heatmap_data[, c(
    intersect(preferred_order, names(heatmap_data)),
    setdiff(names(heatmap_data), preferred_order)
  ), drop = FALSE]

  write.csv(
    heatmap_data,
    file.path(go_dir, "heatmap_filled_with_signed_pvalues.csv"),
    row.names = FALSE, na = ""
  )
} else {
  message("Heatmap input not found; skipped: ", heatmap_file)
}


# 4. Add GO term descriptions to the heatmap table.
# Install GO.db separately if it is not already installed.
description_input <- "D:/RNASequencingwork/Genes and GO Result FILTERED/Heatmap Time/final files/GO_heatmap_host.csv"

if (file.exists(description_input) &&
    requireNamespace("GO.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {

  go_data <- readr::read_csv(description_input, show_col_types = FALSE)

  if (!"GO.Term" %in% names(go_data)) {
    stop("GO_heatmap_host.csv must contain a column named 'GO.Term'.")
  }

  go_data$GO.Term <- trimws(as.character(go_data$GO.Term))
  go_data$Description <- vapply(go_data$GO.Term, function(go_id) {
    term <- GO.db::GOTERM[[go_id]]
    if (is.null(term)) NA_character_ else AnnotationDbi::Term(term)
  }, character(1))

  readr::write_csv(
    go_data,
    file.path(dirname(description_input), "GO_heatmap_host_with_descriptions.csv")
  )

} else if (!file.exists(description_input)) {
  message("GO description input not found; skipped: ", description_input)
} else {
  message("Install GO.db and AnnotationDbi to add GO term descriptions.")
}

message("Finished processing host GO term files.")
