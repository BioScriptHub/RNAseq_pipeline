#!/usr/bin/env Rscript

# DESeq2 differential expression analysis.
# Input:
#   08_counts/gene_counts.txt
#   09_metadata/metadata.tsv
#   09_metadata/contrasts.tsv
# Output:
#   normalized counts, VST matrix, PCA, correlation heatmap, DEG tables, volcano plots, MA plots

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = default)
  if (identical(value, "")) {
    return(default)
  }
  value
}

stop_with_help <- function(error, hint) {
  stop(
    paste0(
      "\n[ERROR] ", error,
      "\n[HINT] ", hint,
      "\n"
    ),
    call. = FALSE
  )
}

runtime_cache_dir <- file.path(getwd(), ".runtime_cache")
dir.create(runtime_cache_dir, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(XDG_CACHE_HOME = runtime_cache_dir)
if (identical(Sys.getenv("TZ"), "")) {
  Sys.setenv(TZ = "Asia/Shanghai")
}

required_packages <- c(
  "DESeq2",
  "tidyverse",
  "pheatmap",
  "ggplot2",
  "ggrepel",
  "RColorBrewer"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop_with_help(
    paste0("Missing R packages: ", paste(missing_packages, collapse = ", ")),
    "Create the conda environment with: mamba env create -f environment.yml, then run conda activate rnaseq_pipeline."
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
})

metadata_file <- get_env("METADATA", "09_metadata/metadata.tsv")
contrasts_file <- get_env("CONTRASTS", "09_metadata/contrasts.tsv")
count_file <- file.path(get_env("COUNT_DIR", "08_counts"), "gene_counts.txt")
out_dir <- get_env("DESEQ2_DIR", "10_deseq2")
fig_dir <- get_env("FIGURE_DIR", "11_figures")
candidate_dir <- get_env("CANDIDATE_DIR", "13_candidate_genes")

min_count <- as.integer(get_env("MIN_COUNT", "10"))
min_samples <- as.integer(get_env("MIN_SAMPLES", "3"))
padj_cutoff <- as.numeric(get_env("PADJ_CUTOFF", "0.05"))
lfc_cutoff <- as.numeric(get_env("LFC_CUTOFF", "1"))
top_n <- as.integer(get_env("TOP_N", "50"))

if (is.na(min_count) || min_count < 0) {
  stop_with_help("MIN_COUNT is invalid.", "Set MIN_COUNT to a non-negative integer in config/config.sh.")
}
if (is.na(min_samples) || min_samples < 1) {
  stop_with_help("MIN_SAMPLES is invalid.", "Set MIN_SAMPLES to a positive integer in config/config.sh.")
}
if (is.na(padj_cutoff) || padj_cutoff <= 0 || padj_cutoff > 1) {
  stop_with_help("PADJ_CUTOFF is invalid.", "Set PADJ_CUTOFF to a number greater than 0 and less than or equal to 1.")
}
if (is.na(lfc_cutoff) || lfc_cutoff < 0) {
  stop_with_help("LFC_CUTOFF is invalid.", "Set LFC_CUTOFF to a non-negative number.")
}
if (is.na(top_n) || top_n < 1) {
  stop_with_help("TOP_N is invalid.", "Set TOP_N to a positive integer.")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "comparisons"), showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(candidate_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(count_file)) {
  stop_with_help(
    paste0("Count file not found: ", count_file),
    "Run bash scripts/04_featureCounts.sh first and confirm 08_counts/gene_counts.txt exists."
  )
}
if (!file.exists(metadata_file)) {
  stop_with_help(
    paste0("Metadata file not found: ", metadata_file),
    "Run cp 09_metadata/metadata.example.tsv 09_metadata/metadata.tsv, then edit it for your samples."
  )
}
if (!file.exists(contrasts_file)) {
  stop_with_help(
    paste0("Contrasts file not found: ", contrasts_file),
    "Run cp 09_metadata/contrasts.example.tsv 09_metadata/contrasts.tsv, then edit numerator and denominator."
  )
}

message("[INFO] Reading featureCounts output")
raw_counts <- read.delim(
  count_file,
  comment.char = "#",
  check.names = FALSE
)

required_count_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (!all(required_count_cols %in% colnames(raw_counts))) {
  stop_with_help(
    "featureCounts output does not contain required columns.",
    "Use the unmodified output from scripts/04_featureCounts.sh. Required columns are Geneid, Chr, Start, End, Strand and Length."
  )
}

count_mat <- raw_counts %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

gene_ids <- count_mat$Geneid
if (anyDuplicated(gene_ids) > 0) {
  duplicated_ids <- unique(gene_ids[duplicated(gene_ids)])
  stop_with_help(
    paste0("Duplicated gene IDs found, for example: ", paste(head(duplicated_ids, 5), collapse = ", ")),
    "Check the gene_id attribute in annotation.gtf. Each gene_id should represent one gene-level feature."
  )
}

count_mat <- count_mat %>%
  dplyr::select(-Geneid) %>%
  as.data.frame()

colnames(count_mat) <- basename(colnames(count_mat))
colnames(count_mat) <- stringr::str_replace(colnames(count_mat), "\\.sorted\\.bam$", "")
rownames(count_mat) <- gene_ids

count_mat <- as.matrix(count_mat)
storage.mode(count_mat) <- "integer"

if (any(is.na(count_mat))) {
  stop_with_help(
    "NA values found in count matrix.",
    "Check 08_counts/gene_counts.txt. It should contain integer raw counts from featureCounts, not edited values."
  )
}

message("[INFO] Reading metadata")
metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(metadata) == 0) {
  stop_with_help(
    "metadata.tsv contains no sample rows.",
    "Add at least two biological replicates per condition. Three replicates per condition are recommended."
  )
}

required_meta_cols <- c("sample_id", "condition", "fastq_1", "fastq_2")
if (!all(required_meta_cols %in% colnames(metadata))) {
  stop_with_help(
    paste0("metadata.tsv has invalid columns: ", paste(colnames(metadata), collapse = ", ")),
    paste0("The header must be exactly: ", paste(required_meta_cols, collapse = "\t"))
  )
}

if ("batch" %in% colnames(metadata)) {
  stop_with_help(
    "metadata.tsv contains a batch column.",
    "This template assumes one experimental batch and uses design = ~ condition. Remove the batch column. Do not mix different batches in this template."
  )
}

if (anyDuplicated(metadata$sample_id) > 0) {
  duplicated_samples <- unique(metadata$sample_id[duplicated(metadata$sample_id)])
  stop_with_help(
    paste0("Duplicated sample_id found: ", paste(duplicated_samples, collapse = ", ")),
    "Edit metadata.tsv so every sample_id is unique."
  )
}

metadata <- metadata %>%
  dplyr::mutate(
    sample_id = as.character(sample_id),
    condition = factor(condition)
  )

if (!all(metadata$sample_id %in% colnames(count_mat))) {
  missing_samples <- setdiff(metadata$sample_id, colnames(count_mat))
  stop_with_help(
    paste0("Samples in metadata but not in count matrix: ", paste(missing_samples, collapse = ", ")),
    "Check that sample_id values match BAM names generated by HISAT2 and columns in featureCounts output."
  )
}

if (!all(colnames(count_mat) %in% metadata$sample_id)) {
  extra_samples <- setdiff(colnames(count_mat), metadata$sample_id)
  stop_with_help(
    paste0("Samples in count matrix but not in metadata: ", paste(extra_samples, collapse = ", ")),
    "Check metadata.tsv. Every count matrix column must have one matching sample_id."
  )
}

count_mat <- count_mat[, metadata$sample_id, drop = FALSE]
rownames(metadata) <- metadata$sample_id

message("[INFO] Reading contrasts")
contrasts <- read.delim(
  contrasts_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (nrow(contrasts) == 0) {
  stop_with_help(
    "contrasts.tsv contains no comparison rows.",
    "Add at least one row such as: Treatment_vs_Control<TAB>Treatment<TAB>Control."
  )
}

required_contrast_cols <- c("comparison", "numerator", "denominator")
if (!all(required_contrast_cols %in% colnames(contrasts))) {
  stop_with_help(
    paste0("contrasts.tsv has invalid columns: ", paste(colnames(contrasts), collapse = ", ")),
    paste0("The header must be exactly: ", paste(required_contrast_cols, collapse = "\t"))
  )
}

if (anyDuplicated(contrasts$comparison) > 0) {
  duplicated_comparisons <- unique(contrasts$comparison[duplicated(contrasts$comparison)])
  stop_with_help(
    paste0("Duplicated comparison names found: ", paste(duplicated_comparisons, collapse = ", ")),
    "Each comparison name must be unique because it is used as the output filename prefix."
  )
}

available_conditions <- as.character(unique(metadata$condition))
contrast_conditions <- unique(c(contrasts$numerator, contrasts$denominator))
missing_conditions <- setdiff(contrast_conditions, available_conditions)
if (length(missing_conditions) > 0) {
  stop_with_help(
    paste0("Conditions in contrasts.tsv but not in metadata: ", paste(missing_conditions, collapse = ", ")),
    "Edit contrasts.tsv so numerator and denominator exactly match condition names in metadata.tsv."
  )
}

condition_counts <- table(metadata$condition)
contrast_condition_counts <- condition_counts[contrast_conditions]
low_rep_conditions <- names(contrast_condition_counts[contrast_condition_counts < 2])
if (length(low_rep_conditions) > 0) {
  stop_with_help(
    paste0("Too few replicates for condition(s): ", paste(low_rep_conditions, collapse = ", ")),
    "Each condition used in contrasts.tsv should have at least 2 biological replicates. Three or more are recommended."
  )
}

message("[INFO] Building DESeq2 model")
design_formula <- as.formula("~ condition")

design_matrix <- model.matrix(design_formula, data = metadata)
if (qr(design_matrix)$rank < ncol(design_matrix)) {
  stop_with_help(
    "Design matrix is not full rank.",
    "Check condition groups in metadata.tsv. Each condition used in contrasts should have valid samples."
  )
}

dds <- tryCatch(
  DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = metadata,
    design = design_formula
  ),
  error = function(e) {
    stop_with_help(
      paste0("Failed to build DESeqDataSet: ", conditionMessage(e)),
      "Check that count matrix contains integer raw counts and metadata sample_id values match count columns."
    )
  }
)

min_samples <- min(min_samples, ncol(dds))
keep <- rowSums(counts(dds) >= min_count) >= min_samples
if (sum(keep) == 0) {
  stop_with_help(
    "No genes retained after low-expression filtering.",
    "Counts may be too low or MIN_COUNT/MIN_SAMPLES may be too strict. Check 08_counts/gene_counts.txt and gene_counts.txt.summary."
  )
}
dds <- dds[keep, ]

dds <- tryCatch(
  DESeq(dds),
  error = function(e) {
    stop_with_help(
      paste0("DESeq2 failed: ", conditionMessage(e)),
      "Check replicate numbers, condition names, count matrix quality, and whether all samples belong to one experimental batch."
    )
  }
)
saveRDS(dds, file.path(out_dir, "dds_rnaseq.rds"))

norm_counts <- counts(dds, normalized = TRUE)
write.csv(
  as.data.frame(norm_counts) %>% tibble::rownames_to_column("gene_id"),
  file.path(out_dir, "normalized_counts.csv"),
  row.names = FALSE
)

vsd <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) {
    message("[WARN] vst failed. Falling back to varianceStabilizingTransformation.")
    varianceStabilizingTransformation(dds, blind = FALSE)
  }
)

vsd_mat <- assay(vsd)
write.csv(
  as.data.frame(vsd_mat) %>% tibble::rownames_to_column("gene_id"),
  file.path(out_dir, "vst_expression_matrix.csv"),
  row.names = FALSE
)

annotation_col <- as.data.frame(colData(vsd)[, "condition", drop = FALSE])

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

message("[INFO] Drawing PCA")
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3, max.overlaps = 100) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(
  filename = file.path(fig_dir, "PCA_condition.pdf"),
  plot = p_pca,
  width = 7.2,
  height = 5.2
)

message("[INFO] Drawing sample correlation heatmap")
sample_cor <- cor(vsd_mat, method = "pearson")

pdf(file.path(fig_dir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100),
  border_color = NA,
  main = "Sample correlation"
)
dev.off()

plot_volcano <- function(res_df, comparison_name, output_file) {
  plot_df <- res_df %>%
    dplyr::mutate(
      neg_log10_padj = -log10(padj),
      regulation = factor(regulation, levels = c("down", "not_significant", "up")),
      change = dplyr::case_when(
        regulation == "up" ~ "Up",
        regulation == "down" ~ "Down",
        TRUE ~ "NoDiff"
      ),
      change = factor(change, levels = c("Up", "Down", "NoDiff"))
    )

  finite_max <- suppressWarnings(max(plot_df$neg_log10_padj[is.finite(plot_df$neg_log10_padj)], na.rm = TRUE))
  if (!is.finite(finite_max)) {
    finite_max <- 0
  }
  plot_df$neg_log10_padj[is.infinite(plot_df$neg_log10_padj)] <- finite_max + 1

  top_labels <- plot_df %>%
    dplyr::filter(regulation != "not_significant") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 10) %>%
    dplyr::pull(gene_id)

  plot_df <- plot_df %>%
    dplyr::mutate(label_gene = ifelse(gene_id %in% top_labels, gene_id, NA_character_))

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = regulation)) +
    geom_point(alpha = 0.7, size = 1.2, na.rm = TRUE) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", linewidth = 0.3) +
    geom_text_repel(aes(label = label_gene), size = 2.5, max.overlaps = 50, na.rm = TRUE) +
    scale_color_manual(
      values = c(
        "up" = "#B2182B",
        "down" = "#2166AC",
        "not_significant" = "grey70"
      )
    ) +
    labs(
      title = comparison_name,
      x = "log2 fold change",
      y = "-log10 adjusted P-value",
      color = "Regulation"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())

  ggsave(output_file, plot = p, width = 7, height = 5)
}

plot_ma <- function(res_df, comparison_name, output_file) {
  plot_df <- res_df %>%
    dplyr::filter(!is.na(baseMean), !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      baseMean_plot = pmax(baseMean, 1),
      significant = regulation != "not_significant"
    )

  p <- ggplot(plot_df, aes(x = baseMean_plot, y = log2FoldChange, color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    scale_x_log10() +
    scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "grey70")) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", linewidth = 0.3) +
    labs(
      title = paste0("MA plot: ", comparison_name),
      x = "mean normalized count",
      y = "log2 fold change",
      color = "Significant"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank())

  ggsave(output_file, plot = p, width = 7, height = 5)
}

all_results <- list()
sig_sets <- list()
up_sets <- list()
down_sets <- list()
summary_rows <- list()

for (i in seq_len(nrow(contrasts))) {
  comparison_name <- contrasts$comparison[i]
  numerator <- contrasts$numerator[i]
  denominator <- contrasts$denominator[i]
  comparison_safe <- safe_name(comparison_name)

  message("[INFO] Contrast: ", comparison_name)

  res <- results(
    dds,
    contrast = c("condition", numerator, denominator),
    alpha = padj_cutoff
  )

  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      comparison = comparison_name,
      regulation = dplyr::case_when(
        !is.na(padj) & padj < padj_cutoff & log2FoldChange >= lfc_cutoff ~ "up",
        !is.na(padj) & padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "down",
        TRUE ~ "not_significant"
      )
    ) %>%
    dplyr::arrange(is.na(padj), padj)

  all_results[[comparison_name]] <- res_df

  sig_df <- res_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)
  up_df <- sig_df %>% dplyr::filter(log2FoldChange >= lfc_cutoff)
  down_df <- sig_df %>% dplyr::filter(log2FoldChange <= -lfc_cutoff)

  comparison_dir <- file.path(out_dir, "comparisons")

  write.csv(res_df, file.path(comparison_dir, paste0(comparison_safe, ".all_results.csv")), row.names = FALSE)
  write.csv(sig_df, file.path(comparison_dir, paste0(comparison_safe, ".significant_DEGs.csv")), row.names = FALSE)
  write.csv(up_df, file.path(comparison_dir, paste0(comparison_safe, ".up_genes.csv")), row.names = FALSE)
  write.csv(down_df, file.path(comparison_dir, paste0(comparison_safe, ".down_genes.csv")), row.names = FALSE)

  plot_volcano(
    res_df,
    comparison_name,
    file.path(fig_dir, paste0(comparison_safe, ".volcano.pdf"))
  )

  plot_ma(
    res_df,
    comparison_name,
    file.path(fig_dir, paste0(comparison_safe, ".MA_plot.pdf"))
  )

  top_candidates <- sig_df %>%
    dplyr::mutate(
      padj_for_score = pmax(padj, .Machine$double.xmin),
      candidate_score = abs(log2FoldChange) * (-log10(padj_for_score))
    ) %>%
    dplyr::select(-padj_for_score) %>%
    dplyr::arrange(dplyr::desc(candidate_score), padj)

  write.csv(
    top_candidates,
    file.path(candidate_dir, paste0(comparison_safe, ".top_candidate_genes_ranked.csv")),
    row.names = FALSE
  )

  write.csv(
    head(top_candidates, top_n),
    file.path(candidate_dir, paste0(comparison_safe, ".top", top_n, "_candidate_genes.csv")),
    row.names = FALSE
  )

  sig_sets[[comparison_name]] <- sig_df$gene_id
  up_sets[[comparison_name]] <- up_df$gene_id
  down_sets[[comparison_name]] <- down_df$gene_id

  summary_rows[[comparison_name]] <- tibble::tibble(
    comparison = comparison_name,
    numerator = numerator,
    denominator = denominator,
    tested_genes = nrow(res_df),
    significant_DEGs = nrow(sig_df),
    up_genes = nrow(up_df),
    down_genes = nrow(down_df)
  )
}

deg_summary <- dplyr::bind_rows(summary_rows)
write.csv(deg_summary, file.path(out_dir, "DEG_summary.csv"), row.names = FALSE)

write_upset_outputs <- function(gene_sets, membership_file, plot_file, plot_title) {
  all_genes <- sort(unique(unlist(gene_sets, use.names = FALSE)))
  if (length(all_genes) > 0) {
    upset_input <- data.frame(gene_id = all_genes, stringsAsFactors = FALSE)
    for (comparison_name in names(gene_sets)) {
      upset_input[[safe_name(comparison_name)]] <- as.integer(upset_input$gene_id %in% gene_sets[[comparison_name]])
    }
  } else {
    upset_input <- data.frame(gene_id = character(), stringsAsFactors = FALSE)
  }

  write.csv(upset_input, membership_file, row.names = FALSE)

  pdf(plot_file, width = 9, height = 6)
  if (length(all_genes) > 0 && length(gene_sets) >= 2 && ncol(upset_input) > 2) {
    membership_mat <- as.matrix(upset_input[, -1, drop = FALSE])
    storage.mode(membership_mat) <- "integer"
    patterns <- apply(membership_mat, 1, paste0, collapse = "")
    pattern_counts <- sort(table(patterns), decreasing = TRUE)
    pattern_counts <- head(pattern_counts, 30)
    pattern_mat <- do.call(rbind, strsplit(names(pattern_counts), split = ""))
    pattern_mat <- apply(pattern_mat, 2, as.integer)
    if (is.null(dim(pattern_mat))) {
      pattern_mat <- matrix(pattern_mat, nrow = 1)
    }
    colnames(pattern_mat) <- colnames(membership_mat)

    old_par <- par(no.readonly = TRUE)
    layout(matrix(c(1, 2), nrow = 2), heights = c(3.3, 2.2))

    par(mar = c(1.5, 4.2, 3.5, 1))
    barplot(
      as.numeric(pattern_counts),
      names.arg = rep("", length(pattern_counts)),
      col = "#3D3D3D",
      border = NA,
      ylab = "Gene count",
      main = plot_title
    )
    grid(nx = NA, ny = NULL, col = "#E5E7EB")

    set_names <- colnames(pattern_mat)
    n_intersections <- nrow(pattern_mat)
    n_sets <- ncol(pattern_mat)
    y_positions <- rev(seq_len(n_sets))

    par(mar = c(4.4, 8.2, 0.5, 1))
    plot(
      NA,
      xlim = c(0.5, n_intersections + 0.5),
      ylim = c(0.5, n_sets + 0.5),
      xaxt = "n",
      yaxt = "n",
      xlab = "Intersection rank",
      ylab = "",
      bty = "n"
    )
    axis(1, at = seq_len(n_intersections), labels = seq_len(n_intersections), las = 2, cex.axis = 0.75)
    axis(2, at = y_positions, labels = set_names, las = 2, cex.axis = 0.85)
    abline(h = y_positions, col = "#E5E7EB", lwd = 0.8)

    for (i in seq_len(n_intersections)) {
      active <- which(pattern_mat[i, ] == 1)
      active_y <- y_positions[active]
      inactive_y <- y_positions[setdiff(seq_len(n_sets), active)]
      if (length(active_y) >= 2) {
        lines(c(i, i), range(active_y), col = "#B9433F", lwd = 1.2)
      }
      if (length(inactive_y) > 0) {
        points(rep(i, length(inactive_y)), inactive_y, pch = 16, col = "#D1D5DB", cex = 1.0)
      }
      points(rep(i, length(active_y)), active_y, pch = 16, col = "#B9433F", cex = 1.35)
    }
    mtext("Rows are comparisons. Red dots indicate that genes in this intersection are present in the comparison.", side = 1, line = 3.4, cex = 0.72, col = "#4B5563")
    par(old_par)
  } else {
    plot.new()
    text(0.5, 0.5, paste0("Not enough genes for ", plot_title))
  }
  dev.off()
}

write_upset_outputs(
  sig_sets,
  file.path(out_dir, "UpSet_input_DEG_membership.csv"),
  file.path(fig_dir, "DEG_UpSet_plot.pdf"),
  "DEG intersections"
)

write_upset_outputs(
  up_sets,
  file.path(out_dir, "UpSet_input_up_genes_membership.csv"),
  file.path(fig_dir, "DEG_UpSet_up_genes.pdf"),
  "Up-regulated gene intersections"
)

write_upset_outputs(
  down_sets,
  file.path(out_dir, "UpSet_input_down_genes_membership.csv"),
  file.path(fig_dir, "DEG_UpSet_down_genes.pdf"),
  "Down-regulated gene intersections"
)

message("[INFO] DESeq2 analysis finished")
