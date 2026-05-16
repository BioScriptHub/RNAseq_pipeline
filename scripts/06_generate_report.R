#!/usr/bin/env Rscript

# Generate a self-contained project-level HTML summary report.
# Input: metadata, contrasts, QC summaries, DESeq2 outputs and figures
# Output: 14_report/RNAseq_analysis_report.html

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

h <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
}

safe_read_delim <- function(file, sep = "\t", ...) {
  if (!file.exists(file)) {
    return(NULL)
  }
  read.delim(file, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

safe_read_csv <- function(file, ...) {
  if (!file.exists(file)) {
    return(NULL)
  }
  read.csv(file, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

table_html <- function(df, max_rows = 20) {
  if (is.null(df) || nrow(df) == 0) {
    return("<p class=\"muted\">未找到可展示的数据。</p>")
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (nrow(df) > max_rows) {
    df <- utils::head(df, max_rows)
  }
  header <- paste0("<tr>", paste0("<th>", h(colnames(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(row) {
    paste0("<tr>", paste0("<td>", h(row), "</td>", collapse = ""), "</tr>")
  })
  paste0("<div class=\"table-wrap\"><table>", header, paste(rows, collapse = "\n"), "</table></div>")
}

href_for_report <- function(path) {
  file.path("..", path)
}

link_if_exists <- function(path, label = basename(path), href = href_for_report(path)) {
  if (file.exists(path)) {
    paste0("<a href=\"", h(href), "\" target=\"_blank\">", h(label), "</a>")
  } else {
    paste0("<span class=\"missing\">未找到：", h(path), "</span>")
  }
}

figure_block <- function(path, title, note = "", href = href_for_report(path)) {
  if (!file.exists(path)) {
    return(paste0(
      "<div class=\"figure-card\"><h4>", h(title), "</h4>",
      "<p class=\"missing\">未找到图文件：", h(path), "</p></div>"
    ))
  }
  paste0(
    "<div class=\"figure-card\">",
    "<h4>", h(title), "</h4>",
    if (nzchar(note)) paste0("<p class=\"muted\">", h(note), "</p>") else "",
    "<iframe src=\"", h(href), "\" title=\"", h(title), "\"></iframe>",
    "<p>", link_if_exists(path, "打开 PDF 图", href = href), "</p>",
    "</div>"
  )
}

metadata_file <- get_env("METADATA", "09_metadata/metadata.tsv")
contrasts_file <- get_env("CONTRASTS", "09_metadata/contrasts.tsv")
count_dir <- get_env("COUNT_DIR", "08_counts")
deseq2_dir <- get_env("DESEQ2_DIR", "10_deseq2")
fig_dir <- get_env("FIGURE_DIR", "11_figures")
candidate_dir <- get_env("CANDIDATE_DIR", "13_candidate_genes")
report_dir <- get_env("REPORT_DIR", "14_report")
log_dir <- get_env("LOG_DIR", "logs")

padj_cutoff <- as.numeric(get_env("PADJ_CUTOFF", "0.05"))
lfc_cutoff <- as.numeric(get_env("LFC_CUTOFF", "1"))
min_count <- as.integer(get_env("MIN_COUNT", "10"))
min_samples <- as.integer(get_env("MIN_SAMPLES", "3"))
strandness <- get_env("STRANDNESS", "0")
feature_type <- get_env("FEATURE_TYPE", "exon")
group_attribute <- get_env("GROUP_ATTRIBUTE", "auto")
hisat2_splice_mode <- get_env("HISAT2_SPLICE_MODE", "auto")
hisat2_dta <- get_env("HISAT2_DTA", "1")

dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
report_file <- file.path(report_dir, "RNAseq_analysis_report.html")

metadata <- safe_read_delim(metadata_file)
if (is.null(metadata)) {
  stop_with_help(
    paste0("Metadata file not found: ", metadata_file),
    "Run bash scripts/00_make_metadata.sh or bash scripts/run_all.sh before generating the report."
  )
}

contrasts <- safe_read_delim(contrasts_file)
if (is.null(contrasts)) {
  stop_with_help(
    paste0("Contrasts file not found: ", contrasts_file),
    "Run bash scripts/00_make_metadata.sh or provide 09_metadata/contrasts.tsv."
  )
}

deg_summary_file <- file.path(deseq2_dir, "DEG_summary.csv")
deg_summary <- safe_read_csv(deg_summary_file)
if (is.null(deg_summary)) {
  stop_with_help(
    paste0("DEG summary file not found: ", deg_summary_file),
    "Run Rscript scripts/05_DESeq2.R before generating the report."
  )
}

sample_n <- nrow(metadata)
condition_summary <- as.data.frame(table(metadata$condition), stringsAsFactors = FALSE)
colnames(condition_summary) <- c("condition", "sample_count")
condition_n <- nrow(condition_summary)
contrast_n <- nrow(contrasts)

mapping_file <- file.path(log_dir, "hisat2_alignment_rate_summary.txt")
mapping_lines <- if (file.exists(mapping_file)) readLines(mapping_file, warn = FALSE) else character()
mapping_df <- data.frame(sample_id = character(), alignment_rate_percent = numeric())
if (length(mapping_lines) > 0) {
  rates <- sub(".*:([0-9.]+)% overall alignment rate.*", "\\1", mapping_lines)
  samples <- basename(sub("\\.hisat2\\.log:.*", "", mapping_lines))
  mapping_df <- data.frame(
    sample_id = samples,
    alignment_rate_percent = suppressWarnings(as.numeric(rates)),
    stringsAsFactors = FALSE
  )
}

count_summary_file <- file.path(count_dir, "gene_counts.txt.summary")
count_summary <- safe_read_delim(count_summary_file)
assigned_df <- NULL
if (!is.null(count_summary) && "Status" %in% colnames(count_summary)) {
  status_col <- count_summary$Status
  value_mat <- count_summary[, setdiff(colnames(count_summary), "Status"), drop = FALSE]
  value_mat[] <- lapply(value_mat, function(x) suppressWarnings(as.numeric(x)))
  assigned <- as.numeric(value_mat[status_col == "Assigned", , drop = TRUE])
  total <- colSums(value_mat, na.rm = TRUE)
  sample_names <- basename(colnames(value_mat))
  sample_names <- sub("\\.sorted\\.bam$", "", sample_names)
  assigned_df <- data.frame(
    sample_id = sample_names,
    assigned_reads = assigned,
    total_counted_reads = total,
    assigned_percent = ifelse(total > 0, assigned / total * 100, NA_real_),
    stringsAsFactors = FALSE
  )
}

tested_genes <- if ("tested_genes" %in% colnames(deg_summary)) unique(deg_summary$tested_genes) else NA
max_deg_row <- if (nrow(deg_summary) > 0) deg_summary[which.max(deg_summary$significant_DEGs), , drop = FALSE] else NULL

qc_text <- if (nrow(mapping_df) > 0) {
  paste0(
    "HISAT2 总比对率范围为 ",
    fmt_num(min(mapping_df$alignment_rate_percent, na.rm = TRUE)), "% 到 ",
    fmt_num(max(mapping_df$alignment_rate_percent, na.rm = TRUE)), "%，平均值为 ",
    fmt_num(mean(mapping_df$alignment_rate_percent, na.rm = TRUE)), "%。"
  )
} else {
  "未找到 HISAT2 比对率汇总文件，报告无法自动评价比对率。"
}

count_text <- if (!is.null(assigned_df) && nrow(assigned_df) > 0) {
  paste0(
    "featureCounts 的 Assigned 比例范围为 ",
    fmt_num(min(assigned_df$assigned_percent, na.rm = TRUE)), "% 到 ",
    fmt_num(max(assigned_df$assigned_percent, na.rm = TRUE)), "%，平均值为 ",
    fmt_num(mean(assigned_df$assigned_percent, na.rm = TRUE)), "%。"
  )
} else {
  "未找到 featureCounts summary，报告无法自动计算 Assigned 比例。"
}

deg_text <- if (!is.null(max_deg_row) && nrow(max_deg_row) == 1) {
  paste0(
    "本次分析共完成 ", contrast_n, " 个差异比较。差异基因数量最多的比较是 ",
    max_deg_row$comparison, "，共检测到 ", max_deg_row$significant_DEGs,
    " 个显著差异基因，其中上调 ", max_deg_row$up_genes,
    " 个，下调 ", max_deg_row$down_genes, " 个。"
  )
} else {
  "未找到差异表达汇总信息。"
}

method_text <- paste0(
  "本流程以双端 RNA-seq 原始测序文件为输入。首先使用 fastp 进行接头识别和低质量 reads 过滤，",
  "随后使用 FastQC 和 MultiQC 汇总测序质量。过滤后的 reads 使用 HISAT2 比对到参考基因组；",
  "当前剪接模式为 ", hisat2_splice_mode, "，HISAT2_DTA=", hisat2_dta, "。比对结果经 samtools sort 和 index 生成排序 BAM 文件。",
  "基因水平计数使用 featureCounts 完成，参数设置为 FEATURE_TYPE=", feature_type,
  "，GROUP_ATTRIBUTE=", group_attribute, "，STRANDNESS=", strandness, "。差异表达分析使用 DESeq2，设计公式为 ~ condition。",
  "低表达过滤条件为至少 ", min_samples, " 个样本 count >= ", min_count, "。显著差异基因阈值为 padj < ",
  padj_cutoff, " 且 |log2FoldChange| >= ", lfc_cutoff, "。"
)

deg_summary_display <- deg_summary
numeric_cols <- intersect(c("tested_genes", "significant_DEGs", "up_genes", "down_genes"), colnames(deg_summary_display))
deg_summary_display[numeric_cols] <- lapply(deg_summary_display[numeric_cols], as.character)

contrast_sections <- character()
for (i in seq_len(nrow(deg_summary))) {
  row <- deg_summary[i, , drop = FALSE]
  comparison <- row$comparison
  comparison_safe <- gsub("[^A-Za-z0-9_.-]+", "_", comparison)
  sig_file <- file.path(deseq2_dir, "comparisons", paste0(comparison_safe, ".significant_DEGs.csv"))
  top_file <- file.path(candidate_dir, paste0(comparison_safe, ".top_candidate_genes_ranked.csv"))
  top_df <- safe_read_csv(top_file)
  if (!is.null(top_df)) {
    keep_cols <- intersect(
      c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "candidate_score"),
      colnames(top_df)
    )
    top_df <- top_df[, keep_cols, drop = FALSE]
  }
  interpretation <- paste0(
    comparison, "：显著差异基因 ", row$significant_DEGs,
    " 个，上调 ", row$up_genes, " 个，下调 ", row$down_genes,
    " 个。若该比较对应明确处理因素，优先关注排名靠前且表达变化方向与生物学假设一致的基因。"
  )
  contrast_sections <- c(
    contrast_sections,
    paste0(
      "<section class=\"subsection\"><h3>", h(comparison), "</h3>",
      "<p>", h(interpretation), "</p>",
      "<p>结果表：", link_if_exists(sig_file, "显著差异基因表"), "；",
      link_if_exists(top_file, "候选基因排序表"), "</p>",
      "<div class=\"figure-grid\">",
      figure_block(file.path(fig_dir, paste0(comparison_safe, ".volcano.pdf")), paste0(comparison, " 火山图"), "横轴为 log2FoldChange，纵轴为 -log10 adjusted P-value。"),
      figure_block(file.path(fig_dir, paste0(comparison_safe, ".MA_plot.pdf")), paste0(comparison, " MA 图"), "用于检查表达量与差异倍数的关系。"),
      "</div>",
      "<h4>Top candidate genes</h4>",
      table_html(top_df, max_rows = 10),
      "</section>"
    )
  )
}

preflight <- safe_read_delim(file.path(log_dir, "preflight_check.tsv"))
reference_stats <- safe_read_delim(file.path(log_dir, "reference_genome.seqkit_stats.txt"), sep = "")
if (!is.null(reference_stats)) {
  reference_stats <- reference_stats[, colSums(!is.na(reference_stats)) > 0, drop = FALSE]
}

html <- paste0(
  "<!doctype html>\n<html lang=\"zh-CN\">\n<head>\n<meta charset=\"utf-8\">\n",
  "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
  "<title>RNA-seq 分析报告</title>\n",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI','Microsoft YaHei',Arial,sans-serif;margin:0;background:#f6f7f9;color:#1f2933;line-height:1.65;}",
  "main{max-width:1180px;margin:0 auto;padding:28px 22px 60px;}",
  "header{background:#111827;color:#fff;padding:34px 40px;border-radius:10px;margin-bottom:22px;}",
  "h1{margin:0 0 8px;font-size:30px;}h2{font-size:22px;margin:30px 0 12px;border-left:5px solid #334155;padding-left:12px;}h3{font-size:18px;margin-top:24px;}h4{margin:16px 0 8px;}",
  ".card{background:#fff;border:1px solid #e5e7eb;border-radius:10px;padding:20px;margin:16px 0;box-shadow:0 1px 2px rgba(15,23,42,.04);}",
  ".summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin-top:16px;}",
  ".metric{background:#f8fafc;border:1px solid #e2e8f0;border-radius:8px;padding:14px;}.metric b{display:block;font-size:24px;color:#0f172a;}.metric span{color:#64748b;font-size:13px;}",
  ".table-wrap{overflow-x:auto;}table{border-collapse:collapse;width:100%;font-size:13px;}th,td{border:1px solid #e5e7eb;padding:7px 9px;text-align:left;}th{background:#f1f5f9;color:#0f172a;}tr:nth-child(even){background:#fafafa;}",
  ".muted{color:#64748b;}.missing{color:#b91c1c;font-weight:600;}.note{background:#fff7ed;border:1px solid #fed7aa;border-radius:8px;padding:12px 14px;}",
  ".figure-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:16px;}.figure-card{background:#fff;border:1px solid #e5e7eb;border-radius:10px;padding:14px;margin:12px 0;}",
  "iframe{width:100%;height:420px;border:1px solid #d1d5db;border-radius:6px;background:#fff;}a{color:#1d4ed8;text-decoration:none;}a:hover{text-decoration:underline;}",
  "code{background:#eef2f7;padding:1px 4px;border-radius:4px;}.footer{color:#64748b;font-size:13px;margin-top:28px;}",
  "</style>\n</head>\n<body>\n<main>\n",
  "<header><h1>RNA-seq 分析报告</h1><p>自动生成时间：", h(format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "</p></header>",
  "<section class=\"card\"><h2>1. 项目概览</h2>",
  "<div class=\"summary\">",
  "<div class=\"metric\"><b>", sample_n, "</b><span>样本数</span></div>",
  "<div class=\"metric\"><b>", condition_n, "</b><span>分组数</span></div>",
  "<div class=\"metric\"><b>", contrast_n, "</b><span>差异比较数</span></div>",
  "<div class=\"metric\"><b>", h(paste(tested_genes, collapse = ", ")), "</b><span>进入差异分析的基因数</span></div>",
  "</div>",
  "<p>", h(qc_text), "</p><p>", h(count_text), "</p><p>", h(deg_text), "</p>",
  "<p class=\"note\">本报告用于整理常规 RNA-seq 分析结果。差异表达表示统计关联，不等同于因果证明；候选基因仍需结合注释、实验设计和后续验证解释。</p>",
  "</section>",
  "<section class=\"card\"><h2>2. 材料与方法</h2>",
  "<h3>2.1 样本与比较设计</h3>",
  "<p>本项目包含 ", sample_n, " 个样本，分属于 ", condition_n, " 个处理组。每个处理组的样本数如下。</p>",
  table_html(condition_summary, max_rows = 50),
  "<h3>2.2 分析流程</h3><p>", h(method_text), "</p>",
  "<h3>2.3 输入文件与软件检查</h3>",
  "<p>预检查结果用于确认样本表、比较表、原始数据、参考文件和软件命令是否满足流程要求。</p>",
  table_html(preflight, max_rows = 30),
  "<h3>2.4 参考基因组概况</h3>",
  table_html(reference_stats, max_rows = 10),
  "</section>",
  "<section class=\"card\"><h2>3. 质控、比对与计数结果</h2>",
  "<p>原始 reads 经过 fastp 过滤后进入 HISAT2 比对。建议先查看 MultiQC 报告，再解释差异表达结果。</p>",
  "<p>FASTQ 质控汇总：", link_if_exists(file.path("03_multiqc", "multiqc_fastq_qc.html"), "multiqc_fastq_qc.html"), "</p>",
  "<p>BAM 质控汇总：", link_if_exists(file.path("03_multiqc", "multiqc_bam_qc.html"), "multiqc_bam_qc.html"), "</p>",
  "<h3>3.1 HISAT2 比对率</h3>",
  table_html(mapping_df, max_rows = 100),
  "<h3>3.2 featureCounts Assigned 比例</h3>",
  table_html(assigned_df, max_rows = 100),
  "</section>",
  "<section class=\"card\"><h2>4. 样本关系</h2>",
  "<p>PCA 和样本相关性用于检查生物学重复是否聚集、是否存在离群样本或明显批次问题。</p>",
  "<div class=\"figure-grid\">",
  figure_block(file.path(fig_dir, "PCA_condition.pdf"), "PCA 图", "同组样本应尽量聚集，明显离群样本需要回查质控和样本来源。"),
  figure_block(file.path(fig_dir, "sample_correlation_heatmap.pdf"), "样本相关性热图", "组内相关性通常应高于组间相关性。"),
  "</div></section>",
  "<section class=\"card\"><h2>5. 差异表达结果</h2>",
  "<p>", h(deg_text), "</p>",
  table_html(deg_summary_display, max_rows = 100),
  "<p>差异结果目录：", link_if_exists(file.path(deseq2_dir, "comparisons"), "10_deseq2/comparisons/"), "</p>",
  paste(contrast_sections, collapse = "\n"),
  "</section>",
  "<section class=\"card\"><h2>6. 差异基因交集</h2>",
  "<p>UpSet 图用于查看多个差异比较之间共享和特异的差异基因集合。</p>",
  figure_block(file.path(fig_dir, "DEG_UpSet_plot.pdf"), "全部 DEG UpSet 图", "交集较大的比较可能共享相似转录响应，特异集合更适合筛选处理特异候选基因。"),
  figure_block(file.path(fig_dir, "DEG_UpSet_up_genes.pdf"), "上调基因 UpSet 图", "只统计各比较中的上调差异基因。"),
  figure_block(file.path(fig_dir, "DEG_UpSet_down_genes.pdf"), "下调基因 UpSet 图", "只统计各比较中的下调差异基因。"),
  "<p>交集输入表：", link_if_exists(file.path(deseq2_dir, "UpSet_input_DEG_membership.csv"), "全部 DEG membership"), "；",
  link_if_exists(file.path(deseq2_dir, "UpSet_input_up_genes_membership.csv"), "上调基因 membership"), "；",
  link_if_exists(file.path(deseq2_dir, "UpSet_input_down_genes_membership.csv"), "下调基因 membership"), "</p>",
  "</section>",
  "<section class=\"card\"><h2>7. 结果解读要点</h2>",
  "<ul>",
  "<li>先确认 FASTQ 质量、比对率和 Assigned 比例，再解释差异基因。</li>",
  "<li>PCA 或相关性热图中若存在明显离群样本，应优先排查样本质量、分组和批次。</li>",
  "<li>显著差异基因数量反映转录扰动强度，但不能单独证明机制。</li>",
  "<li>候选基因应结合功能注释、表达方向、处理特异性和已有生物学知识筛选。</li>",
  "<li>RNA-seq 能提出高质量假说，最终因果关系仍需遗传学、分子实验或表型实验验证。</li>",
  "</ul>",
  "</section>",
  "<p class=\"footer\">报告生成脚本：scripts/06_generate_report.R；报告文件：", h(report_file), "</p>",
  "</main>\n</body>\n</html>\n"
)

writeLines(html, report_file, useBytes = TRUE)
message("[INFO] HTML report written: ", report_file)
