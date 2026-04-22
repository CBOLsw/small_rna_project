#!/usr/bin/env Rscript
# DESeq2差异表达分析脚本
#
# 功能：
# 1. 读取基因计数矩阵和样本信息
# 2. 使用DESeq2进行差异表达分析
# 3. 筛选差异表达基因（DEGs）
# 4. 生成分析结果和统计报告
# 5. 创建可视化图表（火山图、MA图、热图）
#
# 使用方法：
#    Rscript deseq2_analysis.R --counts <计数矩阵文件> --metadata <样本信息文件> --output <输出目录>
#
# 依赖：
#    DESeq2, ggplot2, pheatmap, dplyr, tidyverse
#
# 安装依赖：
#    if (!require("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#    BiocManager::install("DESeq2")
#    install.packages(c("ggplot2", "pheatmap", "dplyr", "tidyverse"))

# 加载必要的库
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(tidyverse)
  library(optparse)
  library(RColorBrewer)
})

# 配置日志
log_message <- function(msg, level = "INFO") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] [%s] %s\n", timestamp, level, msg))
}

# 解析命令行参数
parse_arguments <- function() {
  option_list <- list(
    make_option(c("-c", "--counts"),
                type = "character",
                default = NULL,
                help = "基因计数矩阵文件（CSV/TSV格式）",
                metavar = "FILE"),

    make_option(c("-m", "--metadata"),
                type = "character",
                default = NULL,
                help = "样本信息文件（CSV格式，包含sample和group列）",
                metavar = "FILE"),

    make_option(c("-o", "--output"),
                type = "character",
                default = "results/expression/deseq2",
                help = "输出目录 [默认: %default]",
                metavar = "DIR"),

    make_option(c("--group-col"),
                type = "character",
                default = "group",
                help = "分组列名 [默认: %default]"),

    make_option(c("--control-group"),
                type = "character",
                default = NULL,
                help = "对照组名称（如不指定，使用第一个组）"),

    make_option(c("--treatment-group"),
                type = "character",
                default = NULL,
                help = "处理组名称（如不指定，使用第二个组）"),

    make_option(c("--fc-threshold"),
                type = "numeric",
                default = 1.5,
                help = "倍数变化阈值（log2 fold-change） [默认: %default]",
                metavar = "NUMBER"),

    make_option(c("--padj-threshold"),
                type = "numeric",
                default = 0.05,
                help = "调整后p值阈值 [默认: %default]",
                metavar = "NUMBER"),

    make_option(c("--min-mean-count"),
                type = "numeric",
                default = 10,
                help = "最小平均计数阈值 [默认: %default]",
                metavar = "NUMBER"),

    make_option(c("--filter-low-counts"),
                type = "logical",
                default = TRUE,
                help = "是否过滤低表达基因 [默认: %default]"),

    make_option(c("--alpha"),
                type = "numeric",
                default = 0.05,
                help = "显著性水平 [默认: %default]",
                metavar = "NUMBER"),

    make_option(c("--cooks-cutoff"),
                type = "logical",
                default = TRUE,
                help = "是否使用Cook's距离过滤 [默认: %default]"),

    make_option(c("--independent-filtering"),
                type = "logical",
                default = TRUE,
                help = "是否使用独立过滤 [默认: %default]"),

    make_option(c("--generate-plots"),
                type = "logical",
                default = TRUE,
                help = "是否生成可视化图表 [默认: %default]")
  )

  parser <- OptionParser(option_list = option_list)
  args <- parse_args(parser)

  # 检查必需参数
  if (is.null(args$counts)) {
    print_help(parser)
    stop("--counts 参数是必需的")
  }

  if (is.null(args$metadata)) {
    print_help(parser)
    stop("--metadata 参数是必需的")
  }

  return(args)
}

# 读取计数矩阵
read_count_matrix <- function(count_file) {
  log_message(sprintf("读取计数矩阵: %s", count_file))

  # 根据文件扩展名确定格式
  file_ext <- tools::file_ext(count_file)

  # featureCounts输出通常是TSV格式（即使扩展名是.csv）
  # 先尝试用制表符读取，跳过注释行
  counts <- tryCatch({
    # 读取所有行
    all_lines <- readLines(count_file)

    # 找到数据行（不是#开头的注释）
    data_lines <- all_lines[!grepl("^#", all_lines)]

    # 用制表符分割
    data <- read.delim(textConnection(paste(data_lines, collapse = "\n")),
                       header = TRUE, check.names = FALSE)

    # 第一列是基因ID
    rownames(data) <- data[, 1]
    data <- data[, -1, drop = FALSE]

    data
  }, error = function(e) {
    # 如果失败，尝试其他方法
    if (file_ext %in% c("csv", "CSV")) {
      read.csv(count_file, row.names = 1, check.names = FALSE)
    } else {
      read.delim(count_file, row.names = 1, check.names = FALSE)
    }
  })

  # 从列名中提取样本名（featureCounts输出的是完整路径）
  # 例如: "results/alignment/GAO_1.sorted.bam" -> "GAO_1"
  old_names <- colnames(counts)
  log_message(sprintf("原始列名 (%d): %s", length(old_names), paste(old_names, collapse = ", ")))

  new_names <- gsub(".*/", "", old_names)  # 去掉路径
  new_names <- gsub("\\.sorted\\.bam$", "", new_names)  # 去掉后缀
  colnames(counts) <- new_names

  log_message(sprintf("处理后列名 (%d): %s", length(new_names), paste(new_names, collapse = ", ")))

  # 过滤掉非数值列（如 Chr, Start, End, Strand, Length 等元数据列）
  metadata_cols <- c("Chr", "Start", "End", "Strand", "Length", "Chr.1", "Start.1", "End.1")
  counts <- counts[, !colnames(counts) %in% metadata_cols, drop = FALSE]

  # 确保所有列都是数值型
  for (col in colnames(counts)) {
    counts[[col]] <- as.numeric(counts[[col]])
  }

  log_message(sprintf("计数矩阵维度: %d 基因 x %d 样本",
                      nrow(counts), ncol(counts)))
  log_message(sprintf("总计数: %d", sum(counts, na.rm = TRUE)))

  return(counts)
}

# 读取样本信息
read_metadata <- function(metadata_file, group_col = "group") {
  log_message(sprintf("读取样本信息: %s", metadata_file))

  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

  # 检查必要的列
  required_cols <- c("sample", group_col)
  missing_cols <- setdiff(required_cols, colnames(metadata))

  if (length(missing_cols) > 0) {
    stop(sprintf("样本信息文件缺少必要的列: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  # 确保样本名为字符
  metadata$sample <- as.character(metadata$sample)

  # 设置分组为因子
  metadata[[group_col]] <- factor(metadata[[group_col]])

  log_message(sprintf("样本数: %d", nrow(metadata)))
  log_message(sprintf("分组: %s",
                      paste(levels(metadata[[group_col]]), collapse = " vs ")))

  return(metadata)
}

# 检查计数矩阵和样本信息的一致性
check_consistency <- function(counts, metadata, group_col = "group") {
  log_message("检查计数矩阵和样本信息的一致性")

  # 获取计数矩阵中的样本名
  count_samples <- colnames(counts)
  log_message(sprintf("计数矩阵中的样本名 (%d): %s",
                      length(count_samples), paste(count_samples, collapse = ", ")))

  # 获取样本信息中的样本名
  metadata_samples <- metadata$sample
  log_message(sprintf("样本信息中的样本名 (%d): %s",
                      length(metadata_samples), paste(metadata_samples, collapse = ", ")))

  # 检查样本名是否匹配
  missing_in_metadata <- setdiff(count_samples, metadata_samples)
  missing_in_counts <- setdiff(metadata_samples, count_samples)

  if (length(missing_in_metadata) > 0) {
    warning(sprintf("以下样本在计数矩阵中存在但样本信息中缺失: %s",
                    paste(missing_in_metadata, collapse = ", ")))
  }

  if (length(missing_in_counts) > 0) {
    warning(sprintf("以下样本在样本信息中存在但计数矩阵中缺失: %s",
                    paste(missing_in_counts, collapse = ", ")))
  }

  # 使用交集样本
  common_samples <- intersect(count_samples, metadata_samples)
  log_message(sprintf("共同样本数: %d", length(common_samples)))

  if (length(common_samples) < 2) {
    stop("没有足够的共同样本进行分析")
  }

  log_message(sprintf("使用 %d 个共同样本", length(common_samples)))

  # 过滤计数矩阵和样本信息
  counts_filtered <- counts[, common_samples, drop = FALSE]
  metadata_filtered <- metadata[metadata$sample %in% common_samples, , drop = FALSE]

  # 确保顺序一致
  counts_filtered <- counts_filtered[, metadata_filtered$sample]

  return(list(counts = counts_filtered, metadata = metadata_filtered))
}

# 运行DESeq2分析
run_deseq2 <- function(counts, metadata, group_col = "group",
                       control_group = NULL, treatment_group = NULL,
                       filter_low_counts = TRUE,
                       cooks_cutoff = TRUE,
                       independent_filtering = TRUE,
                       alpha = 0.05) {
  log_message("运行DESeq2差异表达分析")

  # 创建DESeqDataSet
  log_message("创建DESeqDataSet对象")

  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = as.formula(paste("~", group_col))
  )

  # 设置参考水平（对照组）
  if (!is.null(control_group)) {
    dds[[group_col]] <- relevel(dds[[group_col]], ref = control_group)
    log_message(sprintf("设置对照组: %s", control_group))
  }

  # 过滤低表达基因（可选）
  if (filter_low_counts) {
    log_message("过滤低表达基因")
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep, ]
    log_message(sprintf("保留 %d 个基因（原始 %d）", nrow(dds), nrow(counts)))
  }

  # 运行DESeq2
  log_message("运行DESeq2分析")
  dds <- DESeq(dds, quiet = TRUE)

  # 获取对比组信息
  groups <- levels(dds[[group_col]])

  if (is.null(control_group)) {
    control_group <- groups[1]
  }

  if (is.null(treatment_group)) {
    treatment_group <- groups[2]
  }

  if (!control_group %in% groups) {
    stop(sprintf("对照组 '%s' 不在分组水平中: %s",
                 control_group, paste(groups, collapse = ", ")))
  }

  if (!treatment_group %in% groups) {
    stop(sprintf("处理组 '%s' 不在分组水平中: %s",
                 treatment_group, paste(groups, collapse = ", ")))
  }

  contrast <- c(group_col, treatment_group, control_group)

  log_message(sprintf("对比: %s vs %s", treatment_group, control_group))

  # 获取结果
  res <- results(dds,
                 contrast = contrast,
                 alpha = alpha,
                 cooksCutoff = cooks_cutoff,
                 independentFiltering = independent_filtering)

  # 添加基因名
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)

  # 重新排序列
  res_df <- res_df %>%
    select(gene_id, everything())

  # 添加表达水平信息
  res_df$baseMean_control <- rowMeans(counts(dds)[, dds[[group_col]] == control_group])
  res_df$baseMean_treatment <- rowMeans(counts(dds)[, dds[[group_col]] == treatment_group])

  # 计算原始倍数变化（非log2）
  res_df$fold_change <- 2^res_df$log2FoldChange

  log_message(sprintf("分析完成，总基因数: %d", nrow(res_df)))

  return(list(results = res_df,
              dds = dds,
              contrast = contrast,
              control_group = control_group,
              treatment_group = treatment_group))
}

# 筛选差异表达基因
filter_degs <- function(results_df, fc_threshold = 1.5, padj_threshold = 0.05) {
  log_message(sprintf("筛选差异表达基因 (|log2FC| >= %.2f, padj < %.3f)",
                      log2(fc_threshold), padj_threshold))

  # 计算绝对值
  results_df$abs_log2fc <- abs(results_df$log2FoldChange)

  # 筛选DEGs
  degs <- results_df %>%
    filter(!is.na(padj)) %>%
    filter(padj < padj_threshold) %>%
    filter(abs_log2fc >= log2(fc_threshold))

  # 分上下调
  degs_up <- degs %>%
    filter(log2FoldChange > 0) %>%
    arrange(desc(abs_log2fc))

  degs_down <- degs %>%
    filter(log2FoldChange < 0) %>%
    arrange(desc(abs_log2fc))

  log_message(sprintf("差异表达基因总数: %d", nrow(degs)))
  log_message(sprintf("上调基因: %d", nrow(degs_up)))
  log_message(sprintf("下调基因: %d", nrow(degs_down)))

  return(list(all = degs,
              up = degs_up,
              down = degs_down))
}

# 生成分析报告
generate_report <- function(deseq2_results, degs, output_dir,
                            fc_threshold, padj_threshold,
                            control_group, treatment_group) {
  log_message("生成分析报告")

  results_df <- deseq2_results$results

  # 创建报告目录
  report_dir <- file.path(output_dir, "reports")
  dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. 总体统计报告
  overall_stats <- data.frame(
    项目 = c("总基因数", "差异表达基因数", "上调基因数", "下调基因数",
             "对照组", "处理组", "倍数变化阈值", "调整后p值阈值"),
    数值 = c(nrow(results_df), nrow(degs$all), nrow(degs$up), nrow(degs$down),
            control_group, treatment_group, fc_threshold, padj_threshold)
  )

  overall_file <- file.path(report_dir, "deseq2_overall_stats.csv")
  write.csv(overall_stats, overall_file, row.names = FALSE)

  # 2. 完整结果文件
  full_results_file <- file.path(output_dir, "deseq2_full_results.csv")
  write.csv(results_df, full_results_file, row.names = FALSE)

  # 3. 差异表达基因列表
  degs_file <- file.path(output_dir, "deseq2_significant_genes.csv")
  write.csv(degs$all, degs_file, row.names = FALSE)

  # 4. 上调基因列表
  up_degs_file <- file.path(output_dir, "deseq2_upregulated_genes.csv")
  write.csv(degs$up, up_degs_file, row.names = FALSE)

  # 5. 下调基因列表
  down_degs_file <- file.path(output_dir, "deseq2_downregulated_genes.csv")
  write.csv(degs$down, down_degs_file, row.names = FALSE)

  # 6. 生成文本报告
  text_report <- file.path(report_dir, "deseq2_analysis_report.txt")

  sink(text_report)
  cat("=== DESeq2差异表达分析报告 ===\n\n")
  cat(sprintf("生成时间: %s\n\n", Sys.time()))

  cat("1. 分析概览\n")
  cat(sprintf("   总基因数: %d\n", nrow(results_df)))
  cat(sprintf("   对照组: %s (n=%d)\n",
              control_group,
              sum(!is.na(results_df$baseMean_control))))
  cat(sprintf("   处理组: %s (n=%d)\n",
              treatment_group,
              sum(!is.na(results_df$baseMean_treatment))))
  cat(sprintf("   倍数变化阈值: %.2f (|log2FC| >= %.2f)\n",
              fc_threshold, log2(fc_threshold)))
  cat(sprintf("   显著性阈值: padj < %.3f\n\n", padj_threshold))

  cat("2. 差异表达基因统计\n")
  cat(sprintf("   差异表达基因总数: %d\n", nrow(degs$all)))
  cat(sprintf("   上调基因数: %d (%.1f%%)\n",
              nrow(degs$up), nrow(degs$up)/nrow(degs$all)*100))
  cat(sprintf("   下调基因数: %d (%.1f%%)\n",
              nrow(degs$down), nrow(degs$down)/nrow(degs$all)*100))
  cat(sprintf("   非差异表达基因: %d\n\n", nrow(results_df) - nrow(degs$all)))

  cat("3. 表达水平统计\n")
  cat(sprintf("   平均表达水平 (所有基因): %.2f\n", mean(results_df$baseMean, na.rm = TRUE)))
  cat(sprintf("   中位表达水平 (所有基因): %.2f\n", median(results_df$baseMean, na.rm = TRUE)))
  cat(sprintf("   最大表达水平: %.2f\n\n", max(results_df$baseMean, na.rm = TRUE)))

  cat("4. 倍数变化分布\n")
  cat(sprintf("   最大上调倍数: %.2f (log2FC=%.2f)\n",
              max(degs$up$fold_change, na.rm = TRUE),
              max(degs$up$log2FoldChange, na.rm = TRUE)))
  cat(sprintf("   最大下调倍数: %.2f (log2FC=%.2f)\n",
              1/min(degs$down$fold_change, na.rm = TRUE),
              min(degs$down$log2FoldChange, na.rm = TRUE)))
  cat(sprintf("   平均倍数变化: %.2f\n\n", mean(abs(results_df$fold_change), na.rm = TRUE)))

  cat("5. 质量控制指标\n")
  cat(sprintf("   缺失值 (padj=NA): %d\n", sum(is.na(results_df$padj))))
  cat(sprintf("   低表达过滤: %s\n", ifelse(any(grepl("rowSums", deseq2_results$dds@metadata$filterThreshold)), "是", "否")))
  cat(sprintf("   Cook's距离过滤: %s\n", ifelse(any(grepl("cooksCutoff", deseq2_results$dds@metadata$resultsNames)), "是", "否")))

  cat("\n6. 文件清单\n")
  cat("   deseq2_full_results.csv - 完整分析结果\n")
  cat("   deseq2_significant_genes.csv - 差异表达基因列表\n")
  cat("   deseq2_upregulated_genes.csv - 上调基因列表\n")
  cat("   deseq2_downregulated_genes.csv - 下调基因列表\n")

  sink()

  log_message(sprintf("报告已生成: %s", text_report))
}

# 生成可视化图表
generate_plots <- function(deseq2_results, degs, output_dir,
                          control_group, treatment_group,
                          fc_threshold, padj_threshold) {
  log_message("生成可视化图表")

  results_df <- deseq2_results$results
  dds <- deseq2_results$dds

  # 创建图表目录
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. 火山图
  log_message("生成火山图")

  volcano_data <- results_df %>%
    mutate(significant = ifelse(
      !is.na(padj) & padj < padj_threshold & abs(log2FoldChange) >= log2(fc_threshold),
      ifelse(log2FoldChange > 0, "上调", "下调"), "非差异"))

  volcano_plot <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("上调" = "red", "下调" = "blue", "非差异" = "grey")) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_vline(xintercept = c(-log2(fc_threshold), log2(fc_threshold)),
               linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = sprintf("火山图: %s vs %s", treatment_group, control_group),
         x = expression(log[2]~"倍数变化"),
         y = expression(-log[10]~"调整后p值"),
         color = "差异表达") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))

  volcano_file <- file.path(plot_dir, "volcano_plot.png")
  ggsave(volcano_file, volcano_plot, width = 10, height = 8, dpi = 300)

  # 2. MA图
  log_message("生成MA图")

  ma_plot <- ggplot(results_df, aes(x = log10(baseMean), y = log2FoldChange)) +
    geom_point(aes(color = ifelse(!is.na(padj) & padj < padj_threshold &
                                    abs(log2FoldChange) >= log2(fc_threshold),
                                  "差异", "非差异")),
               alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("差异" = "red", "非差异" = "grey")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    geom_hline(yintercept = c(-log2(fc_threshold), log2(fc_threshold)),
               linetype = "dashed", color = "blue", alpha = 0.5) +
    labs(title = sprintf("MA图: %s vs %s", treatment_group, control_group),
         x = expression(log[10]~"平均表达"),
         y = expression(log[2]~"倍数变化"),
         color = "差异表达") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))

  ma_file <- file.path(plot_dir, "ma_plot.png")
  ggsave(ma_file, ma_plot, width = 10, height = 8, dpi = 300)

  # 3. 样本间相关性热图
  log_message("生成样本间相关性热图")

  # 计算标准化计数
  vsd <- vst(dds, blind = FALSE)
  sample_cor <- cor(assay(vsd))

  # 准备注释数据
  annotation_df <- as.data.frame(colData(vsd)[, deseq2_results$contrast[1], drop = FALSE])
  colnames(annotation_df) <- "Group"

  # 创建热图
  heatmap_file <- file.path(plot_dir, "sample_correlation_heatmap.png")

  png(heatmap_file, width = 800, height = 800)
  pheatmap(sample_cor,
           annotation_col = annotation_df,
           annotation_row = annotation_df,
           main = "样本间相关性热图",
           color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           fontsize = 10)
  dev.off()

  # 4. 差异表达基因热图
  log_message("生成差异表达基因热图")

  if (nrow(degs$all) > 0) {
    # 选择top差异表达基因
    top_n <- min(50, nrow(degs$all))
    top_genes <- degs$all %>%
      arrange(padj) %>%
      head(top_n) %>%
      pull(gene_id)

    # 获取这些基因的表达数据
    top_expr <- assay(vsd)[top_genes, , drop = FALSE]

    # 标准化行（基因）
    top_expr_scaled <- t(scale(t(top_expr)))

    # 创建热图
    deg_heatmap_file <- file.path(plot_dir, "deg_expression_heatmap.png")

    png(deg_heatmap_file, width = 800, height = 800)
    pheatmap(top_expr_scaled,
             annotation_col = annotation_df,
             main = sprintf("Top %d 差异表达基因热图", top_n),
             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100),
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             show_rownames = ifelse(top_n <= 30, TRUE, FALSE),
             fontsize = 10)
    dev.off()
  }

  log_message(sprintf("图表已保存到: %s", plot_dir))
}

# 主函数
main <- function() {
  # 解析参数
  args <- parse_arguments()

  # 创建输出目录
  dir.create(args$output, showWarnings = FALSE, recursive = TRUE)

  log_message(sprintf("输出目录: %s", args$output))

  # 1. 读取数据
  counts <- read_count_matrix(args$counts)
  metadata <- read_metadata(args$metadata, args$group_col)

  # 2. 检查一致性
  checked <- check_consistency(counts, metadata, args$group_col)
  counts_filtered <- checked$counts
  metadata_filtered <- checked$metadata

  # 3. 运行DESeq2分析
  deseq2_results <- run_deseq2(
    counts = counts_filtered,
    metadata = metadata_filtered,
    group_col = args$group_col,
    control_group = args$control_group,
    treatment_group = args$treatment_group,
    filter_low_counts = args$filter_low_counts,
    cooks_cutoff = args$cooks_cutoff,
    independent_filtering = args$independent_filtering,
    alpha = args$alpha
  )

  # 4. 筛选差异表达基因
  degs <- filter_degs(deseq2_results$results,
                      fc_threshold = args$fc_threshold,
                      padj_threshold = args$padj_threshold)

  # 5. 生成报告
  generate_report(deseq2_results,
                  degs,
                  args$output,
                  fc_threshold = args$fc_threshold,
                  padj_threshold = args$padj_threshold,
                  control_group = deseq2_results$control_group,
                  treatment_group = deseq2_results$treatment_group)

  # 6. 生成可视化图表
  if (args$generate_plots) {
    generate_plots(deseq2_results,
                   degs,
                   args$output,
                   control_group = deseq2_results$control_group,
                   treatment_group = deseq2_results$treatment_group,
                   fc_threshold = args$fc_threshold,
                   padj_threshold = args$padj_threshold)
  }

  log_message("DESeq2差异表达分析完成")
}

# 运行主函数
if (interactive()) {
  # 交互模式测试
  log_message("交互模式")
} else {
  # 命令行模式
  tryCatch({
    main()
  }, error = function(e) {
    log_message(sprintf("分析过程中出错: %s", e$message), level = "ERROR")
    quit(status = 1)
  })
}