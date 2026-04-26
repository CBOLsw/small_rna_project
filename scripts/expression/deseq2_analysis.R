suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(dplyr)
  library(optparse)
  library(RColorBrewer)
})

log_message <- function(msg, level="INFO"){
  cat(sprintf("[%s] [%s] %s\n", format(Sys.time(), "%F %T"), level, msg))
}

# ---------------- 参数 ----------------
parse_arguments <- function(){
  option_list <- list(
    make_option(c("-c","--counts"), type="character"),
    make_option(c("-m","--metadata"), type="character"),
    make_option(c("-o","--output"), type="character", default="results/differential_expression"),
    make_option(c("--group-col"), type="character", default="group"),
    make_option(c("--control-group"), type="character"),
    make_option(c("--treatment-group"), type="character"),
    make_option(c("--lfc-threshold"), type="numeric", default=1.0),
    make_option(c("--padj-threshold"), type="numeric", default=0.05),
    make_option(c("--min-count"), type="numeric", default=10),
    make_option(c("--min-samples"), type="numeric", default=2),
    make_option(c("--generate-plots"), type="logical", default=TRUE)
  )
  args <- parse_args(OptionParser(option_list=option_list))

  if (is.null(args$counts) || is.null(args$metadata)) {
    stop("必须提供 counts 和 metadata")
  }

  return(args)
}

# ---------------- 读取 ----------------
read_counts <- function(file){

  log_message("读取counts")

  lines <- readLines(file)
  lines <- lines[!grepl("^\\s*#", lines)]

  df <- read.delim(text=paste(lines, collapse="\n"), check.names=FALSE)

  rownames(df) <- df[,1]
  df <- df[,-1]

  # 删除featureCounts统计行
  df <- df[!grepl("^__", rownames(df)), ]

  # 删除注释列
  metadata_cols <- c("Chr","Start","End","Strand","Length")
  df <- df[, !(colnames(df) %in% metadata_cols), drop=FALSE]

  # 清理列名
  colnames(df) <- sub("\\.sorted\\.bam$", "", basename(colnames(df)))

  # 强制数值化
  df[] <- lapply(df, function(x){
    x_num <- suppressWarnings(as.numeric(x))
    if (any(is.na(x_num))) {
      stop("counts中仍包含非数值列，请检查文件格式")
    }
    x_num
  })

  log_message(sprintf("counts维度: %d x %d", nrow(df), ncol(df)))

  return(df)
}

read_metadata <- function(file, group_col){
  log_message("读取metadata")

  meta <- tryCatch(read.csv(file), error=function(e) read.delim(file))

  if (is.null(group_col) || !nzchar(group_col)) group_col <- "group"

  if (!"sample" %in% colnames(meta)) stop("缺少sample列")
  if (!(group_col %in% colnames(meta))) stop("缺少group列")

  rownames(meta) <- meta$sample
  meta[[group_col]] <- factor(meta[[group_col]])

  if (nlevels(meta[[group_col]]) < 2) stop("分组不足")

  return(meta)
}

# ---------------- DESeq2 ----------------
run_deseq2 <- function(counts, meta, group_col, control, treatment, min_count, min_samples){

  dds <- DESeqDataSetFromMatrix(counts, meta, design=as.formula(paste("~", group_col)))

  groups <- levels(meta[[group_col]])
  if (is.null(control)) control <- groups[1]
  if (is.null(treatment)) treatment <- groups[2]

  if (!(control %in% groups)) stop("control错误")
  if (!(treatment %in% groups)) stop("treatment错误")

  dds[[group_col]] <- relevel(dds[[group_col]], ref=control)

  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds <- dds[keep,]

  dds <- DESeq(dds, quiet=TRUE)

  contrast <- c(group_col, treatment, control)

  if (requireNamespace("apeglm", quietly = TRUE)) {
  log_message("使用apeglm shrinkage")
  res <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")
  } else {
  log_message("apeglm未安装，使用normal shrinkage", "WARNING")
  res <- lfcShrink(dds, contrast = contrast, type = "normal")
  }

  df <- as.data.frame(res)
  df$gene_id <- rownames(df)

  df$padj[is.na(df$padj)] <- 1
  df$log2FoldChange[is.na(df$log2FoldChange)] <- 0
  df$log2FoldChange <- pmax(pmin(df$log2FoldChange,20),-20)

  df$baseMean <- rowMeans(counts(dds, normalized=TRUE))

  return(list(df=df, dds=dds, group_col=group_col,
              control=control, treatment=treatment))
}

# ---------------- DEG ----------------
get_degs <- function(df, lfc, padj_val){
  df %>% filter(padj < padj_val & abs(log2FoldChange) >= lfc)
}

# ---------------- 画图 ----------------
plot_all <- function(res, degs, outdir, lfc, padj_val){

  df <- res$df
  dds <- res$dds
  group_col <- res$group_col

  df$padj[df$padj==0] <- 1e-300

  # 火山图
  df$sig <- "NotSig"
  df$sig[df$padj < padj_val & df$log2FoldChange >= lfc] <- "Up"
  df$sig[df$padj < padj_val & df$log2FoldChange <= -lfc] <- "Down"
  deg_count <- sum(df$sig != "NotSig")

  p <- ggplot(df, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(color=sig), alpha=0.6, size=1.5) +
    scale_color_manual(values=c("Up"="red","Down"="blue","NotSig"="grey")) +
    geom_hline(yintercept=-log10(padj_val), linetype="dashed", alpha=0.5) +
    geom_vline(xintercept=c(-lfc, lfc), linetype="dashed", alpha=0.5) +
    labs(title=paste0("Volcano Plot (Up:", sum(df$sig=="Up"),
         " Down:", sum(df$sig=="Down"), ")"),
         x=expression(log[2]~FoldChange), y=expression(-log[10]~Padj)) +
    theme_minimal() + theme(legend.position="bottom")

  # 标注top10 DEG基因（需ggrepel包）
  top_genes <- df[df$sig != "NotSig", ]
  if (nrow(top_genes) > 0 && requireNamespace("ggrepel", quietly=TRUE)){
    top_genes <- top_genes[order(top_genes$padj), ]
    top_genes <- head(top_genes, 10)
    p <- p + ggrepel::geom_text_repel(data=top_genes, aes(label=gene_id),
                                       max.overlaps=15, size=3, segment.alpha=0.3)
  }

  ggsave(file.path(outdir,"volcano_plot.png"), p, width=10, height=8, dpi=300)

  # PCA
  vsd <- vst(dds, blind=FALSE)
  pca <- prcomp(t(assay(vsd)))
  pca_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  pca_df <- as.data.frame(pca$x)
  pca_df$group <- colData(vsd)[[group_col]]
  pca_df$sample <- rownames(colData(vsd))

  has_ggrepel <- requireNamespace("ggrepel", quietly=TRUE)
  p2 <- ggplot(pca_df, aes(PC1, PC2, color=group)) +
    geom_point(size=3) +
    labs(title="PCA Plot",
         x=paste0("PC1 (", pca_var[1], "%)"),
         y=paste0("PC2 (", pca_var[2], "%)")) +
    theme_minimal() + theme(legend.position="bottom")

  # 添加95%置信椭圆（每组≥4个样本时才计算）
  if (nrow(pca_df) >= 8 && all(table(pca_df$group) >= 4)) {
    p2 <- p2 + stat_ellipse(aes(fill=group), geom="polygon", alpha=0.1, show.legend=FALSE)
  }

  # 添加样本标签（如果ggrepel可用）
  if (has_ggrepel) {
    p2 <- p2 + ggrepel::geom_text_repel(aes(label=sample), size=3,
                                         show.legend=FALSE, max.overlaps=20)
  }

  ggsave(file.path(outdir,"pca_plot.png"), p2, width=8, height=6, dpi=300)

  # MA图
  p3 <- ggplot(df, aes(log10(baseMean+1), log2FoldChange)) +
    geom_point(aes(color=sig), alpha=0.5, size=1) +
    scale_color_manual(values=c("Up"="red","Down"="blue","NotSig"="grey")) +
    geom_hline(yintercept=0, linetype="solid", alpha=0.3) +
    labs(title=paste0("MA Plot (", deg_count, " DEGs)"),
         x=expression(log[10]~baseMean), y=expression(log[2]~FoldChange)) +
    theme_minimal() + theme(legend.position="bottom")

  ggsave(file.path(outdir,"ma_plot.png"), p3, width=10, height=8, dpi=300)

  # 热图
  if (nrow(degs) > 0){
    genes <- head(degs$gene_id, min(100, nrow(degs)))
    mat <- assay(vsd)[genes, , drop=FALSE]
    mat <- t(scale(t(mat)))
    mat[is.na(mat)] <- 0

    ann <- data.frame(Group=colData(vsd)[[group_col]])
    rownames(ann) <- colnames(vsd)

    png(file.path(outdir,"heatmap.png"), width=800, height=1000)
    pheatmap(mat, annotation_col=ann,
             main=paste0("Top ", nrow(genes), " DEGs"),
             color=colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100),
             show_rownames=FALSE,
             clustering_distance_rows="euclidean",
             clustering_distance_cols="euclidean")
    dev.off()
  }
}

# ---------------- 主函数 ----------------
main <- function(){

  args <- parse_arguments()
  dir.create(args$output, recursive=TRUE, showWarnings=FALSE)

  counts <- read_counts(args$counts)
  meta <- read_metadata(args$metadata, args$`group-col`)

  common <- intersect(colnames(counts), rownames(meta))
  counts <- counts[,common]
  meta <- meta[common,,drop=FALSE]

  res <- run_deseq2(counts, meta, args$`group-col`,
                    args$`control-group`,
                    args$`treatment-group`,
                    args$`min-count`,
                    args$`min-samples`)

  degs <- get_degs(res$df, args$`lfc-threshold`, args$`padj-threshold`)

  write.csv(res$df, file.path(args$output,"deseq2_results.csv"), row.names=FALSE)
  write.csv(degs, file.path(args$output,"filtered_degs.csv"), row.names=FALSE)

  if (args$`generate-plots`) {
    plot_all(res, degs, args$output, args$`lfc-threshold`, args$`padj-threshold`)
  }

  writeLines(capture.output(sessionInfo()), file.path(args$output,"sessionInfo.txt"))

  log_message("完成")
}

if (!interactive()){
  tryCatch(main(), error=function(e){
    log_message(e$message,"ERROR")
    quit(status=1)
  })
}