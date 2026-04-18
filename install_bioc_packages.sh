#!/bin/bash

# Bioconductor包快速安装脚本
# 从清华镜像源下载，避免官方源慢或卡住的问题

echo "=========================================="
echo "Bioconductor包快速安装"
echo "=========================================="

# 检查conda环境
if [[ -z "$CONDA_DEFAULT_ENV" || "$CONDA_DEFAULT_ENV" != "small_rna_analysis" ]]; then
    echo "请先激活small_rna_analysis环境："
    echo "  conda activate small_rna_analysis"
    exit 1
fi

# 下载并安装Bioconductor包
echo "正在安装Bioconductor包..."

# 使用清华镜像安装GenomeInfoDbData
R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');
# 强制使用清华镜像源
options(repos = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'));
BiocManager::install('GenomeInfoDbData', ask=FALSE, force=TRUE)" 2>&1

# 安装DESeq2和其他核心包
R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');
if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');
BiocManager::install(c('DESeq2', 'GenomicRanges', 'IRanges', 'SummarizedExperiment'), ask=FALSE)" 2>&1

echo "=========================================="
echo "Bioconductor包安装完成！"
echo "=========================================="
