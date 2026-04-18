#!/bin/bash

# Bioconductor包快速安装脚本
# 从清华镜像源下载，避免官方源慢或卡住的问题

# 获取当前时间戳
get_timestamp() {
    date "+%H:%M:%S"
}

print_header() {
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════╗"
    echo "║  Bioconductor包快速安装"
    echo "║  从清华镜像源下载，避免官方源慢或卡住的问题"
    echo "╚═══════════════════════════════════════════════════════════════╝"
}

print_success() {
    echo "  [$(get_timestamp)] ✓ 成功: $1"
}

print_info() {
    echo "  [$(get_timestamp)] ℹ 信息: $1"
}

print_error() {
    echo "  [$(get_timestamp)] ✗ 错误: $1"
}

print_progress() {
    local current=$1
    local total=$2
    local message=$3
    local percent=$((current * 100 / total))
    echo "  [$(get_timestamp)] [$current/$total] $message ($percent%)"
}

print_header

# 检查conda环境
if [[ -z "$CONDA_DEFAULT_ENV" || "$CONDA_DEFAULT_ENV" != "small_rna_analysis" ]]; then
    print_error "请先激活small_rna_analysis环境："
    print_error "  conda activate small_rna_analysis"
    exit 1
fi

print_info "当前环境: $CONDA_DEFAULT_ENV"

# 安装Bioconductor包
packages=('GenomeInfoDbData' 'DESeq2' 'GenomicRanges' 'IRanges' 'SummarizedExperiment')
total_packages=${#packages[@]}

print_info "开始安装 $total_packages 个Bioconductor包..."

for i in "${!packages[@]}"; do
    package=${packages[$i]}
    current_step=$((i + 1))

    print_progress $current_step $total_packages "正在安装 $package"

    # 直接使用清华镜像安装
    if R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor');
              if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager');
              options(repos = c(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'));
              BiocManager::install('$package', ask=FALSE, force=TRUE, update=FALSE)" 2>&1 > /dev/null; then
        print_success "$package 安装成功"
    else
        print_error "$package 安装失败"
        exit 1
    fi
done

print_info "所有Bioconductor包安装完成！"
echo ""
echo "=========================================="
echo "安装成功！"
echo "=========================================="
