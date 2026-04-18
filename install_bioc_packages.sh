#!/bin/bash
# Bioconductor大包快速安装脚本
# 使用清华镜像源加速下载

echo "=========================================="
echo "Bioconductor大包快速安装"
echo "=========================================="

# 检查是否激活了conda环境
if [[ -z "$CONDA_DEFAULT_ENV" || "$CONDA_DEFAULT_ENV" != "small_rna_analysis" ]]; then
    echo "请先激活small_rna_analysis环境："
    echo "  conda activate small_rna_analysis"
    exit 1
fi

# 创建临时下载目录
DOWNLOAD_DIR="/tmp/bioc_packages_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$DOWNLOAD_DIR"
cd "$DOWNLOAD_DIR" || exit 1

echo "下载目录: $DOWNLOAD_DIR"
echo ""

# 定义需要下载的包和它们的清华镜像URL
declare -A BIOC_PACKAGES
BIOC_PACKAGES=(
    ["GenomeInfoDbData_1.2.11.tar.gz"]="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"
)

# 下载并安装包
for pkg in "${!BIOC_PACKAGES[@]}"; do
    url="${BIOC_PACKAGES[$pkg]}"
    echo "正在下载 $pkg ..."

    if wget -q "$url" -O "$pkg"; then
        echo "  ✓ 下载成功"
        echo "  正在安装 $pkg ..."
        R -e "install.packages('$DOWNLOAD_DIR/$pkg', repos = NULL, type = 'source')" 2>&1 | tail -10
        if [ $? -eq 0 ]; then
            echo "  ✓ $pkg 安装成功"
        else
            echo "  ✗ $pkg 安装失败"
        fi
    else
        echo "  ✗ 下载失败: $url"
        echo "  将尝试通过Bioconductor安装..."
        pkg_name=$(echo "$pkg" | cut -d'_' -f1)
        R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor'); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('$pkg_name', ask=FALSE, force=TRUE)"
    fi
    echo ""
done

# 清理
cd - > /dev/null || exit 1
rm -rf "$DOWNLOAD_DIR"

echo "=========================================="
echo "Bioconductor包安装完成！"
echo "=========================================="
