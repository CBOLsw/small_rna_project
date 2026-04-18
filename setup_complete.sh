#!/bin/bash

# Small RNA项目完整环境安装脚本
# 适用于WSL2/Linux系统
# 使用mamba替代conda，速度提升10倍
# 配置国内镜像源加速下载
# 解决Bioconductor数据包（如GenomeInfoDbData）下载慢和卡住的问题

echo "=========================================="
echo "Small RNA项目环境安装脚本"
echo "=========================================="

PROJECT_DIR=$(pwd)
echo "项目目录: $PROJECT_DIR"
echo ""

# 1. 检测系统
if [ -f /etc/os-release ]; then
    . /etc/os-release
    echo "操作系统: $PRETTY_NAME"
else
    echo "操作系统: 未知"
fi

# 2. 检查Miniconda
if [ -f "$HOME/miniconda3/bin/conda" ]; then
    echo "Miniconda已安装"
    export PATH="$HOME/miniconda3/bin:$PATH"
else
    echo "Miniconda未找到，正在下载安装..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash ~/miniconda.sh -b -p ~/miniconda3
    export PATH="$HOME/miniconda3/bin:$PATH"
    conda init bash
    echo "Miniconda安装完成"
fi

# 3. 配置Conda镜像源
cat > ~/.condarc << 'EOF'
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
EOF
echo "conda镜像源配置完成（清华源）"

# 4. 配置R镜像源
mkdir -p ~/.R
cat > ~/.Rprofile << 'EOFR'
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
message("已配置使用清华镜像源:")
message("  - Bioconductor: https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
message("  - CRAN: https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
EOFR
echo "R镜像源配置完成（清华源）"

# 5. 安装mamba
if ! command -v mamba &> /dev/null; then
    echo "正在安装mamba..."
    conda install -y -c conda-forge mamba > /dev/null 2>&1
fi

# 6. 删除旧环境
if conda env list | grep -q "small_rna_analysis"; then
    echo "删除旧环境..."
    mamba env remove -n small_rna_analysis -y > /dev/null 2>&1
fi

# 7. 预下载Bioconductor大包（关键：使用清华镜像）
echo "预下载Bioconductor数据包（避免卡住）..."
mkdir -p /tmp/bioc_packages
cd /tmp/bioc_packages
if [ ! -f "GenomeInfoDbData_1.2.11.tar.gz" ]; then
    echo "正在下载GenomeInfoDbData（清华镜像）..."
    wget -q "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"
fi
cd -

# 8. 创建新环境
echo "正在创建conda环境（使用mamba加速）..."
mamba env create -f envs/small_rna_analysis.yaml

if [ $? -eq 0 ]; then
    echo "环境创建成功"

    eval "$(conda shell.bash hook)"
    conda activate small_rna_analysis

    # 手动安装预下载的Bioconductor包
    if [ -f "/tmp/bioc_packages/GenomeInfoDbData_1.2.11.tar.gz" ]; then
        echo "手动安装GenomeInfoDbData..."
        R -e "install.packages('/tmp/bioc_packages/GenomeInfoDbData_1.2.11.tar.gz', repos = NULL, type = 'source')" 2>&1
    fi

    # 在conda环境中创建.Rprofile
    R_PROFILE_PATH="$CONDA_PREFIX/.Rprofile"
    cat > "$R_PROFILE_PATH" << 'EOFR'
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
EOFR

    # 安装DESeq2和其他Bioconductor包
    echo "正在安装Bioconductor包..."
    R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor'); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DESeq2', ask=FALSE)"

else
    echo "环境创建失败，尝试使用conda..."
    conda env create -f envs/small_rna_analysis.yaml
fi

# 9. 检查apt-get包（可选）
if command -v apt-get &> /dev/null; then
    echo "安装系统依赖..."
    sudo apt-get update -qq && sudo apt-get install -qq -y fastqc samtools bowtie2 trimmomatic > /dev/null 2>&1
fi

# 10. 检查项目文件
echo ""
echo "检查项目文件..."
if [ -f "references/hg38.fa" ]; then
    echo "✓ 参考基因组已存在"
else
    echo "⚠️  参考基因组未找到，正在下载..."
    python download_references.py
fi

if [ -f "data/metadata/sample_info.csv" ]; then
    echo "✓ 样本信息已存在"
else
    echo "⚠️  样本信息未找到"
fi

# 11. 清理
rm -rf /tmp/bioc_packages

echo "=========================================="
echo "环境安装完成！"
echo "=========================================="
echo ""
echo "使用方法："
echo "1. 激活环境: conda activate small_rna_analysis"
echo "2. 查看流程: snakemake -n --configfile config/config.yaml"
echo "3. 运行流程: snakemake --cores 4 --configfile config/config.yaml"
