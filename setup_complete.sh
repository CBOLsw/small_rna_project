#!/bin/bash

# Small RNA项目完整环境安装脚本
# 适用于WSL2/Linux系统
# 使用mamba替代conda，速度提升10倍
# 配置国内镜像源加速下载

echo "=========================================="
echo "Small RNA项目环境安装脚本"
echo "=========================================="

# 项目目录
PROJECT_DIR=$(pwd)
echo "项目目录: $PROJECT_DIR"

# 1. 检测系统
echo ""
echo "1. 检测系统..."
if [ -f /etc/os-release ]; then
    . /etc/os-release
    echo "操作系统: $PRETTY_NAME"
else
    echo "操作系统: 未知"
fi

# 2. 检查Miniconda
echo ""
echo "2. 检查Miniconda..."
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
echo ""
echo "3. 配置Conda镜像源..."
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
echo "镜像源配置完成（使用清华源加速下载）"

# 4. 配置R镜像源（用于Bioconductor包下载）
echo ""
echo "4. 配置R镜像源..."
mkdir -p ~/.R
cat > ~/.Rprofile << 'EOFR'
# R配置文件 - 使用国内镜像源加速Bioconductor包下载

# 设置Bioconductor镜像（清华大学）
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 设置CRAN镜像（清华大学）
options(repos = c(
  CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
))

# 提示信息
message("已配置使用清华镜像源:")
message("  - Bioconductor: https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
message("  - CRAN: https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
EOFR
echo "R镜像源配置完成（使用清华源加速Bioconductor包下载）"

# 5. 安装mamba（快速包管理器）
echo ""
echo "5. 安装mamba（快速包管理器）..."
if ! command -v mamba &> /dev/null; then
    echo "正在安装mamba..."
    conda install -y -c conda-forge mamba > /dev/null 2>&1
    echo "mamba安装完成"
else
    echo "mamba已安装"
fi

# 6. 创建项目环境
echo ""
echo "6. 创建项目环境..."
if conda env list | grep -q "small_rna_analysis"; then
    echo "环境已存在，正在删除..."
    mamba env remove -n small_rna_analysis -y > /dev/null 2>&1
fi

echo "正在创建新环境（使用mamba加速）..."
echo "这可能需要5-10分钟，取决于网络速度..."

# 设置 R 环境变量，确保 Bioconductor 镜像源生效
export R_PROFILE_USER="$HOME/.Rprofile"

# 设置环境变量，确保 R 安装时使用正确的镜像源
export BIOCONDUCTOR_MIRROR="https://mirrors.tuna.tsinghua.edu.cn/bioconductor"
export R_REPOSITORIES="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"

# 执行环境创建，显示完整输出
mamba env create -f envs/small_rna_analysis.yaml

if [ $? -eq 0 ]; then
    echo "环境创建完成"

    # 激活环境并配置R镜像源
    eval "$(conda shell.bash hook)"
    conda activate small_rna_analysis

    # 在conda环境中创建.Rprofile文件
    R_PROFILE_PATH="$CONDA_PREFIX/.Rprofile"
    cat > "$R_PROFILE_PATH" << 'EOFR'
# R配置文件 - 使用国内镜像源加速Bioconductor包下载

# 设置Bioconductor镜像（清华大学）
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

# 设置CRAN镜像（清华大学）
options(repos = c(
  CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
))

# 提示信息
message("已配置使用清华镜像源:")
message("  - Bioconductor: https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
message("  - CRAN: https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
EOFR

    echo "R镜像源配置已写入conda环境"
else
    echo "环境创建失败，尝试使用conda..."
    conda env create -f envs/small_rna_analysis.yaml
fi

# 7. 检查apt-get包（可选）
echo ""
echo "7. 检查生物信息学工具..."
if command -v apt-get &> /dev/null; then
    echo "检测到Debian/Ubuntu系统"
    echo "正在安装生物信息学工具..."
    sudo apt-get update -qq && sudo apt-get install -qq -y fastqc samtools bowtie2 trimmomatic > /dev/null 2>&1
    echo "工具安装完成"
else
    echo "跳过apt-get包安装（需手动安装生物信息学工具）"
fi

# 8. 检查项目文件
echo ""
echo "8. 检查项目文件..."
if [ -f "references/hg38.fa" ]; then
    echo "✅ 参考基因组已存在"
else
    echo "⚠️  参考基因组未找到，正在下载..."
    python download_references.py
fi

if [ -f "data/metadata/sample_info.csv" ]; then
    echo "✅ 样本信息已存在"
else
    echo "⚠️  样本信息未找到"
fi

# 9. 完成提示
echo ""
echo "=========================================="
echo "环境安装完成！"
echo "=========================================="
echo ""
echo "使用方法："
echo "1. 激活环境: conda activate small_rna_analysis"
echo "2. 查看流程: snakemake -n --configfile config/config.yaml"
echo "3. 运行流程: snakemake --cores 8 --configfile config/config.yaml"
echo ""
echo "项目文档请查看: docs/environment_setup.md"
echo ""
