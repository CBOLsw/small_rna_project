#!/bin/bash
# 安装small RNA分析所需的所有生物信息学工具

set -e  # 遇到错误时退出

echo "=== Small RNA分析工具安装脚本 ==="
echo "开始时间: $(date)"
echo

# 检查conda是否安装
if ! command -v conda &> /dev/null; then
    echo "错误: conda未安装。请先安装Miniconda或Anaconda。"
    echo "下载地址: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "1. 创建conda环境: small_rna_analysis"
echo "   使用配置文件: envs/small_rna_analysis.yaml"
echo

# 检查环境文件是否存在
if [ ! -f "envs/small_rna_analysis.yaml" ]; then
    echo "错误: 环境配置文件不存在: envs/small_rna_analysis.yaml"
    echo "请确保在项目根目录运行此脚本"
    exit 1
fi

# 创建conda环境
conda env create -f envs/small_rna_analysis.yaml

echo
echo "2. 激活环境并验证工具安装"
echo

# 验证工具是否安装成功
echo "验证工具版本:"
conda activate small_rna_analysis

echo "Python版本: $(python --version 2>&1)"
echo "FastQC版本: $(fastqc --version 2>&1 | head -1)"
echo "Trimmomatic版本: $(java -jar $(which trimmomatic) -version 2>&1)"
echo "Bowtie2版本: $(bowtie2 --version 2>&1 | head -1)"
echo "samtools版本: $(samtools --version 2>&1 | head -1)"
echo "featureCounts版本: $(featureCounts -v 2>&1)"
echo "R版本: $(R --version 2>&1 | head -1)"
echo "MEME Suite版本: $(meme -version 2>&1 | head -1)"
echo "Snakemake版本: $(snakemake --version 2>&1)"

echo
echo "3. 安装Python包"
echo

# 安装额外的Python包
pip install --upgrade pip
pip install biopython pysam pandas numpy matplotlib seaborn jupyter notebook

echo
echo "=== 安装完成 ==="
echo "完成时间: $(date)"
echo
echo "使用以下命令激活环境:"
echo "  conda activate small_rna_analysis"
echo
echo "要更新环境:"
echo "  conda env update -f envs/small_rna_analysis.yaml"
echo
echo "要删除环境:"
echo "  conda env remove -n small_rna_analysis"