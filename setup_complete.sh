#!/bin/bash

# Small RNA项目完整环境安装脚本
# 适用于WSL2/Linux系统
# 使用mamba替代conda，速度提升10倍
# 配置国内镜像源加速下载
# 解决Bioconductor数据包（如GenomeInfoDbData）下载慢和卡住的问题

print_header() {
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════╗"
    echo "║  $1"
    echo "╚═══════════════════════════════════════════════════════════════╝"
}

print_step() {
    echo ""
    echo "▶ 步骤 $1: $2"
}

print_success() {
    echo "  ✓ $1"
}

print_info() {
    echo "  ℹ $1"
}

print_warning() {
    echo "  ⚠  $1"
}

print_header "Small RNA项目环境安装脚本"

PROJECT_DIR=$(pwd)
print_info "项目目录: $PROJECT_DIR"

# 1. 检测系统
print_step "1" "检测系统信息"
if [ -f /etc/os-release ]; then
    . /etc/os-release
    print_success "操作系统: $PRETTY_NAME"
else
    print_warning "操作系统: 未知"
fi

# 2. 检查Miniconda
print_step "2" "检查Miniconda"
if [ -f "$HOME/miniconda3/bin/conda" ]; then
    print_success "Miniconda已安装"
    export PATH="$HOME/miniconda3/bin:$PATH"
else
    print_warning "Miniconda未找到，正在下载安装..."
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    bash ~/miniconda.sh -b -p ~/miniconda3
    export PATH="$HOME/miniconda3/bin:$PATH"
    conda init bash
    print_success "Miniconda安装完成"
fi

# 3. 配置Conda镜像源
print_step "3" "配置Conda镜像源（清华源）"
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
print_success "Conda镜像源配置完成"

# 4. 配置R镜像源
print_step "4" "配置R镜像源（清华源）"
mkdir -p ~/.R
cat > ~/.Rprofile << 'EOFR'
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
message("已配置使用清华镜像源:")
message("  - Bioconductor: https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
message("  - CRAN: https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
EOFR
print_success "R镜像源配置完成"

# 5. 安装mamba
print_step "5" "检查/安装Mamba"
if ! command -v mamba &> /dev/null; then
    print_info "正在安装Mamba..."
    conda install -y -c conda-forge mamba > /dev/null 2>&1
fi
print_success "Mamba已就绪"

# 6. 删除旧环境
print_step "6" "清理旧环境"
if conda env list | grep -q "small_rna_analysis"; then
    print_info "正在删除旧环境..."
    mamba env remove -n small_rna_analysis -y > /dev/null 2>&1
    print_success "旧环境已删除"
else
    print_info "无旧环境需要清理"
fi

# 7. 预下载Bioconductor大包
print_step "7" "预下载Bioconductor数据包"
mkdir -p /tmp/bioc_packages
cd /tmp/bioc_packages
if [ ! -f "GenomeInfoDbData_1.2.11.tar.gz" ]; then
    print_info "正在下载GenomeInfoDbData（从清华镜像）..."
    wget -q "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz"
    if [ -f "GenomeInfoDbData_1.2.11.tar.gz" ]; then
        print_success "GenomeInfoDbData预下载完成"
    fi
else
    print_success "GenomeInfoDbData已存在，跳过下载"
fi
cd -

# 8. 创建新环境
print_step "8" "创建Conda环境（使用Mamba加速）"
print_info "这可能需要15-25分钟，请耐心等待..."
mamba env create -f envs/small_rna_analysis.yaml

if [ $? -eq 0 ]; then
    print_success "Conda环境创建成功"

    eval "$(conda shell.bash hook)"
    conda activate small_rna_analysis

    # 手动安装预下载的Bioconductor包
    print_step "9" "安装Bioconductor数据包"
    if [ -f "/tmp/bioc_packages/GenomeInfoDbData_1.2.11.tar.gz" ]; then
        print_info "正在安装GenomeInfoDbData..."
        R -e "install.packages('/tmp/bioc_packages/GenomeInfoDbData_1.2.11.tar.gz', repos = NULL, type = 'source')" 2>&1 > /dev/null
        print_success "GenomeInfoDbData安装完成"
    fi

    # 在conda环境中创建.Rprofile
    print_info "配置环境内的R镜像..."
    R_PROFILE_PATH="$CONDA_PREFIX/.Rprofile"
    cat > "$R_PROFILE_PATH" << 'EOFR'
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
EOFR

    # 安装DESeq2和其他Bioconductor包
    print_info "正在安装DESeq2..."
    R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor'); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('DESeq2', ask=FALSE)" 2>&1 > /dev/null
    print_success "DESeq2安装完成"

else
    print_warning "环境创建失败，尝试使用conda..."
    conda env create -f envs/small_rna_analysis.yaml
fi

# 9. 检查apt-get包
print_step "10" "安装系统依赖（可选）"
if command -v apt-get &> /dev/null; then
    print_info "检测到Debian/Ubuntu系统，正在安装工具..."
    sudo apt-get update -qq && sudo apt-get install -qq -y fastqc samtools bowtie2 trimmomatic > /dev/null 2>&1
    print_success "系统依赖安装完成"
else
    print_info "跳过apt-get包安装（需手动安装生物信息学工具）"
fi

# 10. 检查项目文件
print_step "11" "检查项目文件"
if [ -f "references/hg38.fa" ]; then
    print_success "参考基因组已存在"
else
    print_warning "参考基因组未找到"
    print_info "可运行: python download_references.py"
fi

if [ -f "data/metadata/sample_info.csv" ]; then
    print_success "样本信息已存在"
else
    print_warning "样本信息未找到"
fi

# 11. 清理
print_step "12" "清理临时文件"
rm -rf /tmp/bioc_packages
print_success "临时文件已清理"

print_header "环境安装完成！"
echo ""
echo "使用方法："
echo "  1. 激活环境: conda activate small_rna_analysis"
echo "  2. 查看流程: snakemake -n --configfile config/config.yaml"
echo "  3. 运行流程: snakemake --cores 4 --configfile config/config.yaml"
echo ""
