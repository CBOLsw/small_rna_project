#!/bin/bash

# Small RNA项目完整环境安装脚本
# 适用于WSL2/Linux系统
# 配置国内镜像源加速下载

# 获取当前时间戳
get_timestamp() {
    date "+%H:%M:%S"
}

print_header() {
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════╗"
    echo "║  $1"
    echo "╚═══════════════════════════════════════════════════════════════╝"
}

print_step() {
    local step_num=$1
    local step_desc=$2
    echo ""
    echo "┌─────────────────────────────────────────────────────────────────┐"
    echo "│  [$(get_timestamp)] 步骤 $step_num: $step_desc"
    echo "└─────────────────────────────────────────────────────────────────┘"
}

print_success() {
    echo "  [$(get_timestamp)] ✓ 成功: $1"
}

print_info() {
    echo "  [$(get_timestamp)] ℹ 信息: $1"
}

print_warning() {
    echo "  [$(get_timestamp)] ⚠ 警告: $1"
}

print_error() {
    echo "  [$(get_timestamp)] ✗ 错误: $1"
}

print_header "Small RNA项目环境安装脚本"

# 获取项目目录
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
PROJECT_DIR=$(cd "$SCRIPT_DIR/../.." && pwd)
cd "$PROJECT_DIR"
print_info "项目目录: $PROJECT_DIR"

# 确保使用Linux版conda
if [ -f "$HOME/miniconda3/bin/conda" ]; then
    export PATH="$HOME/miniconda3/bin:$PATH"
elif [ -f "$HOME/miniconda/bin/conda" ]; then
    export PATH="$HOME/miniconda/bin:$PATH"
fi

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
elif [ -f "$HOME/miniconda/bin/conda" ]; then
    print_success "Miniconda已安装"
else
    print_warning "Miniconda未找到，正在下载安装..."
    if command -v wget &> /dev/null; then
        wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    elif command -v curl &> /dev/null; then
        curl -s -o ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    else
        print_error "wget或curl都不可用，请先安装其中一个工具"
        exit 1
    fi
    bash ~/miniconda.sh -b -p ~/miniconda3
    export PATH="$HOME/miniconda3/bin:$PATH"
    conda init bash
    print_success "Miniconda安装完成"
fi

# 3. 配置Conda镜像源
print_step "3" "配置Conda镜像源（清华源）"
cat > ~/.condarc << 'EOF'
channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
EOF
print_success "Conda镜像源配置完成"

# 4. 配置R镜像源
print_step "4" "配置R镜像源"
mkdir -p ~/.R
cat > ~/.Rprofile << 'EOFR'
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://bioconductor.org")
EOFR
print_success "R镜像源配置完成"

# 5. 删除旧环境
print_step "5" "清理旧环境和缓存"
if conda env list | grep -q "small_rna_analysis"; then
    print_info "删除旧conda环境..."
    conda env remove -n small_rna_analysis -y 2>/dev/null || true
fi

# 清理残留目录
CONDA_ENV_PATH="$HOME/miniconda3/envs/small_rna_analysis"
if [ -d "$CONDA_ENV_PATH" ]; then
    print_info "清理残留环境目录..."
    rm -rf "$CONDA_ENV_PATH"
fi
CONDA_ENV_PATH2="$HOME/miniconda/envs/small_rna_analysis"
if [ -d "$CONDA_ENV_PATH2" ]; then
    print_info "清理残留环境目录..."
    rm -rf "$CONDA_ENV_PATH2"
fi

print_success "旧环境和缓存清理完成"

# 6. 创建新环境
print_step "6" "创建Conda环境（含DESeq2, apeglm等所有依赖）"
print_info "这可能需要15-25分钟，请耐心等待..."
conda env create -f envs/small_rna_analysis.yaml

if [ $? -eq 0 ]; then
    print_success "Conda环境创建成功"

    eval "$(conda shell.bash hook)"
    conda activate small_rna_analysis

    # 在conda环境中创建.Rprofile
    print_info "配置环境内的R镜像..."
    R_PROFILE_PATH="$CONDA_PREFIX/.Rprofile"
    cat > "$R_PROFILE_PATH" << 'EOFR'
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://bioconductor.org")
EOFR
    print_success "R镜像配置完成"
else
    print_error "Conda环境创建失败"
    print_info "可能的原因："
    print_info "1. 网络连接问题（请检查网络）"
    print_info "2. 磁盘空间不足"
    exit 1
fi

# 7. 安装系统依赖
print_step "7" "安装系统依赖（可选）"
if command -v apt-get &> /dev/null; then
    print_info "检测到Debian/Ubuntu系统，正在安装工具..."
    sudo apt-get update -qq && sudo apt-get install -qq -y fastqc samtools bowtie2 trimmomatic > /dev/null 2>&1
    print_success "系统依赖安装完成"
else
    print_info "跳过apt-get包安装（需手动安装生物信息学工具）"
fi

# 8. 检查项目文件
print_step "8" "检查项目文件"
if [ -f "references/hg38.fa" ] || [ -f "references/hg38.fa.gz" ]; then
    print_success "参考基因组已存在"
else
    print_warning "参考基因组未找到"
    print_info "可运行: python scripts/utils/download_references.py"
fi

if [ -f "references/hg38.knownGene.gtf" ] || [ -f "references/hg38.knownGene.gtf.gz" ]; then
    print_success "基因注释文件已存在"
else
    print_warning "基因注释文件未找到"
    print_info "可运行: python scripts/utils/download_references.py"
fi

if [ -f "data/metadata/sample_info.csv" ]; then
    print_success "样本信息已存在"
else
    print_warning "样本信息未找到"
fi

# 9. 清理临时文件
print_step "9" "清理临时文件"
rm -rf /tmp/bioc_packages
print_success "临时文件已清理"

print_header "环境安装完成！"
echo ""
echo "使用方法："
echo "  1. 激活环境: conda activate small_rna_analysis"
echo "  2. 查看流程: snakemake -n --configfile config/config.yaml"
echo "  3. 运行流程: snakemake --cores 4 --configfile config/config.yaml"
echo ""
