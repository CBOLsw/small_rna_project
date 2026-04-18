#!/bin/bash
set -e

# ==========================================
# Small RNA项目 - 一键安装脚本
# 整合所有功能：环境安装 + 数据检查 + 参考基因组下载
# ==========================================

# 颜色定义
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# 时间戳函数
get_timestamp() {
    date "+%H:%M:%S"
}

# 打印函数
print_info() {
    echo -e "${BLUE}[$(get_timestamp)] [INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[$(get_timestamp)] [SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[$(get_timestamp)] [WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[$(get_timestamp)] [ERROR]${NC} $1"
}

print_header() {
    echo ""
    echo -e "${GREEN}╔═══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║${NC}  $1"
    echo -e "${GREEN}╚═══════════════════════════════════════════════════════════════╝${NC}"
}

print_step() {
    echo ""
    echo -e "${BLUE}┌─────────────────────────────────────────────────────────────────┐${NC}"
    echo -e "${BLUE}│${NC}  步骤 $1: $2"
    echo -e "${BLUE}└─────────────────────────────────────────────────────────────────┘${NC}"
}

# 获取项目目录
PROJECT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$PROJECT_DIR"

# ==================== 主程序开始 ====================
clear
print_header "Small RNA项目 - 一键安装脚本"
print_info "项目目录: $PROJECT_DIR"

# 检查系统
print_step "1/5" "检查系统环境"
if [ "$(uname)" = "Linux" ] || [ -d "/mnt/c" ]; then
    print_success "检测到WSL2/Linux系统"
else
    print_error "此脚本需要在WSL2/Linux系统上运行"
    print_info "请在WSL2中运行: cd /mnt/c/Users/24584/PycharmProjects/small_rna_project"
    exit 1
fi

# 创建必要的目录
print_step "2/5" "创建项目目录结构"
mkdir -p data/raw_fastq/fastq_files
mkdir -p data/metadata
mkdir -p references/bowtie2_index
mkdir -p results/qc
mkdir -p results/alignment
mkdir -p results/counts
mkdir -p results/differential_expression
mkdir -p results/motif_analysis
mkdir -p logs
print_success "目录结构创建完成"

# 检查数据
print_step "3/5" "检查数据文件"
SAMPLE_CSV="data/metadata/sample_info.csv"
if [ -f "$SAMPLE_CSV" ]; then
    print_success "样本信息文件存在: $SAMPLE_CSV"
    NUM_SAMPLES=$(tail -n +2 "$SAMPLE_CSV" | wc -l)
    print_info "包含 $NUM_SAMPLES 个样本"
else
    print_warning "样本信息文件不存在"
    print_info "请确保 $SAMPLE_CSV 文件存在"
fi

# 检查FASTQ文件
FASTQ_COUNT=$(find data/raw_fastq/fastq_files -name "*.fastq.gz" -o -name "*.fq.gz" 2>/dev/null | wc -l)
if [ "$FASTQ_COUNT" -gt 0 ]; then
    print_success "找到 $FASTQ_COUNT 个FASTQ文件"
else
    print_warning "未找到FASTQ文件"
    print_info "请将测序数据放在 data/raw_fastq/fastq_files/ 目录下"
fi

# 下载参考基因组
print_step "4/5" "下载参考基因组 (可能需要30-60分钟)"
if [ -f "references/hg38.fa" ] && [ -f "references/hg38.gtf" ]; then
    print_success "参考基因组已存在，跳过下载"
else
    print_info "正在下载参考基因组..."
    if command -v python3 &> /dev/null; then
        python3 scripts/utils/download_references.py
    elif command -v python &> /dev/null; then
        python scripts/utils/download_references.py
    else
        print_error "Python未找到，无法下载参考基因组"
        exit 1
    fi
fi

# 安装conda环境
print_step "5/5" "安装分析环境 (可能需要30-45分钟)"
if conda env list | grep -q "small_rna_analysis"; then
    print_warning "conda环境已存在"
    read -p "是否重新安装？(y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "删除旧环境..."
        conda env remove -n small_rna_analysis -y
        print_info "重新安装环境..."
        ./scripts/setup/setup_complete.sh
    else
        print_info "跳过环境安装"
    fi
else
    print_info "开始安装conda环境..."
    ./scripts/setup/setup_complete.sh
fi

# 完成
print_header "安装完成！"
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  下一步操作：${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "  1. 激活环境:"
echo "     conda activate small_rna_analysis"
echo ""
echo "  2. 检查数据:"
echo "     ./check_data.sh"
echo ""
echo "  3. 运行分析:"
echo "     snakemake --cores 4 --configfile config/config.yaml"
echo ""
echo -e "${GREEN}========================================${NC}"
