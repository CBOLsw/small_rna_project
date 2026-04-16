#!/bin/bash
# 下载参考基因组和注释文件

set -e

REF_DIR="references"
mkdir -p "$REF_DIR"

echo "=== 下载参考基因组数据 (hg38) ==="
echo "开始时间: $(date)"
echo

# 检查文件是否已存在
check_file() {
    if [ -f "$1" ]; then
        echo "✓ 文件已存在: $1"
        return 0
    else
        return 1
    fi
}

# 1. 下载hg38参考基因组 (UCSC版本)
GENOME_FA="$REF_DIR/hg38.fa"
GENOME_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

if check_file "$GENOME_FA"; then
    echo "  跳过下载，文件已存在"
else
    echo "1. 下载hg38参考基因组..."
    echo "   来源: $GENOME_URL"
    wget -O "$REF_DIR/hg38.fa.gz" "$GENOME_URL"
    echo "   解压..."
    gunzip -c "$REF_DIR/hg38.fa.gz" > "$GENOME_FA"
    rm "$REF_DIR/hg38.fa.gz"
    echo "   ✓ 完成: $GENOME_FA"
fi

echo

# 2. 下载基因注释文件 (GENCODE版本)
GTF_FILE="$REF_DIR/gencode.v44.annotation.gtf"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

if check_file "$GTF_FILE"; then
    echo "  跳过下载，文件已存在"
else
    echo "2. 下载GENCODE注释文件..."
    echo "   来源: $GTF_URL"
    wget -O "$REF_DIR/gencode.v44.annotation.gtf.gz" "$GTF_URL"
    echo "   解压..."
    gunzip -c "$REF_DIR/gencode.v44.annotation.gtf.gz" > "$GTF_FILE"
    rm "$REF_DIR/gencode.v44.annotation.gtf.gz"
    echo "   ✓ 完成: $GTF_FILE"
fi

echo

# 3. 下载small RNA注释 (miRBase)
MIRNA_GTF="$REF_DIR/hg38.mirbase.gff3"
MIRNA_URL="https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3"

if check_file "$MIRNA_GTF"; then
    echo "  跳过下载，文件已存在"
else
    echo "3. 下载miRBase small RNA注释..."
    echo "   来源: $MIRNA_URL"
    wget -O "$MIRNA_GTF" "$MIRNA_URL"
    echo "   ✓ 完成: $MIRNA_GTF"
fi

echo

# 4. 创建Bowtie2索引目录
BOWTIE2_INDEX_DIR="$REF_DIR/bowtie2_index"
mkdir -p "$BOWTIE2_INDEX_DIR"

echo "4. 准备Bowtie2索引..."
if [ -f "$BOWTIE2_INDEX_DIR/hg38.1.bt2" ]; then
    echo "  索引文件已存在，跳过构建"
    echo "  索引目录: $BOWTIE2_INDEX_DIR"
else
    echo "  参考基因组文件: $GENOME_FA"
    echo "  索引目录: $BOWTIE2_INDEX_DIR"
    echo "  注意: Bowtie2索引构建需要较长时间和大量内存"
    echo "  请运行以下命令构建索引:"
    echo "    bowtie2-build $GENOME_FA $BOWTIE2_INDEX_DIR/hg38"
    echo "  或使用预构建的索引"
fi

echo
echo "=== 下载完成 ==="
echo "完成时间: $(date)"
echo
echo "文件列表:"
ls -lh "$REF_DIR/"*.fa "$REF_DIR/"*.gtf "$REF_DIR/"*.gff3 2>/dev/null || true
echo
echo "下一步:"
echo "1. 构建Bowtie2索引: bowtie2-build $GENOME_FA $BOWTIE2_INDEX_DIR/hg38"
echo "2. 检查文件完整性"
echo "3. 更新样本信息文件中的参考文件路径"