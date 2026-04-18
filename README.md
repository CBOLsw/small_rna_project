# GAO组与PAL组Small RNA测序数据分析项目

## 项目介绍

这个项目用于分析GAO组和PAL组的HeLa细胞small RNA测序数据，通过生物信息学方法找出两组细胞之间的差异表达基因，并发现可能的基因调控模式。

## 技术栈

### 核心分析工具
- **FastQC**: 测序数据质量评估
- **Trimmomatic**: 测序数据修剪和过滤
- **Bowtie2**: 序列比对工具（专门为small RNA优化）
- **SAMtools**: BAM文件处理
- **featureCounts**: 基因表达量计数
- **DESeq2**: 差异表达分析（R语言包）
- **MEME Suite**: 转录因子结合位点分析（motif发现）

### 流程管理和编程
- **Snakemake**: 工作流程管理
- **Conda/Mamba**: 环境管理和包管理
- **Python 3.9+**: 数据分析和脚本编写
- **R 4.3+**: 统计分析和可视化

## 快速开始（WSL2/Linux）

### 一键安装（推荐）

```bash
# 1. 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 2. 给安装脚本添加执行权限
chmod +x install_everything.sh

# 3. 运行一键安装（预计45-60分钟）
./install_everything.sh
```

**脚本会自动执行：**
1. 检查系统环境
2. 创建项目目录结构
3. 检查数据文件
4. 下载参考基因组（如果不存在）
5. 安装完整分析环境（Conda + R包）

### 手动下载参考资源（如果自动下载慢）

如果自动下载太慢，可以手动下载以下文件并放到指定目录：

#### 1. hg38参考基因组序列
- **下载地址**：https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
- **备用地址**：https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc.fa.gz
- **目标文件**：`references/hg38.fa`
- **操作步骤**：
  ```bash
  # 下载并解压
  cd references
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz
  ```

#### 2. hg38基因注释文件
- **下载地址**：https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
- **备用地址**：https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.gff.gz
- **目标文件**：`references/hg38.gtf`
- **操作步骤**：
  ```bash
  cd references
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
  gunzip hg38.knownGene.gtf.gz
  mv hg38.knownGene.gtf hg38.gtf
  ```

#### 3. miRBase small RNA注释（可选）
- **下载地址**：https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3
- **备用地址**：ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
- **目标文件**：`references/hg38.mirbase.gff3`
- **操作步骤**：
  ```bash
  cd references
  wget https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3
  mv hsa.gff3 hg38.mirbase.gff3
  ```

#### 4. 测序数据文件（必须）
将您的12个fastq.gz测序文件放到：
- **目标目录**：`data/raw_fastq/fastq_files/`

**示例文件名**：
```
HeLa-GAO1_S87_L002_R1_001.fastq.gz
HeLa-GAO1_S87_L002_R2_001.fastq.gz
HeLa-GAO2_S58_L004_R1_001.fastq.gz
HeLa-GAO2_S58_L004_R2_001.fastq.gz
...
```

**文件检查**：
下载完成后，运行以下命令检查文件完整性：
```bash
ls -lh references/
ls -lh data/raw_fastq/fastq_files/
```

## 分析流程

项目会自动执行以下分析步骤：

1. **数据质量检查** - 检查测序数据的质量
2. **序列修剪** - 去除低质量的序列和接头
3. **序列比对** - 将测序序列比对到人类参考基因组
4. **基因计数** - 统计每个基因的序列数量
5. **差异表达分析** - 找出GAO组和PAL组之间表达量有显著差异的基因
6. **Motif发现** - 发现差异表达基因中可能的调控序列模式

## 项目数据

### 原始数据
- 12个fastq.gz测序文件，位置：`data/raw_fastq/fastq_files/`
- 样本信息：`data/metadata/sample_info.csv`

**样本分组：**
- GAO组：GAO_1, GAO_2, GAO_3
- PAL组：PAL_1, PAL_2, PAL_3

### 参考数据
- hg38参考基因组：`references/hg38.fa`
- hg38基因注释：`references/hg38.gtf`
- Bowtie2索引：`references/bowtie2_index/`（运行时生成）

## 项目依赖的配置文件

- `config/config.yaml`: 分析流程参数配置
- `envs/small_rna_analysis.yaml`: Conda环境配置

## 项目目录详细说明

```
small_rna_project/
├── data/
│   ├── raw_fastq/                    # 原始测序数据
│   │   └── fastq_files/              # FASTQ文件目录
│   ├── processed/                    # 处理后的中间数据
│   └── metadata/                     # 样本信息和分组信息
├── references/                       # 参考基因组和注释文件
│   ├── hg38.fa                       # hg38参考基因组
│   ├── hg38.gtf                      # hg38基因注释
│   ├── hg38.mirbase.gff3             # miRBase注释
│   └── bowtie2_index/                # Bowtie2索引目录
├── scripts/                          # 分析脚本
│   ├── qc/                          # 质量控制脚本
│   ├── alignment/                  # 序列比对脚本
│   ├── expression/                 # 基因表达分析脚本
│   └── motif/                      # Motif分析脚本
├── workflow/                        # Snakemake流程定义
├── config/                          # 配置文件
├── envs/                           # Conda环境配置
├── results/                        # 分析结果
│   ├── qc/
│   ├── alignment/
│   ├── counts/
│   ├── differential_expression/
│   └── motif_analysis/
├── logs/                          # 运行日志
└── install_everything.sh         # 一键安装脚本
```

## 结果输出

分析结果会保存在以下目录：

- `results/qc/` - 数据质量检查结果
- `results/alignment/` - 序列比对结果
- `results/counts/` - 基因计数结果
- `results/differential_expression/` - 差异表达分析结果
- `results/motif_analysis/` - Motif发现结果

## 技术实现细节

### 关键参数设置
- Trimmomatic修剪参数：LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
- Bowtie2比对参数：very-sensitive-local，适合small RNA测序
- DESeq2分析参数：padj<0.05，log2FC>1.0
- MEME motif长度：6-12个碱基对

### 分析统计阈值
- 差异表达基因筛选：padj < 0.05，|log2FC| > 1.0
- Motif显著性：E-value < 10

## 常见问题

### 问：如何运行分析？

```bash
conda activate small_rna_analysis
snakemake --cores 4 --configfile config/config.yaml
```

### 问：如何只运行部分分析？

```bash
# 只运行质量控制
snakemake --cores 4 --configfile config/config.yaml results/qc/qc_summary.csv

# 只运行差异表达分析
snakemake --cores 4 --configfile config/config.yaml results/differential_expression/deseq2_results.csv
```

### 问：需要多少磁盘空间？

建议至少准备20GB的可用空间，中间文件会占用较多空间。

### 问：运行需要多长时间？

使用4个CPU核心，完整分析大约需要2-3小时。如果使用8个核心，时间会缩短。

### 问：系统要求是什么？

- **操作系统**：Windows 10/11（使用WSL2）或 Linux
- **内存**：至少8GB RAM（推荐16GB）
- **磁盘空间**：至少20GB可用空间
- **CPU**：至少4个核心（推荐8个或更多）

## 许可证

本项目仅供学术研究使用。
