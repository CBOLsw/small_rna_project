# Small RNA测序分析项目 - 完整使用说明

## 📋 项目概述

本项目用于分析GAO组和PAL组的HeLa细胞small RNA测序数据，通过完整的生物信息学流程进行差异表达分析和motif发现。

**项目名称**: small-rna-seq-analysis-gao-pal-groups

## 🚀 快速开始（推荐WSL2/Linux）

### 1. 环境安装

#### 方案一：一键安装脚本（推荐）

```bash
# 1. 启动WSL（如果尚未启动）
wsl

# 2. 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 3. 运行一键安装脚本
chmod +x setup_complete.sh
./setup_complete.sh
```

#### 方案二：手动配置

```bash
# 1. 配置Conda环境
conda env create -f envs/small_rna_analysis.yaml
conda activate small_rna_analysis

# 2. 安装生物信息学工具（WSL/Linux）
sudo apt-get update && sudo apt-get install -y fastqc samtools bowtie2 trimmomatic

# 3. 下载参考数据
python download_references.py
```

---

### 2. 流程运行

#### 方法一：使用Python脚本（推荐）

```bash
# 激活Conda环境
conda activate small_rna_analysis

# 查看流程概览（试运行）
python scripts/run_pipeline.py --config config/config.yaml --dry-run

# 运行完整流程
python scripts/run_pipeline.py --config config/config.yaml --cores 8

# 仅运行质量控制模块
python scripts/run_pipeline.py --config config/config.yaml --module qc

# 从上次失败处恢复
python scripts/run_pipeline.py --config config/config.yaml --resume

# 检查流程状态
python scripts/run_pipeline.py --config config/config.yaml --status

# 列出可用模块
python scripts/run_pipeline.py --config config/config.yaml --list-modules
```

#### 方法二：直接使用Snakemake

```bash
# 激活Conda环境
conda activate small_rna_analysis

# 查看流程概览
snakemake -n --configfile config/config.yaml

# 运行完整流程
snakemake --cores 8 --configfile config/config.yaml

# 仅运行质量控制模块
snakemake --cores 4 --configfile config/config.yaml results/qc/qc_summary.csv
```

---

## 🔍 详细使用说明

### 1. 数据准备

#### 原始数据
- **位置**: `data/raw_fastq/fastq_files/`
- **格式**: 12个fastq.gz文件
- **样本信息**: `data/metadata/sample_info.csv`

#### 参考数据
- **hg38基因组**: `references/hg38.fa` (881MB)
- **基因注释**: `references/hg38.gtf` (54MB)
- **Bowtie2索引**: `references/bowtie2_index/` (需要运行时生成)

---

### 2. 主要分析模块

#### 质量控制模块 (QC)
```bash
# 运行FastQC和Trimmomatic
python scripts/qc/fastqc_analysis.py --input data/raw_fastq/fastq_files/ --output results/qc
python scripts/qc/trim_fastq.py --input data/raw_fastq/fastq_files/ --output data/processed
python scripts/qc/qc_summary.py --input-dir results/qc --processed-dir data/processed --output results/qc/qc_summary.csv
```

**输出**:
- FastQC报告: `results/qc/*.html`
- 修剪后的序列: `data/processed/*.fastq.gz`
- 质量摘要: `results/qc/qc_summary.csv`

---

#### 序列比对模块 (Alignment)
```bash
# 构建Bowtie2索引
python scripts/alignment/build_bowtie2_index.py --genome references/hg38.fa --output references/bowtie2_index --prefix hg38

# 运行Bowtie2比对
python scripts/alignment/run_bowtie2.py --input data/processed --output results/alignment --index references/bowtie2_index/hg38

# 计算比对统计
python scripts/alignment/alignment_stats.py --input-dir results/alignment --output results/alignment/alignment_summary.csv
```

**输出**:
- 比对结果: `results/alignment/*.sorted.bam`
- 统计信息: `results/alignment/*_alignment_stats.csv`
- 质量评估: `results/alignment/alignment_summary.csv`

---

#### 基因计数模块 (Counts)
```bash
# 使用featureCounts计算基因计数
python scripts/expression/count_features.py --bams results/alignment/*.sorted.bam --annotation references/hg38.gtf --output results/counts/gene_counts.csv

# 生成表达矩阵
python scripts/expression/generate_expression_matrix.py --input results/counts/gene_counts.csv --output results/counts/counts_summary.csv
```

**输出**:
- 基因计数矩阵: `results/counts/gene_counts.csv`
- 表达矩阵: `results/counts/counts_summary.csv`

---

#### 差异表达分析 (Differential Expression)
```bash
# 运行DESeq2分析
Rscript scripts/expression/deseq2_analysis.R --counts results/counts/gene_counts.csv --metadata data/metadata/sample_info.csv --output-dir results/differential_expression --padj-threshold 0.05 --log2fc-threshold 1.0

# 筛选差异表达基因
python scripts/expression/filter_degs.py --input results/differential_expression/deseq2_results.csv --output results/differential_expression/filtered_degs.csv

# 可视化
python scripts/expression/visualize_degs.py --input results/differential_expression --output results/differential_expression
```

**输出**:
- DESeq2结果: `results/differential_expression/deseq2_results.csv`
- 差异基因: `results/differential_expression/filtered_degs.csv`
- 火山图: `results/differential_expression/volcano_plot.png`
- 热图: `results/differential_expression/heatmap.png`

---

#### Motif分析 (Motif Discovery)
```bash
# 提取差异基因序列
python scripts/motif/extract_sequences.py --genes results/differential_expression/filtered_degs.csv --genome references/hg38.fa --annotation references/hg38.gtf --output results/motif_analysis/gene_sequences.fasta

# 运行MEME分析
python scripts/motif/run_meme.py --fasta results/motif_analysis/gene_sequences.fasta --output results/motif_analysis/meme_results

# 过滤motif结果
python scripts/motif/filter_motifs.py --meme results/motif_analysis/meme_results --output results/motif_analysis

# 运行TomTom比较
python scripts/motif/tomtom_comparison.py --meme-output results/motif_analysis/meme_results --output results/motif_analysis/tomtom_results --db JASPAR_vertebrates

# 可视化
python scripts/motif/visualize_motifs.py --meme results/motif_analysis/meme_results --filtered results/motif_analysis --output results/motif_analysis/visualization
```

**输出**:
- 基因序列: `results/motif_analysis/gene_sequences.fasta`
- MEME结果: `results/motif_analysis/meme_results/`
- 筛选结果: `results/motif_analysis/filtered_motifs.csv`
- TomTom比较: `results/motif_analysis/tomtom_results/`
- 可视化: `results/motif_analysis/visualization/`

---

## 📊 结果说明

### 1. 质量控制结果
| 文件 | 内容 |
|------|------|
| `qc_summary.csv` | 所有样本的质量控制统计信息 |
| `*_fastqc.html` | 每个样本的FastQC报告 |
| `*_fastqc.zip` | FastQC原始数据 |

### 2. 比对结果
| 文件 | 内容 |
|------|------|
| `*.sorted.bam` | 排序后的BAM文件 |
| `*_alignment_stats.csv` | 每个样本的比对统计 |
| `alignment_summary.csv` | 整体比对质量评估 |

### 3. 表达分析
| 文件 | 内容 |
|------|------|
| `gene_counts.csv` | 原始基因计数矩阵 |
| `counts_summary.csv` | 基因计数统计 |
| `deseq2_results.csv` | 差异表达分析结果 |
| `filtered_degs.csv` | 筛选后的差异表达基因 |
| `volcano_plot.png` | 差异表达火山图 |
| `heatmap.png` | 基因表达热图 |

### 4. Motif分析
| 文件 | 内容 |
|------|------|
| `gene_sequences.fasta` | 差异基因序列 |
| `meme_results/` | MEME分析结果 |
| `filtered_motifs.csv` | 筛选后的motif |
| `tomtom_results/` | TomTom比较结果 |
| `visualization/motif_overview.png` | Motif可视化 |

---

## ⚙️ 配置说明

### 主要配置文件
- **`config/config.yaml`**: 流程参数配置
- **`envs/small_rna_analysis.yaml`**: Conda环境配置

### 可调整参数
```yaml
# config/config.yaml 关键参数

project_name: "small_rna_analysis_gao_pal"

# 样本信息
samples:
  metadata_file: "data/metadata/sample_info.csv"
  group_column: "group"
  sample_column: "sample_id"

# 参考基因组
reference:
  genome_fasta: "references/hg38.fa"
  gtf_annotation: "references/hg38.gtf"

# 质量控制
quality_control:
  fastqc:
    threads: 4
  trimmomatic:
    threads: 4
    leading: 3
    trailing: 3
    slidingwindow: "4:15"
    minlen: 18

# 序列比对
alignment:
  bowtie2:
    threads: 8
    preset: "very-sensitive-local"

# 差异表达分析
differential_expression:
  deseq2:
    padj_threshold: 0.05
    log2fc_threshold: 1.0

# Motif分析
motif_analysis:
  meme:
    min_width: 6
    max_width: 12
    max_motifs: 10
```

---

## 🔧 常见问题与解决方案

### 1. 环境安装问题

#### Conda镜像源速度慢
```bash
# 配置清华镜像源
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
```

#### 生物信息学工具安装失败
```bash
# 使用apt-get直接安装（WSL/Linux）
sudo apt-get update && sudo apt-get install -y fastqc samtools bowtie2 trimmomatic
```

#### Bioconductor包下载慢（如GenomeInfoDbData）
```bash
# 方案一：使用setup_complete.sh自动配置（推荐）
# setup_complete.sh 脚本会自动在conda环境内配置R镜像源
# 重新运行安装脚本：
./setup_complete.sh

# 方案二：手动在conda环境中配置R镜像源
# 激活conda环境
conda activate small_rna_analysis

# 在conda环境中创建.Rprofile文件
cat > "$CONDA_PREFIX/.Rprofile" << 'EOF'
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
message("已配置使用清华镜像源:")
message("  - Bioconductor: https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
message("  - CRAN: https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
EOF

# 方案三：如果仍然卡住，手动下载并安装大包
# 以GenomeInfoDbData为例：
conda activate small_rna_analysis
R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor'); options(repos=c(CRAN='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('GenomeInfoDbData', ask=FALSE, force=TRUE)"
```

---

### 2. 流程运行问题

#### 依赖库版本冲突
```bash
# 重置Conda环境
conda env remove -n small_rna_analysis -y
conda env create -f envs/small_rna_analysis.yaml
```

#### 流程运行缓慢
```bash
# 调整使用的核心数
python scripts/run_pipeline.py --config config/config.yaml --cores 4
```

#### 缺少参考数据
```bash
# 重新下载参考数据
python download_references.py --force
```

---

## 📈 性能优化建议

### 1. 资源分配
```bash
# 根据系统资源调整核心数
python scripts/run_pipeline.py --config config/config.yaml --cores 16
```

### 2. 并行处理
```bash
# 调整Snakemake参数
snakemake --cores 16 --configfile config/config.yaml --latency-wait 60
```

---

## 📚 项目文档索引

| 文件 | 内容 |
|------|------|
| `README.md` | 项目概述和快速开始 |
| `USAGE_INSTRUCTIONS.md` | **当前文件 - 完整使用说明** |
| `setup_complete.sh` | 一键安装脚本 |
| `docs/environment_setup.md` | 环境配置详解 |
| `docs/technical_manual.md` | 技术实现细节 |
| `docs/user_guide.md` | 详细用户指南 |
| `docs/wsl_setup_guide.md` | WSL配置指南 |
| `config/config.yaml` | 流程参数配置 |

---

## 🔍 项目结构

```
small_rna_project/
├── data/                          # 数据目录
│   ├── raw_fastq/                 # 原始fastq文件
│   ├── processed/                 # 处理后的序列文件
│   └── metadata/                  # 样本信息
├── references/                    # 参考基因组
├── scripts/                       # 分析脚本
│   ├── qc/                        # 质量控制模块
│   ├── alignment/                 # 序列比对模块
│   ├── expression/                # 基因表达分析
│   ├── motif/                     # Motif分析
│   └── run_pipeline.py            # 流程管理脚本
├── workflow/                      # Snakemake流程
├── config/                        # 配置文件
├── envs/                          # Conda环境配置
├── docs/                          # 文档目录
├── reports/                       # 分析报告
├── results/                       # 结果输出
└── setup_complete.sh              # 一键安装脚本
```

---

## 🔄 测试流程

### 数据质量控制测试

```bash
# 运行质量控制模块测试
python scripts/qc/test_qc.py --test-data test_data/
```

### 流程完整性测试

```bash
# 1. 运行质量控制
python scripts/run_pipeline.py --config config/config.yaml --module qc

# 2. 查看结果
ls -la results/qc/

# 3. 检查FastQC报告
ls results/qc/*_fastqc.html
```

---

## 📝 注意事项

### 系统要求
- **Windows 10/11 或 Linux系统**
- **至少8GB内存**
- **至少20GB磁盘空间**
- **网络连接**（用于下载工具和数据）

### 使用限制
- 生物信息学工具（Bowtie2、SAMtools等）在Windows系统上不可用
- 推荐使用WSL2/Linux系统运行完整流程
- 中间文件可能占用大量磁盘空间

### 版本说明
- **Conda环境**: `small_rna_analysis.yaml` (Python 3.9)
- **生物信息学工具**: FastQC 0.12, Bowtie2 2.5.2, SAMtools 1.19, DESeq2 1.40.2

---

## 🎯 技术支持

### 常见问题查询

```bash
# 查看帮助信息
python scripts/run_pipeline.py --help

# 查看具体脚本帮助
python scripts/qc/fastqc_analysis.py --help
```

### 日志文件
- **位置**: `logs/` 目录
- **格式**: 每个步骤都有对应的日志文件

---

## 📄 许可证

本项目仅供学术研究使用。
