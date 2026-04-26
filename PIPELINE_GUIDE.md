# Small RNA测序分析流程 - 运行说明和参数配置

## 项目概述

这是一个用于分析GAO组和PAL组HeLa细胞small RNA测序数据的生物信息学分析流程。项目使用Snakemake进行流程管理，结合FastQC、Trimmomatic、Bowtie2、Samtools、featureCounts和DESeq2等工具完成从原始数据到差异表达分析的完整流程。

### 核心分析步骤

1. **数据质量控制** - FastQC + Trimmomatic 质量评估和序列修剪
2. **序列比对** - Bowtie2 将测序 reads 比对到参考基因组 hg38
3. **基因计数** - featureCounts 统计每个基因的 reads 数目
4. **差异表达分析** - DESeq2 找出 GAO 组和 PAL 组之间的差异表达基因
5. **Motif 发现** - MEME Suite 对 miRNA reads 进行 de novo motif 发现

## 目录结构

```
small_rna_project/
├── config/                  # 配置文件
│   └── config.yaml         # 主配置文件
├── data/                    # 数据目录
│   ├── raw_fastq/          # 原始FASTQ文件
│   ├── processed/          # 修剪后的FASTQ文件
│   └── metadata/           # 样本元数据 (sample_info.csv)
├── references/             # 参考基因组文件
│   ├── hg38.fa(.gz)       # hg38 基因组序列
│   ├── hg38.knownGene.gtf(.gz)  # 基因注释
│   ├── hsa.mature.fa       # 人类成熟miRNA序列 (自动下载)
│   └── bowtie2_index/     # Bowtie2 索引 (自动构建)
├── scripts/                # 分析脚本
│   ├── qc/               # 质量控制脚本
│   ├── alignment/         # 序列比对脚本
│   ├── expression/        # 表达分析脚本
│   ├── motif/            # Motif分析脚本
│   ├── utils/            # 工具脚本
│   └── run_pipeline.py   # 流程执行入口
├── workflow/              # Snakemake 流程定义
├── envs/                  # Conda 环境配置
│   └── small_rna_analysis.yaml  # 环境依赖文件
└── results/               # 分析结果 (自动生成)
```

## 快速开始

### 第一步：一键安装 (推荐)

```bash
cd small_rna_project

# 给脚本添加执行权限
chmod +x scripts/setup/install_everything.sh

# 运行一键安装 (自动完成：环境配置 + 参考基因组下载 + 目录创建)
./scripts/setup/install_everything.sh
```

### 手动安装

```bash
cd small_rna_project

# 1. 创建并激活 conda 环境
conda env create -f envs/small_rna_analysis.yaml
conda activate small_rna_analysis

# 2. 下载参考基因组 (建议手动下载大文件，见下方说明)
python scripts/utils/download_references.py --output references

# 3. 检查项目状态
python scripts/run_pipeline.py --config config/config.yaml --check
```

### 一键安装脚本功能

`install_everything.sh` 自动执行：
1. 检测 WSL2/Linux 系统
2. 安装 Linux 版 Miniconda (如需要)
3. 创建项目目录结构
4. 下载参考基因组
5. 安装 conda 分析环境
6. 安装 R 包依赖

### 第二步：准备数据

1. **放置测序文件**到 `data/raw_fastq/`
2. **配置样本信息**编辑 `data/metadata/sample_info.csv`

样本信息文件格式示例：
```csv
sample,group,fastq_r1
GAO_1,GAO,data/raw_fastq/HeLa-GAO1_S87_L002_R1_001.fastq.gz
GAO_2,GAO,data/raw_fastq/HeLa-GAO2_S58_L004_R1_001.fastq.gz
GAO_3,GAO,data/raw_fastq/HeLa-GAO3_S59_L004_R1_001.fastq.gz
PAL_1,PAL,data/raw_fastq/HeLa-PAL1_S51_L001_R1_001.fastq.gz
PAL_2,PAL,data/raw_fastq/HeLa-PAL2_S52_L001_R1_001.fastq.gz
PAL_3,PAL,data/raw_fastq/HeLa-PAL3_S53_L001_R1_001.fastq.gz
```

### 第三步：运行分析

```bash
# 完整流程运行 (使用 8 个 CPU 核心)
python scripts/run_pipeline.py --config config/config.yaml --cores 8

# 分模块运行
python scripts/run_pipeline.py --config config/config.yaml --module qc --cores 8    # 质量控制
python scripts/run_pipeline.py --config config/config.yaml --module alignment --cores 8  # 序列比对
python scripts/run_pipeline.py --config config/config.yaml --module counts --cores 8    # 基因计数
python scripts/run_pipeline.py --config config/config.yaml --module de --cores 8       # 差异表达
python scripts/run_pipeline.py --config config/config.yaml --module motif --cores 8    # Motif分析

# 查看流程状态
python scripts/run_pipeline.py --config config/config.yaml --status

# 从失败处恢复
python scripts/run_pipeline.py --config config/config.yaml --resume --cores 8
```

## 参考基因组说明

### 建议手动下载的文件

由于 hg38 基因组文件较大 (~3GB)，建议手动下载以加快部署速度：

| 文件 | 大小 | 下载地址 |
|------|------|----------|
| hg38.fa.gz | ~940 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |
| hg38.knownGene.gtf.gz | ~37 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz |

**放置位置：**
```
references/
├── hg38.fa.gz                 # 放入此处
└── hg38.knownGene.gtf.gz     # 放入此处
```

下载后，脚本会自动检测并解压。

### 自动下载的文件

miRBase small RNA 注释文件较小 (~500KB)，脚本会直接下载：
```
references/hsa.mature.fa
```

## run_pipeline.py 参数说明

| 参数 | 功能 | 示例 |
|------|------|------|
| `--config` | 配置文件路径 (必填) | `--config config/config.yaml` |
| `--cores` | CPU 核心数 | `--cores 8` |
| `--module` | 运行特定模块 | `--module qc` |
| `--resume` | 从失败处恢复 | `--resume` |
| `--dry-run` | 预览流程 | `--dry-run` |
| `--status` | 查看执行状态 | `--status` |
| `--check` | 项目状态检查 | `--check` |
| `--list-modules` | 列出可用模块 | `--list-modules` |

**可用模块：** `qc` | `alignment` | `counts` | `de` | `motif`

## 配置文件详解

所有参数集中在 `config/config.yaml` 中修改。

### 主要配置项

```yaml
# 样本配置
samples:
  metadata_file: "data/metadata/sample_info.csv"
  group_column: "group"
  sample_column: "sample"
  groups: ["GAO", "PAL"]

# 质量控制
quality_control:
  trimmomatic:
    threads: 4
    minlen: 18            # small RNA 最小长度
    slidingwindow: "4:15"

# 序列比对
alignment:
  bowtie2:
    threads: 8
    preset: "very-sensitive-local"
    k: 10

# 差异表达
differential_expression:
  deseq2:
    padj_threshold: 0.05
    log2fc_threshold: 1.0
    control_group: "GAO"
    treatment_group: "PAL"

# Small RNA Motif 分析
motif_analysis:
  mirbase:
    min_len: 18
    max_len: 35
    max_mismatches: 1
  meme:
    min_width: 5
    max_width: 8
    max_motifs: 3
    evalue_threshold: 1e-4
    minsites: 10
    maxsites: 100
    searchsize: 100000
```

## 分析结果说明

流程完成后，结果保存在 `results/` 目录：

| 目录 | 内容 |
|------|------|
| `qc/` | FastQC 质量报告、修剪后reads统计 |
| `alignment/` | BAM 文件、比对统计信息 |
| `counts/` | 基因表达矩阵 (gene_counts.csv) |
| `differential_expression/` | DEGs列表、火山图、热图 |
| `small_rna_motif/` | Small RNA motif分析结果（miRBase比对、唯一reads、MEME de novo motifs）|

## 常见问题

### 1. Snakemake 锁定
如果流程中断后无法重启：
```bash
rm -rf .snakemake/lock
```

### 2. 参考基因组下载慢
建议手动下载 `.gz` 文件放到 `references/` 目录，脚本会自动跳过下载直接解压。

### 3. Snakemake 锁定
如果流程中断后无法重启：
```bash
rm -rf .snakemake/lock
```

### 4. 内存不足
减少 `--cores` 参数，或在 `config.yaml` 中减小线程数。

## 系统要求

- **操作系统：** Windows 10/11 (WSL2) 或 Linux
- **内存：** ≥8GB (推荐 16GB)
- **磁盘空间：** ≥20GB
- **CPU：** ≥4 核心 (推荐 8 核心)

## 技术栈

| 类别 | 工具 |
|------|------|
| 质量控制 | FastQC, Trimmomatic |
| 序列比对 | Bowtie2, SAMtools |
| 基因计数 | featureCounts |
| 差异表达 | DESeq2 (R) |
| Motif 发现 | MEME Suite (de novo motif发现) |
| 流程管理 | Snakemake |
| 环境管理 | conda |

---
**创建时间：** 2026年4月
**最后更新：** 2026年4月25日
