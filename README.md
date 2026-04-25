# Small RNA测序数据分析项目

> **环境要求**：本项目仅在 **Linux / WSL2** 环境下运行，不支持 Windows 原生环境。

## 项目介绍

分析 GAO 组和 PAL 组的 HeLa 细胞 small RNA 测序数据，找出两组之间的差异表达基因，并发现调控序列模式 (motif)。

## 快速开始

### 第一步：环境配置

#### 方法一：一键安装 (推荐)

```bash
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 给脚本添加执行权限
chmod +x scripts/setup/install_everything.sh

# 运行一键安装 (自动完成：环境配置 + 参考基因组下载 + 目录创建)
./scripts/setup/install_everything.sh
```

**一键安装脚本自动执行：**
1. 检测 WSL2/Linux 系统
2. 安装 Linux 版 Miniconda (如需要)
3. 创建项目目录结构
4. 下载参考基因组
5. 安装 conda 分析环境
6. 安装 R 包依赖

#### 方法二：手动安装

```bash
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 1. 创建 conda 环境
conda env create -f envs/small_rna_analysis.yaml
conda activate small_rna_analysis

# 2. 下载参考基因组 (建议手动下载大文件，见下方说明)
python scripts/utils/download_references.py

# 3. 检查项目状态
python scripts/run_pipeline.py --config config/config.yaml --check
```

### 第二步：准备数据

- 测序文件放入 `data/raw_fastq/`
- 样本信息编辑 `data/metadata/sample_info.csv`

### 第三步：运行分析

```bash
# 完整流程
python scripts/run_pipeline.py --config config/config.yaml --cores 8

# 分模块运行
python scripts/run_pipeline.py --config config/config.yaml --module qc --cores 8
python scripts/run_pipeline.py --config config/config.yaml --module alignment --cores 8
python scripts/run_pipeline.py --config config/config.yaml --module counts --cores 8
python scripts/run_pipeline.py --config config/config.yaml --module de --cores 8
python scripts/run_pipeline.py --config config/config.yaml --module motif --cores 8
```

## 参考基因组手动下载

建议手动下载大文件以加快部署：

| 文件 | 大小 | 地址 |
|------|------|------|
| hg38.fa.gz | ~940 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |
| hg38.knownGene.gtf.gz | ~37 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz |

放入 `references/` 目录，脚本会自动检测并解压。

## 分析流程

1. **数据质量控制** - FastQC + Trimmomatic
2. **序列比对** - Bowtie2 → hg38
3. **基因计数** - featureCounts
4. **差异表达分析** - DESeq2
5. **Motif 发现** - MEME Suite + TomTom

## 结果目录

| 目录 | 内容 |
|------|------|
| `results/qc/` | 质量控制报告 |
| `results/alignment/` | BAM文件、比对视图 |
| `results/counts/` | 基因表达矩阵 |
| `results/differential_expression/` | DEGs、火山图、热图 |
| `results/motif_analysis/` | motif结果、可视化 |

## 参数配置

所有参数在 `config/config.yaml` 中修改。详细说明见 [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md)。

## 系统要求

- Windows 10/11 (WSL2) 或 Linux
- 内存 ≥8GB (推荐 16GB)
- 磁盘空间 ≥20GB
- CPU ≥4 核心

## 详细文档

完整说明、参数配置和故障排除请查看 [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md)。
