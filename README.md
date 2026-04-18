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

## 快速开始（推荐WSL2/Linux）

### 一键安装（最简单）

```bash
# 1. 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 2. 给安装脚本添加执行权限
chmod +x setup_complete.sh install_bioc_packages.sh

# 3. 运行一键安装（预计30-45分钟）
./setup_complete.sh
```

### 安装完成后

```bash
# 激活分析环境
conda activate small_rna_analysis

# 查看分析流程
snakemake -n --configfile config/config.yaml

# 运行分析（使用4个CPU核心）
snakemake --cores 4 --configfile config/config.yaml
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
│   ├── raw_fastq/            # 原始测序数据
│   ├── processed/            # 处理后的中间数据
│   └── metadata/             # 样本信息和分组信息
├── references/               # 参考基因组和注释文件
├── scripts/                  # 分析脚本
│   ├── qc/                  # 质量控制脚本
│   ├── alignment/          # 序列比对脚本
│   ├── expression/         # 基因表达分析脚本
│   └── motif/              # Motif分析脚本
├── workflow/                 # Snakemake流程定义
├── config/                   # 配置文件
├── envs/                     # Conda环境配置
├── results/                  # 分析结果
│   ├── qc/
│   ├── alignment/
│   ├── counts/
│   ├── differential_expression/
│   └── motif_analysis/
└── logs/                     # 运行日志
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

### 问：安装卡住了怎么办？

如果在安装Bioconductor数据包（如GenomeInfoDbData）时卡住了，可以按Ctrl+C停止安装，然后单独运行：

```bash
conda activate small_rna_analysis
./install_bioc_packages.sh
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
