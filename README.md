# GAO组与PAL组Small RNA测序数据分析项目

## 项目介绍

这个项目用于分析GAO组和PAL组的HeLa细胞small RNA测序数据，通过生物信息学方法找出两组细胞之间的差异表达基因，并发现可能的基因调控模式。

## 详细运行说明和参数配置

**重要提示：** 完整的运行说明、参数配置和故障排除指南请查看 [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md) 文件。该文件包含：
- 完整的流程执行指南
- 所有超参数配置说明
- 详细的运行示例
- 故障排除和性能优化建议

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

## 项目结构

```
small_rna_project/
├── data/                           # 数据目录
│   ├── raw_fastq/                 # 原始fastq测序文件
│   │   └── fastq_files/           # FASTQ文件目录
│   ├── processed/                 # 处理后的中间文件
│   └── metadata/                  # 样本信息文件
├── references/                    # 参考基因组文件
│   ├── *.fa                       # 参考基因组
│   ├── *.gtf                      # 基因注释
│   ├── *.mirbase.gff3             # miRBase注释
│   └── bowtie2_index/             # Bowtie2索引目录（运行时生成）
├── scripts/                       # 分析脚本目录
│   ├── qc/                       # 质量控制脚本
│   ├── alignment/                # 序列比对脚本
│   ├── expression/               # 基因表达分析脚本
│   ├── motif/                    # Motif分析脚本
│   ├── setup/                    # 环境设置和工具下载脚本
│   │   ├── install_everything.sh    # 一键安装脚本
│   │   ├── install_bioc_packages.sh # R包安装脚本
│   │   └── setup_complete.sh        # 完整环境安装脚本
│   └── utils/                    # 工具和辅助脚本
│       ├── download_references.py   # 参考基因组下载脚本
│       └── final_check.py           # 项目状态检查脚本
├── workflow/                       # Snakemake流程定义
├── config/                        # 配置文件
├── envs/                          # Conda环境配置
├── results/                       # 分析结果目录
│   ├── qc/                       # 质量控制结果
│   ├── alignment/                # 序列比对结果
│   ├── counts/                   # 基因计数结果
│   ├── differential_expression/  # 差异表达分析结果
│   └── motif_analysis/           # Motif分析结果
├── logs/                         # 运行日志
└── reports/                      # 分析报告
```

## 快速开始（WSL2/Linux）

### 一键安装（推荐）

如果您是第一次使用：

```bash
# 1. 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 2. 给安装脚本添加执行权限
chmod +x scripts/setup/install_everything.sh

# 3. 运行一键安装（预计45-60分钟）
./scripts/setup/install_everything.sh
```

**脚本会自动执行：**
1. 检查系统环境和conda版本
2. 自动安装Linux版Miniconda（如果需要）
3. 创建项目目录结构
4. 检查并下载参考基因组（如果不存在）
5. 创建并配置conda环境
6. 安装Bioconductor包
7. 检查项目状态

### 手动环境配置

如果一键安装失败或需要手动操作：

```bash
# 1. 安装Miniconda（如果未安装）
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init bash
exec bash

# 2. 创建conda环境
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project
conda env create -f envs/small_rna_analysis.yaml

# 3. 激活环境
conda activate small_rna_analysis

# 4. 安装R包
bash scripts/setup/install_bioc_packages.sh

# 5. 下载参考基因组（如果需要）
python scripts/utils/download_references.py

# 6. 检查项目状态
python scripts/utils/final_check.py
```

---

## 数据准备

### 1. 测序数据文件
将您的fastq.gz测序文件放到：
- **目标目录**：`data/raw_fastq/fastq_files/`

**示例文件名**：
```
HeLa-GAO1_S87_L002_R1_001.fastq.gz
HeLa-GAO2_S58_L004_R1_001.fastq.gz
HeLa-PAL1_S51_L001_R1_001.fastq.gz
HeLa-PAL2_S52_L001_R1_001.fastq.gz
...
```

### 2. 参考基因组
如果需要手动更新参考基因组：

```bash
# hg38参考基因组序列
cd references
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# hg38基因注释文件
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
gunzip hg38.knownGene.gtf.gz
mv hg38.knownGene.gtf hg38.gtf

# 文件检查
ls -lh references/
ls -lh data/raw_fastq/fastq_files/
```

---

## 分析流程操作

### 快速开始（已配置环境）

```bash
# 1. 激活分析环境
conda activate small_rna_analysis

# 2. 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 3. 检查项目状态和配置
python scripts/run_pipeline.py --config config/config.yaml --check

# 4. 查看将要执行的步骤
python scripts/run_pipeline.py --config config/config.yaml --dry-run

# 5. 运行完整分析（使用8个CPU核心）
python scripts/run_pipeline.py --config config/config.yaml --cores 8
```

### 常用操作

```bash
# 分模块运行分析
python scripts/run_pipeline.py --config config/config.yaml --module qc      # 质量控制
python scripts/run_pipeline.py --config config/config.yaml --module alignment  # 序列比对
python scripts/run_pipeline.py --config config/config.yaml --module counts    # 基因计数
python scripts/run_pipeline.py --config config/config.yaml --module de        # 差异表达分析
python scripts/run_pipeline.py --config config/config.yaml --module motif     # Motif分析

# 从失败处恢复执行
python scripts/run_pipeline.py --config config/config.yaml --resume

# 查看分析状态和结果
python scripts/run_pipeline.py --config config/config.yaml --status
ls -lh results/
```

### 可用模块

| 模块名 | 说明 |
|---------|------|
| all | 运行完整分析流程（默认） |
| check | 仅运行项目状态检查 |
| qc | 仅运行数据质量控制模块 |
| alignment | 仅运行序列比对模块 |
| counts | 仅运行基因计数模块 |
| de | 仅运行差异表达分析模块 |
| motif | 仅运行motif分析模块 |

---

## 分析流程

项目会自动执行以下分析步骤：

1. **数据质量检查** - 使用FastQC检查测序数据的质量
2. **序列修剪** - 使用Trimmomatic去除低质量序列和接头
3. **序列比对** - 使用Bowtie2将测序序列比对到参考基因组
4. **基因计数** - 使用featureCounts统计每个基因的序列数量
5. **差异表达分析** - 使用DESeq2找出GAO组和PAL组之间的差异表达基因
6. **Motif发现** - 使用MEME Suite发现差异表达基因中的调控序列模式

## 结果输出

分析结果会保存在以下目录：

- `results/qc/` - 数据质量检查结果（FastQC报告）
- `results/alignment/` - 序列比对结果（BAM文件和统计信息）
- `results/counts/` - 基因计数结果（表达矩阵）
- `results/differential_expression/` - 差异表达分析结果（DEGs列表、火山图、热图）
- `results/motif_analysis/` - Motif发现结果（MEME结果、TomTom比较、可视化）

---

## 参数配置

### 主要配置文件

项目的所有参数都可以在 `config/config.yaml` 文件中进行修改：

#### 1. 样本信息配置
```yaml
samples:
  metadata_file: "data/metadata/sample_info.csv"
  group_column: "group"
  sample_column: "sample"
  groups:
    - "GAO"
    - "PAL"
```

#### 2. 质量控制参数
```yaml
quality_control:
  trimmomatic:
    threads: 4
    leading: 3
    trailing: 3
    slidingwindow: "4:15"
    minlen: 18
    adapter_file: "config/VAHTS-SmallRNA-V2.fa"
    adapter_type: "vahts_small_rna_v2"
```

#### 3. 差异表达分析参数
```yaml
differential_expression:
  deseq2:
    padj_threshold: 0.05
    log2fc_threshold: 1.0
    shrinkage: true
```

#### 4. 执行参数
```yaml
snakemake:
  cores: 8
  latency_wait: 60
  restart_times: 2
  keep_going: true
```

---

## 常见问题

### 1. 如何检查项目状态和配置？

```bash
# 运行项目状态检查，确认配置正确
python scripts/utils/final_check.py

# 或者通过流程执行脚本检查
python scripts/run_pipeline.py --config config/config.yaml --check
```

### 2. 如何运行分析？

```bash
conda activate small_rna_analysis
snakemake --cores 4 --configfile config/config.yaml
```

### 3. 如何只运行部分分析？

```bash
# 只运行质量控制
snakemake --cores 4 --configfile config/config.yaml results/qc/qc_summary.csv

# 只运行差异表达分析
snakemake --cores 4 --configfile config/config.yaml results/differential_expression/deseq2_results.csv
```

### 4. 需要多少磁盘空间？

建议至少准备20GB的可用空间，中间文件会占用较多空间。

### 5. 运行需要多长时间？

使用4个CPU核心，完整分析大约需要2-3小时。如果使用8个核心，时间会缩短。

### 6. 系统要求是什么？

- **操作系统**：Windows 10/11（使用WSL2）或 Linux
- **内存**：至少8GB RAM（推荐16GB）
- **磁盘空间**：至少20GB可用空间
- **CPU**：至少4个核心（推荐8个或更多）

---

## 技术实现细节

### 关键设计决策

1. **比对工具**: Bowtie2（针对small RNA短序列优化）
2. **差异表达分析**: DESeq2（适合小样本量）
3. **motif发现**: MEME Suite（行业标准）
4. **流程管理**: Snakemake（声明式语法，支持并行）
5. **参数管理**: 集中式配置文件 (`config/config.yaml`)

### 关键参数设置

- Trimmomatic修剪参数：LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
- Bowtie2比对参数：very-sensitive-local，适合small RNA测序
- DESeq2分析参数：padj<0.05，log2FC>1.0
- MEME motif长度：6-12个碱基对

---

## 许可证

本项目仅供学术研究使用。
