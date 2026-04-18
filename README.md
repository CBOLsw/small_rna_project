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
chmod +x scripts/setup/install_everything.sh

# 3. 运行一键安装（预计45-60分钟）
./scripts/setup/install_everything.sh
```

**脚本会自动执行：**
1. 检查系统环境
2. 创建项目目录结构
3. 检查数据文件
4. 下载参考基因组（如果不存在）
5. 安装完整分析环境（Conda + R包）

### 手动下载参考资源（推荐使用WSL2下载）

#### 1. hg38参考基因组序列
**可用源：UCSC官方源**
- **下载地址**：https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
- **目标文件**：`references/hg38.fa`
- **操作步骤**：
  ```bash
  cd references
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz
  ```

#### 2. hg38基因注释文件
**可用源：UCSC官方源**
- **下载地址**：https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz
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
将您的fastq.gz测序文件放到：
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
3. **序列比对** - 将测序序列比对到参考基因组
4. **基因计数** - 统计每个基因的序列数量
5. **差异表达分析** - 找出GAO组和PAL组之间表达量有显著差异的基因
6. **Motif发现** - 发现差异表达基因中可能的调控序列模式

## 项目数据

### 原始数据
- fastq.gz测序文件，位置：`data/raw_fastq/fastq_files/`
- 样本信息：`data/metadata/sample_info.csv`

**样本分组：**
- GAO组：GAO_1, GAO_2, GAO_3
- PAL组：PAL_1, PAL_2, PAL_3

### 参考数据
- 参考基因组：`references/*.fa`
- 基因注释：`references/*.gtf`
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
│   ├── *.fa                       # 参考基因组
│   ├── *.gtf                      # 基因注释
│   ├── *.mirbase.gff3             # miRBase注释
│   └── bowtie2_index/                # Bowtie2索引目录
├── scripts/                          # 分析脚本
│   ├── qc/                          # 质量控制脚本
│   ├── alignment/                  # 序列比对脚本
│   ├── expression/                 # 基因表达分析脚本
│   ├── motif/                      # Motif分析脚本
│   ├── setup/                      # 环境设置和工具下载脚本
│   │   ├── install_everything.sh    # 一键安装脚本
│   │   ├── install_bioc_packages.sh # R包安装脚本
│   │   ├── setup_complete.sh        # 完成设置脚本
│   └── utils/                      # 工具和辅助脚本
│       ├── download_references.py   # 参考基因组下载脚本
│       └── final_check.py           # 项目状态检查脚本
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
```

## 结果输出

分析结果会保存在以下目录：

- `results/qc/` - 数据质量检查结果
- `results/alignment/` - 序列比对结果
- `results/counts/` - 基因计数结果
- `results/differential_expression/` - 差异表达分析结果
- `results/motif_analysis/` - Motif发现结果

## 参数配置指南

本项目的所有参数都可以在 `config/config.yaml` 文件中进行修改。以下是详细的参数说明：

### 1. 样本信息配置
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| 元数据文件 | `data/metadata/sample_info.csv` | `data/metadata/new_samples.csv` | 样本信息CSV文件路径 | `samples.metadata_file` |
| 分组列名 | `group` | `condition`, `treatment` | 样本分组的列名 | `samples.group_column` |
| 样本列名 | `sample` | `sample_id`, `id` | 样本ID的列名 | `samples.sample_column` |
| 分组列表 | `["GAO", "PAL"]` | `["Control", "Treatment"]`, `["WT", "KO"]` | 比较的实验组名称 | `samples.groups` |

### 2. 参考基因组配置
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| 参考基因组序列 | `references/hg38.fa` | `references/mm10.fa`, `references/rn6.fa` | hg38参考基因组FASTA文件 | `reference.genome_fasta` |
| Bowtie2索引 | `references/bowtie2_index/hg38` | `references/bowtie2_index/mm10` | Bowtie2索引前缀 | `reference.bowtie2_index` |
| 基因注释文件 | `references/hg38.gtf` | `references/mm10.gtf`, `references/rn6.gtf` | GTF格式基因注释 | `reference.gtf_annotation` |

### 3. 质量控制参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| Trimmomatic线程数 | 4 | `2`, `8`, `16` | 并行处理线程数 | `quality_control.trimmomatic.threads` |
| 前端质量阈值 | 3 | `2`, `5`, `10` | 切除5'端低于此质量的碱基 | `quality_control.trimmomatic.leading` |
| 末端质量阈值 | 3 | `2`, `5`, `10` | 切除3'端低于此质量的碱基 | `quality_control.trimmomatic.trailing` |
| 滑动窗口 | `4:15` | `3:10`, `5:20` | 窗口大小:平均质量阈值 | `quality_control.trimmomatic.slidingwindow` |
| 最小序列长度 | 18 | `15`, `20`, `25` | 保留的最短序列长度 | `quality_control.trimmomatic.minlen` |

### 4. 序列比对参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| Bowtie2线程数 | 8 | `4`, `12`, `16` | 比对线程数 | `alignment.bowtie2.threads` |
| 比对模式 | `very-sensitive-local` | `sensitive-local`, `very-fast-local` | 比对预设模式 | `alignment.bowtie2.preset` |
| 种子长度 | 15 | `10`, `20` | 种子序列长度 | `alignment.bowtie2.seed_length` |
| 最大错配数 | 1 | `0`, `2` | 允许的最大错配数 | `alignment.bowtie2.max_mismatches` |
| 报告比对数 | 10 | `5`, `20` | 最多报告的比对位置数 | `alignment.bowtie2.k` |

### 5. 基因计数参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| featureCounts线程数 | 8 | `4`, `12`, `16` | 计数线程数 | `counting.featureCounts.threads` |
| 链特异性 | 0 | `1`, `2` | 0=无链,1=有链,2=反向链 | `counting.featureCounts.strandness` |
| 最小重叠 | 1 | `5`, `10`, `20` | 最小重叠碱基数 | `counting.featureCounts.min_overlap` |
| 多重比对计数 | false | `true` | 是否计数多重比对序列 | `counting.featureCounts.count_multi_mapping` |

### 6. 差异表达分析参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| 调整p值阈值 | 0.05 | `0.01`, `0.1` | 差异表达显著性阈值 | `differential_expression.deseq2.padj_threshold` |
| log2倍数变化阈值 | 1.0 | `0.5`, `1.5`, `2.0` | 表达变化倍数阈值 | `differential_expression.deseq2.log2fc_threshold` |
| 收缩效应 | true | `false` | 是否使用效应量收缩 | `differential_expression.deseq2.shrinkage` |

### 7. Motif分析参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| 最小motif宽度 | 6 | `4`, `8` | 最短motif长度 | `motif_analysis.meme.min_width` |
| 最大motif宽度 | 12 | `10`, `15`, `20` | 最长motif长度 | `motif_analysis.meme.max_width` |
| 最大motif数 | 10 | `5`, `20`, `50` | 最多发现的motif数量 | `motif_analysis.meme.max_motifs` |
| E-value阈值 | 1e-4 | `1e-3`, `1e-5` | Motif显著性阈值 | `motif_analysis.meme.evalue_threshold` |
| MEME线程数 | 8 | `4`, `12`, `16` | Motif分析线程数 | `motif_analysis.meme.threads` |
| Tomtom数据库 | `JASPAR_vertebrates` | `HOCOMOCOv11_core`, `UNIPROBE` | Motif比对数据库 | `motif_analysis.tomtom.database` |

### 8. 执行参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| Snakemake核心数 | 8 | `4`, `12`, `16`, `all` | 并行执行的最大核心数 | `snakemake.cores` |
| 等待时间 | 60 | `30`, `120`, `300` | 文件生成等待时间(秒) | `snakemake.latency_wait` |
| 重试次数 | 2 | `1`, `3`, `5` | 失败任务重试次数 | `snakemake.restart_times` |
| 出错继续 | true | `false` | 某个任务失败时是否继续 | `snakemake.keep_going` |

### 9. 可视化参数
| 参数 | 默认值 | 取值示例/推荐值 | 说明 | 修改位置 |
|------|--------|-----------------|------|----------|
| 图片DPI | 300 | `150`, `600` | 输出图片分辨率 | `visualization.dpi` |
| 图片格式 | `png` | `pdf`, `svg`, `jpg` | 输出图片格式 | `visualization.figure_format` |
| 颜色方案 | `Set2` | `viridis`, `Paired`, `tab10` | 绘图颜色调色板 | `visualization.color_palette` |

### 配置文件修改示例

如需修改参数，编辑 `config/config.yaml` 文件：

```yaml
# 示例：修改差异表达分析阈值
differential_expression:
  deseq2:
    padj_threshold: 0.01        # 将显著性阈值改为0.01
    log2fc_threshold: 1.5        # 将倍数变化阈值改为1.5
```

## 技术实现细节

### 关键参数设置
- Trimmomatic修剪参数：LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
- Bowtie2比对参数：very-sensitive-local，适合small RNA测序
- DESeq2分析参数：padj<0.05，log2FC>1.0
- MEME motif长度：6-12个碱基对

### 分析统计阈值
- 差异表达基因筛选：padj < 0.05，|log2FC| > 1.0
- Motif显著性：E-value < 1e-4

## 流程执行脚本使用说明

`scripts/run_pipeline.py` 是项目的核心执行脚本，提供了丰富的命令和参数来管理分析流程。

### 基本命令格式

```bash
python scripts/run_pipeline.py --config config/config.yaml [选项]
```

### 必需参数

- `--config` 或 `-c`: 指定配置文件路径（必须）
  - 示例: `--config config/config.yaml`

### 运行模式选择（互斥，一次只能使用一个）

- 不指定模式: 运行完整分析流程
- `--module <模块名>`: 仅运行指定的模块
- `--resume`: 从上次失败处恢复执行
- `--dry-run`: 仅显示将要执行的步骤，不实际运行
- `--status`: 检查当前流程的执行状态
- `--check`: 运行项目状态和配置检查

### 可用模块

使用 `--list-modules` 可以查看所有可用模块：

| 模块名 | 说明 |
|---------|------|
| all | 运行完整分析流程（默认） |
| check | 仅运行项目状态检查 |
| qc | 仅运行数据质量控制模块 |
| alignment | 仅运行序列比对模块 |
| counts | 仅运行基因计数模块 |
| de | 仅运行差异表达分析模块 |
| motif | 仅运行motif分析模块 |

### 执行参数

- `--cores <数字>`: 指定使用的CPU核心数（默认使用配置文件中的设置）
  - 示例: `--cores 8`
- `-v` 或 `--verbose`: 显示详细输出信息
- `--log-file <文件路径>`: 指定日志文件路径（默认输出到标准输出）
  - 示例: `--log-file logs/pipeline.log`

### 使用示例

```bash
# 1. 运行完整分析流程（使用默认配置）
python scripts/run_pipeline.py --config config/config.yaml

# 2. 仅运行质量控制模块
python scripts/run_pipeline.py --config config/config.yaml --module qc

# 3. 从上次失败处恢复执行
python scripts/run_pipeline.py --config config/config.yaml --resume

# 4. 检查流程状态
python scripts/run_pipeline.py --config config/config.yaml --status

# 5. 运行项目状态检查
python scripts/run_pipeline.py --config config/config.yaml --check

# 6. Dry-run模式（查看将要执行的步骤）
python scripts/run_pipeline.py --config config/config.yaml --dry-run

# 7. 指定CPU核心数运行
python scripts/run_pipeline.py --config config/config.yaml --cores 12

# 8. 保存日志到文件
python scripts/run_pipeline.py --config config/config.yaml --log-file logs/analysis.log

# 9. 显示所有可用模块
python scripts/run_pipeline.py --config config/config.yaml --list-modules

# 10. 显示详细输出
python scripts/run_pipeline.py --config config/config.yaml --verbose
```

### 模块执行顺序

完整流程按以下顺序执行各模块：
1. 数据质量控制 (qc)
2. 序列比对 (alignment)
3. 基因计数 (counts)
4. 差异表达分析 (de)
5. Motif分析 (motif)

## 常见问题

### 问：如何检查项目状态和配置？

```bash
# 运行项目状态检查，确认配置正确
python scripts/utils/final_check.py

# 或者通过流程执行脚本检查
python scripts/run_pipeline.py --config config/config.yaml --check
```

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
