# Small RNA测序分析流程 - 运行说明和参数配置

## 项目概述

这是一个用于分析GAO组和PAL组HeLa细胞small RNA测序数据的生物信息学分析流程。项目使用Snakemake进行流程管理，结合FastQC、Trimmomatic、Bowtie2、Samtools、FeatureCounts和DESeq2等工具完成从原始数据到差异表达分析的完整流程。

## 目录结构

```
small_rna_project/
├── config/                  # 配置文件
├── data/                    # 数据目录
│   ├── raw_fastq/          # 原始数据
│   ├── processed/          # 处理后的中间文件
│   └── metadata/           # 元数据
├── references/             # 参考基因组
├── results/                # 分析结果
│   ├── qc/                 # 质量控制结果
│   ├── alignment/          # 比对结果
│   ├── counts/             # 基因计数
│   ├── differential_expression/ # 差异表达分析结果
│   └── motif_analysis/     # 基序分析结果
├── scripts/                # 分析脚本
│   ├── qc/                 # 质量控制脚本
│   ├── alignment/          # 比对脚本
│   ├── expression/         # 表达分析脚本
│   ├── motif/              # 基序分析脚本
│   ├── utils/              # 工具脚本
│   └── setup/              # 环境设置脚本
└── workflow/               # Snakemake流程规则
```

## 运行前准备

### 1. 环境设置

**✅ 请选择以下**二选一**方法：**

#### 方法1：手动环境配置（推荐给有经验的用户）
```bash
# 1. 创建conda环境
conda create -p ./envs/small_rna_analysis python=3.9

# 2. 激活环境
conda activate ./envs/small_rna_analysis

# 3. 安装所需的包（推荐使用mamba加速安装）
conda install -c bioconda snakemake fastqc trimmomatic bowtie2 samtools bedtools seqkit
conda install -c r r-base r-deseq2
```

#### 方法2：使用项目自带的一键安装脚本（推荐给大多数用户）
```bash
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project
chmod +x scripts/setup/install_everything.sh
./scripts/setup/install_everything.sh
```

**📝 说明：**
- 方法2会自动执行方法1的所有步骤，包括创建conda环境、激活环境和安装所有依赖
- 如果您选择方法2，则不需要再运行方法1的命令
- 如果您已经手动配置了环境（方法1），则不需要再运行方法2
- 一键安装脚本会自动处理依赖关系，避免遗漏重要的安装步骤

## 主要运行方式

### 方式1：使用run_pipeline.py（推荐）

这是项目提供的高级接口，支持多种运行模式和参数配置。

#### 1.1 完整流程运行

```bash
python scripts/run_pipeline.py --config config/config.yaml --cores 8
```

**参数说明：**
- `--config`：指定配置文件（必须）
- `--cores`：指定使用的CPU核心数（推荐4-8）
- `--resume`：从上次失败处恢复运行（可选）
- `--dry-run`：仅显示将要执行的步骤，不实际运行（可选）
- `--status`：检查流程状态（可选）
- `--check`：运行项目状态检查（可选）
- `--module`：指定运行特定模块（可选）
- `--log-file`：指定日志文件路径（可选）

#### 1.2 运行特定模块

```bash
# 仅运行质量控制
python scripts/run_pipeline.py --config config/config.yaml --module qc --cores 8

# 仅运行序列比对
python scripts/run_pipeline.py --config config/config.yaml --module alignment --cores 8

# 仅运行基因计数
python scripts/run_pipeline.py --config config/config.yaml --module counts --cores 8

# 仅运行差异表达分析
python scripts/run_pipeline.py --config config/config.yaml --module de --cores 8

# 仅运行基序分析
python scripts/run_pipeline.py --config config/config.yaml --module motif --cores 8
```

#### 1.3 从失败处恢复

```bash
python scripts/run_pipeline.py --config config/config.yaml --resume --cores 8
```

#### 1.4 查看流程状态

```bash
python scripts/run_pipeline.py --config config/config.yaml --status
```

### 完整可选配置实例

**所有可用选项的综合示例：**
```bash
python scripts/run_pipeline.py \
    --config config/config.yaml \
    --cores 8 \
    --module qc \
    --resume \
    --log-file logs/pipeline.log \
    --dry-run
```

**常用组合：**
```bash
# 快速测试流程
python scripts/run_pipeline.py --config config/config.yaml --dry-run

# 运行完整流程并保存详细日志
python scripts/run_pipeline.py --config config/config.yaml --cores 12 --log-file logs/full_run.log

# 仅运行质量控制和序列比对
python scripts/run_pipeline.py --config config/config.yaml --module qc
python scripts/run_pipeline.py --config config/config.yaml --module alignment --cores 8

# 恢复执行并显示详细输出
python scripts/run_pipeline.py --config config/config.yaml --resume --cores 8 --verbose

# 检查项目状态并生成状态报告
python scripts/run_pipeline.py --config config/config.yaml --check
python scripts/run_pipeline.py --config config/config.yaml --status
```

**run_pipeline.py 所有可选参数详细说明：**

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
| :--- | :--- | :--- | :--- | :--- |
| `--config, -c` | 指定配置文件路径 | 文件路径 | 必填参数 | `--config config/config.yaml` |
| `--cores, -n` | 指定使用的CPU核心数 | 整数 | 4 | `--cores 8` |
| `--module, -m` | 指定运行特定模块 | qc/alignment/counts/de/motif | 无（完整流程） | `--module qc` |
| `--resume, -r` | 从上次失败处恢复运行 | 无参数（开关） | 不恢复 | `--resume` |
| `--dry-run, -d` | 预览执行计划，不实际运行 | 无参数（开关） | 实际运行 | `--dry-run` |
| `--check, -k` | 运行项目状态检查 | 无参数（开关） | 不检查 | `--check` |
| `--status, -s` | 查看流程状态和完成进度 | 无参数（开关） | 不检查 | `--status` |
| `--log-file, -l` | 指定日志文件路径 | 文件路径 | 输出到屏幕 | `--log-file logs/pipeline.log` |
| `--verbose, -v` | 显示详细执行信息 | 无参数（开关） | 正常输出 | `--verbose` |
| `--list-modules` | 列出所有可用的模块 | 无参数（开关） | 不列出 | `--list-modules` |

### 方式2：使用Snakemake直接运行（高级用户）

```bash
# 完整流程
snakemake --cores 8 --configfile config/config.yaml --printshellcmds

# 特定模块
snakemake --cores 8 alignment --configfile config/config.yaml

# 查看将要执行的命令
snakemake --cores 8 --configfile config/config.yaml --dry-run

# 从失败处恢复
snakemake --cores 8 --configfile config/config.yaml --printshellcmds --rerun-incomplete
```

## 配置文件详解（config/config.yaml）

### 1. 项目信息配置

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `project_name` | 项目名称，会影响输出文件名前缀 | 字符串 | `"small_rna_analysis_gao_pal"` | `"small_rna_analysis_2026"` |
| `samples.metadata_file` | 样本元数据文件路径（CSV格式） | 文件路径 | `"data/metadata/sample_info.csv"` | `"data/metadata/samples.csv"` |
| `samples.group_column` | 分组信息所在列名 | 列名 | `"group"` | `"condition"` |
| `samples.sample_column` | 样本ID所在列名 | 列名 | `"sample"` | `"SampleID"` |
| `samples.groups` | 指定的分组（与metadata_file中的分组一致） | 字符串列表 | `["GAO", "PAL"]` | `["control", "treated"]` |

**示例配置：**
```yaml
# 项目名称
project_name: "small_rna_analysis_gao_pal"

# 样本配置
samples:
  metadata_file: "data/metadata/sample_info.csv"
  group_column: "group"
  sample_column: "sample"
  groups:
    - "GAO"
    - "PAL"
```

### 2. 目录配置

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `directories.raw_fastq` | 原始FASTQ文件存放路径 | 文件夹路径 | `"data/raw_fastq"` | `"data/raw"` |
| `directories.processed` | 处理后的中间文件路径 | 文件夹路径 | `"data/processed"` | `"tmp/processed"` |
| `directories.references` | 参考基因组文件路径 | 文件夹路径 | `"references"` | `"data/references"` |
| `directories.results` | 分析结果输出路径 | 文件夹路径 | `"results"` | `"outputs"` |
| `directories.logs` | 日志文件存放路径 | 文件夹路径 | `"logs"` | `"outputs/logs"` |
| `directories.reports` | 报告文件存放路径 | 文件夹路径 | `"reports"` | `"outputs/reports"` |

### 3. 参考基因组配置

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `reference.genome_fasta` | 参考基因组序列（支持.fa和.fa.gz格式） | 文件路径 | `"references/hg38.fa.gz"` | `"references/mm10.fa"` |
| `reference.genome_index` | 基因组索引（自动生成） | 文件路径 | `"references/hg38.fa.fai"` | `"references/mm10.fa.fai"` |
| `reference.bowtie2_index` | Bowtie2索引前缀 | 文件路径 | `"references/bowtie2_index/hg38"` | `"references/bowtie2/mm10"` |
| `reference.gtf_annotation` | 基因注释文件（GTF格式） | 文件路径 | `"references/hg38.gtf"` | `"references/mm10.ensembl.gtf"` |

**注意：** 系统会自动检测并解压压缩的基因组文件（如.gz格式）。

### 4. 质量控制参数

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `quality_control.fastqc.threads` | FastQC使用的线程数 | 整数 | 4 | 8 |
| `quality_control.trimmomatic.threads` | Trimmomatic使用的线程数 | 整数 | 4 | 8 |
| `quality_control.trimmomatic.leading` | 前导质量修剪阈值（小于该值的碱基会被修剪） | 整数(0-40) | 3 | 5 |
| `quality_control.trimmomatic.trailing` | 末尾质量修剪阈值（小于该值的碱基会被修剪） | 整数(0-40) | 3 | 5 |
| `quality_control.trimmomatic.slidingwindow` | 滑动窗口参数（窗口大小:最低质量） | "w:q"格式 | `"4:15"` | `"5:20"` |
| `quality_control.trimmomatic.minlen` | 修剪后序列的最小长度 | 整数 | 18 | 20 |
| `quality_control.trimmomatic.adapter_file` | 接头序列文件路径 | 文件路径 | `"config/VAHTS-SmallRNA-V2.fa"` | `"config/adapters.fa"` |
| `quality_control.trimmomatic.adapter_type` | 接头类型 | vahts_small_rna_v2/illumina | `"vahts_small_rna_v2"` | `"illumina"` |

### 5. 序列比对参数

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `alignment.bowtie2.threads` | Bowtie2使用的线程数 | 整数 | 8 | 12 |
| `alignment.bowtie2.preset` | 比对预设模式 | very-sensitive/sensitive/fast | `"very-sensitive"` | `"sensitive-local"` |
| `alignment.bowtie2.seed_length` | 种子匹配长度 | 整数 | 15 | 20 |
| `alignment.bowtie2.max_mismatches` | 最大错配数 | 整数 | 1 | 2 |
| `alignment.bowtie2.k` | 报告最佳k个比对结果 | 整数 | 10 | 5 |
| `alignment.samtools.threads` | Samtools使用的线程数 | 整数 | 4 | 8 |

### 6. 基因计数参数

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `counting.featureCounts.threads` | FeatureCounts使用的线程数 | 整数 | 8 | 12 |
| `counting.featureCounts.strandness` | 链特异性 | 0=无, 1=有, 2=反向 | 0 | 1 |
| `counting.featureCounts.min_overlap` | 最小重叠长度（碱基） | 整数 | 1 | 10 |
| `counting.featureCounts.count_multi_mapping` | 是否统计多次比对的reads | true/false | false | true |

### 7. 差异表达分析参数

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `differential_expression.deseq2.padj_threshold` | 差异表达的调整后p值阈值 | (0,1) | 0.05 | 0.01 |
| `differential_expression.deseq2.log2fc_threshold` | 差异表达的log2倍数变化阈值 | 正数 | 1.0 | 2.0 |
| `differential_expression.deseq2.shrinkage` | 是否使用shrinkage估计 | true/false | true | false |

### 8. 基序分析参数

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `motif_analysis.meme.min_width` | 基序的最小宽度 | 整数 | 6 | 5 |
| `motif_analysis.meme.max_width` | 基序的最大宽度 | 整数 | 12 | 15 |
| `motif_analysis.meme.max_motifs` | 最大发现的基序数 | 整数 | 10 | 20 |
| `motif_analysis.meme.evalue_threshold` | 基序显著性阈值 | (0,1) | 1e-4 | 1e-3 |
| `motif_analysis.meme.threads` | MEME使用的线程数 | 整数 | 8 | 12 |
| `motif_analysis.tomtom.database` | TomTom数据库名称 | 字符串 | `"JASPAR_vertebrates"` | `"HOCOMOCOv11_core"` |
| `motif_analysis.tomtom.evalue_threshold` | TomTom比对显著性阈值 | (0,1) | 0.05 | 0.01 |
| `motif_analysis.tomtom.min_overlap` | 最小重叠长度 | 整数 | 5 | 4 |
| `motif_analysis.filtering.evalue_threshold` | 基序筛选的E值阈值 | (0,1) | 1e-4 | 1e-3 |
| `motif_analysis.filtering.min_sites` | 基序的最小出现次数 | 整数 | 5 | 10 |
| `motif_analysis.filtering.min_width` | 基序的最小宽度 | 整数 | 6 | 5 |
| `motif_analysis.filtering.max_width` | 基序的最大宽度 | 整数 | 12 | 15 |

### 9. Snakemake参数

| 参数 | 功能 | 可选值 | 默认值 | 示例 |
|------|------|--------|--------|------|
| `snakemake.cores` | 默认使用的CPU核心数 | 整数 | 4 | 12 |
| `snakemake.latency_wait` | 等待输入文件生成的时间（秒） | 整数 | 60 | 120 |
| `snakemake.restart_times` | 失败后重新尝试次数 | 整数 | 2 | 3 |
| `snakemake.keep_going` | 遇到错误时是否继续执行其他任务 | true/false | true | false |

## 运行示例

### 1. 快速开始

```bash
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 1. 检查项目状态
python scripts/run_pipeline.py --config config/config.yaml --check

# 2. 查看流程将要执行的任务
python scripts/run_pipeline.py --config config/config.yaml --dry-run

# 3. 开始完整分析流程（使用8个核心）
python scripts/run_pipeline.py --config config/config.yaml --cores 8
```

### 2. 检查流程状态

```bash
python scripts/run_pipeline.py --config config/config.yaml --status
```

**预期输出示例：**
```
流程状态检查:
  [✓] qc - 完成于 2026-04-19 12:41:02
  [✓] alignment - 完成于 2026-04-19 12:45:30
  [ ] counts - 待执行
  [ ] differential_expression - 待执行
  [ ] motif_analysis - 待执行

已完成: 2/5 个模块
```

### 3. 查看已安装的工具版本

```bash
# 检查核心工具版本
fastqc --version
trimmomatic -version 2>&1 | head -1
bowtie2 --version
samtools --version
```

## 数据分析结果

### 1. 质量控制结果

```
results/qc/
├── fastqc/              # 原始数据质量报告
├── trim/               # 修剪后的质量报告
└── qc_summary.csv      # 质量控制汇总
```

**查看FASTQC报告：**
```bash
ls results/qc/fastqc/*.html
```

### 2. 序列比对结果

```
results/alignment/
├── *.bam               # 比对后的BAM文件
├── *.bai               # BAM索引文件
├── *.sorted.bam        # 排序后的BAM文件
└── alignment_summary.csv # 比对统计
```

**查看比对统计：**
```bash
cat results/alignment/alignment_summary.csv
```

### 3. 基因计数结果

```
results/counts/
├── gene_counts.txt     # 基因计数矩阵
├── *.counts.txt        # 样本计数文件
└── counts_summary.csv
```

### 4. 差异表达分析结果

```
results/differential_expression/
├── deseq2_results.csv  # 完整的差异表达分析结果
├── filtered_degs.csv   # 筛选后的差异表达基因
├── volcano_plot.png    # 火山图
└── heatmap.png         # 热图
```

**查看差异表达基因列表：**
```bash
head -20 results/differential_expression/filtered_degs.csv
```

### 5. 基序分析结果

```
results/motif_analysis/
├── meme_results/       # MEME基序发现结果
├── tomtom_results/     # TomTom基序比对结果
└── motif_analysis.csv  # 基序分析结果
```

## 常见问题和故障排除

### 1. 系统资源问题

**问题：** 流程运行缓慢或卡住

**解决方案：**
- 减少使用的CPU核心数
- 检查内存使用（bowtie2和samtools可能需要大量内存）
- 确保磁盘有足够空间

### 2. 工具版本问题

**问题：** 工具未找到或版本不兼容

**解决方案：**
- 检查conda环境是否正确激活
- 使用conda重新安装所需工具
- 查看`logs/`目录下的错误日志

### 3. 文件权限问题

**问题：** 权限被拒绝

**解决方案：**
- 确保对项目目录有读写权限
- 使用`chmod`命令添加必要的权限

### 4. Snakemake锁问题

**问题：** Snakemake流程被锁定

**解决方案：**
```bash
# 删除锁定文件
rm -f workflow/.snakemake/lock
rm -f .snakemake/lock
```

## 性能优化建议

### 1. 并行化

- 使用`--cores`参数指定更多核心
- 根据系统资源调整`threads`参数

### 2. 内存管理

- Bowtie2索引构建需要大量内存，建议系统至少有8GB内存
- 对于大型数据集，增加`latency_wait`参数

### 3. 存储管理

- 定期清理中间文件（如未使用的BAM文件）
- 使用压缩格式存储大文件

## 更新和维护

### 1. 配置文件修改

**添加新的样本：**
1. 将新的FASTQ文件添加到`data/raw_fastq/`
2. 更新`data/metadata/sample_info.csv`
3. 重新运行流程

**修改参考基因组：**
1. 将新的FASTA文件放置在`references/`目录
2. 更新`config/config.yaml`中的路径
3. 删除现有的索引文件
4. 重新运行流程

### 2. 依赖更新

```bash
# 更新conda环境
conda update --all

# 重新安装工具
conda install -c bioconda snakemake fastqc trimmomatic bowtie2 samtools
```

## 注意事项

1. 原始数据大小可能会影响运行时间
2. 某些步骤（如索引构建）可能需要较长时间
3. 确保有足够的磁盘空间存储中间文件
4. 使用--dry-run参数预览流程步骤
5. 运行过程中会产生大量临时文件

---

**创建时间：** 2026年4月19日
**项目版本：** 1.0
**更新日志：**
- 初始版本创建
- 添加了详细的参数配置说明
- 包含了所有模块的运行示例
- 补充了故障排除指南
