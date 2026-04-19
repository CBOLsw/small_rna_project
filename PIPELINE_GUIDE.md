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

**推荐使用conda创建专门的分析环境：**
```bash
# 创建conda环境
conda create -p ./envs/small_rna_analysis python=3.9

# 激活环境
conda activate ./envs/small_rna_analysis

# 安装所需的包（推荐使用mamba加速安装）
conda install -c bioconda snakemake fastqc trimmomatic bowtie2 samtools bedtools seqkit
conda install -c r r-base r-deseq2
```

**使用项目自带的安装脚本：**
```bash
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project
chmod +x scripts/setup/install_everything.sh
./scripts/setup/install_everything.sh
```

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

**配置说明：**
- `project_name`：项目名称，会影响输出文件名
- `metadata_file`：样本元数据文件路径（CSV格式）
- `group_column`：分组信息列名（在metadata_file中）
- `sample_column`：样本名列名（在metadata_file中）
- `groups`：指定的分组（与metadata_file中的分组一致）

### 2. 目录配置

```yaml
directories:
  raw_fastq: "data/raw_fastq"
  processed: "data/processed"
  references: "references"
  results: "results"
  logs: "logs"
  reports: "reports"
```

### 3. 参考基因组配置

```yaml
reference:
  genome_fasta: "references/hg38.fa.gz"          # 基因组序列（支持压缩格式）
  genome_index: "references/hg38.fa.fai"         # 基因组索引（自动生成）
  bowtie2_index: "references/bowtie2_index/hg38" # Bowtie2索引前缀
  gtf_annotation: "references/hg38.gtf"          # 基因注释文件
```

**注意：** 系统会自动检测并解压压缩的基因组文件（如.gz格式）。

### 4. 质量控制参数

```yaml
quality_control:
  fastqc:
    threads: 4
  trimmomatic:
    threads: 4
    leading: 3                # 前导质量修剪（小于3的碱基）
    trailing: 3               # 末尾质量修剪（小于3的碱基）
    slidingwindow: "4:15"     # 滑动窗口质量控制（窗口大小:最低质量）
    minlen: 18                # 最小序列长度（小于18的会被丢弃）
    adapter_file: "config/VAHTS-SmallRNA-V2.fa"
    adapter_type: "vahts_small_rna_v2"
```

**Trimmomatic参数说明：**
- `ILLUMINACLIP`：接头去除参数（`文件:最大误配:最小分数:保留最短序列`）
- `SLIDINGWINDOW`：滑动窗口（`窗口大小:最低质量`）
- `LEADING`：前导质量修剪
- `TRAILING`：末尾质量修剪
- `MINLEN`：最低长度过滤

### 5. 序列比对参数

```yaml
alignment:
  bowtie2:
    threads: 8
    preset: "very-sensitive"
    seed_length: 15
    max_mismatches: 1
    k: 10                     # 报告最佳k个比对结果
  samtools:
    threads: 4
```

**Bowtie2参数说明：**
- `preset`：比对预设模式（very-sensitive, sensitive, fast等）
- `seed_length`：种子匹配长度
- `max_mismatches`：最大错配数
- `k`：报告的最佳比对数

### 6. 基因计数参数

```yaml
counting:
  featureCounts:
    threads: 8
    strandness: 0             # 链特异性（0=无，1=有，2=反向）
    min_overlap: 1
    count_multi_mapping: false
```

### 7. 差异表达分析参数

```yaml
differential_expression:
  deseq2:
    padj_threshold: 0.05
    log2fc_threshold: 1.0
    shrinkage: true
```

**DESeq2参数说明：**
- `padj_threshold`：调整后p值阈值（用于筛选差异表达基因）
- `log2fc_threshold`：log2倍数变化阈值
- `shrinkage`：是否使用shrinkage估计（推荐使用）

### 8. 基序分析参数

```yaml
motif_analysis:
  meme:
    min_width: 6
    max_width: 12
    max_motifs: 10
    evalue_threshold: 1e-4
    threads: 8
  tomtom:
    database: "JASPAR_vertebrates"
    evalue_threshold: 0.05
    min_overlap: 5
  filtering:
    evalue_threshold: 1e-4
    min_sites: 5
    min_width: 6
    max_width: 12
```

### 9. Snakemake参数

```yaml
snakemake:
  cores: 4
  latency_wait: 60
  restart_times: 2
  keep_going: true
```

**Snakemake参数说明：**
- `cores`：默认核心数
- `latency_wait`：等待输入文件生成的时间（秒）
- `restart_times`：失败后重新尝试次数
- `keep_going`：遇到错误时是否继续执行其他任务

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
