# GAO组与PAL组Small RNA测序数据分析项目

## 项目概述
本项目用于分析GAO组和PAL组的HeLa细胞small RNA测序数据，通过完整的生物信息学流程进行差异表达分析和motif发现。

**变更名称**: small-rna-seq-analysis-gao-pal-groups

## 目录结构说明

```
.
├── data/                           # 数据目录
│   ├── raw_fastq/                 # 原始fastq测序文件
│   ├── processed/                 # 处理后的中间文件
│   └── metadata/                  # 样本信息文件
├── references/                    # 参考基因组文件
│   ├── hg38.fa                    # hg38参考基因组序列
│   ├── hg38.gtf                   # hg38基因注释
│   └── bowtie2_index/             # Bowtie2索引目录
├── scripts/                       # 分析脚本目录
│   ├── qc/                       # 数据质量控制
│   ├── alignment/                # 序列比对
│   ├── expression/               # 基因表达和差异分析
│   ├── motif/                    # Motif分析
│   ├── run_pipeline.py           # 主流程执行脚本
│   └── results_integration.py    # 结果整合
├── workflow/                      # Snakemake流程定义
├── config/                        # 配置文件
├── envs/                          # Conda环境配置
├── docs/                          # 文档目录
│   ├── user_guide.md             # 用户指南
│   └── technical_manual.md       # 技术手册
├── reports/                       # 分析报告
├── results/                       # 分析结果目录
│   ├── qc/                       # 质量控制结果
│   ├── alignment/                # 序列比对结果
│   ├── counts/                   # 基因计数结果
│   ├── differential_expression/  # 差异表达分析结果
│   └── motif_analysis/           # Motif分析结果
├── logs/                         # 运行日志
└── notebooks/                    # Jupyter分析笔记本
```

## 数据准备状态

### 原始数据
- ✅ 已准备13个fastq.gz文件，位置：`data/raw_fastq/fastq_files/`
- ✅ 样本信息已配置：`data/metadata/sample_info.csv`

样本列表：
- GAO_1, GAO_2, GAO_3 (GAO组)
- PAL_1, PAL_2, PAL_3 (PAL组)

### 参考基因组
- ✅ hg38参考基因组 (840.4 MB)：`references/hg38.fa`
- ✅ hg38基因注释 (51.7 MB)：`references/hg38.gtf`
- ⚠️ miRBase注释文件下载失败（404错误）

## 分析流程

本项目执行以下分析步骤：

1. **测序数据质量控制** - FastQC + Trimmomatic
2. **序列比对到参考基因组** - Bowtie2（针对small RNA优化）
3. **基因计数** - featureCounts
4. **差异表达分析** - DESeq2 (GAO组 vs PAL组)
5. **de novo motif发现** - MEME Suite（基于差异表达基因）

## 当前问题和解决方案

### Windows平台限制

**问题**: Bioconda中的生物信息学工具（Bowtie2、SAMtools等）在Windows系统上不可用

**影响**: 
- 无法直接创建conda环境
- 无法构建Bowtie2索引
- 无法运行完整分析流程

**解决方案（按优先级排序）**:

#### 1. 使用WSL2（推荐）
- 安装Windows Subsystem for Linux 2
- 在WSL2中运行完整流程
- 优点: 原生Linux环境，性能好
- 缺点: 需要额外配置

#### 2. 使用Docker
- 使用项目提供的Dockerfile
- 在容器中运行分析
- 优点: 环境一致性好
- 缺点: 需要Docker Desktop

#### 3. 手动安装Windows版本
- 部分工具提供Windows版本
- 需要手动配置路径
- 优点: 不需要额外环境
- 缺点: 配置复杂，部分工具不可用

## 快速开始（Linux环境）

### 1. 环境设置

```bash
# 创建conda环境
conda env create -f envs/small_rna_analysis.yaml
conda activate small_rna_analysis
```

### 2. 构建Bowtie2索引

```bash
python scripts/alignment/build_bowtie2_index.py \
  --genome references/hg38.fa \
  --output references/bowtie2_index/hg38 \
  --threads 4
```

### 3. 运行完整流程

```bash
# 使用Python脚本
python scripts/run_pipeline.py --config config/config.yaml --cores 8

# 或使用Snakemake
snakemake --cores 8 --configfile config/config.yaml
```

## 已完成工作

### 项目清理和准备
- ✅ 清理非核心项目文件
- ✅ 验证项目结构完整性
- ✅ 检查原始数据（13个fastq.gz文件）

### 参考基因组下载
- ✅ 成功下载hg38参考基因组 (840.4 MB)
- ✅ 成功下载hg38基因注释文件 (51.7 MB)

### 模块测试
- ✅ QC模块测试通过（FastQC、Trimmomatic功能验证）
- ✅ 所有Python模块导入测试通过

### 代码完整性
- ✅ 所有分析模块脚本已完成
- ✅ 技术文档和用户指南已创建
- ✅ 配置文件和流程定义已准备
- ✅ Dockerfile已创建

## 文档索引

- `docs/user_guide.md` - 详细使用指南
- `docs/technical_manual.md` - 技术实现细节
- `config/config.yaml` - 分析参数配置
- `envs/small_rna_analysis.yaml` - Conda环境配置
- `openspec/changes/small-rna-seq-analysis-gao-pal-groups/tasks.md` - 完整任务列表

## 注意事项

- 请根据实际数据情况调整分析参数
- small RNA测序数据通常较短（18-35nt），比对参数需要相应调整
- 确保有足够的差异表达基因进行motif分析
- 中间文件可能占用20-30GB空间
- 比对和motif分析需要多核CPU

## 技术栈

- **核心工具**: FastQC, Trimmomatic, Bowtie2, SAMtools, featureCounts, DESeq2, MEME Suite
- **流程管理**: Snakemake
- **编程语言**: Python 3.9+, R 4.3+
- **环境管理**: Conda / Docker

## 许可证

本项目仅供学术研究使用。
