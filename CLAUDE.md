# CLAUDE.md

> **重要提示：这个文件主要是给Claude Code使用的项目指导文档。如果您是用户，请先查看README.md文件了解如何使用本项目。**

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## 项目概述

这是一个生物信息学分析项目，用于分析GAO组和PAL组的HeLa细胞small RNA测序数据。项目使用OpenSpec进行变更管理和工作流程规范。

**核心分析流程**：
1. 数据质量控制（FastQC, Trimmomatic）
2. 序列比对到参考基因组（Bowtie2）
3. 基因计数（featureCounts）
4. 差异表达分析（DESeq2）
5. de novo motif发现（MEME Suite）
6. 结果可视化和报告生成

**技术栈**：
- 主要语言：Python, R, Shell
- 流程管理：Snakemake
- 环境管理：conda
- 项目规范：OpenSpec

## 项目结构

```
├── data/                           # 数据目录
│   ├── raw_fastq/                 # 原始fastq测序文件
│   ├── processed/                 # 处理后的中间文件
│   └── metadata/                  # 样本信息文件
├── references/                    # 参考基因组文件
├── scripts/                       # 分析脚本目录
│   ├── qc/                       # 质量控制脚本
│   ├── alignment/                # 序列比对脚本
│   ├── expression/               # 基因表达分析脚本
│   ├── motif/                    # Motif分析脚本
│   ├── setup/                    # 环境设置和工具下载脚本
│   └── utils/                    # 工具和辅助脚本
├── results/                       # 分析结果目录
│   ├── qc/                       # 质量控制结果
│   ├── alignment/                # 序列比对结果
│   ├── counts/                   # 基因计数结果
│   ├── differential_expression/  # 差异表达分析结果
│   └── motif_analysis/           # Motif分析结果
├── logs/                         # 运行日志
├── reports/                      # 分析报告
├── envs/                         # 环境配置文件
├── notebooks/                    # Jupyter分析笔记本
└── openspec/                     # OpenSpec变更管理
    ├── config.yaml              # OpenSpec配置
    └── changes/                 # 变更规范
        └── small-rna-seq-analysis-gao-pal-groups/
            ├── proposal.md      # 变更提案
            ├── design.md        # 设计文档
            ├── tasks.md         # 任务清单
            └── specs/           # 详细规范
                ├── data-quality-control/
                ├── sequence-alignment/
                ├── differential-expression-analysis/
                ├── motif-discovery/
                └── analysis-pipeline-orchestration/
```

## OpenSpec集成

项目使用OpenSpec进行变更管理。当前活动变更：`small-rna-seq-analysis-gao-pal-groups`

**重要文件**：
- `openspec/changes/small-rna-seq-analysis-gao-pal-groups/tasks.md` - 完整的任务清单
- `openspec/changes/small-rna-seq-analysis-gao-pal-groups/design.md` - 设计决策和工具选择

**使用OpenSpec技能**：
- `/openspec-apply-change` - 开始实现变更任务
- `/openspec-archive-change` - 归档已完成的变更
- `/openspec-explore` - 探索模式，用于调查问题和澄清需求
- `/openspec-propose` - 提出新变更建议

## 常用开发命令

### 环境设置
```bash
# 创建conda环境（基于设计文档中的工具选择）
conda env create -f envs/small_rna_analysis.yaml

# 激活环境
conda activate small_rna_analysis
```

### 流程执行（推荐使用run_pipeline.py）

`run_pipeline.py` 是项目的核心执行脚本，提供了丰富的功能和参数。

```bash
# 基本使用格式
python scripts/run_pipeline.py --config config/config.yaml [选项]

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
python scripts/run_pipeline.py --config config/config.yaml --log-file logs/pipeline.log

# 9. 显示所有可用模块
python scripts/run_pipeline.py --config config/config.yaml --list-modules

# 10. 显示详细输出
python scripts/run_pipeline.py --config config/config.yaml --verbose
```

### run_pipeline.py 参数详细说明

#### 必需参数
| 参数 | 说明 |
|------|------|
| `--config` 或 `-c` | 指定配置文件路径（必须） |

#### 运行模式（互斥，一次只能使用一个）
| 参数 | 说明 |
|------|------|
| 不指定模式 | 运行完整分析流程 |
| `--module <模块名>` | 仅运行指定的模块 |
| `--resume` | 从上次失败处恢复执行 |
| `--dry-run` | 仅显示将要执行的步骤，不实际运行 |
| `--status` | 检查当前流程的执行状态 |
| `--check` | 运行项目状态和配置检查 |

#### 可用模块
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

#### 执行参数
| 参数 | 说明 |
|------|------|
| `--cores <数字>` | 指定使用的CPU核心数（默认使用配置文件中的设置） |
| `-v` 或 `--verbose` | 显示详细输出信息 |
| `--log-file <文件路径>` | 指定日志文件路径（默认输出到标准输出） |

#### 传统Snakemake命令（备选方案）
```bash
# 运行完整分析流程（Snakemake）
snakemake --cores 4 --configfile config/config.yaml

# 运行特定模块
snakemake --cores 4 alignment

# 生成流程图
snakemake --dag | dot -Tpng > workflow.png
```

### 测试命令
```bash
# 运行数据质控测试
python scripts/qc/fastqc_analysis.py --test

# 运行差异表达分析测试
Rscript scripts/expression/deseq2_analysis.R --test

# 运行项目状态检查
python scripts/utils/final_check.py
```

## 架构要点

### 模块化设计
分析流程采用模块化设计，每个模块独立可测试：
1. **数据质控模块** (`scripts/qc/`) - FastQC, Trimmomatic
2. **序列比对模块** (`scripts/alignment/`) - Bowtie2, samtools
3. **基因计数模块** (`scripts/expression/`) - featureCounts
4. **差异表达模块** (`scripts/expression/`) - DESeq2 (R)
5. **motif分析模块** (`scripts/motif/`) - MEME Suite
6. **流程编排模块** (`workflow/`) - Snakemake

### 关键设计决策
1. **比对工具**: Bowtie2（针对small RNA短序列优化）
2. **差异表达分析**: DESeq2（适合小样本量）
3. **motif发现**: MEME Suite（行业标准）
4. **流程管理**: Snakemake（声明式语法，支持并行）
5. **参数管理**: 集中式配置文件 (`config/config.yaml`)

### 数据流
```
原始fastq → 质控 → clean reads → 比对 → BAM文件 → 基因计数 → 表达矩阵 → 差异表达 → 差异基因 → motif发现 → 结果报告
```

## 开发指南

### 实现新模块
1. 参考`openspec/changes/small-rna-seq-analysis-gao-pal-groups/specs/`中的详细规范
2. 按照`tasks.md`中的任务清单逐步实现
3. 确保脚本支持命令行参数和配置文件
4. 添加适当的日志记录和错误处理
5. 实现单元测试和集成测试

### 代码规范
- Python脚本使用snake_case命名
- R脚本使用camelCase命名
- 所有脚本支持`--help`参数显示使用说明
- 关键参数通过配置文件管理，避免硬编码
- 添加详细的文档字符串和注释（使用中文）

### 测试策略
1. **单元测试**: 测试单个函数/模块
2. **集成测试**: 测试模块间数据流
3. **端到端测试**: 使用示例数据运行完整流程
4. **回归测试**: 确保修改不影响现有功能

## 注意事项

1. **数据存储**: 中间文件可能占用20-30GB空间
2. **计算资源**: 比对和motif分析需要多核CPU
3. **参考数据**: 需要hg38参考基因组和基因注释文件
4. **small RNA特性**: 序列短（18-35nt），需要特殊比对参数
5. **可重复性**: 使用conda/Docker确保环境一致性

## 快速开始

1. 设置conda环境：`conda env create -f envs/small_rna_analysis.yaml`
2. 准备参考数据：下载hg38到`references/`
3. 更新样本信息：编辑`data/metadata/sample_info.csv`
4. 运行分析流程：`snakemake --cores 4 --configfile config/config.yaml`
5. 查看结果：`results/` 和 `reports/`