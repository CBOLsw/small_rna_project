# GAO组与PAL组Small RNA测序数据分析项目

## 项目概述
本项目用于分析GAO组和PAL组的small RNA测序数据，通过完整的生物信息学流程进行差异表达分析和motif发现。

## 目录结构说明

```
.
├── data/                           # 数据目录
│   ├── raw_fastq/                 # 原始fastq测序文件
│   ├── processed/                 # 处理后的中间文件
│   └── metadata/                  # 样本信息文件
├── references/                    # 参考基因组文件
├── scripts/                       # 分析脚本目录
├── results/                       # 分析结果目录
│   ├── qc/                       # 质量控制结果
│   ├── alignment/                # 序列比对结果
│   ├── counts/                   # 基因计数结果
│   ├── differential_expression/  # 差异表达分析结果
│   └── motif_analysis/           # Motif分析结果
├── logs/                         # 运行日志
├── reports/                      # 分析报告
├── envs/                         # 环境配置文件
└── notebooks/                    # Jupyter分析笔记本
```

## 数据准备指南

### 1. 原始数据
将fastq测序文件放入 `data/raw_fastq/` 目录：
- GAO_rep1.fastq.gz
- GAO_rep2.fastq.gz  
- PAL_rep1.fastq.gz
- PAL_rep2.fastq.gz

### 2. 样本信息
在 `data/metadata/` 目录中创建样本信息文件，建议命名为 `sample_info.csv`，包含以下列：
- sample: 样本名称
- group: 组别 (GAO/PAL)
- replicate: 重复编号
- fastq_file: fastq文件名

### 3. 参考数据
将以下文件放入 `references/` 目录：
- 参考基因组序列文件 (如 genome.fa)
- 基因注释文件 (如 genome.gtf)
- Bowtie2索引文件 (在 bowtie2_index/ 子目录中)

## 分析流程
本项目计划执行以下分析步骤：
1. 测序数据质量控制
2. 序列比对到参考基因组
3. 基因计数
4. 差异表达分析 (GAO组 vs PAL组)
5. de novo motif发现 (基于差异表达基因)

## 注意事项
- 请根据实际数据情况调整分析参数
- small RNA测序数据通常较短，比对参数需要相应调整
- 确保有足够的差异表达基因进行motif分析