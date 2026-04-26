# motif-discovery Specification

## Purpose
TBD - created by archiving change small-rna-seq-analysis-gao-pal-groups. Update Purpose after archive.
## Requirements
### Requirement: miRNA reads提取
系统应将trimmed reads比对到miRBase参考序列，提取比对上的miRNA reads用于motif发现。

#### Scenario: miRBase比对
- **WHEN** 输入trimmed fastq文件
- **THEN** 系统使用bowtie2将reads比对到miRBase mature miRNA序列，提取比对上的序列

#### Scenario: 长度过滤
- **WHEN** 提取miRNA序列
- **THEN** 系统过滤长度在18-35nt范围内的序列

#### Scenario: 序列去重
- **WHEN** 合并所有样本的miRNA reads
- **THEN** 系统对序列进行精确去重，保留唯一序列用于MEME分析

### Requirement: de novo motif发现
系统应使用MEME Suite在miRNA reads序列中发现富集的序列motif。

#### Scenario: MEME motif发现
- **WHEN** 输入miRNA reads序列
- **THEN** 系统运行MEME分析，发现序列中显著富集的motif，设置motif宽度为5-8nt，最大motif数为3

#### Scenario: 统计显著性评估
- **WHEN** motif发现完成
- **THEN** 系统评估每个motif的统计显著性，默认E-value阈值<1e-4

### Requirement: motif结果可视化和解释
系统应生成motif分析的可视化图表和结果报告。

#### Scenario: motif序列logo生成
- **WHEN** 发现显著motif
- **THEN** MEME生成序列logo图（meme.html），显示motif各位置碱基频率和信息含量

#### Scenario: MEME结果导出
- **WHEN** motif分析完成
- **THEN** 系统导出MEME结果JSON，包含motif序列、宽度、E-value等信息

