# differential-expression-analysis Specification

## Purpose
TBD - created by archiving change small-rna-seq-analysis-gao-pal-groups. Update Purpose after archive.
## Requirements
### Requirement: 基因计数矩阵生成
系统应基于比对结果生成基因表达计数矩阵。

#### Scenario: 基因注释文件处理
- **WHEN** 输入GTF格式基因注释文件
- **THEN** 系统解析注释信息，提取small RNA基因特征（如miRNA、piRNA、snoRNA等）

#### Scenario: 特征计数
- **WHEN** 输入排序后的BAM文件和基因注释
- **THEN** 系统使用featureCounts对每个基因进行计数，统计唯一比对到每个基因的reads数

#### Scenario: 计数矩阵生成
- **WHEN** 所有样本计数完成
- **THEN** 系统生成基因表达计数矩阵，行为基因，列为样本，值为原始reads计数

### Requirement: 表达量标准化
系统应对原始计数进行标准化处理，消除技术变异。

#### Scenario: 文库大小标准化
- **WHEN** 输入原始计数矩阵
- **THEN** 系统计算每个样本的总计数，进行文库大小标准化

#### Scenario: 方差稳定变换
- **WHEN** 标准化后数据
- **THEN** 系统应用DESeq2的方差稳定变换，减少高表达基因的方差膨胀

### Requirement: 差异表达分析
系统应识别GAO组和PAL组间的差异表达small RNA。

#### Scenario: DESeq2差异分析
- **WHEN** 输入计数矩阵和实验设计信息
- **THEN** 系统使用DESeq2进行负二项分布检验，计算每个基因的log2 fold change、p-value和adjusted p-value

#### Scenario: 多重检验校正
- **WHEN** 进行多个基因的假设检验
- **THEN** 系统应用Benjamini-Hochberg方法进行错误发现率校正，生成adjusted p-value

#### Scenario: 差异基因筛选
- **WHEN** 分析结果生成
- **THEN** 系统筛选显著差异表达基因，默认阈值：|log2FC| > 1且adjusted p-value < 0.05

### Requirement: 结果可视化和解释
系统应生成差异表达分析的可视化图表和结果报告。

#### Scenario: 火山图生成
- **WHEN** 差异分析完成
- **THEN** 系统生成火山图，x轴为log2 fold change，y轴为-log10(p-value)，显著基因高亮显示

#### Scenario: 热图生成
- **WHEN** 筛选出差异表达基因
- **THEN** 系统生成表达量热图，显示差异基因在样本间的表达模式

#### Scenario: MA图生成
- **WHEN** 分析完成
- **THEN** 系统生成MA图，显示基因平均表达量与fold change的关系

#### Scenario: 结果表格导出
- **WHEN** 差异分析完成
- **THEN** 系统导出完整结果表格，包含基因ID、基础表达量、log2FC、p-value、adjusted p-value等信息

