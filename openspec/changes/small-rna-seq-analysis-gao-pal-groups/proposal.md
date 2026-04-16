## Why

本研究旨在对GAO组和PAL组的HeLa细胞small RNA测序数据进行系统分析，揭示两组间small RNA表达差异及其调控机制。通过完整的生物信息学分析流程，包括数据质控、序列比对、差异表达分析和de novo motif发现，可以识别两组间的关键差异表达small RNA，并探索其潜在的序列motif特征，为后续功能验证提供重要线索。

## What Changes

本项目将构建一个完整的small RNA测序数据分析流程：

- **数据质量评估**：对原始fastq文件进行质量评估，包括测序质量分布、接头污染检测、长度分布分析
- **序列比对**：将clean reads比对到参考基因组，统计比对率、唯一比对率等指标
- **基因计数**：基于比对结果对small RNA进行计数，生成基因表达矩阵
- **差异表达分析**：比较GAO组和PAL组的small RNA表达差异，识别显著上调和下调的small RNA
- **de novo motif发现**：分别在GAO组和PAL组差异表达的small RNA序列中寻找富集的序列motif
- **结果可视化**：生成分析报告，包括质量控制图、差异表达火山图、热图、motif序列logo等
- **可重复分析框架**：建立模块化的分析脚本，确保分析过程的可重复性

## Capabilities

### New Capabilities

- **data-quality-control**: 对small RNA测序数据进行全面的质量控制，包括FastQC分析、序列质量修剪、接头去除
- **sequence-alignment**: 将质控后的reads比对到参考基因组，支持Bowtie2比对工具，处理small RNA短序列特性
- **differential-expression-analysis**: 基于基因计数矩阵进行统计检验，识别GAO组与PAL组间的差异表达small RNA
- **motif-discovery**: 对差异表达small RNA的序列进行de novo motif分析，发现富集的序列模式
- **analysis-pipeline-orchestration**: 协调整个分析流程，自动化执行各个分析步骤并生成最终报告

### Modified Capabilities

<!-- 本项目为全新分析流程，不修改现有功能 -->

## Impact

- **数据存储**：需要约20-30GB的临时存储空间用于中间文件处理
- **计算资源**：比对和分析步骤需要中等计算资源，建议使用多核CPU
- **软件依赖**：需要安装FastQC、Trimmomatic、Bowtie2、featureCounts、DESeq2、MEME Suite等生物信息学工具
- **参考数据**：需要人类参考基因组（如hg38）和相应的基因注释文件
- **结果输出**：将生成质量控制报告、比对统计、差异表达基因列表、motif分析结果等