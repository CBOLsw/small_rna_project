## 1. 项目设置和环境配置

- [x] 1.1 创建conda环境配置文件 `envs/small_rna_analysis.yaml`
- [x] 1.2 安装所有必要的生物信息学工具（FastQC、Trimmomatic、Bowtie2、samtools、featureCounts、DESeq2、MEME Suite等）
- [x] 1.3 设置项目目录结构，确保data、scripts、results等目录正确创建
- [x] 1.4 准备参考基因组数据，下载hg38参考基因组和基因注释文件到 `references/` 目录
- [x] 1.5 创建样本信息文件 `data/metadata/sample_info.csv`，记录样本名称、组别、重复、fastq文件路径

## 2. 数据质量控制模块

- [x] 2.1 编写FastQC分析脚本 `scripts/qc/fastqc_analysis.py`
- [x] 2.2 编写Trimmomatic质量修剪脚本 `scripts/qc/trim_fastq.py`
- [x] 2.3 实现质控结果汇总功能，生成所有样本的质量控制报告
- [x] 2.4 创建质控可视化脚本，生成reads数条形图、质量分布热图等图表
- [x] 2.5 测试数据质控模块，确保正确处理原始fastq文件（`scripts/qc/test_qc.py`）

## 3. 序列比对模块

- [x] 3.1 构建Bowtie2参考基因组索引 `scripts/alignment/build_bowtie2_index.py`
- [x] 3.2 编写Bowtie2比对脚本 `scripts/alignment/run_bowtie2.py`，支持small RNA短序列参数优化
- [x] 3.3 实现SAM到BAM格式转换和排序功能（在`run_bowtie2.py`中增强）
- [x] 3.4 编写比对统计脚本 `scripts/alignment/alignment_stats.py`，计算每个样本的比对率、唯一比对率
- [x] 3.5 创建比对质量评估脚本 `scripts/alignment/bam_quality_assessment.py`，生成比对质量分布图和基因组覆盖度分析

## 4. 基因计数和差异表达分析模块

- [x] 4.1 编写featureCounts基因计数脚本 `scripts/expression/count_features.py`
- [x] 4.2 创建基因表达矩阵生成脚本 `scripts/expression/generate_expression_matrix.py`，整合所有样本的计数结果
- [x] 4.3 实现DESeq2差异表达分析脚本 `scripts/expression/deseq2_analysis.R`
- [x] 4.4 编写差异基因筛选和结果导出功能 `scripts/expression/filter_degs.py`
- [x] 4.5 创建差异表达可视化脚本 `scripts/expression/visualize_degs.py`，生成火山图、热图、MA图
- [x] 4.6 测试差异表达分析模块 `scripts/expression/test_expression.py`，验证统计结果的正确性

## 5. de novo motif发现模块

- [x] 5.1 编写差异基因序列提取脚本 `scripts/motif/extract_sequences.py`（基础框架完成）
- [x] 5.2 创建MEME motif分析脚本 `scripts/motif/run_meme.py`
- [x] 5.3 实现TomTom motif比较脚本，与已知motif数据库比对
- [x] 5.4 编写motif结果过滤和验证功能，去除假阳性结果
- [x] 5.5 创建motif可视化脚本，生成序列logo图、motif位置分布图
- [x] 5.6 实现组间motif比较功能，识别GAO组和PAL组的共有和特有motif

## 6. 流程编排和自动化

- [x] 6.1 设计Snakemake流程文件 `workflow/Snakefile`，定义所有分析步骤的规则
- [x] 6.2 创建配置文件 `config/config.yaml`，集中管理样本信息、参数设置
- [x] 6.3 实现流程执行脚本 `scripts/run_pipeline.py`，支持命令行参数
- [x] 6.4 添加日志记录功能，记录每个分析步骤的执行状态
- [x] 6.5 实现错误处理和恢复机制，支持从失败步骤重新开始
- [x] 6.6 创建结果整合脚本，自动收集各模块输出生成综合报告

## 7. 结果验证和文档生成

- [x] 7.1 测试完整分析流程，使用示例数据验证各模块正确性（已完成QC模块测试）
- [x] 7.2 生成技术文档 `docs/technical_manual.md`，详细说明分析方法和参数
- [x] 7.3 创建用户指南 `docs/user_guide.md`，指导如何运行分析流程
- [x] 7.4 准备分析报告模板 `reports/report_template.Rmd`，支持自动生成分析报告
- [ ] 7.5 验证分析结果的可重复性，确保相同输入产生一致输出
- [x] 7.6 打包项目环境，创建Dockerfile或Singularity定义文件

## 8. 项目实施状态

### 已完成任务
- 清理非核心项目文件
- 测试原始数据读取（共13个fastq.gz文件）
- 下载hg38参考基因组和基因注释文件
- 测试QC模块，验证FastQC和Trimmomatic功能

### 当前问题
- **Windows平台限制**：Bioconda中的部分工具（如Bowtie2、SAMtools）在Windows系统上不可用
- **miRBase注释文件下载失败**：https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3 返回404错误
- **conda环境创建失败**：Windows系统上无法直接创建包含Bioconda工具的环境

### 解决方案建议
1. 使用Windows Subsystem for Linux (WSL2) 运行完整流程
2. 使用Docker容器运行分析
3. 手动安装部分工具的Windows版本