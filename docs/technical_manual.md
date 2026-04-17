# Small RNA测序分析技术手册

## 概述

本项目实现了一个完整的small RNA测序数据分析流程，专门针对GAO组和PAL组的HeLa细胞small RNA测序数据进行分析。流程采用模块化设计，确保各个分析步骤的独立性和可重复性。

## 分析流程

### 整体架构

```
原始数据 → 数据质控 → 序列比对 → 基因计数 → 差异表达分析 → motif发现 → 结果可视化
```

### 模块设计

1. **数据质控模块** - 对原始测序数据进行质量评估和过滤
2. **序列比对模块** - 将clean reads比对到参考基因组
3. **基因计数模块** - 统计基因表达量
4. **差异表达分析模块** - 识别差异表达的small RNA
5. **motif发现模块** - 寻找差异表达基因中的序列motif
6. **结果整合与报告模块** - 生成综合分析报告

## 工具选择

### 质量控制

- **FastQC** - 质量评估工具
- **Trimmomatic** - 质量修剪和接头去除工具

### 序列比对

- **Bowtie2** - 短序列比对工具，对small RNA优化
- **Samtools** - SAM/BAM格式处理工具

### 基因计数

- **featureCounts** - 高效的基因计数工具

### 差异表达分析

- **DESeq2** - 专门为RNA-seq数据设计的差异表达分析工具

### Motif分析

- **MEME Suite** - 集成的motif发现和分析工具套件
- **JASPAR数据库** - 脊椎动物转录因子结合位点数据库

### 流程管理

- **Snakemake** - 声明式流程管理工具
- **Python** - 主要脚本语言

## 配置参数

### 质量控制参数

```yaml
quality_control:
  fastqc:
    threads: 4
  trimmomatic:
    threads: 4
    leading: 3          # 切除开头碱基质量小于3的序列
    trailing: 3         # 切除结尾碱基质量小于3的序列
    slidingwindow: "4:15" # 滑动窗口质量控制（4个碱基平均质量<15）
    minlen: 18          # 最小序列长度
```

### 序列比对参数

```yaml
alignment:
  bowtie2:
    threads: 8
    preset: "very-sensitive-local"
    seed_length: 15
    max_mismatches: 1
    k: 10
```

### Motif分析参数

```yaml
motif_analysis:
  meme:
    min_width: 6
    max_width: 12
    max_motifs: 10
    evalue_threshold: 1e-4
    threads: 8
```

## 输入数据格式

### 原始数据

- **格式**: fastq.gz（压缩的FASTQ格式）
- **位置**: `data/raw_fastq/` 目录
- **文件名格式**: `<sample_id>_R1.fastq.gz`

### 参考数据

- **人类基因组**: hg38.fa（UCSC版本）
- **基因注释**: hg38.gtf（RefSeq注释）
- **Bowtie2索引**: references/bowtie2_index/ 目录

## 输出结果格式

### 数据质控结果

```
results/qc/
├── <sample_id>_fastqc.html - 质量评估HTML报告
├── <sample_id>_fastqc.zip - 质量评估原始数据
└── qc_summary.csv - 所有样本的质量统计
```

### 序列比对结果

```
results/alignment/
├── <sample_id>.sorted.bam - 比对结果（BAM格式）
├── <sample_id>_alignment_stats.csv - 比对统计信息
└── alignment_summary.csv - 综合比对统计
```

### 基因计数结果

```
results/counts/
├── gene_counts.csv - 基因表达矩阵
└── counts_summary.csv - 表达量统计
```

### 差异表达分析结果

```
results/differential_expression/
├── deseq2_results.csv - 完整差异表达结果
├── filtered_degs.csv - 差异基因筛选结果
├── volcano_plot.png - 火山图
└── heatmap.png - 热图
```

### Motif分析结果

```
results/motif_analysis/
├── gene_sequences.fasta - 差异基因序列
├── filtered_motifs.csv - 过滤后的motif
├── meme_results/ - MEME软件输出
├── tomtom_results/ - TomTom比较结果
└── visualization/ - 可视化图表
```

## 使用方法

### 环境准备

```bash
# 创建conda环境
conda env create -f envs/small_rna_analysis.yaml

# 激活环境
conda activate small_rna_analysis
```

### 流程运行

```bash
# 运行完整分析流程
python scripts/run_pipeline.py --config config/config.yaml --cores 8

# 运行特定模块
python scripts/run_pipeline.py --config config/config.yaml --module qc

# 查看运行状态
python scripts/run_pipeline.py --config config/config.yaml --status
```

### 测试运行

```bash
# 运行质量控制模块测试
python scripts/qc/test_qc.py

# 运行序列比对模块测试
python scripts/alignment/test_alignment.py

# 运行差异表达分析模块测试
python scripts/expression/test_expression.py
```

## 报告生成

### 自动报告

流程完成后会自动生成综合分析报告:
- HTML版本报告在 `reports/` 目录中
- 包含所有分析结果的可视化和统计

### 结果查看

```bash
# 查看报告（使用浏览器打开）
start reports/report_*.html
```

## 可重复性保障

### 环境配置

使用conda管理所有依赖包，确保环境一致性。

### 版本控制

所有分析脚本和配置文件使用Git进行版本控制。

### 流程日志

- 详细的执行日志在 `logs/` 目录中
- 包含每个步骤的开始/结束时间、状态和输出信息

## 故障排除

### 常见问题

1. **FastQC报告显示质量低**
   - 检查序列长度分布和碱基质量分布
   - 调整Trimmomatic参数

2. **比对率过低**
   - 检查基因组索引是否正确
   - 调整Bowtie2参数（如减少seed length）

3. **差异表达分析失败**
   - 检查基因计数矩阵是否包含所有样本
   - 验证样本信息元数据文件的格式

4. **Motif分析无结果**
   - 检查差异基因数量是否足够（至少需要10个基因）
   - 调整MEME参数（如最小宽度、E值阈值）

### 恢复机制

```bash
# 从失败的步骤恢复
python scripts/run_pipeline.py --config config/config.yaml --resume
```

## 性能优化

### 并行执行

```bash
# 使用更多核心运行
python scripts/run_pipeline.py --config config/config.yaml --cores 16
```

### 内存管理

- 确保系统有足够的内存（>16GB建议）
- 对于大文件处理，调整缓存参数

### 存储优化

- 使用压缩格式（如gzip）存储中间结果
- 定期清理不需要的临时文件

## 更新记录

### 版本1.0（2024年）

- 初始版本发布
- 支持GAO和PAL两个实验组的分析
- 包含完整的质量控制、比对、差异表达分析和motif发现

---

**项目链接**: [https://github.com/your-repo/small_rna_analysis_gao_pal](https://github.com/your-repo/small_rna_analysis_gao_pal)

**作者**: 生物信息学小组

**日期**: 2024年4月17日
