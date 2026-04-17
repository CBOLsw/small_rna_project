# Small RNA测序分析项目用户指南

## 快速开始

### 1. 环境准备

首先，确保您的系统满足以下要求：
- 操作系统：Linux / macOS / Windows（推荐使用WSL）
- 内存：至少16GB RAM（推荐32GB）
- 存储：至少100GB可用空间
- Python 3.8+
- R 4.0+

### 2. 安装conda环境

```bash
# 从项目根目录开始
cd small_rna_project

# 创建conda环境（约需要10-20分钟）
conda env create -f envs/small_rna_analysis.yaml

# 激活环境
conda activate small_rna_analysis
```

### 3. 准备数据

#### 3.1 下载参考基因组数据

```bash
# 创建参考基因组目录
mkdir -p references

# 下载hg38参考基因组（来自UCSC）
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
mv hg38.fa references/

# 下载hg38基因注释
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
gzip -d hg38.refGene.gtf.gz
mv hg38.refGene.gtf references/hg38.gtf
```

#### 3.2 准备测序数据

将您的测序数据文件放置在 `data/raw_fastq/` 目录中：

```
data/raw_fastq/
├── GAO_1_R1.fastq.gz
├── GAO_2_R1.fastq.gz
├── GAO_3_R1.fastq.gz
├── PAL_1_R1.fastq.gz
├── PAL_2_R1.fastq.gz
└── PAL_3_R1.fastq.gz
```

#### 3.3 更新样本信息

编辑 `data/metadata/sample_info.csv` 文件，记录样本信息：

```csv
sample_id,group,replicate,fastq_path
GAO_1,GAO,1,data/raw_fastq/GAO_1_R1.fastq.gz
GAO_2,GAO,2,data/raw_fastq/GAO_2_R1.fastq.gz
GAO_3,GAO,3,data/raw_fastq/GAO_3_R1.fastq.gz
PAL_1,PAL,1,data/raw_fastq/PAL_1_R1.fastq.gz
PAL_2,PAL,2,data/raw_fastq/PAL_2_R1.fastq.gz
PAL_3,PAL,3,data/raw_fastq/PAL_3_R1.fastq.gz
```

### 4. 运行分析流程

#### 4.1 构建参考基因组索引

```bash
# 构建Bowtie2索引（约需要30分钟）
python scripts/alignment/build_bowtie2_index.py \
    --genome references/hg38.fa \
    --output references/bowtie2_index \
    --prefix hg38
```

#### 4.2 运行完整分析流程

```bash
# 使用8个核心运行完整流程（约需要2-3小时）
python scripts/run_pipeline.py --config config/config.yaml --cores 8
```

### 5. 查看结果

结果将保存在以下位置：
- `results/` - 分析结果
- `reports/` - 分析报告
- `logs/` - 运行日志

```bash
# 查看生成的报告（浏览器中打开）
start reports/report_*.html
```

## 详细使用方法

### 自定义配置

编辑 `config/config.yaml` 文件来自定义参数：

```yaml
# 修改质量控制参数
quality_control:
  trimmomatic:
    leading: 5          # 更严格的质量阈值
    minlen: 20          # 提高最小序列长度要求

# 修改差异表达分析阈值
differential_expression:
  deseq2:
    padj_threshold: 0.01    # 更严格的统计显著性
    log2fc_threshold: 1.5    # 要求更大的变化幅度
```

### 模块级运行

#### 仅运行质量控制

```bash
python scripts/run_pipeline.py --config config/config.yaml --module qc
```

#### 仅运行序列比对

```bash
python scripts/run_pipeline.py --config config/config.yaml --module alignment
```

#### 仅运行差异表达分析

```bash
python scripts/run_pipeline.py --config config/config.yaml --module de
```

### 检查流程状态

```bash
# 查看当前流程运行状态
python scripts/run_pipeline.py --config config/config.yaml --status
```

### 从失败处恢复

如果流程在执行过程中失败，可以从失败处继续：

```bash
# 从上次失败的步骤恢复
python scripts/run_pipeline.py --config config/config.yaml --resume
```

### 查看可用模块

```bash
# 列出所有可单独运行的模块
python scripts/run_pipeline.py --list-modules
```

## 输出结果说明

### 1. 质量控制结果

- **FastQC报告**: HTML格式的详细质量评估报告
- **QC摘要**: 所有样本的质量统计汇总

### 2. 序列比对结果

- **BAM文件**: 排序后的比对结果
- **比对统计**: 每个样本的比对率、唯一比对率等信息

### 3. 基因计数结果

- **表达矩阵**: 所有样本的基因表达量矩阵
- **统计摘要**: 表达量分布统计

### 4. 差异表达分析结果

- **差异基因完整列表**: 包含所有基因的统计结果
- **筛选后的差异基因**: 显著差异表达的基因
- **可视化图表**: 火山图、热图等

### 5. Motif分析结果

- **发现的motif**: 序列motif及其统计信息
- **数据库比对结果**: 与已知motif数据库的比对结果
- **可视化**: 序列logo图等

## 常见问题

### Q: 如何添加新的样本？

A: 编辑 `data/metadata/sample_info.csv` 文件，添加新的样本信息，然后重新运行流程。Snakemake会自动检测新样本并只运行新样本的分析。

### Q: 如何调整Motif分析参数？

A: 修改 `config/config.yaml` 中的 `motif_analysis` 部分：

```yaml
motif_analysis:
  meme:
    min_width: 4       # 更短的最小motif宽度
    max_width: 16      # 更长的最大motif宽度
    max_motifs: 20     # 发现更多的motif
```

### Q: 如何在服务器上运行？

A: 使用以下命令在后台运行：

```bash
# 使用nohup在后台运行
nohup python scripts/run_pipeline.py --config config/config.yaml --cores 16 > run.log 2>&1 &

# 或使用screen
screen -S small_rna
python scripts/run_pipeline.py --config config/config.yaml --cores 16
# 按 Ctrl+A 然后 D 分离
```

### Q: 分析需要多长时间？

A: 取决于数据量和计算资源，典型运行时间：
- 质量控制：~30分钟
- 序列比对：~1小时
- 差异表达：~30分钟
- Motif分析：~1小时
- 完整流程：~2-3小时（使用8核CPU）

### Q: 如何释放磁盘空间？

A: 可以删除中间结果来释放空间：

```bash
# 删除处理中的中间文件（保留最终结果）
rm -rf data/processed/*.fastq.gz
```

### Q: 如何导出结果用于论文？

A: 使用结果整合脚本：

```bash
python scripts/results_integration.py --config config/config.yaml
```

这将生成包含所有结果的综合报告，包含高质量的图表和统计数据。

## 进阶使用

### 1. 使用SLURM集群提交作业

```bash
# 创建SLURM提交脚本
cat > submit_job.sh << 'EOF'
#!/bin/bash
#SBATCH --job-name=small_rna
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

conda activate small_rna_analysis
python scripts/run_pipeline.py --config config/config.yaml --cores 16
EOF

# 提交作业
sbatch submit_job.sh
```

### 2. 使用自定义参考基因组

如果需要使用其他参考基因组（如hg19），需要：

1. 下载对应的fasta和注释文件
2. 修改 `config/config.yaml` 中的路径
3. 重新构建Bowtie2索引

### 3. 分析其他物种

修改 `config/config.yaml`，更新参考基因组和注释文件的路径即可。

## 技术支持

如果您遇到问题：

1. 查看日志文件：`logs/` 目录
2. 检查系统资源：磁盘空间、内存使用情况
3. 验证输入数据：确保所有文件路径正确
4. 查阅技术文档：`docs/technical_manual.md`

---

**最后更新**: 2024年4月17日

**项目版本**: 1.0
