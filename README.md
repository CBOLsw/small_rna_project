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

## 完整环境配置流程（推荐WSL2/Linux）

### 方案一：一键安装（快速推荐）

```bash
# 1. 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 2. 给所有脚本添加执行权限
chmod +x setup_complete.sh install_bioc_packages.sh

# 3. 运行一键安装（预计30-45分钟）
./setup_complete.sh
```

### 方案二：手动分步安装（用于调试）

#### 第一步：配置镜像源
```bash
# 配置conda镜像源（清华源）
cat > ~/.condarc << 'EOF'
channels:
  - defaults
show_channel_urls: true
default_channels:
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r
  - https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2
custom_channels:
  conda-forge: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  bioconda: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  menpo: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  pytorch: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
  simpleitk: https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud
EOF

# 配置R镜像源
mkdir -p ~/.R
cat > ~/.Rprofile << 'EOFR'
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
message("已配置使用清华镜像源:")
message("  - Bioconductor: https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
message("  - CRAN: https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
EOFR
```

#### 第二步：创建并激活环境
```bash
# 安装mamba（如果没有）
conda install -y -c conda-forge mamba

# 创建conda环境
mamba env create -f envs/small_rna_analysis.yaml

# 激活环境
conda activate small_rna_analysis
```

#### 第三步：安装Bioconductor数据包
```bash
# 运行专门的Bioconductor包安装脚本（解决GenomeInfoDbData下载慢问题）
./install_bioc_packages.sh
```

#### 第四步：安装系统依赖（WSL2）
```bash
sudo apt-get update && sudo apt-get install -y fastqc samtools bowtie2 trimmomatic
```

### 验证安装

```bash
# 检查环境是否激活
conda env list
which python
python --version
which R
R --version

# 检查关键工具
which fastqc
which bowtie2
which samtools
which meme
which snakemake

# 检查R包
R -e "library(DESeq2); print('DESeq2 loaded')"
R -e "library(GenomeInfoDbData); print('GenomeInfoDbData loaded')"

# 检查Python包
python -c "import Bio; print('Biopython loaded')"
python -c "import pysam; print('pysam loaded')"
python -c "import pandas as pd; print('pandas loaded')"
```

### 运行分析流程

#### 步骤一：构建Bowtie2索引

```bash
python scripts/alignment/build_bowtie2_index.py \
  --genome references/hg38.fa \
  --output references/bowtie2_index/hg38 \
  --threads 4
```

#### 步骤二：运行完整流程

```bash
# 使用Python脚本（推荐）
python scripts/run_pipeline.py --config config/config.yaml --cores 8

# 或使用Snakemake
snakemake --cores 8 --configfile config/config.yaml
```

### 模块级运行（可选）

```bash
# 仅运行质量控制
python scripts/run_pipeline.py --config config/config.yaml --module qc

# 仅运行序列比对
python scripts/run_pipeline.py --config config/config.yaml --module alignment

# 仅运行差异表达分析
python scripts/run_pipeline.py --config config/config.yaml --module de

# 查看流程状态
python scripts/run_pipeline.py --config config/config.yaml --status
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

## 常见问题解决

### 问题1：Bioconductor数据包下载慢或卡住

**症状**：安装GenomeInfoDbData等大包时，下载速度极慢或直接卡住

**原因**：bioconductor-data-packages的安装脚本使用硬编码的官方URL，绕过了镜像源配置

**解决方案**：
```bash
# 方法一：使用专门的Bioconductor安装脚本（推荐）
conda activate small_rna_analysis
./install_bioc_packages.sh

# 方法二：直接在R中使用清华镜像安装
R -e "options(BioC_mirror='https://mirrors.tuna.tsinghua.edu.cn/bioconductor'); if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('GenomeInfoDbData', ask=FALSE, force=TRUE)"

# 方法三：手动下载并安装
cd /tmp
wget https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.18/data/annotation/src/contrib/GenomeInfoDbData_1.2.11.tar.gz
R -e "install.packages('/tmp/GenomeInfoDbData_1.2.11.tar.gz', repos = NULL, type = 'source')"
```

### 问题2：conda环境未激活

**症状**：运行脚本时提示找不到命令

**解决方案**：
```bash
# 检查环境列表
conda env list

# 激活环境
conda activate small_rna_analysis

# 如果激活失败，检查是否已初始化
conda init bash
```

### 问题3：系统依赖未安装

**症状**：fastqc、bowtie2等工具找不到

**解决方案**：
```bash
# 在WSL2中安装
sudo apt-get update && sudo apt-get install -y fastqc samtools bowtie2 trimmomatic
```

### 问题4：内存不足

**症状**：比对或DESeq2分析时崩溃

**解决方案**：
- 减少使用的核心数：`--cores 4` 而不是 `--cores 8`
- 关闭其他内存占用大的程序
- 在WSL2中调整内存限制：编辑`.wslconfig`文件

## 文档索引

- `USAGE_INSTRUCTIONS.md` - **完整使用说明**（必读）
- `setup_complete.sh` - **一键环境安装脚本**（推荐）
- `install_bioc_packages.sh` - **Bioconductor数据包快速安装脚本**
- `docs/environment_setup.md` - 完整环境配置详解
- `docs/wsl_setup_guide.md` - WSL2详细配置指南
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
