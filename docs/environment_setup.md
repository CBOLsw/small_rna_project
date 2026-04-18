# 项目环境配置指南

## 概述

本项目是一个small RNA测序数据分析项目，使用Snakemake进行流程管理，Conda进行环境管理。

## 环境配置方案

### 方案一：Windows系统（推荐用于数据分析和可视化）

```bash
# 在Windows PowerShell中执行
conda env create -f envs/small_rna_analysis_windows.yaml
conda activate small_rna_analysis_windows
```

### 方案二：WSL2/Linux系统（推荐用于完整流程运行）

**方法1：一键安装脚本（推荐，使用mamba加速）**
```bash
# 在WSL/Linux终端中执行
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project
chmod +x setup_complete.sh
./setup_complete.sh
```

**方法2：手动配置（使用mamba加速）**
```bash
# 在WSL/Linux终端中执行
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p ~/miniconda3
export PATH="$HOME/miniconda3/bin:$PATH"
conda init bash
exec $SHELL

# 配置清华镜像源
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
EOF

# 安装mamba（快速包管理器）
conda install -y -c conda-forge mamba

# 使用mamba创建环境（速度是conda的10倍）
mamba env create -f envs/small_rna_analysis.yaml
conda activate small_rna_analysis
```

## 生物信息学工具安装

### Windows系统（WSL2/Linux中已包含）

```bash
# 使用apt-get安装（WSL/Linux）
sudo apt-get update && sudo apt-get install -y fastqc samtools bowtie2 trimmomatic
```

### 使用Conda安装

```bash
# 使用Conda安装所有工具
conda install -c bioconda fastqc samtools bowtie2 trimmomatic
```

## 参考数据下载

```bash
python download_references.py
```

## 项目运行

```bash
# 查看流程概览
snakemake -n --configfile config/config.yaml

# 运行完整流程
snakemake --cores 8 --configfile config/config.yaml
```

## 常见问题

### Conda安装速度慢

```bash
# 配置清华镜像源
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/msys2/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
```

### 包安装失败

```bash
# 使用pip安装
pip install pandas numpy matplotlib seaborn jupyter biopython pysam
```

## 环境重置

```bash
# 删除旧环境
conda env remove -n small_rna_analysis -y

# 重新创建
conda env create -f envs/small_rna_analysis.yaml
```

## 系统要求

- Windows 10/11 或 Linux系统
- 至少8GB内存
- 至少20GB磁盘空间
- 网络连接（用于下载工具和数据）
