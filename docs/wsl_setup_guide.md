# WSL环境详细配置指南

## 第一步：启动WSL并检查状态

在PowerShell或命令提示符中运行：
```powershell
wsl
```

如果遇到网络问题（localhost相关错误），尝试：
```powershell
wsl --shutdown
wsl
```

## 第二步：在WSL中安装Miniconda

进入WSL后，执行以下命令：

```bash
# 下载Miniconda（使用清华镜像源，速度更快）
cd ~
wget https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-latest-Linux-x86_64.sh

# 安装Miniconda（静默安装）
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3

# 初始化conda
~/miniconda3/bin/conda init bash

# 重新加载shell配置
source ~/.bashrc

# 验证安装
conda --version
```

## 第三步：创建项目环境

```bash
# 进入项目目录
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 创建conda环境（这需要10-20分钟）
conda env create -f envs/small_rna_analysis.yaml

# 激活环境
conda activate small_rna_analysis

# 验证环境
python --version
which bowtie2
which samtools
which fastqc
```

## 第四步：运行分析流程

```bash
# 确保在项目目录中
cd /mnt/c/Users/24584/PycharmProjects/small_rna_project

# 激活环境
conda activate small_rna_analysis

# 查看Snakemake流程（试运行）
snakemake -n --configfile config/config.yaml

# 运行完整流程（使用8个CPU核心）
snakemake --cores 8 --configfile config/config.yaml
```

## 常见问题

### 1. WSL网络问题
如果WSL无法连接网络：
```powershell
# 在Windows PowerShell中
wsl --shutdown
# 然后重新启动WSL
wsl
```

### 2. conda下载慢
配置清华镜像源：
```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
```

### 3. 文件权限问题
在WSL中访问Windows文件时可能遇到权限问题：
```bash
# 修改文件权限
chmod +x scripts/*.py
chmod +x workflow/Snakefile
```

## 当前项目状态

- Windows环境：已创建 `small_rna_analysis_windows`（有限功能）
- WSL环境：需要手动按上述步骤配置
- 参考基因组：已下载到 `references/`
- 分析脚本：已就绪

## 推荐操作流程

1. 在Windows PowerShell中运行 `wsl` 启动Linux子系统
2. 按照本指南第二步安装Miniconda
3. 创建项目环境
4. 开始分析！
