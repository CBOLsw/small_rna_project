# Small RNA测序分析项目 Dockerfile
# 使用mambaorg/micromamba作为基础镜像，方便conda环境管理
FROM mambaorg/micromamba:jammy

LABEL maintainer="CBOLS"
LABEL description="Small RNA sequencing analysis pipeline"
LABEL version="1.0"

# 设置工作目录
WORKDIR /small_rna_project

# 切换到root用户安装系统依赖
USER root

# 安装系统级依赖
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    curl \
    git \
    unzip \
    bzip2 \
    ca-certificates \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    build-essential \
    procps \
    && rm -rf /var/lib/apt/lists/*

# 复制conda环境配置文件
COPY --chown=$MAMBA_USER:$MAMBA_USER envs/small_rna_analysis.yaml /tmp/env.yaml

# 创建conda环境
RUN micromamba create -y -f /tmp/env.yaml \
    && micromamba clean --all --yes

# 配置环境变量
ENV PATH /opt/conda/envs/small_rna_analysis/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE 1

# 切换回mamba用户
USER $MAMBA_USER

# 复制项目文件
COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# 创建必要的目录
RUN mkdir -p data/raw_fastq \
             data/processed \
             data/metadata \
             references \
             results \
             logs \
             reports

# 设置执行权限
RUN chmod +x scripts/*.py \
    && chmod +x scripts/qc/*.py \
    && chmod +x scripts/alignment/*.py \
    && chmod +x scripts/expression/*.py \
    && chmod +x scripts/motif/*.py

# 设置工作目录
WORKDIR /small_rna_project

# 测试环境是否正确设置
RUN micromamba run -n small_rna_analysis python --version \
    && micromamba run -n small_rna_analysis R --version

# 默认命令：显示帮助信息
CMD ["python", "scripts/run_pipeline.py", "--help"]
