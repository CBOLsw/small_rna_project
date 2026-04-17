# 完整流程执行计划

## 概述
本计划详细描述了从参考基因组下载到完整分析执行的所有步骤，确保所有任务按正确顺序完成。

---

## 📋 待办任务列表

### 第一阶段：环境准备
- [x] ✅ 清理非核心项目文件
- [ ] 📥 下载hg38参考基因组（约3GB）
- [ ] 📥 下载hg38基因注释文件（约125MB）
- [ ] 📥 下载small RNA注释（约2MB）
- [ ] 🏗️ 构建Bowtie2基因组索引（约20GB，可能需要2-3小时）
- [ ] ✅ 创建conda环境

### 第二阶段：数据预处理
- [ ] 🧪 运行数据质量控制（QC）
- [ ] ✂️ 序列质量修剪
- [ ] 📊 生成质量控制摘要

### 第三阶段：序列比对
- [ ] 🎯 运行Bowtie2序列比对
- [ ] 📈 比对结果统计
- [ ] 🏷️ 读取计数（featureCounts）
- [ ] 📋 生成基因表达矩阵

### 第四阶段：差异表达分析
- [ ] 🔬 运行DESeq2差异表达分析
- [ ] 🎨 可视化差异表达结果
- [ ] 📝 筛选显著差异表达基因（FDR < 0.05，|log2FC| > 1）

### 第五阶段：Motif分析
- [ ] 🧬 提取差异表达基因的序列
- [ ] 🔍 运行MEME motif分析
- [ ] 📊 motif结果过滤
- [ ] 🔗 运行TomTom motif比较
- [ ] 🎨 motif可视化

### 第六阶段：结果整合
- [ ] 📚 收集和整理所有结果
- [ ] 📈 生成分析报告
- [ ] 📊 结果可视化整合

---

## 🚀 立即执行步骤

### 1. 快速开始 - 下载参考基因组
```bash
python download_references.py --threads 4 --force 2>&1 | tee download_references.log
```

### 2. 构建基因组索引
```bash
python scripts/alignment/build_bowtie2_index.py \
  --genome references/hg38.fa \
  --output references/bowtie2_index/hg38 \
  --threads 4 --small-rna 2>&1 | tee bowtie2_index_build.log
```

### 3. 快速运行QC模块（预览模式）
```bash
python scripts/run_pipeline.py \
  --config config/config.yaml \
  --module qc \
  --preview 2>&1 | tee qc_preview.log
```

---

## 📊 当前状态检查

### 项目结构完整性
✅ **所有脚本文件** - 完整  
✅ **所有配置文件** - 完整  
✅ **目录结构** - 已创建  
✅ **数据准备** - 已准备好13个fastq.gz文件  

### 已下载文件（检查）
```bash
ls -lh references/
```

### 待下载文件
- `hg38.fa` - 参考基因组
- `hg38.gtf` - 基因注释
- `hg38.mirbase.gff3` - small RNA注释

---

## 📈 预期资源使用情况

### 磁盘空间需求
- 参考基因组: ~3GB
- 索引文件: ~20GB
- 结果文件: ~10GB (取决于样本数量)
- **总需求: 约40GB**

### 内存使用
- 下载: ~2GB
- 索引构建: 8GB+
- 比对: 8GB+
- **建议配置: 16GB内存**

### CPU需求
- 建议: 8-16个CPU核心
- 最小: 4个CPU核心

---

## ⏱️ 时间预估（16核CPU）

1. 参考基因组下载: 30-60分钟
2. 索引构建: 2-3小时
3. 数据预处理: 30-60分钟
4. 序列比对: 1-2小时
5. 基因计数: 30分钟
6. 差异表达分析: 20分钟
7. Motif分析: 30-60分钟
8. 结果整合: 15分钟

**总时间: 6-8小时**

---

## 🔧 应急措施

### 如果下载失败
1. 检查网络连接
2. 使用VPN或镜像站点
3. 手动下载后放置到正确位置

### 如果索引构建失败
1. 增加系统资源
2. 使用预构建索引文件
3. 减少并发线程数

### 如果分析流程出错
1. 查看日志文件
2. 检查脚本输出
3. 使用预览模式调试

---

## 📞 支持与帮助

### 运行脚本帮助
```bash
python scripts/run_pipeline.py --help
python scripts/alignment/build_bowtie2_index.py --help
python download_references.py --help
```

### 查看状态
```bash
python scripts/run_pipeline.py --config config/config.yaml --status
```

---

## 📋 执行日志

### 下载过程跟踪
```bash
tail -f download_references.log
```

### 索引构建跟踪
```bash
tail -f bowtie2_index_build.log
```

### 流程执行跟踪
```bash
python scripts/run_pipeline.py --config config/config.yaml --status
```
