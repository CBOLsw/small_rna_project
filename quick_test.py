#!/usr/bin/env python3
"""
快速测试脚本 - 用于验证分析管道的基本功能
"""

import os
import sys
from pathlib import Path
import gzip
import pandas as pd

# 添加项目路径
sys.path.append(str(Path(__file__).parent))

print("=" * 60)
print("Small RNA分析管道快速测试")
print("=" * 60)

print("\n1. 测试数据读取...")
data_dir = Path("data/raw_fastq/fastq_files")
if data_dir.exists():
    fastq_files = list(data_dir.glob("*.fastq.gz"))
    print(f"   [OK] 找到 {len(fastq_files)} 个fastq.gz文件")

    # 读取第一个文件的前几行
    if fastq_files:
        first_file = fastq_files[0]
        try:
            with gzip.open(first_file, 'rt') as f:
                lines_read = 0
                for line in f:
                    lines_read += 1
                    if lines_read >= 40:  # 读取前10个reads
                        break
            print(f"   [OK] 成功读取 {lines_read} 行数据")
        except Exception as e:
            print(f"   [FAIL] 数据读取失败: {e}")

print("\n2. 测试Python模块导入...")
modules_to_test = [
    ("QC模块", "scripts.qc.fastqc_analysis"),
    ("QC模块", "scripts.qc.trim_fastq"),
    ("比对模块", "scripts.alignment.run_bowtie2"),
    ("表达模块", "scripts.expression.count_features"),
    ("Motif模块", "scripts.motif.extract_sequences"),
]

all_modules_ok = True
for module_desc, module_name in modules_to_test:
    try:
        __import__(module_name.replace('.', '/').replace('/', '.'))
        print(f"   [OK] {module_desc}导入成功")
    except Exception as e:
        print(f"   [FAIL] {module_desc}导入失败: {e}")
        all_modules_ok = False

print("\n3. 检查样本信息...")
sample_file = Path("data/metadata/sample_info.csv")
if sample_file.exists():
    try:
        df = pd.read_csv(sample_file)
        print(f"   [OK] 样本信息文件包含 {len(df)} 个样本")
        print(f"   样本列表:")
        for _, row in df.iterrows():
            print(f"     - {row['sample']} ({row['group']})")
    except Exception as e:
        print(f"   [FAIL] 样本信息读取失败: {e}")

print("\n4. 检查配置和流程文件...")
config_files = [
    ("配置文件", "config/config.yaml"),
    ("流程文件", "workflow/Snakefile"),
    ("QC脚本", "scripts/qc/test_qc.py"),
]

for config_desc, config_path in config_files:
    if Path(config_path).exists():
        print(f"   [OK] {config_desc}找到: {config_path}")
    else:
        print(f"   [FAIL] {config_desc}未找到")

print("\n5. 检查QC测试脚本...")
qc_test_script = Path("scripts/qc/test_qc.py")
if qc_test_script.exists():
    print(f"   [OK] QC测试脚本找到")
    print(f"\n提示: 您可以运行以下命令来测试QC模块:")
    print(f"  cd scripts/qc/")
    print(f"  python test_qc.py")

print("\n" + "=" * 60)
print("快速测试完成")
print("=" * 60)

print("\n下一步:")
print("1. 下载并准备参考基因组数据")
print("2. 构建Bowtie2索引")
print("3. 运行完整分析流程")
print("\n详细执行计划请查看: complete_pipeline_execution_plan.md")
