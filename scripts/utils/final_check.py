#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
最终项目状态检查脚本
"""
import sys
import yaml
from pathlib import Path
import subprocess

# 修复Windows编码问题
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

def check_config():
    """检查配置文件"""
    print("=== 检查配置文件 ===")
    config_path = Path("config/config.yaml")
    if not config_path.exists():
        print("❌ config.yaml 不存在")
        return False

    try:
        with open(config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)

        print("✅ config.yaml 解析成功")

        # 检查VAHTS配置
        if config.get('quality_control', {}).get('trimmomatic', {}).get('adapter_type') == 'vahts_small_rna_v2':
            print("✅ adapter_type 配置正确")
        else:
            print("❌ adapter_type 配置不正确")
            return False

        if config.get('quality_control', {}).get('trimmomatic', {}).get('adapter_file') == 'VAHTS-SmallRNA-V2.fa':
            print("✅ adapter_file 配置正确")
        else:
            print("❌ adapter_file 配置不正确")
            return False

    except Exception as e:
        print(f"❌ 解析 config.yaml 出错: {e}")
        return False

    return True


def check_vahts_file():
    """检查VAHTS接头序列文件"""
    print("\n=== 检查VAHTS接头序列文件 ===")
    vahts_path = Path("config/VAHTS-SmallRNA-V2.fa")
    if not vahts_path.exists():
        print("❌ VAHTS-SmallRNA-V2.fa 不存在")
        return False

    try:
        with open(vahts_path, 'r', encoding='utf-8') as f:
            content = f.read()

        if ">VAHTS_SmallRNA_V2_3prime" in content and ">VAHTS_SmallRNA_V2_5prime" in content:
            print("✅ VAHTS接头序列文件内容正确")
        else:
            print("❌ VAHTS接头序列文件内容不正确")
            return False

    except Exception as e:
        print(f"❌ 读取 VAHTS-SmallRNA-V2.fa 出错: {e}")
        return False

    return True


def check_trim_script():
    """检查trim_fastq.py脚本"""
    print("\n=== 检查trim_fastq.py脚本 ===")
    script_path = Path("scripts/qc/trim_fastq.py")
    if not script_path.exists():
        print("❌ trim_fastq.py 不存在")
        return False

    try:
        with open(script_path, 'r', encoding='utf-8') as f:
            content = f.read()

        if 'vahts_small_rna_v2' in content:
            print("✅ trim_fastq.py 包含VAHTS接头类型")
        else:
            print("❌ trim_fastq.py 缺少VAHTS接头类型")
            return False

        if 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC' in content:
            print("✅ VAHTS接头序列正确")
        else:
            print("❌ VAHTS接头序列不正确")
            return False

    except Exception as e:
        print(f"❌ 读取 trim_fastq.py 出错: {e}")
        return False

    return True


def check_project_structure():
    """检查项目目录结构"""
    print("\n=== 检查项目目录结构 ===")
    required_dirs = [
        "data", "data/raw_fastq", "data/processed", "data/metadata",
        "references", "references/bowtie2_index",
        "results", "results/qc", "results/alignment", "results/counts",
        "results/differential_expression", "results/motif_analysis",
        "logs", "reports", "scripts", "scripts/qc", "scripts/alignment",
        "scripts/expression", "scripts/motif", "scripts/setup", "scripts/utils",
        "workflow", "config", "envs", "openspec", "openspec/specs"
    ]

    all_exists = True
    for dir_path in required_dirs:
        if not Path(dir_path).exists():
            print(f"⚠️  目录不存在: {dir_path}")
            all_exists = False

    if all_exists:
        print("✅ 所有必需目录均存在")
    else:
        print("⚠️  部分目录不存在（可能需要用户手动创建）")

    return True


def check_git_ignore():
    """检查.gitignore文件"""
    print("\n=== 检查.gitignore配置 ===")
    gitignore_path = Path(".gitignore")
    if not gitignore_path.exists():
        print("❌ .gitignore 不存在")
        return False

    try:
        with open(gitignore_path, 'r', encoding='utf-8') as f:
            content = f.read()

        if 'config/VAHTS-SmallRNA-V2.fa' in content:
            print("✅ .gitignore 正确配置了VAHTS接头序列文件")
        else:
            print("❌ .gitignore 缺少VAHTS接头序列文件配置")
            return False

        if 'openspec/changes/archive/' in content:
            print("✅ .gitignore 正确配置了OpenSpec归档目录")
        else:
            print("❌ .gitignore 缺少OpenSpec归档目录配置")
            return False

    except Exception as e:
        print(f"❌ 读取 .gitignore 出错: {e}")
        return False

    return True


def check_environment_config():
    """检查Conda环境配置"""
    print("\n=== 检查Conda环境配置 ===")
    env_path = Path("envs/small_rna_analysis.yaml")
    if not env_path.exists():
        print("❌ small_rna_analysis.yaml 不存在")
        return False

    print("✅ Conda环境配置文件存在")
    return True


def run_all_checks():
    """运行所有检查"""
    print("=" * 50)
    print("项目最终状态检查")
    print("=" * 50)

    checks = [
        check_config,
        check_vahts_file,
        check_trim_script,
        check_project_structure,
        check_git_ignore,
        check_environment_config
    ]

    all_passed = True
    for check_func in checks:
        if not check_func():
            all_passed = False

    print("\n" + "=" * 50)
    if all_passed:
        print("✅ 项目状态检查全部通过！")
    else:
        print("❌ 项目状态检查存在问题！")

    return all_passed

if __name__ == "__main__":
    success = run_all_checks()
    sys.exit(0 if success else 1)
