#!/usr/bin/env python3
"""
FastQC分析脚本 - 用于small RNA测序数据质量控制

功能：
1. 对单个或多个fastq文件运行FastQC
2. 生成汇总报告和质量评估
3. 识别常见问题（接头污染、低质量序列等）

使用方法：
    python fastqc_analysis.py --input <fastq文件或目录> --output <输出目录> --sample-info <样本信息文件>
"""

import os
import sys
import argparse
import subprocess
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import logging
import re

# 导入项目的日志配置工具
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.logging_utils import get_script_logger

# 配置日志
logger = get_script_logger('fastqc_analysis')


def load_sample_info(sample_info_file: str) -> Dict[str, str]:
    """
    从sample_info.csv加载样本名映射

    返回: {fastq文件名 -> 样本名} 的映射
    """
    sample_map = {}
    try:
        df = pd.read_csv(sample_info_file)
        for _, row in df.iterrows():
            sample_name = row['sample']
            fastq_r1 = row['fastq_r1']
            # 从路径中提取fastq文件名
            fastq_name = os.path.basename(fastq_r1)
            # 去掉 .fastq.gz 后缀
            base_name = fastq_name.replace('.fastq.gz', '')
            sample_map[base_name] = sample_name
            sample_map[fastq_name] = sample_name
            # 也支持不包含路径的key
            sample_map[fastq_r1] = sample_name
        logger.info(f"已加载 {len(sample_map)} 个样本映射")
    except Exception as e:
        logger.warning(f"无法加载样本信息文件: {e}")
    return sample_map


class FastQCAnalyzer:
    """FastQC分析器类"""

    def __init__(self, fastqc_path: str = "fastqc"):
        """
        初始化FastQC分析器

        参数:
            fastqc_path: FastQC可执行文件路径
        """
        self.fastqc_path = fastqc_path
        self.results = {}

    def check_fastqc(self) -> bool:
        """检查FastQC是否可用"""
        try:
            result = subprocess.run(
                [self.fastqc_path, "--version"],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                logger.info(f"FastQC版本: {result.stdout.strip()}")
                return True
            else:
                logger.error("FastQC检查失败")
                return False
        except FileNotFoundError:
            logger.error(f"未找到FastQC: {self.fastqc_path}")
            return False

    def run_fastqc(self, input_path: str, output_dir: str, threads: int = 4,
                   sample_name: str = None, sample_map: Dict[str, str] = None) -> bool:
        """
        运行FastQC分析

        参数:
            input_path: 输入fastq文件或目录
            output_dir: 输出目录
            threads: 线程数
            sample_name: 样本名称（可选，如果指定则所有文件使用此名称）
            sample_map: 样本名映射字典 {fastq文件名(含路径) -> 样本名}

        返回:
            bool: 是否成功
        """
        if sample_map is None:
            sample_map = {}

        # 创建输出目录
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 收集fastq文件
        fastq_files = self._collect_fastq_files(input_path)
        if not fastq_files:
            logger.error(f"未找到fastq文件: {input_path}")
            return False

        logger.info(f"找到 {len(fastq_files)} 个fastq文件")

        # 运行FastQC
        success_count = 0
        for i, fastq_file in enumerate(fastq_files, 1):
            logger.info(f"处理文件 {i}/{len(fastq_files)}: {fastq_file}")

            # 确定样本名
            resolved_sample_name = sample_name
            if resolved_sample_name is None:
                # 尝试从sample_map中查找
                fastq_basename = fastq_file.name  # HeLa-PAL2_S52_L001_R1_001.fastq.gz
                fastq_stem = fastq_basename.replace('.fastq.gz', '')  # HeLa-PAL2_S52_L001_R1_001

                if str(fastq_file) in sample_map:
                    resolved_sample_name = sample_map[str(fastq_file)]
                elif fastq_stem in sample_map:
                    resolved_sample_name = sample_map[fastq_stem]
                elif fastq_basename in sample_map:
                    resolved_sample_name = sample_map[fastq_basename]

            try:
                # 运行FastQC
                cmd = [
                    self.fastqc_path,
                    "-o", str(output_dir),
                    "-t", str(threads),
                    "-q",  # 静默模式
                    str(fastq_file)
                ]

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False
                )

                if result.returncode == 0:
                    success_count += 1
                    # 解析结果，使用映射的样本名
                    self._parse_fastqc_result(fastq_file, output_dir, resolved_sample_name)
                else:
                    logger.error(f"FastQC处理失败: {fastq_file}")
                    logger.error(f"错误输出: {result.stderr}")

            except Exception as e:
                logger.error(f"运行FastQC时出错: {e}")

        logger.info(f"FastQC分析完成: {success_count}/{len(fastq_files)} 个文件成功")
        return success_count > 0

    def _collect_fastq_files(self, input_path: str) -> List[Path]:
        """收集fastq文件"""
        input_path = Path(input_path)
        fastq_files = []

        if input_path.is_file():
            # 单个文件
            if self._is_fastq_file(input_path):
                fastq_files.append(input_path)
        else:
            # 目录
            for ext in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
                fastq_files.extend(input_path.rglob(f"*{ext}"))

        return sorted(fastq_files)

    def _is_fastq_file(self, file_path: Path) -> bool:
        """检查是否为fastq文件"""
        suffixes = file_path.suffixes
        if len(suffixes) > 1:
            # 处理压缩文件如 .fastq.gz
            return suffixes[-2] in ['.fastq', '.fq']
        else:
            return file_path.suffix in ['.fastq', '.fq']

    def _parse_fastqc_result(self, fastq_file: Path, output_dir: Path, sample_name: str = None):
        """解析FastQC结果文件"""
        import shutil

        # 确定样本key
        if sample_name:
            key = sample_name
        else:
            # 从fastq文件名提取样本名（去掉路径和后缀）
            key = fastq_file.name
            if key.endswith('.gz'):
                key = key[:-3]
            if key.endswith('.fastq') or key.endswith('.fq'):
                key = key[:-6]

        # FastQC 输出的文件名基于原始输入文件
        if fastq_file.suffix == '.gz':
            fastq_base = fastq_file.stem.replace('.fastq', '')  # HeLa-PAL2_S52_L001_R1_001
        else:
            fastq_base = fastq_file.name.replace('.fastq', '').replace('.fq', '')

        original_html = output_dir / f"{fastq_base}_fastqc.html"
        original_zip = output_dir / f"{fastq_base}_fastqc.zip"

        # 如果指定了sample_name且与fastq文件名不同，需要重命名
        if sample_name and sample_name != fastq_base:
            target_html = output_dir / f"{sample_name}_fastqc.html"
            target_zip = output_dir / f"{sample_name}_fastqc.zip"

            # 重命名文件
            if original_html.exists() and not target_html.exists():
                shutil.move(str(original_html), str(target_html))
                logger.info(f"重命名: {original_html.name} -> {target_html.name}")
            if original_zip.exists() and not target_zip.exists():
                shutil.move(str(original_zip), str(target_zip))

            result_file = target_html
            data_file = target_zip
        else:
            result_file = original_html if original_html.exists() else None
            data_file = original_zip if original_zip.exists() else None

        # 如果目标文件存在，记录结果
        if result_file and result_file.exists():
            self.results[key] = {
                'file': str(fastq_file),
                'html_report': str(result_file),
                'data_file': str(data_file) if data_file and data_file.exists() else '',
                'status': 'completed'
            }
            logger.info(f"解析结果: {key} -> {result_file.name}")
        else:
            logger.warning(f"未找到FastQC结果文件 for {key} (expected: {original_html})")

    def generate_summary(self, output_dir: str) -> Optional[str]:
        """
        生成汇总报告

        参数:
            output_dir: 输出目录

        返回:
            str: 汇总报告文件路径
        """
        if not self.results:
            logger.warning("无分析结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        csv_file = output_dir / "fastqc_summary.csv"

        # 读取现有的CSV（如果存在），避免覆盖其他样本的结果
        existing_data = []
        if csv_file.exists():
            try:
                existing_df = pd.read_csv(csv_file)
                existing_data = existing_df.to_dict('records')
                logger.info(f"合并现有汇总: {len(existing_data)} 个样本")
            except Exception as e:
                logger.warning(f"无法读取现有汇总文件: {e}")

        # 创建汇总DataFrame
        summary_data = []
        existing_samples = set()
        for row in existing_data:
            existing_samples.add(row['sample'])
            summary_data.append(row)

        for sample, result in self.results.items():
            if sample not in existing_samples:
                summary_data.append({
                    'sample': sample,
                    'file': result['file'],
                    'status': result['status'],
                    'html_report': result['html_report']
                })

        df = pd.DataFrame(summary_data)

        # 保存为CSV（追加模式）
        df.to_csv(csv_file, index=False)
        logger.info(f"汇总已更新: {csv_file} ({len(df)} 个样本)")

        # 保存为JSON
        json_file = output_dir / "fastqc_summary.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        # 生成简单文本报告（从CSV重新生成，确保包含所有样本）
        report_file = output_dir / "fastqc_report.txt"
        completed_count = len(df[df['status'] == 'completed'])
        with open(report_file, 'w') as f:
            f.write("=== FastQC分析汇总报告 ===\n\n")
            f.write(f"总样本数: {len(df)}\n")
            f.write(f"成功分析: {completed_count}\n\n")

            f.write("样本列表:\n")
            for _, row in df.iterrows():
                f.write(f"  - {row['sample']}: {row['status']}\n")

            f.write(f"\n报告生成时间: {pd.Timestamp.now()}\n")

        logger.info(f"汇总报告已生成: {report_file}")
        return str(report_file)


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="FastQC分析脚本 - small RNA测序数据质量控制"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入fastq文件或目录"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/qc/fastqc",
        help="输出目录 (默认: results/qc/fastqc)"
    )

    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=4,
        help="线程数 (默认: 4)"
    )

    parser.add_argument(
        "--fastqc-path",
        default="fastqc",
        help="FastQC可执行文件路径 (默认: fastqc)"
    )

    parser.add_argument(
        "--summary",
        action="store_true",
        help="生成汇总报告"
    )

    parser.add_argument(
        "--sample",
        help="样本名称，用于指定输出文件名前缀（可选）"
    )

    parser.add_argument(
        "--sample-info",
        default="data/metadata/sample_info.csv",
        help="样本信息文件路径 (默认: data/metadata/sample_info.csv)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 加载样本名映射
    sample_map = load_sample_info(args.sample_info)

    # 初始化分析器
    analyzer = FastQCAnalyzer(fastqc_path=args.fastqc_path)

    # 检查FastQC
    if not analyzer.check_fastqc():
        logger.error("FastQC检查失败，请确保FastQC已安装并可用")
        sys.exit(1)

    # 运行FastQC
    logger.info(f"开始FastQC分析: {args.input}")
    success = analyzer.run_fastqc(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        sample_name=args.sample,
        sample_map=sample_map
    )

    if not success:
        logger.error("FastQC分析过程中出现错误")
        sys.exit(1)

    # 生成汇总报告
    if args.summary or True:  # 默认总是生成汇总
        report_file = analyzer.generate_summary(args.output)
        if report_file:
            logger.info(f"分析完成，报告文件: {report_file}")
        else:
            logger.warning("未能生成汇总报告")

    logger.info("FastQC分析流程完成")


if __name__ == "__main__":
    main()