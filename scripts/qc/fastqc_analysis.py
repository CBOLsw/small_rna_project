#!/usr/bin/env python3
"""
FastQC分析脚本 - 用于small RNA测序数据质量控制

功能：
1. 对单个或多个fastq文件运行FastQC
2. 生成汇总报告和质量评估
3. 识别常见问题（接头污染、低质量序列等）

使用方法：
    python fastqc_analysis.py --input <fastq文件或目录> --output <输出目录>
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

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


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

    def run_fastqc(self, input_path: str, output_dir: str, threads: int = 4, sample_name: str = None) -> bool:
        """
        运行FastQC分析

        参数:
            input_path: 输入fastq文件或目录
            output_dir: 输出目录
            threads: 线程数
            sample_name: 样本名称，用于指定输出文件名前缀（可选）

        返回:
            bool: 是否成功
        """
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
                    # 解析结果
                    self._parse_fastqc_result(fastq_file, output_dir, sample_name)
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
        # FastQC输出文件名模式
        # 首先尝试使用原始文件名
        base_name = fastq_file.name
        if base_name.endswith('.gz'):
            base_name = base_name[:-3]
        if base_name.endswith('.fastq') or base_name.endswith('.fq'):
            base_name = base_name[:-6]

        result_file = output_dir / f"{base_name}_fastqc.html"
        data_file = output_dir / f"{base_name}_fastqc.zip"

        # 如果提供了sample_name，并且原始文件名的结果不存在，则尝试重命名
        if sample_name and not result_file.exists():
            # 使用sample_name来查找结果文件
            sample_result_file = output_dir / f"{sample_name}_fastqc.html"
            sample_data_file = output_dir / f"{sample_name}_fastqc.zip"

            # 如果sample_name的结果文件存在，则使用它们
            if sample_result_file.exists():
                result_file = sample_result_file
                data_file = sample_data_file
            else:
                # 否则，尝试将原始文件名的结果重命名为sample_name
                original_result_file = output_dir / f"{base_name}_fastqc.html"
                original_data_file = output_dir / f"{base_name}_fastqc.zip"

                if original_result_file.exists():
                    import shutil
                    shutil.move(str(original_result_file), str(sample_result_file))
                    shutil.move(str(original_data_file), str(sample_data_file))
                    result_file = sample_result_file
                    data_file = sample_data_file

        if result_file.exists():
            # 提取基本信息
            sample_name = base_name
            self.results[sample_name] = {
                'file': str(fastq_file),
                'html_report': str(result_file),
                'data_file': str(data_file),
                'status': 'completed'
            }
            logger.info(f"解析结果: {sample_name}")
        else:
            logger.warning(f"未找到结果文件: {result_file}")

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

        # 创建汇总DataFrame
        summary_data = []
        for sample, result in self.results.items():
            summary_data.append({
                'sample': sample,
                'file': result['file'],
                'status': result['status'],
                'html_report': result['html_report']
            })

        df = pd.DataFrame(summary_data)

        # 保存为CSV
        csv_file = output_dir / "fastqc_summary.csv"
        df.to_csv(csv_file, index=False)

        # 保存为JSON
        json_file = output_dir / "fastqc_summary.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        # 生成简单文本报告
        report_file = output_dir / "fastqc_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== FastQC分析汇总报告 ===\n\n")
            f.write(f"总样本数: {len(self.results)}\n")
            f.write(f"成功分析: {len([r for r in self.results.values() if r['status'] == 'completed'])}\n\n")

            f.write("样本列表:\n")
            for sample in sorted(self.results.keys()):
                f.write(f"  - {sample}: {self.results[sample]['status']}\n")

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

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

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
        sample_name=args.sample
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