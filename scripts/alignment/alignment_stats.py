#!/usr/bin/env python3
"""
比对统计脚本

功能：
1. 从Bowtie2日志文件或BAM文件计算比对统计
2. 计算比对率、唯一比对率、多重比对率
3. 生成比对统计汇总报告
4. 支持批量处理多个样本

使用方法：
    python alignment_stats.py --input <日志文件目录或BAM文件目录> --output <输出目录>
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Any
import logging
import json
import yaml

# 导入项目的日志配置工具
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.logging_utils import get_script_logger

# 配置日志
logger = get_script_logger('alignment_stats')


class AlignmentStatsCalculator:
    """比对统计计算器"""

    def __init__(self, samtools_path: str = "samtools"):
        """
        初始化统计计算器

        参数:
            samtools_path: samtools可执行文件路径
        """
        self.samtools_path = samtools_path
        self.results = {}

    def calculate_from_log(self, log_file: str, sample_name: Optional[str] = None) -> Dict[str, Any]:
        """
        从Bowtie2日志文件计算统计

        参数:
            log_file: Bowtie2日志文件路径
            sample_name: 样本名称

        返回:
            Dict: 比对统计
        """
        if sample_name is None:
            sample_name = Path(log_file).stem.replace('_bowtie2', '')

        logger.info(f"从日志文件计算统计: {sample_name}")

        stats = {
            'sample': sample_name,
            'log_file': log_file,
            'source': 'bowtie2_log'
        }

        try:
            with open(log_file, 'r') as f:
                content = f.read()

            lines = content.split('\n')
            for line in lines:
                line = line.strip()
                # 解析Bowtie2统计信息
                if 'reads; of these:' in line:
                    stats['total_reads'] = int(line.split()[0])
                elif 'aligned 0 times' in line:
                    stats['failed_to_align'] = int(line.split()[0])
                elif 'aligned exactly 1 time' in line:
                    stats['unique_alignments'] = int(line.split()[0])
                elif 'aligned >1 times' in line:
                    stats['multiple_alignments'] = int(line.split()[0])
                elif 'overall alignment rate' in line:
                    rate_str = line.split()[0].replace('%', '')
                    stats['alignment_rate'] = float(rate_str) / 100.0

            # 计算总比对reads
            stats['aligned_reads'] = stats.get('unique_alignments', 0) + stats.get('multiple_alignments', 0)

            # 计算唯一比对率
            if stats.get('aligned_reads', 0) > 0:
                stats['unique_alignment_rate'] = stats.get('unique_alignments', 0) / stats['aligned_reads']
                stats['multiple_alignment_rate'] = stats.get('multiple_alignments', 0) / stats['aligned_reads']
            else:
                stats['unique_alignment_rate'] = 0.0
                stats['multiple_alignment_rate'] = 0.0

            stats['success'] = True

        except Exception as e:
            logger.error(f"解析日志文件时出错: {e}")
            stats.update({
                'success': False,
                'error': str(e)
            })

        self.results[sample_name] = stats
        return stats

    def calculate_from_bam(self, bam_file: str, sample_name: Optional[str] = None) -> Dict[str, Any]:
        """
        从BAM文件计算统计

        参数:
            bam_file: BAM文件路径
            sample_name: 样本名称

        返回:
            Dict: 比对统计
        """
        if sample_name is None:
            sample_name = Path(bam_file).stem
            if sample_name.endswith('_sorted'):
                sample_name = sample_name[:-7]

        logger.info(f"从BAM文件计算统计: {sample_name}")

        stats = {
            'sample': sample_name,
            'bam_file': bam_file,
            'source': 'bam_file'
        }

        try:
            # 使用samtools flagstat获取基本统计
            cmd = [self.samtools_path, "flagstat", bam_file]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            lines = result.stdout.strip().split('\n')
            for line in lines:
                if 'in total' in line:
                    stats['total_reads'] = int(line.split()[0])
                elif 'mapped (' in line:
                    stats['mapped_reads'] = int(line.split()[0])
                elif 'paired in sequencing' in line:
                    stats['paired_reads'] = int(line.split()[0])
                elif 'properly paired' in line:
                    stats['properly_paired'] = int(line.split()[0])
                elif 'duplicates' in line:
                    stats['duplicate_reads'] = int(line.split()[0])

            # 计算比对率
            if stats.get('total_reads', 0) > 0:
                stats['alignment_rate'] = stats.get('mapped_reads', 0) / stats['total_reads']
                if stats.get('paired_reads', 0) > 0:
                    stats['proper_pair_rate'] = stats.get('properly_paired', 0) / stats['paired_reads']
                stats['duplicate_rate'] = stats.get('duplicate_reads', 0) / stats['total_reads']

            # 获取唯一/多重比对统计（需要额外处理）
            # 这里使用简单估计：通过MAPQ分布
            mapq_stats = self._estimate_unique_multiple_from_bam(bam_file)
            stats.update(mapq_stats)

            stats['success'] = True

        except Exception as e:
            logger.error(f"分析BAM文件时出错: {e}")
            stats.update({
                'success': False,
                'error': str(e)
            })

        self.results[sample_name] = stats
        return stats

    def _estimate_unique_multiple_from_bam(self, bam_file: str) -> Dict[str, Any]:
        """从BAM文件估计唯一和多重比对"""
        stats = {}

        try:
            # 提取MAPQ值，高MAPQ通常表示唯一比对
            cmd = f"{self.samtools_path} view {bam_file} | awk '{{print $5}}' | head -10000"
            result = subprocess.run(
                cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode == 0 and result.stdout.strip():
                mapq_values = [int(x) for x in result.stdout.strip().split() if x]
                if mapq_values:
                    # 简单启发式：MAPQ >= 30 认为是唯一比对
                    high_mapq = sum(1 for m in mapq_values if m >= 30)
                    low_mapq = len(mapq_values) - high_mapq

                    stats['estimated_unique_alignments'] = high_mapq
                    stats['estimated_multiple_alignments'] = low_mapq

                    if len(mapq_values) > 0:
                        stats['estimated_unique_rate'] = high_mapq / len(mapq_values)
                        stats['estimated_multiple_rate'] = low_mapq / len(mapq_values)

        except Exception as e:
            logger.warning(f"估计唯一/多重比对时出错: {e}")

        return stats

    def save_stats(self, stats: Dict[str, Any], output_dir: Path):
        """保存统计信息"""
        sample_name = stats.get('sample', 'unknown')
        output_dir.mkdir(parents=True, exist_ok=True)

        # 文本报告
        txt_file = output_dir / f"{sample_name}_alignment_stats.txt"
        with open(txt_file, 'w') as f:
            f.write(f"=== 比对统计报告 ===\n\n")
            f.write(f"样本: {sample_name}\n")
            f.write(f"数据源: {stats.get('source', 'unknown')}\n\n")

            f.write("基本统计:\n")
            f.write(f"  总reads数: {stats.get('total_reads', 0):,}\n")
            f.write(f"  比对reads数: {stats.get('aligned_reads', stats.get('mapped_reads', 0)):,}\n")
            f.write(f"  比对率: {stats.get('alignment_rate', 0)*100:.2f}%\n\n")

            if 'unique_alignments' in stats:
                f.write("比对类型统计:\n")
                f.write(f"  唯一比对: {stats.get('unique_alignments', 0):,}\n")
                f.write(f"  多重比对: {stats.get('multiple_alignments', 0):,}\n")
                f.write(f"  未比对: {stats.get('failed_to_align', 0):,}\n")
                f.write(f"  唯一比对率: {stats.get('unique_alignment_rate', 0)*100:.2f}%\n")
                f.write(f"  多重比对率: {stats.get('multiple_alignment_rate', 0)*100:.2f}%\n\n")

            if 'paired_reads' in stats:
                f.write("双端统计:\n")
                f.write(f"  双端reads数: {stats.get('paired_reads', 0):,}\n")
                f.write(f"  正确配对: {stats.get('properly_paired', 0):,}\n")
                f.write(f"  正确配对率: {stats.get('proper_pair_rate', 0)*100:.2f}%\n\n")

            f.write(f"重复率: {stats.get('duplicate_rate', 0)*100:.2f}%\n\n")

            if 'estimated_unique_alignments' in stats:
                f.write("估计比对类型（基于MAPQ）:\n")
                f.write(f"  估计唯一比对: {stats.get('estimated_unique_alignments', 0):,}\n")
                f.write(f"  估计多重比对: {stats.get('estimated_multiple_alignments', 0):,}\n")
                f.write(f"  估计唯一比对率: {stats.get('estimated_unique_rate', 0)*100:.2f}%\n\n")

            f.write(f"生成时间: {pd.Timestamp.now()}\n")

        # JSON格式
        json_file = output_dir / f"{sample_name}_alignment_stats.json"
        with open(json_file, 'w') as f:
            json.dump(stats, f, indent=2)

        logger.info(f"比对统计已保存: {txt_file}")
        return str(txt_file)

    def generate_summary_report(self, output_dir: str) -> Optional[str]:
        """
        生成比对统计汇总报告

        参数:
            output_dir: 输出目录

        返回:
            str: 报告文件路径
        """
        if not self.results:
            logger.warning("无比对统计结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建DataFrame
        summary_data = []
        for sample, stats in self.results.items():
            if stats.get('success'):
                summary_data.append({
                    'sample': sample,
                    'total_reads': stats.get('total_reads', 0),
                    'aligned_reads': stats.get('aligned_reads', stats.get('mapped_reads', 0)),
                    'alignment_rate': stats.get('alignment_rate', 0),
                    'unique_alignments': stats.get('unique_alignments', stats.get('estimated_unique_alignments', 0)),
                    'multiple_alignments': stats.get('multiple_alignments', stats.get('estimated_multiple_alignments', 0)),
                    'unique_rate': stats.get('unique_alignment_rate', stats.get('estimated_unique_rate', 0)),
                    'duplicate_rate': stats.get('duplicate_rate', 0),
                    'source': stats.get('source', 'unknown')
                })

        if not summary_data:
            return None

        df = pd.DataFrame(summary_data)

        # 保存CSV
        csv_file = output_dir / "alignment_stats_summary.csv"
        df.to_csv(csv_file, index=False)

        # 生成文本报告
        report_file = output_dir / "alignment_stats_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== 比对统计汇总报告 ===\n\n")
            f.write(f"总样本数: {len(summary_data)}\n\n")

            f.write("样本比对统计:\n")
            for _, row in df.iterrows():
                f.write(f"  - {row['sample']}:\n")
                f.write(f"    总reads: {row['total_reads']:,}\n")
                f.write(f"    比对reads: {row['aligned_reads']:,}\n")
                f.write(f"    比对率: {row['alignment_rate']*100:.1f}%\n")
                if not pd.isna(row['unique_rate']):
                    f.write(f"    唯一比对率: {row['unique_rate']*100:.1f}%\n")
                if not pd.isna(row['duplicate_rate']):
                    f.write(f"    重复率: {row['duplicate_rate']*100:.1f}%\n")
                f.write("\n")

            f.write("总体统计:\n")
            f.write(f"  平均比对率: {df['alignment_rate'].mean()*100:.1f}%\n")
            if 'unique_rate' in df.columns and not df['unique_rate'].isna().all():
                f.write(f"  平均唯一比对率: {df['unique_rate'].mean()*100:.1f}%\n")
            if 'duplicate_rate' in df.columns and not df['duplicate_rate'].isna().all():
                f.write(f"  平均重复率: {df['duplicate_rate'].mean()*100:.1f}%\n")

            f.write(f"\n报告生成时间: {pd.Timestamp.now()}\n")

        logger.info(f"比对统计汇总报告已生成: {report_file}")
        return str(report_file)


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """加载配置文件"""
    default_config = {
        'samtools_path': 'samtools',
        'threads': 4
    }

    if config_file and Path(config_file).exists():
        try:
            with open(config_file, 'r') as f:
                user_config = yaml.safe_load(f)
                default_config.update(user_config)
                logger.info(f"已加载配置文件: {config_file}")
        except Exception as e:
            logger.warning(f"加载配置文件失败: {e}，使用默认配置")

    return default_config


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="比对统计脚本"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入文件或目录（Bowtie2日志文件、BAM文件或包含这些文件的目录）"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/alignment/stats",
        help="输出目录或CSV文件路径 (默认: results/alignment/stats)"
    )

    parser.add_argument(
        "--type", "-t",
        choices=['log', 'bam', 'auto'],
        default='auto',
        help="输入文件类型: log(Bowtie2日志), bam(BAM文件), auto(自动检测)"
    )

    parser.add_argument(
        "--config", "-c",
        help="配置文件 (YAML格式)"
    )

    parser.add_argument(
        "--sample-info",
        help="样本信息CSV文件，包含sample和file_path列"
    )

    parser.add_argument(
        "--summary",
        action="store_true",
        default=True,
        help="生成汇总报告 (默认: True)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 加载配置
    config = load_config(args.config)

    # 初始化统计计算器
    calculator = AlignmentStatsCalculator(samtools_path=config.get('samtools_path', 'samtools'))

    # 处理输入
    input_path = Path(args.input)
    output_path = Path(args.output)

    # 判断输出是文件还是目录
    is_file_output = False
    if output_path.suffix == '.csv':
        is_file_output = True
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = output_path
        output_dir.mkdir(parents=True, exist_ok=True)

    # 根据输入类型处理
    if args.sample_info:
        # 使用样本信息CSV文件
        logger.info(f"使用样本信息文件: {args.sample_info}")
        df = pd.read_csv(args.sample_info)

        for _, row in df.iterrows():
            sample = row.get('sample', f"sample_{_}")
            file_path = row.get('file_path')
            file_type = row.get('type', args.type)

            if pd.notna(file_path):
                if file_type == 'log' or (file_type == 'auto' and Path(file_path).suffix in ['.log', '.txt']):
                    calculator.calculate_from_log(file_path, sample)
                elif file_type == 'bam' or (file_type == 'auto' and Path(file_path).suffix == '.bam'):
                    calculator.calculate_from_bam(file_path, sample)

    elif input_path.is_file():
        # 单个文件
        if args.type == 'log' or (args.type == 'auto' and input_path.suffix in ['.log', '.txt']):
            calculator.calculate_from_log(str(input_path))
        elif args.type == 'bam' or (args.type == 'auto' and input_path.suffix == '.bam'):
            calculator.calculate_from_bam(str(input_path))

    elif input_path.is_dir():
        # 目录下的所有文件
        logger.info(f"扫描目录: {input_path}")

        if args.type in ['log', 'auto']:
            # 查找日志文件
            log_files = list(input_path.glob("*.log")) + list(input_path.glob("*_bowtie2.log")) + \
                       list(input_path.glob("*.txt"))
            for log_file in log_files:
                calculator.calculate_from_log(str(log_file))

        if args.type in ['bam', 'auto']:
            # 查找BAM文件
            bam_files = list(input_path.glob("*.bam")) + list(input_path.glob("*.BAM"))
            for bam_file in bam_files:
                calculator.calculate_from_bam(str(bam_file))

    else:
        logger.error("不支持的输入类型")
        sys.exit(1)

    # 处理输出
    if is_file_output:
        # 直接输出到CSV文件
        logger.info(f"直接保存到CSV文件: {output_path}")
        if len(calculator.results) > 0:
            # 获取第一个样本的统计结果
            first_sample = next(iter(calculator.results.values()))
            if first_sample.get('success'):
                # 转换为DataFrame并保存到CSV
                df = pd.DataFrame([first_sample])
                df.to_csv(output_path, index=False)
                logger.info(f"统计结果已保存到: {output_path}")
            else:
                logger.error("统计计算失败，无法保存到CSV文件")
                sys.exit(1)
    else:
        # 输出到目录
        for sample, stats in calculator.results.items():
            if stats.get('success'):
                calculator.save_stats(stats, output_dir)

        # 生成汇总报告
        if args.summary:
            report_file = calculator.generate_summary_report(output_dir)
            if report_file:
                logger.info(f"比对统计完成，报告文件: {report_file}")
            else:
                logger.warning("未能生成汇总报告")

    logger.info("比对统计流程完成")


if __name__ == "__main__":
    main()