#!/usr/bin/env python3
"""
质量控制汇总脚本 - 整合FastQC和Trimmomatic结果

功能：
1. 收集所有样本的FastQC和Trimmomatic结果
2. 生成综合质量控制报告
3. 可视化质量指标
4. 识别问题样本

使用方法：
    python qc_summary.py --fastqc-dir <fastqc结果目录> --trim-dir <trimmomatic结果目录> --output <输出目录>
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Optional, Any
import logging
import matplotlib.pyplot as plt
import seaborn as sns

# 导入项目的日志配置工具
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.logging_utils import get_script_logger

# 配置日志
logger = get_script_logger('qc_summary')


class QCSummary:
    """质量控制汇总类"""

    def __init__(self):
        self.fastqc_data = {}
        self.trim_data = {}
        self.combined_data = {}

    def load_fastqc_results(self, fastqc_dir: str) -> bool:
        """
        加载FastQC结果

        参数:
            fastqc_dir: FastQC结果目录

        返回:
            bool: 是否成功
        """
        fastqc_dir = Path(fastqc_dir)
        if not fastqc_dir.exists():
            logger.error(f"FastQC目录不存在: {fastqc_dir}")
            return False

        # 查找FastQC汇总文件
        summary_files = list(fastqc_dir.glob("*fastqc_summary.csv")) + \
                       list(fastqc_dir.glob("*fastqc_summary.json"))

        if not summary_files:
            logger.warning(f"未找到FastQC汇总文件，尝试解析原始结果")
            return self._parse_raw_fastqc_results(fastqc_dir)

        # 加载CSV汇总文件
        csv_files = [f for f in summary_files if f.suffix == '.csv']
        if csv_files:
            try:
                df = pd.read_csv(csv_files[0])
                for _, row in df.iterrows():
                    sample = row.get('sample')
                    if sample:
                        self.fastqc_data[sample] = row.to_dict()
                logger.info(f"加载FastQC数据: {len(self.fastqc_data)} 个样本")
                return True
            except Exception as e:
                logger.error(f"加载FastQC CSV文件失败: {e}")

        # 加载JSON汇总文件
        json_files = [f for f in summary_files if f.suffix == '.json']
        if json_files:
            try:
                with open(json_files[0], 'r') as f:
                    self.fastqc_data = json.load(f)
                logger.info(f"加载FastQC数据: {len(self.fastqc_data)} 个样本")
                return True
            except Exception as e:
                logger.error(f"加载FastQC JSON文件失败: {e}")

        return False

    def _parse_raw_fastqc_results(self, fastqc_dir: Path) -> bool:
        """解析原始FastQC结果文件"""
        # 查找FastQC HTML报告
        html_files = list(fastqc_dir.glob("*_fastqc.html"))
        if not html_files:
            logger.error("未找到FastQC结果文件")
            return False

        for html_file in html_files:
            sample_name = html_file.name.replace('_fastqc.html', '')
            self.fastqc_data[sample_name] = {
                'sample': sample_name,
                'html_report': str(html_file),
                'status': 'parsed'
            }

        logger.info(f"解析FastQC数据: {len(self.fastqc_data)} 个样本")
        return True

    def load_trimmomatic_results(self, trim_dir: str) -> bool:
        """
        加载Trimmomatic结果

        参数:
            trim_dir: Trimmomatic结果目录

        返回:
            bool: 是否成功
        """
        trim_dir = Path(trim_dir)
        if not trim_dir.exists():
            logger.error(f"Trimmomatic目录不存在: {trim_dir}")
            return False

        # 查找Trimmomatic汇总文件
        summary_files = list(trim_dir.glob("*trimmomatic_summary.csv")) + \
                       list(trim_dir.glob("*trimmomatic_summary.json"))

        if not summary_files:
            logger.warning(f"未找到Trimmomatic汇总文件")
            return False

        # 加载CSV汇总文件
        csv_files = [f for f in summary_files if f.suffix == '.csv']
        if csv_files:
            try:
                df = pd.read_csv(csv_files[0])
                for _, row in df.iterrows():
                    sample = row.get('sample')
                    if sample:
                        self.trim_data[sample] = row.to_dict()
                logger.info(f"加载Trimmomatic数据: {len(self.trim_data)} 个样本")
                return True
            except Exception as e:
                logger.error(f"加载Trimmomatic CSV文件失败: {e}")

        # 加载JSON汇总文件
        json_files = [f for f in summary_files if f.suffix == '.json']
        if json_files:
            try:
                with open(json_files[0], 'r') as f:
                    self.trim_data = json.load(f)
                logger.info(f"加载Trimmomatic数据: {len(self.trim_data)} 个样本")
                return True
            except Exception as e:
                logger.error(f"加载Trimmomatic JSON文件失败: {e}")

        return False

    def combine_results(self):
        """合并FastQC和Trimmomatic结果"""
        all_samples = set(list(self.fastqc_data.keys()) + list(self.trim_data.keys()))

        for sample in all_samples:
            combined = {
                'sample': sample,
                'has_fastqc': sample in self.fastqc_data,
                'has_trimmomatic': sample in self.trim_data
            }

            # 添加FastQC数据
            if sample in self.fastqc_data:
                fastqc_info = self.fastqc_data[sample]
                if isinstance(fastqc_info, dict):
                    combined.update({f'fastqc_{k}': v for k, v in fastqc_info.items()})

            # 添加Trimmomatic数据
            if sample in self.trim_data:
                trim_info = self.trim_data[sample]
                if isinstance(trim_info, dict):
                    combined.update({f'trim_{k}': v for k, v in trim_info.items()})

            self.combined_data[sample] = combined

        logger.info(f"合并数据: {len(self.combined_data)} 个样本")

    def generate_report(self, output_dir: str) -> Optional[str]:
        """
        生成综合报告

        参数:
            output_dir: 输出目录

        返回:
            str: 报告文件路径
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        if not self.combined_data:
            logger.warning("无数据可生成报告")
            return None

        # 创建DataFrame
        df = pd.DataFrame(list(self.combined_data.values()))

        # 1. 保存为CSV
        csv_file = output_dir / "qc_combined_summary.csv"
        df.to_csv(csv_file, index=False)
        logger.info(f"保存CSV汇总: {csv_file}")

        # 2. 保存为JSON
        json_file = output_dir / "qc_combined_summary.json"
        with open(json_file, 'w') as f:
            json.dump(self.combined_data, f, indent=2)
        logger.info(f"保存JSON汇总: {json_file}")

        # 3. 生成文本报告
        report_file = output_dir / "qc_summary_report.txt"
        self._generate_text_report(report_file, df)
        logger.info(f"生成文本报告: {report_file}")

        # 4. 生成可视化
        try:
            self._generate_visualizations(df, output_dir)
        except Exception as e:
            logger.warning(f"生成可视化时出错: {e}")

        return str(report_file)

    def _generate_text_report(self, report_file: Path, df: pd.DataFrame):
        """生成文本报告"""
        with open(report_file, 'w') as f:
            f.write("=== 质量控制综合报告 ===\n\n")

            # 基本统计
            total_samples = len(df)
            fastqc_samples = df['has_fastqc'].sum()
            trim_samples = df['has_trimmomatic'].sum()

            f.write(f"总样本数: {total_samples}\n")
            f.write(f"有FastQC结果的样本: {fastqc_samples}\n")
            f.write(f"有Trimmomatic结果的样本: {trim_samples}\n\n")

            # Trimmomatic统计
            if trim_samples > 0:
                f.write("Trimmomatic修剪统计:\n")
                if 'trim_input_reads' in df.columns:
                    total_input = df['trim_input_reads'].sum()
                    total_surviving = df['trim_surviving_reads'].sum()
                    if total_input > 0:
                        survival_rate = total_surviving / total_input * 100
                        f.write(f"  总输入reads: {total_input:,}\n")
                        f.write(f"  总保留reads: {total_surviving:,}\n")
                        f.write(f"  总保留率: {survival_rate:.1f}%\n\n")

                # 各样本保留率
                f.write("各样本保留率:\n")
                for _, row in df.iterrows():
                    if row['has_trimmomatic'] and 'trim_survival_rate' in row:
                        rate = row['trim_survival_rate']
                        if isinstance(rate, (int, float)):
                            f.write(f"  - {row['sample']}: {rate*100:.1f}%\n")

            # 问题样本检测
            f.write("\n问题样本检测:\n")
            problem_samples = []

            # 低保留率样本
            if 'trim_survival_rate' in df.columns:
                low_survival = df[df['trim_survival_rate'] < 0.5]
                for _, row in low_survival.iterrows():
                    problem_samples.append(f"{row['sample']} (保留率: {row['trim_survival_rate']*100:.1f}%)")

            # 无FastQC结果样本
            no_fastqc = df[~df['has_fastqc']]
            for _, row in no_fastqc.iterrows():
                problem_samples.append(f"{row['sample']} (缺少FastQC结果)")

            # 无Trimmomatic结果样本
            no_trim = df[~df['has_trimmomatic']]
            for _, row in no_trim.iterrows():
                problem_samples.append(f"{row['sample']} (缺少Trimmomatic结果)")

            if problem_samples:
                for problem in problem_samples:
                    f.write(f"  - {problem}\n")
            else:
                f.write("  未检测到明显问题\n")

            f.write(f"\n报告生成时间: {pd.Timestamp.now()}\n")

    def _generate_visualizations(self, df: pd.DataFrame, output_dir: Path):
        """生成可视化图表"""
        plt.style.use('seaborn-v0_8')

        # 1. Trimmomatic保留率柱状图
        if 'trim_survival_rate' in df.columns and df['has_trimmomatic'].any():
            plt.figure(figsize=(10, 6))
            trim_df = df[df['has_trimmomatic']].copy()
            trim_df = trim_df.sort_values('trim_survival_rate', ascending=False)

            bars = plt.bar(range(len(trim_df)), trim_df['trim_survival_rate'] * 100)
            plt.xticks(range(len(trim_df)), trim_df['sample'], rotation=45, ha='right')
            plt.ylabel('保留率 (%)')
            plt.title('Trimmomatic reads保留率')
            plt.tight_layout()

            # 添加数值标签
            for bar in bars:
                height = bar.get_height()
                plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                        f'{height:.1f}%', ha='center', va='bottom', fontsize=9)

            plt.savefig(output_dir / 'trim_survival_rate.png', dpi=300)
            plt.close()

        # 2. 输入reads vs 保留reads散点图
        if 'trim_input_reads' in df.columns and 'trim_surviving_reads' in df.columns:
            plt.figure(figsize=(8, 6))
            plt.scatter(df['trim_input_reads'], df['trim_surviving_reads'], alpha=0.7)
            plt.xlabel('输入reads数')
            plt.ylabel('保留reads数')
            plt.title('输入 vs 保留 reads')

            # 添加对角线
            max_val = max(df['trim_input_reads'].max(), df['trim_surviving_reads'].max())
            plt.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)

            plt.tight_layout()
            plt.savefig(output_dir / 'input_vs_surviving.png', dpi=300)
            plt.close()

        # 3. 样本状态热图
        status_data = []
        for sample, data in self.combined_data.items():
            status_data.append({
                'sample': sample,
                'FastQC': 1 if data.get('has_fastqc') else 0,
                'Trimmomatic': 1 if data.get('has_trimmomatic') else 0
            })

        status_df = pd.DataFrame(status_data).set_index('sample')
        if len(status_df) > 0:
            plt.figure(figsize=(8, max(4, len(status_df) * 0.3)))
            sns.heatmap(status_df, annot=True, cmap='RdYlGn', cbar=False,
                       linewidths=1, linecolor='gray')
            plt.title('样本分析状态')
            plt.tight_layout()
            plt.savefig(output_dir / 'sample_status_heatmap.png', dpi=300)
            plt.close()

        logger.info(f"生成可视化图表到: {output_dir}")


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="质量控制汇总脚本 - 整合FastQC和Trimmomatic结果"
    )

    parser.add_argument(
        "--fastqc-dir",
        required=True,
        help="FastQC结果目录"
    )

    parser.add_argument(
        "--trim-dir",
        required=True,
        help="Trimmomatic结果目录"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/qc/summary",
        help="输出目录 (默认: results/qc/summary)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 初始化汇总器
    qc_summary = QCSummary()

    # 加载数据
    logger.info("加载FastQC结果...")
    if not qc_summary.load_fastqc_results(args.fastqc_dir):
        logger.warning("FastQC数据加载失败，继续处理")

    logger.info("加载Trimmomatic结果...")
    if not qc_summary.load_trimmomatic_results(args.trim_dir):
        logger.warning("Trimmomatic数据加载失败，继续处理")

    # 合并结果
    qc_summary.combine_results()

    # 生成报告
    report_file = qc_summary.generate_report(args.output)
    if report_file:
        logger.info(f"质量控制报告已生成: {report_file}")
    else:
        logger.error("未能生成质量控制报告")
        sys.exit(1)

    logger.info("质量控制汇总流程完成")


if __name__ == "__main__":
    main()