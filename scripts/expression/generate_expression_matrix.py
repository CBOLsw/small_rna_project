#!/usr/bin/env python3
"""
基因表达矩阵生成脚本

功能：
1. 整合所有样本的featureCounts输出文件
2. 生成基因表达矩阵（CSV、TSV、Excel格式）
3. 过滤低表达基因
4. 标准化计数（可选）
5. 生成表达矩阵统计报告

使用方法：
    python generate_expression_matrix.py --input <计数文件目录> --output <输出目录>
"""

import os
import sys
import argparse
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
logger = get_script_logger('generate_expression_matrix')


class ExpressionMatrixGenerator:
    """基因表达矩阵生成器"""

    def __init__(self):
        self.results = {}
        self.count_matrix = None

    def load_count_files(self, input_dir: str, pattern: str = "*_counts.txt") -> Dict[str, pd.DataFrame]:
        """
        加载所有计数文件

        参数:
            input_dir: 输入目录
            pattern: 文件匹配模式

        返回:
            Dict: 样本名 -> 计数DataFrame
        """
        input_path = Path(input_dir)
        if not input_path.exists():
            logger.error(f"输入目录不存在: {input_dir}")
            return {}

        # 首先检查是否有 featureCounts 批量输出的完整计数矩阵
        count_matrix_file = input_path / "gene_counts.csv"
        if count_matrix_file.exists():
            logger.info(f"找到 featureCounts 批量输出文件: {count_matrix_file}")
            return self._load_featurecounts_matrix(count_matrix_file)

        count_files = list(input_path.glob(pattern))
        if not count_files:
            logger.warning(f"未找到匹配 {pattern} 的计数文件")
            # 尝试其他常见模式
            count_files = list(input_path.glob("*.txt")) + list(input_path.glob("*.counts"))

        logger.info(f"找到 {len(count_files)} 个计数文件")

        count_data = {}
        for count_file in count_files:
            try:
                # 从文件名提取样本名
                sample_name = count_file.stem
                for suffix in ['_counts', '.counts', '_featurecounts']:
                    if sample_name.endswith(suffix):
                        sample_name = sample_name[:-len(suffix)]

                # 读取计数文件
                df = self._read_count_file(count_file)
                if df is not None:
                    count_data[sample_name] = df
                    logger.info(f"加载成功: {sample_name} ({count_file.name})")
                else:
                    logger.warning(f"无法读取文件: {count_file}")

            except Exception as e:
                logger.warning(f"加载文件失败 {count_file}: {e}")

        return count_data

    def _load_featurecounts_matrix(self, count_file: Path) -> Dict[str, pd.DataFrame]:
        """
        读取 featureCounts 批量输出的完整计数矩阵

        参数:
            count_file: 计数矩阵文件路径

        返回:
            Dict: 样本名 -> 计数DataFrame（单个基因列）
        """
        try:
            # 读取 featureCounts 输出文件（跳过注释行）
            df = pd.read_csv(count_file, sep='\t', comment='#', header=0)

            logger.info(f"featureCounts 矩阵维度: {df.shape}")
            logger.info(f"列名: {list(df.columns)}")

            # 识别基因ID列
            gene_col = None
            for col in df.columns:
                if 'gene' in col.lower() or col == 'Geneid':
                    gene_col = col
                    break

            if gene_col is None:
                logger.error("未找到基因ID列")
                return {}

            # 识别样本列（通常是 BAM 文件名或样本名）
            sample_cols = []
            for col in df.columns:
                if col != gene_col and col not in ['Chr', 'Start', 'End', 'Strand', 'Length']:
                    # 提取样本名（从完整路径或文件名）
                    sample_name = Path(col).stem.replace('.sorted', '').replace('.bam', '')
                    sample_cols.append((col, sample_name))

            logger.info(f"识别到 {len(sample_cols)} 个样本列")

            # 转换为 {样本名: DataFrame} 格式
            count_data = {}
            for col, sample_name in sample_cols:
                df_sample = df[[gene_col, col]].copy()
                df_sample.columns = ['gene_id', 'count']
                df_sample['count'] = df_sample['count'].astype(int)
                count_data[sample_name] = df_sample
                logger.info(f"加载样本: {sample_name}")

            return count_data

        except Exception as e:
            logger.error(f"读取 featureCounts 矩阵失败: {e}")
            return {}

    def _read_count_file(self, count_file: Path) -> Optional[pd.DataFrame]:
        """读取单个计数文件"""
        try:
            # 根据文件扩展名选择分隔符
            if count_file.suffix == '.csv' or count_file.suffix == '.CSV':
                df = pd.read_csv(count_file, sep=',', comment='#', header=0)
            else:
                # 尝试读取featureCounts输出格式（TSV）
                df = pd.read_csv(count_file, sep='\t', comment='#', header=0)

            # 检查必要的列
            required_cols = []
            for col in df.columns:
                if 'gene' in col.lower():
                    required_cols.append(col)
                elif 'feature' in col.lower():
                    required_cols.append(col)

            if not required_cols:
                logger.warning(f"未找到基因ID列: {count_file}")
                return None

            gene_col = required_cols[0]

            # 查找计数列
            count_cols = []
            for col in df.columns:
                if col != gene_col and ('count' in col.lower() or col.endswith('.bam') or col.isdigit()):
                    count_cols.append(col)

            if not count_cols:
                # 如果没有明确的计数列，使用最后一列
                count_cols = [df.columns[-1]]

            # 提取基因ID和计数
            result = df[[gene_col, count_cols[0]]].copy()
            result.columns = ['gene_id', 'count']

            # 确保计数为整数
            result['count'] = result['count'].astype(int)

            return result

        except Exception as e:
            logger.warning(f"读取计数文件失败 {count_file}: {e}")
            return None

    def create_count_matrix(self, count_data: Dict[str, pd.DataFrame],
                           min_count: int = 1, min_samples: int = 1) -> pd.DataFrame:
        """
        创建计数矩阵

        参数:
            count_data: 样本计数数据字典
            min_count: 最小计数阈值
            min_samples: 最小样本数阈值

        返回:
            pd.DataFrame: 计数矩阵
        """
        if not count_data:
            logger.error("无计数数据可处理")
            return pd.DataFrame()

        logger.info(f"创建计数矩阵，样本数: {len(count_data)}")

        # 合并所有样本的计数
        count_frames = []
        for sample, df in count_data.items():
            # 设置基因ID为索引
            df_sample = df.set_index('gene_id')[['count']].copy()
            df_sample.columns = [sample]
            count_frames.append(df_sample)

        # 合并所有样本
        try:
            # 从第一个数据框开始
            count_matrix = count_frames[0]

            for i in range(1, len(count_frames)):
                count_matrix = count_matrix.join(count_frames[i], how='outer')

            # 填充NaN值为0
            count_matrix = count_matrix.fillna(0).astype(int)

            # 过滤低表达基因
            if min_count > 0 or min_samples > 0:
                count_matrix = self._filter_low_expression(count_matrix, min_count, min_samples)

            self.count_matrix = count_matrix
            logger.info(f"计数矩阵创建完成，维度: {count_matrix.shape}")

            return count_matrix

        except Exception as e:
            logger.error(f"创建计数矩阵时出错: {e}")
            return pd.DataFrame()

    def _filter_low_expression(self, count_matrix: pd.DataFrame,
                              min_count: int, min_samples: int) -> pd.DataFrame:
        """过滤低表达基因"""
        original_genes = count_matrix.shape[0]

        # 计算每个基因在多少样本中达到最小计数
        if min_count > 0:
            # 至少达到min_count的样本数
            expressed_samples = (count_matrix >= min_count).sum(axis=1)
            keep_genes = expressed_samples >= min_samples
        else:
            # 只需要在min_samples个样本中有表达（>0）
            expressed_samples = (count_matrix > 0).sum(axis=1)
            keep_genes = expressed_samples >= min_samples

        filtered_matrix = count_matrix[keep_genes]

        filtered_genes = filtered_matrix.shape[0]
        logger.info(f"过滤低表达基因: {original_genes} -> {filtered_genes} "
                   f"(保留 {filtered_genes/original_genes*100:.1f}%)")

        return filtered_matrix

    def normalize_counts(self, method: str = 'cpm') -> pd.DataFrame:
        """
        标准化计数

        参数:
            method: 标准化方法 ('cpm', 'tpm', 'rpkm', 'none')

        返回:
            pd.DataFrame: 标准化后的矩阵
        """
        if self.count_matrix is None or self.count_matrix.empty:
            logger.error("无计数矩阵可标准化")
            return pd.DataFrame()

        logger.info(f"使用 {method} 方法标准化计数")

        if method == 'none' or method == 'raw':
            return self.count_matrix.copy()

        normalized_matrix = self.count_matrix.copy().astype(float)

        if method == 'cpm':
            # Counts Per Million
            for sample in normalized_matrix.columns:
                total_counts = normalized_matrix[sample].sum()
                if total_counts > 0:
                    normalized_matrix[sample] = (normalized_matrix[sample] / total_counts) * 1e6

        elif method == 'tpm':
            # Transcripts Per Million (需要基因长度信息)
            logger.warning("TPM标准化需要基因长度信息，暂不支持，使用CPM代替")
            return self.normalize_counts('cpm')

        elif method == 'rpkm':
            # Reads Per Kilobase per Million (需要基因长度信息)
            logger.warning("RPKM标准化需要基因长度信息，暂不支持，使用CPM代替")
            return self.normalize_counts('cpm')

        else:
            logger.warning(f"未知标准化方法: {method}，使用CPM")
            return self.normalize_counts('cpm')

        return normalized_matrix

    def save_matrix(self, matrix: pd.DataFrame, output_dir: Path,
                   prefix: str = "gene_expression"):
        """
        保存表达矩阵

        参数:
            matrix: 表达矩阵
            output_dir: 输出目录
            prefix: 文件名前缀
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # 保存为CSV
        csv_file = output_dir / f"{prefix}_matrix.csv"
        matrix.to_csv(csv_file)
        logger.info(f"CSV格式矩阵已保存: {csv_file}")

        # 保存为TSV
        tsv_file = output_dir / f"{prefix}_matrix.tsv"
        matrix.to_csv(tsv_file, sep='\t')
        logger.info(f"TSV格式矩阵已保存: {tsv_file}")

        # 尝试保存为Excel
        try:
            excel_file = output_dir / f"{prefix}_matrix.xlsx"
            matrix.to_excel(excel_file)
            logger.info(f"Excel格式矩阵已保存: {excel_file}")
        except ImportError:
            logger.warning("未安装openpyxl，跳过Excel格式保存")
        except Exception as e:
            logger.warning(f"保存Excel格式失败: {e}")

        # 保存基因列表
        gene_list_file = output_dir / f"{prefix}_genes.txt"
        with open(gene_list_file, 'w') as f:
            for gene in matrix.index:
                f.write(f"{gene}\n")
        logger.info(f"基因列表已保存: {gene_list_file}")

    def generate_report(self, count_matrix: pd.DataFrame,
                       normalized_matrix: Optional[pd.DataFrame],
                       output_dir: Path):
        """
        生成表达矩阵报告

        参数:
            count_matrix: 原始计数矩阵
            normalized_matrix: 标准化矩阵
            output_dir: 输出目录
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        report_file = output_dir / "expression_matrix_report.txt"

        try:
            with open(report_file, 'w') as f:
                f.write("=== 基因表达矩阵分析报告 ===\n\n")
                f.write(f"生成时间: {pd.Timestamp.now()}\n\n")

                f.write("1. 数据概览\n")
                f.write(f"   基因数: {count_matrix.shape[0]:,}\n")
                f.write(f"   样本数: {count_matrix.shape[1]:,}\n")
                f.write(f"   总计数: {count_matrix.sum().sum():,}\n\n")

                f.write("2. 样本统计\n")
                sample_stats = count_matrix.sum().sort_values(ascending=False)
                f.write("   样本总计数:\n")
                for sample, total in sample_stats.items():
                    f.write(f"     {sample}: {total:,}\n")

                f.write("\n   样本平均计数:\n")
                sample_means = count_matrix.mean()
                for sample, mean in sample_means.sort_values(ascending=False).items():
                    f.write(f"     {sample}: {mean:.1f}\n")

                f.write("\n3. 基因表达水平分布\n")
                gene_means = count_matrix.mean(axis=1)

                # 表达水平分布
                f.write("   按平均表达水平分布:\n")
                bins = [0, 0.1, 1, 10, 100, 1000, float('inf')]
                bin_labels = ['0', '0.1-1', '1-10', '10-100', '100-1000', '>1000']

                for i in range(len(bins)-1):
                    low = bins[i]
                    high = bins[i+1]
                    if high == float('inf'):
                        count = (gene_means >= low).sum()
                    else:
                        count = ((gene_means >= low) & (gene_means < high)).sum()

                    proportion = count / len(gene_means) * 100
                    f.write(f"     {bin_labels[i]}: {count:,} genes ({proportion:.1f}%)\n")

                f.write(f"\n   零表达基因: {(gene_means == 0).sum():,} "
                       f"({(gene_means == 0).mean()*100:.1f}%)\n")

                f.write(f"   中位表达水平: {gene_means.median():.2f}\n")
                f.write(f"   平均表达水平: {gene_means.mean():.2f}\n")
                f.write(f"   最大表达水平: {gene_means.max():.2f}\n\n")

                f.write("4. 表达相关性（样本间）\n")
                if count_matrix.shape[1] > 1:
                    correlations = count_matrix.corr()
                    f.write("   样本间相关性矩阵:\n")
                    f.write(correlations.to_string())
                    f.write("\n\n")

                if normalized_matrix is not None:
                    f.write("5. 标准化矩阵统计\n")
                    f.write(f"   标准化方法: CPM (Counts Per Million)\n")
                    f.write(f"   标准化后总计数（每百万）: {normalized_matrix.sum().sum()/1e6:.1f} million\n\n")

                f.write("6. 质量控制指标\n")
                # 计算质量控制指标
                f.write(f"   基因检出率（在至少一个样本中表达）: "
                       f"{(count_matrix.sum(axis=1) > 0).mean()*100:.1f}%\n")

                if count_matrix.shape[1] > 1:
                    # 计算样本间相似性
                    sample_similarity = count_matrix.corr().mean().mean()
                    f.write(f"   平均样本间相关性: {sample_similarity:.3f}\n")

                f.write("\n7. 文件列表\n")
                f.write("   gene_expression_matrix.csv - CSV格式表达矩阵\n")
                f.write("   gene_expression_matrix.tsv - TSV格式表达矩阵\n")
                f.write("   gene_expression_genes.txt - 基因列表\n")
                if Path(output_dir / "gene_expression_matrix.xlsx").exists():
                    f.write("   gene_expression_matrix.xlsx - Excel格式表达矩阵\n")

            logger.info(f"表达矩阵报告已生成: {report_file}")

        except Exception as e:
            logger.warning(f"生成报告时出错: {e}")

    def run_analysis(self, input_dir: str, output_dir: str,
                    min_count: int = 1, min_samples: int = 1,
                    normalize: bool = True):
        """
        运行完整分析流程

        参数:
            input_dir: 输入目录或文件（featureCounts批量输出的gene_counts.csv）
            output_dir: 输出目录或CSV文件
            min_count: 最小计数阈值
            min_samples: 最小样本数阈值
            normalize: 是否标准化
        """
        output_path = Path(output_dir)

        # 检查输出是文件还是目录
        is_file_output = output_path.suffix == '.csv'
        if is_file_output:
            # 确保父目录存在
            output_path.parent.mkdir(parents=True, exist_ok=True)
            actual_output_dir = output_path.parent
        else:
            actual_output_dir = output_path
            actual_output_dir.mkdir(parents=True, exist_ok=True)

        # 1. 加载计数文件
        logger.info("步骤1: 加载计数文件")

        # 检查输入是文件还是目录
        input_path = Path(input_dir)
        if input_path.is_file() and input_path.name == "gene_counts.csv":
            # 直接处理文件
            logger.info(f"检测到gene_counts.csv文件: {input_path}")
            count_data = self._load_featurecounts_matrix(input_path)
        else:
            # 按目录处理
            count_data = self.load_count_files(input_dir)

        if not count_data:
            logger.error("未找到有效的计数文件")
            return

        # 2. 创建计数矩阵
        logger.info("步骤2: 创建计数矩阵")
        count_matrix = self.create_count_matrix(count_data, min_count, min_samples)

        if count_matrix.empty:
            logger.error("无法创建计数矩阵")
            return

        # 3. 根据输出类型保存
        logger.info("步骤3: 保存计数矩阵")
        if is_file_output:
            # 直接保存矩阵到指定文件
            count_matrix.to_csv(output_path)
            logger.info(f"计数矩阵已保存: {output_path}")
        else:
            # 保存原始计数矩阵到目录
            self.save_matrix(count_matrix, actual_output_dir, "raw_counts")

            # 4. 标准化计数
            normalized_matrix = None
            if normalize:
                logger.info("步骤4: 标准化计数")
                normalized_matrix = self.normalize_counts('cpm')
                self.save_matrix(normalized_matrix, actual_output_dir, "normalized_counts")

            # 5. 生成报告
                logger.info("步骤5: 生成分析报告")
                self.generate_report(count_matrix, normalized_matrix, actual_output_dir)

                logger.info("表达矩阵分析完成")


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="基因表达矩阵生成脚本"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入目录或文件（featureCounts输出的gene_counts.csv）"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/expression/matrix",
        help="输出目录 (默认: results/expression/matrix)"
    )

    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="最小计数阈值 (默认: 1)"
    )

    parser.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="最小样本数阈值 (默认: 1)"
    )

    parser.add_argument(
        "--no-normalize",
        action="store_true",
        help="不进行计数标准化"
    )

    parser.add_argument(
        "--pattern",
        default="*_counts.txt",
        help="计数文件名模式 (默认: *_counts.txt)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 初始化生成器
    generator = ExpressionMatrixGenerator()

    # 运行分析
    generator.run_analysis(
        input_dir=args.input,
        output_dir=args.output,
        min_count=args.min_count,
        min_samples=args.min_samples,
        normalize=not args.no_normalize
    )

    logger.info("基因表达矩阵生成流程完成")


if __name__ == "__main__":
    main()