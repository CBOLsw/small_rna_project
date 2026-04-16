#!/usr/bin/env python3
"""
差异表达基因筛选和结果导出脚本

功能：
1. 从DESeq2结果中筛选差异表达基因
2. 应用多种筛选标准（fold-change, p-value, 表达水平等）
3. 导出不同格式的结果文件
4. 生成基因列表用于下游分析
5. 创建筛选统计报告

使用方法：
    python filter_degs.py --input <DESeq2结果文件> --output <输出目录>
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

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DEGFilter:
    """差异表达基因筛选器"""

    def __init__(self):
        self.results = None
        self.filtered_results = {}
        self.stats = {}

    def load_results(self, input_file: str) -> pd.DataFrame:
        """
        加载DESeq2结果文件

        参数:
            input_file: 输入文件路径

        返回:
            pd.DataFrame: 结果数据框
        """
        logger.info(f"加载DESeq2结果: {input_file}")

        input_path = Path(input_file)
        if not input_path.exists():
            logger.error(f"输入文件不存在: {input_file}")
            return pd.DataFrame()

        # 根据文件扩展名读取
        file_ext = input_path.suffix.lower()

        try:
            if file_ext == '.csv':
                df = pd.read_csv(input_file)
            elif file_ext in ['.tsv', '.txt']:
                df = pd.read_csv(input_file, sep='\t')
            elif file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(input_file)
            else:
                logger.error(f"不支持的文件格式: {file_ext}")
                return pd.DataFrame()

            logger.info(f"加载成功，数据维度: {df.shape}")
            logger.info(f"列名: {list(df.columns)}")

            # 检查必要的列
            required_cols = ['gene_id', 'log2FoldChange', 'pvalue', 'padj']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                logger.warning(f"缺少必要的列: {missing_cols}")
                logger.warning("尝试识别类似的列...")

                # 尝试识别类似的列名
                col_mapping = {}
                for req_col in missing_cols:
                    for actual_col in df.columns:
                        if req_col.lower() in actual_col.lower():
                            col_mapping[req_col] = actual_col
                            logger.info(f"  映射 {req_col} -> {actual_col}")
                            break

                # 重命名列
                df = df.rename(columns=col_mapping)

            self.results = df
            return df

        except Exception as e:
            logger.error(f"加载结果文件时出错: {e}")
            return pd.DataFrame()

    def filter_degs(self, df: pd.DataFrame,
                   log2fc_threshold: float = 1.0,
                   padj_threshold: float = 0.05,
                   pvalue_threshold: Optional[float] = None,
                   base_mean_threshold: float = 10.0,
                   min_fold_change: float = 1.5,
                   direction: str = 'both') -> pd.DataFrame:
        """
        筛选差异表达基因

        参数:
            df: 输入数据框
            log2fc_threshold: log2倍数变化阈值（绝对值）
            padj_threshold: 调整后p值阈值
            pvalue_threshold: 原始p值阈值（可选）
            base_mean_threshold: 基础表达水平阈值
            min_fold_change: 最小倍数变化阈值（原始值）
            direction: 筛选方向 ('both', 'up', 'down')

        返回:
            pd.DataFrame: 筛选后的数据框
        """
        if df.empty:
            logger.error("输入数据框为空")
            return pd.DataFrame()

        logger.info(f"筛选差异表达基因 (log2FC≥{log2fc_threshold}, padj<{padj_threshold})")

        # 创建副本
        filtered = df.copy()

        # 计算绝对值
        filtered['abs_log2fc'] = filtered['log2FoldChange'].abs()

        # 应用筛选条件
        conditions = []

        # 1. 倍数变化筛选
        if log2fc_threshold > 0:
            conditions.append(filtered['abs_log2fc'] >= log2fc_threshold)

        # 2. 显著性筛选
        if padj_threshold:
            conditions.append(filtered['padj'] < padj_threshold)

        # 3. 原始p值筛选（可选）
        if pvalue_threshold and 'pvalue' in filtered.columns:
            conditions.append(filtered['pvalue'] < pvalue_threshold)

        # 4. 表达水平筛选
        if base_mean_threshold > 0 and 'baseMean' in filtered.columns:
            conditions.append(filtered['baseMean'] >= base_mean_threshold)

        # 5. 原始倍数变化筛选
        if min_fold_change > 1:
            # 计算原始倍数变化
            filtered['fold_change'] = 2 ** filtered['log2FoldChange'].abs()
            conditions.append(filtered['fold_change'] >= min_fold_change)

        # 组合所有条件
        if conditions:
            mask = conditions[0]
            for condition in conditions[1:]:
                mask = mask & condition
        else:
            mask = pd.Series([True] * len(filtered), index=filtered.index)

        # 按方向筛选
        if direction == 'up':
            mask = mask & (filtered['log2FoldChange'] > 0)
        elif direction == 'down':
            mask = mask & (filtered['log2FoldChange'] < 0)
        # 'both' 不额外筛选

        # 应用筛选
        degs = filtered[mask].copy()

        # 添加上下调标签
        degs['regulation'] = degs['log2FoldChange'].apply(
            lambda x: 'up' if x > 0 else 'down'
        )

        # 排序：先按显著性，再按倍数变化
        sort_columns = []
        if 'padj' in degs.columns:
            sort_columns.append('padj')
        elif 'pvalue' in degs.columns:
            sort_columns.append('pvalue')
        sort_columns.append('abs_log2fc')

        degs = degs.sort_values(by=sort_columns, ascending=[True, False])

        # 统计信息
        total_genes = len(df)
        deg_count = len(degs)
        up_count = (degs['regulation'] == 'up').sum() if not degs.empty else 0
        down_count = (degs['regulation'] == 'down').sum() if not degs.empty else 0

        logger.info(f"筛选结果: {deg_count}/{total_genes} 基因 ({deg_count/total_genes*100:.1f}%)")
        logger.info(f"  上调基因: {up_count} ({up_count/deg_count*100:.1f}% of DEGs)")
        logger.info(f"  下调基因: {down_count} ({down_count/deg_count*100:.1f}% of DEGs)")

        # 保存统计信息
        self.stats['filtering'] = {
            'total_genes': total_genes,
            'deg_count': deg_count,
            'up_count': up_count,
            'down_count': down_count,
            'deg_percentage': deg_count / total_genes * 100 if total_genes > 0 else 0,
            'parameters': {
                'log2fc_threshold': log2fc_threshold,
                'padj_threshold': padj_threshold,
                'pvalue_threshold': pvalue_threshold,
                'base_mean_threshold': base_mean_threshold,
                'min_fold_change': min_fold_change,
                'direction': direction
            }
        }

        return degs

    def apply_multiple_criteria(self, df: pd.DataFrame,
                              criteria_list: List[Dict[str, Any]]) -> Dict[str, pd.DataFrame]:
        """
        应用多组筛选标准

        参数:
            df: 输入数据框
            criteria_list: 筛选标准列表

        返回:
            Dict: 标准名称 -> 筛选结果
        """
        results = {}

        for i, criteria in enumerate(criteria_list):
            name = criteria.get('name', f'criteria_{i+1}')

            logger.info(f"应用筛选标准: {name}")

            # 提取参数
            log2fc = criteria.get('log2fc_threshold', 1.0)
            padj = criteria.get('padj_threshold', 0.05)
            pvalue = criteria.get('pvalue_threshold', None)
            base_mean = criteria.get('base_mean_threshold', 10.0)
            min_fc = criteria.get('min_fold_change', 1.5)
            direction = criteria.get('direction', 'both')

            # 筛选
            filtered = self.filter_degs(
                df=df,
                log2fc_threshold=log2fc,
                padj_threshold=padj,
                pvalue_threshold=pvalue,
                base_mean_threshold=base_mean,
                min_fold_change=min_fc,
                direction=direction
            )

            results[name] = filtered

        self.filtered_results = results
        return results

    def export_results(self, results: Dict[str, pd.DataFrame],
                      output_dir: str, prefix: str = "degs"):
        """
        导出筛选结果

        参数:
            results: 筛选结果字典
            output_dir: 输出目录
            prefix: 文件名前缀
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        exported_files = []

        for name, df in results.items():
            if df.empty:
                logger.warning(f"无结果可导出: {name}")
                continue

            # 清理名称用于文件名
            safe_name = name.replace(' ', '_').replace('/', '_').lower()

            # 导出完整结果
            full_file = output_path / f"{prefix}_{safe_name}_full.csv"
            df.to_csv(full_file, index=False)
            exported_files.append(str(full_file))

            # 导出简化版本（仅关键列）
            simple_cols = ['gene_id', 'log2FoldChange', 'padj', 'regulation']
            available_cols = [col for col in simple_cols if col in df.columns]

            if available_cols:
                simple_df = df[available_cols].copy()
                simple_file = output_path / f"{prefix}_{safe_name}_simple.csv"
                simple_df.to_csv(simple_file, index=False)
                exported_files.append(str(simple_file))

            # 导出基因列表
            genes_file = output_path / f"{prefix}_{safe_name}_genes.txt"
            with open(genes_file, 'w') as f:
                for gene in df['gene_id']:
                    f.write(f"{gene}\n")
            exported_files.append(str(genes_file))

            # 分上下调导出
            if 'regulation' in df.columns:
                # 上调基因
                up_df = df[df['regulation'] == 'up']
                if not up_df.empty:
                    up_file = output_path / f"{prefix}_{safe_name}_up.csv"
                    up_df.to_csv(up_file, index=False)
                    exported_files.append(str(up_file))

                    up_genes_file = output_path / f"{prefix}_{safe_name}_up_genes.txt"
                    with open(up_genes_file, 'w') as f:
                        for gene in up_df['gene_id']:
                            f.write(f"{gene}\n")
                    exported_files.append(str(up_genes_file))

                # 下调基因
                down_df = df[df['regulation'] == 'down']
                if not down_df.empty:
                    down_file = output_path / f"{prefix}_{safe_name}_down.csv"
                    down_df.to_csv(down_file, index=False)
                    exported_files.append(str(down_file))

                    down_genes_file = output_path / f"{prefix}_{safe_name}_down_genes.txt"
                    with open(down_genes_file, 'w') as f:
                        for gene in down_df['gene_id']:
                            f.write(f"{gene}\n")
                    exported_files.append(str(down_genes_file))

            logger.info(f"导出结果: {name} ({len(df)} 基因)")

        # 生成文件清单
        manifest_file = output_path / f"{prefix}_file_manifest.txt"
        with open(manifest_file, 'w') as f:
            f.write("=== 差异表达基因文件清单 ===\n\n")
            f.write(f"生成时间: {pd.Timestamp.now()}\n\n")
            f.write("文件列表:\n")
            for file_path in exported_files:
                f.write(f"  {Path(file_path).name}\n")
        exported_files.append(str(manifest_file))

        logger.info(f"所有结果已导出到: {output_dir}")
        return exported_files

    def generate_summary_report(self, results: Dict[str, pd.DataFrame],
                               output_dir: str, prefix: str = "degs"):
        """
        生成筛选结果汇总报告

        参数:
            results: 筛选结果字典
            output_dir: 输出目录
            prefix: 文件名前缀
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        report_file = output_path / f"{prefix}_filtering_report.txt"

        try:
            with open(report_file, 'w') as f:
                f.write("=== 差异表达基因筛选报告 ===\n\n")
                f.write(f"生成时间: {pd.Timestamp.now()}\n\n")

                f.write("1. 筛选标准汇总\n")
                for name, df in results.items():
                    f.write(f"\n  标准: {name}\n")
                    f.write(f"    总基因数: {len(self.results) if self.results is not None else 'N/A'}\n")
                    f.write(f"    差异基因数: {len(df)}\n")
                    if len(df) > 0:
                        up_count = (df['regulation'] == 'up').sum()
                        down_count = (df['regulation'] == 'down').sum()
                        f.write(f"      上调基因: {up_count}\n")
                        f.write(f"      下调基因: {down_count}\n")
                        f.write(f"    差异基因比例: {len(df)/len(self.results)*100:.1f}%\n")

                f.write("\n2. 基因重叠分析\n")
                if len(results) > 1:
                    # 计算不同标准间的重叠
                    gene_sets = {}
                    for name, df in results.items():
                        gene_sets[name] = set(df['gene_id'])

                    f.write("   不同标准间的基因重叠:\n")
                    names = list(gene_sets.keys())
                    for i, name1 in enumerate(names):
                        for name2 in names[i+1:]:
                            overlap = gene_sets[name1] & gene_sets[name2]
                            union = gene_sets[name1] | gene_sets[name2]
                            overlap_pct = len(overlap) / len(union) * 100 if union else 0

                            f.write(f"     {name1} ∩ {name2}: {len(overlap)} 基因 ({overlap_pct:.1f}%)\n")

                f.write("\n3. 表达水平统计\n")
                if self.results is not None and not self.results.empty:
                    f.write("   所有基因:\n")
                    if 'baseMean' in self.results.columns:
                        f.write(f"     平均表达: {self.results['baseMean'].mean():.2f}\n")
                        f.write(f"     中位表达: {self.results['baseMean'].median():.2f}\n")

                    if 'log2FoldChange' in self.results.columns:
                        f.write(f"     平均log2FC: {self.results['log2FoldChange'].mean():.2f}\n")
                        f.write(f"     log2FC范围: [{self.results['log2FoldChange'].min():.2f}, "
                               f"{self.results['log2FoldChange'].max():.2f}]\n")

                f.write("\n4. 筛选参数\n")
                if 'filtering' in self.stats:
                    params = self.stats['filtering'].get('parameters', {})
                    for key, value in params.items():
                        f.write(f"    {key}: {value}\n")

                f.write("\n5. 文件说明\n")
                f.write("   *_full.csv - 完整筛选结果（包含所有列）\n")
                f.write("   *_simple.csv - 简化结果（仅关键列）\n")
                f.write("   *_genes.txt - 基因列表（每行一个基因）\n")
                f.write("   *_up.csv - 上调基因完整结果\n")
                f.write("   *_up_genes.txt - 上调基因列表\n")
                f.write("   *_down.csv - 下调基因完整结果\n")
                f.write("   *_down_genes.txt - 下调基因列表\n")
                f.write("   *_file_manifest.txt - 文件清单\n")

            logger.info(f"筛选报告已生成: {report_file}")

        except Exception as e:
            logger.warning(f"生成报告时出错: {e}")

    def run_analysis(self, input_file: str, output_dir: str,
                    criteria_list: Optional[List[Dict[str, Any]]] = None):
        """
        运行完整筛选分析

        参数:
            input_file: 输入文件路径
            output_dir: 输出目录
            criteria_list: 筛选标准列表
        """
        # 1. 加载结果
        df = self.load_results(input_file)
        if df.empty:
            logger.error("无法加载结果文件")
            return

        # 2. 应用筛选标准
        if criteria_list is None:
            # 默认筛选标准
            criteria_list = [
                {
                    'name': 'stringent',
                    'log2fc_threshold': 1.0,  # |log2FC| >= 1
                    'padj_threshold': 0.01,
                    'base_mean_threshold': 10.0
                },
                {
                    'name': 'moderate',
                    'log2fc_threshold': 0.5,  # |log2FC| >= 0.5
                    'padj_threshold': 0.05,
                    'base_mean_threshold': 5.0
                },
                {
                    'name': 'lenient',
                    'log2fc_threshold': 0.0,  # 仅基于显著性
                    'padj_threshold': 0.1,
                    'base_mean_threshold': 0.0
                }
            ]

        filtered_results = self.apply_multiple_criteria(df, criteria_list)

        # 3. 导出结果
        exported_files = self.export_results(filtered_results, output_dir)

        # 4. 生成报告
        self.generate_summary_report(filtered_results, output_dir)

        logger.info(f"差异表达基因筛选完成，文件数: {len(exported_files)}")


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """加载配置文件"""
    default_config = {
        'default_criteria': [
            {
                'name': 'stringent',
                'log2fc_threshold': 1.0,
                'padj_threshold': 0.01,
                'base_mean_threshold': 10.0
            },
            {
                'name': 'moderate',
                'log2fc_threshold': 0.5,
                'padj_threshold': 0.05,
                'base_mean_threshold': 5.0
            }
        ],
        'output_prefix': 'degs',
        'export_formats': ['csv', 'txt', 'simple']
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
        description="差异表达基因筛选和结果导出脚本"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="DESeq2结果文件（CSV/TSV格式）"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/expression/degs",
        help="输出目录 (默认: results/expression/degs)"
    )

    parser.add_argument(
        "--config", "-c",
        help="配置文件 (YAML格式)"
    )

    parser.add_argument(
        "--log2fc-threshold",
        type=float,
        default=1.0,
        help="log2倍数变化阈值 (默认: 1.0)"
    )

    parser.add_argument(
        "--padj-threshold",
        type=float,
        default=0.05,
        help="调整后p值阈值 (默认: 0.05)"
    )

    parser.add_argument(
        "--base-mean-threshold",
        type=float,
        default=10.0,
        help="基础表达水平阈值 (默认: 10.0)"
    )

    parser.add_argument(
        "--direction",
        choices=['both', 'up', 'down'],
        default='both',
        help="筛选方向 (默认: both)"
    )

    parser.add_argument(
        "--single-criteria",
        action="store_true",
        help="使用单一筛选标准（而非多组标准）"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 加载配置
    config = load_config(args.config)

    # 初始化筛选器
    filter_tool = DEGFilter()

    if args.single_criteria:
        # 使用单一筛选标准
        logger.info("使用单一筛选标准")

        # 加载结果
        df = filter_tool.load_results(args.input)
        if df.empty:
            logger.error("无法加载结果文件")
            sys.exit(1)

        # 筛选
        degs = filter_tool.filter_degs(
            df=df,
            log2fc_threshold=args.log2fc_threshold,
            padj_threshold=args.padj_threshold,
            base_mean_threshold=args.base_mean_threshold,
            direction=args.direction
        )

        # 导出
        results = {'single_criteria': degs}
        filter_tool.export_results(results, args.output)
        filter_tool.generate_summary_report(results, args.output)

    else:
        # 使用多组筛选标准
        logger.info("使用多组筛选标准")

        # 运行完整分析
        filter_tool.run_analysis(
            input_file=args.input,
            output_dir=args.output,
            criteria_list=config.get('default_criteria')
        )

    logger.info("差异表达基因筛选流程完成")


if __name__ == "__main__":
    main()