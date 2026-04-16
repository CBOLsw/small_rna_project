#!/usr/bin/env python3
"""
featureCounts基因计数脚本

功能：
1. 使用featureCounts从BAM文件计算基因计数
2. 支持small RNA特征计数（miRNA, snoRNA等）
3. 生成每个样本的计数文件
4. 汇总所有样本的计数矩阵

使用方法：
    python count_features.py --input <BAM文件或目录> --annotation <GTF/GFF文件> --output <输出目录>
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

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class FeatureCounter:
    """featureCounts基因计数类"""

    def __init__(self, featurecounts_path: str = "featureCounts"):
        """
        初始化基因计数器

        参数:
            featurecounts_path: featureCounts可执行文件路径
        """
        self.featurecounts_path = featurecounts_path
        self.results = {}

    def check_featurecounts(self) -> bool:
        """检查featureCounts是否可用"""
        try:
            result = subprocess.run(
                [self.featurecounts_path, "-v"],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                version_line = result.stdout.strip()
                logger.info(f"featureCounts版本: {version_line}")
                return True
            else:
                logger.error(f"featureCounts检查失败: {result.stderr}")
                return False
        except FileNotFoundError:
            logger.error(f"未找到featureCounts: {self.featurecounts_path}")
            return False

    def count_single_bam(self, bam_file: str, annotation_file: str,
                        output_dir: str, sample_name: Optional[str] = None,
                        config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        对单个BAM文件进行基因计数

        参数:
            bam_file: BAM文件路径
            annotation_file: 基因注释文件（GTF/GFF）
            output_dir: 输出目录
            sample_name: 样本名称
            config: 配置参数

        返回:
            Dict: 计数结果
        """
        if config is None:
            config = {}

        bam_path = Path(bam_file)
        if not bam_path.exists():
            logger.error(f"BAM文件不存在: {bam_file}")
            return {'success': False, 'error': '文件不存在'}

        annotation_path = Path(annotation_file)
        if not annotation_path.exists():
            logger.error(f"注释文件不存在: {annotation_file}")
            return {'success': False, 'error': '注释文件不存在'}

        if sample_name is None:
            sample_name = bam_path.stem
            if sample_name.endswith('_sorted'):
                sample_name = sample_name[:-7]

        logger.info(f"基因计数: {sample_name}")

        # 准备输出目录
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 输出文件
        count_file = output_dir / f"{sample_name}_counts.txt"
        summary_file = output_dir / f"{sample_name}_counts_summary.txt"
        log_file = output_dir / f"{sample_name}_featurecounts.log"

        # 构建featureCounts命令
        cmd = self._build_featurecounts_command(
            bam_file, annotation_file, count_file, summary_file, config
        )

        # 运行featureCounts
        success, count_stats = self._run_featurecounts(cmd, log_file, sample_name)

        result = {
            'sample': sample_name,
            'bam_file': bam_file,
            'annotation_file': annotation_file,
            'count_file': str(count_file) if success else None,
            'summary_file': str(summary_file) if success else None,
            'log_file': str(log_file),
            'success': success
        }

        if success:
            # 解析计数结果
            count_data = self._parse_count_file(count_file, sample_name)
            summary_data = self._parse_summary_file(summary_file, sample_name)

            result.update({
                'total_features': count_data.get('total_features', 0),
                'assigned_reads': summary_data.get('assigned', 0),
                'unassigned_reads': summary_data.get('unassigned', 0),
                'assignment_rate': summary_data.get('assignment_rate', 0),
                'count_data': count_data
            })

            logger.info(f"基因计数完成: {sample_name}, "
                       f"分配reads: {summary_data.get('assigned', 0):,}, "
                       f"分配率: {summary_data.get('assignment_rate', 0)*100:.1f}%")

        self.results[sample_name] = result
        return result

    def _build_featurecounts_command(self, bam_file: str, annotation_file: str,
                                   count_file: Path, summary_file: Path,
                                   config: Dict[str, Any]) -> List[str]:
        """构建featureCounts命令"""
        cmd = [
            self.featurecounts_path,
            "-a", annotation_file,
            "-o", str(count_file),
            "-T", str(config.get('threads', 4)),
        ]

        # 输入文件
        cmd.append(bam_file)

        # 计数参数
        if config.get('feature_type'):
            cmd.extend(["-t", config['feature_type']])
        else:
            cmd.extend(["-t", "exon"])  # 默认计数exon

        if config.get('attribute'):
            cmd.extend(["-g", config['attribute']])
        else:
            cmd.extend(["-g", "gene_id"])  # 默认使用gene_id

        # small RNA特定参数
        if config.get('small_rna_mode', True):
            # small RNA通常较短，需要调整参数
            cmd.extend([
                "-O",  # 允许reads分配到多个特征
                "-M",  # 允许多重比对reads
                "--primary",  # 只计数主要比对
                "-s", str(config.get('strand_specific', 0)),  # 链特异性
            ])

            # 对于small RNA，通常计数整个基因而不是exon
            cmd.extend(["-t", "gene"])
            cmd.extend(["-g", "gene_name"])

        # 其他参数
        if config.get('min_overlap'):
            cmd.extend(["--minOverlap", str(config['min_overlap'])])

        if config.get('frac_overlap'):
            cmd.extend(["--fracOverlap", str(config['frac_overlap'])])

        if config.get('frac_overlap_feature'):
            cmd.extend(["--fracOverlapFeature", str(config['frac_overlap_feature'])])

        if config.get('largest_overlap'):
            cmd.append("--largestOverlap")

        if config.get('read_extension'):
            cmd.extend(["-E", str(config['read_extension'])])

        # 输出选项
        cmd.extend([
            "-R", "BAM",  # 输出每个read的分配信息到BAM文件
            "--verbose",
        ])

        return cmd

    def _run_featurecounts(self, cmd: List[str], log_file: Path,
                          sample_name: str) -> Tuple[bool, str]:
        """运行featureCounts"""
        logger.info(f"运行featureCounts: {sample_name}")
        logger.debug(f"命令: {' '.join(cmd)}")

        try:
            with open(log_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )

            if result.returncode == 0:
                logger.info(f"featureCounts运行成功: {sample_name}")
                # 读取日志文件获取统计信息
                with open(log_file, 'r') as f:
                    log_content = f.read()
                return True, log_content
            else:
                logger.error(f"featureCounts运行失败: {sample_name}")
                logger.error(f"返回码: {result.returncode}")
                # 读取日志文件获取错误信息
                with open(log_file, 'r') as f:
                    log_content = f.read()
                logger.error(f"错误信息: {log_content[:500]}...")
                return False, log_content

        except Exception as e:
            logger.error(f"运行featureCounts时出错: {e}")
            return False, ""

    def _parse_count_file(self, count_file: Path, sample_name: str) -> Dict[str, Any]:
        """解析计数文件"""
        try:
            # featureCounts输出格式：第一行是header，然后是数据
            df = pd.read_csv(count_file, sep='\t', comment='#', header=0)

            # 获取计数列（通常是最后一列）
            count_columns = [col for col in df.columns if 'count' in col.lower() or col.endswith('.bam')]
            if count_columns:
                count_col = count_columns[0]
                counts = df[count_col].astype(int)
            else:
                # 如果没有找到计数列，使用最后一列
                count_col = df.columns[-1]
                counts = df[count_col].astype(int)

            # 计算统计
            total_features = len(df)
            total_counts = counts.sum()
            expressed_features = (counts > 0).sum()
            zero_count_features = (counts == 0).sum()

            return {
                'total_features': total_features,
                'total_counts': total_counts,
                'expressed_features': expressed_features,
                'expression_rate': expressed_features / total_features if total_features > 0 else 0,
                'zero_count_features': zero_count_features,
                'mean_count': counts.mean() if total_features > 0 else 0,
                'median_count': counts.median() if total_features > 0 else 0,
                'max_count': counts.max() if total_features > 0 else 0,
                'dataframe_shape': df.shape,
                'count_column': count_col
            }

        except Exception as e:
            logger.warning(f"解析计数文件时出错: {e}")
            return {}

    def _parse_summary_file(self, summary_file: Path, sample_name: str) -> Dict[str, Any]:
        """解析汇总文件"""
        try:
            if not summary_file.exists():
                # 如果没有单独的汇总文件，尝试从计数文件中提取
                return self._parse_summary_from_count_file(summary_file.parent / f"{sample_name}_counts.txt.summary")

            df = pd.read_csv(summary_file, sep='\t', header=0, index_col=0)

            # 获取统计信息
            assigned = df.loc['Assigned', 'Count'] if 'Assigned' in df.index else 0
            unassigned = df.loc['Unassigned', 'Count'] if 'Unassigned' in df.index else 0
            total = assigned + unassigned if assigned + unassigned > 0 else 1

            return {
                'assigned': int(assigned),
                'unassigned': int(unassigned),
                'total': int(total),
                'assignment_rate': assigned / total if total > 0 else 0,
                'summary_data': df.to_dict()
            }

        except Exception as e:
            logger.warning(f"解析汇总文件时出错: {e}")
            return {}

    def _parse_summary_from_count_file(self, count_file: Path) -> Dict[str, Any]:
        """从计数文件解析汇总信息（备用方法）"""
        try:
            # 读取计数文件的header部分
            with open(count_file, 'r') as f:
                lines = f.readlines()

            assigned = 0
            unassigned = 0

            for line in lines:
                if line.startswith('#') and 'Assigned' in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            assigned = int(parts[-1])
                        except ValueError:
                            pass
                elif line.startswith('#') and 'Unassigned' in line:
                    parts = line.strip().split()
                    if len(parts) >= 2:
                        try:
                            unassigned = int(parts[-1])
                        except ValueError:
                            pass

            total = assigned + unassigned if assigned + unassigned > 0 else 1

            return {
                'assigned': assigned,
                'unassigned': unassigned,
                'total': total,
                'assignment_rate': assigned / total if total > 0 else 0
            }

        except Exception as e:
            logger.warning(f"从计数文件解析汇总信息时出错: {e}")
            return {'assigned': 0, 'unassigned': 0, 'total': 0, 'assignment_rate': 0}

    def generate_count_matrix(self, output_dir: str) -> Optional[str]:
        """
        生成计数矩阵（所有样本的合并计数表）

        参数:
            output_dir: 输出目录

        返回:
            str: 计数矩阵文件路径
        """
        if not self.results:
            logger.warning("无计数结果可合并")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 收集所有样本的计数数据
        count_dataframes = []
        sample_names = []

        for sample, result in self.results.items():
            if result.get('success') and result.get('count_file'):
                count_file = result['count_file']
                try:
                    df = pd.read_csv(count_file, sep='\t', comment='#', header=0)

                    # 获取基因ID和计数列
                    gene_id_col = None
                    count_col = None

                    # 寻找基因ID列
                    for col in df.columns:
                        if 'gene' in col.lower() or 'feature' in col.lower():
                            gene_id_col = col
                            break

                    # 寻找计数列
                    for col in df.columns:
                        if 'count' in col.lower() or col.endswith('.bam'):
                            count_col = col
                            break

                    if gene_id_col and count_col:
                        # 提取基因ID和计数
                        counts = df[[gene_id_col, count_col]].copy()
                        counts.columns = ['gene_id', sample]

                        count_dataframes.append(counts.set_index('gene_id'))
                        sample_names.append(sample)

                except Exception as e:
                    logger.warning(f"读取计数文件失败 ({sample}): {e}")

        if not count_dataframes:
            logger.error("无法读取任何计数文件")
            return None

        # 合并所有样本的计数
        try:
            # 从第一个数据框开始合并
            count_matrix = count_dataframes[0]

            for i in range(1, len(count_dataframes)):
                count_matrix = count_matrix.join(count_dataframes[i], how='outer')

            # 填充NaN值为0
            count_matrix = count_matrix.fillna(0).astype(int)

            # 保存计数矩阵
            matrix_file = output_dir / "gene_count_matrix.csv"
            count_matrix.to_csv(matrix_file)

            # 保存为TSV格式
            tsv_file = output_dir / "gene_count_matrix.tsv"
            count_matrix.to_csv(tsv_file, sep='\t')

            # 保存为Excel格式（可选）
            try:
                excel_file = output_dir / "gene_count_matrix.xlsx"
                count_matrix.to_excel(excel_file)
            except ImportError:
                logger.warning("未安装openpyxl，跳过Excel格式导出")

            # 生成统计报告
            self._generate_matrix_report(count_matrix, output_dir)

            logger.info(f"计数矩阵已生成: {matrix_file}")
            logger.info(f"矩阵维度: {count_matrix.shape}")
            logger.info(f"包含样本: {', '.join(sample_names)}")

            return str(matrix_file)

        except Exception as e:
            logger.error(f"生成计数矩阵时出错: {e}")
            return None

    def _generate_matrix_report(self, count_matrix: pd.DataFrame, output_dir: Path):
        """生成计数矩阵统计报告"""
        try:
            report_file = output_dir / "count_matrix_report.txt"

            with open(report_file, 'w') as f:
                f.write("=== 基因计数矩阵统计报告 ===\n\n")
                f.write(f"生成时间: {pd.Timestamp.now()}\n\n")

                f.write("基本统计:\n")
                f.write(f"  基因数: {count_matrix.shape[0]:,}\n")
                f.write(f"  样本数: {count_matrix.shape[1]:,}\n")
                f.write(f"  总计数: {count_matrix.sum().sum():,}\n\n")

                f.write("样本统计:\n")
                sample_stats = count_matrix.sum().sort_values(ascending=False)
                for sample, total in sample_stats.items():
                    f.write(f"  {sample}: {total:,} counts\n")

                f.write("\n基因表达水平统计:\n")
                # 计算每个基因在所有样本中的平均表达
                gene_means = count_matrix.mean(axis=1)
                f.write(f"  平均表达 (所有基因): {gene_means.mean():.2f}\n")
                f.write(f"  中位表达 (所有基因): {gene_means.median():.2f}\n")
                f.write(f"  最大表达: {gene_means.max():.2f}\n")
                f.write(f"  零表达基因数: {(gene_means == 0).sum():,}\n")
                f.write(f"  零表达比例: {(gene_means == 0).mean()*100:.1f}%\n\n")

                f.write("表达水平分布:\n")
                bins = [0, 1, 10, 100, 1000, 10000, float('inf')]
                bin_labels = ['0', '1-10', '11-100', '101-1000', '1001-10000', '>10000']

                for i in range(len(bins)-1):
                    low = bins[i]
                    high = bins[i+1]
                    if high == float('inf'):
                        count = (gene_means >= low).sum()
                    else:
                        count = ((gene_means >= low) & (gene_means < high)).sum()

                    proportion = count / len(gene_means) * 100
                    f.write(f"  {bin_labels[i]}: {count:,} genes ({proportion:.1f}%)\n")

        except Exception as e:
            logger.warning(f"生成矩阵报告时出错: {e}")

    def generate_summary_report(self, output_dir: str) -> Optional[str]:
        """
        生成基因计数汇总报告

        参数:
            output_dir: 输出目录

        返回:
            str: 报告文件路径
        """
        if not self.results:
            logger.warning("无基因计数结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建汇总数据
        summary_data = []
        for sample, result in self.results.items():
            if result.get('success'):
                summary_data.append({
                    'sample': sample,
                    'bam_file': result.get('bam_file', ''),
                    'total_features': result.get('total_features', 0),
                    'assigned_reads': result.get('assigned_reads', 0),
                    'unassigned_reads': result.get('unassigned_reads', 0),
                    'assignment_rate': result.get('assignment_rate', 0),
                    'success': result.get('success', False)
                })

        df = pd.DataFrame(summary_data)

        # 保存CSV
        csv_file = output_dir / "featurecounts_summary.csv"
        df.to_csv(csv_file, index=False)

        # 生成文本报告
        report_file = output_dir / "featurecounts_summary_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== featureCounts基因计数汇总报告 ===\n\n")
            f.write(f"总样本数: {len(summary_data)}\n")
            f.write(f"成功计数: {len([r for r in summary_data if r['success']])}\n\n")

            f.write("样本计数统计:\n")
            for _, row in df.iterrows():
                if row['success']:
                    f.write(f"  - {row['sample']}:\n")
                    f.write(f"    分配reads: {row['assigned_reads']:,}\n")
                    f.write(f"    分配率: {row['assignment_rate']*100:.1f}%\n")
                    f.write(f"    特征数: {row['total_features']:,}\n\n")

            if 'assignment_rate' in df.columns and not df[df['success']].empty:
                avg_rate = df[df['success']]['assignment_rate'].mean() * 100
                f.write(f"平均分配率: {avg_rate:.1f}%\n")

            f.write(f"\n报告生成时间: {pd.Timestamp.now()}\n")

        logger.info(f"基因计数汇总报告已生成: {report_file}")
        return str(report_file)


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """加载配置文件"""
    default_config = {
        'featurecounts_path': 'featureCounts',
        'threads': 4,
        'small_rna_mode': True,
        'feature_type': 'gene',  # small RNA通常计数基因
        'attribute': 'gene_name',
        'strand_specific': 0,  # 0=无链特异性, 1=正链, 2=负链
        'min_overlap': 1,
        'frac_overlap': 0.0,
        'frac_overlap_feature': 0.0,
        'largest_overlap': True,
        'read_extension': 0
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
        description="featureCounts基因计数脚本"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入BAM文件或包含BAM文件的目录"
    )

    parser.add_argument(
        "--annotation", "-a",
        required=True,
        help="基因注释文件（GTF/GFF格式）"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/expression/counts",
        help="输出目录 (默认: results/expression/counts)"
    )

    parser.add_argument(
        "--config", "-c",
        help="配置文件 (YAML格式)"
    )

    parser.add_argument(
        "--sample-info",
        help="样本信息CSV文件，包含sample和bam_file列"
    )

    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=4,
        help="线程数 (默认: 4)"
    )

    parser.add_argument(
        "--matrix",
        action="store_true",
        default=True,
        help="生成计数矩阵 (默认: True)"
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
    config['threads'] = args.threads

    # 初始化基因计数器
    counter = FeatureCounter(featurecounts_path=config.get('featurecounts_path', 'featureCounts'))

    # 检查featureCounts
    if not counter.check_featurecounts():
        logger.error("featureCounts检查失败，请确保featureCounts已安装")
        sys.exit(1)

    # 处理输入
    input_path = Path(args.input)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 根据输入类型处理
    if args.sample_info:
        # 使用样本信息CSV文件
        logger.info(f"使用样本信息文件: {args.sample_info}")
        df = pd.read_csv(args.sample_info)

        for _, row in df.iterrows():
            sample = row.get('sample', f"sample_{_}")
            bam_file = row.get('bam_file')

            if pd.notna(bam_file):
                counter.count_single_bam(
                    bam_file=bam_file,
                    annotation_file=args.annotation,
                    output_dir=output_dir,
                    sample_name=sample,
                    config=config
                )

    elif input_path.is_file() and input_path.suffix == '.bam':
        # 单个BAM文件
        counter.count_single_bam(
            bam_file=str(input_path),
            annotation_file=args.annotation,
            output_dir=output_dir,
            sample_name=None,
            config=config
        )

    elif input_path.is_dir():
        # 目录下的所有BAM文件
        logger.info(f"扫描目录: {input_path}")
        bam_files = list(input_path.glob("*.bam")) + list(input_path.glob("*.BAM"))

        for bam_file in bam_files:
            counter.count_single_bam(
                bam_file=str(bam_file),
                annotation_file=args.annotation,
                output_dir=output_dir,
                sample_name=None,
                config=config
            )

    else:
        logger.error("不支持的输入类型，请提供BAM文件、包含BAM文件的目录或样本信息CSV")
        sys.exit(1)

    # 生成计数矩阵
    if args.matrix:
        matrix_file = counter.generate_count_matrix(output_dir)
        if matrix_file:
            logger.info(f"计数矩阵已生成: {matrix_file}")
        else:
            logger.warning("未能生成计数矩阵")

    # 生成汇总报告
    if args.summary:
        report_file = counter.generate_summary_report(output_dir)
        if report_file:
            logger.info(f"基因计数完成，报告文件: {report_file}")
        else:
            logger.warning("未能生成汇总报告")

    logger.info("featureCounts基因计数流程完成")


if __name__ == "__main__":
    main()