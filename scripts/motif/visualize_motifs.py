#!/usr/bin/env python3
"""
motif可视化脚本
用于生成序列logo图、motif位置分布图等可视化图表

功能：
1. 生成序列logo图
2. 创建motif位置分布图
3. 生成motif热图
4. 可视化motif比较结果
5. 导出高质量图片

使用方法：
    python visualize_motifs.py --meme <MEME结果目录> --output <输出目录> [--format <图片格式>]
"""

import os
import sys
import argparse
import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional, Any, Tuple
import logging

# 尝试导入可视化库
try:
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.patches import Rectangle
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("警告: matplotlib未安装，无法生成可视化图表")

try:
    import seaborn as sns
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class MotifVisualizer:
    """motif可视化器"""

    def __init__(self, output_dir: str, config: Dict[str, Any] = None):
        """
        初始化可视化器

        参数:
            output_dir: 输出目录
            config: 配置参数
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if config is None:
            config = {}

        # 默认配置
        defaults = {
            'dpi': 300,
            'format': 'png',
            'color_palette': 'Set2',
            'logo_style': 'classic',
            'max_motifs_per_plot': 10,
            'figsize': (10, 6),
        }

        self.config = {**defaults, **config}
        self.motifs = []
        self.comparisons = []

    def load_data(self, meme_dir: str, filtered_dir: str = None) -> bool:
        """
        加载motif数据

        参数:
            meme_dir: MEME结果目录
            filtered_dir: 过滤后结果目录（可选）

        返回:
            bool: 是否成功加载
        """
        # 优先使用过滤后的结果
        if filtered_dir and Path(filtered_dir).exists():
            logger.info(f"加载过滤后结果: {filtered_dir}")
            if self._load_filtered_results(filtered_dir):
                return True

        # 否则使用原始MEME结果
        logger.info(f"加载MEME结果: {meme_dir}")
        return self._load_meme_results(meme_dir)

    def _load_meme_results(self, meme_dir: str) -> bool:
        """加载MEME结果"""
        meme_path = Path(meme_dir)

        # 尝试加载JSON摘要
        summary_file = meme_path / 'meme_summary.json'
        if summary_file.exists():
            try:
                with open(summary_file, 'r') as f:
                    meme_summary = json.load(f)
                self.motifs = meme_summary.get('motifs', [])
                logger.info(f"加载 {len(self.motifs)} 个motif")
                return True
            except Exception as e:
                logger.warning(f"加载JSON摘要失败: {e}")

        # 尝试加载CSV
        csv_file = meme_path / 'motifs_summary.csv'
        if csv_file.exists():
            try:
                df = pd.read_csv(csv_file)
                self.motifs = df.to_dict('records')
                logger.info(f"从CSV加载 {len(self.motifs)} 个motif")
                return True
            except Exception as e:
                logger.warning(f"加载CSV失败: {e}")

        return False

    def _load_filtered_results(self, filtered_dir: str) -> bool:
        """加载过滤后结果"""
        filtered_path = Path(filtered_dir)

        # 尝试加载质量评估文件
        quality_file = filtered_path / 'motifs_with_quality.csv'
        if quality_file.exists():
            try:
                df = pd.read_csv(quality_file)
                self.motifs = df.to_dict('records')
                logger.info(f"从质量文件加载 {len(self.motifs)} 个motif")
                return True
            except Exception as e:
                logger.warning(f"加载质量文件失败: {e}")

        # 尝试加载基本motif文件
        motifs_file = filtered_path / 'filtered_motifs.csv'
        if motifs_file.exists():
            try:
                df = pd.read_csv(motifs_file)
                self.motifs = df.to_dict('records')
                logger.info(f"从过滤文件加载 {len(self.motifs)} 个motif")
                return True
            except Exception as e:
                logger.warning(f"加载过滤文件失败: {e}")

        return False

    def generate_all_visualizations(self) -> Dict[str, str]:
        """
        生成所有可视化图表

        返回:
            Dict: 生成的图片文件路径
        """
        if not HAS_MATPLOTLIB:
            logger.error("matplotlib未安装，无法生成可视化图表")
            return {}

        output_files = {}

        # 1. motif概览图
        if self.motifs:
            overview_file = self._create_motif_overview()
            if overview_file:
                output_files['overview'] = str(overview_file)

            # 2. 序列logo图（如果有频率矩阵）
            logo_files = self._create_sequence_logos()
            output_files.update(logo_files)

            # 3. motif质量分布图
            quality_file = self._create_quality_distribution()
            if quality_file:
                output_files['quality'] = str(quality_file)

            # 4. motif参数关系图
            param_file = self._create_parameter_relationships()
            if param_file:
                output_files['parameters'] = str(param_file)

        return output_files

    def _create_motif_overview(self) -> Optional[Path]:
        """创建motif概览图"""
        if not self.motifs:
            return None

        try:
            # 准备数据
            df = pd.DataFrame(self.motifs)
            max_motifs = min(len(df), self.config['max_motifs_per_plot'])
            df_top = df.head(max_motifs).copy()

            # 设置图形
            fig, axes = plt.subplots(2, 2, figsize=self.config['figsize'])
            fig.suptitle('Motif Analysis Overview', fontsize=16)

            # 1. E-value分布（子图1）
            if 'evalue' in df_top.columns:
                ax1 = axes[0, 0]
                df_top['log_evalue'] = -np.log10(df_top['evalue'].clip(lower=1e-10))
                df_top_sorted = df_top.sort_values('log_evalue', ascending=False)

                colors = plt.cm.Set2(np.arange(len(df_top_sorted)) / len(df_top_sorted))
                bars = ax1.bar(range(len(df_top_sorted)), df_top_sorted['log_evalue'], color=colors)

                ax1.set_xlabel('Motif Rank')
                ax1.set_ylabel('-log10(E-value)')
                ax1.set_title('Motif Significance')
                ax1.set_xticks(range(len(df_top_sorted)))
                ax1.set_xticklabels(df_top_sorted['id'], rotation=45, ha='right')

                # 添加数值标签
                for i, (bar, evalue) in enumerate(zip(bars, df_top_sorted['evalue'])):
                    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                            f'{evalue:.1e}', ha='center', va='bottom', fontsize=8)

            # 2. 位点数分布（子图2）
            if 'sites' in df_top.columns:
                ax2 = axes[0, 1]
                df_top_sorted = df_top.sort_values('sites', ascending=False)

                colors = plt.cm.Set2(np.arange(len(df_top_sorted)) / len(df_top_sorted))
                bars = ax2.bar(range(len(df_top_sorted)), df_top_sorted['sites'], color=colors)

                ax2.set_xlabel('Motif Rank')
                ax2.set_ylabel('Number of Sites')
                ax2.set_title('Motif Site Count')
                ax2.set_xticks(range(len(df_top_sorted)))
                ax2.set_xticklabels(df_top_sorted['id'], rotation=45, ha='right')

            # 3. 宽度分布（子图3）
            if 'width' in df_top.columns:
                ax3 = axes[1, 0]
                width_counts = df_top['width'].value_counts().sort_index()

                colors = plt.cm.Set2(np.arange(len(width_counts)) / len(width_counts))
                bars = ax3.bar(width_counts.index.astype(str), width_counts.values, color=colors)

                ax3.set_xlabel('Motif Width (nt)')
                ax3.set_ylabel('Count')
                ax3.set_title('Motif Width Distribution')

                # 添加数值标签
                for bar in bars:
                    height = bar.get_height()
                    ax3.text(bar.get_x() + bar.get_width()/2, height,
                            str(int(height)), ha='center', va='bottom')

            # 4. 质量评分分布（如果有）
            if 'quality_score' in df_top.columns:
                ax4 = axes[1, 1]
                df_top_sorted = df_top.sort_values('quality_score', ascending=False)

                colors = plt.cm.Set2(np.arange(len(df_top_sorted)) / len(df_top_sorted))
                bars = ax4.bar(range(len(df_top_sorted)), df_top_sorted['quality_score'], color=colors)

                ax4.set_xlabel('Motif Rank')
                ax4.set_ylabel('Quality Score')
                ax4.set_title('Motif Quality Assessment')
                ax4.set_xticks(range(len(df_top_sorted)))
                ax4.set_xticklabels(df_top_sorted['id'], rotation=45, ha='right')

                # 添加阈值线
                ax4.axhline(y=50, color='r', linestyle='--', alpha=0.5, label='Threshold')
                ax4.legend()

            plt.tight_layout()

            # 保存图片
            output_file = self.output_dir / f"motif_overview.{self.config['format']}"
            plt.savefig(output_file, dpi=self.config['dpi'], bbox_inches='tight')
            plt.close()

            logger.info(f"motif概览图已保存: {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"创建motif概览图时出错: {e}")
            return None

    def _create_sequence_logos(self) -> Dict[str, str]:
        """创建序列logo图"""
        logo_files = {}

        if not self.motifs:
            return logo_files

        try:
            # 检查是否有频率矩阵数据
            has_freq_matrix = any('frequency_matrix' in motif for motif in self.motifs)

            if has_freq_matrix:
                # 如果有频率矩阵，可以创建更精确的logo
                for i, motif in enumerate(self.motifs[:5]):  # 只创建前5个
                    if 'frequency_matrix' in motif:
                        logo_file = self._create_detailed_logo(motif, i)
                        if logo_file:
                            logo_files[f'logo_{i+1}'] = str(logo_file)
            else:
                # 否则创建简单的共识序列图
                for i, motif in enumerate(self.motifs[:5]):
                    logo_file = self._create_simple_logo(motif, i)
                    if logo_file:
                        logo_files[f'logo_{i+1}'] = str(logo_file)

        except Exception as e:
            logger.error(f"创建序列logo图时出错: {e}")

        return logo_files

    def _create_simple_logo(self, motif: Dict[str, Any], index: int) -> Optional[Path]:
        """创建简单的共识序列logo"""
        try:
            consensus = motif.get('consensus', '')
            if not consensus:
                # 尝试从ID或其他字段提取
                motif_id = motif.get('id', f'motif_{index+1}')
                consensus = motif_id.split('_')[-1] if '_' in motif_id else 'NNNNNN'

            width = motif.get('width', len(consensus))
            evalue = motif.get('evalue', 1.0)
            sites = motif.get('sites', 0)

            fig, ax = plt.subplots(figsize=(max(width/2, 6), 4))

            # 创建简单的序列表示
            positions = range(1, width + 1)
            for pos in positions:
                if pos - 1 < len(consensus):
                    base = consensus[pos - 1]
                    # 为不同碱基使用不同颜色
                    if base == 'A':
                        color = 'red'
                    elif base == 'C':
                        color = 'blue'
                    elif base == 'G':
                        color = 'orange'
                    elif base == 'T':
                        color = 'green'
                    else:
                        color = 'gray'

                    rect = Rectangle((pos - 0.4, 0.1), 0.8, 0.8,
                                   facecolor=color, edgecolor='black')
                    ax.add_patch(rect)

                    # 添加碱基标签
                    ax.text(pos, 0.5, base, ha='center', va='center',
                           fontsize=12, fontweight='bold', color='white')

            ax.set_xlim(0, width + 1)
            ax.set_ylim(0, 1)
            ax.set_xlabel('Position')
            ax.set_title(f'Motif {motif.get("id", f"motif_{index+1}")}\n'
                        f'Width: {width}nt, Sites: {sites}, E-value: {evalue:.1e}',
                        fontsize=10)
            ax.set_yticks([])
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

            output_file = self.output_dir / f"logo_simple_{index+1}.{self.config['format']}"
            plt.savefig(output_file, dpi=self.config['dpi'], bbox_inches='tight')
            plt.close()

            logger.info(f"简单logo图已保存: {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"创建简单logo图时出错: {e}")
            return None

    def _create_detailed_logo(self, motif: Dict[str, Any], index: int) -> Optional[Path]:
        """创建详细的序列logo（需要频率矩阵）"""
        try:
            freq_matrix = motif.get('frequency_matrix', [])
            if not freq_matrix:
                return None

            width = len(freq_matrix)
            evalue = motif.get('evalue', 1.0)
            sites = motif.get('sites', 0)

            fig, ax = plt.subplots(figsize=(max(width/2, 8), 6))

            # 计算信息含量
            information_content = []
            for pos in freq_matrix:
                total = sum(pos.values())
                if total > 0:
                    # 计算相对频率
                    freqs = {base: count/total for base, count in pos.items()}
                    # 计算信息含量（简化版本）
                    ic = 2 + sum(f * np.log2(f) if f > 0 else 0 for f in freqs.values())
                    information_content.append(ic)
                else:
                    information_content.append(0)

            # 创建堆叠条形图
            bottom = np.zeros(width)
            colors = {'A': 'red', 'C': 'blue', 'G': 'orange', 'T': 'green'}

            for base in ['A', 'C', 'G', 'T']:
                heights = []
                for pos in freq_matrix:
                    total = sum(pos.values())
                    if total > 0:
                        height = pos.get(base, 0) / total * information_content[len(heights)]
                    else:
                        height = 0
                    heights.append(height)

                ax.bar(range(1, width + 1), heights, bottom=bottom,
                      color=colors[base], label=base, edgecolor='black')
                bottom += heights

            ax.set_xlim(0.5, width + 0.5)
            ax.set_ylim(0, max(information_content) * 1.1)
            ax.set_xlabel('Position')
            ax.set_ylabel('Information Content (bits)')
            ax.set_title(f'Sequence Logo: {motif.get("id", f"motif_{index+1}")}\n'
                        f'Width: {width}nt, Sites: {sites}, E-value: {evalue:.1e}')
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3)

            plt.tight_layout()

            output_file = self.output_dir / f"logo_detailed_{index+1}.{self.config['format']}"
            plt.savefig(output_file, dpi=self.config['dpi'], bbox_inches='tight')
            plt.close()

            logger.info(f"详细logo图已保存: {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"创建详细logo图时出错: {e}")
            return None

    def _create_quality_distribution(self) -> Optional[Path]:
        """创建motif质量分布图"""
        if not self.motifs:
            return None

        try:
            df = pd.DataFrame(self.motifs)

            # 检查是否有质量评分
            if 'quality_score' not in df.columns:
                # 如果没有质量评分，创建基于其他参数的质量评估
                quality_scores = []
                for motif in self.motifs:
                    score = 0
                    evalue = motif.get('evalue', 1.0)
                    sites = motif.get('sites', 0)
                    llr = motif.get('llr', 0)

                    # 简单评分规则
                    if evalue <= 1e-5:
                        score += 30
                    elif evalue <= 1e-3:
                        score += 20
                    elif evalue <= 0.05:
                        score += 10

                    if sites >= 10:
                        score += 30
                    elif sites >= 5:
                        score += 20

                    if llr >= 20:
                        score += 40
                    elif llr >= 10:
                        score += 30
                    elif llr >= 5:
                        score += 20

                    quality_scores.append(score)

                df['quality_score'] = quality_scores

            fig, axes = plt.subplots(1, 3, figsize=(15, 5))

            # 1. 质量评分分布
            ax1 = axes[0]
            ax1.hist(df['quality_score'], bins=20, color='skyblue', edgecolor='black')
            ax1.set_xlabel('Quality Score')
            ax1.set_ylabel('Count')
            ax1.set_title('Motif Quality Score Distribution')
            ax1.axvline(x=50, color='red', linestyle='--', label='Threshold')
            ax1.legend()

            # 2. 质量 vs E-value
            ax2 = axes[1]
            if 'evalue' in df.columns:
                scatter = ax2.scatter(-np.log10(df['evalue'].clip(lower=1e-10)),
                                    df['quality_score'],
                                    c=df['sites'] if 'sites' in df.columns else 'blue',
                                    cmap='viridis', alpha=0.6)
                ax2.set_xlabel('-log10(E-value)')
                ax2.set_ylabel('Quality Score')
                ax2.set_title('Quality vs Significance')
                if 'sites' in df.columns:
                    plt.colorbar(scatter, ax=ax2, label='Site Count')

            # 3. 质量分类
            ax3 = axes[2]
            quality_categories = []
            for score in df['quality_score']:
                if score >= 70:
                    quality_categories.append('Excellent')
                elif score >= 50:
                    quality_categories.append('Good')
                elif score >= 30:
                    quality_categories.append('Fair')
                else:
                    quality_categories.append('Poor')

            category_counts = pd.Series(quality_categories).value_counts()
            colors = ['green', 'lightgreen', 'orange', 'red']
            ax3.pie(category_counts.values, labels=category_counts.index,
                   autopct='%1.1f%%', colors=colors[:len(category_counts)])
            ax3.set_title('Motif Quality Categories')

            plt.suptitle('Motif Quality Assessment', fontsize=16)
            plt.tight_layout()

            output_file = self.output_dir / f"quality_distribution.{self.config['format']}"
            plt.savefig(output_file, dpi=self.config['dpi'], bbox_inches='tight')
            plt.close()

            logger.info(f"质量分布图已保存: {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"创建质量分布图时出错: {e}")
            return None

    def _create_parameter_relationships(self) -> Optional[Path]:
        """创建motif参数关系图"""
        if not self.motifs or len(self.motifs) < 3:
            return None

        try:
            df = pd.DataFrame(self.motifs)
            required_columns = ['evalue', 'sites', 'width']
            if not all(col in df.columns for col in required_columns):
                return None

            # 设置Seaborn样式（如果可用）
            if HAS_SEABORN:
                sns.set_style("whitegrid")
                palette = sns.color_palette(self.config['color_palette'])
            else:
                palette = plt.cm.Set2(np.arange(len(df)) / len(df))

            fig, axes = plt.subplots(2, 2, figsize=(12, 10))

            # 1. E-value vs 位点数
            ax1 = axes[0, 0]
            scatter1 = ax1.scatter(-np.log10(df['evalue'].clip(lower=1e-10)),
                                 df['sites'],
                                 c=df['width'] if 'width' in df.columns else 'blue',
                                 cmap='viridis', alpha=0.6, s=50)
            ax1.set_xlabel('-log10(E-value)')
            ax1.set_ylabel('Number of Sites')
            ax1.set_title('Significance vs Site Count')
            if 'width' in df.columns:
                plt.colorbar(scatter1, ax=ax1, label='Motif Width')

            # 2. 宽度 vs 位点数
            ax2 = axes[0, 1]
            scatter2 = ax2.scatter(df['width'], df['sites'],
                                 c=-np.log10(df['evalue'].clip(lower=1e-10)),
                                 cmap='plasma', alpha=0.6, s=50)
            ax2.set_xlabel('Motif Width (nt)')
            ax2.set_ylabel('Number of Sites')
            ax2.set_title('Width vs Site Count')
            plt.colorbar(scatter2, ax=ax2, label='-log10(E-value)')

            # 3. 相关矩阵热图
            ax3 = axes[1, 0]
            numeric_cols = ['evalue', 'sites', 'width']
            if 'llr' in df.columns:
                numeric_cols.append('llr')
            if 'quality_score' in df.columns:
                numeric_cols.append('quality_score')

            numeric_df = df[numeric_cols].copy()
            numeric_df['log_evalue'] = -np.log10(numeric_df['evalue'].clip(lower=1e-10))
            numeric_df = numeric_df.drop('evalue', axis=1)

            # 计算相关性
            corr_matrix = numeric_df.corr()

            im = ax3.imshow(corr_matrix, cmap='coolwarm', vmin=-1, vmax=1)
            ax3.set_xticks(range(len(corr_matrix.columns)))
            ax3.set_yticks(range(len(corr_matrix.columns)))
            ax3.set_xticklabels(corr_matrix.columns, rotation=45, ha='right')
            ax3.set_yticklabels(corr_matrix.columns)
            ax3.set_title('Parameter Correlation Matrix')

            # 添加相关性数值
            for i in range(len(corr_matrix.columns)):
                for j in range(len(corr_matrix.columns)):
                    text = ax3.text(j, i, f'{corr_matrix.iloc[i, j]:.2f}',
                                   ha="center", va="center", color="black")

            plt.colorbar(im, ax=ax3)

            # 4. 3D散点图（如果可用）
            ax4 = axes[1, 1]
            if 'llr' in df.columns:
                scatter3 = ax4.scatter(df['width'], df['sites'],
                                      s=df['llr']/5,  # 大小表示LLR
                                      c=-np.log10(df['evalue'].clip(lower=1e-10)),
                                      cmap='rainbow', alpha=0.6)
                ax4.set_xlabel('Width')
                ax4.set_ylabel('Sites')
                ax4.set_title('Width vs Sites (size=LLR, color=E-value)')
                plt.colorbar(scatter3, ax=ax4, label='-log10(E-value)')
            else:
                ax4.axis('off')
                ax4.text(0.5, 0.5, 'LLR data not available\nfor 3D visualization',
                        ha='center', va='center', fontsize=12)

            plt.suptitle('Motif Parameter Relationships', fontsize=16)
            plt.tight_layout()

            output_file = self.output_dir / f"parameter_relationships.{self.config['format']}"
            plt.savefig(output_file, dpi=self.config['dpi'], bbox_inches='tight')
            plt.close()

            logger.info(f"参数关系图已保存: {output_file}")
            return output_file

        except Exception as e:
            logger.error(f"创建参数关系图时出错: {e}")
            return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="motif可视化脚本 - 生成序列logo和分布图",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法：使用MEME结果
  python visualize_motifs.py --meme meme_results/ --output visualization/

  # 使用过滤后结果
  python visualize_motifs.py --meme meme_results/ --filtered filtered_results/ --output viz/

  # 自定义格式和DPI
  python visualize_motifs.py --meme meme_results/ --output results/ --format pdf --dpi 600
        """
    )

    # 输入参数
    parser.add_argument('--meme', required=True, help='MEME结果目录')
    parser.add_argument('--filtered', help='过滤后结果目录（可选，优先使用）')
    parser.add_argument('--output', required=True, help='输出目录')

    # 可视化参数
    parser.add_argument('--format', choices=['png', 'pdf', 'svg', 'jpg'], default='png',
                       help='输出图片格式（默认: png）')
    parser.add_argument('--dpi', type=int, default=300,
                       help='图片DPI（默认: 300）')
    parser.add_argument('--max-motifs', type=int, default=10,
                       help='每张图最多显示motif数（默认: 10）')

    args = parser.parse_args()

    if not HAS_MATPLOTLIB:
        logger.error("错误: matplotlib库未安装")
        logger.error("请使用以下命令安装: pip install matplotlib seaborn")
        return 1

    logger.info("=== motif可视化开始 ===")
    logger.info(f"MEME结果目录: {args.meme}")
    if args.filtered:
        logger.info(f"过滤后结果目录: {args.filtered}")
    logger.info(f"输出目录: {args.output}")
    logger.info(f"图片格式: {args.format}, DPI: {args.dpi}")

    # 创建可视化器
    config = {
        'format': args.format,
        'dpi': args.dpi,
        'max_motifs_per_plot': args.max_motifs,
    }

    visualizer = MotifVisualizer(args.output, config)

    # 加载数据
    if not visualizer.load_data(args.meme, args.filtered):
        logger.error("加载motif数据失败")
        return 1

    # 生成可视化图表
    output_files = visualizer.generate_all_visualizations()

    logger.info("=== 可视化完成 ===")
    logger.info(f"输出目录: {args.output}")

    if output_files:
        logger.info("生成的图片文件:")
        for file_type, file_path in output_files.items():
            logger.info(f"  {file_type}: {Path(file_path).name}")
    else:
        logger.warning("未生成任何可视化图表")

    return 0


if __name__ == '__main__':
    sys.exit(main())