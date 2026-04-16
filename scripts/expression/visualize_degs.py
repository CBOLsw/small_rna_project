#!/usr/bin/env python3
"""
差异表达可视化脚本

功能：
1. 火山图（Volcano plot）
2. MA图（MA plot）
3. 热图（Heatmap）
4. 表达分布图
5. 维恩图（多组比较）
6. 通路富集可视化（如果可用）

使用方法：
    python visualize_degs.py --input <差异表达结果文件> --counts <计数矩阵文件> --output <输出目录>
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

# 尝试导入可视化库
try:
    import matplotlib
    matplotlib.use('Agg')  # 非交互模式
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("警告: matplotlib/seaborn未安装，将跳过图表生成")

try:
    from scipy import stats, cluster
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class DEGVisualizer:
    """差异表达可视化器"""

    def __init__(self, style: str = 'default'):
        """
        初始化可视化器

        参数:
            style: 绘图风格 ('default', 'seaborn', 'ggplot')
        """
        self.style = style
        self.results = None
        self.count_matrix = None
        self.metadata = None

        # 设置绘图风格
        if HAS_MATPLOTLIB:
            self._set_plot_style(style)

    def _set_plot_style(self, style: str):
        """设置绘图风格"""
        if style == 'seaborn':
            plt.style.use('seaborn-v0_8-whitegrid')
            sns.set_palette("husl")
        elif style == 'ggplot':
            plt.style.use('ggplot')
        else:
            plt.style.use('default')

        # 设置中文字体（如果需要）
        try:
            plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial']
            plt.rcParams['axes.unicode_minus'] = False
        except:
            pass

    def load_data(self, deg_file: str, count_file: Optional[str] = None,
                 metadata_file: Optional[str] = None):
        """
        加载数据

        参数:
            deg_file: 差异表达结果文件
            count_file: 计数矩阵文件（可选）
            metadata_file: 样本信息文件（可选）
        """
        # 1. 加载差异表达结果
        logger.info(f"加载差异表达结果: {deg_file}")
        self.results = self._load_dataframe(deg_file)

        if self.results is None or self.results.empty:
            logger.error("无法加载差异表达结果")
            return False

        logger.info(f"差异表达结果维度: {self.results.shape}")
        logger.info(f"列名: {list(self.results.columns)}")

        # 2. 加载计数矩阵（如果提供）
        if count_file:
            logger.info(f"加载计数矩阵: {count_file}")
            self.count_matrix = self._load_dataframe(count_file, index_col=0)

            if self.count_matrix is not None:
                logger.info(f"计数矩阵维度: {self.count_matrix.shape}")
            else:
                logger.warning("无法加载计数矩阵")

        # 3. 加载样本信息（如果提供）
        if metadata_file:
            logger.info(f"加载样本信息: {metadata_file}")
            self.metadata = self._load_dataframe(metadata_file)

            if self.metadata is not None:
                logger.info(f"样本信息维度: {self.metadata.shape}")
            else:
                logger.warning("无法加载样本信息")

        return True

    def _load_dataframe(self, file_path: str, index_col: Optional[int] = None) -> Optional[pd.DataFrame]:
        """加载数据框"""
        try:
            path = Path(file_path)
            if not path.exists():
                logger.error(f"文件不存在: {file_path}")
                return None

            file_ext = path.suffix.lower()

            if file_ext == '.csv':
                df = pd.read_csv(file_path, index_col=index_col)
            elif file_ext in ['.tsv', '.txt', '.tab']:
                df = pd.read_csv(file_path, sep='\t', index_col=index_col)
            elif file_ext in ['.xlsx', '.xls']:
                df = pd.read_excel(file_path, index_col=index_col)
            else:
                logger.error(f"不支持的文件格式: {file_ext}")
                return None

            return df

        except Exception as e:
            logger.error(f"加载文件失败 {file_path}: {e}")
            return None

    def create_volcano_plot(self, output_dir: Path,
                           log2fc_threshold: float = 1.0,
                           padj_threshold: float = 0.05,
                           title: str = "火山图",
                           figsize: Tuple[int, int] = (10, 8)):
        """
        创建火山图

        参数:
            output_dir: 输出目录
            log2fc_threshold: log2倍数变化阈值
            padj_threshold: 调整后p值阈值
            title: 图表标题
            figsize: 图形大小
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib未安装，跳过火山图生成")
            return None

        if self.results is None or self.results.empty:
            logger.error("无差异表达结果数据")
            return None

        # 检查必要的列
        required_cols = ['log2FoldChange']
        if 'padj' in self.results.columns:
            pval_col = 'padj'
        elif 'pvalue' in self.results.columns:
            pval_col = 'pvalue'
        else:
            logger.error("结果中缺少p值列 (padj或pvalue)")
            return None

        required_cols.append(pval_col)

        missing_cols = [col for col in required_cols if col not in self.results.columns]
        if missing_cols:
            logger.error(f"结果中缺少必要的列: {missing_cols}")
            return None

        logger.info(f"创建火山图 (log2FC≥{log2fc_threshold}, {pval_col}<{padj_threshold})")

        try:
            # 准备数据
            df = self.results.copy()
            df = df.dropna(subset=[pval_col, 'log2FoldChange'])

            # 计算-log10(p值)
            df['neg_log10_pval'] = -np.log10(df[pval_col])

            # 标记差异表达基因
            conditions = [
                (df['log2FoldChange'].abs() >= log2fc_threshold) & (df[pval_col] < padj_threshold) & (df['log2FoldChange'] > 0),
                (df['log2FoldChange'].abs() >= log2fc_threshold) & (df[pval_col] < padj_threshold) & (df['log2FoldChange'] < 0),
                True  # 默认值
            ]
            choices = ['上调', '下调', '非差异']
            df['regulation'] = np.select(conditions, choices)

            # 创建图形
            plt.figure(figsize=figsize)

            # 定义颜色
            colors = {'上调': '#E74C3C', '下调': '#3498DB', '非差异': '#95A5A6'}
            color_vector = df['regulation'].map(colors)

            # 散点图
            scatter = plt.scatter(df['log2FoldChange'], df['neg_log10_pval'],
                                 c=color_vector, alpha=0.6, s=20, edgecolors='none')

            # 添加阈值线
            plt.axhline(y=-np.log10(padj_threshold), color='gray', linestyle='--', alpha=0.5, linewidth=1)
            plt.axvline(x=log2fc_threshold, color='gray', linestyle='--', alpha=0.5, linewidth=1)
            plt.axvline(x=-log2fc_threshold, color='gray', linestyle='--', alpha=0.5, linewidth=1)

            # 添加标签和标题
            plt.xlabel('log₂ 倍数变化', fontsize=12)
            plt.ylabel(f'-log₁₀ {pval_col}', fontsize=12)
            plt.title(title, fontsize=14, fontweight='bold')

            # 添加图例
            legend_elements = [
                mpatches.Patch(color=colors['上调'], label=f'上调 (n={(df["regulation"] == "上调").sum()})'),
                mpatches.Patch(color=colors['下调'], label=f'下调 (n={(df["regulation"] == "下调").sum()})'),
                mpatches.Patch(color=colors['非差异'], label=f'非差异 (n={(df["regulation"] == "非差异").sum()})')
            ]
            plt.legend(handles=legend_elements, loc='best', frameon=True)

            # 添加网格
            plt.grid(True, alpha=0.3, linestyle='--')

            # 调整布局
            plt.tight_layout()

            # 保存图形
            output_file = output_dir / "volcano_plot.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"火山图已保存: {output_file}")

            # 保存绘图数据
            data_file = output_dir / "volcano_plot_data.csv"
            df.to_csv(data_file)
            logger.info(f"火山图数据已保存: {data_file}")

            return str(output_file)

        except Exception as e:
            logger.error(f"创建火山图时出错: {e}")
            return None

    def create_ma_plot(self, output_dir: Path,
                      padj_threshold: float = 0.05,
                      log2fc_threshold: float = 1.0,
                      title: str = "MA图",
                      figsize: Tuple[int, int] = (10, 8)):
        """
        创建MA图

        参数:
            output_dir: 输出目录
            padj_threshold: 调整后p值阈值
            log2fc_threshold: log2倍数变化阈值
            title: 图表标题
            figsize: 图形大小
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib未安装，跳过MA图生成")
            return None

        if self.results is None or self.results.empty:
            logger.error("无差异表达结果数据")
            return None

        # 检查必要的列
        required_cols = ['log2FoldChange']
        if 'baseMean' in self.results.columns:
            mean_col = 'baseMean'
        elif 'baseMean_control' in self.results.columns and 'baseMean_treatment' in self.results.columns:
            # 计算平均表达
            self.results['baseMean'] = (self.results['baseMean_control'] + self.results['baseMean_treatment']) / 2
            mean_col = 'baseMean'
        else:
            logger.error("结果中缺少基础表达水平列 (baseMean)")
            return None

        if 'padj' in self.results.columns:
            pval_col = 'padj'
        elif 'pvalue' in self.results.columns:
            pval_col = 'pvalue'
        else:
            logger.error("结果中缺少p值列 (padj或pvalue)")
            return None

        logger.info(f"创建MA图 ({pval_col}<{padj_threshold}, |log2FC|≥{log2fc_threshold})")

        try:
            # 准备数据
            df = self.results.copy()
            df = df.dropna(subset=[pval_col, 'log2FoldChange', mean_col])

            # 标记差异表达基因
            conditions = [
                (df['log2FoldChange'].abs() >= log2fc_threshold) & (df[pval_col] < padj_threshold) & (df['log2FoldChange'] > 0),
                (df['log2FoldChange'].abs() >= log2fc_threshold) & (df[pval_col] < padj_threshold) & (df['log2FoldChange'] < 0),
                True  # 默认值
            ]
            choices = ['上调', '下调', '非差异']
            df['regulation'] = np.select(conditions, choices)

            # 创建图形
            plt.figure(figsize=figsize)

            # 定义颜色
            colors = {'上调': '#E74C3C', '下调': '#3498DB', '非差异': '#95A5A6'}
            color_vector = df['regulation'].map(colors)

            # 散点图
            plt.scatter(np.log10(df[mean_col]), df['log2FoldChange'],
                       c=color_vector, alpha=0.6, s=20, edgecolors='none')

            # 添加阈值线
            plt.axhline(y=0, color='black', linestyle='-', alpha=0.3, linewidth=1)
            plt.axhline(y=log2fc_threshold, color='gray', linestyle='--', alpha=0.5, linewidth=1)
            plt.axhline(y=-log2fc_threshold, color='gray', linestyle='--', alpha=0.5, linewidth=1)

            # 添加标签和标题
            plt.xlabel('log₁₀ 平均表达', fontsize=12)
            plt.ylabel('log₂ 倍数变化', fontsize=12)
            plt.title(title, fontsize=14, fontweight='bold')

            # 添加图例
            legend_elements = [
                mpatches.Patch(color=colors['上调'], label=f'上调 (n={(df["regulation"] == "上调").sum()})'),
                mpatches.Patch(color=colors['下调'], label=f'下调 (n={(df["regulation"] == "下调").sum()})'),
                mpatches.Patch(color=colors['非差异'], label=f'非差异 (n={(df["regulation"] == "非差异").sum()})')
            ]
            plt.legend(handles=legend_elements, loc='best', frameon=True)

            # 添加网格
            plt.grid(True, alpha=0.3, linestyle='--')

            # 调整布局
            plt.tight_layout()

            # 保存图形
            output_file = output_dir / "ma_plot.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"MA图已保存: {output_file}")

            # 保存绘图数据
            data_file = output_dir / "ma_plot_data.csv"
            df.to_csv(data_file)
            logger.info(f"MA图数据已保存: {data_file}")

            return str(output_file)

        except Exception as e:
            logger.error(f"创建MA图时出错: {e}")
            return None

    def create_heatmap(self, output_dir: Path,
                      top_n: int = 50,
                      group_col: Optional[str] = None,
                      title: str = "差异表达基因热图",
                      figsize: Tuple[int, int] = (12, 10)):
        """
        创建差异表达基因热图

        参数:
            output_dir: 输出目录
            top_n: 显示top N差异表达基因
            group_col: 分组列名（用于注释）
            title: 图表标题
            figsize: 图形大小
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib未安装，跳过热图生成")
            return None

        if self.results is None or self.results.empty:
            logger.error("无差异表达结果数据")
            return None

        if self.count_matrix is None or self.count_matrix.empty:
            logger.error("无计数矩阵数据")
            return None

        logger.info(f"创建热图 (top {top_n} 差异表达基因)")

        try:
            # 准备差异表达基因数据
            df = self.results.copy()

            # 排序：按显著性和倍数变化
            sort_cols = []
            if 'padj' in df.columns:
                sort_cols.append('padj')
            elif 'pvalue' in df.columns:
                sort_cols.append('pvalue')

            if 'log2FoldChange' in df.columns:
                sort_cols.append('log2FoldChange')

            if sort_cols:
                df = df.sort_values(by=sort_cols, ascending=[True, False])
            else:
                df = df.sort_index()

            # 获取top基因
            top_genes = df.head(top_n).index.tolist()

            # 从计数矩阵中提取这些基因的表达数据
            common_genes = [gene for gene in top_genes if gene in self.count_matrix.index]
            if not common_genes:
                logger.error("计数矩阵中未找到差异表达基因")
                return None

            expr_data = self.count_matrix.loc[common_genes]

            # 标准化表达数据（Z-score标准化）
            expr_data_z = expr_data.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0, axis=1)

            # 创建图形
            plt.figure(figsize=figsize)

            # 创建热图
            if HAS_SCIPY:
                # 使用层次聚类
                import seaborn as sns

                # 创建列注释（如果提供了分组信息）
                col_colors = None
                if group_col and self.metadata is not None and group_col in self.metadata.columns:
                    # 创建分组颜色映射
                    groups = self.metadata[group_col].astype('category')
                    group_palette = sns.color_palette("husl", len(groups.cat.categories))
                    group_lut = dict(zip(groups.cat.categories, group_palette))

                    # 创建列颜色
                    col_colors = groups.map(group_lut)

                # 绘制热图
                g = sns.clustermap(expr_data_z,
                                  cmap='RdBu_r',
                                  center=0,
                                  figsize=figsize,
                                  row_cluster=True,
                                  col_cluster=True,
                                  yticklabels=True,
                                  xticklabels=True,
                                  cbar_kws={"label": "Z-score"},
                                  col_colors=col_colors)

                plt.title(title, fontsize=14, fontweight='bold', pad=20)

                # 调整布局
                g.fig.subplots_adjust(right=0.85)

                # 保存图形
                output_file = output_dir / "deg_heatmap.png"
                g.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()

            else:
                # 简单热图（无聚类）
                import seaborn as sns

                # 创建热图
                plt.figure(figsize=figsize)
                sns.heatmap(expr_data_z,
                           cmap='RdBu_r',
                           center=0,
                           yticklabels=True,
                           xticklabels=True,
                           cbar_kws={"label": "Z-score"})

                plt.title(title, fontsize=14, fontweight='bold')
                plt.xlabel('样本', fontsize=12)
                plt.ylabel('基因', fontsize=12)

                # 调整布局
                plt.tight_layout()

                # 保存图形
                output_file = output_dir / "deg_heatmap.png"
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
                plt.close()

            logger.info(f"热图已保存: {output_file}")

            # 保存热图数据
            data_file = output_dir / "heatmap_data.csv"
            expr_data_z.to_csv(data_file)
            logger.info(f"热图数据已保存: {data_file}")

            return str(output_file)

        except Exception as e:
            logger.error(f"创建热图时出错: {e}")
            return None

    def create_expression_distribution(self, output_dir: Path,
                                      group_col: Optional[str] = None,
                                      title: str = "表达分布图",
                                      figsize: Tuple[int, int] = (12, 8)):
        """
        创建表达分布图

        参数:
            output_dir: 输出目录
            group_col: 分组列名
            title: 图表标题
            figsize: 图形大小
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib未安装，跳过表达分布图生成")
            return None

        if self.results is None or self.results.empty:
            logger.error("无差异表达结果数据")
            return None

        logger.info("创建表达分布图")

        try:
            # 准备数据
            df = self.results.copy()

            # 检查是否有分组信息
            if group_col and self.metadata is not None and group_col in self.metadata.columns:
                # 如果有分组，创建分组表达分布
                groups = self.metadata[group_col].unique()
                n_groups = len(groups)

                if n_groups == 2:
                    # 两个组：创建对比图
                    fig, axes = plt.subplots(1, 3, figsize=figsize)

                    # 获取组名
                    group1, group2 = groups

                    # 检查是否有分组表达数据
                    if f'baseMean_{group1}' in df.columns and f'baseMean_{group2}' in df.columns:
                        # 图1: 组1表达分布
                        axes[0].hist(np.log10(df[f'baseMean_{group1}'].dropna() + 1),
                                    bins=50, alpha=0.7, color='blue', edgecolor='black')
                        axes[0].set_title(f'{group1} 表达分布', fontsize=12)
                        axes[0].set_xlabel('log₁₀(表达+1)', fontsize=10)
                        axes[0].set_ylabel('基因数', fontsize=10)

                        # 图2: 组2表达分布
                        axes[1].hist(np.log10(df[f'baseMean_{group2}'].dropna() + 1),
                                    bins=50, alpha=0.7, color='red', edgecolor='black')
                        axes[1].set_title(f'{group2} 表达分布', fontsize=12)
                        axes[1].set_xlabel('log₁₀(表达+1)', fontsize=10)
                        axes[1].set_ylabel('基因数', fontsize=10)

                        # 图3: 两组对比
                        axes[2].hist(np.log10(df[f'baseMean_{group1}'].dropna() + 1),
                                    bins=50, alpha=0.5, color='blue', edgecolor='black', label=group1)
                        axes[2].hist(np.log10(df[f'baseMean_{group2}'].dropna() + 1),
                                    bins=50, alpha=0.5, color='red', edgecolor='black', label=group2)
                        axes[2].set_title('表达分布对比', fontsize=12)
                        axes[2].set_xlabel('log₁₀(表达+1)', fontsize=10)
                        axes[2].set_ylabel('基因数', fontsize=10)
                        axes[2].legend()

                    else:
                        # 没有分组表达数据，使用平均表达
                        if 'baseMean' in df.columns:
                            axes[0].hist(np.log10(df['baseMean'].dropna() + 1),
                                        bins=50, alpha=0.7, color='blue', edgecolor='black')
                            axes[0].set_title('平均表达分布', fontsize=12)
                            axes[0].set_xlabel('log₁₀(表达+1)', fontsize=10)
                            axes[0].set_ylabel('基因数', fontsize=10)

                            # 隐藏其他子图
                            axes[1].axis('off')
                            axes[2].axis('off')

                else:
                    # 多个组或多个子图
                    pass

            else:
                # 无分组信息，创建简单表达分布
                fig, axes = plt.subplots(1, 2, figsize=figsize)

                # 图1: 平均表达分布
                if 'baseMean' in df.columns:
                    axes[0].hist(np.log10(df['baseMean'].dropna() + 1),
                                bins=50, alpha=0.7, color='steelblue', edgecolor='black')
                    axes[0].set_title('平均表达分布', fontsize=12)
                    axes[0].set_xlabel('log₁₀(表达+1)', fontsize=10)
                    axes[0].set_ylabel('基因数', fontsize=10)

                # 图2: log2FC分布
                if 'log2FoldChange' in df.columns:
                    axes[1].hist(df['log2FoldChange'].dropna(),
                                bins=50, alpha=0.7, color='forestgreen', edgecolor='black')
                    axes[1].set_title('log₂倍数变化分布', fontsize=12)
                    axes[1].set_xlabel('log₂倍数变化', fontsize=10)
                    axes[1].set_ylabel('基因数', fontsize=10)

            # 设置总标题
            fig.suptitle(title, fontsize=14, fontweight='bold')

            # 调整布局
            plt.tight_layout()

            # 保存图形
            output_file = output_dir / "expression_distribution.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()

            logger.info(f"表达分布图已保存: {output_file}")

            return str(output_file)

        except Exception as e:
            logger.error(f"创建表达分布图时出错: {e}")
            return None

    def generate_all_plots(self, output_dir: str,
                          log2fc_threshold: float = 1.0,
                          padj_threshold: float = 0.05,
                          top_n: int = 50,
                          group_col: Optional[str] = None):
        """
        生成所有可视化图表

        参数:
            output_dir: 输出目录
            log2fc_threshold: log2倍数变化阈值
            padj_threshold: 调整后p值阈值
            top_n: 热图显示的top基因数
            group_col: 分组列名
        """
        if not HAS_MATPLOTLIB:
            logger.warning("matplotlib未安装，跳过所有图表生成")
            return []

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        generated_files = []

        # 1. 火山图
        volcano_file = self.create_volcano_plot(
            output_path,
            log2fc_threshold=log2fc_threshold,
            padj_threshold=padj_threshold,
            title=f"火山图 (|log2FC|≥{log2fc_threshold}, padj<{padj_threshold})"
        )
        if volcano_file:
            generated_files.append(volcano_file)

        # 2. MA图
        ma_file = self.create_ma_plot(
            output_path,
            padj_threshold=padj_threshold,
            log2fc_threshold=log2fc_threshold,
            title=f"MA图 (padj<{padj_threshold}, |log2FC|≥{log2fc_threshold})"
        )
        if ma_file:
            generated_files.append(ma_file)

        # 3. 热图（需要计数矩阵）
        if self.count_matrix is not None and not self.count_matrix.empty:
            heatmap_file = self.create_heatmap(
                output_path,
                top_n=top_n,
                group_col=group_col,
                title=f"Top {top_n} 差异表达基因热图"
            )
            if heatmap_file:
                generated_files.append(heatmap_file)

        # 4. 表达分布图
        dist_file = self.create_expression_distribution(
            output_path,
            group_col=group_col,
            title="基因表达分布"
        )
        if dist_file:
            generated_files.append(dist_file)

        # 生成文件清单
        if generated_files:
            manifest_file = output_path / "visualization_manifest.txt"
            with open(manifest_file, 'w') as f:
                f.write("=== 差异表达可视化文件清单 ===\n\n")
                f.write(f"生成时间: {pd.Timestamp.now()}\n\n")
                f.write("生成的文件:\n")
                for file_path in generated_files:
                    f.write(f"  {Path(file_path).name}\n")

            generated_files.append(str(manifest_file))

        logger.info(f"生成 {len(generated_files)} 个可视化文件")
        return generated_files


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="差异表达可视化脚本"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="差异表达结果文件（CSV/TSV格式）"
    )

    parser.add_argument(
        "--counts", "-c",
        help="计数矩阵文件（用于热图）"
    )

    parser.add_argument(
        "--metadata", "-m",
        help="样本信息文件（用于分组）"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/expression/visualization",
        help="输出目录 (默认: results/expression/visualization)"
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
        "--top-n",
        type=int,
        default=50,
        help="热图显示的top基因数 (默认: 50)"
    )

    parser.add_argument(
        "--group-col",
        help="分组列名（用于热图注释和表达分布）"
    )

    parser.add_argument(
        "--style",
        choices=['default', 'seaborn', 'ggplot'],
        default='default',
        help="绘图风格 (默认: default)"
    )

    parser.add_argument(
        "--all-plots",
        action="store_true",
        default=True,
        help="生成所有图表 (默认: True)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 检查matplotlib是否可用
    if not HAS_MATPLOTLIB:
        logger.error("matplotlib未安装，无法生成可视化图表")
        logger.error("请安装: pip install matplotlib seaborn")
        sys.exit(1)

    # 初始化可视化器
    visualizer = DEGVisualizer(style=args.style)

    # 加载数据
    success = visualizer.load_data(
        deg_file=args.input,
        count_file=args.counts,
        metadata_file=args.metadata
    )

    if not success:
        logger.error("数据加载失败")
        sys.exit(1)

    # 生成图表
    if args.all_plots:
        generated_files = visualizer.generate_all_plots(
            output_dir=args.output,
            log2fc_threshold=args.log2fc_threshold,
            padj_threshold=args.padj_threshold,
            top_n=args.top_n,
            group_col=args.group_col
        )

        logger.info(f"差异表达可视化完成，生成 {len(generated_files)} 个文件")
    else:
        logger.info("未选择生成图表，使用 --all-plots 生成所有图表")

    logger.info("差异表达可视化流程完成")


if __name__ == "__main__":
    main()