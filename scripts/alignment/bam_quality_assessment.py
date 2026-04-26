#!/usr/bin/env python3
"""
BAM文件质量评估脚本

功能：
1. 比对质量分布分析
2. 基因组覆盖度分析
3. 插入片段大小分布（双端数据）
4. 生成质量评估报告和可视化图表

使用方法：
    python bam_quality_assessment.py --input <BAM文件或目录> --output <输出目录>
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

# 尝试导入可视化库
try:
    import matplotlib
    matplotlib.use('Agg')  # 非交互模式

    # 配置matplotlib字体，避免中文字体警告
    matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans']
    matplotlib.rcParams['font.family'] = 'sans-serif'
    matplotlib.rcParams['axes.unicode_minus'] = False

    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_VIS = True
except ImportError:
    HAS_VIS = False
    print("警告: matplotlib/seaborn未安装，将跳过图表生成")

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 压制matplotlib的冗余INFO日志
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('matplotlib.category').setLevel(logging.WARNING)


class BAMQualityAssessor:
    """BAM文件质量评估类"""

    def __init__(self, samtools_path: str = "samtools"):
        """
        初始化质量评估器

        参数:
            samtools_path: samtools可执行文件路径
        """
        self.samtools_path = samtools_path
        self.results = {}

    def check_samtools(self) -> bool:
        """检查samtools是否可用"""
        try:
            result = subprocess.run(
                [self.samtools_path, "--version"],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                version_line = result.stdout.strip().split('\n')[0]
                logger.info(f"samtools版本: {version_line}")
                return True
            else:
                logger.error(f"samtools检查失败: {result.stderr}")
                return False
        except FileNotFoundError:
            logger.error(f"未找到samtools: {self.samtools_path}")
            return False

    def assess_single_bam(self, bam_file: str, output_dir: str,
                         sample_name: Optional[str] = None,
                         config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        评估单个BAM文件质量

        参数:
            bam_file: BAM文件路径
            output_dir: 输出目录
            sample_name: 样本名称
            config: 配置参数

        返回:
            Dict: 质量评估结果
        """
        if config is None:
            config = {}

        bam_path = Path(bam_file)
        if not bam_path.exists():
            logger.error(f"BAM文件不存在: {bam_file}")
            return {'success': False, 'error': '文件不存在'}

        if sample_name is None:
            sample_name = bam_path.stem
            if sample_name.endswith('_sorted'):
                sample_name = sample_name[:-7]

        logger.info(f"评估BAM文件质量: {sample_name}")

        # 准备输出目录
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 输出文件
        stats_file = output_dir / f"{sample_name}_quality_stats.txt"
        json_file = output_dir / f"{sample_name}_quality_metrics.json"

        # 收集所有质量指标
        metrics = {
            'sample': sample_name,
            'bam_file': str(bam_path),
            'success': True
        }

        try:
            # 1. 基本统计
            basic_stats = self._get_basic_stats(bam_file, sample_name)
            metrics.update(basic_stats)

            # 2. 比对质量分布
            mapq_stats = self._analyze_mapping_quality(bam_file, sample_name, output_dir)
            metrics.update(mapq_stats)

            # 3. 基因组覆盖度
            coverage_stats = self._analyze_genome_coverage(bam_file, sample_name, output_dir)
            metrics.update(coverage_stats)

            # 4. 如果是双端数据，分析插入片段大小
            if metrics.get('paired_reads', 0) > 0:
                insert_stats = self._analyze_insert_size(bam_file, sample_name, output_dir)
                metrics.update(insert_stats)

            # 5. 重复率分析
            duplicate_stats = self._analyze_duplicates(bam_file, sample_name)
            metrics.update(duplicate_stats)

            # 保存结果
            self._save_quality_metrics(metrics, stats_file, json_file)

            logger.info(f"BAM质量评估完成: {sample_name}")

        except Exception as e:
            logger.error(f"评估BAM文件时出错: {e}")
            metrics.update({
                'success': False,
                'error': str(e)
            })

        self.results[sample_name] = metrics
        return metrics

    def _get_basic_stats(self, bam_file: str, sample_name: str) -> Dict[str, Any]:
        """获取基本统计信息"""
        stats = {}

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
                elif 'singletons' in line:
                    stats['singletons'] = int(line.split()[0])
                elif 'duplicates' in line:
                    stats['duplicate_reads'] = int(line.split()[0])

            # 计算比率
            if stats.get('total_reads', 0) > 0:
                stats['mapping_rate'] = stats.get('mapped_reads', 0) / stats['total_reads']
                stats['proper_pair_rate'] = stats.get('properly_paired', 0) / max(stats.get('paired_reads', 1), 1)
                stats['duplicate_rate'] = stats.get('duplicate_reads', 0) / stats['total_reads']

        except Exception as e:
            logger.warning(f"获取基本统计时出错: {e}")

        return stats

    def _analyze_mapping_quality(self, bam_file: str, sample_name: str,
                               output_dir: Path) -> Dict[str, Any]:
        """分析比对质量分布"""
        stats = {}

        try:
            # 提取MAPQ值
            cmd = [
                self.samtools_path, "view", bam_file,
                "|", "awk", "'{print $5}'",
                "|", "sort", "-n", "|", "uniq", "-c"
            ]

            # 使用shell执行管道
            full_cmd = f"{self.samtools_path} view {bam_file} | awk '{{print $5}}' | sort -n | uniq -c"
            result = subprocess.run(
                full_cmd,
                shell=True,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                mapq_counts = {}
                total = 0

                for line in lines:
                    if line.strip():
                        count, mapq = line.strip().split()
                        mapq_counts[int(mapq)] = int(count)
                        total += int(count)

                if total > 0:
                    # 计算统计量
                    mapq_values = []
                    for mapq, count in mapq_counts.items():
                        mapq_values.extend([mapq] * count)

                    mapq_array = np.array(mapq_values)
                    stats.update({
                        'mapq_mean': float(np.mean(mapq_array)),
                        'mapq_median': float(np.median(mapq_array)),
                        'mapq_std': float(np.std(mapq_array)),
                        'mapq_min': int(np.min(mapq_array)),
                        'mapq_max': int(np.max(mapq_array)),
                        'high_quality_reads': int(np.sum(mapq_array >= 30)),
                        'high_quality_rate': float(np.sum(mapq_array >= 30) / len(mapq_array))
                    })

                    # 保存原始分布
                    dist_file = output_dir / f"{sample_name}_mapq_distribution.txt"
                    with open(dist_file, 'w') as f:
                        f.write("MAPQ\tCount\tProportion\n")
                        for mapq in sorted(mapq_counts.keys()):
                            prop = mapq_counts[mapq] / total
                            f.write(f"{mapq}\t{mapq_counts[mapq]}\t{prop:.4f}\n")

                    # 生成可视化（如果可用）
                    if HAS_VIS:
                        self._plot_mapq_distribution(mapq_counts, sample_name, output_dir)

        except Exception as e:
            logger.warning(f"分析比对质量时出错: {e}")

        return stats

    def _analyze_genome_coverage(self, bam_file: str, sample_name: str,
                               output_dir: Path) -> Dict[str, Any]:
        """分析基因组覆盖度（使用mosdepth，比samtools depth快10-50倍）"""
        stats = {}

        try:
            prefix = output_dir / f"{sample_name}_mosdepth"
            cmd = ["mosdepth", "-n", str(prefix), bam_file]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False,
                timeout=600
            )

            if result.returncode != 0:
                logger.warning(f"mosdepth运行失败: {result.stderr[:200]}")
                return stats

            # 解析全局覆盖度分布
            dist_file = Path(f"{prefix}.mosdepth.global.dist.txt")
            if dist_file.exists():
                depths = []
                props = []
                with open(dist_file) as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) == 3 and parts[0] == "total":
                            depths.append(int(parts[1]))
                            props.append(float(parts[2]))

                if len(props) > 1:
                    # 基尼系数近似：面积法
                    stats['mean_coverage'] = float(np.trapz(props, depths))
                    stats['median_coverage'] = 0.0
                    for d, p in zip(depths, props):
                        if p >= 0.5:
                            stats['median_coverage'] = float(d)
                            break

                    stats['max_coverage'] = max(depths)

                    # 累积分布中提取关键阈值
                    for d, p in zip(depths, props):
                        if d == 1 and p >= 0:
                            stats['covered_bases_rate'] = p
                        if d >= 10 and 'coverage_10x' not in stats:
                            stats['coverage_10x'] = p
                        if d >= 30 and 'coverage_30x' not in stats:
                            stats['coverage_30x'] = p

                    # 保存分布数据
                    dist_out = output_dir / f"{sample_name}_coverage_dist.txt"
                    with open(dist_out, 'w') as f:
                        f.write("Depth\tProportion≥Depth\n")
                        for d, p in zip(depths, props):
                            f.write(f"{d}\t{p:.6f}\n")

                    # 可视化
                    if HAS_VIS:
                        self._plot_coverage_distribution_mosdepth(
                            depths, props, sample_name, output_dir
                        )

            # 解析染色体级别汇总
            summary_file = Path(f"{prefix}.mosdepth.summary.txt")
            if summary_file.exists():
                chrom_depths = []
                with open(summary_file) as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 4 and parts[0] != "chrom":
                            try:
                                chrom_depths.append(float(parts[3]))
                            except ValueError:
                                pass
                if chrom_depths:
                    stats['mean_coverage'] = float(np.mean(chrom_depths))

            # 清理mosdepth临时文件
            for f in Path(output_dir).glob(f"{sample_name}_mosdepth*"):
                f.unlink()

        except subprocess.TimeoutExpired:
            logger.warning(f"mosdepth超时（{sample_name}），跳过覆盖度分析")
        except FileNotFoundError:
            logger.warning("mosdepth未安装，跳过覆盖度分析")
        except Exception as e:
            logger.warning(f"分析基因组覆盖度时出错: {e}")

        return stats


    def _analyze_insert_size(self, bam_file: str, sample_name: str,
                           output_dir: Path) -> Dict[str, Any]:
        """分析插入片段大小分布（仅适用于双端数据）"""
        stats = {}

        try:
            # 使用samtools stats获取插入片段大小
            cmd = [self.samtools_path, "stats", bam_file]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode == 0:
                lines = result.stdout.strip().split('\n')
                insert_sizes = []

                for line in lines:
                    if line.startswith('ID'):
                        continue
                    if 'insert size' in line and 'average' in line:
                        parts = line.split()
                        for i, part in enumerate(parts):
                            if part == 'average':
                                if i + 1 < len(parts):
                                    stats['insert_size_mean'] = float(parts[i+1])
                                    break
                    elif 'insert size' in line and 'standard deviation' in line:
                        parts = line.split()
                        for i, part in enumerate(parts):
                            if part == 'deviation':
                                if i + 1 < len(parts):
                                    stats['insert_size_std'] = float(parts[i+1])
                                    break

                # 如果samtools stats没有提供，尝试其他方法
                if 'insert_size_mean' not in stats:
                    # 简单方法：提取TLEN字段
                    cmd = f"{self.samtools_path} view {bam_file} | awk '{{if($9>0) print $9}}' | head -1000"
                    result = subprocess.run(
                        cmd,
                        shell=True,
                        capture_output=True,
                        text=True,
                        check=False
                    )

                    if result.returncode == 0 and result.stdout.strip():
                        values = [int(x) for x in result.stdout.strip().split() if x]
                        if values:
                            ins_array = np.array(values)
                            stats.update({
                                'insert_size_mean': float(np.mean(ins_array)),
                                'insert_size_median': float(np.median(ins_array)),
                                'insert_size_std': float(np.std(ins_array)),
                                'insert_size_min': int(np.min(ins_array)),
                                'insert_size_max': int(np.max(ins_array))
                            })

                # 保存插入片段统计
                if 'insert_size_mean' in stats:
                    ins_file = output_dir / f"{sample_name}_insert_size_stats.txt"
                    with open(ins_file, 'w') as f:
                        f.write(f"平均插入片段大小: {stats['insert_size_mean']:.1f} bp\n")
                        f.write(f"插入片段大小标准差: {stats['insert_size_std']:.1f} bp\n")
                        if 'insert_size_median' in stats:
                            f.write(f"中位插入片段大小: {stats['insert_size_median']:.1f} bp\n")

                    # 生成可视化
                    if HAS_VIS:
                        self._plot_insert_size_distribution(values, sample_name, output_dir)

        except Exception as e:
            logger.warning(f"分析插入片段大小时出错: {e}")

        return stats

    def _analyze_duplicates(self, bam_file: str, sample_name: str) -> Dict[str, Any]:
        """分析重复reads"""
        stats = {}

        try:
            # 使用samtools flagstat已经提供了重复reads数
            # 这里可以添加更详细的分析
            pass

        except Exception as e:
            logger.warning(f"分析重复reads时出错: {e}")

        return stats

    def _save_quality_metrics(self, metrics: Dict[str, Any],
                            stats_file: Path, json_file: Path):
        """保存质量指标"""
        try:
            # 保存为文本
            with open(stats_file, 'w') as f:
                f.write(f"=== BAM Quality Assessment Report ===\n\n")
                f.write(f"Sample: {metrics.get('sample', 'unknown')}\n")
                f.write(f"BAM File: {metrics.get('bam_file', 'unknown')}\n\n")

                f.write("Basic Statistics:\n")
                f.write(f"  Total reads: {metrics.get('total_reads', 0):,}\n")
                f.write(f"  Mapped reads: {metrics.get('mapped_reads', 0):,}\n")
                f.write(f"  Mapping rate: {metrics.get('mapping_rate', 0)*100:.2f}%\n")

                if 'paired_reads' in metrics:
                    f.write(f"  Paired reads: {metrics.get('paired_reads', 0):,}\n")
                    f.write(f"  Properly paired rate: {metrics.get('proper_pair_rate', 0)*100:.2f}%\n")

                f.write(f"  Duplicate rate: {metrics.get('duplicate_rate', 0)*100:.2f}%\n\n")

                f.write("Mapping Quality Statistics:\n")
                f.write(f"  MAPQ mean: {metrics.get('mapq_mean', 0):.2f}\n")
                f.write(f"  MAPQ median: {metrics.get('mapq_median', 0):.2f}\n")
                f.write(f"  High quality rate (MAPQ≥30): {metrics.get('high_quality_rate', 0)*100:.2f}%\n\n")

                f.write("Genome Coverage Statistics:\n")
                f.write(f"  Mean coverage: {metrics.get('mean_coverage', 0):.2f}\n")
                f.write(f"  Median coverage: {metrics.get('median_coverage', 0):.2f}\n")
                f.write(f"  Genome covered (≥1x): {metrics.get('covered_bases_rate', 0)*100:.2f}%\n")
                f.write(f"  ≥10x coverage: {metrics.get('coverage_10x', 0)*100:.2f}%\n")
                f.write(f"  ≥30x coverage: {metrics.get('coverage_30x', 0)*100:.2f}%\n")

            # 保存为JSON
            with open(json_file, 'w') as f:
                json.dump(metrics, f, indent=2)

            logger.info(f"质量评估报告已保存: {stats_file}")

        except Exception as e:
            logger.warning(f"保存质量指标时出错: {e}")

    def _plot_mapq_distribution(self, mapq_counts: Dict[int, int],
                              sample_name: str, output_dir: Path):
        """绘制MAPQ分布图"""
        try:
            plt.figure(figsize=(10, 6))

            # 准备数据
            mapqs = sorted(mapq_counts.keys())
            counts = [mapq_counts[m] for m in mapqs]
            proportions = [c / sum(counts) for c in counts]

            # 条形图
            bars = plt.bar(mapqs, proportions, alpha=0.7, color='steelblue')

            # 添加数值标签
            for bar, prop in zip(bars, proportions):
                if prop > 0.05:  # 只显示较大的比例
                    height = bar.get_height()
                    plt.text(bar.get_x() + bar.get_width()/2., height,
                            f'{prop:.1%}', ha='center', va='bottom')

            plt.xlabel('Mapping Quality (MAPQ)')
            plt.ylabel('Proportion')
            plt.title(f'{sample_name} - Mapping Quality Distribution')
            plt.grid(True, alpha=0.3)

            # 保存图片
            plot_file = output_dir / f"{sample_name}_mapq_distribution.png"
            plt.tight_layout()
            plt.savefig(plot_file, dpi=150)
            plt.close()

            logger.info(f"MAPQ分布图已生成: {plot_file}")

        except Exception as e:
            logger.warning(f"绘制MAPQ分布图时出错: {e}")

    def _plot_coverage_distribution(self, depths: np.ndarray,
                                  sample_name: str, output_dir: Path):
        """绘制覆盖度分布图"""
        try:
            plt.figure(figsize=(12, 5))

            # 子图1: 覆盖度直方图
            plt.subplot(1, 2, 1)
            max_depth = min(np.max(depths), 100)  # 限制最大显示深度
            bins = min(50, max_depth)

            plt.hist(depths, bins=bins, alpha=0.7, color='forestgreen',
                    edgecolor='black')
            plt.xlabel('Coverage Depth')
            plt.ylabel('Number of Positions')
            plt.title(f'{sample_name} - Coverage Depth Distribution')
            plt.grid(True, alpha=0.3)

            # 子图2: 累积覆盖度
            plt.subplot(1, 2, 2)
            sorted_depths = np.sort(depths)
            cumulative = np.arange(1, len(sorted_depths) + 1) / len(sorted_depths)

            plt.plot(sorted_depths, cumulative, 'b-', linewidth=2)
            plt.xlabel('Coverage Depth')
            plt.ylabel('Cumulative Proportion')
            plt.title(f'{sample_name} - Cumulative Coverage')
            plt.grid(True, alpha=0.3)

            # 添加关键点标记
            for threshold in [1, 5, 10, 30]:
                idx = np.searchsorted(sorted_depths, threshold, side='right')
                if idx < len(cumulative):
                    plt.plot(threshold, cumulative[idx], 'ro')
                    plt.text(threshold, cumulative[idx] + 0.02,
                            f'{cumulative[idx]:.1%}', ha='center')

            # 保存图片
            plot_file = output_dir / f"{sample_name}_coverage_distribution.png"
            plt.tight_layout()
            plt.savefig(plot_file, dpi=150)
            plt.close()

            logger.info(f"覆盖度分布图已生成: {plot_file}")

        except Exception as e:
            logger.warning(f"绘制覆盖度分布图时出错: {e}")

    def _plot_coverage_distribution_mosdepth(self, depths: List[int],
                                            props: List[float],
                                            sample_name: str, output_dir: Path):
        """绘制mosdepth覆盖度累积分布图"""
        try:
            plt.figure(figsize=(12, 5))

            # 子图1: 覆盖率柱状图（关键阈值）
            plt.subplot(1, 2, 1)
            thresholds = [1, 5, 10, 15, 20, 30, 50]
            t_props = []
            for t in thresholds:
                for d, p in zip(depths, props):
                    if d >= t:
                        t_props.append(p)
                        break
                else:
                    t_props.append(0)

            colors = ['#fee5d9','#fcae91','#fb6a4a','#de2d26','#a50f15','#67000d','#300000']
            bars = plt.bar([str(t) for t in thresholds], t_props,
                          color=colors[:len(thresholds)], alpha=0.8)
            for bar, p in zip(bars, t_props):
                if p > 0.05:
                    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                            f'{p:.0%}', ha='center', va='bottom', fontsize=9)

            plt.xlabel('Coverage Depth Threshold')
            plt.ylabel('Proportion of Genome ≥ Threshold')
            plt.title(f'{sample_name} - Coverage at Thresholds')
            plt.grid(True, alpha=0.3, axis='y')

            # 子图2: 累积覆盖度曲线
            plt.subplot(1, 2, 2)
            plt.plot(depths, props, 'b-', linewidth=2)
            plt.fill_between(depths, props, alpha=0.15, color='blue')
            for t in [1, 5, 10, 30]:
                for d, p in zip(depths, props):
                    if d == t:
                        plt.plot(t, p, 'ro')
                        plt.text(t, p + 0.03, f'{p:.1%}', ha='center')
                    elif d < t < depths[props.index(p) + 1] if props.index(p) + 1 < len(props) else True:
                        pass

            plt.xlabel('Coverage Depth')
            plt.ylabel('Proportion of Genome ≥ Depth')
            plt.title(f'{sample_name} - Cumulative Coverage')
            plt.grid(True, alpha=0.3)
            plt.xlim(0, min(max(depths), 100))

            plot_file = output_dir / f"{sample_name}_coverage_distribution.png"
            plt.tight_layout()
            plt.savefig(plot_file, dpi=150)
            plt.close()
            logger.info(f"覆盖度分布图已生成: {plot_file}")

        except Exception as e:
            logger.warning(f"绘制mosdepth覆盖度分布图时出错: {e}")

    def _plot_insert_size_distribution(self, insert_sizes: List[int],
                                     sample_name: str, output_dir: Path):
        """绘制插入片段大小分布图"""
        try:
            if len(insert_sizes) < 10:
                return

            plt.figure(figsize=(10, 6))

            # 直方图
            plt.hist(insert_sizes, bins=50, alpha=0.7, color='purple',
                    edgecolor='black')

            # 添加统计信息
            mean_size = np.mean(insert_sizes)
            median_size = np.median(insert_sizes)

            plt.axvline(mean_size, color='red', linestyle='--', linewidth=2,
                       label=f'平均值: {mean_size:.1f} bp')
            plt.axvline(median_size, color='blue', linestyle=':', linewidth=2,
                       label=f'中位数: {median_size:.1f} bp')

            plt.xlabel('Insert Size (bp)')
            plt.ylabel('Read Count')
            plt.title(f'{sample_name} - Insert Size Distribution')
            plt.legend()
            plt.grid(True, alpha=0.3)

            # 保存图片
            plot_file = output_dir / f"{sample_name}_insert_size_distribution.png"
            plt.tight_layout()
            plt.savefig(plot_file, dpi=150)
            plt.close()

            logger.info(f"插入片段大小分布图已生成: {plot_file}")

        except Exception as e:
            logger.warning(f"绘制插入片段大小分布图时出错: {e}")

    def generate_summary_report(self, output_dir: str, output_file: str = None) -> Optional[str]:
        """
        生成质量评估汇总报告

        参数:
            output_dir: 输出目录
            output_file: 可选的输出文件路径

        返回:
            str: 报告文件路径
        """
        if not self.results:
            logger.warning("无质量评估结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建DataFrame
        summary_data = []
        for sample, metrics in self.results.items():
            if metrics.get('success'):
                summary_data.append({
                    'sample': sample,
                    'total_reads': metrics.get('total_reads', 0),
                    'mapped_reads': metrics.get('mapped_reads', 0),
                    'mapping_rate': metrics.get('mapping_rate', 0),
                    'mapq_mean': metrics.get('mapq_mean', 0),
                    'high_quality_rate': metrics.get('high_quality_rate', 0),
                    'mean_coverage': metrics.get('mean_coverage', 0),
                    'coverage_10x': metrics.get('coverage_10x', 0),
                    'duplicate_rate': metrics.get('duplicate_rate', 0)
                })

        if not summary_data:
            return None

        df = pd.DataFrame(summary_data)

        if output_file:
            # 直接保存到指定文件
            csv_file = Path(output_file)
            df.to_csv(csv_file, index=False)
            logger.info(f"质量评估汇总已保存: {csv_file}")
            return str(csv_file)

        # 保存CSV
        csv_file = output_dir / "bam_quality_summary.csv"
        df.to_csv(csv_file, index=False)

        # 生成文本报告
        report_file = output_dir / "bam_quality_assessment_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== BAM Quality Assessment Summary Report ===\n\n")
            f.write(f"Number of samples assessed: {len(summary_data)}\n\n")

            f.write("Sample quality metrics:\n")
            for _, row in df.iterrows():
                f.write(f"  - {row['sample']}:\n")
                f.write(f"    Mapping rate: {row['mapping_rate']*100:.1f}%\n")
                f.write(f"    High quality rate: {row['high_quality_rate']*100:.1f}%\n")
                f.write(f"    Mean coverage: {row['mean_coverage']:.2f}\n")
                f.write(f"    ≥10x coverage: {row['coverage_10x']*100:.1f}%\n\n")

            f.write("Overall statistics:\n")
            f.write(f"  Mean mapping rate: {df['mapping_rate'].mean()*100:.1f}%\n")
            f.write(f"  Mean high quality rate: {df['high_quality_rate'].mean()*100:.1f}%\n")
            f.write(f"  Mean coverage: {df['mean_coverage'].mean():.2f}\n")
            f.write(f"  Mean ≥10x coverage: {df['coverage_10x'].mean()*100:.1f}%\n\n")

            f.write(f"报告生成时间: {pd.Timestamp.now()}\n")

        logger.info(f"质量评估汇总报告已生成: {report_file}")
        return str(report_file)


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """加载配置文件"""
    default_config = {
        'samtools_path': 'samtools',
        'threads': 4,
        'min_coverage_depth': 1,
        'max_plot_depth': 100,
        'generate_plots': True
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
        description="BAM文件质量评估脚本"
    )

    parser.add_argument(
        "--input", "-i",
        help="输入BAM文件或包含BAM文件的目录"
    )

    parser.add_argument(
        "--input-dir",
        help="输入目录（与--input功能相同）"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/alignment/quality_assessment",
        help="输出目录 (默认: results/alignment/quality_assessment)"
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

    # 初始化评估器
    assessor = BAMQualityAssessor(samtools_path=config.get('samtools_path', 'samtools'))

    # 检查samtools
    if not assessor.check_samtools():
        logger.error("samtools检查失败，请确保samtools已安装")
        sys.exit(1)

    # 处理输入
    if args.input_dir:
        input_path = Path(args.input_dir)
    elif args.input:
        input_path = Path(args.input)
    else:
        logger.error("必须提供--input或--input-dir参数")
        sys.exit(1)

    # 处理输出路径
    output_path = Path(args.output)
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
            bam_file = row.get('bam_file')

            if pd.notna(bam_file):
                assessor.assess_single_bam(
                    bam_file=bam_file,
                    output_dir=output_dir,
                    sample_name=sample,
                    config=config
                )

    elif input_path.is_file() and input_path.suffix == '.bam':
        # 单个BAM文件
        assessor.assess_single_bam(
            bam_file=str(input_path),
            output_dir=output_dir,
            sample_name=None,
            config=config
        )

    elif input_path.is_dir():
        # 目录下的所有BAM文件
        logger.info(f"扫描目录: {input_path}")
        bam_files = list(input_path.glob("*.bam")) + list(input_path.glob("*.BAM"))

        for bam_file in bam_files:
            assessor.assess_single_bam(
                bam_file=str(bam_file),
                output_dir=output_dir,
                sample_name=None,
                config=config
            )

    else:
        logger.error("不支持的输入类型，请提供BAM文件、包含BAM文件的目录或样本信息CSV")
        sys.exit(1)

    # 生成汇总报告
    if args.summary:
        if is_file_output:
            # 直接保存到指定的CSV文件
            report_file = assessor.generate_summary_report(output_dir, str(output_path))
            if report_file:
                logger.info(f"质量评估完成，报告文件: {report_file}")
        else:
            report_file = assessor.generate_summary_report(output_dir)
            if report_file:
                logger.info(f"质量评估完成，报告文件: {report_file}")
            else:
                logger.warning("未能生成汇总报告")

    logger.info("BAM文件质量评估流程完成")


if __name__ == "__main__":
    main()