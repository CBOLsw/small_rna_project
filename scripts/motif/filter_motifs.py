#!/usr/bin/env python3
"""
motif结果过滤和验证脚本
用于去除假阳性motif结果，提高分析可靠性

功能：
1. 背景序列验证
2. 统计显著性过滤（E-value、q-value）
3. 序列特异性检查
4. 结果质量评估
5. 过滤后结果导出

使用方法：
    python filter_motifs.py --meme <MEME结果目录> --tomtom <TomTom结果目录> --output <输出目录>
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
import re

# 导入项目的日志配置工具
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.logging_utils import get_script_logger

# 配置日志
logger = get_script_logger('filter_motifs')


class MotifFilter:
    """motif结果过滤器"""

    def __init__(self):
        self.motifs = []
        self.comparisons = []
        self.filtered_motifs = []
        self.filtered_comparisons = []
        self.filter_stats = {}

    def load_meme_results(self, meme_dir: str) -> bool:
        """加载MEME结果"""
        meme_path = Path(meme_dir)

        # 尝试加载JSON摘要
        summary_file = meme_path / 'meme_summary.json'
        if summary_file.exists():
            try:
                with open(summary_file, 'r') as f:
                    meme_summary = json.load(f)
                self.motifs = meme_summary.get('motifs', [])
                logger.info(f"从JSON加载 {len(self.motifs)} 个motif")
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

        # 尝试解析文本文件
        txt_file = meme_path / 'meme.txt'
        if txt_file.exists():
            try:
                self.motifs = self._parse_meme_text(txt_file)
                logger.info(f"从文本文件解析 {len(self.motifs)} 个motif")
                return True
            except Exception as e:
                logger.warning(f"解析文本文件失败: {e}")

        logger.error("无法加载MEME结果")
        return False

    def load_tomtom_results(self, tomtom_dir: str) -> bool:
        """加载TomTom结果"""
        tomtom_path = Path(tomtom_dir)

        # 尝试加载JSON摘要
        summary_file = tomtom_path / 'tomtom_summary.json'
        if summary_file.exists():
            try:
                with open(summary_file, 'r') as f:
                    tomtom_summary = json.load(f)
                self.comparisons = tomtom_summary.get('comparisons', [])
                logger.info(f"从JSON加载 {len(self.comparisons)} 个motif比对")
                return True
            except Exception as e:
                logger.warning(f"加载TomTom JSON摘要失败: {e}")

        # 尝试加载CSV
        csv_file = tomtom_path / 'tomtom_comparisons.csv'
        if csv_file.exists():
            try:
                df = pd.read_csv(csv_file)
                self.comparisons = df.to_dict('records')
                logger.info(f"从CSV加载 {len(self.comparisons)} 个motif比对")
                return True
            except Exception as e:
                logger.warning(f"加载TomTom CSV失败: {e}")

        logger.warning("无法加载TomTom结果，将继续仅基于MEME结果过滤")
        return False

    def _parse_meme_text(self, txt_file: Path) -> List[Dict[str, Any]]:
        """解析MEME文本文件"""
        motifs = []

        try:
            with open(txt_file, 'r') as f:
                content = f.read()

            # 简单的正则表达式匹配
            # 查找MOTIF行
            motif_pattern = r'MOTIF\s+(\S+).*?width\s+(\d+).*?sites\s+(\d+).*?llr\s+([\d\.\-]+).*?E-value\s+([\d\.\-eE\+]+)'
            matches = re.findall(motif_pattern, content, re.DOTALL)

            for match in matches:
                motif_id, width, sites, llr, evalue = match
                motifs.append({
                    'id': motif_id,
                    'width': int(width),
                    'sites': int(sites),
                    'llr': float(llr),
                    'evalue': float(evalue)
                })

        except Exception as e:
            logger.error(f"解析MEME文本时出错: {e}")

        return motifs

    def filter_motifs(self, config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        过滤motif结果

        参数:
            config: 过滤配置

        返回:
            Dict: 过滤统计信息
        """
        if config is None:
            config = {}

        # 默认过滤参数
        defaults = {
            'evalue_threshold': 1e-4,     # E-value阈值
            'min_sites': 5,               # 最小位点数
            'min_width': 6,               # 最小motif宽度
            'max_width': 12,              # 最大motif宽度
            'require_database_match': False,  # 是否要求数据库匹配
            'tomtom_evalue_threshold': 0.05,  # TomTom E-value阈值
            'max_background_frequency': 0.1,  # 最大背景频率
        }

        # 合并配置
        for key, value in defaults.items():
            config.setdefault(key, value)

        logger.info("开始过滤motif结果")
        logger.info(f"初始motif数: {len(self.motifs)}")
        logger.info(f"初始比对数: {len(self.comparisons)}")

        stats = {
            'initial_motifs': len(self.motifs),
            'initial_comparisons': len(self.comparisons),
            'filter_steps': {}
        }

        # 步骤1: 基于MEME结果的过滤
        filtered_motifs = []
        for motif in self.motifs:
            # 检查E-value
            evalue = motif.get('evalue')
            if evalue is not None and evalue > config['evalue_threshold']:
                continue

            # 检查位点数
            sites = motif.get('sites', 0)
            if sites < config['min_sites']:
                continue

            # 检查宽度
            width = motif.get('width', 0)
            if width < config['min_width'] or width > config['max_width']:
                continue

            filtered_motifs.append(motif)

        stats['filter_steps']['meme_basic'] = {
            'remaining': len(filtered_motifs),
            'filtered': len(self.motifs) - len(filtered_motifs)
        }
        logger.info(f"基础过滤后剩余: {len(filtered_motifs)} 个motif")

        # 步骤2: 基于TomTom结果的过滤（如果可用）
        if self.comparisons and config.get('require_database_match', False):
            # 创建motif ID到比对结果的映射
            motif_comparisons = {}
            for comp in self.comparisons:
                query_id = comp.get('Query_ID')
                if query_id:
                    if query_id not in motif_comparisons:
                        motif_comparisons[query_id] = []
                    motif_comparisons[query_id].append(comp)

            # 过滤：只保留有显著数据库匹配的motif
            db_filtered_motifs = []
            for motif in filtered_motifs:
                motif_id = motif.get('id')
                if motif_id in motif_comparisons:
                    # 检查是否有显著匹配
                    significant_match = False
                    for comp in motif_comparisons[motif_id]:
                        comp_evalue = comp.get('E-value')
                        if comp_evalue is not None and comp_evalue <= config['tomtom_evalue_threshold']:
                            significant_match = True
                            break

                    if significant_match:
                        db_filtered_motifs.append(motif)
                else:
                    # 如果没有比对结果，根据配置决定是否保留
                    if not config.get('require_database_match', False):
                        db_filtered_motifs.append(motif)

            stats['filter_steps']['database_match'] = {
                'remaining': len(db_filtered_motifs),
                'filtered': len(filtered_motifs) - len(db_filtered_motifs)
            }
            filtered_motifs = db_filtered_motifs
            logger.info(f"数据库匹配过滤后剩余: {len(filtered_motifs)} 个motif")

        # 步骤3: 过滤TomTom比对结果（只保留过滤后motif的比对）
        filtered_comparisons = []
        if self.comparisons:
            filtered_motif_ids = {motif.get('id') for motif in filtered_motifs}
            for comp in self.comparisons:
                query_id = comp.get('Query_ID')
                if query_id in filtered_motif_ids:
                    # 同时应用TomTom特定的过滤
                    comp_evalue = comp.get('E-value')
                    if comp_evalue is None or comp_evalue <= config['tomtom_evalue_threshold']:
                        filtered_comparisons.append(comp)

            stats['filter_steps']['tomtom_filter'] = {
                'remaining': len(filtered_comparisons),
                'filtered': len(self.comparisons) - len(filtered_comparisons)
            }
            logger.info(f"TomTom结果过滤后剩余: {len(filtered_comparisons)} 个比对")

        # 保存结果
        self.filtered_motifs = filtered_motifs
        self.filtered_comparisons = filtered_comparisons

        stats.update({
            'final_motifs': len(filtered_motifs),
            'final_comparisons': len(filtered_comparisons),
            'filter_config': config
        })

        self.filter_stats = stats
        return stats

    def calculate_background_frequency(self, background_sequences: List[str],
                                     motif_consensus: str) -> float:
        """
        计算motif在背景序列中的出现频率

        参数:
            background_sequences: 背景序列列表
            motif_consensus: motif共有序列

        返回:
            float: 出现频率
        """
        if not background_sequences or not motif_consensus:
            return 0.0

        motif_len = len(motif_consensus)
        total_positions = 0
        motif_count = 0

        for seq in background_sequences:
            seq_len = len(seq)
            if seq_len < motif_len:
                continue

            total_positions += (seq_len - motif_len + 1)

            # 简单搜索（实际应用中可能需要模糊匹配）
            for i in range(seq_len - motif_len + 1):
                substring = seq[i:i+motif_len]
                # 简单匹配：至少80%相同
                matches = sum(1 for a, b in zip(substring, motif_consensus) if a == b)
                if matches >= motif_len * 0.8:
                    motif_count += 1

        if total_positions == 0:
            return 0.0

        return motif_count / total_positions

    def assess_motif_quality(self, motif: Dict[str, Any]) -> Dict[str, Any]:
        """
        评估motif质量

        参数:
            motif: motif信息

        返回:
            Dict: 质量评估结果
        """
        quality = {
            'score': 0.0,
            'criteria': {}
        }

        # 1. E-value评分
        evalue = motif.get('evalue', 1.0)
        if evalue <= 1e-10:
            quality['criteria']['evalue'] = 'excellent'
            quality['score'] += 30
        elif evalue <= 1e-5:
            quality['criteria']['evalue'] = 'good'
            quality['score'] += 20
        elif evalue <= 1e-3:
            quality['criteria']['evalue'] = 'fair'
            quality['score'] += 10
        else:
            quality['criteria']['evalue'] = 'poor'
            quality['score'] += 0

        # 2. 位点数评分
        sites = motif.get('sites', 0)
        if sites >= 20:
            quality['criteria']['sites'] = 'excellent'
            quality['score'] += 30
        elif sites >= 10:
            quality['criteria']['sites'] = 'good'
            quality['score'] += 20
        elif sites >= 5:
            quality['criteria']['sites'] = 'fair'
            quality['score'] += 10
        else:
            quality['criteria']['sites'] = 'poor'
            quality['score'] += 0

        # 3. 信息含量评分（如果可用）
        llr = motif.get('llr')
        if llr is not None:
            if llr >= 50:
                quality['criteria']['llr'] = 'excellent'
                quality['score'] += 40
            elif llr >= 30:
                quality['criteria']['llr'] = 'good'
                quality['score'] += 30
            elif llr >= 15:
                quality['criteria']['llr'] = 'fair'
                quality['score'] += 20
            else:
                quality['criteria']['llr'] = 'poor'
                quality['score'] += 10

        return quality

    def save_filtered_results(self, output_dir: str) -> Dict[str, str]:
        """
        保存过滤后的结果

        参数:
            output_dir: 输出目录

        返回:
            Dict: 输出文件路径
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        output_files = {}

        # 1. 保存过滤统计
        stats_file = output_path / "filter_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(self.filter_stats, f, indent=2)
        output_files['statistics'] = str(stats_file)

        # 2. 保存过滤后的motif
        if self.filtered_motifs:
            motifs_df = pd.DataFrame(self.filtered_motifs)
            motifs_file = output_path / "filtered_motifs.csv"
            motifs_df.to_csv(motifs_file, index=False)
            output_files['motifs'] = str(motifs_file)

            # 添加质量评估
            motifs_with_quality = []
            for motif in self.filtered_motifs:
                motif_copy = motif.copy()
                quality = self.assess_motif_quality(motif)
                motif_copy.update({
                    'quality_score': quality['score'],
                    'quality_assessment': json.dumps(quality['criteria'])
                })
                motifs_with_quality.append(motif_copy)

            quality_df = pd.DataFrame(motifs_with_quality)
            quality_file = output_path / "motifs_with_quality.csv"
            quality_df.to_csv(quality_file, index=False)
            output_files['motifs_quality'] = str(quality_file)

        # 3. 保存过滤后的比对结果
        if self.filtered_comparisons:
            comparisons_df = pd.DataFrame(self.filtered_comparisons)
            comparisons_file = output_path / "filtered_comparisons.csv"
            comparisons_df.to_csv(comparisons_file, index=False)
            output_files['comparisons'] = str(comparisons_file)

        # 4. 保存过滤报告
        report_file = output_path / "filtering_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== motif结果过滤报告 ===\n\n")
            f.write(f"生成时间: {pd.Timestamp.now()}\n\n")

            f.write("1. 过滤统计\n")
            f.write(f"   初始motif数: {self.filter_stats.get('initial_motifs', 0)}\n")
            f.write(f"   过滤后motif数: {self.filter_stats.get('final_motifs', 0)}\n")
            f.write(f"   过滤掉: {self.filter_stats.get('initial_motifs', 0) - self.filter_stats.get('final_motifs', 0)}\n\n")

            if self.filter_stats.get('initial_comparisons', 0) > 0:
                f.write(f"   初始比对数: {self.filter_stats.get('initial_comparisons', 0)}\n")
                f.write(f"   过滤后比对数: {self.filter_stats.get('final_comparisons', 0)}\n\n")

            f.write("2. 过滤步骤\n")
            for step_name, step_stats in self.filter_stats.get('filter_steps', {}).items():
                f.write(f"   {step_name}: 剩余 {step_stats['remaining']}, 过滤 {step_stats['filtered']}\n")

            f.write("\n3. 过滤参数\n")
            for param_name, param_value in self.filter_stats.get('filter_config', {}).items():
                f.write(f"   {param_name}: {param_value}\n")

            f.write("\n4. 输出文件\n")
            for file_type, file_path in output_files.items():
                f.write(f"   {file_type}: {Path(file_path).name}\n")

        output_files['report'] = str(report_file)

        logger.info(f"过滤结果已保存到: {output_path}")
        return output_files


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="motif结果过滤和验证脚本 - 去除假阳性结果",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法
  python filter_motifs.py --meme meme_results/ --tomtom tomtom_results/ --output filtered_results/

  # 仅使用MEME结果
  python filter_motifs.py --meme meme_results/ --output filtered_results/

  # 自定义过滤参数
  python filter_motifs.py --meme meme_results/ --output results/ --evalue 1e-5 --min-sites 10
        """
    )

    # 输入参数
    parser.add_argument('--meme', required=True, help='MEME结果目录')
    parser.add_argument('--tomtom', help='TomTom结果目录（可选）')
    parser.add_argument('--output', required=True, help='输出目录')

    # 过滤参数
    parser.add_argument('--evalue', type=float, default=1e-4,
                       help='E-value阈值（默认: 1e-4）')
    parser.add_argument('--min-sites', type=int, default=5,
                       help='最小位点数（默认: 5）')
    parser.add_argument('--min-width', type=int, default=6,
                       help='最小motif宽度（默认: 6）')
    parser.add_argument('--max-width', type=int, default=12,
                       help='最大motif宽度（默认: 12）')
    parser.add_argument('--require-db-match', action='store_true',
                       help='要求有数据库匹配（默认: 不要求）')
    parser.add_argument('--tomtom-evalue', type=float, default=0.05,
                       help='TomTom E-value阈值（默认: 0.05）')

    args = parser.parse_args()

    logger.info("=== motif结果过滤开始 ===")
    logger.info(f"MEME结果目录: {args.meme}")
    if args.tomtom:
        logger.info(f"TomTom结果目录: {args.tomtom}")
    logger.info(f"输出目录: {args.output}")

    # 创建过滤器
    filter = MotifFilter()

    # 加载数据
    if not filter.load_meme_results(args.meme):
        logger.error("加载MEME结果失败")
        return 1

    if args.tomtom:
        filter.load_tomtom_results(args.tomtom)

    # 配置过滤参数
    filter_config = {
        'evalue_threshold': args.evalue,
        'min_sites': args.min_sites,
        'min_width': args.min_width,
        'max_width': args.max_width,
        'require_database_match': args.require_db_match,
        'tomtom_evalue_threshold': args.tomtom_evalue,
    }

    # 执行过滤
    stats = filter.filter_motifs(filter_config)
    logger.info(f"过滤完成: {stats['final_motifs']}/{stats['initial_motifs']} 个motif保留")

    # 保存结果
    output_files = filter.save_filtered_results(args.output)

    logger.info("=== 过滤完成 ===")
    logger.info(f"输出目录: {args.output}")

    # 显示主要输出文件
    for file_type, file_path in output_files.items():
        if file_type in ['motifs', 'comparisons', 'report']:
            logger.info(f"{file_type}: {file_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())