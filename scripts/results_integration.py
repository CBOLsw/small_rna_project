#!/usr/bin/env python3
"""
结果整合脚本
用于自动收集各模块分析结果，生成综合报告

功能：
1. 收集所有模块的分析结果
2. 生成结果摘要统计
3. 创建综合报告目录结构
4. 整合图表和数据
5. 生成HTML报告

使用方法：
    python scripts/results_integration.py --config config/config.yaml
"""

import os
import sys
import argparse
import logging
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List
import json

# 配置日志
def setup_logging() -> None:
    """配置日志系统"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def load_config(config_file: str) -> Dict[str, Any]:
    """加载配置文件"""
    logger = logging.getLogger(__name__)
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        logger.info(f"配置文件已加载: {config_file}")
        return config
    except Exception as e:
        logger.error(f"加载配置文件失败: {e}")
        sys.exit(1)


def collect_qc_results(config: Dict[str, Any]) -> pd.DataFrame:
    """收集质量控制结果"""
    logger = logging.getLogger(__name__)
    qc_dir = os.path.join(config['directories']['results'], 'qc')

    qc_results = []

    # 查找所有 fastqc 结果
    for file in Path(qc_dir).glob('*.html'):
        if '_fastqc.html' in str(file):
            sample = str(file.name).replace('_fastqc.html', '')
            qc_results.append({
                'sample': sample,
                'fastqc_html': str(file),
                'status': 'completed'
            })

    # 加载 qc_summary.csv
    summary_file = os.path.join(qc_dir, 'qc_summary.csv')
    if os.path.exists(summary_file):
        try:
            qc_summary = pd.read_csv(summary_file)
            logger.info(f"质量控制结果已收集: {len(qc_summary)} 个样本")
            return qc_summary
        except Exception as e:
            logger.warning(f"无法读取QC摘要: {e}")

    return pd.DataFrame(qc_results)


def collect_alignment_results(config: Dict[str, Any]) -> pd.DataFrame:
    """收集序列比对结果"""
    logger = logging.getLogger(__name__)
    alignment_dir = os.path.join(config['directories']['results'], 'alignment')

    alignment_files = []

    # 查找所有比对统计文件
    for file in Path(alignment_dir).glob('*_alignment_stats.csv'):
        if '_alignment_stats.csv' in str(file):
            sample = str(file.name).replace('_alignment_stats.csv', '')
            alignment_files.append({
                'sample': sample,
                'stats_file': str(file),
                'bam_file': os.path.join(alignment_dir, f"{sample}.sorted.bam"),
                'bai_file': os.path.join(alignment_dir, f"{sample}.sorted.bam.bai")
            })

    # 加载 alignment_summary.csv
    summary_file = os.path.join(alignment_dir, 'alignment_summary.csv')
    if os.path.exists(summary_file):
        try:
            alignment_summary = pd.read_csv(summary_file)
            logger.info(f"序列比对结果已收集: {len(alignment_summary)} 个样本")
            return alignment_summary
        except Exception as e:
            logger.warning(f"无法读取比对摘要: {e}")

    return pd.DataFrame(alignment_files)


def collect_expression_results(config: Dict[str, Any]) -> Dict[str, Any]:
    """收集基因计数和差异表达结果"""
    logger = logging.getLogger(__name__)
    counts_dir = os.path.join(config['directories']['results'], 'counts')
    de_dir = os.path.join(config['directories']['results'], 'differential_expression')

    results = {}

    # 基因计数结果
    counts_file = os.path.join(counts_dir, 'gene_counts.csv')
    if os.path.exists(counts_file):
        try:
            counts_df = pd.read_csv(counts_file)
            results['gene_counts'] = counts_df
            logger.info(f"基因计数结果已收集: {counts_df.shape[0]} 个基因")
        except Exception as e:
            logger.warning(f"无法读取基因计数: {e}")

    # 差异表达结果
    de_file = os.path.join(de_dir, 'deseq2_results.csv')
    if os.path.exists(de_file):
        try:
            de_df = pd.read_csv(de_file)
            results['deseq2_results'] = de_df
            logger.info(f"差异表达结果已收集: {de_df.shape[0]} 个基因")
        except Exception as e:
            logger.warning(f"无法读取差异表达结果: {e}")

    # 过滤后的差异表达基因
    filtered_file = os.path.join(de_dir, 'filtered_degs.csv')
    if os.path.exists(filtered_file):
        try:
            filtered_df = pd.read_csv(filtered_file)
            results['filtered_degs'] = filtered_df
            logger.info(f"过滤后的差异表达基因已收集: {filtered_df.shape[0]} 个基因")
        except Exception as e:
            logger.warning(f"无法读取过滤后的DEGs: {e}")

    return results


def collect_motif_results(config: Dict[str, Any]) -> Dict[str, Any]:
    """收集motif分析结果"""
    logger = logging.getLogger(__name__)
    motif_dir = os.path.join(config['directories']['results'], 'small_rna_motif')

    results = {}

    # MEME结果
    meme_dir = os.path.join(motif_dir, 'meme_results')
    if os.path.exists(meme_dir):
        results['meme_dir'] = meme_dir
        summary_file = os.path.join(meme_dir, 'meme_summary.json')
        if os.path.exists(summary_file):
            try:
                with open(summary_file, 'r') as f:
                    meme_summary = json.load(f)
                results['meme_summary'] = meme_summary
                results['filtered_motifs'] = meme_summary.get('motifs', [])
                logger.info(f"MEME结果已收集: {meme_summary.get('motifs_found', 0)} 个motif")
            except Exception as e:
                logger.warning(f"无法读取MEME摘要: {e}")

    # miRNA reads FASTA
    mirna_fasta = os.path.join(motif_dir, 'mirna_reads.fasta')
    if os.path.exists(mirna_fasta):
        results['mirna_fasta'] = mirna_fasta

    return results


def generate_comprehensive_report(config: Dict[str, Any],
                                 report_dir: str,
                                 qc_data: pd.DataFrame,
                                 alignment_data: pd.DataFrame,
                                 expression_data: Dict[str, Any],
                                 motif_data: Dict[str, Any]) -> str:
    """生成综合报告"""
    logger = logging.getLogger(__name__)

    # 生成HTML报告
    html_content = generate_html_report(config, qc_data, alignment_data, expression_data, motif_data)

    report_file = os.path.join(report_dir, 'comprehensive_report.html')
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

    logger.info(f"综合报告已生成: {report_file}")

    # 生成报告摘要
    generate_report_summary(config, report_dir, qc_data, alignment_data, expression_data, motif_data)

    return report_file


def generate_html_report(config: Dict[str, Any],
                        qc_data: pd.DataFrame,
                        alignment_data: pd.DataFrame,
                        expression_data: Dict[str, Any],
                        motif_data: Dict[str, Any]) -> str:
    """生成HTML报告内容"""
    # 简单的HTML模板
    html = f"""
    <!DOCTYPE html>
    <html lang="zh-CN">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Small RNA测序分析报告 - {config['project_name']}</title>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 0 10px rgba(0,0,0,0.1);
            }}
            h1 {{
                color: #2c3e50;
                text-align: center;
            }}
            .section {{
                margin: 30px 0;
                padding: 20px;
                border: 1px solid #eee;
                border-radius: 8px;
                background-color: #fafafa;
            }}
            h2 {{
                color: #34495e;
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
            }}
            .stats {{
                display: flex;
                justify-content: space-around;
                margin: 20px 0;
            }}
            .stat-card {{
                background-color: white;
                padding: 20px;
                border-radius: 8px;
                box-shadow: 0 2px 5px rgba(0,0,0,0.1);
                text-align: center;
                width: 200px;
            }}
            .stat-value {{
                font-size: 24px;
                font-weight: bold;
                color: #3498db;
            }}
            .stat-label {{
                margin-top: 5px;
                color: #666;
            }}
            table {{
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
            }}
            table, th, td {{
                border: 1px solid #ddd;
            }}
            th, td {{
                padding: 12px;
                text-align: left;
            }}
            th {{
                background-color: #f2f2f2;
                font-weight: bold;
            }}
            .status {{
                padding: 5px 10px;
                border-radius: 4px;
                font-size: 12px;
                font-weight: bold;
            }}
            .status.completed {{
                background-color: #d4edda;
                color: #155724;
            }}
            .status.pending {{
                background-color: #fff3cd;
                color: #856404;
            }}
            .file-link {{
                color: #007bff;
                text-decoration: none;
            }}
            .file-link:hover {{
                text-decoration: underline;
            }}
            .visualization {{
                display: flex;
                flex-wrap: wrap;
                gap: 20px;
            }}
            .visualization img {{
                max-width: 45%;
                border: 1px solid #ddd;
                border-radius: 4px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>Small RNA测序分析报告</h1>
            <p style="text-align: center; color: #666;">
                项目: {config['project_name']} |
                生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </p>

            <!-- 项目概览 -->
            <div class="section">
                <h2>项目概览</h2>
                <div class="stats">
                    <div class="stat-card">
                        <div class="stat-value">{len(qc_data)}</div>
                        <div class="stat-label">样本数量</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">
                            {len(expression_data.get('gene_counts', [])) if 'gene_counts' in expression_data else 0}
                        </div>
                        <div class="stat-label">基因数量</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">
                            {len(expression_data.get('filtered_degs', [])) if 'filtered_degs' in expression_data else 0}
                        </div>
                        <div class="stat-label">差异基因</div>
                    </div>
                    <div class="stat-card">
                        <div class="stat-value">
                            {len(motif_data.get('filtered_motifs', [])) if 'filtered_motifs' in motif_data else 0}
                        </div>
                        <div class="stat-label">Motif数量</div>
                    </div>
                </div>
            </div>

            <!-- 数据质量控制 -->
            <div class="section">
                <h2>数据质量控制</h2>
                <table>
                    <tr>
                        <th>样本</th>
                        <th>状态</th>
                        <th>详细报告</th>
                    </tr>
    """

    for _, row in qc_data.head(10).iterrows():  # 只显示前10个样本
        sample = row.get('sample', str(row.get('sample_name', '未知')))
        html += f"""
                    <tr>
                        <td>{sample}</td>
                        <td><span class="status completed">完成</span></td>
                        <td><a href="{row.get('fastqc_html', '#')}" class="file-link" target="_blank">FastQC报告</a></td>
                    </tr>
        """

    if len(qc_data) > 10:
        html += f"<tr><td colspan='3'>还有 {len(qc_data) - 10} 个样本...</td></tr>"

    html += """
                </table>
            </div>

            <!-- 序列比对 -->
            <div class="section">
                <h2>序列比对</h2>
                <table>
                    <tr>
                        <th>样本</th>
                        <th>比对率</th>
                        <th>唯一比对</th>
                        <th>多重比对</th>
                    </tr>
    """

    for _, row in alignment_data.head(10).iterrows():  # 只显示前10个样本
        sample = row.get('sample', '未知')
        mapping_rate = row.get('mapping_rate', 'N/A')
        unique_mapping = row.get('unique_mapping', 'N/A')
        multiple_mapping = row.get('multiple_mapping', 'N/A')

        html += f"""
                    <tr>
                        <td>{sample}</td>
                        <td>{mapping_rate}</td>
                        <td>{unique_mapping}</td>
                        <td>{multiple_mapping}</td>
                    </tr>
        """

    if len(alignment_data) > 10:
        html += f"<tr><td colspan='4'>还有 {len(alignment_data) - 10} 个样本...</td></tr>"

    html += """
                </table>
            </div>

            <!-- 差异表达分析 -->
            <div class="section">
                <h2>差异表达分析</h2>
                <div class="visualization">
    """

    de_dir = os.path.join(config['directories']['results'], 'differential_expression')
    volcano_plot = os.path.join(de_dir, 'volcano_plot.png')
    if os.path.exists(volcano_plot):
        html += f"<img src='{volcano_plot}' alt='火山图'>"

    heatmap = os.path.join(de_dir, 'heatmap.png')
    if os.path.exists(heatmap):
        html += f"<img src='{heatmap}' alt='热图'>"

    html += """
                </div>
                <p>差异表达基因数量:
                <strong>
    """

    html += f"{len(expression_data.get('filtered_degs', [])) if 'filtered_degs' in expression_data else 0}"

    html += """
                </strong>
                </p>
            </div>

            <!-- Motif分析 -->
            <div class="section">
                <h2>Motif分析</h2>
                <p>发现的motif数量:
                <strong>
    """

    html += f"{len(motif_data.get('filtered_motifs', [])) if 'filtered_motifs' in motif_data else 0}"

    html += """
                </strong>
                </p>

                <div class="visualization">
    """

    motif_dir = os.path.join(config['directories']['results'], 'small_rna_motif')
    meme_html = os.path.join(motif_dir, 'meme_results', 'meme.html')
    if os.path.exists(meme_html):
        html += f"<p>MEME结果: <a href='{meme_html}' target='_blank'>meme.html</a></p>"

    html += """
                </div>
            </div>

            <!-- 报告信息 -->
            <div class="section">
                <h2>报告信息</h2>
                <table>
                    <tr>
                        <th>项目</th>
                        <td><a href="https://github.com/your-repo" target="_blank">
                        small_rna_analysis_gao_pal</a></td>
                    </tr>
                    <tr>
                        <th>分析流程</th>
                        <td>Snakemake</td>
                    </tr>
                    <tr>
                        <th>参考基因组</th>
                        <td>hg38</td>
                    </tr>
                    <tr>
                        <th>motif发现工具</th>
                        <td>MEME Suite (de novo)</td>
                    </tr>
                </table>
            </div>
        </div>
    </body>
    </html>
    """

    return html


def generate_report_summary(config: Dict[str, Any],
                          report_dir: str,
                          qc_data: pd.DataFrame,
                          alignment_data: pd.DataFrame,
                          expression_data: Dict[str, Any],
                          motif_data: Dict[str, Any]) -> str:
    """生成报告摘要"""
    logger = logging.getLogger(__name__)

    summary = {
        'report_generation_time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'project_name': config['project_name'],
        'samples_analyzed': len(qc_data),
        'qc_results': len(qc_data),
        'alignment_results': len(alignment_data),
        'genes_expressed': len(expression_data.get('gene_counts', [])) if 'gene_counts' in expression_data else 0,
        'differential_genes': len(expression_data.get('filtered_degs', [])) if 'filtered_degs' in expression_data else 0,
        'motifs_found': len(motif_data.get('filtered_motifs', [])) if 'filtered_motifs' in motif_data else 0,
        'config': {
            'sample_groups': config['samples']['groups'],
            'reference_genome': str(os.path.basename(config['reference']['genome_fasta'])),
            'alignment_tool': 'Bowtie2',
            'expression_analysis': 'DESeq2',
            'motif_discovery': 'MEME Suite'
        }
    }

    summary_file = os.path.join(report_dir, 'report_summary.json')
    with open(summary_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    logger.info(f"报告摘要已生成: {summary_file}")

    return summary_file


def copy_visualization_files(config: Dict[str, Any], report_dir: str) -> None:
    """复制可视化文件到报告目录"""
    logger = logging.getLogger(__name__)
    results_dir = config['directories']['results']
    visualization_dir = os.path.join(report_dir, 'visualizations')

    Path(visualization_dir).mkdir(exist_ok=True)

    # 复制差异表达分析图表
    de_dir = os.path.join(results_dir, 'differential_expression')
    for file in Path(de_dir).glob('*.png'):
        target_path = os.path.join(visualization_dir, os.path.basename(file))
        if not os.path.exists(target_path):
            try:
                import shutil
                shutil.copy2(str(file), str(target_path))
            except Exception as e:
                logger.warning(f"无法复制文件 {file}: {e}")

    # 复制motif可视化图表
    motif_dir = os.path.join(results_dir, 'small_rna_motif', 'meme_results')
    if os.path.exists(motif_dir):
        for file in Path(motif_dir).glob('*.html'):
            target_path = os.path.join(visualization_dir, os.path.basename(file))
            if not os.path.exists(target_path):
                try:
                    import shutil
                    shutil.copy2(str(file), str(target_path))
                except Exception as e:
                    logger.warning(f"无法复制文件 {file}: {e}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="结果整合脚本 - 用于自动收集各模块分析结果，生成综合报告"
    )

    parser.add_argument('--config', required=True,
                       help='配置文件路径 (config/config.yaml)')
    parser.add_argument('--output',
                       help='输出报告目录（默认使用 timestamp-based 目录）')

    args = parser.parse_args()

    setup_logging()
    logger = logging.getLogger(__name__)

    logger.info("=" * 60)
    logger.info("结果整合脚本开始")
    logger.info("=" * 60)

    # 加载配置
    config = load_config(args.config)

    # 创建报告目录（使用results/summary作为默认输出）
    report_dir = args.output if args.output else os.path.join(config['directories']['results'], 'summary')
    Path(report_dir).mkdir(parents=True, exist_ok=True)

    # 收集各模块结果
    logger.info("开始收集分析结果...")

    qc_data = collect_qc_results(config)
    alignment_data = collect_alignment_results(config)
    expression_data = collect_expression_results(config)
    motif_data = collect_motif_results(config)

    logger.info("分析结果收集完成")

    # 生成报告
    logger.info("开始生成报告...")
    report_file = generate_comprehensive_report(
        config, report_dir, qc_data, alignment_data,
        expression_data, motif_data
    )

    # 复制可视化文件
    copy_visualization_files(config, report_dir)

    logger.info("=" * 60)
    logger.info("结果整合完成！")
    logger.info(f"报告位置: {report_file}")
    logger.info(f"报告目录: {report_dir}")
    logger.info("=" * 60)

    return 0


if __name__ == '__main__':
    sys.exit(main())
