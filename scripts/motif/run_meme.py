#!/usr/bin/env python3
"""
MEME motif分析脚本
用于在差异表达small RNA序列中发现富集的序列motif

功能：
1. 运行MEME de novo motif发现
2. 针对small RNA优化参数（短序列、高GC含量）
3. 统计显著性评估（E-value < 1e-4）
4. 结果过滤和导出

使用方法：
    python run_meme.py --fasta <输入序列文件> --output <输出目录>
"""

import os
import sys
import argparse
import subprocess
import logging
import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Any

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def check_meme_installed(meme_path: str = "meme") -> bool:
    """检查MEME是否安装"""
    try:
        result = subprocess.run(
            [meme_path, '-version'],
            capture_output=True,
            text=True,
            check=False
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def run_meme_analysis(fasta_file: str, output_dir: str,
                     width_min: int = 6, width_max: int = 12,
                     max_motifs: int = 10, evalue_threshold: float = 1e-4) -> Dict[str, Any]:
    """
    运行MEME motif分析

    参数:
        fasta_file: 输入FASTA序列文件
        output_dir: 输出目录
        width_min: motif最小宽度
        width_max: motif最大宽度
        max_motifs: 最大motif数
        evalue_threshold: E-value阈值

    返回:
        Dict: 分析结果
    """
    # 检查输入文件
    if not Path(fasta_file).exists():
        logger.error(f"输入文件不存在: {fasta_file}")
        return {'success': False, 'error': '输入文件不存在'}

    # 检查MEME是否安装
    if not check_meme_installed():
        logger.error("MEME未安装或不在PATH中")
        return {'success': False, 'error': 'MEME未安装'}

    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger.info(f"开始MEME分析，输入文件: {fasta_file}")
    logger.info(f"motif宽度: {width_min}-{width_max} nt")
    logger.info(f"最大motif数: {max_motifs}")
    logger.info(f"E-value阈值: {evalue_threshold}")

    try:
        # 构建MEME命令（针对small RNA优化）
        cmd = [
            'meme',
            fasta_file,
            '-o', str(output_path),
            '-dna',              # DNA序列
            '-mod', 'zoops',     # 零或一次出现模型（适合motif发现）
            '-nmotifs', str(max_motifs),
            '-minw', str(width_min),
            '-maxw', str(width_max),
            '-minsites', '5',    # 每个motif最少位点数
            '-maxsites', '100',  # 每个motif最多位点数
            '-evt', str(evalue_threshold),
            '-revcomp',          # 考虑反向互补链
        ]

        logger.info(f"运行MEME命令: {' '.join(cmd)}")

        # 执行MEME命令
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        logger.info("MEME分析完成")

        # 检查输出文件
        result_files = {}
        expected_files = ['meme.html', 'meme.xml', 'meme.txt']

        for file_name in expected_files:
            file_path = output_path / file_name
            if file_path.exists():
                result_files[file_name] = str(file_path)

        # 解析主要结果
        txt_file = output_path / 'meme.txt'
        motifs = parse_meme_text_result(txt_file) if txt_file.exists() else []

        result = {
            'success': True,
            'output_dir': str(output_path),
            'files': result_files,
            'motifs_found': len(motifs),
            'motifs': motifs,
            'parameters': {
                'width_min': width_min,
                'width_max': width_max,
                'max_motifs': max_motifs,
                'evalue_threshold': evalue_threshold
            }
        }

        # 保存结果摘要
        save_result_summary(result, output_path)

        logger.info(f"发现 {len(motifs)} 个显著motif")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"MEME执行失败: {e}")
        logger.error(f"STDERR: {e.stderr[:500]}")
        return {'success': False, 'error': f'MEME执行失败: {e.stderr[:200]}'}
    except Exception as e:
        logger.error(f"运行MEME时出错: {e}")
        return {'success': False, 'error': str(e)}


def parse_meme_text_result(txt_file: Path) -> List[Dict[str, Any]]:
    """解析MEME文本格式结果"""
    motifs = []

    try:
        with open(txt_file, 'r') as f:
            content = f.read()

        # 查找motif部分
        lines = content.split('\n')
        in_motif_section = False

        for i, line in enumerate(lines):
            line = line.strip()

            # 查找MOTIF行
            if line.startswith('MOTIF'):
                parts = line.split()
                if len(parts) >= 7:
                    motif_id = parts[1].strip()

                    # 查找宽度、位点数等信息
                    width_idx = -1
                    sites_idx = -1
                    llr_idx = -1
                    evalue_idx = -1

                    for j, part in enumerate(parts):
                        if part == 'width':
                            width_idx = j + 1
                        elif part == 'sites':
                            sites_idx = j + 1
                        elif part == 'llr':
                            llr_idx = j + 1
                        elif part == 'E-value':
                            evalue_idx = j + 1

                    motif = {
                        'id': motif_id,
                        'width': int(parts[width_idx]) if width_idx > 0 else 0,
                        'sites': int(parts[sites_idx]) if sites_idx > 0 else 0,
                        'llr': float(parts[llr_idx]) if llr_idx > 0 else 0,
                        'evalue': float(parts[evalue_idx]) if evalue_idx > 0 else 1.0,
                    }

                    # 查找后续行的共有序列
                    for j in range(i+1, min(i+10, len(lines))):
                        next_line = lines[j].strip()
                        if 'regular expression' in next_line and 'sites' in next_line:
                            # 提取共有序列
                            if ':' in next_line:
                                consensus = next_line.split(':', 1)[1].strip()
                                motif['consensus'] = consensus
                            break

                    motifs.append(motif)

    except Exception as e:
        logger.warning(f"解析MEME结果时出错: {e}")

    return motifs


def save_result_summary(result: Dict[str, Any], output_dir: Path):
    """保存结果摘要"""
    # JSON摘要
    summary_file = output_dir / "meme_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(result, f, indent=2)

    # CSV表格
    if result.get('motifs'):
        df_data = []
        for motif in result['motifs']:
            row = {
                'motif_id': motif.get('id', ''),
                'width': motif.get('width', 0),
                'sites': motif.get('sites', 0),
                'llr': motif.get('llr', 0),
                'evalue': motif.get('evalue', 0),
                'consensus': motif.get('consensus', '')
            }
            df_data.append(row)

        df = pd.DataFrame(df_data)
        csv_file = output_dir / "motifs_summary.csv"
        df.to_csv(csv_file, index=False)

    logger.info(f"结果摘要已保存: {summary_file}")


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="MEME motif分析脚本 - 用于small RNA序列motif发现",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  python run_meme.py --fasta gene_sequences.fasta --output meme_results/
  python run_meme.py --fasta sequences.fasta --output results/ --width 6,12 --nmotifs 10 --evalth 1e-4
        """
    )

    # 必需参数
    parser.add_argument('--fasta', required=True,
                       help='输入FASTA序列文件（来自extract_sequences.py的输出）')
    parser.add_argument('--output', required=True,
                       help='输出目录（将在此目录下创建MEME结果）')

    # 可选参数
    parser.add_argument('--width', type=str, default='6,12',
                       help='motif宽度范围，格式: min,max（默认: 6,12）')
    parser.add_argument('--nmotifs', type=int, default=10,
                       help='最大motif数（默认: 10）')
    parser.add_argument('--evalth', type=float, default=1e-4,
                       help='E-value阈值（默认: 1e-4）')

    parser.add_argument('--minsites', type=int, default=5,
                       help='每个motif最少位点数（默认: 5）')
    parser.add_argument('--maxsites', type=int, default=100,
                       help='每个motif最多位点数（默认: 100）')

    args = parser.parse_args()

    # 解析宽度参数
    try:
        width_parts = args.width.split(',')
        width_min = int(width_parts[0])
        width_max = int(width_parts[1]) if len(width_parts) > 1 else width_min
    except:
        logger.warning(f"无法解析宽度参数: {args.width}，使用默认值6,12")
        width_min, width_max = 6, 12

    logger.info("=== MEME motif分析开始 ===")
    logger.info(f"输入文件: {args.fasta}")
    logger.info(f"输出目录: {args.output}")

    # 运行分析
    result = run_meme_analysis(
        fasta_file=args.fasta,
        output_dir=args.output,
        width_min=width_min,
        width_max=width_max,
        max_motifs=args.nmotifs,
        evalue_threshold=args.evalth
    )

    # 输出结果
    if result.get('success'):
        logger.info("=== MEME分析成功 ===")
        logger.info(f"输出目录: {result.get('output_dir')}")

        motifs = result.get('motifs', [])
        if motifs:
            logger.info(f"发现 {len(motifs)} 个显著motif:")
            for i, motif in enumerate(motifs):
                logger.info(f"  Motif {i+1}: ID={motif.get('id')}, "
                          f"宽度={motif.get('width')}nt, "
                          f"位点数={motif.get('sites')}, "
                          f"E-value={motif.get('evalue', 0):.2e}")
                if 'consensus' in motif:
                    logger.info(f"    共有序列: {motif['consensus']}")

        # 检查重要文件
        files = result.get('files', {})
        if 'meme.html' in files:
            logger.info(f"HTML报告: {files['meme.html']}")
        if 'meme.txt' in files:
            logger.info(f"文本结果: {files['meme.txt']}")

        logger.info("分析完成！")
        return 0
    else:
        logger.error("=== MEME分析失败 ===")
        logger.error(f"错误: {result.get('error', '未知错误')}")
        return 1


if __name__ == '__main__':
    sys.exit(main())