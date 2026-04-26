#!/usr/bin/env python3
"""
TomTom motif比较脚本
用于与已知motif数据库（JASPAR、CIS-BP等）比对，识别相似motif

功能：
1. 运行TomTom motif比较
2. 支持多种motif数据库格式
3. 统计显著性评估（E-value阈值）
4. 结果过滤和导出

使用方法：
    python tomtom_comparison.py --motif <输入motif文件> --output <输出目录> [--db <数据库路径>]
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


def check_tomtom_installed(tomtom_path: str = "tomtom") -> bool:
    """检查TomTom是否安装"""
    try:
        result = subprocess.run(
            [tomtom_path, '-version'],
            capture_output=True,
            text=True,
            check=False
        )
        return result.returncode == 0
    except FileNotFoundError:
        return False


def ensure_motif_database(database: Optional[str]) -> Optional[str]:
    """
    确保motif数据库可用，必要时自动下载JASPAR
    """
    if database and Path(database).exists():
        return database

    # 检查项目本地路径
    local_db = Path("references/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant.meme")
    if local_db.exists():
        logger.info(f"使用本地数据库: {local_db}")
        return str(local_db)

    old_db = Path("references/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme")
    if old_db.exists():
        logger.info(f"使用本地数据库: {old_db}")
        return str(old_db)

    # 检查MEME默认安装路径
    for prefix in ["/usr/local", "/usr", os.path.expanduser("~/miniconda3/envs/small_rna_analysis")]:
        for db_path in [
            f"{prefix}/share/meme/db/motif_databases/JASPAR/JASPAR2024_CORE_vertebrates_non-redundant.meme",
            f"{prefix}/share/meme/db/motif_databases/JASPAR/JASPAR2022_CORE_vertebrates_non-redundant.meme",
        ]:
            if Path(db_path).exists():
                logger.info(f"使用系统数据库: {db_path}")
                return db_path

    # 自动下载JASPAR数据库
    logger.info("未找到JASPAR数据库，开始自动下载...")
    db_dir = Path("references/motif_databases/JASPAR")
    db_dir.mkdir(parents=True, exist_ok=True)

    import urllib.request
    import tarfile

    url = "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.25.tgz"
    tgz_path = Path("references/motif_databases/motif_databases.tgz")

    try:
        logger.info(f"下载motif数据库: {url}")
        urllib.request.urlretrieve(url, tgz_path)
        logger.info("下载完成，解压中...")

        with tarfile.open(tgz_path, "r:gz") as tar:
            tar.extractall(path="references/motif_databases/")
        tgz_path.unlink()

        # 查找JASPAR文件
        for p in db_dir.parent.rglob("JASPAR*_CORE_vertebrates_non-redundant.meme"):
            logger.info(f"数据库就绪: {p}")
            return str(p)

        logger.error("下载完成但未找到JASPAR数据库文件")
    except Exception as e:
        logger.error(f"下载数据库失败: {e}")
        if tgz_path.exists():
            tgz_path.unlink()

    return database


def run_tomtom_comparison(motif_file: str, output_dir: str,
                         database: str = None,
                         evalue_threshold: float = 0.05,
                         min_overlap: int = 5) -> Dict[str, Any]:
    """
    运行TomTom motif比较分析

    参数:
        motif_file: 输入motif文件（MEME格式）
        output_dir: 输出目录
        database: motif数据库路径（默认使用JASPAR脊椎动物核心数据库）
        evalue_threshold: E-value阈值
        min_overlap: 最小重叠长度

    返回:
        Dict: 分析结果
    """
    # 检查输入文件
    if not Path(motif_file).exists():
        logger.error(f"输入文件不存在: {motif_file}")
        return {'success': False, 'error': '输入文件不存在'}

    # 检查TomTom是否安装
    if not check_tomtom_installed():
        logger.error("TomTom未安装或不在PATH中")
        return {'success': False, 'error': 'TomTom未安装'}

    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 确保数据库可用（本地查找或自动下载）
    database = ensure_motif_database(database)
    if database is None:
        logger.error("未找到motif数据库且自动下载失败")
        return {'success': False, 'error': 'motif数据库未找到'}

    logger.info(f"开始TomTom motif比较分析")
    logger.info(f"输入motif文件: {motif_file}")
    logger.info(f"数据库: {database}")
    logger.info(f"E-value阈值: {evalue_threshold}")
    logger.info(f"最小重叠长度: {min_overlap}")

    try:
        # 构建TomTom命令
        cmd = [
            'tomtom',
            '-o', str(output_path),
            '-evalue', str(evalue_threshold),
            '-min-overlap', str(min_overlap),
            '-text',  # 输出文本格式便于解析
            motif_file,
            database
        ]

        logger.info(f"运行TomTom命令: {' '.join(cmd)}")

        # 执行TomTom命令
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )

        logger.info("TomTom分析完成")

        # 检查输出文件
        result_files = {}
        expected_files = ['tomtom.html', 'tomtom.xml', 'tomtom.txt']

        for file_name in expected_files:
            file_path = output_path / file_name
            if file_path.exists():
                result_files[file_name] = str(file_path)

        # 解析主要结果
        txt_file = output_path / 'tomtom.txt'
        comparisons = parse_tomtom_text_result(txt_file) if txt_file.exists() else []

        result = {
            'success': True,
            'output_dir': str(output_path),
            'files': result_files,
            'comparisons_found': len(comparisons),
            'comparisons': comparisons,
            'parameters': {
                'database': database,
                'evalue_threshold': evalue_threshold,
                'min_overlap': min_overlap
            }
        }

        # 保存结果摘要
        save_tomtom_summary(result, output_path)

        logger.info(f"发现 {len(comparisons)} 个显著motif比对")
        return result

    except subprocess.CalledProcessError as e:
        logger.error(f"TomTom执行失败: {e}")
        logger.error(f"STDERR: {e.stderr[:500]}")
        return {'success': False, 'error': f'TomTom执行失败: {e.stderr[:200]}'}
    except Exception as e:
        logger.error(f"运行TomTom时出错: {e}")
        return {'success': False, 'error': str(e)}


def parse_tomtom_text_result(txt_file: Path) -> List[Dict[str, Any]]:
    """解析TomTom文本格式结果"""
    comparisons = []

    try:
        with open(txt_file, 'r') as f:
            lines = f.readlines()

        # TomTom文本格式：制表符分隔
        # 跳过标题行
        header_found = False
        for line in lines:
            line = line.strip()

            if not line:
                continue

            # 查找标题行
            if line.startswith('#Query_ID'):
                header_found = True
                headers = line[1:].split('\t')  # 去掉#号
                continue

            if header_found and line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    comparison = {}
                    for i, header in enumerate(headers):
                        if i < len(parts):
                            value = parts[i].strip()

                            # 转换数值类型
                            if header in ['E-value', 'q-value', 'Overlap', 'Offset']:
                                try:
                                    if value == 'N/A' or value == '':
                                        value = None
                                    else:
                                        value = float(value)
                                except:
                                    value = None
                            elif header in ['Orientation']:
                                value = value

                            comparison[header] = value

                    comparisons.append(comparison)

    except Exception as e:
        logger.warning(f"解析TomTom结果时出错: {e}")

    return comparisons


def save_tomtom_summary(result: Dict[str, Any], output_dir: Path):
    """保存TomTom结果摘要"""
    # JSON摘要
    summary_file = output_dir / "tomtom_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(result, f, indent=2)

    # CSV表格
    if result.get('comparisons'):
        df = pd.DataFrame(result['comparisons'])
        csv_file = output_dir / "tomtom_comparisons.csv"
        df.to_csv(csv_file, index=False)

        # 创建简化的Top N结果
        top_n = 20
        if len(df) > 0:
            # 按E-value排序
            df_sorted = df.copy()
            if 'E-value' in df_sorted.columns:
                df_sorted = df_sorted.sort_values('E-value')
                top_file = output_dir / f"tomtom_top{top_n}.csv"
                df_sorted.head(top_n).to_csv(top_file, index=False)

    logger.info(f"TomTom结果摘要已保存: {summary_file}")


def prepare_motif_input(meme_output_dir: str, output_file: str = None) -> str:
    """
    从MEME输出目录准备TomTom输入文件

    参数:
        meme_output_dir: MEME输出目录
        output_file: 输出文件路径（如未指定则自动生成）

    返回:
        str: 准备好的motif文件路径
    """
    meme_dir = Path(meme_output_dir)

    # 检查可能的motif文件
    motif_files = [
        meme_dir / 'meme.xml',      # XML格式
        meme_dir / 'meme.txt',      # 文本格式
        meme_dir / 'meme.html',     # HTML格式（通常不是最佳）
    ]

    for file_path in motif_files:
        if file_path.exists():
            logger.info(f"找到motif文件: {file_path}")

            # 如果未指定输出文件，使用默认名称
            if output_file is None:
                output_file = str(meme_dir / 'motifs_for_tomtom.meme')

            # 如果是XML或文本格式，可能需要转换
            # 这里简化处理：直接使用原始文件
            # 实际应用中可能需要格式转换
            if file_path.suffix == '.xml':
                # 可以在这里添加XML到MEME格式的转换
                pass

            return str(file_path)

    logger.error(f"在目录 {meme_dir} 中未找到motif文件")
    return None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="TomTom motif比较脚本 - 用于与已知数据库比对motif",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 基本用法：使用默认数据库
  python tomtom_comparison.py --motif meme_output/meme.xml --output tomtom_results/

  # 指定数据库
  python tomtom_comparison.py --motif meme.xml --output results/ --db /path/to/database.meme

  # 从MEME输出目录准备输入
  python tomtom_comparison.py --meme-output meme_results/ --output tomtom_results/
        """
    )

    # 输入参数组
    input_group = parser.add_argument_group('输入参数')
    input_mutex = input_group.add_mutually_exclusive_group(required=True)
    input_mutex.add_argument('--motif', help='输入motif文件（MEME格式）')
    input_mutex.add_argument('--meme-output', help='MEME输出目录（自动提取motif文件）')

    # 输出参数
    parser.add_argument('--output', required=True, help='输出目录')

    # TomTom参数
    parser.add_argument('--db', help='motif数据库路径（默认：JASPAR脊椎动物核心数据库）')
    parser.add_argument('--evalth', type=float, default=0.05,
                       help='E-value阈值（默认: 0.05）')
    parser.add_argument('--min-overlap', type=int, default=5,
                       help='最小重叠长度（默认: 5）')

    # 其他参数
    parser.add_argument('--top-n', type=int, default=20,
                       help='保存前N个最显著结果（默认: 20）')

    args = parser.parse_args()

    # 确定输入文件
    motif_file = args.motif
    if args.meme_output:
        logger.info(f"从MEME输出目录准备输入: {args.meme_output}")
        motif_file = prepare_motif_input(args.meme_output)
        if motif_file is None:
            logger.error("无法从MEME输出目录准备输入文件")
            return 1

    logger.info("=== TomTom motif比较分析开始 ===")
    logger.info(f"输入文件: {motif_file}")
    logger.info(f"输出目录: {args.output}")

    # 运行分析
    result = run_tomtom_comparison(
        motif_file=motif_file,
        output_dir=args.output,
        database=args.db,
        evalue_threshold=args.evalth,
        min_overlap=args.min_overlap
    )

    # 输出结果
    if result.get('success'):
        logger.info("=== TomTom分析成功 ===")
        logger.info(f"输出目录: {result.get('output_dir')}")

        comparisons = result.get('comparisons', [])
        if comparisons:
            logger.info(f"发现 {len(comparisons)} 个显著motif比对")

            # 显示前5个最显著结果
            sorted_comparisons = sorted(
                comparisons,
                key=lambda x: x.get('E-value', float('inf')) if x.get('E-value') is not None else float('inf')
            )

            top_n = min(5, len(sorted_comparisons))
            logger.info(f"前 {top_n} 个最显著结果:")
            for i, comp in enumerate(sorted_comparisons[:top_n]):
                query_id = comp.get('Query_ID', 'N/A')
                target_id = comp.get('Target_ID', 'N/A')
                evalue = comp.get('E-value', 'N/A')
                qvalue = comp.get('q-value', 'N/A')

                evalue_str = f"{evalue:.2e}" if isinstance(evalue, float) else evalue
                qvalue_str = f"{qvalue:.2e}" if isinstance(qvalue, float) else qvalue

                logger.info(f"  {i+1}. {query_id} -> {target_id} "
                          f"(E-value: {evalue_str}, q-value: {qvalue_str})")

        # 检查重要文件
        files = result.get('files', {})
        if 'tomtom.html' in files:
            logger.info(f"HTML报告: {files['tomtom.html']}")
        if 'tomtom.txt' in files:
            logger.info(f"文本结果: {files['tomtom.txt']}")
        if 'tomtom_comparisons.csv' in files:
            logger.info(f"CSV结果: {files.get('tomtom_comparisons.csv', 'tomtom_comparisons.csv')}")

        logger.info("分析完成！")
        return 0
    else:
        logger.error("=== TomTom分析失败 ===")
        logger.error(f"错误: {result.get('error', '未知错误')}")
        return 1


if __name__ == '__main__':
    sys.exit(main())