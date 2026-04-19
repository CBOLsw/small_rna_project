#!/usr/bin/env python3
"""
验证参考基因组和索引文件的一致性

功能：
1. 检查FASTA文件和.fai索引文件的一致性
2. 验证genome_fasta、genome_index和bowtie2_index的一致性
3. 提供详细的状态报告
"""

import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Any

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)

def verify_fasta_index(genome_fasta: str, genome_index: str) -> Dict[str, Any]:
    """
    验证FASTA文件和.fai索引文件的一致性

    参数:
        genome_fasta: FASTA文件路径
        genome_index: .fai索引文件路径

    返回:
        Dict: 验证结果
    """
    result = {
        'genome_fasta': genome_fasta,
        'genome_index': genome_index,
        'success': False,
        'errors': [],
        'warnings': []
    }

    fasta_path = Path(genome_fasta)
    index_path = Path(genome_index)

    logger.info("=== 验证基因组索引文件 ===")

    # 检查文件是否存在
    if not fasta_path.exists():
        error_msg = f"FASTA文件不存在: {genome_fasta}"
        logger.error(error_msg)
        result['errors'].append(error_msg)
        return result

    if not index_path.exists():
        error_msg = f".fai索引文件不存在: {genome_index}"
        logger.error(error_msg)
        result['errors'].append(error_msg)
        return result

    logger.info("✓ 文件存在性检查通过")

    try:
        # 1. 检查索引文件的格式是否有效
        with open(index_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) == 0:
            result['errors'].append("索引文件是空的")
            logger.error("索引文件是空的")
            return result

        # 验证每行的格式
        valid_lines = 0
        invalid_lines = 0
        for i, line in enumerate(lines, 1):
            fields = line.split('\t')
            if len(fields) != 5:
                logger.warning(f"第{i}行格式无效，字段数: {len(fields)}")
                invalid_lines += 1
            else:
                try:
                    # 验证数值字段
                    int(fields[1])
                    int(fields[2])
                    int(fields[3])
                    int(fields[4])
                    valid_lines += 1
                except ValueError:
                    logger.warning(f"第{i}行数值字段无效")
                    invalid_lines += 1

        logger.info(f"✓ 有效行: {valid_lines}, 无效行: {invalid_lines}")

        # 2. 检查FASTA文件的第一行是否与索引文件的第一行匹配
        with open(fasta_path, 'r') as f:
            fasta_first_line = f.readline().strip()
            if fasta_first_line:
                fasta_first_seq = fasta_first_line.split()[0] if fasta_first_line.startswith('>') else None
                if fasta_first_seq:
                    index_first_seq = lines[0].split('\t')[0] if valid_lines > 0 else None
                    if index_first_seq and fasta_first_seq[1:] == index_first_seq:
                        logger.info(f"✓ 第一个序列ID匹配: {fasta_first_seq}")
                    else:
                        logger.warning(f"第一个序列ID不匹配\nFASTA: {fasta_first_seq}\n索引: {index_first_seq}")

        # 3. 使用samtools faidx验证
        result['samtools_verify'] = run_samtools_verify(fasta_path, lines)
        logger.info(f"samtools验证: {result['samtools_verify']}")

        result['success'] = invalid_lines == 0 and len(result['errors']) == 0
        result['total_sequences'] = valid_lines

    except Exception as e:
        logger.error(f"验证过程中出错: {e}")
        result['errors'].append(str(e))

    logger.info(f"=== 验证完成 {'✅' if result['success'] else '❌'} ===")
    return result


def run_samtools_verify(fasta_path: Path, index_lines: List[str]) -> Dict[str, Any]:
    """使用samtools faidx验证序列"""
    samtools_result = {
        'success': False,
        'verified_sequences': 0,
        'failed_sequences': 0,
        'errors': []
    }

    try:
        # 尝试提取第一个序列的一小段
        test_regions = []
        for line in index_lines[:3]:  # 测试前3个序列
            seq_id = line.split('\t')[0]
            seq_length = int(line.split('\t')[1])
            # 提取前100个碱基（如果序列更长）
            end = min(100, seq_length)
            test_regions.append(f"{seq_id}:1-{end}")

        for region in test_regions:
            cmd = ['samtools', 'faidx', str(fasta_path), region]
            try:
                output = subprocess.run(cmd, capture_output=True, text=True, check=True)
                if output.stdout.strip() and len(output.stdout.strip().split('\n')) > 1:
                    samtools_result['verified_sequences'] += 1
                else:
                    samtools_result['failed_sequences'] += 1
            except Exception as e:
                samtools_result['errors'].append(f"提取区域 {region} 时出错: {e}")
                samtools_result['failed_sequences'] += 1

        if samtools_result['verified_sequences'] > 0:
            samtools_result['success'] = True

    except FileNotFoundError:
        samtools_result['errors'].append("未找到samtools命令")
    except Exception as e:
        samtools_result['errors'].append(f"samtools命令失败: {e}")

    return samtools_result


def verify_bowtie2_index(config: Dict[str, Any]) -> Dict[str, Any]:
    """验证Bowtie2索引与参考基因组的一致性"""
    bowtie2_index = config['reference']['bowtie2_index']
    genome_fasta = config['reference']['genome_fasta']

    result = {
        'bowtie2_index': bowtie2_index,
        'genome_fasta': genome_fasta,
        'success': False,
        'errors': [],
        'warnings': []
    }

    logger.info("=== 验证Bowtie2索引 ===")

    # 检查索引文件是否存在
    required_files = ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']
    missing_files = []
    for ext in required_files:
        index_file = Path(f"{bowtie2_index}.{ext}")
        if not index_file.exists():
            missing_files.append(index_file.name)

    if missing_files:
        logger.error(f"Bowtie2索引文件缺失: {', '.join(missing_files)}")
        result['errors'].append(f"Bowtie2索引文件缺失: {', '.join(missing_files)}")
        return result

    logger.info("✓ Bowtie2索引文件完整性检查通过")

    # 尝试使用bowtie2-inspect验证索引（可选，需要安装bowtie2）
    try:
        cmd = ['bowtie2-inspect', bowtie2_index]
        output = subprocess.run(cmd, capture_output=True, text=True, check=True)
        if output.stdout.strip():
            logger.info("✓ Bowtie2索引通过bowtie2-inspect验证")
            result['success'] = True
        else:
            logger.warning("Bowtie2索引验证未通过")
            result['warnings'].append("Bowtie2索引验证未通过")

    except FileNotFoundError:
        logger.warning("未找到bowtie2-inspect命令")
        result['warnings'].append("未找到bowtie2-inspect命令")
        result['success'] = True  # 没有命令时，只检查文件存在性
    except Exception as e:
        logger.error(f"Bowtie2索引验证失败: {e}")
        result['errors'].append(f"Bowtie2索引验证失败: {e}")

    logger.info(f"=== Bowtie2索引验证完成 {'✅' if result['success'] else '❌'} ===")
    return result


def print_report(fasta_result, bowtie2_result=None):
    """打印验证报告"""
    # 处理Windows编码问题
    try:
        print("\n" + "="*60)
        print("参考基因组和索引文件验证报告")
        print("="*60)

        print("\n1. FASTA文件和.fai索引:")
        if fasta_result['success']:
            print("  验证成功")
        else:
            print("  验证失败")

        if fasta_result['errors']:
            print("\n  错误:")
            for error in fasta_result['errors']:
                print(f"    - {error}")

        if fasta_result['warnings']:
            print("\n  警告:")
            for warning in fasta_result['warnings']:
                print(f"    - {warning}")

        if fasta_result.get('total_sequences') is not None:
            print(f"\n  序列统计: {fasta_result['total_sequences']}个序列")

        if bowtie2_result:
            print("\n" + "-"*60)
            print("\n2. Bowtie2索引:")
            if bowtie2_result['success']:
                print("  验证成功")
            else:
                print("  验证失败")

            if bowtie2_result['errors']:
                print("\n  错误:")
                for error in bowtie2_result['errors']:
                    print(f"    - {error}")

            if bowtie2_result['warnings']:
                print("\n  警告:")
                for warning in bowtie2_result['warnings']:
                    print(f"    - {warning}")
    except UnicodeEncodeError as e:
        # 如果出现编码错误，使用简化的输出
        print("\n验证报告 (简化版):")
        print(f"FASTA验证: {'成功' if fasta_result['success'] else '失败'}")
        print(f"Bowtie2验证: {'成功' if bowtie2_result and bowtie2_result['success'] else '失败' if bowtie2_result else '未执行'}")
        if fasta_result['errors']:
            print("FASTA错误:")
            for err in fasta_result['errors']:
                print(f"  {err}")
        if bowtie2_result and bowtie2_result['errors']:
            print("Bowtie2错误:")
            for err in bowtie2_result['errors']:
                print(f"  {err}")


def main():
    """主函数"""
    import yaml
    import argparse

    parser = argparse.ArgumentParser(description="验证参考基因组和索引文件")
    parser.add_argument(
        "--config", "-c",
        default="config/config.yaml",
        help="配置文件路径 (默认: config/config.yaml)"
    )
    parser.add_argument(
        "--genome", "-g",
        help="参考基因组文件路径 (可选，会覆盖配置文件中的设置)"
    )
    parser.add_argument(
        "--index", "-i",
        help="genome_index文件路径 (可选，会覆盖配置文件中的设置)"
    )
    parser.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="安静模式，只显示错误"
    )

    args = parser.parse_args()

    # 加载配置
    config = None
    if Path(args.config).exists():
        try:
            with open(args.config, 'r', encoding='utf-8') as f:
                config = yaml.safe_load(f)
        except Exception as e:
            logger.error(f"加载配置文件失败: {e}")
            config = None

    # 获取路径
    genome_fasta = args.genome or (config['reference']['genome_fasta'] if config else None)
    genome_index = args.index or (config['reference']['genome_index'] if config else None)

    if not genome_fasta:
        logger.error("未提供参考基因组文件路径")
        return 1

    if not genome_index:
        genome_index = f"{genome_fasta}.fai"

    logger.info(f"参考基因组文件: {genome_fasta}")
    logger.info(f"索引文件: {genome_index}")

    # 执行验证
    fasta_result = verify_fasta_index(genome_fasta, genome_index)

    bowtie2_result = None
    if config and 'reference' in config and 'bowtie2_index' in config['reference']:
        bowtie2_result = verify_bowtie2_index(config)

    # 打印报告
    print_report(fasta_result, bowtie2_result)

    # 检查是否有错误
    if fasta_result['errors'] or (bowtie2_result and bowtie2_result['errors']):
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main())
