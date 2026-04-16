#!/usr/bin/env python3
"""
Bowtie2参考基因组索引构建脚本

功能：
1. 构建Bowtie2参考基因组索引
2. 支持small RNA测序的优化参数
3. 检查现有索引，避免重复构建
4. 生成索引构建报告

使用方法：
    python build_bowtie2_index.py --genome <参考基因组文件> --output <输出目录>
"""

import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple
import logging
import json

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class Bowtie2IndexBuilder:
    """Bowtie2索引构建器类"""

    # Bowtie2索引文件扩展名
    INDEX_EXTENSIONS = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']

    def __init__(self, bowtie2_build_path: str = "bowtie2-build"):
        """
        初始化Bowtie2索引构建器

        参数:
            bowtie2_build_path: bowtie2-build可执行文件路径
        """
        self.bowtie2_build_path = bowtie2_build_path
        self.results = {}

    def check_bowtie2(self) -> bool:
        """检查Bowtie2是否可用"""
        try:
            result = subprocess.run(
                [self.bowtie2_build_path, "--version"],
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                version_line = result.stdout.strip().split('\n')[0]
                logger.info(f"Bowtie2版本: {version_line}")
                return True
            else:
                logger.error(f"Bowtie2检查失败: {result.stderr}")
                return False
        except FileNotFoundError:
            logger.error(f"未找到Bowtie2: {self.bowtie2_build_path}")
            return False

    def check_index_exists(self, index_prefix: str) -> bool:
        """
        检查索引文件是否已存在

        参数:
            index_prefix: 索引文件前缀

        返回:
            bool: 索引是否完整存在
        """
        index_prefix = Path(index_prefix)
        missing_files = []

        for ext in self.INDEX_EXTENSIONS:
            index_file = index_prefix.with_suffix(ext)
            if not index_file.exists():
                missing_files.append(index_file.name)

        if missing_files:
            logger.info(f"索引文件不完整，缺失 {len(missing_files)} 个文件")
            return False
        else:
            logger.info(f"索引文件已存在: {index_prefix}")
            return True

    def build_index(self, genome_file: str, output_prefix: str,
                    threads: int = 4, small_rna_mode: bool = True) -> Dict[str, any]:
        """
        构建Bowtie2索引

        参数:
            genome_file: 参考基因组文件路径
            output_prefix: 输出索引文件前缀
            threads: 线程数
            small_rna_mode: 是否使用small RNA优化参数

        返回:
            Dict: 构建结果
        """
        genome_path = Path(genome_file)
        output_path = Path(output_prefix)

        # 检查输入文件
        if not genome_path.exists():
            logger.error(f"参考基因组文件不存在: {genome_path}")
            return {'success': False, 'error': '输入文件不存在'}

        # 检查输出目录
        output_dir = output_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)

        # 检查是否已存在索引
        if self.check_index_exists(output_prefix):
            logger.info(f"索引已存在，跳过构建: {output_prefix}")
            return {
                'success': True,
                'skipped': True,
                'index_prefix': str(output_prefix),
                'message': '索引已存在'
            }

        # 准备构建命令
        cmd = self._build_index_command(
            genome_file, output_prefix, threads, small_rna_mode
        )

        logger.info(f"开始构建Bowtie2索引: {genome_path.name}")
        logger.info(f"输出前缀: {output_prefix}")
        logger.info(f"线程数: {threads}")
        logger.info(f"small RNA模式: {small_rna_mode}")

        # 执行构建命令
        success, stats = self._run_bowtie2_build(cmd, genome_path, output_prefix)

        result = {
            'success': success,
            'skipped': False,
            'genome_file': str(genome_path),
            'index_prefix': str(output_prefix),
            'threads': threads,
            'small_rna_mode': small_rna_mode,
            'stats': stats
        }

        self.results[str(output_prefix)] = result
        return result

    def _build_index_command(self, genome_file: str, output_prefix: str,
                             threads: int, small_rna_mode: bool) -> List[str]:
        """构建Bowtie2索引命令"""
        cmd = [
            self.bowtie2_build_path,
            "--threads", str(threads),
        ]

        # small RNA优化参数
        if small_rna_mode:
            # 对于small RNA，使用较小的seed length和更宽松的参数
            cmd.extend([
                "--seed", "16",  # 较小的seed length适应短序列
                "--bmax", "100",  # 减少内存使用
                "--dcv", "1024",  # 减少内存使用
                "--nodc",  # 禁用深度压缩
            ])

        cmd.extend([
            genome_file,
            str(output_prefix)
        ])

        return cmd

    def _run_bowtie2_build(self, cmd: List[str], genome_path: Path,
                           output_prefix: str) -> Tuple[bool, Dict[str, any]]:
        """运行Bowtie2索引构建命令"""
        logger.info(f"运行命令: {' '.join(cmd)}")

        try:
            # 记录开始时间
            import time
            start_time = time.time()

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            end_time = time.time()
            elapsed_time = end_time - start_time

            if result.returncode == 0:
                logger.info(f"索引构建成功: {output_prefix}")
                logger.info(f"构建时间: {elapsed_time:.1f} 秒")

                # 检查生成的索引文件
                index_files = []
                for ext in self.INDEX_EXTENSIONS:
                    index_file = Path(output_prefix).with_suffix(ext)
                    if index_file.exists():
                        index_files.append(str(index_file))
                        file_size = index_file.stat().st_size / (1024**2)  # MB
                        logger.debug(f"  索引文件: {index_file.name} ({file_size:.1f} MB)")

                stats = {
                    'elapsed_time': elapsed_time,
                    'index_files': index_files,
                    'num_files': len(index_files),
                    'stdout': result.stdout.strip(),
                    'stderr': result.stderr.strip()
                }

                return True, stats
            else:
                logger.error(f"索引构建失败: {output_prefix}")
                logger.error(f"错误输出: {result.stderr}")
                logger.error(f"返回码: {result.returncode}")

                stats = {
                    'elapsed_time': elapsed_time,
                    'stdout': result.stdout.strip(),
                    'stderr': result.stderr.strip(),
                    'returncode': result.returncode
                }

                return False, stats

        except Exception as e:
            logger.error(f"运行Bowtie2构建时出错: {e}")
            return False, {'error': str(e)}

    def generate_report(self, output_dir: str) -> Optional[str]:
        """
        生成索引构建报告

        参数:
            output_dir: 输出目录

        返回:
            str: 报告文件路径
        """
        if not self.results:
            logger.warning("无构建结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建汇总数据
        summary_data = []
        for index_prefix, result in self.results.items():
            stats = result.get('stats', {})
            summary_data.append({
                'index_prefix': index_prefix,
                'genome_file': result.get('genome_file'),
                'success': result.get('success', False),
                'skipped': result.get('skipped', False),
                'threads': result.get('threads', 0),
                'small_rna_mode': result.get('small_rna_mode', False),
                'elapsed_time': stats.get('elapsed_time', 0),
                'num_index_files': stats.get('num_files', 0),
                'error': stats.get('error', '') if not result.get('success') else ''
            })

        # 保存为JSON
        json_file = output_dir / "bowtie2_index_build.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        # 生成文本报告
        report_file = output_dir / "bowtie2_index_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== Bowtie2索引构建报告 ===\n\n")

            total = len(self.results)
            successful = len([r for r in self.results.values() if r.get('success')])
            skipped = len([r for r in self.results.values() if r.get('skipped')])

            f.write(f"总构建任务: {total}\n")
            f.write(f"成功构建: {successful}\n")
            f.write(f"跳过构建: {skipped}\n")
            f.write(f"失败构建: {total - successful - skipped}\n\n")

            for index_prefix, result in self.results.items():
                f.write(f"索引: {index_prefix}\n")
                f.write(f"  状态: {'成功' if result.get('success') else '失败'}\n")
                if result.get('skipped'):
                    f.write(f"  操作: 跳过（索引已存在）\n")
                else:
                    stats = result.get('stats', {})
                    f.write(f"  构建时间: {stats.get('elapsed_time', 0):.1f} 秒\n")
                    f.write(f"  索引文件数: {stats.get('num_files', 0)}\n")
                f.write("\n")

            f.write(f"报告生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

        logger.info(f"索引构建报告已生成: {report_file}")
        return str(report_file)


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="Bowtie2参考基因组索引构建脚本"
    )

    parser.add_argument(
        "--genome", "-g",
        required=True,
        help="参考基因组文件路径 (FASTA格式)"
    )

    parser.add_argument(
        "--output", "-o",
        required=True,
        help="输出索引文件前缀 (如: references/bowtie2_index/hg38)"
    )

    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=4,
        help="线程数 (默认: 4)"
    )

    parser.add_argument(
        "--no-small-rna",
        action="store_true",
        help="禁用small RNA优化参数"
    )

    parser.add_argument(
        "--bowtie2-path",
        default="bowtie2-build",
        help="bowtie2-build可执行文件路径 (默认: bowtie2-build)"
    )

    parser.add_argument(
        "--report-dir",
        default="results/alignment/index_build",
        help="报告输出目录 (默认: results/alignment/index_build)"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 初始化构建器
    builder = Bowtie2IndexBuilder(bowtie2_build_path=args.bowtie2_path)

    # 检查Bowtie2
    if not builder.check_bowtie2():
        logger.error("Bowtie2检查失败，请确保Bowtie2已安装")
        sys.exit(1)

    # 构建索引
    small_rna_mode = not args.no_small_rna
    result = builder.build_index(
        genome_file=args.genome,
        output_prefix=args.output,
        threads=args.threads,
        small_rna_mode=small_rna_mode
    )

    if not result.get('success'):
        logger.error("索引构建失败")
        sys.exit(1)

    # 生成报告
    report_file = builder.generate_report(args.report_dir)
    if report_file:
        logger.info(f"索引构建完成，报告文件: {report_file}")
    else:
        logger.warning("未能生成构建报告")

    logger.info("Bowtie2索引构建流程完成")


if __name__ == "__main__":
    main()