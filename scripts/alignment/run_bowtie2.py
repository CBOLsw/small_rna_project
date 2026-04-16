#!/usr/bin/env python3
"""
Bowtie2序列比对脚本 - 针对small RNA测序数据优化

功能：
1. 将fastq reads比对到参考基因组
2. 优化small RNA短序列比对参数
3. 支持单端和双端测序数据
4. 生成比对统计和质量报告

使用方法：
    python run_bowtie2.py --input <fastq文件或样本信息> --index <Bowtie2索引前缀> --output <输出目录>
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


class Bowtie2Aligner:
    """Bowtie2比对器类"""

    def __init__(self, bowtie2_path: str = "bowtie2"):
        """
        初始化Bowtie2比对器

        参数:
            bowtie2_path: bowtie2可执行文件路径
        """
        self.bowtie2_path = bowtie2_path
        self.results = {}

    def check_bowtie2(self) -> bool:
        """检查Bowtie2是否可用"""
        try:
            result = subprocess.run(
                [self.bowtie2_path, "--version"],
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
            logger.error(f"未找到Bowtie2: {self.bowtie2_path}")
            return False

    def align_single_end(self, fastq_file: str, index_prefix: str,
                         output_dir: str, sample_name: str,
                         config: Dict[str, Any]) -> Dict[str, Any]:
        """
        单端测序数据比对

        参数:
            fastq_file: fastq文件路径
            index_prefix: Bowtie2索引前缀
            output_dir: 输出目录
            sample_name: 样本名称
            config: 配置参数

        返回:
            Dict: 比对结果
        """
        logger.info(f"单端比对: {sample_name}")

        # 准备输出文件
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        sam_file = output_dir / f"{sample_name}.sam"
        bam_file = output_dir / f"{sample_name}.bam"
        sorted_bam = output_dir / f"{sample_name}_sorted.bam"
        log_file = output_dir / f"{sample_name}_bowtie2.log"
        stats_file = output_dir / f"{sample_name}_alignment_stats.txt"

        # 构建Bowtie2命令
        cmd = self._build_single_end_command(
            fastq_file, index_prefix, sam_file, config
        )

        # 运行比对
        success, align_stats = self._run_bowtie2(cmd, log_file, sample_name)

        if success:
            # 转换SAM到BAM并排序
            threads = config.get('threads', 4)
            bam_success = self._sam_to_bam(
                str(sam_file), str(sorted_bam), sample_name,
                threads=threads, keep_sam=config.get('keep_sam', False)
            )
            if bam_success:
                # 索引BAM文件
                self._index_bam(str(sorted_bam), sample_name)

            # 解析比对统计
            stats = self._parse_alignment_stats(log_file, align_stats)
            stats.update({
                'bam_file': str(sorted_bam) if bam_success else None,
                'bam_indexed': bam_success
            })

            # 保存统计信息
            self._save_alignment_stats(stats, stats_file)
        else:
            stats = {'success': False}

        result = {
            'sample': sample_name,
            'type': 'single_end',
            'fastq_file': fastq_file,
            'index_prefix': index_prefix,
            'sam_file': str(sam_file),
            'bam_file': str(sorted_bam) if success else None,
            'log_file': str(log_file),
            'stats_file': str(stats_file) if success else None,
            'success': success,
            'stats': stats
        }

        self.results[sample_name] = result
        return result

    def align_paired_end(self, fastq_r1: str, fastq_r2: str, index_prefix: str,
                         output_dir: str, sample_name: str,
                         config: Dict[str, Any]) -> Dict[str, Any]:
        """
        双端测序数据比对

        参数:
            fastq_r1: R1 fastq文件路径
            fastq_r2: R2 fastq文件路径
            index_prefix: Bowtie2索引前缀
            output_dir: 输出目录
            sample_name: 样本名称
            config: 配置参数

        返回:
            Dict: 比对结果
        """
        logger.info(f"双端比对: {sample_name}")

        # 准备输出文件
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        sam_file = output_dir / f"{sample_name}.sam"
        bam_file = output_dir / f"{sample_name}.bam"
        sorted_bam = output_dir / f"{sample_name}_sorted.bam"
        log_file = output_dir / f"{sample_name}_bowtie2.log"
        stats_file = output_dir / f"{sample_name}_alignment_stats.txt"

        # 构建Bowtie2命令
        cmd = self._build_paired_end_command(
            fastq_r1, fastq_r2, index_prefix, sam_file, config
        )

        # 运行比对
        success, align_stats = self._run_bowtie2(cmd, log_file, sample_name)

        if success:
            # 转换SAM到BAM并排序
            threads = config.get('threads', 4)
            bam_success = self._sam_to_bam(
                str(sam_file), str(sorted_bam), sample_name,
                threads=threads, keep_sam=config.get('keep_sam', False)
            )
            if bam_success:
                # 索引BAM文件
                self._index_bam(str(sorted_bam), sample_name)

            # 解析比对统计
            stats = self._parse_alignment_stats(log_file, align_stats)
            stats.update({
                'bam_file': str(sorted_bam) if bam_success else None,
                'bam_indexed': bam_success
            })

            # 保存统计信息
            self._save_alignment_stats(stats, stats_file)
        else:
            stats = {'success': False}

        result = {
            'sample': sample_name,
            'type': 'paired_end',
            'fastq_r1': fastq_r1,
            'fastq_r2': fastq_r2,
            'index_prefix': index_prefix,
            'sam_file': str(sam_file),
            'bam_file': str(sorted_bam) if success else None,
            'log_file': str(log_file),
            'stats_file': str(stats_file) if success else None,
            'success': success,
            'stats': stats
        }

        self.results[sample_name] = result
        return result

    def _build_single_end_command(self, fastq_file: str, index_prefix: str,
                                  sam_file: Path, config: Dict[str, Any]) -> List[str]:
        """构建单端Bowtie2命令"""
        cmd = [
            self.bowtie2_path,
            "-x", index_prefix,
            "-U", fastq_file,
            "-S", str(sam_file),
            "--threads", str(config.get('threads', 4)),
            "--quiet",  # 减少输出
        ]

        # 添加small RNA优化参数
        cmd.extend(self._get_small_rna_params(config))

        # 添加其他参数
        cmd.extend(self._get_common_params(config))

        return cmd

    def _build_paired_end_command(self, fastq_r1: str, fastq_r2: str,
                                  index_prefix: str, sam_file: Path,
                                  config: Dict[str, Any]) -> List[str]:
        """构建双端Bowtie2命令"""
        cmd = [
            self.bowtie2_path,
            "-x", index_prefix,
            "-1", fastq_r1,
            "-2", fastq_r2,
            "-S", str(sam_file),
            "--threads", str(config.get('threads', 4)),
            "--quiet",
        ]

        # 添加small RNA优化参数
        cmd.extend(self._get_small_rna_params(config))

        # 添加双端特定参数
        if config.get('pe_max_insert'):
            cmd.extend(["--maxins", str(config['pe_max_insert'])])

        # 添加其他参数
        cmd.extend(self._get_common_params(config))

        return cmd

    def _get_small_rna_params(self, config: Dict[str, Any]) -> List[str]:
        """获取small RNA优化参数"""
        params = []

        if config.get('small_rna_mode', True):
            # small RNA特定优化
            params.extend([
                "-L", "16",  # seed length (small RNA通常18-35nt)
                "-N", "1",   # 允许1个错配
                "-i", "L,0,0.50",  # 间隔函数
                "--score-min", "L,0,-0.6",  # 最小得分
                "--no-unal",  # 不输出未比对序列
                "--no-head",  # 不在SAM中输出header
                "--no-sq",    # 不在SAM中输出参考序列信息
            ])

            # 如果指定了长度范围
            if config.get('min_len') and config.get('max_len'):
                # 注意：Bowtie2没有直接的长度过滤，需要在比对前或后处理
                pass

        return params

    def _get_common_params(self, config: Dict[str, Any]) -> List[str]:
        """获取通用参数"""
        params = []

        # 质量参数
        if config.get('phred'):
            params.extend(["--phred", str(config['phred'])])

        # 比对参数
        if config.get('sensitive'):
            params.append("--sensitive")

        if config.get('very_sensitive'):
            params.append("--very-sensitive")

        # 报告参数
        if config.get('report_all'):
            params.append("-a")

        if config.get('report_k'):
            params.extend(["-k", str(config['report_k'])])

        return params

    def _run_bowtie2(self, cmd: List[str], log_file: Path,
                     sample_name: str) -> Tuple[bool, str]:
        """运行Bowtie2命令"""
        logger.info(f"运行Bowtie2: {sample_name}")
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
                logger.info(f"Bowtie2比对成功: {sample_name}")
                # 读取日志文件获取统计信息
                with open(log_file, 'r') as f:
                    log_content = f.read()
                return True, log_content
            else:
                logger.error(f"Bowtie2比对失败: {sample_name}")
                logger.error(f"返回码: {result.returncode}")
                return False, ""

        except Exception as e:
            logger.error(f"运行Bowtie2时出错: {e}")
            return False, ""

    def _sam_to_bam(self, sam_file: str, bam_file: str, sample_name: str,
                    threads: int = 4, keep_sam: bool = False) -> bool:
        """
        转换SAM到BAM并排序

        参数:
            sam_file: 输入SAM文件路径
            bam_file: 输出BAM文件路径
            sample_name: 样本名称
            threads: 排序使用的线程数
            keep_sam: 是否保留中间SAM文件

        返回:
            bool: 转换是否成功
        """
        try:
            # 检查samtools是否可用
            samtools_check = subprocess.run(
                ["samtools", "--version"],
                capture_output=True,
                text=True,
                check=False
            )
            if samtools_check.returncode != 0:
                logger.warning("samtools不可用，跳过SAM到BAM转换")
                return False

            # 检查SAM文件是否存在
            sam_path = Path(sam_file)
            if not sam_path.exists():
                logger.error(f"SAM文件不存在: {sam_file}")
                return False

            # 检查SAM文件大小
            sam_size = sam_path.stat().st_size
            if sam_size == 0:
                logger.warning(f"SAM文件为空: {sam_file}")
                if not keep_sam:
                    sam_path.unlink(missing_ok=True)
                return False

            logger.info(f"转换SAM到BAM: {sample_name} (文件大小: {sam_size/1024/1024:.2f} MB)")

            # 步骤1: SAM到BAM转换
            logger.debug(f"步骤1: SAM到BAM转换")
            view_cmd = ["samtools", "view", "-bS", sam_file]

            # 步骤2: 排序BAM
            logger.debug(f"步骤2: 排序BAM (线程数: {threads})")
            sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam_file]

            # 使用管道连接两个步骤
            try:
                view_process = subprocess.Popen(view_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                sort_process = subprocess.Popen(sort_cmd, stdin=view_process.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                view_process.stdout.close()
                sort_stdout, sort_stderr = sort_process.communicate()

                # 检查错误
                if view_process.returncode != 0:
                    logger.error(f"SAM到BAM转换步骤失败，返回码: {view_process.returncode}")
                    if view_process.stderr:
                        logger.error(f"错误信息: {view_process.stderr.read().decode() if hasattr(view_process.stderr, 'read') else view_process.stderr}")
                    return False

                if sort_process.returncode != 0:
                    logger.error(f"BAM排序步骤失败，返回码: {sort_process.returncode}")
                    if sort_stderr:
                        logger.error(f"错误信息: {sort_stderr.decode()}")
                    return False

            except Exception as pipe_error:
                logger.error(f"管道执行失败: {pipe_error}")
                # 回退到分步执行
                return self._sam_to_bam_stepwise(sam_file, bam_file, sample_name, threads, keep_sam)

            # 验证输出文件
            bam_path = Path(bam_file)
            if bam_path.exists() and bam_path.stat().st_size > 0:
                logger.info(f"BAM文件已创建: {bam_file} (大小: {bam_path.stat().st_size/1024/1024:.2f} MB)")

                # 检查BAM文件是否有效
                check_cmd = ["samtools", "quickcheck", bam_file]
                check_result = subprocess.run(check_cmd, capture_output=True, text=True)
                if check_result.returncode == 0:
                    logger.info(f"BAM文件验证通过: {bam_file}")
                else:
                    logger.warning(f"BAM文件验证失败: {check_result.stderr}")

                # 删除临时SAM文件（如果不需要保留）
                if not keep_sam:
                    try:
                        sam_path.unlink(missing_ok=True)
                        logger.debug(f"临时SAM文件已删除: {sam_file}")
                    except Exception as del_error:
                        logger.warning(f"删除SAM文件失败: {del_error}")

                return True
            else:
                logger.error(f"BAM文件创建失败或为空: {bam_file}")
                return False

        except Exception as e:
            logger.error(f"SAM到BAM转换时出错: {e}")
            return False

    def _sam_to_bam_stepwise(self, sam_file: str, bam_file: str, sample_name: str,
                           threads: int = 4, keep_sam: bool = False) -> bool:
        """分步执行SAM到BAM转换（管道失败时的回退方法）"""
        try:
            logger.info(f"使用分步方法转换SAM到BAM: {sample_name}")

            # 中间BAM文件
            temp_bam = bam_file.replace('.bam', '_temp.bam')

            # 步骤1: SAM到BAM
            logger.debug(f"步骤1: SAM到中间BAM")
            view_cmd = ["samtools", "view", "-bS", "-o", temp_bam, sam_file]
            view_result = subprocess.run(view_cmd, capture_output=True, text=True)

            if view_result.returncode != 0:
                logger.error(f"SAM到BAM转换失败: {view_result.stderr}")
                Path(temp_bam).unlink(missing_ok=True)
                return False

            # 步骤2: 排序
            logger.debug(f"步骤2: 排序BAM")
            sort_cmd = ["samtools", "sort", "-@", str(threads), "-o", bam_file, temp_bam]
            sort_result = subprocess.run(sort_cmd, capture_output=True, text=True)

            # 删除中间文件
            Path(temp_bam).unlink(missing_ok=True)

            if sort_result.returncode != 0:
                logger.error(f"BAM排序失败: {sort_result.stderr}")
                return False

            # 删除临时SAM文件
            if not keep_sam:
                Path(sam_file).unlink(missing_ok=True)

            return True

        except Exception as e:
            logger.error(f"分步转换失败: {e}")
            return False

    def _index_bam(self, bam_file: str, sample_name: str) -> bool:
        """创建BAM索引"""
        try:
            logger.info(f"创建BAM索引: {sample_name}")
            result = subprocess.run(
                ["samtools", "index", bam_file],
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode == 0:
                logger.info(f"BAM索引已创建: {bam_file}.bai")
                return True
            else:
                logger.warning(f"BAM索引创建失败: {result.stderr}")
                return False

        except Exception as e:
            logger.error(f"创建BAM索引时出错: {e}")
            return False

    def _parse_alignment_stats(self, log_file: Path, log_content: str) -> Dict[str, Any]:
        """解析比对统计信息"""
        stats = {
            'total_reads': 0,
            'aligned_reads': 0,
            'alignment_rate': 0.0,
            'unique_alignments': 0,
            'multiple_alignments': 0,
            'failed_to_align': 0
        }

        try:
            if log_content:
                lines = log_content.split('\n')
            else:
                with open(log_file, 'r') as f:
                    lines = f.readlines()

            for line in lines:
                line = line.strip()
                # 解析Bowtie2统计信息
                if 'reads; of these:' in line:
                    stats['total_reads'] = int(line.split()[0])
                elif 'aligned 0 times' in line:
                    stats['failed_to_align'] = int(line.split()[0])
                elif 'aligned exactly 1 time' in line:
                    stats['unique_alignments'] = int(line.split()[0])
                elif 'aligned >1 times' in line:
                    stats['multiple_alignments'] = int(line.split()[0])
                elif 'overall alignment rate' in line:
                    rate_str = line.split()[0].replace('%', '')
                    stats['alignment_rate'] = float(rate_str) / 100.0

            # 计算总比对reads
            stats['aligned_reads'] = stats['unique_alignments'] + stats['multiple_alignments']

        except Exception as e:
            logger.warning(f"解析比对统计时出错: {e}")

        return stats

    def _save_alignment_stats(self, stats: Dict[str, Any], stats_file: Path):
        """保存比对统计信息"""
        try:
            with open(stats_file, 'w') as f:
                f.write("=== Bowtie2比对统计 ===\n\n")
                f.write(f"总reads数: {stats.get('total_reads', 0):,}\n")
                f.write(f"比对reads数: {stats.get('aligned_reads', 0):,}\n")
                f.write(f"比对率: {stats.get('alignment_rate', 0.0)*100:.1f}%\n")
                f.write(f"唯一比对: {stats.get('unique_alignments', 0):,}\n")
                f.write(f"多重比对: {stats.get('multiple_alignments', 0):,}\n")
                f.write(f"未比对: {stats.get('failed_to_align', 0):,}\n")
        except Exception as e:
            logger.warning(f"保存统计信息时出错: {e}")

    def generate_summary(self, output_dir: str) -> Optional[str]:
        """
        生成比对汇总报告

        参数:
            output_dir: 输出目录

        返回:
            str: 报告文件路径
        """
        if not self.results:
            logger.warning("无比对结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建汇总数据
        summary_data = []
        for sample, result in self.results.items():
            stats = result.get('stats', {})
            summary_data.append({
                'sample': sample,
                'type': result.get('type', 'unknown'),
                'success': result.get('success', False),
                'total_reads': stats.get('total_reads', 0),
                'aligned_reads': stats.get('aligned_reads', 0),
                'alignment_rate': stats.get('alignment_rate', 0.0),
                'unique_alignments': stats.get('unique_alignments', 0),
                'multiple_alignments': stats.get('multiple_alignments', 0),
                'bam_file': result.get('bam_file', ''),
                'bam_indexed': result.get('stats', {}).get('bam_indexed', False)
            })

        # 保存为DataFrame
        df = pd.DataFrame(summary_data)

        # 保存为CSV
        csv_file = output_dir / "bowtie2_alignment_summary.csv"
        df.to_csv(csv_file, index=False)

        # 保存为JSON
        json_file = output_dir / "bowtie2_alignment_summary.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        # 生成文本报告
        report_file = output_dir / "bowtie2_alignment_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== Bowtie2比对汇总报告 ===\n\n")
            f.write(f"总样本数: {len(self.results)}\n")
            f.write(f"成功比对: {len([r for r in self.results.values() if r.get('success')])}\n\n")

            f.write("样本比对统计:\n")
            for _, row in df.iterrows():
                if row['success']:
                    f.write(f"  - {row['sample']}: {row['alignment_rate']*100:.1f}% "
                           f"({row['aligned_reads']:,}/{row['total_reads']:,} reads)\n")

            if 'alignment_rate' in df.columns:
                avg_rate = df[df['success']]['alignment_rate'].mean() * 100
                f.write(f"\n平均比对率: {avg_rate:.1f}%\n")

            f.write(f"\n报告生成时间: {pd.Timestamp.now()}\n")

        logger.info(f"比对汇总报告已生成: {report_file}")
        return str(report_file)


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """加载配置文件"""
    default_config = {
        'threads': 4,
        'small_rna_mode': True,
        'phred': 33,
        'sensitive': False,
        'very_sensitive': True,
        'report_all': False,
        'report_k': 1,
        'pe_max_insert': 1000,
        'bowtie2_path': 'bowtie2',
        'min_len': 18,
        'max_len': 35,
        'keep_sam': False
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
        description="Bowtie2序列比对脚本 - 针对small RNA测序数据优化"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入文件或目录，或样本信息CSV文件"
    )

    parser.add_argument(
        "--index", "-x",
        required=True,
        help="Bowtie2索引文件前缀"
    )

    parser.add_argument(
        "--output", "-o",
        default="results/alignment/bam",
        help="输出目录 (默认: results/alignment/bam)"
    )

    parser.add_argument(
        "--config", "-c",
        help="配置文件 (YAML格式)"
    )

    parser.add_argument(
        "--sample-info",
        help="样本信息CSV文件，包含fastq_r1和fastq_r2列"
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
        help="生成汇总报告"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    # 加载配置
    config = load_config(args.config)
    config['threads'] = args.threads

    # 初始化比对器
    aligner = Bowtie2Aligner(bowtie2_path=config.get('bowtie2_path', 'bowtie2'))

    # 检查Bowtie2
    if not aligner.check_bowtie2():
        logger.error("Bowtie2检查失败，请确保Bowtie2已安装")
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
            fastq_r1 = row.get('fastq_r1')
            fastq_r2 = row.get('fastq_r2')

            if pd.notna(fastq_r1) and pd.notna(fastq_r2):
                # 双端数据
                aligner.align_paired_end(
                    fastq_r1=fastq_r1,
                    fastq_r2=fastq_r2,
                    index_prefix=args.index,
                    output_dir=output_dir,
                    sample_name=sample,
                    config=config
                )
            elif pd.notna(fastq_r1):
                # 单端数据
                aligner.align_single_end(
                    fastq_file=fastq_r1,
                    index_prefix=args.index,
                    output_dir=output_dir,
                    sample_name=sample,
                    config=config
                )

    elif input_path.is_file() and input_path.suffix in ['.fastq', '.fq', '.fastq.gz', '.fq.gz']:
        # 单个fastq文件
        sample_name = input_path.stem
        if sample_name.endswith('.fastq') or sample_name.endswith('.fq'):
            sample_name = sample_name[:-6]
        if sample_name.endswith('_R1') or sample_name.endswith('_R2'):
            sample_name = sample_name[:-3]

        aligner.align_single_end(
            fastq_file=str(input_path),
            index_prefix=args.index,
            output_dir=output_dir,
            sample_name=sample_name,
            config=config
        )

    else:
        logger.error("不支持的输入类型，请提供样本信息CSV或fastq文件")
        sys.exit(1)

    # 生成汇总报告
    if args.summary or True:  # 默认总是生成汇总
        report_file = aligner.generate_summary(output_dir)
        if report_file:
            logger.info(f"比对完成，报告文件: {report_file}")
        else:
            logger.warning("未能生成汇总报告")

    logger.info("Bowtie2比对流程完成")


if __name__ == "__main__":
    main()