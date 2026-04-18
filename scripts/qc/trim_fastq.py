#!/usr/bin/env python3
"""
Trimmomatic质量修剪脚本 - 用于small RNA测序数据预处理

功能：
1. 对fastq文件进行质量修剪和接头去除
2. 特别优化small RNA测序数据参数
3. 处理单端和双端测序数据
4. 生成修剪统计报告

使用方法：
    python trim_fastq.py --input <fastq文件或样本信息> --output <输出目录>
"""

import os
import sys
import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any
import logging
import json
import yaml

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class TrimmomaticProcessor:
    """Trimmomatic处理器类"""

    # small RNA测序常用接头序列
    ADAPTERS = {
        'illumina_universal': 'AGATCGGAAGAG',
        'illumina_small_rna': 'TGGAATTCTCGG',
        'nextera': 'CTGTCTCTTATA',
        'truseq': 'AGATCGGAAGAGC',
        'vahts_small_rna_v2': 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'  # VAHTS Small RNA Library Prep Kit for Illumina V2
    }

    def __init__(self, trimmomatic_path: str = "trimmomatic",
                 java_path: str = "java"):
        """
        初始化Trimmomatic处理器

        参数:
            trimmomatic_path: Trimmomatic可执行jar路径
            java_path: Java可执行文件路径
        """
        self.trimmomatic_path = trimmomatic_path
        self.java_path = java_path
        self.results = {}

    def check_trimmomatic(self) -> bool:
        """检查Trimmomatic是否可用"""
        try:
            # 检查Java
            java_check = subprocess.run(
                [self.java_path, "-version"],
                capture_output=True,
                text=True,
                check=False
            )
            if java_check.returncode != 0:
                logger.error("Java检查失败")
                return False

            # 检查Trimmomatic
            cmd = [
                self.java_path, "-jar", self.trimmomatic_path,
                "-version"
            ]
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )
            if result.returncode == 0:
                logger.info(f"Trimmomatic版本: {result.stdout.strip()}")
                return True
            else:
                logger.error(f"Trimmomatic检查失败: {result.stderr}")
                return False
        except FileNotFoundError as e:
            logger.error(f"未找到可执行文件: {e}")
            return False

    def trim_single_end(self, input_file: str, output_dir: str,
                        sample_name: str, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        处理单端测序数据

        参数:
            input_file: 输入fastq文件
            output_dir: 输出目录
            sample_name: 样本名称
            config: 配置参数

        返回:
            Dict: 处理结果
        """
        logger.info(f"处理单端样本: {sample_name}")

        # 构建输出文件路径
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        output_file = output_dir / f"{sample_name}_trimmed.fastq.gz"
        log_file = output_dir / f"{sample_name}_trimmomatic.log"

        # 构建Trimmomatic命令
        cmd = self._build_single_end_command(
            input_file, output_file, log_file, config
        )

        # 执行命令
        success, stats = self._run_trimmomatic(cmd, log_file, sample_name)

        result = {
            'sample': sample_name,
            'input_file': input_file,
            'output_file': str(output_file),
            'log_file': str(log_file),
            'success': success,
            'type': 'single_end',
            'stats': stats
        }

        self.results[sample_name] = result
        return result

    def trim_paired_end(self, input_r1: str, input_r2: str, output_dir: str,
                        sample_name: str, config: Dict[str, Any]) -> Dict[str, Any]:
        """
        处理双端测序数据

        参数:
            input_r1: R1输入文件
            input_r2: R2输入文件
            output_dir: 输出目录
            sample_name: 样本名称
            config: 配置参数

        返回:
            Dict: 处理结果
        """
        logger.info(f"处理双端样本: {sample_name}")

        # 构建输出文件路径
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 输出文件命名
        output_r1_paired = output_dir / f"{sample_name}_R1_paired.fastq.gz"
        output_r1_unpaired = output_dir / f"{sample_name}_R1_unpaired.fastq.gz"
        output_r2_paired = output_dir / f"{sample_name}_R2_paired.fastq.gz"
        output_r2_unpaired = output_dir / f"{sample_name}_R2_unpaired.fastq.gz"
        log_file = output_dir / f"{sample_name}_trimmomatic.log"

        # 构建Trimmomatic命令
        cmd = self._build_paired_end_command(
            input_r1, input_r2,
            output_r1_paired, output_r1_unpaired,
            output_r2_paired, output_r2_unpaired,
            log_file, config
        )

        # 执行命令
        success, stats = self._run_trimmomatic(cmd, log_file, sample_name)

        result = {
            'sample': sample_name,
            'input_r1': input_r1,
            'input_r2': input_r2,
            'output_r1_paired': str(output_r1_paired),
            'output_r1_unpaired': str(output_r1_unpaired),
            'output_r2_paired': str(output_r2_paired),
            'output_r2_unpaired': str(output_r2_unpaired),
            'log_file': str(log_file),
            'success': success,
            'type': 'paired_end',
            'stats': stats
        }

        self.results[sample_name] = result
        return result

    def _build_single_end_command(self, input_file: Path, output_file: Path,
                                  log_file: Path, config: Dict[str, Any]) -> List[str]:
        """构建单端Trimmomatic命令"""
        cmd = [
            self.java_path, "-jar", self.trimmomatic_path,
            "SE",
            "-threads", str(config.get('threads', 4)),
            "-phred33",  # small RNA通常使用phred33
            str(input_file),
            str(output_file)
        ]

        # 添加Trimmomatic步骤
        cmd.extend(self._build_trimmomatic_steps(config))

        # 添加日志重定向
        cmd.extend(["-trimlog", str(log_file)])

        return cmd

    def _build_paired_end_command(self, input_r1: Path, input_r2: Path,
                                  output_r1_paired: Path, output_r1_unpaired: Path,
                                  output_r2_paired: Path, output_r2_unpaired: Path,
                                  log_file: Path, config: Dict[str, Any]) -> List[str]:
        """构建双端Trimmomatic命令"""
        cmd = [
            self.java_path, "-jar", self.trimmomatic_path,
            "PE",
            "-threads", str(config.get('threads', 4)),
            "-phred33",
            str(input_r1),
            str(input_r2),
            str(output_r1_paired),
            str(output_r1_unpaired),
            str(output_r2_paired),
            str(output_r2_unpaired)
        ]

        # 添加Trimmomatic步骤
        cmd.extend(self._build_trimmomatic_steps(config))

        # 添加日志重定向
        cmd.extend(["-trimlog", str(log_file)])

        return cmd

    def _build_trimmomatic_steps(self, config: Dict[str, Any]) -> List[str]:
        """构建Trimmomatic处理步骤"""
        steps = []

        # 1. 接头去除 (针对small RNA优化)
        adapter = config.get('adapter', 'illumina_small_rna')
        adapter_seq = self.ADAPTERS.get(adapter, adapter)
        steps.extend(["ILLUMINACLIP", f"{adapter_seq}:2:30:10"])

        # 2. 滑动窗口质量修剪
        window_size = config.get('window_size', 4)  # small RNA窗口较小
        required_quality = config.get('required_quality', 20)
        steps.extend(["SLIDINGWINDOW", f"{window_size}:{required_quality}"])

        # 3. 前导质量修剪
        leading_quality = config.get('leading_quality', 20)
        steps.extend(["LEADING", str(leading_quality)])

        # 4. 末尾质量修剪
        trailing_quality = config.get('trailing_quality', 20)
        steps.extend(["TRAILING", str(trailing_quality)])

        # 5. 最小长度筛选 (small RNA通常18-35nt)
        min_length = config.get('min_length', 18)
        steps.extend(["MINLEN", str(min_length)])

        # 6. 裁剪到固定长度 (可选，用于标准化)
        if config.get('crop_length'):
            steps.extend(["CROP", str(config['crop_length'])])

        return steps

    def _run_trimmomatic(self, cmd: List[str], log_file: Path,
                         sample_name: str) -> Tuple[bool, Dict[str, Any]]:
        """运行Trimmomatic命令并解析结果"""
        logger.info(f"运行命令: {' '.join(cmd)}")

        try:
            with open(log_file, 'w') as log_f:
                result = subprocess.run(
                    cmd,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    text=True,
                    check=False
                )

            if result.returncode == 0:
                logger.info(f"Trimmomatic处理成功: {sample_name}")
                stats = self._parse_trimmomatic_log(log_file, sample_name)
                return True, stats
            else:
                logger.error(f"Trimmomatic处理失败: {sample_name}")
                return False, {}

        except Exception as e:
            logger.error(f"运行Trimmomatic时出错: {e}")
            return False, {}

    def _parse_trimmomatic_log(self, log_file: Path, sample_name: str) -> Dict[str, Any]:
        """解析Trimmomatic日志文件"""
        stats = {
            'input_reads': 0,
            'surviving_reads': 0,
            'dropped_reads': 0,
            'survival_rate': 0.0
        }

        try:
            with open(log_file, 'r') as f:
                lines = f.readlines()

            # 查找关键统计信息
            for line in lines:
                line = line.strip()
                if 'Input Read Pairs:' in line:
                    parts = line.split()
                    stats['input_reads'] = int(parts[3])
                    stats['surviving_reads'] = int(parts[6])
                    stats['dropped_reads'] = int(parts[9])
                    if stats['input_reads'] > 0:
                        stats['survival_rate'] = stats['surviving_reads'] / stats['input_reads']

        except Exception as e:
            logger.warning(f"解析日志文件失败: {e}")

        return stats

    def generate_summary(self, output_dir: str) -> Optional[str]:
        """
        生成汇总报告

        参数:
            output_dir: 输出目录

        返回:
            str: 汇总报告文件路径
        """
        if not self.results:
            logger.warning("无处理结果可汇总")
            return None

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # 创建汇总DataFrame
        summary_data = []
        for sample, result in self.results.items():
            stats = result.get('stats', {})
            summary_data.append({
                'sample': sample,
                'type': result.get('type', 'unknown'),
                'success': result.get('success', False),
                'input_reads': stats.get('input_reads', 0),
                'surviving_reads': stats.get('surviving_reads', 0),
                'dropped_reads': stats.get('dropped_reads', 0),
                'survival_rate': stats.get('survival_rate', 0.0),
                'output_file': result.get('output_file', '') if result.get('type') == 'single_end'
                else result.get('output_r1_paired', '')
            })

        df = pd.DataFrame(summary_data)

        # 保存为CSV
        csv_file = output_dir / "trimmomatic_summary.csv"
        df.to_csv(csv_file, index=False)

        # 保存为JSON
        json_file = output_dir / "trimmomatic_summary.json"
        with open(json_file, 'w') as f:
            json.dump(self.results, f, indent=2)

        # 生成简单文本报告
        report_file = output_dir / "trimmomatic_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== Trimmomatic质量修剪汇总报告 ===\n\n")
            f.write(f"总样本数: {len(self.results)}\n")
            f.write(f"成功处理: {len([r for r in self.results.values() if r.get('success')])}\n\n")

            f.write("样本统计:\n")
            for sample, result in self.results.items():
                if result.get('success'):
                    stats = result.get('stats', {})
                    f.write(f"  - {sample}: {stats.get('input_reads', 0)} reads, "
                           f"保留 {stats.get('surviving_reads', 0)} "
                           f"({stats.get('survival_rate', 0.0)*100:.1f}%)\n")

            f.write(f"\n报告生成时间: {pd.Timestamp.now()}\n")

        logger.info(f"汇总报告已生成: {report_file}")
        return str(report_file)


def load_config(config_file: Optional[str] = None) -> Dict[str, Any]:
    """加载配置文件"""
    default_config = {
        'threads': 4,
        'adapter': 'vahts_small_rna_v2',  # 默认使用VAHTS接头
        'window_size': 4,
        'required_quality': 20,
        'leading_quality': 20,
        'trailing_quality': 20,
        'min_length': 18,
        'max_length': 35,
        'crop_length': None,
        'trimmomatic_path': 'trimmomatic',
        'java_path': 'java'
    }

    if config_file and Path(config_file).exists():
        try:
            with open(config_file, 'r') as f:
                user_config = yaml.safe_load(f)
                # 从 config.yaml 的 quality_control.trimmomatic 部分读取配置
                if 'quality_control' in user_config and 'trimmomatic' in user_config['quality_control']:
                    trimmomatic_config = user_config['quality_control']['trimmomatic']
                    if 'leading' in trimmomatic_config:
                        default_config['leading_quality'] = trimmomatic_config['leading']
                    if 'trailing' in trimmomatic_config:
                        default_config['trailing_quality'] = trimmomatic_config['trailing']
                    if 'slidingwindow' in trimmomatic_config:
                        # slidingwindow 格式为 "4:15"，需要解析
                        sw = trimmomatic_config['slidingwindow']
                        if ':' in sw:
                            window, quality = sw.split(':')
                            default_config['window_size'] = int(window)
                            default_config['required_quality'] = int(quality)
                    if 'minlen' in trimmomatic_config:
                        default_config['min_length'] = trimmomatic_config['minlen']
                    if 'adapter_type' in trimmomatic_config:
                        default_config['adapter'] = trimmomatic_config['adapter_type']
                logger.info(f"已加载配置文件: {config_file}")
        except Exception as e:
            logger.warning(f"加载配置文件失败: {e}，使用默认配置")

    return default_config


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="Trimmomatic质量修剪脚本 - small RNA测序数据预处理"
    )

    parser.add_argument(
        "--input", "-i",
        required=True,
        help="输入文件或目录，或样本信息CSV文件"
    )

    parser.add_argument(
        "--output", "-o",
        default="data/processed/trimmed",
        help="输出目录 (默认: data/processed/trimmed)"
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

    # 初始化处理器
    processor = TrimmomaticProcessor(
        trimmomatic_path=config['trimmomatic_path'],
        java_path=config['java_path']
    )

    # 检查Trimmomatic
    if not processor.check_trimmomatic():
        logger.error("Trimmomatic检查失败，请确保Trimmomatic和Java已安装")
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
                processor.trim_paired_end(
                    input_r1=fastq_r1,
                    input_r2=fastq_r2,
                    output_dir=output_dir,
                    sample_name=sample,
                    config=config
                )
            elif pd.notna(fastq_r1):
                # 单端数据
                processor.trim_single_end(
                    input_file=fastq_r1,
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

        processor.trim_single_end(
            input_file=str(input_path),
            output_dir=output_dir,
            sample_name=sample_name,
            config=config
        )

    else:
        logger.error("不支持的输入类型，请提供样本信息CSV或fastq文件")
        sys.exit(1)

    # 生成汇总报告
    if args.summary or True:  # 默认总是生成汇总
        report_file = processor.generate_summary(output_dir)
        if report_file:
            logger.info(f"处理完成，报告文件: {report_file}")
        else:
            logger.warning("未能生成汇总报告")

    logger.info("Trimmomatic处理流程完成")


if __name__ == "__main__":
    main()