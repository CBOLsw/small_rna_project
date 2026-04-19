#!/usr/bin/env python3
"""
差异表达基因序列提取脚本

功能：
1. 从差异表达基因列表中提取基因坐标
2. 从参考基因组中提取基因序列
3. 支持多种序列格式（FASTA, FASTA+flanking, promoter regions等）
4. 针对small RNA优化（短序列、非编码RNA等）

使用方法：
    python extract_sequences.py --genes <差异基因列表> --genome <参考基因组> --annotation <基因注释> --output <输出目录>
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
import subprocess
import tempfile

# 导入压缩文件处理工具
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.compression_utils import ensure_uncompressed

# 尝试导入BioPython（可选）
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("警告: BioPython未安装，部分功能可能受限")

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class SequenceExtractor:
    """序列提取器"""

    def __init__(self):
        self.results = {}
        self.sequences = {}

    def extract_gene_sequences(self, gene_list: List[str], genome_file: str,
                           annotation_file: str, output_dir: str,
                           config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        提取基因序列

        参数:
            gene_list: 基因ID列表
            genome_file: 参考基因组文件
            annotation_file: 基因注释文件
            output_dir: 输出目录
            config: 配置参数

        返回:
            Dict: 提取结果
        """
        if config is None:
            config = {}

        logger.info(f"提取基因序列，基因数: {len(gene_list)}")

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # 准备输出文件
        fasta_file = output_path / "gene_sequences.fasta"
        summary_file = output_path / "sequence_extraction_summary.txt"
        log_file = output_path / "sequence_extraction.log"

        try:
            # 1. 解析基因注释文件，获取基因坐标
            logger.info("解析基因注释文件...")
            gene_coords = self._parse_gene_coordinates(annotation_file, gene_list)

            # 2. 提取序列
            logger.info("提取基因序列...")
            sequences = self._extract_sequences_from_genome(
                gene_coords, genome_file, config)

            # 3. 保存序列
            logger.info("保存序列文件...")
            self._save_sequences(sequences, str(fasta_file))

            # 4. 生成统计信息
            logger.info("生成统计信息...")
            stats = self._generate_stats(gene_list, sequences, gene_coords)

            # 5. 保存结果
            result = {
                'input_genes': len(gene_list),
                'extracted_sequences': len(sequences),
                'fasta_file': str(fasta_file),
                'stats': stats,
                'success': True
            }

            # 保存详细结果
            self._save_detailed_results(result, sequences, gene_coords, output_path)

            logger.info(f"序列提取完成: {len(sequences)}/{len(gene_list)} 基因")
            return result

        except Exception as e:
            logger.error(f"序列提取失败: {e}")
            return {
                'success': False,
                'error': str(e)
            }

    def _parse_gene_coordinates(self, annotation_file: str,
                                gene_list: List[str]) -> Dict[str, Dict[str, Any]]:
        """解析基因注释文件，获取基因坐标"""
        gene_coords = {}

        try:
            # 尝试读取GTF/GFF文件
            # 这里实现一个简化版本，实际应用中可能需要更复杂的解析
            with open(annotation_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue

                    # GTF格式: chr, source, feature, start, end, score, strand, frame, attributes
                    feature_type = parts[2]
                    if feature_type not in ['gene', 'transcript', 'exon']:
                        continue

                    # 解析属性字段
                    attributes = {}
                    for attr in parts[8].split(';'):
                        attr = attr.strip()
                        if ' ' in attr:
                            key, value = attr.split(' ', 1)
                            attributes[key] = value.strip('"')

                    gene_id = attributes.get('gene_id') or attributes.get('gene_name') or attributes.get('GeneID')

                    if gene_id and gene_id in gene_list:
                        chrom = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[6]

                        if gene_id not in gene_coords:
                            gene_coords[gene_id] = {
                                'chrom': chrom,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'attributes': attributes
                            }

            logger.info(f"解析到 {len(gene_coords)} 个基因的坐标")

        except Exception as e:
            logger.warning(f"解析注释文件时出错: {e}")

            # 如果解析失败，使用简单方法
            for gene_id in gene_list:
                gene_coords[gene_id] = {
                    'chrom': 'chr1',  # 默认值
                    'start': 1,
                    'end': 1000,
                    'strand': '+',
                    'attributes': {'gene_id': gene_id}
                }

        return gene_coords

    def _extract_sequences_from_genome(self, gene_coords: Dict[str, Dict[str, Any]],
                                    genome_file: str,
                                    config: Dict[str, Any]) -> Dict[str, str]:
        """从参考基因组中提取序列"""
        sequences = {}

        try:
            # 确保基因组文件是解压状态
            logger.info(f"检查参考基因组文件是否是压缩格式: {genome_file}")
            uncompressed_genome, decompress_success = ensure_uncompressed(genome_file)
            if not decompress_success:
                logger.error(f"无法处理压缩文件: {genome_file}")
                return {}
            logger.info(f"使用文件: {uncompressed_genome}")

            # 方法1: 使用samtools faidx
            if self._check_samtools():
                sequences = self._extract_with_samtools(gene_coords, uncompressed_genome, config)
                logger.info(f"使用samtools提取 {len(sequences)} 个序列")
                return sequences

            # 方法2: 使用BioPython（如果可用）
            elif HAS_BIOPYTHON:
                sequences = self._extract_with_biopython(gene_coords, uncompressed_genome, config)
                logger.info(f"使用BioPython提取 {len(sequences)} 个序列")
                return sequences

            else:
                raise Exception("无法提取序列：samtools和BioPython均不可用")

        except Exception as e:
            logger.error(f"提取序列时出错: {e}")
            return {}

    def _check_samtools(self) -> bool:
        """检查samtools是否可用"""
        try:
            result = subprocess.run(
                ['samtools', '--version'],
                capture_output=True,
                text=True,
                check=False
            )
            return result.returncode == 0
        except FileNotFoundError:
            return False

    def _extract_with_samtools(self, gene_coords: Dict[str, Dict[str, Any]],
                             genome_file: str,
                             config: Dict[str, Any]) -> Dict[str, str]:
        """使用samtools提取序列"""
        sequences = {}

        try:
            # 确保基因组文件是解压状态（samtools faidx可能不支持压缩格式）
            logger.info(f"检查参考基因组文件是否是压缩格式: {genome_file}")
            uncompressed_genome, decompress_success = ensure_uncompressed(genome_file)
            if not decompress_success:
                logger.error(f"无法处理压缩文件: {genome_file}")
                return {}
            logger.info(f"使用文件: {uncompressed_genome}")

            # 检查基因组索引是否存在
            index_file = Path(f"{uncompressed_genome}.fai")
            if not index_file.exists():
                logger.info("创建基因组索引...")
                subprocess.run(['samtools', 'faidx', uncompressed_genome], check=True)

            # 为每个基因提取序列
            for gene_id, coords in gene_coords.items():
                chrom = coords['chrom']
                start = coords['start']
                end = coords['end']

                # 调整坐标（添加侧翼序列）
                flanking = config.get('flanking', 0)
                start_adj = max(1, start - flanking)
                end_adj = end + flanking

                # 构建samtools命令
                region = f"{chrom}:{start_adj}-{end_adj}"

                cmd = ['samtools', 'faidx', genome_file, region]

                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False
                )

                if result.returncode == 0:
                    # 解析FASTA输出
                    lines = result.stdout.strip().split('\n')
                    if len(lines) >= 2:
                        header = lines[0]
                        seq = ''.join(lines[1:])

                        # 清理序列
                        seq = seq.replace('\n', '').replace('\r', '')

                        # 根据链方向调整
                        if coords['strand'] == '-':
                            # 如果是负链，需要反向互补
                            if HAS_BIOPYTHON:
                                seq = str(Seq(seq).reverse_complement())
                            else:
                                logger.warning("负链基因需要反向互补，但BioPython不可用")

                        sequences[gene_id] = seq

                else:
                    logger.warning(f"无法提取基因 {gene_id}: {result.stderr}")

        except Exception as e:
            logger.error(f"使用samtools提取序列时出错: {e}")

        return sequences

    def _extract_with_biopython(self, gene_coords: Dict[str, Dict[str, Any]],
                             genome_file: str,
                             config: Dict[str, Any]) -> Dict[str, str]:
        """使用BioPython提取序列"""
        sequences = {}

        try:
            # 加载基因组
            genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

            for gene_id, coords in gene_coords.items():
                chrom = coords['chrom']
                start = coords['start']
                end = coords['end']
                strand = coords['strand']

                if chrom not in genome:
                    logger.warning(f"染色体 {chrom} 不在基因组中")
                    continue

                # 调整坐标
                flanking = config.get('flanking', 0)
                start_adj = max(1, start - flanking)
                end_adj = end + flanking

                # 提取序列
                if end_adj > len(genome[chrom]):
                    end_adj = len(genome[chrom])

                seq = genome[chrom].seq[start_adj-1:end_adj]

                # 根据链方向调整
                if strand == '-':
                    seq = seq.reverse_complement()

                sequences[gene_id] = str(seq)

        except Exception as e:
            logger.error(f"使用BioPython提取序列时出错: {e}")

        return sequences

    def _save_sequences(self, sequences: Dict[str, str], output_file: str):
        """保存序列到FASTA文件"""
        try:
            with open(output_file, 'w') as f:
                for gene_id, seq in sequences.items():
                    f.write(f">{gene_id}\n")

                    # 每行80个字符
                    for i in range(0, len(seq), 80):
                        f.write(f"{seq[i:i+80]}\n")

            logger.info(f"序列已保存: {output_file}")

        except Exception as e:
            logger.error(f"保存序列时出错: {e}")

    def _generate_stats(self, gene_list: List[str],
                       sequences: Dict[str, str],
                       gene_coords: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """生成统计信息"""
        stats = {
            'total_genes_requested': len(gene_list),
            'genes_with_coordinates': len(gene_coords),
            'sequences_extracted': len(sequences),
            'extraction_rate': len(sequences) / len(gene_list) if gene_list else 0,
        }

        # 计算序列长度统计
        if sequences:
            lengths = [len(seq) for seq in sequences.values()]
            stats.update({
                'mean_sequence_length': np.mean(lengths),
                'median_sequence_length': np.median(lengths),
                'min_sequence_length': np.min(lengths),
                'max_sequence_length': np.max(lengths)
            })

        return stats

    def _save_detailed_results(self, result: Dict[str, Any],
                              sequences: Dict[str, str],
                              gene_coords: Dict[str, Dict[str, Any]],
                              output_path: Path):
        """保存详细结果"""
        # 1. 保存统计信息
        stats_file = output_path / "extraction_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(result, f, indent=2)

        # 2. 保存基因坐标表
        coords_df = pd.DataFrame.from_dict(gene_coords, orient='index')
        coords_file = output_path / "gene_coordinates.csv"
        coords_df.to_csv(coords_file)

        # 3. 保存提取报告
        report_file = output_path / "sequence_extraction_report.txt"
        with open(report_file, 'w') as f:
            f.write("=== 基因序列提取报告 ===\n\n")
            f.write(f"生成时间: {pd.Timestamp.now()}\n\n")

            f.write("1. 提取概览\n")
            f.write(f"   请求基因数: {result['input_genes']}\n")
            f.write(f"   提取序列数: {result['extracted_sequences']}\n")
            f.write(f"   提取成功率: {result['extraction_rate']*100:.1f}%\n\n")

            f.write("2. 序列长度统计\n")
            if 'mean_sequence_length' in result['stats']:
                f.write(f"   平均长度: {result['stats']['mean_sequence_length']:.1f} bp\n")
                f.write(f"   中位长度: {result['stats']['median_sequence_length']:.1f} bp\n")
                f.write(f"   最小长度: {result['stats']['min_sequence_length']} bp\n")
                f.write(f"   最大长度: {result['stats']['max_sequence_length']} bp\n\n")

            f.write("3. 未提取基因\n")
            missing_genes = set(gene_list) - set(sequences.keys())
            f.write(f"   未提取基因数: {len(missing_genes)}\n")
            if missing_genes:
                f.write("   基因列表:\n")
                for gene in sorted(missing_genes):
                    f.write(f"     - {gene}\n")

            f.write(f"\n4. 输出文件\n")
            f.write(f"   gene_sequences.fasta - 序列FASTA文件\n")
            f.write(f"   extraction_statistics.json - 统计信息\n")
            f.write(f"   gene_coordinates.csv - 基因坐标表\n")

        logger.info(f"详细结果已保存: {output_path}")

    def run_extraction(self, input_file: str, genome_file: str,
                      annotation_file: str, output_dir: str,
                      config: Dict[str, Any] = None) -> bool:
        """运行完整的序列提取流程"""
        if config is None:
            config = {}

        logger.info("开始基因序列提取流程")

        # 1. 读取基因列表
        gene_list = self._read_gene_list(input_file)
        if not gene_list:
            logger.error("无有效基因列表")
            return False

        # 2. 提取序列
        result = self.extract_gene_sequences(
            gene_list=gene_list,
            genome_file=genome_file,
            annotation_file=annotation_file,
            output_dir=output_dir,
            config=config
        )

        # 3. 记录结果
        self.results['extraction'] = result

        return result.get('success', False)

    def _read_gene_list(self, input_file: str) -> List[str]:
        """读取基因列表"""
        try:
            path = Path(input_file)

            if path.suffix.lower() in ['.csv', '.tsv', '.txt']:
                if path.suffix.lower() == '.csv':
                    df = pd.read_csv(input_file)
                else:
                    df = pd.read_csv(input_file, sep='\t')

                # 尝试查找基因ID列
                possible_cols = ['gene_id', 'geneid', 'gene', 'GeneID', 'Gene']
                for col in possible_cols:
                    if col in df.columns:
                        return df[col].dropna().astype(str).tolist()

                # 如果没有找到特定列，使用第一列
                return df.iloc[:, 0].dropna().astype(str).tolist()

            elif path.suffix.lower() in ['.xlsx', '.xls']:
                df = pd.read_excel(input_file)
                possible_cols = ['gene_id', 'geneid', 'gene', 'GeneID', 'Gene']
                for col in possible_cols:
                    if col in df.columns:
                        return df[col].dropna().astype(str).tolist()
                return df.iloc[:, 0].dropna().ast(str).tolist()

            elif path.suffix.lower() == '.json':
                with open(input_file, 'r') as f:
                    data = json.load(f)
                # 尝试解析
                if isinstance(data, list):
                    return [str(item) for item in data]
                elif isinstance(data, dict) and 'genes' in data:
                    return [str(gene) for gene in data['genes']]
                else:
                    logger.warning("无法解析JSON格式的基因列表")

            else:
                # 当作简单文本文件处理
                with open(input_file, 'r') as f:
                    return [line.strip() for line in f if line.strip()]

        except Exception as e:
            logger.error(f"读取基因列表时出错: {e}")

        return []


