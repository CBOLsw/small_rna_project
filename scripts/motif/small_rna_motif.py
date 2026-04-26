#!/usr/bin/env python3
"""
Small RNA motif分析脚本
用于从small RNA测序数据中发现富集的序列motif

正确流程：
1. 下载miRBase mature miRNA序列（如不存在）
2. 构建miRBase Bowtie2索引（如不存在）
3. 将trimmed reads比对到miRBase
4. 提取比对上的small RNA reads
5. 合并所有样本的reads，运行MEME motif发现

使用方法：
    python small_rna_motif.py --config config/config.yaml
"""

import os
import sys
import argparse
import subprocess
import logging
import json
import time
from pathlib import Path
from typing import List, Dict, Optional, Any, Tuple

import yaml
import pandas as pd

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.logging_utils import get_script_logger

logger = get_script_logger('small_rna_motif')


def ensure_mirbase_fasta(config: Dict[str, Any]) -> str:
    """确保miRBase mature miRNA序列文件存在，必要时下载"""
    mirbase_fasta = config.get('reference', {}).get('mirbase_fasta', 'references/hsa.mature.fa')

    if Path(mirbase_fasta).exists():
        logger.info(f"miRBase序列已存在: {mirbase_fasta}")
        return mirbase_fasta

    # 自动下载miRBase mature miRNA序列
    logger.info("未找到miRBase序列，开始下载...")
    os.makedirs(os.path.dirname(mirbase_fasta), exist_ok=True)

    import urllib.request

    # 优先使用预过滤的人类miRNA文件（GitHub镜像）
    # 备用：官方miRBase全物种文件（需过滤）
    urls = [
        "https://raw.githubusercontent.com/hank95179/biologic-data-analysis/master/hsa_mature.fa",
        "https://www.mirbase.org/download/mature.fa",
    ]

    downloaded = False
    for url in urls:
        try:
            logger.info(f"下载miRBase: {url}")
            response = urllib.request.urlopen(url, timeout=60)
            content = response.read().decode('utf-8')

            # 过滤人类miRNA（>hsa-*）
            lines = []
            for line in content.split('\n'):
                if line.startswith('>'):
                    if line.lower().startswith('>hsa-'):
                        lines.append(line)
                else:
                    lines.append(line)

            if len(lines) == 0:
                logger.warning(f"未从 {url} 找到人类miRNA，尝试下一个源")
                continue

            with open(mirbase_fasta, 'w') as f_out:
                f_out.write('\n'.join(lines))

            logger.info(f"miRBase人类成熟miRNA序列已保存: {mirbase_fasta} (共{sum(1 for l in lines if l.startswith('>'))}条)")
            downloaded = True
            break
        except Exception as e:
            logger.warning(f"下载失败 ({url}): {e}")
            continue

    if not downloaded:
        logger.error("所有miRBase下载源均失败")
        return mirbase_fasta

    return mirbase_fasta


def build_mirbase_index(mirbase_fasta: str, index_path: str, threads: int = 4) -> bool:
    """构建miRBase Bowtie2索引"""
    if Path(index_path + ".1.bt2").exists():
        logger.info(f"miRBase索引已存在: {index_path}")
        return True

    logger.info(f"构建miRBase Bowtie2索引: {mirbase_fasta}")
    try:
        cmd = [
            'bowtie2-build',
            '--threads', str(threads),
            mirbase_fasta,
            index_path
        ]
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("索引构建完成")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"索引构建失败: {e.stderr}")
        return False


def get_sample_fastq_r1(sample: str, config: Dict[str, Any]) -> str:
    """获取样本的R1 fastq文件路径"""
    metadata_file = config.get('samples', {}).get('metadata_file', 'data/metadata/sample_info.csv')
    if os.path.exists(metadata_file):
        df = pd.read_csv(metadata_file)
        sample_row = df[df[config['samples']['sample_column']] == sample]
        if not sample_row.empty and 'fastq_r1' in df.columns:
            return str(sample_row.iloc[0]['fastq_r1'])
    # 默认路径
    raw_fastq = config.get('directories', {}).get('raw_fastq', 'data/raw_fastq')
    return os.path.join(raw_fastq, f"{sample}_R1.fastq.gz")


def get_sample_trimmed_fastq(sample: str, config: Dict[str, Any]) -> str:
    """获取样本的trimmed fastq文件路径"""
    processed = config.get('directories', {}).get('processed', 'data/processed')
    return os.path.join(processed, f"{sample}_trimmed.fastq.gz")


def map_to_mirbase(trimmed_fastq: str, index_path: str, output_sam: str,
                   min_len: int = 18, max_len: int = 35,
                   max_mismatches: int = 1, threads: int = 4) -> Tuple[bool, int]:
    """
    将small RNA reads比对到miRBase

    返回: (成功标志, 比对上的read数)
    """
    if not Path(trimmed_fastq).exists():
        logger.warning(f"Trimmed文件不存在，跳过: {trimmed_fastq}")
        return False, 0

    logger.info(f"比对到miRBase: {trimmed_fastq}")

    # small RNA比对参数
    cmd = [
        'bowtie2',
        '--threads', str(threads),
        '-x', index_path,
        '-U', trimmed_fastq,
        '-S', output_sam,
        '--norc',  # 不比对反向互补链（miRNA有方向性）
        '-L', str(max(5, min_len // 2)),  # seed长度
        '-D', '20',  # 最大扩展次数
        '-R', '3',   # 重试次数
        '--score-min', f'C,{max_mismatches * -6},0',  # 允许的最大错配惩罚
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        # 解析比对结果
        mapped_reads = 0
        for line in result.stdout.split('\n'):
            if line.startswith('('):
                parts = line.split()
                for p in parts:
                    if p.endswith('reads;') and p[0].isdigit():
                        try:
                            num = int(p.rstrip('reads;'))
                            if 'unaligned' not in line and 'aligned concordantly' not in line:
                                mapped_reads = num
                        except:
                            pass

        # 直接统计SAM文件中的比对数
        if Path(output_sam).exists():
            mapped_reads = sum(1 for line in open(output_sam) if not line.startswith('@'))

        logger.info(f"比对完成，{mapped_reads}条reads比对到miRBase")
        return True, mapped_reads

    except subprocess.CalledProcessError as e:
        logger.error(f"miRBase比对失败: {e.stderr[:300] if e.stderr else ''}")
        return False, 0


def extract_mirna_reads(sam_file: str, min_len: int = 18, max_len: int = 35) -> List[Tuple[str, str]]:
    """
    从SAM文件中提取比对上的miRNA reads

    返回: [(read_id, sequence), ...]
    """
    reads = []
    if not Path(sam_file).exists():
        return reads

    # 从fastq文件读取序列（trimmed文件）
    # SAM只存储比对信息，需要从原始fastq提取序列
    # 这里用简化方法：假设序列在fastq文件中，需要读取原始fastq

    try:
        import gzip
        # 获取对应的fastq文件
        fastq_path = sam_file.replace('_mirbase.sam', '_trimmed.fastq.gz')
        if not Path(fastq_path).exists():
            fastq_path = sam_file.replace('.sam', '_trimmed.fastq.gz')

        # 建立read id到序列的映射
        read_seqs = {}
        if Path(fastq_path).exists():
            with gzip.open(fastq_path, 'rt') as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:  # header
                        read_id = line.strip().split()[0][1:]
                    elif i % 4 == 1:  # sequence
                        read_seqs[read_id] = line.strip()

        # 解析SAM，提取比对的reads
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                flag = int(parts[1])
                # 0x4: 未比对，跳过
                if flag & 0x4:
                    continue
                read_id = parts[0]
                seq = parts[9]
                # 长度过滤
                if min_len <= len(seq) <= max_len:
                    if read_id in read_seqs:
                        reads.append((read_id, read_seqs[read_id]))
    except Exception as e:
        logger.error(f"提取miRNA reads失败: {e}")

    return reads


def extract_mirna_reads_direct(trimmed_fastq: str, sam_file: str,
                               min_len: int = 18, max_len: int = 35) -> List[str]:
    """
    直接从trimmed fastq中提取与miRBase比对的reads序列
    使用minimap2或bowtie2的unaligned输出辅助
    """
    reads = []
    if not Path(trimmed_fastq).exists():
        return reads

    # 建立SAM中比对的read id集合
    mapped_ids = set()
    if Path(sam_file).exists():
        with open(sam_file, 'r') as f:
            for line in f:
                if line.startswith('@'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                flag = int(parts[1])
                if flag & 0x4:  # 未比对
                    continue
                mapped_ids.add(parts[0])

    # 从fastq提取这些reads的序列
    try:
        import gzip
        with gzip.open(trimmed_fastq, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 0:  # header
                    read_id = line.strip().split()[0][1:]
                elif i % 4 == 1:  # sequence
                    seq = line.strip()
                    if read_id in mapped_ids and min_len <= len(seq) <= max_len:
                        reads.append(seq)
    except Exception as e:
        logger.error(f"从fastq提取reads失败: {e}")

    return reads


def save_reads_to_fasta(reads: List[str], output_fasta: str, max_seqs: int = 0):
    """将reads保存为FASTA格式，带去重"""
    seen = set()
    unique_reads = []
    for seq in reads:
        if seq not in seen:
            seen.add(seq)
            unique_reads.append(seq)

    if max_seqs > 0 and len(unique_reads) > max_seqs:
        logger.info(f"序列过多({len(unique_reads)})，随机取{max_seqs}条")
        import random
        random.seed(42)
        unique_reads = random.sample(unique_reads, max_seqs)

    with open(output_fasta, 'w') as f:
        for i, seq in enumerate(unique_reads):
            f.write(f">read_{i}\n{seq}\n")

    logger.info(f"保存{len(unique_reads)}条唯一序列到 {output_fasta}")


def run_meme_on_small_rna(fasta_file: str, output_dir: str,
                          width_min: int = 5, width_max: int = 8,
                          max_motifs: int = 3, evalue_threshold: float = 1e-4,
                          min_sites: int = 10, max_sites: int = 100,
                          searchsize: int = 100000) -> Dict[str, Any]:
    """在small RNA reads上运行MEME"""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    cmd = [
        'meme',
        fasta_file,
        '-oc', str(output_path),
        '-dna',
        '-mod', 'zoops',
        '-nmotifs', str(max_motifs),
        '-minw', str(width_min),
        '-maxw', str(width_max),
        '-minsites', str(min_sites),
        '-maxsites', str(max_sites),
        '-evt', str(evalue_threshold),
        '-searchsize', str(searchsize),
        '-revcomp',
    ]

    logger.info(f"运行MEME: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=600)
        logger.info("MEME分析完成")

        # 解析结果
        motifs = parse_meme_result(output_path / 'meme.txt') if (output_path / 'meme.txt').exists() else []

        # 去重（去除正向/反向互补重复的motif）
        unique_motifs = deduplicate_motifs(motifs)

        return {
            'success': True,
            'motifs_found': len(unique_motifs),
            'motifs': unique_motifs,
            'all_motifs': motifs,  # 保留原始列表用于参考
            'output_dir': str(output_path)
        }
    except subprocess.TimeoutExpired:
        logger.error("MEME运行超时(10分钟)")
        return {'success': False, 'error': 'MEME运行超时'}
    except subprocess.CalledProcessError as e:
        logger.error(f"MEME运行失败: {e.stderr[:300] if e.stderr else ''}")
        return {'success': False, 'error': str(e)}


def parse_meme_result(meme_txt: Path) -> List[Dict[str, Any]]:
    """解析MEME文本结果"""
    motifs = []
    if not meme_txt.exists():
        return motifs

    with open(meme_txt, 'r') as f:
        content = f.read()

    for line in content.split('\n'):
        if line.startswith('MOTIF'):
            parts = line.split()
            if len(parts) >= 3:
                motif = {'id': parts[1]}
                # 尝试提取宽度
                for i, p in enumerate(parts):
                    if p == 'width' and i + 1 < len(parts):
                        try:
                            motif['width'] = int(parts[i + 1])
                        except:
                            pass
                    elif p == 'E-value' and i + 1 < len(parts):
                        try:
                            motif['evalue'] = float(parts[i + 1])
                        except:
                            pass
                motifs.append(motif)

    return motifs


def reverse_complement(seq: str) -> str:
    """返回序列的反向互补"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                  'U': 'A', 'R': 'Y', 'Y': 'R', 'S': 'S',
                  'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V',
                  'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N'}
    return ''.join(complement.get(b, b) for b in reversed(seq.upper()))


def normalize_motif(seq: str) -> str:
    """规范化motif序列（取正向和反向互补中字典序较小的）"""
    rc = reverse_complement(seq)
    return min(seq.upper(), rc.upper())


def deduplicate_motifs(motifs: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    合并重复的motif
    规则：
    1. 规范化motif序列（正向/反向互补取字典序较小的）
    2. 按E-value排序
    3. 相同规范化序列只保留E-value最好的一个
    4. 不限制数量，保留所有唯一motif
    """
    # 按E-value排序
    sorted_motifs = sorted(motifs, key=lambda x: x.get('evalue', float('inf')))

    seen = {}  # {normalized_seq: best_motif}
    for motif in sorted_motifs:
        # 提取序列（motif id格式：SEQUENCE MEME-N）
        seq = motif.get('id', '').split()[0] if ' ' in motif.get('id', '') else motif.get('id', '')
        if not seq:
            continue

        norm_seq = normalize_motif(seq)
        if norm_seq not in seen:
            motif['normalized_seq'] = norm_seq
            motif['original_seq'] = seq
            seen[norm_seq] = motif

    unique_motifs = list(seen.values())
    logger.info(f"Motif去重：原始{len(motifs)}个 → 去重后{len(unique_motifs)}个")
    return unique_motifs


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


def ensure_motif_database(database: Optional[str] = None) -> Optional[str]:
    """确保motif数据库可用，优先使用本地JASPAR数据库，必要时自动下载"""
    import urllib.request

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

    # 检查用户指定的数据库
    if database and Path(database).exists():
        return database

    # 自动下载JASPAR数据库
    logger.info("未找到JASPAR数据库，开始自动下载...")
    db_dir = Path("references/motif_databases/JASPAR")
    db_dir.mkdir(parents=True, exist_ok=True)

    # JASPAR MEME格式数据库下载地址
    # 优先尝试JASPAR2024vertebrates
    urls = [
        "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.25.tgz",
        "https://meme.nbwr.com/meme/meme-software/Databases/motifs/motif_databases.12.23.tgz",
    ]

    tgz_path = db_dir.parent / "motif_databases.tgz"

    for url in urls:
        try:
            logger.info(f"下载motif数据库: {url}")
            urllib.request.urlretrieve(url, tgz_path)
            logger.info("下载完成，解压中...")

            import tarfile
            with tarfile.open(tgz_path, "r:gz") as tar:
                tar.extractall(path=db_dir.parent)
            tgz_path.unlink()

            # 查找JASPAR文件
            for jaspar_db in db_dir.glob("JASPAR*_CORE_vertebrates_non-redundant.meme"):
                logger.info(f"数据库就绪: {jaspar_db}")
                return str(jaspar_db)

            break
        except Exception as e:
            logger.warning(f"下载失败 ({url}): {e}")
            if tgz_path.exists():
                tgz_path.unlink()
            continue

    # 备选：直接下载JASPAR MEME文件（如果完整包太大）
    meme_url = "https://jaspar2024.genereg.net/static/latest/JASPAR2024_CORE_vertebrates_non-redundant.meme.txt"
    meme_path = db_dir / "JASPAR2024_CORE_vertebrates_non-redundant.meme"

    if not meme_path.exists():
        try:
            logger.info(f"尝试直接下载JASPAR: {meme_url}")
            urllib.request.urlretrieve(meme_url, meme_path)
            logger.info(f"数据库已保存: {meme_path}")
            return str(meme_path)
        except Exception as e:
            logger.warning(f"直接下载失败: {e}")

    logger.warning("未找到JASPAR数据库，TomTom分析将跳过")
    return None


def run_tomtom_analysis(motif_file: str, output_dir: str,
                       database: str = None,
                       evalue_threshold: float = 0.05,
                       min_overlap: int = 5) -> Dict[str, Any]:
    """运行TomTom motif比较分析"""
    # 检查输入文件
    if not Path(motif_file).exists():
        logger.error(f"输入文件不存在: {motif_file}")
        return {'success': False, 'error': '输入文件不存在'}

    # 检查TomTom是否安装
    if not check_tomtom_installed():
        logger.error("TomTom未安装或不在PATH中，跳过TomTom分析")
        return {'success': False, 'error': 'TomTom未安装'}

    # 创建输出目录
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # 确保数据库可用
    database = ensure_motif_database(database)
    if database is None:
        logger.warning("未找到motif数据库，TomTom分析跳过")
        return {'success': False, 'error': 'motif数据库未找到'}

    logger.info(f"开始TomTom motif比较分析")
    logger.info(f"输入motif文件: {motif_file}")
    logger.info(f"数据库: {database}")

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

        logger.info(f"运行TomTom: {' '.join(cmd)}")

        # 执行TomTom命令
        process = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True,
            timeout=300
        )

        logger.info("TomTom分析完成")

        # 解析主要结果
        txt_file = output_path / 'tomtom.txt'
        comparisons = parse_tomtom_result(txt_file) if txt_file.exists() else []

        result = {
            'success': True,
            'output_dir': str(output_path),
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

    except subprocess.TimeoutExpired:
        logger.error("TomTom运行超时(5分钟)")
        return {'success': False, 'error': 'TomTom运行超时'}
    except subprocess.CalledProcessError as e:
        logger.error(f"TomTom执行失败: {e.stderr[:300] if e.stderr else ''}")
        return {'success': False, 'error': f'TomTom执行失败'}
    except Exception as e:
        logger.error(f"运行TomTom时出错: {e}")
        return {'success': False, 'error': str(e)}


def parse_tomtom_result(txt_file: Path) -> List[Dict[str, Any]]:
    """解析TomTom文本格式结果"""
    comparisons = []

    if not txt_file.exists():
        return comparisons

    try:
        with open(txt_file, 'r') as f:
            lines = f.readlines()

        # 查找标题行
        header_found = False
        headers = []
        for line in lines:
            line = line.strip()

            if not line:
                continue

            if line.startswith('#Query_ID'):
                header_found = True
                headers = line[1:].split('\t')
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
        if len(df) > 0 and 'E-value' in df.columns:
            df_sorted = df.sort_values('E-value')
            top_file = output_dir / f"tomtom_top{top_n}.csv"
            df_sorted.head(top_n).to_csv(top_file, index=False)

    logger.info(f"TomTom结果摘要已保存: {summary_file}")


def run_small_rna_motif_analysis(config: Dict[str, Any]) -> Dict[str, Any]:
    """运行完整的small RNA motif分析流程"""
    motif_cfg = config.get('motif_analysis', {})
    meme_cfg = motif_cfg.get('meme', {})
    mirbase_cfg = motif_cfg.get('mirbase', {})
    ref_cfg = config.get('reference', {})

    min_len = mirbase_cfg.get('min_len', 18)
    max_len = mirbase_cfg.get('max_len', 35)
    max_misms = mirbase_cfg.get('max_mismatches', 1)
    threads = meme_cfg.get('threads', 4)

    mirbase_fasta = ref_cfg.get('mirbase_fasta', 'references/hsa.mature.fa')
    index_path = ref_cfg.get('mirbase_bowtie2_index', 'references/bowtie2_index/hsa_mirbase')

    results_dir = Path(config.get('directories', {}).get('results', 'results'))
    motif_results_dir = results_dir / 'small_rna_motif'
    motif_results_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: 确保miRBase序列存在
    mirbase_fasta = ensure_mirbase_fasta(config)

    # Step 2: 构建索引
    if not build_mirbase_index(mirbase_fasta, index_path, threads):
        return {'success': False, 'error': 'miRBase索引构建失败'}

    # Step 3: 获取所有样本
    metadata_file = config.get('samples', {}).get('metadata_file', 'data/metadata/sample_info.csv')
    if os.path.exists(metadata_file):
        df = pd.read_csv(metadata_file)
        samples = df[config['samples']['sample_column']].tolist()
    else:
        samples = ["GAO_1", "GAO_2", "GAO_3", "PAL_1", "PAL_2", "PAL_3"]

    # Step 4: 每个样本比对到miRBase
    all_reads = []
    sample_stats = {}

    for sample in samples:
        trimmed_fastq = get_sample_trimmed_fastq(sample, config)
        sam_file = motif_results_dir / f"{sample}_mirbase.sam"

        success, mapped = map_to_mirbase(
            trimmed_fastq, index_path, str(sam_file),
            min_len, max_len, max_misms, threads
        )
        sample_stats[sample] = {
            'success': success,
            'mapped_reads': mapped
        }

        # 提取reads
        if success and sam_file.exists():
            reads = extract_mirna_reads_direct(str(trimmed_fastq), str(sam_file), min_len, max_len)
            all_reads.extend(reads)
            logger.info(f"  {sample}: 提取{len(reads)}条reads")

    # Step 5: 保存合并的reads
    combined_fasta = motif_results_dir / 'mirna_reads.fasta'
    save_reads_to_fasta(all_reads, str(combined_fasta))

    # Step 6: 运行MEME
    meme_result = run_meme_on_small_rna(
        str(combined_fasta),
        str(motif_results_dir / 'meme_results'),
        width_min=meme_cfg.get('min_width', 5),
        width_max=meme_cfg.get('max_width', 8),
        max_motifs=meme_cfg.get('max_motifs', 3),
        evalue_threshold=meme_cfg.get('evalue_threshold', 1e-4),
        min_sites=meme_cfg.get('minsites', 10),
        max_sites=meme_cfg.get('maxsites', 100),
        searchsize=meme_cfg.get('searchsize', 100000)
    )

    # 保存结果摘要
    summary = {
        'success': meme_result.get('success', False),
        'samples': sample_stats,
        'total_unique_reads': len(set(all_reads)),
        'combined_fasta': str(combined_fasta),
        'meme_result': meme_result
    }

    # 保存主摘要JSON
    with open(motif_results_dir / 'small_rna_motif_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    # 保存MEME摘要JSON（Snakemake需要）
    meme_summary_file = motif_results_dir / 'meme_results' / 'meme_summary.json'
    if meme_summary_file.parent.exists() or meme_summary_file.parent == motif_results_dir / 'meme_results':
        meme_summary_file.parent.mkdir(parents=True, exist_ok=True)
        meme_summary_data = {
            'success': meme_result.get('success', False),
            'motifs_found': meme_result.get('motifs_found', 0),
            'motifs': meme_result.get('motifs', []),
            'output_dir': meme_result.get('output_dir', '')
        }
        with open(meme_summary_file, 'w') as f:
            json.dump(meme_summary_data, f, indent=2)

    # Step 7: 运行TomTom与已知motif数据库比对
    tomtom_cfg = motif_cfg.get('tomtom', {})
    if tomtom_cfg.get('enabled', True) and meme_result.get('success', False):
        tomtom_result = run_tomtom_analysis(
            motif_file=str(motif_results_dir / 'meme_results' / 'meme.xml'),
            output_dir=str(motif_results_dir / 'tomtom_results'),
            database=tomtom_cfg.get('database', None),
            evalue_threshold=tomtom_cfg.get('evalue_threshold', 0.05),
            min_overlap=tomtom_cfg.get('min_overlap', 5)
        )
        summary['tomtom_result'] = tomtom_result

        # 保存TomTom摘要JSON
        if tomtom_result.get('success', False):
            tomtom_summary_file = motif_results_dir / 'tomtom_results' / 'tomtom_summary.json'
            with open(tomtom_summary_file, 'w') as f:
                json.dump(tomtom_result, f, indent=2)

    logger.info(f"=== Small RNA Motif分析完成 ===")
    logger.info(f"总样本数: {len(samples)}")
    logger.info(f"总唯一reads: {len(set(all_reads))}")
    logger.info(f"Motif发现数: {meme_result.get('motifs_found', 0)}")
    logger.info(f"结果目录: {motif_results_dir}")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Small RNA motif分析 - 从测序数据发现miRNA上的富集motif"
    )
    parser.add_argument('--config', default='config/config.yaml',
                       help='配置文件路径')
    parser.add_argument('--sample', help='只分析指定样本')

    args = parser.parse_args()

    # 加载配置
    if not Path(args.config).exists():
        logger.error(f"配置文件不存在: {args.config}")
        return 1

    with open(args.config, 'r') as f:
        config = yaml.safe_load(f)

    logger.info("=== Small RNA Motif分析开始 ===")

    result = run_small_rna_motif_analysis(config)

    if result.get('success'):
        logger.info("分析成功完成!")
        return 0
    else:
        logger.error(f"分析失败: {result.get('error', '未知错误')}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
