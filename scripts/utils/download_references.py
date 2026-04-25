#!/usr/bin/env python3
"""
Python参考基因组下载脚本
跨平台下载hg38参考基因组和注释文件
使用国内镜像源加速下载
"""

import os
import sys
import gzip
import shutil
from pathlib import Path
from urllib import request
from urllib.error import URLError, HTTPError
import logging
from typing import Optional

# 配置日志 - 简化格式，避免Windows编码问题
logging.basicConfig(
    level=logging.INFO,
    format='%(message)s'
)
logger = logging.getLogger(__name__)

# 修复Windows编码问题
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')


class ReferenceDownloader:
    """参考基因组下载器"""

    def __init__(self, ref_dir: str = "references"):
        self.ref_dir = Path(ref_dir)
        self.ref_dir.mkdir(parents=True, exist_ok=True)

    def print_header(self, title: str):
        """打印标题"""
        logger.info("")
        logger.info("+" + "-" * 58 + "+")
        logger.info(f"|  {title:^52}  |")
        logger.info("+" + "-" * 58 + "+")

    def download_file(self, url: str, dest_path: str, description: str = "文件") -> Optional[str]:
        """下载文件"""
        dest_path = Path(dest_path)

        # 1. 检查解压后的文件是否已存在
        if dest_path.exists():
            logger.info(f"  [OK] {description}已存在")
            return str(dest_path)

        # 2. 检查对应的压缩包是否已存在（避免重复下载）
        gz_path = Path(str(dest_path) + '.gz')
        if gz_path.exists():
            logger.info(f"  [OK] {description}压缩包已存在，跳过下载")
            logger.info(f"  正在解压...")
            try:
                with gzip.open(gz_path, 'rb') as f_in:
                    with open(dest_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                if dest_path.exists():
                    size_mb = dest_path.stat().st_size / (1024 * 1024)
                    logger.info(f"  解压完成！文件大小: {size_mb:.1f} MB")
                    return str(dest_path)
            except Exception as e:
                logger.warning(f"  解压失败: {e}，将重新下载")
                gz_path = None  # 标记需要重新下载

        logger.info(f"  正在下载{description}...")

        try:
            # 显示进度条的下载函数
            def report_progress(block_num, block_size, total_size):
                if total_size > 0:
                    percent = min(100, (block_num * block_size * 100) / total_size)
                    bar_length = 40
                    filled = int(percent / 100 * bar_length)
                    bar = "=" * filled + "-" * (bar_length - filled)
                    sys.stdout.write(f"\r  [{bar}] {percent:.1f}%")
                    sys.stdout.flush()

            # 下载文件
            temp_path = dest_path.with_suffix(dest_path.suffix + '.tmp')
            request.urlretrieve(url, temp_path, reporthook=report_progress)

            sys.stdout.write("\r  [OK] 下载完成！\n")
            sys.stdout.flush()

            # 如果是gz文件，解压
            if dest_path.suffix == '.fa' or dest_path.suffix == '.gtf':
                # 检查是否是压缩文件
                if temp_path.suffix == '.gz':
                    logger.info("  正在解压文件...")
                    with gzip.open(temp_path, 'rb') as f_in:
                        with open(dest_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                    # 删除临时压缩文件
                    temp_path.unlink()
                else:
                    # 直接重命名
                    temp_path.rename(dest_path)
            else:
                temp_path.rename(dest_path)

            # 显示文件大小
            if dest_path.exists():
                size_mb = dest_path.stat().st_size / (1024 * 1024)
                logger.info(f"  文件大小: {size_mb:.1f} MB")

            return str(dest_path)

        except HTTPError as e:
            logger.error(f"  [ERROR] HTTP错误 ({e.code}): {e.reason}")
            if 'temp_path' in locals() and temp_path.exists():
                try:
                    temp_path.unlink()
                except Exception as unlink_error:
                    logger.warning(f"  [WARNING] 无法删除临时文件: {unlink_error}")
            return None
        except URLError as e:
            logger.error(f"  [ERROR] URL错误: {e.reason}")
            if 'temp_path' in locals() and temp_path.exists():
                try:
                    temp_path.unlink()
                except Exception as unlink_error:
                    logger.warning(f"  [WARNING] 无法删除临时文件: {unlink_error}")
            return None
        except Exception as e:
            logger.error(f"  [ERROR] 下载失败: {e}")
            if 'temp_path' in locals() and temp_path.exists():
                try:
                    temp_path.unlink()
                except Exception as unlink_error:
                    logger.warning(f"  [WARNING] 无法删除临时文件: {unlink_error}")
            return None

    def download_hg38_genome(self) -> Optional[str]:
        """下载hg38参考基因组 - 使用UCSC官方源"""
        # 提示：建议手动下载（https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz）
        # 放置到 references/hg38.fa.gz，脚本会自动检测并跳过下载
        url_candidates = [
            # 源1：UCSC官方HTTPS（稳定可用）
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
        ]

        dest_path = self.ref_dir / "hg38.fa"

        # 尝试每个URL
        for i, url in enumerate(url_candidates):
            if i > 0:
                logger.info(f"  尝试备用源 ({i+1}/{len(url_candidates)})")
            result = self.download_file(url, dest_path, "hg38参考基因组")
            if result is not None:
                return result

        return None

    def download_hg38_gtf(self) -> Optional[str]:
        """下载hg38基因注释 - 使用UCSC官方源"""
        # 提示：建议手动下载（https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz）
        # 放置到 references/hg38.knownGene.gtf.gz，脚本会自动检测并跳过下载
        url_candidates = [
            # 源1：UCSC官方HTTPS（稳定可用）
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.knownGene.gtf.gz",
        ]

        dest_path = self.ref_dir / "hg38.knownGene.gtf"

        # 尝试每个URL
        for i, url in enumerate(url_candidates):
            if i > 0:
                logger.info(f"  尝试备用源 ({i+1}/{len(url_candidates)})")
            result = self.download_file(url, dest_path, "hg38基因注释")
            if result is not None:
                return result

        return None

    def download_small_rna_annotation(self) -> Optional[str]:
        """下载small RNA注释"""
        # 多个备选URL，按优先级尝试
        url_candidates = [
            "https://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3",
            "http://www.mirbase.org/ftp/CURRENT/genomes/hsa.gff3",
            "ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3",
        ]

        dest_path = self.ref_dir / "hg38.mirbase.gff3"

        # 如果文件已存在，直接返回
        if dest_path.exists():
            logger.info("  [OK] miRBase small RNA注释已存在")
            return str(dest_path)

        # 尝试每个URL
        for i, url in enumerate(url_candidates):
            logger.info(f"  尝试下载miRBase注释 (源{i+1}/{len(url_candidates)})")
            result = self.download_file(url, dest_path, "miRBase small RNA注释")
            if result is not None:
                logger.info("  [OK] 成功下载miRBase注释")
                return result

        # 所有URL都失败
        logger.warning("  [WARNING] 所有miRBase下载源都失败，small RNA注释将不可用")
        return None

    def prepare_index_directory(self) -> Path:
        """准备Bowtie2索引目录"""
        index_dir = self.ref_dir / "bowtie2_index"
        index_dir.mkdir(parents=True, exist_ok=True)
        return index_dir

    def download_all(self) -> dict:
        """下载所有参考数据"""
        results = {}

        self.print_header("开始下载参考基因组数据 (hg38)")

        # 下载hg38基因组
        genome_file = self.download_hg38_genome()
        if genome_file:
            results['genome'] = genome_file

        # 下载基因注释
        gtf_file = self.download_hg38_gtf()
        if gtf_file:
            results['gtf'] = gtf_file

        # 下载small RNA注释
        mirna_file = self.download_small_rna_annotation()
        if mirna_file:
            results['mirna'] = mirna_file

        # 准备索引目录
        index_dir = self.prepare_index_directory()
        results['index_dir'] = str(index_dir)

        self.print_header("下载完成!")

        # 显示文件列表
        logger.info("")
        logger.info("文件统计:")
        for file_type, file_path in results.items():
            if file_type in ['genome', 'gtf', 'mirna']:
                file_path = Path(file_path)
                if file_path.exists():
                    size_mb = file_path.stat().st_size / (1024 * 1024)
                    logger.info(f"  [OK] {file_type}: {file_path.name} ({size_mb:.1f} MB)")

        logger.info(f"\n索引目录: {index_dir}")

        return results


def main():
    """主函数"""
    import argparse

    parser = argparse.ArgumentParser(description="参考基因组下载脚本")
    parser.add_argument("--output", "-o", default="references",
                       help="输出目录 (默认: references)")
    parser.add_argument("--force", "-f", action="store_true",
                       help="强制重新下载，即使文件已存在")

    args = parser.parse_args()

    # 如果强制重新下载，删除现有文件
    if args.force:
        logger.warning("强制模式：将删除现有参考数据文件")
        ref_dir = Path(args.output)
        for file_pattern in ['*.fa', '*.gtf', '*.gff3']:
            for file in ref_dir.glob(file_pattern):
                logger.info(f"  删除文件: {file}")
                file.unlink()

    # 创建下载器并下载
    downloader = ReferenceDownloader(ref_dir=args.output)
    results = downloader.download_all()

    # 生成下一步提示
    if 'genome' in results:
        genome_file = results['genome']
        index_dir = results.get('index_dir', 'references/bowtie2_index')

        print("\n下一步:")
        print("  1. 构建Bowtie2索引:")
        print(f"     python scripts/alignment/build_bowtie2_index.py --genome {genome_file} --output {index_dir}/hg38")
        print("\n  2. 运行完整分析:")
        print("     python scripts/run_pipeline.py --config config/config.yaml --cores 8")

    return 0 if 'genome' in results and 'gtf' in results else 1


if __name__ == "__main__":
    sys.exit(main())
