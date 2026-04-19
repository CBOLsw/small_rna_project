#!/usr/bin/env python3
"""
压缩文件处理工具函数

功能：
1. 检测并处理压缩的FASTA/FASTQ文件
2. 支持gzip (.gz) 格式
3. 提供解压和压缩的辅助函数
"""

import os
import gzip
import shutil
from pathlib import Path
from typing import Optional, Tuple
import logging

logger = logging.getLogger(__name__)


def is_compressed(file_path: str) -> bool:
    """
    检查文件是否是压缩文件

    参数:
        file_path: 文件路径

    返回:
        bool: 是否是压缩文件
    """
    path = Path(file_path)
    if not path.exists():
        return False

    # 检查文件扩展名
    if path.suffix.lower() == '.gz':
        return True

    # 检查文件内容（魔数）
    try:
        with open(file_path, 'rb') as f:
            magic_number = f.read(2)
            if magic_number == b'\x1f\x8b':
                return True
    except Exception:
        pass

    return False


def get_uncompressed_path(file_path: str, output_dir: Optional[str] = None) -> str:
    """
    获取解压后的文件路径

    参数:
        file_path: 原文件路径
        output_dir: 输出目录（可选，默认与原文件同目录）

    返回:
        str: 解压后的文件路径
    """
    path = Path(file_path)

    if not is_compressed(file_path):
        return file_path

    # 生成解压后的文件名
    if path.suffix.lower() == '.gz':
        uncompressed_name = path.stem
    else:
        uncompressed_name = f"{path.name}.uncompressed"

    if output_dir:
        uncompressed_path = Path(output_dir) / uncompressed_name
    else:
        uncompressed_path = path.parent / uncompressed_name

    return str(uncompressed_path)


def decompress_file(file_path: str, output_path: Optional[str] = None,
                    keep_original: bool = True) -> Tuple[str, bool]:
    """
    解压gzip压缩的文件

    参数:
        file_path: 压缩文件路径
        output_path: 输出文件路径（可选）
        keep_original: 是否保留原文件（默认True）

    返回:
        Tuple[str, bool]: (解压后的文件路径, 是否成功)
    """
    if not is_compressed(file_path):
        logger.info(f"文件不是压缩格式，直接返回: {file_path}")
        return file_path, True

    path = Path(file_path)
    if not path.exists():
        logger.error(f"文件不存在: {file_path}")
        return "", False

    if output_path is None:
        output_path = get_uncompressed_path(file_path)

    output_path_obj = Path(output_path)

    # 如果解压后的文件已存在，检查是否需要重新解压
    if output_path_obj.exists():
        src_mtime = path.stat().st_mtime
        dst_mtime = output_path_obj.stat().st_mtime
        if dst_mtime >= src_mtime:
            logger.info(f"解压后的文件已存在且是最新的: {output_path}")
            return output_path, True
        else:
            logger.info(f"解压后的文件已过期，重新解压: {output_path}")

    logger.info(f"正在解压文件: {file_path} -> {output_path}")

    try:
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        logger.info(f"解压完成: {output_path}")
        return output_path, True

    except Exception as e:
        logger.error(f"解压失败: {e}")
        if output_path_obj.exists():
            output_path_obj.unlink()
        return "", False


def compress_file(file_path: str, output_path: Optional[str] = None,
                  keep_original: bool = True) -> Tuple[str, bool]:
    """
    压缩文件为gzip格式

    参数:
        file_path: 原文件路径
        output_path: 输出文件路径（可选）
        keep_original: 是否保留原文件（默认True）

    返回:
        Tuple[str, bool]: (压缩后的文件路径, 是否成功)
    """
    path = Path(file_path)
    if not path.exists():
        logger.error(f"文件不存在: {file_path}")
        return "", False

    if output_path is None:
        output_path = f"{file_path}.gz"

    output_path_obj = Path(output_path)

    logger.info(f"正在压缩文件: {file_path} -> {output_path}")

    try:
        with open(file_path, 'rb') as f_in:
            with gzip.open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        logger.info(f"压缩完成: {output_path}")
        return output_path, True

    except Exception as e:
        logger.error(f"压缩失败: {e}")
        if output_path_obj.exists():
            output_path_obj.unlink()
        return "", False


def ensure_uncompressed(file_path: str, output_dir: Optional[str] = None,
                        keep_original: bool = True) -> Tuple[str, bool]:
    """
    确保文件是解压状态，如果是压缩文件则解压

    参数:
        file_path: 文件路径
        output_dir: 输出目录（可选）
        keep_original: 是否保留原文件（默认True）

    返回:
        Tuple[str, bool]: (解压后的文件路径, 是否成功)
    """
    if not is_compressed(file_path):
        return file_path, True

    return decompress_file(file_path, None, keep_original)
