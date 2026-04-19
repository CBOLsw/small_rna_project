#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
统一日志配置工具
用于为整个项目提供一致的日志格式和配置
"""

import logging
import sys
import os
from pathlib import Path
from typing import Optional, List

# 全局统一日志文件路径
GLOBAL_LOG_FILE = None
# 保存所有已创建的日志记录器
_LOGGERS: List[logging.Logger] = []


def configure_logging(
    name: str = None,
    level: int = logging.INFO,
    log_file: Optional[str] = None,
    show_process: bool = True,
    show_thread: bool = True,
    show_file: bool = False,
    use_global_log: bool = True
) -> logging.Logger:
    """
    配置统一的日志格式

    参数:
        name: 日志记录器名称
        level: 日志级别
        log_file: 可选的日志文件路径
        show_process: 是否显示进程ID
        show_thread: 是否显示线程ID
        show_file: 是否显示文件名和行号
        use_global_log: 是否同时输出到全局统一日志文件

    返回:
        配置好的日志记录器
    """
    global GLOBAL_LOG_FILE

    # 构建日志格式
    format_parts = [
        "%(asctime)s",
        "%(levelname)-8s"
    ]

    if name:
        format_parts.append(f"[{name}]")

    if show_process:
        format_parts.append("%(process)d")

    if show_thread:
        format_parts.append("%(thread)d")

    if show_file:
        format_parts.append("%(filename)s:%(lineno)d")

    format_parts.append("%(message)s")

    log_format = " - ".join(format_parts)

    # 配置基础日志
    logging.basicConfig(
        level=level,
        format=log_format,
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    # 检查根记录器是否已有处理器
    root_logger = logging.getLogger()
    has_root_handlers = len(root_logger.handlers) > 0

    # 如果需要全局统一日志
    if use_global_log and GLOBAL_LOG_FILE is not None:
        log_dir = Path(GLOBAL_LOG_FILE).parent
        log_dir.mkdir(parents=True, exist_ok=True)

        # 检查是否已有全局日志文件处理器
        has_global_handler = False
        for handler in root_logger.handlers:
            if (isinstance(handler, logging.FileHandler) and
                handler.baseFilename == str(Path(GLOBAL_LOG_FILE).resolve())):
                has_global_handler = True
                break

        if not has_global_handler:
            # 添加全局统一日志文件处理器
            global_file_handler = logging.FileHandler(GLOBAL_LOG_FILE, encoding='utf-8')
            global_file_handler.setLevel(level)
            global_file_formatter = logging.Formatter(log_format, datefmt="%Y-%m-%d %H:%M:%S")
            global_file_handler.setFormatter(global_file_formatter)
            root_logger.addHandler(global_file_handler)

    # 如果指定了单独的日志文件，添加文件处理器
    if log_file and log_file != GLOBAL_LOG_FILE:
        log_dir = Path(log_file).parent
        log_dir.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(log_file, encoding='utf-8')
        file_handler.setLevel(level)
        file_formatter = logging.Formatter(log_format, datefmt="%Y-%m-%d %H:%M:%S")
        file_handler.setFormatter(file_formatter)

        root_logger.addHandler(file_handler)

    # 创建或获取日志记录器
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # 确保传播到根记录器
    logger.propagate = True

    return logger


def set_global_log_file(log_file: str):
    """
    设置全局统一日志文件

    参数:
        log_file: 全局日志文件路径
    """
    global GLOBAL_LOG_FILE
    GLOBAL_LOG_FILE = log_file

    # 为所有已创建的日志记录器添加全局日志文件处理器
    log_dir = Path(GLOBAL_LOG_FILE).parent
    log_dir.mkdir(parents=True, exist_ok=True)

    # 创建统一的日志格式
    log_format = ("%(asctime)s - %(levelname)-8s - [%(name)s] - "
                  "%(process)d - %(thread)d - %(message)s")

    # 创建文件处理器
    global_file_handler = logging.FileHandler(GLOBAL_LOG_FILE, encoding='utf-8')
    global_file_handler.setLevel(logging.INFO)
    global_file_formatter = logging.Formatter(
        log_format,
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    global_file_handler.setFormatter(global_file_formatter)

    # 为所有已保存的日志记录器添加处理器
    for logger in _LOGGERS:
        # 检查是否已经有类似的处理器
        has_global_handler = False
        for handler in logger.handlers:
            if (isinstance(handler, logging.FileHandler) and
                handler.baseFilename == str(Path(GLOBAL_LOG_FILE).resolve())):
                has_global_handler = True
                break

        if not has_global_handler:
            logger.addHandler(global_file_handler)
            logger.debug(f"已添加全局日志文件处理器到: {logger.name}")


def get_logger(name: str, **kwargs) -> logging.Logger:
    """
    获取配置好的日志记录器

    参数:
        name: 日志记录器名称
        **kwargs: 传递给 configure_logging 的其他参数

    返回:
        配置好的日志记录器
    """
    logger = configure_logging(name, **kwargs)
    if logger not in _LOGGERS:
        _LOGGERS.append(logger)
    return logger


def get_script_logger(script_name: str, log_dir: str = None, **kwargs) -> logging.Logger:
    """
    获取脚本专用的日志记录器

    参数:
        script_name: 脚本名称（用于生成日志文件名）
        log_dir: 日志文件保存目录
        **kwargs: 传递给 configure_logging 的其他参数

    返回:
        配置好的日志记录器
    """
    log_file = None
    if log_dir:
        log_file = os.path.join(log_dir, f"{script_name}.log")

    return get_logger(script_name, log_file=log_file, **kwargs)


# 测试用例
if __name__ == "__main__":
    # 测试1: 设置全局日志
    set_global_log_file("pipeline.log")

    # 测试2: 基础配置
    logger1 = get_logger("TestLogger")
    logger1.debug("调试信息")
    logger1.info("普通信息")
    logger1.warning("警告信息")
    logger1.error("错误信息")
    logger1.critical("严重错误")

    # 测试3: 包含文件日志
    logger2 = get_logger("FileLogger", log_file="test.log")
    logger2.info("这是一条写入文件的日志")

    # 测试4: 脚本专用配置
    logger3 = get_script_logger("my_script", log_dir="logs")
    logger3.info("脚本运行信息")

    print("日志配置测试完成")
