#!/usr/bin/env python3
"""
Small RNA测序分析流程执行脚本
支持命令行参数，提供完整的分析流程管理功能

功能：
1. 支持完整流程或单个模块运行
2. 提供进度跟踪和日志记录
3. 支持错误处理和恢复机制
4. 提供流程状态查询功能

使用方法：
    python scripts/run_pipeline.py --config config/config.yaml
    python scripts/run_pipeline.py --config config/config.yaml --module qc
    python scripts/run_pipeline.py --config config/config.yaml --resume
"""

import os
import sys
import argparse
import subprocess
import logging
import json
import yaml
import time
from pathlib import Path
from typing import List, Dict, Optional, Any
from datetime import datetime

# 添加项目根目录到Python路径
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# 尝试导入tqdm，如果失败则禁用进度条功能
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False

# 配置日志
def setup_logging(log_file: str = None, verbose: bool = False) -> None:
    """配置日志系统"""
    log_level = logging.DEBUG if verbose else logging.INFO

    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )


def load_config(config_file: str) -> Dict[str, Any]:
    """加载配置文件"""
    logger = logging.getLogger(__name__)
    try:
        with open(config_file, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        logger.info(f"配置文件已加载: {config_file}")
        return config
    except Exception as e:
        logger.error(f"加载配置文件失败: {e}")
        sys.exit(1)


def create_directories(config: Dict[str, Any]) -> None:
    """创建必要的目录"""
    logger = logging.getLogger(__name__)
    directories = [
        config['directories']['processed'],
        config['directories']['results'],
        config['directories']['logs'],
        os.path.join(config['directories']['results'], 'qc'),
        os.path.join(config['directories']['results'], 'alignment'),
        os.path.join(config['directories']['results'], 'counts'),
        os.path.join(config['directories']['results'], 'differential_expression'),
        os.path.join(config['directories']['results'], 'small_rna_motif'),
    ]

    for directory in directories:
        Path(directory).mkdir(parents=True, exist_ok=True)
        logger.debug(f"目录已确认: {directory}")


def run_snakemake(config: Dict[str, Any],
                  target: str = None,
                  dry_run: bool = False,
                  resume: bool = False,
                  cores: int = None,
                  verbose: bool = False,
                  show_progress: bool = True) -> int:
    """
    运行Snakemake流程

    参数:
        config: 配置字典
        target: 目标规则（可选，用于运行特定模块）
        dry_run: 是否执行 dry-run
        resume: 是否从上次失败处恢复
        cores: 使用的核心数
        verbose: 是否详细输出
        show_progress: 是否显示进度条（保留参数，现在总是直接显示Snakemake输出）

    返回:
        int: Snakemake退出码
    """
    logger = logging.getLogger(__name__)

    # 设置全局统一日志文件
    log_dir = config['directories']['logs']
    pipeline_log_file = os.path.join(log_dir, 'pipeline.log')
    from utils.logging_utils import set_global_log_file
    set_global_log_file(pipeline_log_file)
    logger.info(f"全局日志文件已设置为: {pipeline_log_file}")

    # 首先尝试解锁（如果有锁的话）
    unlock_cmd = [
        'snakemake',
        '--snakefile', 'workflow/Snakefile',
        '--configfile', 'config/config.yaml',
        '--unlock'
    ]
    logger.info("检查并清理Snakemake锁...")
    subprocess.run(unlock_cmd, capture_output=True, check=False)

    # 确定核心数
    if cores is None:
        cores = config['snakemake']['cores']

    # 构建Snakemake命令，添加--printshellcmds来显示执行的命令
    cmd = [
        'snakemake',
        '--snakefile', 'workflow/Snakefile',
        '--configfile', 'config/config.yaml',
        '--cores', str(cores),
        '--latency-wait', str(config['snakemake']['latency_wait']),
        '--restart-times', str(config['snakemake']['restart_times']),
        '--printshellcmds',  # 显示正在执行的shell命令
        '--rerun-triggers', 'mtime',  # 仅根据文件修改时间判断是否重新运行
    ]

    if dry_run:
        cmd.append('--dry-run')

    if resume:
        cmd.append('--rerun-incomplete')
        cmd.append('--keep-going')
    elif config['snakemake']['keep_going']:
        cmd.append('--keep-going')

    if verbose:
        cmd.append('--verbose')

    # 添加目标规则
    if target:
        cmd.append(target)

    logger.info(f"运行Snakemake命令: {' '.join(cmd)}")
    logger.info(f"开始执行时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 60)
    logger.info("分析开始，Snakemake将显示实时进度...")
    logger.info("=" * 60)

    try:
        start_time = time.time()

        # 直接运行Snakemake，让它自己显示进度
        result = subprocess.run(cmd, check=False)

        end_time = time.time()

        elapsed_time = end_time - start_time
        hours, remainder = divmod(elapsed_time, 3600)
        minutes, seconds = divmod(remainder, 60)

        logger.info("=" * 60)
        logger.info(f"执行完成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info(f"总执行时间: {int(hours)}小时 {int(minutes)}分钟 {int(seconds)}秒")

        if result.returncode == 0:
            logger.info("Snakemake执行成功完成！")
        else:
            logger.error(f"Snakemake执行失败，退出码: {result.returncode}")

        return result.returncode

    except KeyboardInterrupt:
        logger.warning("执行被用户中断")
        return 130
    except Exception as e:
        logger.error(f"执行Snakemake时出错: {e}")
        return 1


def check_pipeline_status(config: Dict[str, Any]) -> None:
    """检查流程状态"""
    logger = logging.getLogger(__name__)

    results_dir = Path(config['directories']['results'])
    log_dir = Path(config['directories']['logs'])
    if not results_dir.exists():
        logger.info("结果目录不存在，流程尚未开始")
        return

    # 检查各模块输出
    modules = {
        'qc': {
            'output': os.path.join(results_dir, 'qc', 'qc_summary.csv'),
            'log_prefix': 'fastqc_'
        },
        'alignment': {
            'output': os.path.join(results_dir, 'alignment', 'alignment_summary.csv'),
            'log_prefix': 'bowtie2_'
        },
        'counts': {
            'output': os.path.join(results_dir, 'counts', 'counts_summary.csv'),
            'log_prefix': 'count_features_'
        },
        'differential_expression': {
            'output': os.path.join(results_dir, 'differential_expression', 'deseq2_results.csv'),
            'log_prefix': 'deseq2_analysis'
        },
        'motif_analysis': {
            'output': os.path.join(results_dir, 'small_rna_motif', 'meme_results', 'meme_summary.json'),
            'log_prefix': 'meme_'
        },
    }

    logger.info("流程状态检查:")
    completed_modules = []
    running_modules = []
    pending_modules = []

    for module_name, module_info in modules.items():
        output_file = module_info['output']
        log_prefix = module_info['log_prefix']

        if Path(output_file).exists():
            mtime = datetime.fromtimestamp(Path(output_file).stat().st_mtime)
            logger.info(f"  [✓] {module_name:25s} - 完成于 {mtime.strftime('%Y-%m-%d %H:%M:%S')}")
            completed_modules.append(module_name)
        else:
            # 检查是否有正在运行的日志文件
            is_running = False
            if log_dir.exists():
                for log_file in log_dir.glob(f"{log_prefix}*.log"):
                    try:
                        # 检查日志文件的修改时间是否在最近5分钟内
                        mtime = datetime.fromtimestamp(log_file.stat().st_mtime)
                        time_diff = (datetime.now() - mtime).total_seconds()
                        if time_diff < 300:  # 5分钟内更新过，认为正在运行
                            is_running = True
                            break
                    except:
                        pass

            if is_running:
                logger.info(f"  [⚡] {module_name:25s} - 执行中...")
                running_modules.append(module_name)
            else:
                logger.info(f"  [ ] {module_name:25s} - 待执行")
                pending_modules.append(module_name)

    logger.info(f"\n已完成: {len(completed_modules)}/{len(modules)} 个模块")
    if completed_modules:
        logger.info(f"已完成模块: {', '.join(completed_modules)}")
    if running_modules:
        logger.info(f"执行中模块: {', '.join(running_modules)}")
    if pending_modules:
        logger.info(f"待完成模块: {', '.join(pending_modules)}")


def get_available_modules() -> Dict[str, str]:
    """获取可用的模块列表"""
    return {
        'all': '运行完整分析流程',
        'check': '仅运行项目状态检查',
        'qc': '仅运行数据质量控制模块',
        'alignment': '仅运行序列比对模块',
        'counts': '仅运行基因计数模块',
        'de': '仅运行差异表达分析模块',
        'motif': '仅运行motif分析模块（完整流程）',
        'mirbase': '仅运行miRBase比对步骤',
        'meme': '仅运行MEME分析步骤',
    }


def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description="Small RNA测序分析流程执行脚本",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例:
  # 运行完整流程
  python scripts/run_pipeline.py --config config/config.yaml

  # 运行项目状态检查
  python scripts/run_pipeline.py --config config/config.yaml --check
  python scripts/run_pipeline.py --config config/config.yaml --module check

  # 仅运行质量控制模块
  python scripts/run_pipeline.py --config config/config.yaml --module qc

  # 从上次失败处恢复
  python scripts/run_pipeline.py --config config/config.yaml --resume

  # 检查流程状态
  python scripts/run_pipeline.py --config config/config.yaml --status

  # Dry-run (查看将要执行的步骤)
  python scripts/run_pipeline.py --config config/config.yaml --dry-run
        """
    )

    # 必需参数
    parser.add_argument('--config', required=True,
                       help='配置文件路径 (config/config.yaml)')

    # 运行模式
    parser.add_argument('--module', choices=get_available_modules().keys(),
                       help='指定运行的模块 (默认: all)')
    parser.add_argument('--resume', action='store_true',
                       help='从上次失败处恢复执行')
    parser.add_argument('--dry-run', action='store_true',
                       help='仅显示将要执行的步骤，不实际运行')

    # 状态查询
    parser.add_argument('--status', action='store_true',
                       help='检查当前流程状态')
    parser.add_argument('--check', action='store_true',
                       help='运行项目状态检查')

    # 执行参数
    parser.add_argument('--cores', type=int,
                       help='使用的核心数 (默认: 配置文件中的设置)')
    parser.add_argument('--log-file',
                       help='日志文件路径 (默认: 标准输出)')

    # 输出参数
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='显示详细输出')
    parser.add_argument('--list-modules', action='store_true',
                       help='列出所有可用模块')
    parser.add_argument('--no-progress', action='store_true',
                       help='禁用进度条显示')

    args = parser.parse_args()

    # 处理列出模块请求
    if args.list_modules:
        print("可用模块:")
        for module_name, description in get_available_modules().items():
            print(f"  {module_name:12s} - {description}")
        return 0

    # 设置日志
    setup_logging(args.log_file, args.verbose)
    logger = logging.getLogger(__name__)

    logger.info("=" * 60)
    logger.info("Small RNA测序分析流程")
    logger.info("=" * 60)

    # 加载配置
    config = load_config(args.config)

    # 创建必要的目录
    create_directories(config)

    # 处理状态检查请求
    if args.status:
        check_pipeline_status(config)
        return 0

    # 处理项目状态检查请求
    if args.check:
        logger.info("运行项目状态检查")
        check_command = ['python', 'scripts/utils/final_check.py']
        result = subprocess.run(check_command, capture_output=True, text=True, encoding='utf-8', errors='replace')
        if result.returncode == 0:
            logger.info("项目状态检查通过")
            if result.stdout:
                sys.stdout.buffer.write(result.stdout.encode('utf-8', 'replace'))
                print()
        else:
            logger.error("项目状态检查失败")
            if result.stderr:
                sys.stderr.buffer.write(result.stderr.encode('utf-8', 'replace'))
                print()
        return result.returncode

    # 确定目标
    target_map = {
        'check': os.path.join(config['directories']['logs'], '.project_check_complete'),
        'qc': os.path.join(config['directories']['results'], 'qc', 'qc_summary.csv'),
        'alignment': os.path.join(config['directories']['results'], 'alignment', 'alignment_summary.csv'),
        'counts': os.path.join(config['directories']['results'], 'counts', 'counts_summary.csv'),
        'de': os.path.join(config['directories']['results'], 'differential_expression', 'deseq2_results.csv'),
        'motif': os.path.join(config['directories']['results'], 'small_rna_motif', 'meme_results', 'meme_summary.json'),
        # motif子步骤
        'mirbase': os.path.join(config['directories']['results'], 'small_rna_motif', 'mirna_reads.fasta'),
        'meme': os.path.join(config['directories']['results'], 'small_rna_motif', 'meme_results', 'meme_summary.json'),
    }

    target = target_map.get(args.module) if args.module and args.module != 'all' else None

    # 运行流程
    logger.info(f"项目名称: {config['project_name']}")
    logger.info(f"运行模式: {'dry-run' if args.dry_run else '实际执行'}")
    if args.module:
        logger.info(f"指定模块: {args.module}")
    if args.resume:
        logger.info("恢复模式: 开启")

    return_code = run_snakemake(
        config=config,
        target=target,
        dry_run=args.dry_run,
        resume=args.resume,
        cores=args.cores,
        verbose=args.verbose,
        show_progress=not args.no_progress and TQDM_AVAILABLE
    )

    logger.info("=" * 60)
    if return_code == 0:
        logger.info("流程执行成功！")
    else:
        logger.error(f"流程执行失败 (退出码: {return_code})")
    logger.info("=" * 60)

    return return_code


if __name__ == '__main__':
    sys.exit(main())
