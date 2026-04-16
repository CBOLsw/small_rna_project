#!/usr/bin/env python3
"""
数据质控模块测试脚本

功能：
1. 测试FastQC分析脚本导入和基本功能
2. 测试Trimmomatic脚本导入和基本功能
3. 测试汇总脚本导入和基本功能
4. 验证模块完整性

使用方法：
    python test_qc.py [--test-data <测试数据目录>]
"""

import os
import sys
import argparse
import subprocess
import tempfile
from pathlib import Path
import logging

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class QCTester:
    """质量控制测试类"""

    def __init__(self, test_data_dir: str = None):
        self.test_data_dir = Path(test_data_dir) if test_data_dir else None
        self.test_results = {}

    def test_module_imports(self) -> bool:
        """测试模块导入"""
        logger.info("测试模块导入...")

        modules_to_test = [
            ('fastqc_analysis', 'FastQCAnalyzer'),
            ('trim_fastq', 'TrimmomaticProcessor'),
            ('qc_summary', 'QCSummary')
        ]

        all_passed = True
        for module_name, class_name in modules_to_test:
            try:
                module = __import__(module_name)
                if hasattr(module, class_name):
                    logger.info(f"  ✓ {module_name}: 成功导入 {class_name}")
                else:
                    logger.error(f"  ✗ {module_name}: 未找到类 {class_name}")
                    all_passed = False
            except ImportError as e:
                logger.error(f"  ✗ {module_name}: 导入失败 - {e}")
                all_passed = False

        self.test_results['module_imports'] = all_passed
        return all_passed

    def test_script_execution(self) -> bool:
        """测试脚本执行（不实际运行分析）"""
        logger.info("测试脚本执行...")

        scripts_to_test = [
            ('fastqc_analysis.py', ['--help']),
            ('trim_fastq.py', ['--help']),
            ('qc_summary.py', ['--help'])
        ]

        all_passed = True
        for script_name, args in scripts_to_test:
            script_path = Path(__file__).parent / script_name
            if not script_path.exists():
                logger.error(f"  ✗ {script_name}: 文件不存在")
                all_passed = False
                continue

            try:
                cmd = [sys.executable, str(script_path)] + args
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    timeout=10,
                    check=False
                )

                if result.returncode == 0:
                    logger.info(f"  ✓ {script_name}: 成功执行 (返回码: {result.returncode})")
                else:
                    logger.warning(f"  ⚠ {script_name}: 执行完成但返回非零码 (返回码: {result.returncode})")
                    logger.debug(f"    错误输出: {result.stderr[:200]}")

            except subprocess.TimeoutExpired:
                logger.error(f"  ✗ {script_name}: 执行超时")
                all_passed = False
            except Exception as e:
                logger.error(f"  ✗ {script_name}: 执行失败 - {e}")
                all_passed = False

        self.test_results['script_execution'] = all_passed
        return all_passed

    def test_fastqc_analyzer(self) -> bool:
        """测试FastQC分析器类"""
        logger.info("测试FastQC分析器类...")

        try:
            from fastqc_analysis import FastQCAnalyzer

            # 创建分析器实例
            analyzer = FastQCAnalyzer()

            # 测试方法存在性
            required_methods = ['check_fastqc', 'run_fastqc', 'generate_summary']
            for method in required_methods:
                if hasattr(analyzer, method):
                    logger.info(f"  ✓ FastQCAnalyzer 有方法: {method}")
                else:
                    logger.error(f"  ✗ FastQCAnalyzer 缺少方法: {method}")
                    return False

            # 测试文件收集逻辑（模拟）
            test_path = Path(__file__).parent.parent / "data"
            if test_path.exists():
                # 测试内部方法
                files = analyzer._collect_fastq_files(str(test_path))
                logger.info(f"  ✓ 文件收集测试: 找到 {len(files)} 个文件")

            self.test_results['fastqc_analyzer'] = True
            return True

        except Exception as e:
            logger.error(f"  ✗ FastQC分析器测试失败: {e}")
            self.test_results['fastqc_analyzer'] = False
            return False

    def test_trimmomatic_processor(self) -> bool:
        """测试Trimmomatic处理器类"""
        logger.info("测试Trimmomatic处理器类...")

        try:
            from trim_fastq import TrimmomaticProcessor

            # 创建处理器实例
            processor = TrimmomaticProcessor()

            # 测试方法存在性
            required_methods = ['check_trimmomatic', 'trim_single_end', 'trim_paired_end', 'generate_summary']
            for method in required_methods:
                if hasattr(processor, method):
                    logger.info(f"  ✓ TrimmomaticProcessor 有方法: {method}")
                else:
                    logger.error(f"  ✗ TrimmomaticProcessor 缺少方法: {method}")
                    return False

            # 测试适配器序列
            if hasattr(processor, 'ADAPTERS'):
                adapters = processor.ADAPTERS
                logger.info(f"  ✓ 适配器序列: {len(adapters)} 个")
                for name, seq in adapters.items():
                    logger.debug(f"    {name}: {seq}")

            self.test_results['trimmomatic_processor'] = True
            return True

        except Exception as e:
            logger.error(f"  ✗ Trimmomatic处理器测试失败: {e}")
            self.test_results['trimmomatic_processor'] = False
            return False

    def test_qc_summary(self) -> bool:
        """测试质量控制汇总类"""
        logger.info("测试质量控制汇总类...")

        try:
            from qc_summary import QCSummary

            # 创建汇总器实例
            qc_summary = QCSummary()

            # 测试方法存在性
            required_methods = ['load_fastqc_results', 'load_trimmomatic_results',
                               'combine_results', 'generate_report']
            for method in required_methods:
                if hasattr(qc_summary, method):
                    logger.info(f"  ✓ QCSummary 有方法: {method}")
                else:
                    logger.error(f"  ✗ QCSummary 缺少方法: {method}")
                    return False

            self.test_results['qc_summary'] = True
            return True

        except Exception as e:
            logger.error(f"  ✗ 质量控制汇总测试失败: {e}")
            self.test_results['qc_summary'] = False
            return False

    def create_test_data(self) -> bool:
        """创建测试数据（如果需要）"""
        if not self.test_data_dir:
            logger.info("未指定测试数据目录，跳过测试数据创建")
            return True

        self.test_data_dir.mkdir(parents=True, exist_ok=True)

        # 创建简单的测试fastq文件
        test_fastq = self.test_data_dir / "test_sample.fastq"
        if not test_fastq.exists():
            logger.info(f"创建测试fastq文件: {test_fastq}")
            with open(test_fastq, 'w') as f:
                # 简单的fastq格式
                f.write("@TEST:1\n")
                f.write("ACGTACGTACGTACGT\n")
                f.write("+\n")
                f.write("FFFFFFFFFFFFFFFF\n")

        # 创建测试样本信息
        test_csv = self.test_data_dir / "test_samples.csv"
        if not test_csv.exists():
            logger.info(f"创建测试样本信息: {test_csv}")
            import pandas as pd
            df = pd.DataFrame({
                'sample': ['test1', 'test2'],
                'group': ['GAO', 'PAL'],
                'replicate': [1, 1],
                'fastq_r1': [str(test_fastq), str(test_fastq)],
                'fastq_r2': [str(test_fastq), str(test_fastq)]
            })
            df.to_csv(test_csv, index=False)

        return True

    def run_all_tests(self) -> bool:
        """运行所有测试"""
        logger.info("=== 开始数据质控模块测试 ===")

        tests = [
            ('模块导入', self.test_module_imports),
            ('脚本执行', self.test_script_execution),
            ('FastQC分析器', self.test_fastqc_analyzer),
            ('Trimmomatic处理器', self.test_trimmomatic_processor),
            ('质量控制汇总', self.test_qc_summary)
        ]

        if self.test_data_dir:
            tests.append(('测试数据创建', self.create_test_data))

        all_passed = True
        for test_name, test_func in tests:
            logger.info(f"\n测试: {test_name}")
            try:
                if test_func():
                    logger.info(f"  ✓ {test_name}: 通过")
                else:
                    logger.error(f"  ✗ {test_name}: 失败")
                    all_passed = False
            except Exception as e:
                logger.error(f"  ✗ {test_name}: 异常 - {e}")
                all_passed = False

        # 输出测试结果汇总
        logger.info("\n=== 测试结果汇总 ===")
        for test_name, result in self.test_results.items():
            status = "✓ 通过" if result else "✗ 失败"
            logger.info(f"  {test_name}: {status}")

        if all_passed:
            logger.info("所有测试通过!")
        else:
            logger.error("部分测试失败!")

        return all_passed


def parse_arguments():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description="数据质控模块测试脚本"
    )

    parser.add_argument(
        "--test-data",
        help="测试数据目录"
    )

    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="详细输出"
    )

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_arguments()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # 运行测试
    tester = QCTester(test_data_dir=args.test_data)
    success = tester.run_all_tests()

    if success:
        logger.info("\n数据质控模块测试完成: 所有测试通过")
        sys.exit(0)
    else:
        logger.error("\n数据质控模块测试完成: 部分测试失败")
        sys.exit(1)


if __name__ == "__main__":
    main()