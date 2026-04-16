#!/usr/bin/env python3
"""
比对模块测试脚本

测试功能：
1. Bowtie2索引构建
2. 序列比对
3. SAM到BAM转换
4. 比对统计
5. 质量评估

注意：此测试需要Bowtie2和samtools已安装
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path
import subprocess
import logging

# 添加父目录到路径
sys.path.insert(0, str(Path(__file__).parent.parent))

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class AlignmentModuleTest:
    """比对模块测试类"""

    def __init__(self):
        self.test_dir = Path("test_data/alignment")
        self.test_dir.mkdir(parents=True, exist_ok=True)
        self.results = {}

    def check_tools(self) -> bool:
        """检查必要工具是否安装"""
        tools = ['bowtie2', 'samtools']
        all_available = True

        for tool in tools:
            try:
                result = subprocess.run(
                    [tool, "--version"],
                    capture_output=True,
                    text=True,
                    check=False
                )
                if result.returncode == 0:
                    version_line = result.stdout.strip().split('\n')[0]
                    logger.info(f"{tool}版本: {version_line}")
                else:
                    logger.warning(f"{tool}检查失败: {result.stderr}")
                    all_available = False
            except FileNotFoundError:
                logger.warning(f"未找到工具: {tool}")
                all_available = False

        return all_available

    def test_bowtie2_index_build(self) -> bool:
        """测试Bowtie2索引构建"""
        logger.info("测试Bowtie2索引构建...")

        try:
            # 导入模块
            from scripts.alignment.build_bowtie2_index import Bowtie2IndexBuilder

            # 创建测试FASTA文件（小基因组片段）
            test_fasta = self.test_dir / "test_reference.fa"
            with open(test_fasta, 'w') as f:
                f.write(">chr1\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
                f.write(">chr2\n")
                f.write("GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")

            # 构建索引
            builder = Bowtie2IndexBuilder()
            success = builder.build_index(
                fasta_file=str(test_fasta),
                output_prefix=str(self.test_dir / "test_index"),
                threads=1
            )

            if success:
                logger.info("Bowtie2索引构建测试通过")
                self.results['index_build'] = True
                return True
            else:
                logger.error("Bowtie2索引构建测试失败")
                self.results['index_build'] = False
                return False

        except Exception as e:
            logger.error(f"Bowtie2索引构建测试出错: {e}")
            self.results['index_build'] = False
            return False

    def test_bowtie2_alignment(self) -> bool:
        """测试Bowtie2比对"""
        logger.info("测试Bowtie2比对...")

        try:
            # 导入模块
            from scripts.alignment.run_bowtie2 import Bowtie2Aligner

            # 创建测试fastq文件
            test_fastq = self.test_dir / "test_sample.fastq"
            with open(test_fastq, 'w') as f:
                # 生成一些测试reads
                for i in range(10):
                    f.write(f"@read_{i}\n")
                    f.write("ATCGATCGATCGATCGATCG\n")  # 匹配参考基因组
                    f.write("+\n")
                    f.write("IIIIIIIIIIIIIIIIIIII\n")

            # 检查索引是否存在
            index_prefix = self.test_dir / "test_index"
            if not Path(f"{index_prefix}.1.bt2").exists():
                logger.warning("索引文件不存在，先构建索引")
                if not self.test_bowtie2_index_build():
                    logger.error("无法构建索引，跳过比对测试")
                    self.results['alignment'] = False
                    return False

            # 运行比对
            aligner = Bowtie2Aligner()
            config = {
                'threads': 1,
                'small_rna_mode': True,
                'keep_sam': True  # 保留SAM文件用于测试
            }

            result = aligner.align_single_end(
                fastq_file=str(test_fastq),
                index_prefix=str(index_prefix),
                output_dir=str(self.test_dir / "alignment_output"),
                sample_name="test_sample",
                config=config
            )

            if result.get('success'):
                logger.info("Bowtie2比对测试通过")
                logger.info(f"比对统计: {result.get('stats', {})}")
                self.results['alignment'] = True
                return True
            else:
                logger.error("Bowtie2比对测试失败")
                self.results['alignment'] = False
                return False

        except Exception as e:
            logger.error(f"Bowtie2比对测试出错: {e}")
            self.results['alignment'] = False
            return False

    def test_sam_to_bam_conversion(self) -> bool:
        """测试SAM到BAM转换"""
        logger.info("测试SAM到BAM转换...")

        try:
            # 导入模块
            from scripts.alignment.run_bowtie2 import Bowtie2Aligner

            # 创建测试SAM文件
            test_sam = self.test_dir / "test_conversion.sam"
            with open(test_sam, 'w') as f:
                f.write("@HD\tVN:1.0\tSO:unsorted\n")
                f.write("@SQ\tSN:chr1\tLN:1000\n")
                f.write("read1\t0\tchr1\t1\t255\t10M\t*\t0\t0\tATCGATCGAT\tIIIIIIIIII\n")

            # 测试转换
            aligner = Bowtie2Aligner()
            bam_file = self.test_dir / "test_conversion.bam"

            success = aligner._sam_to_bam(
                sam_file=str(test_sam),
                bam_file=str(bam_file),
                sample_name="test_conversion",
                threads=1,
                keep_sam=True
            )

            if success and bam_file.exists() and bam_file.stat().st_size > 0:
                logger.info("SAM到BAM转换测试通过")
                self.results['sam_to_bam'] = True
                return True
            else:
                logger.error("SAM到BAM转换测试失败")
                self.results['sam_to_bam'] = False
                return False

        except Exception as e:
            logger.error(f"SAM到BAM转换测试出错: {e}")
            self.results['sam_to_bam'] = False
            return False

    def test_alignment_stats(self) -> bool:
        """测试比对统计"""
        logger.info("测试比对统计...")

        try:
            # 导入模块
            from scripts.alignment.alignment_stats import AlignmentStatsCalculator

            # 创建测试Bowtie2日志
            test_log = self.test_dir / "test_stats.log"
            with open(test_log, 'w') as f:
                f.write("10000 reads; of these:\n")
                f.write("  10000 (100.00%) were unpaired; of these:\n")
                f.write("    2000 (20.00%) aligned 0 times\n")
                f.write("    7000 (70.00%) aligned exactly 1 time\n")
                f.write("    1000 (10.00%) aligned >1 times\n")
                f.write("80.00% overall alignment rate\n")

            # 计算统计
            calculator = AlignmentStatsCalculator()
            stats = calculator.calculate_from_log(str(test_log), "test_sample")

            if stats.get('success'):
                logger.info("比对统计测试通过")
                logger.info(f"统计结果: {stats}")
                self.results['alignment_stats'] = True
                return True
            else:
                logger.error("比对统计测试失败")
                self.results['alignment_stats'] = False
                return False

        except Exception as e:
            logger.error(f"比对统计测试出错: {e}")
            self.results['alignment_stats'] = False
            return False

    def test_bam_quality_assessment(self) -> bool:
        """测试BAM质量评估"""
        logger.info("测试BAM质量评估...")

        try:
            # 导入模块
            from scripts.alignment.bam_quality_assessment import BAMQualityAssessor

            # 需要先有一个BAM文件
            bam_file = self.test_dir / "test_conversion.bam"
            if not bam_file.exists():
                logger.warning("BAM文件不存在，先测试SAM到BAM转换")
                if not self.test_sam_to_bam_conversion():
                    logger.error("无法创建BAM文件，跳过质量评估测试")
                    self.results['quality_assessment'] = False
                    return False

            # 运行质量评估
            assessor = BAMQualityAssessor()
            metrics = assessor.assess_single_bam(
                bam_file=str(bam_file),
                output_dir=str(self.test_dir / "quality_output"),
                sample_name="test_sample",
                config={}
            )

            if metrics.get('success'):
                logger.info("BAM质量评估测试通过")
                logger.info(f"质量指标: {metrics.get('total_reads', 'N/A')} reads")
                self.results['quality_assessment'] = True
                return True
            else:
                logger.error("BAM质量评估测试失败")
                self.results['quality_assessment'] = False
                return False

        except Exception as e:
            logger.error(f"BAM质量评估测试出错: {e}")
            self.results['quality_assessment'] = False
            return False

    def run_all_tests(self) -> bool:
        """运行所有测试"""
        logger.info("开始比对模块测试...")

        # 检查工具
        if not self.check_tools():
            logger.warning("部分工具未安装，某些测试可能失败")

        # 运行测试（跳过需要实际比对的测试如果工具不可用）
        tests = [
            ('sam_to_bam_conversion', self.test_sam_to_bam_conversion),
            ('alignment_stats', self.test_alignment_stats),
        ]

        # 如果工具可用，添加更多测试
        if self.check_tools():
            tests.insert(0, ('bowtie2_index_build', self.test_bowtie2_index_build))
            tests.insert(1, ('bowtie2_alignment', self.test_bowtie2_alignment))
            tests.append(('bam_quality_assessment', self.test_bam_quality_assessment))

        # 执行测试
        all_passed = True
        for test_name, test_func in tests:
            logger.info(f"运行测试: {test_name}")
            try:
                if not test_func():
                    all_passed = False
                    logger.error(f"测试失败: {test_name}")
            except Exception as e:
                logger.error(f"测试执行出错 ({test_name}): {e}")
                all_passed = False

        # 输出结果
        logger.info("=" * 50)
        logger.info("比对模块测试结果:")
        for test_name, passed in self.results.items():
            status = "通过" if passed else "失败"
            logger.info(f"  {test_name}: {status}")

        if all_passed:
            logger.info("所有测试通过!")
        else:
            logger.warning("部分测试失败")

        return all_passed

    def cleanup(self):
        """清理测试文件"""
        try:
            if self.test_dir.exists():
                # 保留目录，只删除内容
                for item in self.test_dir.iterdir():
                    if item.is_file():
                        item.unlink()
                    elif item.is_dir():
                        shutil.rmtree(item)
                logger.info("测试文件已清理")
        except Exception as e:
            logger.warning(f"清理测试文件时出错: {e}")


def main():
    """主函数"""
    test = AlignmentModuleTest()

    try:
        # 运行测试
        success = test.run_all_tests()

        # 清理（可选）
        # test.cleanup()

        # 返回退出码
        sys.exit(0 if success else 1)

    except KeyboardInterrupt:
        logger.info("测试被用户中断")
        sys.exit(1)
    except Exception as e:
        logger.error(f"测试过程出错: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()