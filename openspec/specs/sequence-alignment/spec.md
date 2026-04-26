# sequence-alignment Specification

## Purpose
TBD - created by archiving change small-rna-seq-analysis-gao-pal-groups. Update Purpose after archive.
## Requirements
### Requirement: 参考基因组索引构建
系统应构建参考基因组的比对索引，支持Bowtie2比对工具。

#### Scenario: 参考基因组准备
- **WHEN** 输入参考基因组FASTA文件
- **THEN** 系统验证文件格式和完整性，检查染色体命名一致性

#### Scenario: Bowtie2索引构建
- **WHEN** 参考基因组准备就绪
- **THEN** 系统使用bowtie2-build命令构建索引文件，支持small RNA短序列比对优化

### Requirement: 序列比对执行
系统应将质控后的clean reads比对到参考基因组，统计比对效率。

#### Scenario: small RNA优化比对
- **WHEN** 比对small RNA测序数据
- **THEN** 系统使用Bowtie2的--local模式，设置seed length为15，允许最多2个错配，报告最佳比对位置

#### Scenario: 配对端数据处理
- **WHEN** 输入配对端fastq文件
- **THEN** 系统正确处理R1和R2文件，确保配对信息一致性，过滤不配对的reads

#### Scenario: SAM格式输出
- **WHEN** 比对完成
- **THEN** 系统输出SAM格式比对结果，包含比对位置、比对质量、CIGAR字符串等完整信息

### Requirement: 比对结果统计和过滤
系统应统计比对效率，过滤低质量比对，生成唯一比对结果。

#### Scenario: 比对率统计
- **WHEN** 比对完成
- **THEN** 系统统计每个样本的总reads数、比对reads数、唯一比对reads数、多重比对reads数，计算比对率、唯一比对率

#### Scenario: 唯一比对筛选
- **WHEN** 存在多重比对reads
- **THEN** 系统保留比对质量最高的唯一比对，过滤低质量多重比对

#### Scenario: BAM格式转换和排序
- **WHEN** SAM文件生成
- **THEN** 系统将SAM转换为BAM格式，按基因组坐标排序，创建索引文件

### Requirement: 比对质量评估
系统应评估比对质量，识别潜在问题。

#### Scenario: 比对质量分布
- **WHEN** 分析比对结果
- **THEN** 系统生成比对质量分数分布图，识别低质量比对区域

#### Scenario: 基因组覆盖度分析
- **WHEN** 比对完成
- **THEN** 系统分析reads在基因组上的覆盖度，生成覆盖深度分布图，识别覆盖偏好性

