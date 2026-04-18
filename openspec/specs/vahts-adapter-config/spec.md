## ADDED Requirements

### Requirement: VAHTS接头序列配置

系统 SHALL 支持配置和使用 VAHTS Small RNA Library Prep Kit for Illumina V2 文库结构的接头序列。

#### Scenario: 加载VAHTS接头序列配置

- **WHEN** 用户在配置文件中设置 `adapter_type: vahts_small_rna_v2`
- **THEN** 系统 SHALL 加载对应的VAHTS V2接头序列

#### Scenario: 获取VAHTS接头序列

- **WHEN** 系统需要VAHTS接头序列时
- **THEN** 系统 SHALL 返回正确的VAHTS V2接头序列：`AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`

#### Scenario: 使用VAHTS接头文件

- **WHEN** 系统使用VAHTS接头配置时
- **THEN** 系统 SHALL 使用 `config/VAHTS-SmallRNA-V2.fa` 文件中的接头序列

### Requirement: 接头序列类型选择

系统 SHALL 支持通过配置文件选择不同类型的接头序列。

#### Scenario: 通过配置文件设置接头类型

- **WHEN** 用户在 config.yaml 中设置 `adapter_type` 参数
- **THEN** 系统 SHALL 使用该配置指定的接头类型

#### Scenario: 默认使用VAHTS接头

- **WHEN** 用户未配置 `adapter_type` 参数
- **THEN** 系统 SHALL 默认使用VAHTS V2接头类型

#### Scenario: 支持多种接头类型

- **WHEN** 用户配置不同的接头类型（如 `illumina_small_rna`, `truseq`, `vahts_small_rna_v2` 等）
- **THEN** 系统 SHALL 正确加载对应的接头序列

## MODIFIED Requirements

### Requirement: Trimmomatic接头去除配置

系统 SHALL 配置Trimmomatic使用正确的接头序列进行质量控制。

#### Scenario: 使用VAHTS接头进行Trimmomatic处理

- **WHEN** 系统使用Trimmomatic进行质量修剪和接头去除
- **THEN** 系统 SHALL 使用VAHTS V2接头序列进行接头去除

#### Scenario: 支持FASTA格式接头文件

- **WHEN** 系统需要加载接头序列
- **THEN** 系统 SHALL 支持从FASTA格式的接头序列文件中读取数据

#### Scenario: 配置Trimmomatic参数

- **WHEN** 系统配置Trimmomatic的ILLUMINACLIP参数
- **THEN** 系统 SHALL 使用正确的接头文件路径和参数设置
