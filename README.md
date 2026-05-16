# 一键式 RNA-seq 分析流程使用教程

这个仓库是一套双端 RNA-seq 标准分析流程。用户把原始数据、参考基因组和比较关系放到指定位置后，可以一条命令完成质控、比对、计数、差异分析、基础作图和 HTML 分析报告。

需要准备三类文件：

1. 原始测序数据，通常是 `.fastq.gz` 或测序公司给的 `.zip`、`.tar.gz`
2. 参考基因组 FASTA，例如 `.fa`、`.fasta`、`.fna`
3. 基因注释文件，例如 `.gtf`、`.gff3`、`.gff`

正常情况下，只需要手动填写 `09_metadata/contrasts.csv`。样本表会由程序根据 FASTQ 文件名自动生成。

本模板默认所有样本来自同一批次。存在不同提取批次、建库批次或上机批次时，应重新设计模型，不要把批次问题塞进这个入门模板。

## 1. 最快使用方法

下载仓库后，进入项目目录：

```bash
cd rnaseq-pipeline-template
```

安装软件环境：

```bash
mamba env create -f environment.yml
conda activate rnaseq_pipeline
```

放入原始数据：

```text
00_rawdata/      # 放 .fastq.gz 或 .fq.gz
00_archives/    # 放测序公司交付的 .zip、.tar.gz 等压缩包
```

放入参考基因组和注释：

```text
04_reference/   # 放一个基因组 FASTA 和一个 GTF/GFF 注释文件
```

填写比较关系：

```bash
cp 09_metadata/contrasts.example.csv 09_metadata/contrasts.csv
```

编辑 `09_metadata/contrasts.csv`：

```text
comparison,numerator,denominator
Treatment_vs_Control,Treatment,Control
```

检查配置文件：

```text
config/config.sh
```

多数用户只需要看线程数、物种类型、链特异性和 featureCounts 注释字段。配置文件内已经写了中文注释。

运行全流程：

```bash
bash scripts/run_all.sh
```

中途停止后，修好报错，再运行同一条命令。默认 `RESUME_MODE=1`，已经完成的样本级结果会自动跳过。

第一次使用时，主要关心这些文件：

```text
00_rawdata/ 或 00_archives/        # 原始数据
04_reference/                      # 参考基因组和注释
09_metadata/contrasts.csv          # 唯一推荐手动填写的表
config/config.sh                   # 少量参数
```

## 2. 这个流程会做什么

流程顺序如下：

```text
原始 FASTQ
  -> 预检查样本表、比较表、参考文件和软件环境
  -> fastp 去接头和低质量 reads
  -> FastQC / MultiQC 质控汇总
  -> HISAT2 比对到参考基因组
  -> samtools 整理 BAM 和统计比对质量
  -> featureCounts 生成基因 raw counts
  -> DESeq2 做差异表达分析
  -> PCA / 样本相关性 / DEG 表 / 火山图 / MA 图 / 候选基因表
```

最终最重要的结果是：

```text
03_multiqc/multiqc_fastq_qc.html
03_multiqc/multiqc_bam_qc.html
logs/pipeline_handoff_plan.tsv
08_counts/gene_counts.txt
10_deseq2/DEG_summary.csv
10_deseq2/comparisons/*.significant_DEGs.csv
11_figures/PCA_condition.pdf
11_figures/sample_correlation_heatmap.pdf
11_figures/*.volcano.pdf
11_figures/*.MA_plot.pdf
11_figures/DEG_UpSet_plot.pdf
11_figures/DEG_UpSet_up_genes.pdf
11_figures/DEG_UpSet_down_genes.pdf
12_enrichment/README.md
13_candidate_genes/*.top_candidate_genes_ranked.csv
14_report/RNAseq_analysis_report.html
```

## 3. 下载后先看目录

不要改目录名。脚本按这些相对路径找文件。

```text
rnaseq-pipeline-template/
├── 00_archives/                 # 放测序公司交付压缩包
├── 00_rawdata/                  # 放最终 FASTQ 文件
├── 01_clean_fastq/              # fastp 输出
├── 02_fastqc/                   # fastp 和 FastQC 报告
├── 03_multiqc/                  # MultiQC 汇总报告
├── 04_reference/                # 放 genome.fa 和 annotation.gtf
├── 05_hisat2_index/             # HISAT2 index
├── 07_bam/                      # sorted BAM
├── 08_counts/                   # featureCounts 计数结果
├── 09_metadata/                 # 样本表和比较表
├── 10_deseq2/                   # DESeq2 结果
├── 11_figures/                  # PCA、热图、火山图、MA 图
├── 12_enrichment/               # 富集分析结果预留目录
├── 13_candidate_genes/          # 候选基因表
├── 14_report/                   # HTML 分析报告
├── config/config.sh             # 参数配置
├── logs/                        # 每一步日志
├── scripts/                     # 所有分析脚本
├── environment.yml              # conda 环境文件
└── README.md
```

这些数据目录和结果目录里都放了一个可见的 `README.md` 占位文件。GitHub 不保存真正的空文件夹；保留这些占位文件，别人下载 zip 后才能看到完整目录结构。

## 4. 安装软件环境

推荐使用 mamba。

```bash
mamba env create -f environment.yml
conda activate rnaseq_pipeline
```

如果以前已经创建过环境，后来仓库更新了依赖，用下面的命令补齐软件包：

```bash
mamba env update -n rnaseq_pipeline -f environment.yml
conda activate rnaseq_pipeline
```

检查软件是否可用：

```bash
fastp --version
fastqc --version
multiqc --version
hisat2 --version
samtools --version
featureCounts -v
Rscript --version
```

这些命令都能输出版本号，说明环境基本可用。某个命令不存在，就先修环境，不要继续跑流程。

## 5. 放原始数据

原始数据有两种常见情况。

### 5.1 已经拿到 FASTQ

如果文件是 `.fastq.gz` 或 `.fq.gz`，直接放到：

```text
00_rawdata/
```

示例：

```text
00_rawdata/
├── Control_1_R1.fastq.gz
├── Control_1_R2.fastq.gz
├── Control_2_R1.fastq.gz
├── Control_2_R2.fastq.gz
├── Treatment_1_R1.fastq.gz
└── Treatment_1_R2.fastq.gz
```

`.fastq.gz` 是标准 FASTQ 输入格式。不要手动解成 `.fastq`。

### 5.2 拿到的是测序公司压缩包

如果测序公司给的是 `.zip`、`.tar`、`.tar.gz`、`.tgz`、`.tar.bz2`、`.tbz2`、`.tar.xz` 或 `.txz`，放到：

```text
00_archives/
```

然后运行：

```bash
bash scripts/00_unpack_rawdata.sh
```

本步做的事：

```text
交付压缩包
  -> 解包到 00_rawdata/
  -> 递归解开内层 zip 或 tar 包
  -> 保留 fastq.gz 和 fq.gz，不再解压
  -> 生成 FASTQ 路径清单
```

本步输出：

```text
logs/raw_fastq_files_after_unpack.tsv
logs/unpacked_delivery_archives.tsv
```

看 `logs/raw_fastq_files_after_unpack.tsv`。后面填写 metadata 时，使用这个文件里列出的 FASTQ 路径。

## 6. 放参考基因组和注释

把参考基因组和注释文件放到：

```text
04_reference/
```

不需要改成固定文件名。脚本按后缀自动识别。

```text
基因组：*.fa、*.fasta、*.fna、*.fas
注释：*.gtf、*.gff3、*.gff
```

NCBI 下载的文件通常可以直接放进去，例如：

```text
04_reference/GCF_000000000.1_genomic.fna
04_reference/genomic.gff
```

要求：

- 同一目录中只能放一个基因组 FASTA 和一个注释文件
- 基因组和注释必须来自同一参考版本
- 如果有多个参考文件，脚本会停止并要求你删除多余文件或在 `config.sh` 中指定精确路径
- 参考和注释不匹配，会导致比对率低或 counts 分配率低

## 7. FASTQ 命名规则

默认不需要填写样本表。脚本会从 FASTQ 文件名自动推断 `sample_id` 和 `condition`。

推荐命名：

```text
00_rawdata/Control_1_R1.fastq.gz
00_rawdata/Control_1_R2.fastq.gz
00_rawdata/Control_2_R1.fastq.gz
00_rawdata/Control_2_R2.fastq.gz
00_rawdata/Treatment_1_R1.fastq.gz
00_rawdata/Treatment_1_R2.fastq.gz
```

脚本会推断为：

```text
sample_id     condition
Control_1     Control
Control_2     Control
Treatment_1   Treatment
```

也支持这类测序公司常见命名：

```text
CGA_1.dedup.R1.fastq.gz
CGA_1.dedup.R2.fastq.gz
WT_1.dedup.R1.fastq.gz
WT_1.dedup.R2.fastq.gz
```

脚本会推断为：

```text
sample_id   condition
CGA_1       CGA
WT_1        WT
```

如果文件名不符合这个规律，再使用兜底表：

```bash
cp 09_metadata/samples.example.csv 09_metadata/samples.csv
```

`samples.csv` 只需要两列：

```text
sample_id,condition
CGA_1,CGA
CGA_2,CGA
WT_1,WT
WT_2,WT
```

多数情况下不需要 `samples.csv`。

## 8. 填写比较关系

这是唯一推荐手动填写的表。CSV 可以用 Excel、WPS 或文本编辑器打开，不需要手动输入 tab。

```bash
cp 09_metadata/contrasts.example.csv 09_metadata/contrasts.csv
```

打开 `09_metadata/contrasts.csv`，写清楚要比较哪些组。

格式：

```text
comparison,numerator,denominator
Treatment_vs_Control,Treatment,Control
```

复制示例文件后，必须把示例行改成自己的比较。不要在示例行下面直接追加新比较后保留原示例行。

`contrasts.csv` 是唯一推荐手动修改的比较文件。`metadata.tsv` 和 `contrasts.tsv` 是程序生成的派生文件，不要手动改。

每次运行 `run_all.sh` 时，只要检测到 `09_metadata/contrasts.csv`，程序都会先删除旧的：

```text
09_metadata/metadata.tsv
09_metadata/contrasts.tsv
logs/metadata_fastq_match.tsv
```

然后从 FASTQ 文件名和 `contrasts.csv` 重新生成。这样旧 TSV 中的残留行不会进入后续分析。

如果不小心在 `contrasts.csv` 保留了这一行：

```text
Treatment_vs_Control,Treatment,Control
```

而你的 FASTQ 中没有 `Treatment` 和 `Control` 分组，脚本会自动忽略这条示例行，并在日志中给出警告。

含义：

| 列名 | 怎么填 |
|---|---|
| `comparison` | 输出文件名前缀 |
| `numerator` | 分子组，log2FC 的上方 |
| `denominator` | 分母组，log2FC 的下方 |

上面的例子表示：

```text
Treatment_vs_Control = Treatment / Control
```

如果某个基因 log2FoldChange 大于 0，表示它在 `Treatment` 中更高。

运行 `run_all.sh` 时，脚本会自动生成程序真正使用的：

```text
09_metadata/metadata.tsv
09_metadata/contrasts.tsv
logs/metadata_fastq_match.tsv
```

先看 `logs/metadata_fastq_match.tsv`。这个文件会列出每个样本推断出的分组和匹配到的 R1/R2。

## 9. 检查配置文件

配置文件是：

```text
config/config.sh
```

多数初学者只需要检查四类参数。

### 9.1 线程与并行

```bash
THREADS=8
FASTP_PARALLEL_JOBS=4
FASTP_THREADS_PER_SAMPLE=1
FASTQC_THREADS=8

ALIGN_PARALLEL_JOBS=1
HISAT2_THREADS_PER_SAMPLE=8
SAMTOOLS_SORT_THREADS_PER_SAMPLE=2
SAMTOOLS_INDEX_THREADS=2

FEATURECOUNTS_THREADS=8
BAM_QC_THREADS=8
RESUME_MODE=1
```

这套默认值按普通家用电脑设置，约 24 逻辑 CPU、16GB 内存。它不追求把 CPU 打满，而是优先避免内存、磁盘 I/O 和排序步骤互相拖慢。

`THREADS` 是重计算步骤的保守线程预算。fastp 主要按样本并行提速，默认是 `4` 个样本并行、每个样本 `1` 线程。

如果使用 SSD，并且想提高 fastp 速度，可以改大：

```bash
FASTP_PARALLEL_JOBS=8
FASTP_THREADS_PER_SAMPLE=1
```

如果数据在机械硬盘、移动硬盘或网络盘上，fastp 同时跑太多样本可能被磁盘 I/O 拖慢，保持默认或继续降低：

```bash
FASTP_PARALLEL_JOBS=2
FASTP_THREADS_PER_SAMPLE=1
```

`RESUME_MODE=1` 表示默认断点续跑。已经完整生成的 clean FASTQ、fastp 报告、BAM 和索引会被复用。只有明确想全部重算时，才改成：

```bash
RESUME_MODE=0
```

HISAT2 默认一次只比对一个样本，因为比对阶段更吃内存和磁盘。如果机器配置较高，可以改成：

```bash
ALIGN_PARALLEL_JOBS=2
HISAT2_THREADS_PER_SAMPLE=4
SAMTOOLS_SORT_THREADS_PER_SAMPLE=2
```

不要让并行任务总线程数长期明显超过机器实际 CPU 核心数。内存不足、机械硬盘或移动硬盘运行时，不要盲目提高并行数。

### 9.2 参考文件识别

默认自动识别 `04_reference/` 中的参考文件：

```bash
REFERENCE_DIR=04_reference
GENOME_SUFFIXES="*.fa *.fasta *.fna *.fas"
ANNOTATION_SUFFIXES="*.gtf *.gff3 *.gff"
```

如果目录里只有一个 `.fna` 和一个 `.gff`，不需要改名。只有在同一目录有多个候选文件时，才需要在 `config/config.sh` 中手动指定：

```bash
GENOME_FA=04_reference/your_genome.fna
ANNOTATION_GTF=04_reference/your_annotation.gff
```

### 9.3 物种类型

HISAT2 的主流程可以用于动物、植物和真菌，但剪接参数不应机械相同。

推荐起点：

| 数据类型 | 推荐设置 |
|---|---|
| 动物真核转录组 | `HISAT2_SPLICE_MODE=auto`，`HISAT2_MAX_INTRONLEN=500000`，`HISAT2_DTA=1` |
| 高等植物转录组 | `HISAT2_SPLICE_MODE=auto`，`HISAT2_MAX_INTRONLEN=100000` 到 `500000`，`HISAT2_DTA=1` |
| 真菌转录组 | `HISAT2_SPLICE_MODE=auto`，`HISAT2_MAX_INTRONLEN=20000` 到 `50000`，`HISAT2_DTA=1` |
| 细菌或古菌转录组 | `HISAT2_SPLICE_MODE=no`，`HISAT2_DTA=0` |

原核样品没有真核剪接结构，应关闭剪接比对。

### 9.4 链特异性

链特异性由建库方式决定，不由物种决定。拿不准时先看测序公司报告。

| 文库类型 | HISAT2 设置 | featureCounts 设置 |
|---|---|---|
| 非链特异性 paired-end | `HISAT2_RNA_STRANDNESS=` | `STRANDNESS=0` |
| R1 与转录本同向 | `HISAT2_RNA_STRANDNESS=FR` | `STRANDNESS=1` |
| R1 与转录本反向 | `HISAT2_RNA_STRANDNESS=RF` | `STRANDNESS=2` |

链特异性设置错，会明显降低 `featureCounts` 的 Assigned 比例。

## 10. 一键运行

确认下面这些文件已经存在：

```text
04_reference/ 中的基因组 FASTA
04_reference/ 中的 GTF/GFF 注释
09_metadata/contrasts.csv
```

如果存在 `09_metadata/contrasts.csv`，`run_all.sh` 会先调用 `00_make_metadata.sh`，从 FASTQ 文件名和 `contrasts.csv` 重建 TSV。旧的 `metadata.tsv` 和 `contrasts.tsv` 会被删除。若要完全手动使用 `metadata.tsv` 和 `contrasts.tsv`，不要保留 `contrasts.csv`。

原始数据至少满足以下一种情况：

```text
00_rawdata/*.fastq.gz
00_archives/*.zip
00_archives/*.tar.gz
```

运行完整流程：

```bash
bash scripts/run_all.sh
```

如果中途停止，修正问题后再次运行同一条命令。默认 `RESUME_MODE=1`，已经完成的样本级结果会跳过，不会把 fastp 和 HISAT2 从头全部重跑。

临时指定线程运行：

```bash
THREADS=16 bash scripts/run_all.sh
```

这只会覆盖 `THREADS`。若要改变 fastp 或 HISAT2 的样本并行数，直接编辑 `config/config.sh` 更清楚。

本步会依次运行：

```text
00_unpack_rawdata.sh
00_make_metadata.sh
00_preflight_check.sh
00_prepare_reference.sh
01_fastp_qc.sh
02_hisat2_align.sh
03_bam_qc.sh
04_featureCounts.sh
05_DESeq2.R
06_generate_report.R
```

失败时流程会停止。先看终端里的 `[ERROR]` 和 `[HINT]`，再打开 `logs/` 中对应步骤的日志。

## 11. 分步运行

第一次使用建议分步跑。每一步成功后再进入下一步。

### 11.1 解包原始数据

用 `00_unpack_rawdata.sh` 把测序交付包整理成可用于 metadata 的 FASTQ 路径。

```bash
bash scripts/00_unpack_rawdata.sh
```

输入：

```text
00_archives/ 中的交付压缩包
00_rawdata/ 中已有的 FASTQ
```

输出：

```text
logs/raw_fastq_files_after_unpack.tsv
logs/unpacked_delivery_archives.tsv
```

解读：

- `raw_fastq_files_after_unpack.tsv` 列出最终 FASTQ 路径
- `.fastq.gz` 不会继续解压
- 如果没有 FASTQ，脚本会停止

### 11.2 预检查

用 `00_preflight_check.sh` 在正式分析前检查样本表、比较表、参考文件、FASTQ 路径和软件命令，提前发现会导致下游卡住的问题。

```bash
bash scripts/00_preflight_check.sh
```

输入：

```text
config/config.sh
09_metadata/metadata.tsv
09_metadata/contrasts.tsv
04_reference/genome.fa
04_reference/annotation.gtf
metadata 中 fastq_1 和 fastq_2 指向的 FASTQ
```

输出：

```text
logs/preflight_check.tsv
logs/pipeline_handoff_plan.tsv
```

解读：

- `preflight_check.tsv` 全部为 `PASS`，说明基本输入可进入正式分析
- `pipeline_handoff_plan.tsv` 列出每个样本从 raw FASTQ 到 clean FASTQ 再到 BAM 的预期文件名
- 这一步失败时，先修 metadata、contrasts、参考文件或软件环境，不要继续跑

### 11.3 准备参考

用 HISAT2 建立基因组 index，并从 GTF 中提取剪接位点信息。

```bash
bash scripts/00_prepare_reference.sh
```

输入：

```text
04_reference/genome.fa
04_reference/annotation.gtf
```

输出：

```text
05_hisat2_index/genome*.ht2
05_hisat2_index/splice_sites.txt
05_hisat2_index/exons.txt
logs/hisat2_build.log
logs/reference_genome.seqkit_stats.txt
logs/annotation_gtf_record_count.tsv
```

解读：

- `genome*.ht2` 存在，说明 index 建好
- `annotation_gtf_record_count.tsv` 可粗略检查 GTF 中 gene 和 exon 数量
- 原核模式下 `splice_sites.txt` 和 `exons.txt` 为空是正常的

### 11.4 FASTQ 质控和过滤

用 fastp 去接头、过滤低质量 reads，并用 FastQC 和 MultiQC 汇总质控。

```bash
bash scripts/01_fastp_qc.sh
```

输入：

```text
09_metadata/metadata.tsv
metadata 中 fastq_1 和 fastq_2 指向的 FASTQ
```

输出：

```text
01_clean_fastq/*_clean_R1.fastq.gz
01_clean_fastq/*_clean_R2.fastq.gz
02_fastqc/*.fastp.html
02_fastqc/*_fastqc.html
03_multiqc/multiqc_fastq_qc.html
```

解读：

- 先看 `multiqc_fastq_qc.html`
- 重点看 Q30、接头污染、GC 分布和过滤比例
- 某个样本明显异常，应先查原始数据和样本表

### 11.5 reads 比对

用 HISAT2 将 clean reads 比对到参考基因组，并生成排序后的 BAM。

```bash
bash scripts/02_hisat2_align.sh
```

输入：

```text
01_clean_fastq/*_clean_R1.fastq.gz
01_clean_fastq/*_clean_R2.fastq.gz
05_hisat2_index/genome*.ht2
```

输出：

```text
07_bam/*.sorted.bam
07_bam/*.sorted.bam.bai
logs/*.hisat2.log
logs/hisat2_alignment_rate_summary.txt
```

解读：

- 先看 `hisat2_alignment_rate_summary.txt`
- 同一项目内样本比对率应大体接近
- 比对率低时，优先检查参考版本、污染、宿主 RNA 混入和 FASTQ 质量

### 11.6 BAM 质控

用 samtools 检查 BAM 文件的比对质量，并用 MultiQC 汇总。

```bash
bash scripts/03_bam_qc.sh
```

输入：

```text
07_bam/*.sorted.bam
```

输出：

```text
07_bam/qc/*.flagstat.txt
07_bam/qc/*.stats.txt
07_bam/qc/mapping_rate_check.tsv
03_multiqc/multiqc_bam_qc.html
```

解读：

- `mapping_rate_check.tsv` 会标记低比对率样本
- `flagstat.txt` 中 mapped 和 properly paired 过低时，不应继续解释 DEG

### 11.7 基因计数

用 featureCounts 把 BAM 中的 reads 汇总到基因层面，生成 DESeq2 需要的 raw counts。

```bash
bash scripts/04_featureCounts.sh
```

输入：

```text
07_bam/*.sorted.bam
04_reference/annotation.gtf
09_metadata/metadata.tsv
```

输出：

```text
08_counts/gene_counts.txt
08_counts/gene_counts.txt.summary
logs/featureCounts.log
```

解读：

- `gene_counts.txt` 是 DESeq2 的输入
- `gene_counts.txt.summary` 中 `Assigned` 越高越好
- `Unassigned_NoFeatures` 高，常见原因是参考基因组和 GTF 不匹配
- `Assigned` 很低，也要检查 `STRANDNESS`

### 11.8 差异表达分析

用 DESeq2 基于 raw counts 做差异表达分析，并输出图表和 DEG 表。

```bash
Rscript scripts/05_DESeq2.R
```

输入：

```text
08_counts/gene_counts.txt
09_metadata/metadata.tsv
09_metadata/contrasts.tsv
```

输出：

```text
10_deseq2/normalized_counts.csv
10_deseq2/vst_expression_matrix.csv
10_deseq2/DEG_summary.csv
10_deseq2/comparisons/*.all_results.csv
10_deseq2/comparisons/*.significant_DEGs.csv
11_figures/PCA_condition.pdf
11_figures/sample_correlation_heatmap.pdf
11_figures/*.volcano.pdf
11_figures/*.MA_plot.pdf
11_figures/DEG_UpSet_plot.pdf
11_figures/DEG_UpSet_up_genes.pdf
11_figures/DEG_UpSet_down_genes.pdf
12_enrichment/README.md
13_candidate_genes/*.top_candidate_genes_ranked.csv
14_report/RNAseq_analysis_report.html
```

解读：

- 先看 PCA 和样本相关性
- 样本关系正常后，再看 DEG 数量和 DEG 表
- DEG 是统计关联，不是因果证明

### 11.9 生成 HTML 分析报告

用 `06_generate_report.R` 汇总材料方法、质控结果、差异表达结果、候选基因和结果解读，生成一个可直接打开的网页报告。

```bash
Rscript scripts/06_generate_report.R
```

输入：

```text
09_metadata/metadata.tsv
09_metadata/contrasts.tsv
logs/*.tsv
03_multiqc/*.html
08_counts/gene_counts.txt.summary
10_deseq2/DEG_summary.csv
10_deseq2/comparisons/*.csv
11_figures/*.pdf
13_candidate_genes/*.csv
```

输出：

```text
14_report/RNAseq_analysis_report.html
```

解读：

- 这是最终给用户查看的总报告
- 报告内包含材料方法、项目概览、质控结果、差异表达结果和候选基因表
- 图表以 PDF 形式嵌入或链接，若浏览器不直接显示 PDF，可点击链接单独打开

## 12. 结果先看什么

按这个顺序检查：

1. `03_multiqc/multiqc_fastq_qc.html`
2. `logs/hisat2_alignment_rate_summary.txt`
3. `03_multiqc/multiqc_bam_qc.html`
4. `08_counts/gene_counts.txt.summary`
5. `11_figures/PCA_condition.pdf`
6. `11_figures/sample_correlation_heatmap.pdf`
7. `10_deseq2/DEG_summary.csv`
8. `10_deseq2/comparisons/*.significant_DEGs.csv`
9. `13_candidate_genes/*.top_candidate_genes_ranked.csv`
10. `14_report/RNAseq_analysis_report.html`

判断原则：

- 质控不过，不解释差异基因
- 样本聚类异常，不解释差异基因
- 只看 DEG 数量没有意义，必须结合样本关系和研究设计
- 候选基因表用于后续注释、富集和实验验证

## 13. 常见报错怎么处理

脚本停止时会输出两类信息：

```text
[ERROR] 停止原因
[HINT] 处理建议
```

常见情况：

| 报错 | 常见原因 | 处理 |
|---|---|---|
| `Contrast file not found` | 没有创建 `contrasts.csv` | 运行 `cp 09_metadata/contrasts.example.csv 09_metadata/contrasts.csv` 后编辑 |
| `Removed leftover example contrast` | 比较表中存在默认示例行，但样本里没有对应分组 | 脚本已自动忽略这行；最终使用的文件以重新生成的 `contrasts.tsv` 为准 |
| `Cannot detect R1/R2` | FASTQ 文件名没有明确 R1/R2 | 改成 `Sample_1_R1.fastq.gz` 和 `Sample_1_R2.fastq.gz` |
| `Ambiguous FASTQ pairing` | 一个样本匹配到多个 R1 或 R2 | 改 FASTQ 文件名，或提供 `samples.csv` 兜底 |
| `Conditions in contrasts are not present` | 自动推断出的分组和 contrasts 不一致 | 检查 FASTQ 命名，或提供 `samples.csv` 明确定义分组 |
| `Metadata file not found` | 没有创建 `metadata.tsv` | 优先运行 `bash scripts/00_make_metadata.sh` 自动生成 |
| `Metadata header is invalid` | 表头不是规定格式 | 第一行必须是 `sample_id condition fastq_1 fastq_2`，用 tab 分隔 |
| `R1 FASTQ not found` | FASTQ 路径写错 | 按 `logs/raw_fastq_files_after_unpack.tsv` 改 metadata |
| `No FASTQ files found in 00_rawdata` | 没有放 FASTQ 或交付包 | 把 `.fastq.gz` 放入 `00_rawdata/`，或把交付压缩包放入 `00_archives/` |
| `Nested archive unpacking stopped` | 内层压缩包嵌套过深 | 调高 `UNPACK_MAX_ROUNDS`，或手动检查 `00_rawdata/` |
| `HISAT2 index not found` | 没有建 index | 先运行 `bash scripts/00_prepare_reference.sh` |
| `Clean FASTQ not found` | 没有运行 fastp | 先运行 `bash scripts/01_fastp_qc.sh` |
| `BAM not found` | 没有完成比对 | 先运行 `bash scripts/02_hisat2_align.sh` |
| `Invalid STRANDNESS` | 链特异性参数写错 | 使用 `0`、`1` 或 `2` |
| `metadata.tsv contains a batch column` | 加入了 batch 列 | 删除 batch 列，本模板只处理同一批次样本 |
| `Conditions in contrasts.tsv but not in metadata` | 分组名不一致 | 保证 contrasts 中的分组名出现在 metadata 的 `condition` 列 |

不要跳过报错继续运行。报错说明当前输入不满足后续分析要求。

## 14. 输出文件含义

| 文件 | 含义 | 怎么看 |
|---|---|---|
| `multiqc_fastq_qc.html` | FASTQ 质控汇总 | 看测序质量、接头污染、GC 分布 |
| `hisat2_alignment_rate_summary.txt` | HISAT2 比对率 | 样本间应大体一致 |
| `multiqc_bam_qc.html` | BAM 质控汇总 | 看 mapped、properly paired |
| `gene_counts.txt.summary` | reads 分配统计 | `Assigned` 越高越好 |
| `normalized_counts.csv` | DESeq2 标准化 counts | 用于展示，不用于重新做差异分析 |
| `vst_expression_matrix.csv` | 方差稳定化矩阵 | 用于 PCA、聚类、热图 |
| `PCA_condition.pdf` | PCA 图 | 看重复是否聚在一起 |
| `sample_correlation_heatmap.pdf` | 样本相关性热图 | 看组内一致性 |
| `DEG_summary.csv` | DEG 数量汇总 | 比较不同 contrast 的扰动强度 |
| `*.all_results.csv` | 所有基因差异结果 | 保留完整统计结果 |
| `*.significant_DEGs.csv` | 显著差异基因 | 默认 `padj < 0.05` 且 `abs(log2FC) >= 1` |
| `*.volcano.pdf` | 火山图 | 看显著性和变化倍数 |
| `*.MA_plot.pdf` | MA 图 | 看差异是否集中在低表达基因 |
| `DEG_UpSet_plot.pdf` | 全部显著 DEG 交集图 | 不区分上下调 |
| `DEG_UpSet_up_genes.pdf` | 上调基因交集图 | 只统计上调 DEG |
| `DEG_UpSet_down_genes.pdf` | 下调基因交集图 | 只统计下调 DEG |
| `*.top_candidate_genes_ranked.csv` | 候选基因排序 | 用于后续注释和验证 |
| `RNAseq_analysis_report.html` | 最终 HTML 报告 | 汇总材料方法、结果和解读 |
