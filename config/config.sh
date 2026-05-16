#!/usr/bin/env bash

# ============================================================
# 易用转录组分析流程配置文件
# ============================================================
# 一般用户只需要改少数几项：
#   1. 线程和并行数
#   2. 物种相关的 HISAT2 剪接参数
#   3. 文库链特异性
#   4. featureCounts 的注释字段
#
# 所有路径都是相对于项目根目录的相对路径。
# 不建议把这里改成绝对路径，否则别人下载仓库后不容易复现。
#
# 写法说明：
#   THREADS="${THREADS:-8}"
# 表示默认值是 8，也允许临时在命令行覆盖：
#   THREADS=16 bash scripts/run_all.sh

# ============================================================
# 1. 计算资源
# ============================================================
# 默认按家用工作站设置：约 24 逻辑 CPU、16 GB 内存。
# THREADS 是多数单步骤使用的保守线程数，不追求把 CPU 打满。
THREADS="${THREADS:-8}"

# 所有脚本日志统一写入 logs/。
LOG_DIR="${LOG_DIR:-logs}"

# RESUME_MODE=1 表示断点续跑。
# 已经完整生成的过滤后测序文件、比对文件和比对索引会被复用。
# 只有明确要全部重算时，才改为 0。
RESUME_MODE="${RESUME_MODE:-1}"

# ============================================================
# 2. fastp / FastQC 并行设置
# ============================================================
# fastp 单样本多线程收益有限，更适合多个样本同时跑。
# 默认同时处理 4 个样本，每个样本 1 个线程。
# 总 CPU 约为：
#   FASTP_PARALLEL_JOBS * FASTP_THREADS_PER_SAMPLE
FASTP_PARALLEL_JOBS="${FASTP_PARALLEL_JOBS:-4}"
FASTP_THREADS_PER_SAMPLE="${FASTP_THREADS_PER_SAMPLE:-1}"

# FastQC 可以内部使用多线程。默认使用 THREADS。
FASTQC_THREADS="${FASTQC_THREADS:-${THREADS}}"

# ============================================================
# 3. HISAT2 / samtools 并行设置
# ============================================================
# HISAT2 比对更吃 CPU、内存和磁盘 I/O。
# 默认一次比对 1 个样本，每个样本给 8 个 HISAT2 线程。
# 如果机器内存较小、参考基因组很大、数据在移动硬盘上，保持默认更稳。
ALIGN_PARALLEL_JOBS="${ALIGN_PARALLEL_JOBS:-1}"
HISAT2_THREADS_PER_SAMPLE="${HISAT2_THREADS_PER_SAMPLE:-${THREADS}}"

# samtools sort 会额外占用 CPU 和内存。
# 不要把 HISAT2 和 samtools sort 的线程数都开得很高。
SAMTOOLS_SORT_THREADS_PER_SAMPLE="${SAMTOOLS_SORT_THREADS_PER_SAMPLE:-2}"
SAMTOOLS_INDEX_THREADS="${SAMTOOLS_INDEX_THREADS:-2}"

# featureCounts 是一次性读取全部 BAM 计数，使用线程级并行。
FEATURECOUNTS_THREADS="${FEATURECOUNTS_THREADS:-${THREADS}}"

# BAM 质控相对轻量，默认使用 THREADS。
BAM_QC_THREADS="${BAM_QC_THREADS:-${THREADS}}"

# ============================================================
# 4. 输入表格
# ============================================================
# metadata.tsv 和 contrasts.tsv 通常由 scripts/00_make_metadata.sh 自动生成。
# 新手推荐只手动填写 09_metadata/contrasts.csv。
METADATA="${METADATA:-09_metadata/metadata.tsv}"
CONTRASTS="${CONTRASTS:-09_metadata/contrasts.tsv}"

# ============================================================
# 5. 参考基因组和注释
# ============================================================
# 把基因组 FASTA 和注释 GTF/GFF 放入 04_reference/。
# 脚本会优先使用 GENOME_FA 和 ANNOTATION_GTF 指定的文件。
# 如果这两个默认文件不存在，会按后缀自动搜索 REFERENCE_DIR。
REFERENCE_DIR="${REFERENCE_DIR:-04_reference}"
GENOME_FA="${GENOME_FA:-04_reference/genome.fa}"
ANNOTATION_GTF="${ANNOTATION_GTF:-04_reference/annotation.gtf}"

# 自动识别参考文件时使用的后缀。
# NCBI 常见基因组是 .fna，注释常见是 .gff 或 .gff3。
GENOME_SUFFIXES="${GENOME_SUFFIXES:-*.fa *.fasta *.fna *.fas}"
ANNOTATION_SUFFIXES="${ANNOTATION_SUFFIXES:-*.gtf *.gff3 *.gff}"

# HISAT2 index 输出前缀。通常不用改。
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_PREFIX:-05_hisat2_index/genome}"

# FORCE_REBUILD_INDEX=1 会强制重建 HISAT2 index。
# 只有更换参考基因组或怀疑 index 损坏时才需要。
FORCE_REBUILD_INDEX="${FORCE_REBUILD_INDEX:-0}"

# ============================================================
# 6. HISAT2 比对参数
# ============================================================
# HISAT2_SPLICE_MODE 控制是否使用真核剪接比对：
#   auto = 自动使用剪接位点和 exon 提示，适合多数动物、植物、真菌
#   yes  = 强制剪接比对
#   no   = 关闭剪接比对，适合细菌和古菌
HISAT2_SPLICE_MODE="${HISAT2_SPLICE_MODE:-auto}"

# HISAT2_DTA=1 会输出适合转录本组装或下游计数的比对结果。
# 原核 RNA-seq 可设为 0。
HISAT2_DTA="${HISAT2_DTA:-1}"

# 内含子长度范围。
# 动物可保留较大值；植物可适当降低；真菌通常建议更低。
# 示例：
#   动物：HISAT2_MAX_INTRONLEN=500000
#   植物：HISAT2_MAX_INTRONLEN=100000
#   真菌：HISAT2_MAX_INTRONLEN=20000
HISAT2_MIN_INTRONLEN="${HISAT2_MIN_INTRONLEN:-20}"
HISAT2_MAX_INTRONLEN="${HISAT2_MAX_INTRONLEN:-500000}"

# HISAT2_RNA_STRANDNESS 由建库方式决定，不由物种决定。
# 非链特异性：留空
# R1 与转录本同向：FR
# R1 与转录本反向：RF
HISAT2_RNA_STRANDNESS="${HISAT2_RNA_STRANDNESS:-}"

# 额外 HISAT2 参数。普通用户留空。
# 示例：HISAT2_EXTRA_ARGS="--no-mixed --no-discordant"
HISAT2_EXTRA_ARGS="${HISAT2_EXTRA_ARGS:-}"

# ============================================================
# 7. fastp 过滤参数
# ============================================================
# 自动识别双端测序接头。
FASTP_DETECT_ADAPTER_FOR_PE="${FASTP_DETECT_ADAPTER_FOR_PE:-1}"

# 低质量碱基阈值。默认 Q15。
FASTP_QUALIFIED_QUALITY_PHRED="${FASTP_QUALIFIED_QUALITY_PHRED:-15}"

# 单条测序序列中允许低质量碱基的最高比例，单位是百分比。
FASTP_UNQUALIFIED_PERCENT_LIMIT="${FASTP_UNQUALIFIED_PERCENT_LIMIT:-40}"

# 过滤后测序序列的最短长度。
FASTP_LENGTH_REQUIRED="${FASTP_LENGTH_REQUIRED:-30}"

# 额外 fastp 参数。普通用户留空。
# 示例：FASTP_EXTRA_ARGS="--cut_front --cut_tail"
FASTP_EXTRA_ARGS="${FASTP_EXTRA_ARGS:-}"

# ============================================================
# 8. 目录设置
# ============================================================
# 除非你知道自己在改什么，否则不要改这些目录。
ARCHIVE_DIR="${ARCHIVE_DIR:-00_archives}"
RAW_DIR="${RAW_DIR:-00_rawdata}"
CLEAN_DIR="${CLEAN_DIR:-01_clean_fastq}"
FASTQC_DIR="${FASTQC_DIR:-02_fastqc}"
MULTIQC_DIR="${MULTIQC_DIR:-03_multiqc}"
BAM_DIR="${BAM_DIR:-07_bam}"
COUNT_DIR="${COUNT_DIR:-08_counts}"
DESEQ2_DIR="${DESEQ2_DIR:-10_deseq2}"
FIGURE_DIR="${FIGURE_DIR:-11_figures}"
CANDIDATE_DIR="${CANDIDATE_DIR:-13_candidate_genes}"
REPORT_DIR="${REPORT_DIR:-14_report}"

# ============================================================
# 9. 原始数据解包
# ============================================================
# 测序公司交付的 zip、tar.gz 等压缩包会被解包。
# fastq.gz 和 fq.gz 是最终输入，不会继续解压。
UNPACK_NESTED_ARCHIVES="${UNPACK_NESTED_ARCHIVES:-1}"

# 递归解包最多轮数，防止异常压缩包导致无限循环。
UNPACK_MAX_ROUNDS="${UNPACK_MAX_ROUNDS:-5}"

# ============================================================
# 10. featureCounts 参数
# ============================================================
# STRANDNESS:
#   0 = 非链特异性
#   1 = 正向链特异性
#   2 = 反向链特异性
# 必须与建库方式匹配。
STRANDNESS="${STRANDNESS:-0}"

# FEATURE_TYPE 是 GTF/GFF 第三列的 feature 类型。
# 真核基因组通常用 exon。
# 原核 GFF 有时需要改为 CDS 或 gene。
FEATURE_TYPE="${FEATURE_TYPE:-exon}"

# GROUP_ATTRIBUTE 是 GTF/GFF 第 9 列中用于归并测序片段的属性。
# auto 会依次尝试 gene_id、gene、locus_tag、Parent、ID。
# 如果计数结果 Assigned 很低，优先检查这里。
GROUP_ATTRIBUTE="${GROUP_ATTRIBUTE:-auto}"

# 双端转录组测序通常保持下面三个参数为 1。
FEATURECOUNTS_COUNT_READ_PAIRS="${FEATURECOUNTS_COUNT_READ_PAIRS:-1}"
FEATURECOUNTS_REQUIRE_BOTH_ENDS="${FEATURECOUNTS_REQUIRE_BOTH_ENDS:-1}"
FEATURECOUNTS_CHECK_CHIMERA="${FEATURECOUNTS_CHECK_CHIMERA:-1}"

# 额外 featureCounts 参数。普通用户留空。
FEATURECOUNTS_EXTRA_ARGS="${FEATURECOUNTS_EXTRA_ARGS:-}"

# ============================================================
# 11. 质控和差异分析阈值
# ============================================================
# 比对率低于该值时，BAM QC 会提示风险。单位是百分比。
ALIGNMENT_RATE_MIN="${ALIGNMENT_RATE_MIN:-70}"

# DESeq2 低表达过滤：
# 至少 MIN_SAMPLES 个样本中 count >= MIN_COUNT 的基因才保留。
MIN_COUNT="${MIN_COUNT:-10}"
MIN_SAMPLES="${MIN_SAMPLES:-3}"

# 差异基因阈值。
PADJ_CUTOFF="${PADJ_CUTOFF:-0.05}"
LFC_CUTOFF="${LFC_CUTOFF:-1}"

# 每个比较输出前 TOP_N 个候选基因。
TOP_N="${TOP_N:-50}"
