#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] run_all.sh stopped unexpectedly at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check the latest log file in logs/ and verify input files are in the expected folders." >&2' ERR

# Run the full RNA-seq pipeline.
# Input:  config, metadata, contrasts, raw FASTQ, genome FASTA, annotation GTF
# Output: all pipeline results

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

if [[ ! -f "config/config.sh" ]]; then
  echo "[ERROR] Missing config/config.sh" >&2
  echo "[HINT] Run this command from the project folder or restore config/config.sh from the GitHub template." >&2
  exit 1
fi

source config/config.sh

mkdir -p "${LOG_DIR}"

export THREADS LOG_DIR RESUME_MODE
export FASTP_PARALLEL_JOBS FASTP_THREADS_PER_SAMPLE FASTQC_THREADS
export ALIGN_PARALLEL_JOBS HISAT2_THREADS_PER_SAMPLE SAMTOOLS_SORT_THREADS_PER_SAMPLE SAMTOOLS_INDEX_THREADS
export FEATURECOUNTS_THREADS BAM_QC_THREADS
export METADATA CONTRASTS
export REFERENCE_DIR GENOME_FA ANNOTATION_GTF GENOME_SUFFIXES ANNOTATION_SUFFIXES HISAT2_INDEX_PREFIX FORCE_REBUILD_INDEX
export HISAT2_SPLICE_MODE HISAT2_DTA HISAT2_MIN_INTRONLEN HISAT2_MAX_INTRONLEN HISAT2_RNA_STRANDNESS HISAT2_EXTRA_ARGS
export FASTP_DETECT_ADAPTER_FOR_PE FASTP_QUALIFIED_QUALITY_PHRED FASTP_UNQUALIFIED_PERCENT_LIMIT FASTP_LENGTH_REQUIRED FASTP_EXTRA_ARGS
export ARCHIVE_DIR RAW_DIR CLEAN_DIR FASTQC_DIR MULTIQC_DIR BAM_DIR COUNT_DIR DESEQ2_DIR FIGURE_DIR CANDIDATE_DIR REPORT_DIR
export UNPACK_NESTED_ARCHIVES UNPACK_MAX_ROUNDS
export STRANDNESS FEATURE_TYPE GROUP_ATTRIBUTE FEATURECOUNTS_COUNT_READ_PAIRS FEATURECOUNTS_REQUIRE_BOTH_ENDS FEATURECOUNTS_CHECK_CHIMERA FEATURECOUNTS_EXTRA_ARGS
export ALIGNMENT_RATE_MIN MIN_COUNT MIN_SAMPLES PADJ_CUTOFF LFC_CUTOFF TOP_N

run_step() {
  local step_name="$1"
  local step_cmd="$2"
  local log_file="${LOG_DIR}/${step_name}.log"

  echo "[INFO] Starting ${step_name}"
  echo "[INFO] Command: ${step_cmd}"
  echo "[INFO] Log: ${log_file}"

  if ! bash -c "${step_cmd}" > "${log_file}" 2>&1; then
    echo "[ERROR] Step failed: ${step_name}" >&2
    echo "[HINT] Open ${log_file} to see the exact error and suggested fix." >&2
    exit 1
  fi

  echo "[INFO] Finished ${step_name}"
}

run_step "00_unpack_rawdata" "bash scripts/00_unpack_rawdata.sh"
if [[ -f "09_metadata/contrasts.csv" || ! -f "${METADATA}" || ! -f "${CONTRASTS}" ]]; then
  run_step "00_make_metadata" "bash scripts/00_make_metadata.sh"
fi
run_step "00_preflight_check" "bash scripts/00_preflight_check.sh"
run_step "00_prepare_reference" "bash scripts/00_prepare_reference.sh"
run_step "01_fastp_qc" "bash scripts/01_fastp_qc.sh"
run_step "02_hisat2_align" "bash scripts/02_hisat2_align.sh"
run_step "03_bam_qc" "bash scripts/03_bam_qc.sh"
run_step "04_featureCounts" "bash scripts/04_featureCounts.sh"
run_step "05_DESeq2" "Rscript scripts/05_DESeq2.R"
run_step "06_generate_report" "Rscript scripts/06_generate_report.R"

echo "[INFO] Full RNA-seq pipeline finished"
