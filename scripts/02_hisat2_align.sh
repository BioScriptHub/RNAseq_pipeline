#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check HISAT2/samtools installation, index files, clean FASTQ files, disk space, and the log above." >&2' ERR

# Align clean paired-end reads to reference genome using HISAT2.
# Input:  clean FASTQ files and HISAT2 index
# Output: sorted BAM files, BAM index files, HISAT2 logs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

source config/config.sh

die() {
  echo "[ERROR] $1" >&2
  echo "[HINT] $2" >&2
  exit 1
}

is_nonnegative_int() {
  [[ "$1" =~ ^[0-9]+$ ]]
}

is_positive_int() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

mkdir -p "${BAM_DIR}" "${LOG_DIR}"

if [[ ! -f "${METADATA}" ]]; then
  die "Metadata file not found: ${METADATA}" \
      "Create it from 09_metadata/metadata.example.tsv before running alignment."
fi

if [[ ! -f "${HISAT2_INDEX_PREFIX}.1.ht2" && ! -f "${HISAT2_INDEX_PREFIX}.1.ht2l" ]]; then
  die "HISAT2 index not found with prefix: ${HISAT2_INDEX_PREFIX}" \
      "Run bash scripts/00_prepare_reference.sh first, or check HISAT2_INDEX_PREFIX in config/config.sh."
fi

case "${HISAT2_SPLICE_MODE}" in
  auto|yes|no)
    ;;
  *)
    die "Invalid HISAT2_SPLICE_MODE: ${HISAT2_SPLICE_MODE}" \
        "Use auto, yes, or no in config/config.sh. Use no for bacteria and archaea."
    ;;
esac

case "${HISAT2_DTA}" in
  0|1)
    ;;
  *)
    die "Invalid HISAT2_DTA: ${HISAT2_DTA}" \
        "Use HISAT2_DTA=1 to enable --dta for eukaryotic RNA-seq, or HISAT2_DTA=0 to disable it."
    ;;
esac

if ! is_nonnegative_int "${HISAT2_MIN_INTRONLEN}" || ! is_nonnegative_int "${HISAT2_MAX_INTRONLEN}"; then
  die "Invalid intron length setting: min=${HISAT2_MIN_INTRONLEN}, max=${HISAT2_MAX_INTRONLEN}" \
      "Set HISAT2_MIN_INTRONLEN and HISAT2_MAX_INTRONLEN to non-negative integers."
fi

if (( HISAT2_MIN_INTRONLEN > HISAT2_MAX_INTRONLEN )); then
  die "HISAT2_MIN_INTRONLEN is greater than HISAT2_MAX_INTRONLEN" \
      "Check the organism-specific intron length settings in config/config.sh."
fi

for value_name in ALIGN_PARALLEL_JOBS HISAT2_THREADS_PER_SAMPLE SAMTOOLS_SORT_THREADS_PER_SAMPLE SAMTOOLS_INDEX_THREADS; do
  value="${!value_name}"
  if ! is_positive_int "${value}"; then
    die "Invalid ${value_name}: ${value}" \
        "Set ${value_name} to a positive integer in config/config.sh."
  fi
done

case "${RESUME_MODE}" in
  0|1)
    ;;
  *)
    die "Invalid RESUME_MODE: ${RESUME_MODE}" \
        "Use RESUME_MODE=1 to reuse completed outputs, or RESUME_MODE=0 to recompute."
    ;;
esac

if (( ALIGN_PARALLEL_JOBS * (HISAT2_THREADS_PER_SAMPLE + SAMTOOLS_SORT_THREADS_PER_SAMPLE) > THREADS )); then
  echo "[WARN] alignment configured CPU use exceeds THREADS: ${ALIGN_PARALLEL_JOBS} jobs * (${HISAT2_THREADS_PER_SAMPLE} HISAT2 + ${SAMTOOLS_SORT_THREADS_PER_SAMPLE} sort threads) > THREADS=${THREADS}" >&2
  echo "[WARN] This can be valid on large machines, but may cause slow alignment or high memory pressure." >&2
fi

case "${HISAT2_RNA_STRANDNESS}" in
  ""|FR|RF)
    ;;
  *)
    die "Invalid HISAT2_RNA_STRANDNESS for paired-end data: ${HISAT2_RNA_STRANDNESS}" \
        "Use an empty value for unstranded libraries, FR when R1 is sense to transcripts, or RF when R1 is antisense to transcripts."
    ;;
esac

splice_sites="$(dirname "${HISAT2_INDEX_PREFIX}")/splice_sites.txt"
hisat2_common_args=(
  -p "${HISAT2_THREADS_PER_SAMPLE}"
  -x "${HISAT2_INDEX_PREFIX}"
)

if [[ "${HISAT2_SPLICE_MODE}" == "no" ]]; then
  hisat2_common_args+=(--no-spliced-alignment)
else
  hisat2_common_args+=(--min-intronlen "${HISAT2_MIN_INTRONLEN}")
  hisat2_common_args+=(--max-intronlen "${HISAT2_MAX_INTRONLEN}")
fi

if [[ "${HISAT2_DTA}" == "1" && "${HISAT2_SPLICE_MODE}" != "no" ]]; then
  hisat2_common_args+=(--dta)
fi

if [[ -n "${HISAT2_RNA_STRANDNESS}" ]]; then
  hisat2_common_args+=(--rna-strandness "${HISAT2_RNA_STRANDNESS}")
fi

if [[ "${HISAT2_SPLICE_MODE}" != "no" && -s "${splice_sites}" ]]; then
  known_splice_arg=(--known-splicesite-infile "${splice_sites}")
  hisat2_common_args+=("${known_splice_arg[@]}")
elif [[ "${HISAT2_SPLICE_MODE}" == "yes" ]]; then
  echo "[WARN] Known splice-site file is empty or missing: ${splice_sites}" >&2
  echo "[WARN] HISAT2 will still run splice-aware alignment without known splice-site hints." >&2
fi

if [[ -n "${HISAT2_EXTRA_ARGS}" ]]; then
  read -r -a hisat2_extra_args <<< "${HISAT2_EXTRA_ARGS}"
  hisat2_common_args+=("${hisat2_extra_args[@]}")
fi

sample_count=0
job_count=0
failed_jobs=0
sample_table="$(mktemp "${LOG_DIR}/hisat2_samples.XXXXXX.tsv")"
trap 'rm -f "${sample_table}"' EXIT

while IFS=$'\t' read -r sample_id condition fastq_1 fastq_2 || [[ -n "${sample_id:-}" ]]; do
  if [[ -z "${sample_id:-}" ]]; then
    continue
  fi

  sample_id="${sample_id%$'\r'}"
  condition="${condition%$'\r'}"
  fastq_1="${fastq_1%$'\r'}"
  fastq_2="${fastq_2%$'\r'}"

  clean_r1="${CLEAN_DIR}/${sample_id}_clean_R1.fastq.gz"
  clean_r2="${CLEAN_DIR}/${sample_id}_clean_R2.fastq.gz"
  sorted_bam="${BAM_DIR}/${sample_id}.sorted.bam"
  hisat2_log="${LOG_DIR}/${sample_id}.hisat2.log"

  if [[ ! -f "${clean_r1}" ]]; then
    die "Clean R1 FASTQ not found for sample ${sample_id}: ${clean_r1}" \
        "Run bash scripts/01_fastp_qc.sh first, or check sample_id in metadata."
  fi

  if [[ ! -f "${clean_r2}" ]]; then
    die "Clean R2 FASTQ not found for sample ${sample_id}: ${clean_r2}" \
        "Run bash scripts/01_fastp_qc.sh first, or check sample_id in metadata."
  fi

  sample_count=$((sample_count + 1))
  printf "%s\t%s\t%s\n" "${sample_id}" "${clean_r1}" "${clean_r2}" >> "${sample_table}"
done < <(tail -n +2 "${METADATA}")

if (( sample_count == 0 )); then
  die "No samples found in metadata" \
      "Add at least one sample row below the header in ${METADATA}."
fi

run_hisat2_sample() {
  local sample_id="$1"
  local clean_r1="$2"
  local clean_r2="$3"
  local sorted_bam hisat2_log

  sorted_bam="${BAM_DIR}/${sample_id}.sorted.bam"
  hisat2_log="${LOG_DIR}/${sample_id}.hisat2.log"

  if [[ "${RESUME_MODE}" == "1" &&
        -s "${sorted_bam}" &&
        -s "${sorted_bam}.bai" &&
        -s "${hisat2_log}" ]]; then
    echo "[INFO] HISAT2 skip completed sample: ${sample_id}"
    return 0
  fi

  echo "[INFO] HISAT2: ${sample_id}"
  hisat2 \
    "${hisat2_common_args[@]}" \
    -1 "${clean_r1}" \
    -2 "${clean_r2}" \
    2> "${hisat2_log}" \
    | samtools sort \
        -@ "${SAMTOOLS_SORT_THREADS_PER_SAMPLE}" \
        -o "${sorted_bam}" \
        -

  samtools index -@ "${SAMTOOLS_INDEX_THREADS}" "${sorted_bam}"
}

echo "[INFO] HISAT2 parallel jobs: ${ALIGN_PARALLEL_JOBS}"
echo "[INFO] HISAT2 threads per sample: ${HISAT2_THREADS_PER_SAMPLE}"
echo "[INFO] samtools sort threads per sample: ${SAMTOOLS_SORT_THREADS_PER_SAMPLE}"

while IFS=$'\t' read -r sample_id clean_r1 clean_r2 || [[ -n "${sample_id:-}" ]]; do
  run_hisat2_sample "${sample_id}" "${clean_r1}" "${clean_r2}" &
  job_count=$((job_count + 1))

  if (( job_count >= ALIGN_PARALLEL_JOBS )); then
    if ! wait -n; then
      failed_jobs=$((failed_jobs + 1))
    fi
    job_count=$((job_count - 1))
  fi
done < "${sample_table}"

while (( job_count > 0 )); do
  if ! wait -n; then
    failed_jobs=$((failed_jobs + 1))
  fi
  job_count=$((job_count - 1))
done

if (( failed_jobs > 0 )); then
  die "HISAT2 failed for ${failed_jobs} sample job(s)" \
      "Open logs/*.hisat2.log to find the failed sample and exact HISAT2 error."
fi

grep -H "overall alignment rate" "${LOG_DIR}"/*.hisat2.log \
  > "${LOG_DIR}/hisat2_alignment_rate_summary.txt" || true

echo "[INFO] HISAT2 alignment finished"
