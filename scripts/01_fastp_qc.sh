#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check fastp/FastQC/MultiQC installation, FASTQ paths, file permissions, disk space, and the log above." >&2' ERR

# Filter raw paired-end FASTQ files using fastp, then run FastQC and MultiQC.
# Input:  metadata.tsv and raw FASTQ files
# Output: clean FASTQ files, fastp reports, FastQC reports, MultiQC report

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

mkdir -p "${CLEAN_DIR}" "${FASTQC_DIR}" "${MULTIQC_DIR}" "${LOG_DIR}"

if [[ ! -f "${METADATA}" ]]; then
  die "Metadata file not found: ${METADATA}" \
      "Run: cp 09_metadata/metadata.example.tsv 09_metadata/metadata.tsv, then edit FASTQ paths."
fi

expected_header=$'sample_id\tcondition\tfastq_1\tfastq_2'
metadata_header="$(head -n 1 "${METADATA}" | tr -d '\r')"
if [[ "${metadata_header}" != "${expected_header}" ]]; then
  die "Metadata header is invalid: ${metadata_header}" \
      "The first line must be exactly: ${expected_header}"
fi

case "${FASTP_DETECT_ADAPTER_FOR_PE}" in
  0|1)
    ;;
  *)
    die "Invalid FASTP_DETECT_ADAPTER_FOR_PE: ${FASTP_DETECT_ADAPTER_FOR_PE}" \
        "Use FASTP_DETECT_ADAPTER_FOR_PE=1 to enable paired-end adapter detection, or 0 to disable it."
    ;;
esac

for value_name in FASTP_QUALIFIED_QUALITY_PHRED FASTP_UNQUALIFIED_PERCENT_LIMIT FASTP_LENGTH_REQUIRED; do
  value="${!value_name}"
  if ! is_nonnegative_int "${value}"; then
    die "Invalid ${value_name}: ${value}" \
        "Set ${value_name} to a non-negative integer in config/config.sh."
  fi
done

for value_name in FASTP_PARALLEL_JOBS FASTP_THREADS_PER_SAMPLE FASTQC_THREADS; do
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

if (( FASTP_PARALLEL_JOBS * FASTP_THREADS_PER_SAMPLE > THREADS )); then
  echo "[WARN] fastp configured CPU use exceeds THREADS: ${FASTP_PARALLEL_JOBS} jobs * ${FASTP_THREADS_PER_SAMPLE} threads > THREADS=${THREADS}" >&2
  echo "[WARN] This can be valid on large machines, but may slow down laptops or external disks." >&2
fi

fastp_common_args=(
  --thread "${FASTP_THREADS_PER_SAMPLE}"
  --qualified_quality_phred "${FASTP_QUALIFIED_QUALITY_PHRED}"
  --unqualified_percent_limit "${FASTP_UNQUALIFIED_PERCENT_LIMIT}"
  --length_required "${FASTP_LENGTH_REQUIRED}"
)

if [[ "${FASTP_DETECT_ADAPTER_FOR_PE}" == "1" ]]; then
  fastp_common_args+=(--detect_adapter_for_pe)
fi

if [[ -n "${FASTP_EXTRA_ARGS}" ]]; then
  read -r -a fastp_extra_args <<< "${FASTP_EXTRA_ARGS}"
  fastp_common_args+=("${fastp_extra_args[@]}")
fi

duplicate_samples="$(tail -n +2 "${METADATA}" | cut -f1 | sort | uniq -d || true)"
if [[ -n "${duplicate_samples}" ]]; then
  die "Duplicate sample_id found: ${duplicate_samples}" \
      "Each sample_id must be unique. Edit the first column of ${METADATA}."
fi

sample_count=0
job_count=0
failed_jobs=0
sample_table="$(mktemp "${LOG_DIR}/fastp_samples.XXXXXX.tsv")"
trap 'rm -f "${sample_table}"' EXIT

while IFS=$'\t' read -r sample_id condition fastq_1 fastq_2 || [[ -n "${sample_id:-}" ]]; do
  if [[ -z "${sample_id:-}" ]]; then
    continue
  fi

  sample_id="${sample_id%$'\r'}"
  condition="${condition%$'\r'}"
  fastq_1="${fastq_1%$'\r'}"
  fastq_2="${fastq_2%$'\r'}"

  sample_count=$((sample_count + 1))

  if [[ ! -f "${fastq_1}" ]]; then
    die "R1 FASTQ not found for sample ${sample_id}: ${fastq_1}" \
        "Check the fastq_1 path in ${METADATA}. Paths are relative to the project root."
  fi

  if [[ ! -f "${fastq_2}" ]]; then
    die "R2 FASTQ not found for sample ${sample_id}: ${fastq_2}" \
        "Check the fastq_2 path in ${METADATA}. Paths are relative to the project root."
  fi

  printf "%s\t%s\t%s\t%s\n" "${sample_id}" "${condition}" "${fastq_1}" "${fastq_2}" >> "${sample_table}"
done < <(tail -n +2 "${METADATA}")

if (( sample_count == 0 )); then
  die "No samples found in metadata" \
      "Add at least one sample row below the header in ${METADATA}."
fi

run_fastp_sample() {
  local sample_id="$1"
  local fastq_1="$2"
  local fastq_2="$3"
  local out_r1 out_r2 html_report json_report log_file

  out_r1="${CLEAN_DIR}/${sample_id}_clean_R1.fastq.gz"
  out_r2="${CLEAN_DIR}/${sample_id}_clean_R2.fastq.gz"
  html_report="${FASTQC_DIR}/${sample_id}.fastp.html"
  json_report="${FASTQC_DIR}/${sample_id}.fastp.json"
  log_file="${LOG_DIR}/${sample_id}.fastp.log"

  if [[ "${RESUME_MODE}" == "1" &&
        -s "${out_r1}" &&
        -s "${out_r2}" &&
        -s "${html_report}" &&
        -s "${json_report}" ]]; then
    echo "[INFO] fastp skip completed sample: ${sample_id}"
    return 0
  fi

  echo "[INFO] fastp: ${sample_id}"
  fastp \
    --in1 "${fastq_1}" \
    --in2 "${fastq_2}" \
    --out1 "${out_r1}" \
    --out2 "${out_r2}" \
    "${fastp_common_args[@]}" \
    --html "${html_report}" \
    --json "${json_report}" \
    > "${log_file}" 2>&1
}

echo "[INFO] fastp parallel jobs: ${FASTP_PARALLEL_JOBS}"
echo "[INFO] fastp threads per sample: ${FASTP_THREADS_PER_SAMPLE}"

while IFS=$'\t' read -r sample_id condition fastq_1 fastq_2 || [[ -n "${sample_id:-}" ]]; do
  run_fastp_sample "${sample_id}" "${fastq_1}" "${fastq_2}" &
  job_count=$((job_count + 1))

  if (( job_count >= FASTP_PARALLEL_JOBS )); then
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
  die "fastp failed for ${failed_jobs} sample job(s)" \
      "Open logs/*.fastp.log to find the failed sample and exact fastp error."
fi

if ! compgen -G "${CLEAN_DIR}/*_clean_R1.fastq.gz" > /dev/null || \
   ! compgen -G "${CLEAN_DIR}/*_clean_R2.fastq.gz" > /dev/null; then
  die "No clean FASTQ files were produced by fastp" \
      "Check logs/*.fastp.log and available disk space."
fi

echo "[INFO] Running FastQC on clean FASTQ files"
fastqc \
  --threads "${FASTQC_THREADS}" \
  --outdir "${FASTQC_DIR}" \
  "${CLEAN_DIR}"/*_clean_R1.fastq.gz \
  "${CLEAN_DIR}"/*_clean_R2.fastq.gz

echo "[INFO] Running MultiQC for FASTQ QC"
multiqc \
  --force \
  --outdir "${MULTIQC_DIR}" \
  --filename "multiqc_fastq_qc.html" \
  "${FASTQC_DIR}" "${LOG_DIR}"

echo "[INFO] FASTQ QC finished"
