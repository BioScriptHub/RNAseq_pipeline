#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check samtools/MultiQC installation, BAM files, file permissions, disk space, and the log above." >&2' ERR

# Generate BAM QC reports using samtools flagstat, samtools stats and MultiQC.
# Input:  sorted BAM files
# Output: flagstat, stats, mapping rate table, MultiQC report

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

source config/config.sh

die() {
  echo "[ERROR] $1" >&2
  echo "[HINT] $2" >&2
  exit 1
}

is_positive_int() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

QC_DIR="${BAM_DIR}/qc"
mkdir -p "${QC_DIR}" "${MULTIQC_DIR}" "${LOG_DIR}"

if ! is_positive_int "${BAM_QC_THREADS}"; then
  die "Invalid BAM_QC_THREADS: ${BAM_QC_THREADS}" \
      "Set BAM_QC_THREADS to a positive integer in config/config.sh."
fi

shopt -s nullglob
bam_files=("${BAM_DIR}"/*.sorted.bam)
shopt -u nullglob

if (( ${#bam_files[@]} == 0 )); then
  die "No sorted BAM files found in ${BAM_DIR}" \
      "Run bash scripts/02_hisat2_align.sh first. Expected files like 07_bam/sample.sorted.bam."
fi

: > "${QC_DIR}/mapping_rate_check.body.tsv"

for bam in "${bam_files[@]}"; do
  sample_id="$(basename "${bam}" .sorted.bam)"
  echo "[INFO] BAM QC: ${sample_id}"

  samtools flagstat \
    -@ "${BAM_QC_THREADS}" \
    "${bam}" \
    > "${QC_DIR}/${sample_id}.flagstat.txt"

  samtools stats \
    -@ "${BAM_QC_THREADS}" \
    "${bam}" \
    > "${QC_DIR}/${sample_id}.stats.txt"

  mapped_rate="$(awk -F'[()%]' '/ mapped \(/ {print $2; exit}' "${QC_DIR}/${sample_id}.flagstat.txt")"
  if [[ -z "${mapped_rate}" ]]; then
    mapped_rate="0"
  fi

  status="LOW_MAPPING_RATE"
  if awk -v rate="${mapped_rate}" -v min="${ALIGNMENT_RATE_MIN}" 'BEGIN {exit !(rate + 0 >= min + 0)}'; then
    status="PASS"
  fi

  printf "%s\t%.2f\t%s\n" "${sample_id}" "${mapped_rate}" "${status}" \
    >> "${QC_DIR}/mapping_rate_check.body.tsv"
done

{
  echo -e "sample_id\tmapped_rate_percent\tstatus"
  sort "${QC_DIR}/mapping_rate_check.body.tsv"
} > "${QC_DIR}/mapping_rate_check.tsv"

rm -f "${QC_DIR}/mapping_rate_check.body.tsv"

multiqc \
  --force \
  --outdir "${MULTIQC_DIR}" \
  --filename "multiqc_bam_qc.html" \
  "${QC_DIR}" "${LOG_DIR}"

echo "[INFO] BAM QC finished"
