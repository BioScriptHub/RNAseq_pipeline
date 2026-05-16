#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check metadata, contrasts, reference files, FASTQ paths, software environment, and the log above." >&2' ERR

# Validate required inputs before the pipeline starts expensive steps.
# Input:  config, metadata, contrasts, raw FASTQ, reference FASTA, annotation GTF
# Output: logs/preflight_check.tsv and logs/pipeline_handoff_plan.tsv

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

source config/config.sh
source scripts/lib_reference.sh

die() {
  echo "[ERROR] $1" >&2
  echo "[HINT] $2" >&2
  exit 1
}

require_file() {
  local file="$1"
  local hint="$2"

  if [[ ! -f "${file}" ]]; then
    die "Required file not found: ${file}" "${hint}"
  fi
}

require_command() {
  local command_name="$1"

  if ! command -v "${command_name}" >/dev/null 2>&1; then
    die "Required command not found: ${command_name}" \
        "Activate the conda environment with: conda activate rnaseq_pipeline."
  fi
}

check_tsv_field_count() {
  local file="$1"
  local expected_fields="$2"
  local bad_rows

  bad_rows="$(awk -F'\t' -v expected="${expected_fields}" 'NR > 1 && NF != expected {print NR ":" NF}' "${file}" | paste -sd "," - || true)"
  if [[ -n "${bad_rows}" ]]; then
    die "Invalid column count in ${file}: ${bad_rows}" \
        "Use tab-separated TSV files. Do not use comma-separated CSV files or add extra columns."
  fi
}

mkdir -p "${LOG_DIR}"

require_file "config/config.sh" "Restore config/config.sh from the template."
require_file "${METADATA}" "Create it from 09_metadata/metadata.example.tsv, then edit sample paths."
require_file "${CONTRASTS}" "Create it from 09_metadata/contrasts.example.tsv, then edit comparisons."
resolve_reference_files
echo "[INFO] Genome FASTA: ${GENOME_FA}"
echo "[INFO] Annotation file: ${ANNOTATION_GTF}"

for command_name in fastp fastqc multiqc hisat2 hisat2-build samtools featureCounts Rscript; do
  require_command "${command_name}"
done

metadata_header="$(head -n 1 "${METADATA}" | tr -d '\r')"
expected_metadata_header=$'sample_id\tcondition\tfastq_1\tfastq_2'
if [[ "${metadata_header}" != "${expected_metadata_header}" ]]; then
  die "Metadata header is invalid: ${metadata_header}" \
      "The first line must be exactly: sample_id<TAB>condition<TAB>fastq_1<TAB>fastq_2."
fi
check_tsv_field_count "${METADATA}" 4

contrasts_header="$(head -n 1 "${CONTRASTS}" | tr -d '\r')"
expected_contrasts_header=$'comparison\tnumerator\tdenominator'
if [[ "${contrasts_header}" != "${expected_contrasts_header}" ]]; then
  die "Contrasts header is invalid: ${contrasts_header}" \
      "The first line must be exactly: comparison<TAB>numerator<TAB>denominator."
fi
check_tsv_field_count "${CONTRASTS}" 3

sample_rows="$(tail -n +2 "${METADATA}" | awk 'NF > 0 {count++} END {print count + 0}')"
if (( sample_rows == 0 )); then
  die "metadata.tsv contains no sample rows" \
      "Add at least two biological replicates per condition. Three or more are recommended."
fi

contrast_rows="$(tail -n +2 "${CONTRASTS}" | awk 'NF > 0 {count++} END {print count + 0}')"
if (( contrast_rows == 0 )); then
  die "contrasts.tsv contains no comparison rows" \
      "Add at least one row such as Treatment_vs_Control<TAB>Treatment<TAB>Control."
fi

duplicate_samples="$(tail -n +2 "${METADATA}" | cut -f1 | tr -d '\r' | sort | uniq -d | paste -sd "," - || true)"
if [[ -n "${duplicate_samples}" ]]; then
  die "Duplicate sample_id found: ${duplicate_samples}" \
      "Each sample_id must be unique."
fi

duplicate_comparisons="$(tail -n +2 "${CONTRASTS}" | cut -f1 | tr -d '\r' | sort | uniq -d | paste -sd "," - || true)"
if [[ -n "${duplicate_comparisons}" ]]; then
  die "Duplicate comparison name found: ${duplicate_comparisons}" \
      "Each comparison name must be unique because it is used as an output filename prefix."
fi

metadata_conditions="$(mktemp)"
contrast_conditions="$(mktemp)"
trap 'rm -f "${metadata_conditions}" "${contrast_conditions}"' EXIT

{
  echo -e "sample_id\tcondition\traw_r1\traw_r2\tclean_r1\tclean_r2\tsorted_bam"

  while IFS=$'\t' read -r sample_id condition fastq_1 fastq_2 || [[ -n "${sample_id:-}" ]]; do
    [[ -z "${sample_id:-}" ]] && continue

    sample_id="${sample_id%$'\r'}"
    condition="${condition%$'\r'}"
    fastq_1="${fastq_1%$'\r'}"
    fastq_2="${fastq_2%$'\r'}"

    if [[ -z "${sample_id}" || -z "${condition}" || -z "${fastq_1}" || -z "${fastq_2}" ]]; then
      die "Empty field found in metadata.tsv for sample ${sample_id:-UNKNOWN}" \
          "Every row must contain sample_id, condition, fastq_1 and fastq_2."
    fi

    if [[ ! "${sample_id}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
      die "Invalid sample_id: ${sample_id}" \
          "Use only letters, numbers, underscore, dot and hyphen. sample_id is used in output filenames."
    fi

    if [[ ! -f "${fastq_1}" ]]; then
      die "R1 FASTQ not found for sample ${sample_id}: ${fastq_1}" \
          "Check fastq_1 in ${METADATA}. Paths must be relative to the project root."
    fi

    if [[ ! -f "${fastq_2}" ]]; then
      die "R2 FASTQ not found for sample ${sample_id}: ${fastq_2}" \
          "Check fastq_2 in ${METADATA}. Paths must be relative to the project root."
    fi

    if [[ "${fastq_1}" == "${fastq_2}" ]]; then
      die "R1 and R2 are identical for sample ${sample_id}" \
          "Paired-end metadata must point to two different FASTQ files."
    fi

    printf "%s\n" "${condition}" >> "${metadata_conditions}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${sample_id}" \
      "${condition}" \
      "${fastq_1}" \
      "${fastq_2}" \
      "${CLEAN_DIR}/${sample_id}_clean_R1.fastq.gz" \
      "${CLEAN_DIR}/${sample_id}_clean_R2.fastq.gz" \
      "${BAM_DIR}/${sample_id}.sorted.bam"
  done < <(tail -n +2 "${METADATA}")
} > "${LOG_DIR}/pipeline_handoff_plan.tsv"

while IFS=$'\t' read -r comparison numerator denominator || [[ -n "${comparison:-}" ]]; do
  [[ -z "${comparison:-}" ]] && continue

  comparison="${comparison%$'\r'}"
  numerator="${numerator%$'\r'}"
  denominator="${denominator%$'\r'}"

  if [[ -z "${comparison}" || -z "${numerator}" || -z "${denominator}" ]]; then
    die "Empty field found in contrasts.tsv for comparison ${comparison:-UNKNOWN}" \
        "Every row must contain comparison, numerator and denominator."
  fi

  if [[ ! "${comparison}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
    die "Invalid comparison name: ${comparison}" \
        "Use only letters, numbers, underscore, dot and hyphen. comparison is used in output filenames."
  fi

  printf "%s\n%s\n" "${numerator}" "${denominator}" >> "${contrast_conditions}"
done < <(tail -n +2 "${CONTRASTS}")

missing_conditions="$(
  comm -13 \
    <(sort -u "${metadata_conditions}") \
    <(sort -u "${contrast_conditions}") \
    | paste -sd "," - || true
)"
if [[ -n "${missing_conditions}" ]]; then
  die "Conditions in contrasts.tsv but not in metadata.tsv: ${missing_conditions}" \
      "Make numerator and denominator exactly match the condition values in metadata.tsv."
fi

low_rep_conditions="$(
  sort "${metadata_conditions}" \
    | uniq -c \
    | awk '$1 < 2 {print $2}' \
    | paste -sd "," - || true
)"
if [[ -n "${low_rep_conditions}" ]]; then
  die "Too few replicates for condition(s): ${low_rep_conditions}" \
      "Each condition should have at least two biological replicates. Three or more are recommended."
fi

{
  echo -e "item\tstatus"
  echo -e "metadata\tPASS"
  echo -e "contrasts\tPASS"
  echo -e "raw_fastq\tPASS"
  echo -e "reference_genome\tPASS"
  echo -e "annotation_gtf\tPASS"
  echo -e "software_commands\tPASS"
} > "${LOG_DIR}/preflight_check.tsv"

echo "[INFO] Preflight check finished"
echo "[INFO] Handoff plan: ${LOG_DIR}/pipeline_handoff_plan.tsv"
