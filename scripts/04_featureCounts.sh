#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check featureCounts installation, BAM files, annotation GTF, STRANDNESS setting, disk space, and the log above." >&2' ERR

# Count reads at gene level using featureCounts.
# Input:  sorted BAM files and annotation GTF
# Output: gene-level raw count matrix

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

is_positive_int() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

detect_group_attribute() {
  local annotation_file="$1"
  local feature_type="$2"
  local candidate
  local feature_rows

  feature_rows="$(awk -F'\t' -v t="${feature_type}" '$0 !~ /^#/ && $3 == t {print $9; count++; if (count >= 1000) exit}' "${annotation_file}")"
  if [[ -z "${feature_rows}" ]]; then
    die "No ${feature_type} records found in annotation file: ${annotation_file}" \
        "Check FEATURE_TYPE in config/config.sh. Common values are exon for eukaryotes and CDS for prokaryotes."
  fi

  for candidate in gene_id gene locus_tag Parent ID; do
    if awk -v key="${candidate}" '
      {
        n = split($0, fields, ";")
        for (i = 1; i <= n; i++) {
          field = fields[i]
          sub(/^[[:space:]]+/, "", field)
          if (field ~ ("^" key "[ =]")) {
            found = 1
            exit
          }
        }
      }
      END {exit(found ? 0 : 1)}
    ' <<< "${feature_rows}"; then
      printf "%s" "${candidate}"
      return 0
    fi
  done

  die "Cannot auto-detect a gene grouping attribute from annotation file" \
      "Set GROUP_ATTRIBUTE manually in config/config.sh. For GTF this is often gene_id. For NCBI GFF3, try gene or locus_tag."
}

mkdir -p "${COUNT_DIR}" "${LOG_DIR}"

resolve_reference_files
echo "[INFO] Annotation file: ${ANNOTATION_GTF}"

if [[ "${GROUP_ATTRIBUTE}" == "auto" ]]; then
  GROUP_ATTRIBUTE="$(detect_group_attribute "${ANNOTATION_GTF}" "${FEATURE_TYPE}")"
  echo "[INFO] Auto-detected GROUP_ATTRIBUTE: ${GROUP_ATTRIBUTE}"
fi

if ! command -v featureCounts >/dev/null 2>&1; then
  die "featureCounts command not found" \
      "Activate the conda environment or install subread before running this step."
fi

case "${STRANDNESS}" in
  0|1|2)
    ;;
  *)
    die "Invalid STRANDNESS: ${STRANDNESS}" \
        "Use 0 for unstranded, 1 for stranded, or 2 for reversely stranded libraries."
    ;;
esac

for flag_name in FEATURECOUNTS_COUNT_READ_PAIRS FEATURECOUNTS_REQUIRE_BOTH_ENDS FEATURECOUNTS_CHECK_CHIMERA; do
  flag_value="${!flag_name}"
  case "${flag_value}" in
    0|1)
      ;;
    *)
      die "Invalid ${flag_name}: ${flag_value}" \
          "Use 1 to enable this featureCounts option or 0 to disable it."
      ;;
  esac
done

if ! is_positive_int "${FEATURECOUNTS_THREADS}"; then
  die "Invalid FEATURECOUNTS_THREADS: ${FEATURECOUNTS_THREADS}" \
      "Set FEATURECOUNTS_THREADS to a positive integer in config/config.sh."
fi

bam_files=()
while IFS=$'\t' read -r sample_id condition fastq_1 fastq_2 || [[ -n "${sample_id:-}" ]]; do
  if [[ -z "${sample_id:-}" ]]; then
    continue
  fi

  sample_id="${sample_id%$'\r'}"
  condition="${condition%$'\r'}"
  fastq_1="${fastq_1%$'\r'}"
  fastq_2="${fastq_2%$'\r'}"

  bam="${BAM_DIR}/${sample_id}.sorted.bam"
  if [[ ! -f "${bam}" ]]; then
    die "BAM not found for sample ${sample_id}: ${bam}" \
        "Run bash scripts/02_hisat2_align.sh first, or check sample_id consistency in metadata."
  fi
  bam_files+=("${bam}")
done < <(tail -n +2 "${METADATA}")

if (( ${#bam_files[@]} == 0 )); then
  die "No BAM files collected from metadata" \
      "Check ${METADATA}. It must contain at least one sample row."
fi

out_counts="${COUNT_DIR}/gene_counts.txt"

featurecounts_common_args=(
  -T "${FEATURECOUNTS_THREADS}"
  -p
)

if [[ "${FEATURECOUNTS_COUNT_READ_PAIRS}" == "1" ]]; then
  featurecounts_help="$(featureCounts -h 2>&1 || true)"
  if printf "%s\n" "${featurecounts_help}" | grep -q -- "--countReadPairs"; then
    featurecounts_common_args+=(--countReadPairs)
  else
    echo "[WARN] featureCounts does not list --countReadPairs in help output." >&2
    echo "[WARN] Continuing with -p only. For Subread 2.x, consider upgrading if fragment counting is required." >&2
  fi
fi

if [[ "${FEATURECOUNTS_REQUIRE_BOTH_ENDS}" == "1" ]]; then
  featurecounts_common_args+=(-B)
fi

if [[ "${FEATURECOUNTS_CHECK_CHIMERA}" == "1" ]]; then
  featurecounts_common_args+=(-C)
fi

featurecounts_common_args+=(
  -s "${STRANDNESS}"
  -t "${FEATURE_TYPE}"
  -g "${GROUP_ATTRIBUTE}"
)

if [[ -n "${FEATURECOUNTS_EXTRA_ARGS}" ]]; then
  read -r -a featurecounts_extra_args <<< "${FEATURECOUNTS_EXTRA_ARGS}"
  featurecounts_common_args+=("${featurecounts_extra_args[@]}")
fi

featureCounts \
  "${featurecounts_common_args[@]}" \
  -a "${ANNOTATION_GTF}" \
  -o "${out_counts}" \
  "${bam_files[@]}" \
  > "${LOG_DIR}/featureCounts.log" 2>&1

echo "[INFO] featureCounts finished: ${out_counts}"
