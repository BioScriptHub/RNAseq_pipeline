#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check FASTQ filenames, 09_metadata/contrasts.csv, and the log above." >&2' ERR

# Generate metadata.tsv and contrasts.tsv.
# Beginner mode:
#   User only edits 09_metadata/contrasts.csv.
#   Sample IDs and conditions are inferred from FASTQ filenames.
# Fallback mode:
#   If filenames cannot encode conditions clearly, user may provide 09_metadata/samples.csv.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

source config/config.sh

SAMPLES_CSV="${SAMPLES_CSV:-09_metadata/samples.csv}"
CONTRASTS_CSV="${CONTRASTS_CSV:-09_metadata/contrasts.csv}"
METADATA_OUT="${METADATA:-09_metadata/metadata.tsv}"
CONTRASTS_OUT="${CONTRASTS:-09_metadata/contrasts.tsv}"
MATCH_LOG="${LOG_DIR}/metadata_fastq_match.tsv"

die() {
  echo "[ERROR] $1" >&2
  echo "[HINT] $2" >&2
  exit 1
}

trim_cr() {
  tr -d '\r'
}

strip_fastq_ext() {
  local name="$1"
  case "${name}" in
    *.fastq.gz) name="${name%.fastq.gz}" ;;
    *.fq.gz) name="${name%.fq.gz}" ;;
    *.fastq) name="${name%.fastq}" ;;
    *.fq) name="${name%.fq}" ;;
  esac
  printf "%s" "${name}"
}

strip_technical_suffix() {
  local sample="$1"
  local changed=1

  shopt -s nocasematch
  while (( changed == 1 )); do
    changed=0
    if [[ "${sample}" =~ ^(.+)[._-](dedup|clean|trimmed|trim|filtered|filter)$ ]]; then
      sample="${BASH_REMATCH[1]}"
      changed=1
    fi
  done
  shopt -u nocasematch

  printf "%s" "${sample}"
}

parse_fastq_name() {
  local fastq="$1"
  local name stem sample read suffix

  name="$(basename "${fastq}")"
  stem="$(strip_fastq_ext "${name}")"

  shopt -s nocasematch
  if [[ "${stem}" =~ ^(.+)[._-](R1|READ1|READ_1)([._-].*)?$ ]]; then
    sample="${BASH_REMATCH[1]}"
    read="R1"
  elif [[ "${stem}" =~ ^(.+)[._-](R2|READ2|READ_2)([._-].*)?$ ]]; then
    sample="${BASH_REMATCH[1]}"
    read="R2"
  elif [[ "${stem}" =~ ^(.+)[._-]1([._-].*)?$ ]]; then
    sample="${BASH_REMATCH[1]}"
    read="R1"
  elif [[ "${stem}" =~ ^(.+)[._-]2([._-].*)?$ ]]; then
    sample="${BASH_REMATCH[1]}"
    read="R2"
  else
    shopt -u nocasematch
    die "Cannot detect R1/R2 from FASTQ filename: ${name}" \
        "Use filenames such as Sample_1_R1.fastq.gz and Sample_1_R2.fastq.gz, or provide 09_metadata/samples.csv."
  fi
  shopt -u nocasematch

  sample="$(strip_technical_suffix "${sample}")"
  if [[ -z "${sample}" ]]; then
    die "Empty sample_id inferred from FASTQ filename: ${name}" \
        "Rename FASTQ files so the sample name appears before R1/R2."
  fi

  if [[ ! "${sample}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
    die "Invalid sample_id inferred from FASTQ filename: ${sample}" \
        "Use only letters, numbers, underscore, dot and hyphen in FASTQ sample names."
  fi

  printf "%s\t%s\t%s\n" "${sample}" "${read}" "${fastq}"
}

infer_condition() {
  local sample="$1"
  local condition

  condition="${sample}"
  if [[ "${sample}" =~ ^(.+)[._-](rep)?[0-9]+$ ]]; then
    condition="${BASH_REMATCH[1]}"
  fi

  if [[ -z "${condition}" || "${condition}" == "${sample}" ]]; then
    printf "%s" "${condition}"
    return 0
  fi

  printf "%s" "${condition}"
}

build_fastq_index() {
  local index_file="$1"
  local sample read fastq key parsed
  declare -A seen=()

  : > "${index_file}"

  while IFS= read -r -d '' fastq; do
    parsed="$(parse_fastq_name "${fastq}")"
    sample="$(printf "%s" "${parsed}" | cut -f1)"
    read="$(printf "%s" "${parsed}" | cut -f2)"
    key="${sample}:${read}"

    if [[ -n "${seen[${key}]:-}" ]]; then
      {
        echo "Duplicate ${read} FASTQ inferred for sample ${sample}:"
        echo "  ${seen[${key}]}"
        echo "  ${fastq}"
      } >&2
      die "Ambiguous FASTQ pairing for sample ${sample}" \
          "Rename FASTQ files so each sample has exactly one R1 and one R2, or provide 09_metadata/samples.csv."
    fi

    seen["${key}"]="${fastq}"
    printf "%s\n" "${parsed}" >> "${index_file}"
  done < <(
    find "${RAW_DIR}" -type f \( \
      -name "*.fastq" -o \
      -name "*.fq" -o \
      -name "*.fastq.gz" -o \
      -name "*.fq.gz" \
    \) -print0
  )

  if [[ ! -s "${index_file}" ]]; then
    die "No FASTQ files found under ${RAW_DIR}" \
        "Put FASTQ files into ${RAW_DIR}, or run bash scripts/00_unpack_rawdata.sh first."
  fi
}

lookup_fastq() {
  local index_file="$1"
  local sample="$2"
  local read="$3"
  local match

  match="$(awk -F'\t' -v s="${sample}" -v r="${read}" '$1 == s && $2 == r {print $3}' "${index_file}")"
  if [[ -z "${match}" ]]; then
    die "No ${read} FASTQ found for sample ${sample}" \
        "Check FASTQ filenames under ${RAW_DIR}, or provide 09_metadata/samples.csv."
  fi

  printf "%s" "${match}"
}

generate_metadata_from_samples_csv() {
  local index_file="$1"
  local sample_id condition extra fastq_1 fastq_2
  local tmp_out

  samples_header="$(head -n 1 "${SAMPLES_CSV}" | trim_cr)"
  if [[ "${samples_header}" != "sample_id,condition" ]]; then
    die "Invalid samples.csv header: ${samples_header}" \
        "The first line must be exactly: sample_id,condition."
  fi

  tmp_out="$(mktemp "${METADATA_OUT}.tmp.XXXXXX")"
  {
    echo -e "sample_id\tcondition\tfastq_1\tfastq_2"
    tail -n +2 "${SAMPLES_CSV}" | trim_cr | while IFS=',' read -r sample_id condition extra || [[ -n "${sample_id:-}" ]]; do
      [[ -z "${sample_id:-}" ]] && continue

      if [[ -n "${extra:-}" ]]; then
        die "Too many columns in samples.csv for sample ${sample_id}" \
            "samples.csv must have exactly two columns: sample_id,condition."
      fi
      if [[ -z "${sample_id}" || -z "${condition}" ]]; then
        die "Empty field found in samples.csv" \
            "Every row must contain sample_id and condition."
      fi

      fastq_1="$(lookup_fastq "${index_file}" "${sample_id}" "R1")"
      fastq_2="$(lookup_fastq "${index_file}" "${sample_id}" "R2")"
      printf "%s\t%s\t%s\t%s\n" "${sample_id}" "${condition}" "${fastq_1}" "${fastq_2}"
    done
  } > "${tmp_out}"
  mv "${tmp_out}" "${METADATA_OUT}"
}

generate_metadata_automatically() {
  local index_file="$1"
  local sample fastq_1 fastq_2 condition
  local tmp_out

  tmp_out="$(mktemp "${METADATA_OUT}.tmp.XXXXXX")"
  {
    echo -e "sample_id\tcondition\tfastq_1\tfastq_2"
    cut -f1 "${index_file}" | sort -u | while IFS= read -r sample; do
      fastq_1="$(lookup_fastq "${index_file}" "${sample}" "R1")"
      fastq_2="$(lookup_fastq "${index_file}" "${sample}" "R2")"
      condition="$(infer_condition "${sample}")"

      printf "%s\t%s\t%s\t%s\n" "${sample}" "${condition}" "${fastq_1}" "${fastq_2}"
    done
  } > "${tmp_out}"
  mv "${tmp_out}" "${METADATA_OUT}"
}

prepare_contrasts() {
  local tmp_out

  if [[ -f "${CONTRASTS_CSV}" ]]; then
    contrasts_header="$(head -n 1 "${CONTRASTS_CSV}" | trim_cr)"
    if [[ "${contrasts_header}" != "comparison,numerator,denominator" ]]; then
      die "Invalid contrasts.csv header: ${contrasts_header}" \
          "The first line must be exactly: comparison,numerator,denominator."
    fi

    tmp_out="$(mktemp "${CONTRASTS_OUT}.tmp.XXXXXX")"
    {
      echo -e "comparison\tnumerator\tdenominator"
      tail -n +2 "${CONTRASTS_CSV}" | trim_cr | while IFS=',' read -r comparison numerator denominator extra || [[ -n "${comparison:-}" ]]; do
        [[ -z "${comparison:-}" ]] && continue

        if [[ -n "${extra:-}" ]]; then
          die "Too many columns in contrasts.csv for comparison ${comparison}" \
              "contrasts.csv must have exactly three columns: comparison,numerator,denominator."
        fi
        if [[ -z "${comparison}" || -z "${numerator}" || -z "${denominator}" ]]; then
          die "Empty field found in contrasts.csv" \
              "Every row must contain comparison, numerator and denominator."
        fi
        if [[ ! "${comparison}" =~ ^[A-Za-z0-9_.-]+$ ]]; then
          die "Invalid comparison name in contrasts.csv: ${comparison}" \
              "Use only letters, numbers, underscore, dot and hyphen."
        fi

        printf "%s\t%s\t%s\n" "${comparison}" "${numerator}" "${denominator}"
      done
    } > "${tmp_out}"
    mv "${tmp_out}" "${CONTRASTS_OUT}"
  elif [[ -f "${CONTRASTS_OUT}" ]]; then
    true
  else
    die "Contrast file not found" \
        "Copy 09_metadata/contrasts.example.csv to 09_metadata/contrasts.csv, then edit the comparison rows."
  fi

  contrast_count="$(tail -n +2 "${CONTRASTS_OUT}" | awk 'NF > 0 {count++} END {print count + 0}')"
  if (( contrast_count == 0 )); then
    die "No comparison rows were found in ${CONTRASTS_OUT}" \
        "Add at least one comparison row to 09_metadata/contrasts.csv."
  fi
}

validate_metadata_against_contrasts() {
  local metadata_conditions contrast_conditions missing_conditions contrast_count filtered_count

  metadata_conditions="$(mktemp)"
  contrast_conditions="$(mktemp)"

  tail -n +2 "${METADATA_OUT}" | cut -f2 | sort -u > "${metadata_conditions}"

  contrast_count="$(tail -n +2 "${CONTRASTS_OUT}" | awk 'NF > 0 {count++} END {print count + 0}')"
  if (( contrast_count > 1 )); then
    if grep -Fxq $'Treatment_vs_Control\tTreatment\tControl' "${CONTRASTS_OUT}" &&
       ! grep -Fxq "Treatment" "${metadata_conditions}" &&
       ! grep -Fxq "Control" "${metadata_conditions}"; then
      awk -F'\t' 'BEGIN {OFS = FS} NR == 1 || !($1 == "Treatment_vs_Control" && $2 == "Treatment" && $3 == "Control")' \
        "${CONTRASTS_OUT}" > "${CONTRASTS_OUT}.tmp"
      mv "${CONTRASTS_OUT}.tmp" "${CONTRASTS_OUT}"
      echo "[WARN] Removed leftover example contrast: Treatment_vs_Control,Treatment,Control" >&2
      echo "[WARN] This row is ignored because no matching Treatment/Control samples exist." >&2
    fi
  fi

  filtered_count="$(tail -n +2 "${CONTRASTS_OUT}" | awk 'NF > 0 {count++} END {print count + 0}')"
  if (( filtered_count == 0 )); then
    rm -f "${metadata_conditions}" "${contrast_conditions}"
    die "No valid comparison rows remain after removing example rows" \
        "Edit 09_metadata/contrasts.csv. Keep the header and add real comparisons such as CGA_vs_WT,CGA,WT."
  fi

  tail -n +2 "${CONTRASTS_OUT}" | awk -F'\t' '{print $2 "\n" $3}' | sort -u > "${contrast_conditions}"

  missing_conditions="$(
    comm -13 "${metadata_conditions}" "${contrast_conditions}" | paste -sd "," - || true
  )"

  rm -f "${metadata_conditions}" "${contrast_conditions}"

  if [[ -n "${missing_conditions}" ]]; then
    die "Conditions in contrasts are not present in inferred metadata: ${missing_conditions}" \
        "Use FASTQ names like Condition_1_R1.fastq.gz, or provide 09_metadata/samples.csv to define sample_id and condition explicitly."
  fi
}

write_match_log() {
  local tmp_out

  tmp_out="$(mktemp "${MATCH_LOG}.tmp.XXXXXX")"
  {
    echo -e "sample_id\tcondition\tfastq_1\tfastq_2"
    tail -n +2 "${METADATA_OUT}"
  } > "${tmp_out}"
  mv "${tmp_out}" "${MATCH_LOG}"
}

clear_derived_files_when_csv_is_source() {
  if [[ -f "${CONTRASTS_CSV}" ]]; then
    rm -f "${METADATA_OUT}" "${CONTRASTS_OUT}" "${MATCH_LOG}"
    echo "[INFO] ${CONTRASTS_CSV} found. Rebuilding derived TSV files from scratch."
  fi
}

main() {
  local sample_count

  mkdir -p "$(dirname "${METADATA_OUT}")" "${LOG_DIR}"
  FASTQ_INDEX_FILE="$(mktemp)"
  trap 'rm -f "${FASTQ_INDEX_FILE}"' EXIT

  clear_derived_files_when_csv_is_source
  build_fastq_index "${FASTQ_INDEX_FILE}"
  prepare_contrasts

  if [[ -f "${SAMPLES_CSV}" ]]; then
    generate_metadata_from_samples_csv "${FASTQ_INDEX_FILE}"
    echo "[INFO] Generated ${METADATA_OUT} from ${SAMPLES_CSV}"
  else
    generate_metadata_automatically "${FASTQ_INDEX_FILE}"
    echo "[INFO] Generated ${METADATA_OUT} automatically from FASTQ filenames"
  fi

  sample_count="$(tail -n +2 "${METADATA_OUT}" | awk 'NF > 0 {count++} END {print count + 0}')"
  if (( sample_count == 0 )); then
    die "No sample rows were generated in ${METADATA_OUT}" \
        "Check FASTQ filenames under ${RAW_DIR}."
  fi

  validate_metadata_against_contrasts
  write_match_log

  echo "[INFO] Generated ${CONTRASTS_OUT}"
  echo "[INFO] FASTQ matching table: ${MATCH_LOG}"
}

main "$@"
