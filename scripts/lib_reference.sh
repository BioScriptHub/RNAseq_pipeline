#!/usr/bin/env bash

# Resolve reference genome and annotation files by suffix.
# User may keep NCBI filenames such as GCF_xxx_genomic.fna and genomic.gff.

resolve_one_by_suffix() {
  local configured_file="$1"
  local search_dir="$2"
  local suffix_string="$3"
  local label="$4"
  local -a matches=()
  local -a suffixes=()
  local suffix file

  if [[ -f "${configured_file}" ]]; then
    printf "%s" "${configured_file}"
    return 0
  fi

  if [[ ! -d "${search_dir}" ]]; then
    echo "[ERROR] Reference directory not found: ${search_dir}" >&2
    echo "[HINT] Create ${search_dir} and put your ${label} file there." >&2
    return 1
  fi

  read -r -a suffixes <<< "${suffix_string}"
  shopt -s nullglob
  for suffix in "${suffixes[@]}"; do
    for file in "${search_dir}"/${suffix}; do
      [[ -f "${file}" ]] && matches+=("${file}")
    done
  done
  shopt -u nullglob

  if (( ${#matches[@]} == 0 )); then
    echo "[ERROR] No ${label} file found in ${search_dir}" >&2
    echo "[HINT] Expected suffix: ${suffix_string}. Put one matching file in ${search_dir}, or set the exact path in config/config.sh." >&2
    return 1
  fi

  if (( ${#matches[@]} > 1 )); then
    echo "[ERROR] Multiple ${label} files found in ${search_dir}" >&2
    printf "  %s\n" "${matches[@]}" >&2
    echo "[HINT] Keep only one ${label} file in ${search_dir}, or set the exact path in config/config.sh." >&2
    return 1
  fi

  printf "%s" "${matches[0]}"
}

resolve_reference_files() {
  GENOME_FA="$(resolve_one_by_suffix "${GENOME_FA}" "${REFERENCE_DIR}" "${GENOME_SUFFIXES}" "genome FASTA")" || return 1
  ANNOTATION_GTF="$(resolve_one_by_suffix "${ANNOTATION_GTF}" "${REFERENCE_DIR}" "${ANNOTATION_SUFFIXES}" "annotation GTF/GFF")" || return 1
  export GENOME_FA ANNOTATION_GTF
}
