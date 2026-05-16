#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check archive format, file permissions, disk space, and the log above." >&2' ERR

# Prepare raw RNA-seq files.
# Delivery archives: zip, tar, tar.gz, tgz, tar.bz2, tbz2, tar.xz, txz.
# Final FASTQ files: fastq, fq, fastq.gz, fq.gz. These are never decompressed.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_DIR}"

source config/config.sh

ARCHIVE_DIR="${ARCHIVE_DIR:-00_archives}"
RAW_DIR="${RAW_DIR:-00_rawdata}"
UNPACK_NESTED_ARCHIVES="${UNPACK_NESTED_ARCHIVES:-1}"
UNPACK_MAX_ROUNDS="${UNPACK_MAX_ROUNDS:-5}"

die() {
  echo "[ERROR] $1" >&2
  echo "[HINT] $2" >&2
  exit 1
}

is_positive_int() {
  [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

validate_config() {
  case "${UNPACK_NESTED_ARCHIVES}" in
    0|1)
      ;;
    *)
      die "Invalid UNPACK_NESTED_ARCHIVES: ${UNPACK_NESTED_ARCHIVES}" \
          "Use UNPACK_NESTED_ARCHIVES=1 to unpack nested delivery archives, or 0 to disable it."
      ;;
  esac

  if ! is_positive_int "${UNPACK_MAX_ROUNDS}"; then
    die "Invalid UNPACK_MAX_ROUNDS: ${UNPACK_MAX_ROUNDS}" \
        "Set UNPACK_MAX_ROUNDS to a positive integer, for example 5."
  fi
}

find_delivery_archives() {
  local root="$1"

  [[ -d "${root}" ]] || return 0

  find "${root}" -type f \( \
    -name "*.zip" -o \
    -name "*.tar" -o \
    -name "*.tar.gz" -o \
    -name "*.tgz" -o \
    -name "*.tar.bz2" -o \
    -name "*.tbz2" -o \
    -name "*.tar.xz" -o \
    -name "*.txz" \
  \) -print0
}

find_fastq_files() {
  local root="$1"

  [[ -d "${root}" ]] || return 0

  find "${root}" -type f \( \
    -name "*.fastq" -o \
    -name "*.fq" -o \
    -name "*.fastq.gz" -o \
    -name "*.fq.gz" \
  \) -print0
}

copy_fastq_from_archive_dir() {
  local fastq rel_path out_path copied=0

  while IFS= read -r -d '' fastq; do
    rel_path="${fastq#${ARCHIVE_DIR}/}"
    out_path="${RAW_DIR}/${rel_path}"
    mkdir -p "$(dirname "${out_path}")"
    cp -n "${fastq}" "${out_path}"
    echo "[INFO] FASTQ kept compressed: ${out_path}"
    copied=$((copied + 1))
  done < <(find_fastq_files "${ARCHIVE_DIR}")

  echo "[INFO] FASTQ files copied from ${ARCHIVE_DIR}: ${copied}"
}

unpack_archive() {
  local archive="$1"
  local destination="$2"

  mkdir -p "${destination}"
  echo "[INFO] Unpacking delivery archive: ${archive}"

  case "${archive}" in
    *.zip)
      command -v unzip >/dev/null 2>&1 || \
        die "unzip command not found, cannot unpack ${archive}" \
            "Install unzip or manually extract this archive into ${destination}."
      unzip -n "${archive}" -d "${destination}"
      ;;
    *.tar)
      tar -xf "${archive}" -C "${destination}"
      ;;
    *.tar.gz|*.tgz)
      tar -xzf "${archive}" -C "${destination}"
      ;;
    *.tar.bz2|*.tbz2)
      tar -xjf "${archive}" -C "${destination}"
      ;;
    *.tar.xz|*.txz)
      tar -xJf "${archive}" -C "${destination}"
      ;;
    *)
      die "Unsupported delivery archive format: ${archive}" \
          "Supported delivery archives are zip, tar, tar.gz, tgz, tar.bz2, tbz2, tar.xz, and txz. FASTQ gzip files must end with fastq.gz or fq.gz and are not unpacked."
      ;;
  esac

  printf "%s\t%s\n" "${archive}" "${destination}" >> "${LOG_DIR}/unpacked_delivery_archives.tsv"
}

unpack_archives_from_archive_dir() {
  local archive unpacked=0

  while IFS= read -r -d '' archive; do
    unpack_archive "${archive}" "${RAW_DIR}"
    unpacked=$((unpacked + 1))
  done < <(find_delivery_archives "${ARCHIVE_DIR}")

  echo "[INFO] Delivery archives unpacked from ${ARCHIVE_DIR}: ${unpacked}"
}

unpack_nested_archives_from_raw_dir() {
  local round archive unpacked_this_round total_unpacked=0 remaining_unprocessed=0
  declare -A processed=()

  for ((round = 1; round <= UNPACK_MAX_ROUNDS; round++)); do
    unpacked_this_round=0

    while IFS= read -r -d '' archive; do
      [[ -n "${processed[${archive}]:-}" ]] && continue
      processed["${archive}"]=1
      unpack_archive "${archive}" "$(dirname "${archive}")"
      unpacked_this_round=$((unpacked_this_round + 1))
      total_unpacked=$((total_unpacked + 1))
    done < <(find_delivery_archives "${RAW_DIR}")

    (( unpacked_this_round == 0 )) && break
    echo "[INFO] Nested archive round ${round}: ${unpacked_this_round} archive(s)"
  done

  if (( total_unpacked > 0 )); then
    echo "[INFO] Nested delivery archives unpacked from ${RAW_DIR}: ${total_unpacked}"
  fi

  while IFS= read -r -d '' archive; do
    [[ -n "${processed[${archive}]:-}" ]] && continue
    remaining_unprocessed=$((remaining_unprocessed + 1))
  done < <(find_delivery_archives "${RAW_DIR}")

  if (( remaining_unprocessed > 0 )); then
    die "Nested archive unpacking stopped before all delivery archives were processed" \
        "Increase UNPACK_MAX_ROUNDS in config/config.sh, or manually inspect nested archives under ${RAW_DIR}."
  fi
}

warn_if_nested_archives_remain_disabled() {
  local archive_count=0

  while IFS= read -r -d '' _archive; do
    archive_count=$((archive_count + 1))
  done < <(find_delivery_archives "${RAW_DIR}")

  if (( archive_count > 0 )); then
    echo "[WARN] ${archive_count} delivery archive(s) found under ${RAW_DIR}, but UNPACK_NESTED_ARCHIVES=0." >&2
    echo "[WARN] FASTQ files inside those archives will not be detected." >&2
  fi
}

write_fastq_manifest() {
  local fastq_count=0

  {
    echo -e "fastq_path"
    while IFS= read -r -d '' fastq; do
      printf "%s\n" "${fastq}"
      fastq_count=$((fastq_count + 1))
    done < <(find_fastq_files "${RAW_DIR}" | sort -z)
  } > "${LOG_DIR}/raw_fastq_files_after_unpack.tsv"

  if (( fastq_count == 0 )); then
    die "No FASTQ files found in ${RAW_DIR}" \
        "Put fastq.gz/fq.gz files in ${RAW_DIR}, or put delivery archives in ${ARCHIVE_DIR}. Do not decompress fastq.gz into fastq unless you have a specific reason."
  fi

  echo "[INFO] FASTQ files detected: ${fastq_count}"
  echo "[INFO] FASTQ manifest: ${LOG_DIR}/raw_fastq_files_after_unpack.tsv"
}

main() {
  validate_config
  mkdir -p "${ARCHIVE_DIR}" "${RAW_DIR}" "${LOG_DIR}"

  echo -e "archive_path\tdestination" > "${LOG_DIR}/unpacked_delivery_archives.tsv"
  echo "[INFO] Delivery archives will be unpacked."
  echo "[INFO] FASTQ gzip files (*.fastq.gz, *.fq.gz) are final inputs and will not be decompressed."

  copy_fastq_from_archive_dir
  unpack_archives_from_archive_dir

  if [[ "${UNPACK_NESTED_ARCHIVES}" == "1" ]]; then
    unpack_nested_archives_from_raw_dir
  else
    warn_if_nested_archives_remain_disabled
  fi

  write_fastq_manifest
  echo "[INFO] Delivery archive log: ${LOG_DIR}/unpacked_delivery_archives.tsv"
}

main "$@"
