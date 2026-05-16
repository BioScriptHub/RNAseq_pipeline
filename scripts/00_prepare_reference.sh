#!/usr/bin/env bash
set -euo pipefail
trap 'echo "[ERROR] Command failed in $(basename "$0") at line ${LINENO}: ${BASH_COMMAND}" >&2; echo "[HINT] Check tool installation, input paths, file permissions, disk space, and the log above." >&2' ERR

# Build HISAT2 index and extract splice site and exon information.
# Input:  config/config.sh, genome FASTA, annotation GTF
# Output: HISAT2 index, splice_sites.txt, exons.txt

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

mkdir -p "$(dirname "${HISAT2_INDEX_PREFIX}")" "${LOG_DIR}"

resolve_reference_files
echo "[INFO] Genome FASTA: ${GENOME_FA}"
echo "[INFO] Annotation file: ${ANNOTATION_GTF}"

case "${HISAT2_SPLICE_MODE}" in
  auto|yes|no)
    ;;
  *)
    die "Invalid HISAT2_SPLICE_MODE: ${HISAT2_SPLICE_MODE}" \
        "Use auto, yes, or no in config/config.sh. Use no for bacteria and archaea."
    ;;
esac

if [[ "${FORCE_REBUILD_INDEX}" != "1" ]]; then
  if [[ -f "${HISAT2_INDEX_PREFIX}.1.ht2" || -f "${HISAT2_INDEX_PREFIX}.1.ht2l" ]]; then
    echo "[INFO] HISAT2 index already exists: ${HISAT2_INDEX_PREFIX}"
    echo "[INFO] Set FORCE_REBUILD_INDEX=1 to rebuild"
    exit 0
  fi
fi

if command -v seqkit >/dev/null 2>&1; then
  seqkit stats "${GENOME_FA}" > "${LOG_DIR}/reference_genome.seqkit_stats.txt"
fi

awk '$3 == "gene" {gene++} $3 == "exon" {exon++} END {print "gene_records\t" gene; print "exon_records\t" exon}' \
  "${ANNOTATION_GTF}" > "${LOG_DIR}/annotation_gtf_record_count.tsv"

splice_out="$(dirname "${HISAT2_INDEX_PREFIX}")/splice_sites.txt"
exon_out="$(dirname "${HISAT2_INDEX_PREFIX}")/exons.txt"

splice_script="$(command -v hisat2_extract_splice_sites.py || true)"
exon_script="$(command -v hisat2_extract_exons.py || true)"

if [[ "${HISAT2_SPLICE_MODE}" == "no" ]]; then
  echo "[INFO] HISAT2_SPLICE_MODE=no. Skipping splice-site and exon extraction."
  : > "${splice_out}"
  : > "${exon_out}"
else
  if [[ -n "${splice_script}" ]]; then
    if ! "${splice_script}" "${ANNOTATION_GTF}" > "${splice_out}"; then
      if [[ "${HISAT2_SPLICE_MODE}" == "yes" ]]; then
        die "Failed to extract splice sites from annotation file: ${ANNOTATION_GTF}" \
            "Use a valid GTF file for strict splice-aware mode, or set HISAT2_SPLICE_MODE=auto to continue without splice-site hints."
      fi
      echo "[WARN] Failed to extract splice sites. Continuing without known splice-site hints." >&2
      : > "${splice_out}"
    fi
  elif [[ "${HISAT2_SPLICE_MODE}" == "yes" ]]; then
    die "hisat2_extract_splice_sites.py not found" \
        "Install HISAT2 correctly or set HISAT2_SPLICE_MODE=auto to continue without known splice-site hints."
  else
    echo "[WARN] hisat2_extract_splice_sites.py not found. Creating empty splice site file." >&2
    : > "${splice_out}"
  fi

  if [[ -n "${exon_script}" ]]; then
    if ! "${exon_script}" "${ANNOTATION_GTF}" > "${exon_out}"; then
      if [[ "${HISAT2_SPLICE_MODE}" == "yes" ]]; then
        die "Failed to extract exons from annotation file: ${ANNOTATION_GTF}" \
            "Use a valid GTF file for strict splice-aware mode, or set HISAT2_SPLICE_MODE=auto to continue without exon hints."
      fi
      echo "[WARN] Failed to extract exons. Continuing without exon hints." >&2
      : > "${exon_out}"
    fi
  elif [[ "${HISAT2_SPLICE_MODE}" == "yes" ]]; then
    die "hisat2_extract_exons.py not found" \
        "Install HISAT2 correctly or set HISAT2_SPLICE_MODE=auto to continue without exon hints."
  else
    echo "[WARN] hisat2_extract_exons.py not found. Creating empty exon file." >&2
    : > "${exon_out}"
  fi
fi

build_args=()
if [[ -s "${splice_out}" ]]; then
  build_args+=(--ss "${splice_out}")
fi
if [[ -s "${exon_out}" ]]; then
  build_args+=(--exon "${exon_out}")
fi

echo "[INFO] Building HISAT2 index"
hisat2-build \
  "${build_args[@]}" \
  "${GENOME_FA}" \
  "${HISAT2_INDEX_PREFIX}" \
  > "${LOG_DIR}/hisat2_build.log" 2>&1

echo "[INFO] Reference preparation finished"
