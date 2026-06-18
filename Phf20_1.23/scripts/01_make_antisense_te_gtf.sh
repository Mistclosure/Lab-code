#!/usr/bin/env bash
set -euo pipefail

WORKDIR="/mnt/disk1/qiuzerui/expriments/Phf20_1.23"
ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"

ANALYSIS_DIR="${WORKDIR}/sense_antisense_ERV"
OUT_ANNO_DIR="${ANNO_DIR}/TE_antisense_strand_flipped"
LOG_DIR="${ANALYSIS_DIR}/logs"
ANTISENSE_TE_GTF="${OUT_ANNO_DIR}/m39_TE.antisense.strand_flipped.gtf"
README_FILE="${OUT_ANNO_DIR}/README.txt"
LOG_FILE="${LOG_DIR}/01_make_antisense_te_gtf.log"

mkdir -p "${OUT_ANNO_DIR}" "${LOG_DIR}"

{
  echo "[$(date)] Start making antisense TE GTF"
  echo "Input GTF_TE: ${GTF_TE}"
  echo "Antisense GTF: ${ANTISENSE_TE_GTF}"

  if [ ! -s "${GTF_TE}" ]; then
    echo "ERROR: TE GTF not found or empty: ${GTF_TE}" >&2
    exit 1
  fi

  awk 'BEGIN{FS=OFS="\t"}
       /^#/ {print; next}
       NF >= 7 {
         if ($7 == "+") $7 = "-";
         else if ($7 == "-") $7 = "+";
         print;
         next
       }
       {print}' "${GTF_TE}" > "${ANTISENSE_TE_GTF}.tmp"

  mv "${ANTISENSE_TE_GTF}.tmp" "${ANTISENSE_TE_GTF}"

  {
    echo "# TE antisense strand-flipped annotation"
    echo
    echo "Source sense TE GTF:"
    echo "${GTF_TE}"
    echo
    echo "Generated file:"
    echo "${ANTISENSE_TE_GTF}"
    echo
    echo "Rule:"
    echo "GTF column 7 strand was flipped: + -> -, - -> +. All other columns were kept unchanged."
    echo
    echo "Intended usage for reverse-stranded RNA-seq TEcount:"
    echo "sense: original m39_TE.gtf + --stranded reverse"
    echo "antisense: this strand-flipped GTF + --stranded reverse"
    echo
    echo "Generated on: $(date)"
  } > "${README_FILE}"

  echo "Original strand counts:"
  awk 'BEGIN{FS="\t"} !/^#/ {s[$7]++} END{for (k in s) print k, s[k]}' "${GTF_TE}" | sort
  echo "Flipped strand counts:"
  awk 'BEGIN{FS="\t"} !/^#/ {s[$7]++} END{for (k in s) print k, s[k]}' "${ANTISENSE_TE_GTF}" | sort
  echo "Flipped GTF size:"
  ls -lh "${ANTISENSE_TE_GTF}"
  echo "[$(date)] Done"
} 2>&1 | tee "${LOG_FILE}"
