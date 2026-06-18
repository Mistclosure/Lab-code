#!/usr/bin/env bash
set -eo pipefail

# By default this script is a dry run. Use DRY_RUN=0 to actually run TEcount.
DRY_RUN="${DRY_RUN:-1}"
MAX_JOBS="${MAX_JOBS:-4}"
FORCE_RERUN_SENSE="${FORCE_RERUN_SENSE:-0}"
FORCE_RERUN_ANTISENSE="${FORCE_RERUN_ANTISENSE:-0}"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  # shellcheck disable=SC1091
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi
conda activate te_env

set -u

ulimit -n 65535

WORKDIR="/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
ALIGN_DIR="${WORKDIR}/alignments"
EXISTING_SENSE_COUNTS_DIR="${WORKDIR}/counts"

ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"
ANTISENSE_TE_GTF="${ANNO_DIR}/TE_antisense_strand_flipped/m39_TE.antisense.strand_flipped.gtf"

ANALYSIS_DIR="${WORKDIR}/sense_antisense_ERV"
SENSE_COUNTS_DIR="${ANALYSIS_DIR}/counts/TEcount_sense"
ANTISENSE_COUNTS_DIR="${ANALYSIS_DIR}/counts/TEcount_antisense"
LOG_DIR="${ANALYSIS_DIR}/logs"
if [ "${DRY_RUN}" != "1" ]; then
  mkdir -p "${SENSE_COUNTS_DIR}" "${ANTISENSE_COUNTS_DIR}" "${LOG_DIR}"
fi

valid_bam() {
  local bam_file="$1"
  [ -s "$bam_file" ] && samtools quickcheck "$bam_file" >/dev/null 2>&1
}

run_or_print() {
  local log_file="$1"
  shift
  if [ "${DRY_RUN}" = "1" ]; then
    printf '[DRY_RUN] '
    printf '%q ' "$@"
    printf '\n'
  else
    "$@" > "${log_file}" 2>&1
  fi
}

echo "=== TEcount sense/antisense ERV workflow ==="
echo "DRY_RUN=${DRY_RUN}"
echo "WORKDIR=${WORKDIR}"
echo "GTF_GENE=${GTF_GENE}"
echo "GTF_TE sense=${GTF_TE}"
echo "GTF_TE antisense=${ANTISENSE_TE_GTF}"
echo "Important: both sense and antisense use --stranded reverse"

if ! command -v TEcount >/dev/null 2>&1; then
  echo "ERROR: TEcount not found in te_env." >&2
  exit 1
fi
if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not found in te_env." >&2
  exit 1
fi
if [ ! -s "${GTF_GENE}" ] || [ ! -s "${GTF_TE}" ]; then
  echo "ERROR: GTF_GENE or GTF_TE is missing." >&2
  exit 1
fi
if [ ! -s "${ANTISENSE_TE_GTF}" ] && [ "${DRY_RUN}" = "1" ]; then
  echo "WARN: Antisense TE GTF is not present yet; dry-run will still print commands."
  echo "      Run scripts/01_make_antisense_te_gtf.sh before the real TEcount run."
elif [ ! -s "${ANTISENSE_TE_GTF}" ]; then
  echo "ERROR: Antisense TE GTF is missing. Run scripts/01_make_antisense_te_gtf.sh first." >&2
  exit 1
fi

shopt -s nullglob
bam_files=("${ALIGN_DIR}"/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

if [ ${#bam_files[@]} -eq 0 ]; then
  echo "ERROR: No BAM files found in ${ALIGN_DIR}" >&2
  exit 1
fi

for bam_file in "${bam_files[@]}"; do
  sample_name=$(basename "${bam_file}" .Aligned.sortedByCoord.out.bam)
  if ! valid_bam "${bam_file}"; then
    echo "WARN: skip invalid BAM: ${bam_file}" >&2
    continue
  fi

  sense_out="${SENSE_COUNTS_DIR}/${sample_name}.cntTable"
  existing_sense="${EXISTING_SENSE_COUNTS_DIR}/${sample_name}.cntTable"
  if [ "${FORCE_RERUN_SENSE}" != "1" ] && [ -s "${existing_sense}" ]; then
    if [ "${DRY_RUN}" = "1" ]; then
      echo "[DRY_RUN] ln -sf '${existing_sense}' '${sense_out}'"
    else
      ln -sf "${existing_sense}" "${sense_out}"
    fi
    echo "Reuse completed sense TEcount: ${sample_name}"
  elif [ -s "${sense_out}" ] && [ "${FORCE_RERUN_SENSE}" != "1" ]; then
    echo "Skip existing sense TEcount: ${sample_name}"
  else
    (
      run_or_print "${LOG_DIR}/${sample_name}.sense.TEcount.log" \
        TEcount --sortByPos --format BAM --mode multi \
        --GTF "${GTF_GENE}" --TE "${GTF_TE}" \
        --project "${SENSE_COUNTS_DIR}/${sample_name}" \
        --stranded reverse -b "${bam_file}"
      echo "Done sense TEcount: ${sample_name}"
    ) &
  fi

  antisense_out="${ANTISENSE_COUNTS_DIR}/${sample_name}.cntTable"
  if [ -s "${antisense_out}" ] && [ "${FORCE_RERUN_ANTISENSE}" != "1" ]; then
    echo "Skip existing antisense TEcount: ${sample_name}"
  else
    (
      run_or_print "${LOG_DIR}/${sample_name}.antisense.TEcount.log" \
        TEcount --sortByPos --format BAM --mode multi \
        --GTF "${GTF_GENE}" --TE "${ANTISENSE_TE_GTF}" \
        --project "${ANTISENSE_COUNTS_DIR}/${sample_name}" \
        --stranded reverse -b "${bam_file}"
      echo "Done antisense TEcount: ${sample_name}"
    ) &
  fi

  while (( $(jobs -r -p | wc -l) >= MAX_JOBS )); do
    sleep 10
  done
done

wait
echo "=== Finished TEcount sense/antisense driver ==="
