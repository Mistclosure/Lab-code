#!/usr/bin/env bash
set -euo pipefail

# Run from anywhere. Long-running steps are listed explicitly.

PROJECT_CODE="/home/zerui/code/Phf20_1.23"

echo "# 1. Build or refresh the shared strand-flipped TE annotation"
echo "bash \"${PROJECT_CODE}/scripts/01_make_antisense_te_gtf.sh\""

echo
echo "# 2. Run fastp + bowtie2 rRNA/mtDNA removal + STAR + sense TEcount"
echo "bash \"${PROJECT_CODE}/run_step1_4_TEcount_with_rRNA_removal.sh\""

echo
echo "# 3. Run sense/antisense TEcount driver"
echo "# Sense counts from counts/ are linked/reused; antisense uses the strand-flipped TE GTF."
echo "DRY_RUN=0 MAX_JOBS=4 bash \"${PROJECT_CODE}/scripts/02_run_TEcount_sense_antisense.sh\""

echo
echo "# 4. Merge TEcount matrices and infer metadata if needed"
echo "Rscript \"${PROJECT_CODE}/scripts/03_merge_TEcount_matrix.R\""

echo
echo "# 5. edgeR differential analysis and volcano plots for LINE/SINE/LTR"
echo "Rscript \"${PROJECT_CODE}/scripts/04_edgeR_sense_antisense_LINE_SINE_LTR_volcano.R\""

echo
echo "# 6. Plot D heatmaps"
echo "Rscript \"${PROJECT_CODE}/scripts/05_plot_D_heatmap.R\""

echo
echo "# 7. Plot E summary"
echo "Rscript \"${PROJECT_CODE}/scripts/06_plot_E_summary.R\""
