#!/usr/bin/env bash
set -euo pipefail

CODE_DIR="/home/zerui/code/Phf20 CT26 ribosomal depletion RNA-seq"

echo "# 1. Make antisense TE GTF. This is short and only flips column 7 of the TE GTF."
echo "bash \"${CODE_DIR}/scripts/01_make_antisense_te_gtf.sh\""
echo
echo "# 2. Dry-run TEcount commands. This does not start long TEcount jobs."
echo "DRY_RUN=1 bash \"${CODE_DIR}/scripts/02_run_TEcount_sense_antisense.sh\""
echo
echo "# 3. Actual TEcount run, after you confirm the dry-run output."
echo "DRY_RUN=0 bash \"${CODE_DIR}/scripts/02_run_TEcount_sense_antisense.sh\""
echo
echo "# 4. Merge matrices after TEcount is complete."
echo "Rscript \"${CODE_DIR}/scripts/03_merge_TEcount_matrix.R\""
echo
echo "# 5. edgeR analysis and volcano plots for LINE/SINE/LTR."
echo "Rscript \"${CODE_DIR}/scripts/04_edgeR_sense_antisense_LINE_SINE_LTR_volcano.R\""
echo
echo "# 6. Plot D-style heatmaps for LINE/SINE/LTR."
echo "Rscript \"${CODE_DIR}/scripts/05_plot_D_heatmap.R\""
echo
echo "# 7. Plot E/F-style bidirectional LINE/SINE/LTR fold-change summary."
echo "Rscript \"${CODE_DIR}/scripts/06_plot_E_summary.R\""
