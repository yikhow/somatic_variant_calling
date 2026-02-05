
#!/bin/bash
# ==============================================================================
# Script: run_pipeline.sh
# Description: Master orchestrator.
# ==============================================================================

set -euo pipefail

# Resolve paths
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )/scr"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

echo "----------------------------------------------------------------"
echo "[MASTER] Starting Pipeline in: $PROJECT_DIR"
echo "----------------------------------------------------------------"

# Make executable
chmod +x "${SCRIPT_DIR}"/*.sh

# 1. Preprocessing
"${SCRIPT_DIR}/preprocessing.sh"

# 2. Variant Calling
"${SCRIPT_DIR}/variant_calling.sh"

# 3. Post Processing
"${SCRIPT_DIR}/post_processing.sh"

# 4. R Visualization
# Pass project directory to R script
if [ -f "${SCRIPT_DIR}/somatic_visualization.r" ]; then
    echo "[MASTER] Running R Visualization..."
    Rscript --vanilla "${SCRIPT_DIR}/somatic_visualization.r" "$PROJECT_DIR"
else
    echo "[MASTER] R script missing."
fi

echo "[MASTER] Pipeline Complete."
