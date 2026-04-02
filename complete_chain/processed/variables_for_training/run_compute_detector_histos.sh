#!/bin/bash
# run_compute_detector_histos.sh

# ─── User‐configurable defaults ────────────────────────────────────────────────
DELPHES_PATH="/LUSTRE/home/acota/delphes/Delphes_default/Delphes-3.5.0"
DEFAULT_IN="/LUSTRE/home/acota/.../ap_pa_ttbar_NLO_100k.root"
DEFAULT_OUT="/LUSTRE/home/acota/.../output_processed_root_files/processed_ap_pa_ttbar_NLO_100k.root"
DEFAULT_MAX=10000
DEFAULT_OFF=0

# ─── Override with positional arguments ───────────────────────────────────────
INPUT_FILE=${1:-$DEFAULT_IN}
OUTPUT_FILE=${2:-$DEFAULT_OUT}
MAX_EVENTS=${3:-$DEFAULT_MAX}
OFFSET=${4:-$DEFAULT_OFF}

# ─── Create output directory if needed ────────────────────────────────────────
mkdir -p "$(dirname "$OUTPUT_FILE")"

# support “all” shorthand
if [[ "$MAX_EVENTS" == "all" ]]; then
  MAX_EVENTS=-1
fi

echo "INPUT_FILE  = $INPUT_FILE"
echo "OUTPUT_FILE = $OUTPUT_FILE"
echo "MAX_EVENTS  = $MAX_EVENTS"
echo "OFFSET      = $OFFSET"

# ─── 1) Activate ROOT/Delphes environment ────────────────────────────────────
export CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate latest_root

# ─── 2) Make Delphes visible ─────────────────────────────────────────────────
export LD_LIBRARY_PATH="$DELPHES_PATH:$LD_LIBRARY_PATH"
export ROOT_INCLUDE_PATH="$DELPHES_PATH:$DELPHES_PATH/external:$ROOT_INCLUDE_PATH"

# ─── 3) Run the macro ────────────────────────────────────────────────────────
root -l -b -q "compute_detector_histos.C(\"$INPUT_FILE\",\"$OUTPUT_FILE\",$MAX_EVENTS,$OFFSET)"

