#!/bin/bash
# run_compute_bjorken.sh

# ─── User‐configurable defaults ────────────────────────────────────────────────
DELPHES_PATH="/lustre/home/acota/ttbar_semi_exclusive_pheno_paper/delphes/Delphes-3.5.0"
DEFAULT_IN="/LUSTRE/home/acota/.../ap_pa_ttbar_NLO_100k.root"
DEFAULT_OUT="/LUSTRE/home/acota/.../output_processed_root_files/processed_ap_pa_ttbar_NLO_100k.root"
DEFAULT_MAX=10000

# ─── Override with positional arguments ───────────────────────────────────────
INPUT_FILE=${1:-$DEFAULT_IN}
OUTPUT_FILE=${2:-$DEFAULT_OUT}
MAX_EVENTS=${3:-$DEFAULT_MAX}

# ─── Create output directory if needed ────────────────────────────────────────
mkdir -p "$(dirname "$OUTPUT_FILE")"

if [[ "$MAX_EVENTS" == "all" ]]; then
  MAX_EVENTS=-1
fi

echo "INPUT_FILE  = $INPUT_FILE"
echo "OUTPUT_FILE = $OUTPUT_FILE"
echo "MAX_EVENTS  = $MAX_EVENTS"

# ─── 1) Activate ROOT/Delphes environment ────────────────────────────────────
export CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate latest_root

# ─── 2) Make Delphes visible ─────────────────────────────────────────────────
export LD_LIBRARY_PATH="$DELPHES_PATH:$LD_LIBRARY_PATH"
export ROOT_INCLUDE_PATH="$DELPHES_PATH/classes:$DELPHES_PATH/external:$ROOT_INCLUDE_PATH"

# ─── 3) Run the macro ────────────────────────────────────────────────────────
root -l -b -q "compute_bjorken_pomeron.C(\"$INPUT_FILE\",\"$OUTPUT_FILE\",$MAX_EVENTS)"

