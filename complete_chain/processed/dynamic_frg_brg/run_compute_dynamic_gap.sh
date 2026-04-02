#!/bin/bash

# Usage: bash run_compute_dynamic_gap.sh sampleType input.root output.root N_events offset

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <sampleType> <input.root> <output.root> <N_events> <offset>"
  exit 1
fi

sampleType=$1
inputFile=$2
outputFile=$3
maxEvents=$4
offset=$5

# Environment
DELPHES_PATH="/LUSTRE/home/acota/delphes/Delphes_default/Delphes-3.5.0"
export CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate latest_root

export LD_LIBRARY_PATH="$DELPHES_PATH:$LD_LIBRARY_PATH"
export ROOT_INCLUDE_PATH="$DELPHES_PATH:$DELPHES_PATH/external:$ROOT_INCLUDE_PATH"

mkdir -p "$(dirname "$outputFile")"

# Run ROOT macro
echo "Processing chunk: offset=$offset, maxEvents=$maxEvents"
root -l -b -q "compute_dynamic_gap.C(\"$inputFile\",\"$outputFile\",$sampleType,$maxEvents,$offset)"

