#!/bin/bash

# Usage:
# bash run_compute_general_sample.sh sampleType input.root output.root N_events offset
# sampleType: 0 = photon, 1 = pomeron, 2 = background

if [[ $# -ne 5 ]]; then
  echo "Usage: $0 <sampleType> <input.root> <output.root> <N_events> <offset>"
  echo "sampleType: 0 = photon, 1 = pomeron, 2 = background"
  exit 1
fi

sampleType=$1
inputFile=$2
outputFile=$3
maxEvents=$4
offset=$5

DELPHES_PATH="/lustre/home/amendieta/software/Delphes-3.5.0"

export CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate latest_root

export LD_LIBRARY_PATH="$DELPHES_PATH:$LD_LIBRARY_PATH"
export ROOT_INCLUDE_PATH="$DELPHES_PATH:$DELPHES_PATH/external:$ROOT_INCLUDE_PATH"

mkdir -p "$(dirname "$outputFile")"

echo "INPUT_FILE   = $inputFile"
echo "OUTPUT_FILE  = $outputFile"
echo "SAMPLE_TYPE  = $sampleType"
echo "MAX_EVENTS   = $maxEvents"
echo "OFFSET       = $offset"

root -l -b -q "compute_general_sample.C(\"$inputFile\",\"$outputFile\",$sampleType,$maxEvents,$offset)"
