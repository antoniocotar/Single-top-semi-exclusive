#!/bin/bash

# Edit these once per sample submission.
#0=photon, 1=pomeron, 2=qcd
SAMPLE_TYPE=2
#Pomeron:
#INPUT_FILE="/lustre/proyectos/p018/single-top_semi-exclusive/samples/delphes/top_pomeron_pomset6_pomflux8.root"
#Photon:
#INPUT_FILE="/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/twphoto.root"
#Background:
INPUT_FILE="/lustre/proyectos/p018/ttbar_semi-exclusive_pheno/ocotillo/delphes_samples/delphes_ttbar_central_NLO_pp13TeV.root"
#Pomeron:
#OUTPUT_DIR="/lustre/home/amendieta/guille/Pp_Alessio/singletop/samples_processed/pomeron/top_pomeron_pomset6_pomflux8"
#Photon:
#OUTPUT_DIR="/lustre/home/amendieta/guille/Pp_Alessio/singletop/samples_processed/photon"
#Background:
OUTPUT_DIR="/lustre/home/amendieta/guille/Pp_Alessio/singletop/samples_processed/background"
#Pomeron:
#OUTPUT_BASENAME="top_pomeron_pomset6_pomflux8_processed"
#Photon:
#OUTPUT_BASENAME="twphoto_processed"
#Background:
OUTPUT_BASENAME="delphes_ttbar_central_NLO_pp13TeV_processed"
FINAL_OUTPUT_DIR="/lustre/proyectos/p018/single-top_semi-exclusive/samples/processed"
FINAL_OUTPUT_NAME="${OUTPUT_BASENAME}.root"
CHUNK_SIZE=10000

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SAMPLE_SLURM="${SCRIPT_DIR}/run_general_sample.slrm"
RUN_HADD_SLURM="${SCRIPT_DIR}/run_general_hadd.slrm"

if [ ! -f "$RUN_SAMPLE_SLURM" ]; then
  echo "ERROR: $RUN_SAMPLE_SLURM not found"
  exit 1
fi

if [ ! -f "$RUN_HADD_SLURM" ]; then
  echo "ERROR: $RUN_HADD_SLURM not found"
  exit 1
fi

module load conda
source activate latest_root

TOTAL_ENTRIES=$(root -l -b -q -e "
TFile f(\"${INPUT_FILE}\");
TTree *t = (TTree*)f.Get(\"Delphes\");
if (!t) { std::cerr << \"NO_TREE\" << std::endl; gSystem->Exit(1); }
std::cout << t->GetEntries() << std::endl;
gSystem->Exit(0);
" 2>/dev/null | tail -n 1)


if ! [[ "$TOTAL_ENTRIES" =~ ^[0-9]+$ ]]; then
  echo "ERROR: could not determine total entries for $INPUT_FILE"
  exit 1
fi

TOTAL_CHUNKS=$(( (TOTAL_ENTRIES + CHUNK_SIZE - 1) / CHUNK_SIZE ))
if [ "$TOTAL_CHUNKS" -le 0 ]; then
  echo "ERROR: total chunks computed as $TOTAL_CHUNKS"
  exit 1
fi

LAST_TASK=$((TOTAL_CHUNKS - 1))

echo "TOTAL_ENTRIES = $TOTAL_ENTRIES"
echo "CHUNK_SIZE    = $CHUNK_SIZE"
echo "TOTAL_CHUNKS  = $TOTAL_CHUNKS"

EXPORT_VARS="ALL,SAMPLE_TYPE=${SAMPLE_TYPE},INPUT_FILE=${INPUT_FILE},OUTPUT_DIR=${OUTPUT_DIR},OUTPUT_BASENAME=${OUTPUT_BASENAME},CHUNK_SIZE=${CHUNK_SIZE}"
ARRAY_SUBMISSION=$(sbatch --array="0-${LAST_TASK}" --export="$EXPORT_VARS" "$RUN_SAMPLE_SLURM")
ARRAY_JOB_ID=$(echo "$ARRAY_SUBMISSION" | awk '{print $4}')

if [ -z "$ARRAY_JOB_ID" ]; then
  echo "ERROR: could not parse array job id from: $ARRAY_SUBMISSION"
  exit 1
fi

echo "Submitted chunk array job: $ARRAY_JOB_ID"

HADD_EXPORT_VARS="ALL,OUTPUT_DIR=${OUTPUT_DIR},OUTPUT_BASENAME=${OUTPUT_BASENAME},FINAL_OUTPUT_DIR=${FINAL_OUTPUT_DIR},FINAL_OUTPUT_NAME=${FINAL_OUTPUT_NAME}"
HADD_SUBMISSION=$(sbatch --dependency=afterok:${ARRAY_JOB_ID} --export="$HADD_EXPORT_VARS" "$RUN_HADD_SLURM")
HADD_JOB_ID=$(echo "$HADD_SUBMISSION" | awk '{print $4}')

if [ -z "$HADD_JOB_ID" ]; then
  echo "ERROR: could not parse hadd job id from: $HADD_SUBMISSION"
  exit 1
fi

echo "Submitted hadd job: $HADD_JOB_ID"
echo "Pipeline ready:"
echo "  array job = $ARRAY_JOB_ID"
echo "  hadd job  = $HADD_JOB_ID"
