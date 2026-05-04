// singletop_compute_detector_histos.C
#ifdef __CLING__
  R__LOAD_LIBRARY(libDelphes)
  #include "classes/DelphesClasses.h"
  #include "ExRootAnalysis/ExRootTreeReader.h"
#else
  class ExRootTreeReader;
  class GenParticle;
  class Jet;
  class Tower;
  class Track;
  class Electron;
  class Muon;
  class MissingET;
#endif

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <Riostream.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

void singletop_compute_detector_histos(const char* inputFile,
                             const char* outputFile,
                             Long64_t    maxEvents = -1,
                             Long64_t    offset    = 0,
                             float       Emin      = 100.0f,
                             float       Emax      = 300.0f)
{
  // Load Delphes
  gSystem->Load("libDelphes");

  // Open input chain
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create a tree reader
  ExRootTreeReader* reader = new ExRootTreeReader(&chain);
  Long64_t nEntries = reader->GetEntries();

  // Validate offset
  if (offset < 0) offset = 0;
  if (offset >= nEntries) {
    std::cerr << "Offset ("<<offset<<") >= entries ("<<nEntries<<"). Aborting.\n";
    return;
  }

  // Compute number of events to loop
  Long64_t nLoop = (maxEvents>0 && maxEvents < (nEntries-offset))
                   ? maxEvents
                   : (nEntries - offset);
  std::cout << "Processing "<<nLoop<<" events starting at "<<offset<<"\n";

  // Set up branches
  TClonesArray* branchJet       = reader->UseBranch("Jet");
  TClonesArray* branchTrack     = reader->UseBranch("Track");
  TClonesArray* branchTower     = reader->UseBranch("Tower");
  TClonesArray* branchElectron  = reader->UseBranch("Electron");
  TClonesArray* branchMuon      = reader->UseBranch("Muon");
  TClonesArray* branchMissingET = reader->UseBranch("MissingET");

  // Prepare output file & tree
  TFile* outFile = TFile::Open(outputFile, "RECREATE");
  TTree* tree = new TTree("Events","Reco observables per event");

  // Declare variables and branch them
  float rapidityGapPlus  = -1;
  float rapidityGapMinus = -1;
  float hfSumPlus        =  0;
  float hfSumMinus       =  0;
  float nchPV            =  0;
  float nForwardJets     =  0;
  float nTotalJets       =  0;
  float leadLightJetPt   =  0;
  float maxHFEnergy      =  std::numeric_limits<float>::lowest();
  float minHFEnergy      =  std::numeric_limits<float>::max();
  float METval           =  0;
  float maxHFEnergyPlus   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus  = std::numeric_limits<float>::lowest();
  float minOfMaxHF        = -1.0f;
  float minOfHFSums       = -1.0f;

  // newnew floats
  float maxHFEnergy_30_31       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_31       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_31   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_31  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_31        = -1.0f;

  float maxHFEnergy_30_32       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_32       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_32   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_32  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_32        = -1.0f;

  float maxHFEnergy_30_33       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_33       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_33   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_33  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_33        = -1.0f;

  float maxHFEnergy_30_34       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_34       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_34   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_34  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_34        = -1.0f;

  float maxHFEnergy_30_35       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_35       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_35   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_35  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_35        = -1.0f;

  float maxHFEnergy_30_36       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_36       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_36   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_36  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_36        = -1.0f;

  float maxHFEnergy_30_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_37        = -1.0f;

  float maxHFEnergy_30_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_38        = -1.0f;

  float maxHFEnergy_30_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_39        = -1.0f;

  float maxHFEnergy_30_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_40        = -1.0f;

  float maxHFEnergy_30_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_41        = -1.0f;

  float maxHFEnergy_30_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_42        = -1.0f;

  float maxHFEnergy_30_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_43        = -1.0f;

  float maxHFEnergy_30_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_44        = -1.0f;

  float maxHFEnergy_30_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_45        = -1.0f;

  float maxHFEnergy_30_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_46        = -1.0f;

  float maxHFEnergy_30_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_47        = -1.0f;

  float maxHFEnergy_30_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_48        = -1.0f;

  float maxHFEnergy_30_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_49        = -1.0f;

  float maxHFEnergy_30_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_50        = -1.0f;

  float maxHFEnergy_30_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_51        = -1.0f;

  float maxHFEnergy_30_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_30_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_30_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_30_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_30_52        = -1.0f;



  float maxHFEnergy_31_32       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_32       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_32   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_32  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_32        = -1.0f;

  float maxHFEnergy_31_33       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_33       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_33   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_33  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_33        = -1.0f;

  float maxHFEnergy_31_34       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_34       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_34   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_34  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_34        = -1.0f;

  float maxHFEnergy_31_35       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_35       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_35   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_35  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_35        = -1.0f;

  float maxHFEnergy_31_36       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_36       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_36   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_36  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_36        = -1.0f;

  float maxHFEnergy_31_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_37        = -1.0f;

  float maxHFEnergy_31_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_38        = -1.0f;

  float maxHFEnergy_31_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_39        = -1.0f;

  float maxHFEnergy_31_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_40        = -1.0f;

  float maxHFEnergy_31_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_41        = -1.0f;

  float maxHFEnergy_31_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_42        = -1.0f;

  float maxHFEnergy_31_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_43        = -1.0f;

  float maxHFEnergy_31_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_44        = -1.0f;

  float maxHFEnergy_31_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_45        = -1.0f;

  float maxHFEnergy_31_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_46        = -1.0f;

  float maxHFEnergy_31_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_47        = -1.0f;

  float maxHFEnergy_31_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_48        = -1.0f;

  float maxHFEnergy_31_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_49        = -1.0f;

  float maxHFEnergy_31_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_50        = -1.0f;

  float maxHFEnergy_31_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_51        = -1.0f;

  float maxHFEnergy_31_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_31_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_31_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_31_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_31_52        = -1.0f;



  float maxHFEnergy_32_33       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_33       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_33   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_33  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_33        = -1.0f;

  float maxHFEnergy_32_34       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_34       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_34   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_34  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_34        = -1.0f;

  float maxHFEnergy_32_35       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_35       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_35   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_35  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_35        = -1.0f;

  float maxHFEnergy_32_36       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_36       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_36   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_36  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_36        = -1.0f;

  float maxHFEnergy_32_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_37        = -1.0f;

  float maxHFEnergy_32_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_38        = -1.0f;

  float maxHFEnergy_32_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_39        = -1.0f;

  float maxHFEnergy_32_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_40        = -1.0f;

  float maxHFEnergy_32_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_41        = -1.0f;

  float maxHFEnergy_32_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_42        = -1.0f;

  float maxHFEnergy_32_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_43        = -1.0f;

  float maxHFEnergy_32_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_44        = -1.0f;

  float maxHFEnergy_32_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_45        = -1.0f;

  float maxHFEnergy_32_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_46        = -1.0f;

  float maxHFEnergy_32_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_47        = -1.0f;

  float maxHFEnergy_32_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_48        = -1.0f;

  float maxHFEnergy_32_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_49        = -1.0f;

  float maxHFEnergy_32_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_50        = -1.0f;

  float maxHFEnergy_32_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_51        = -1.0f;

  float maxHFEnergy_32_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_32_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_32_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_32_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_32_52        = -1.0f;



  float maxHFEnergy_33_34       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_34       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_34   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_34  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_34        = -1.0f;

  float maxHFEnergy_33_35       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_35       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_35   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_35  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_35        = -1.0f;

  float maxHFEnergy_33_36       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_36       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_36   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_36  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_36        = -1.0f;

  float maxHFEnergy_33_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_37        = -1.0f;

  float maxHFEnergy_33_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_38        = -1.0f;

  float maxHFEnergy_33_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_39        = -1.0f;

  float maxHFEnergy_33_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_40        = -1.0f;

  float maxHFEnergy_33_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_41        = -1.0f;

  float maxHFEnergy_33_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_42        = -1.0f;

  float maxHFEnergy_33_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_43        = -1.0f;

  float maxHFEnergy_33_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_44        = -1.0f;

  float maxHFEnergy_33_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_45        = -1.0f;

  float maxHFEnergy_33_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_46        = -1.0f;

  float maxHFEnergy_33_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_47        = -1.0f;

  float maxHFEnergy_33_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_48        = -1.0f;

  float maxHFEnergy_33_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_49        = -1.0f;

  float maxHFEnergy_33_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_50        = -1.0f;

  float maxHFEnergy_33_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_51        = -1.0f;

  float maxHFEnergy_33_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_33_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_33_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_33_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_33_52        = -1.0f;



  float maxHFEnergy_34_35       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_35       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_35   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_35  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_35        = -1.0f;

  float maxHFEnergy_34_36       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_36       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_36   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_36  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_36        = -1.0f;

  float maxHFEnergy_34_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_37        = -1.0f;

  float maxHFEnergy_34_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_38        = -1.0f;

  float maxHFEnergy_34_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_39        = -1.0f;

  float maxHFEnergy_34_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_40        = -1.0f;

  float maxHFEnergy_34_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_41        = -1.0f;

  float maxHFEnergy_34_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_42        = -1.0f;

  float maxHFEnergy_34_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_43        = -1.0f;

  float maxHFEnergy_34_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_44        = -1.0f;

  float maxHFEnergy_34_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_45        = -1.0f;

  float maxHFEnergy_34_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_46        = -1.0f;

  float maxHFEnergy_34_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_47        = -1.0f;

  float maxHFEnergy_34_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_48        = -1.0f;

  float maxHFEnergy_34_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_49        = -1.0f;

  float maxHFEnergy_34_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_50        = -1.0f;

  float maxHFEnergy_34_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_51        = -1.0f;

  float maxHFEnergy_34_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_34_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_34_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_34_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_34_52        = -1.0f;



  float maxHFEnergy_35_36       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_36       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_36   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_36  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_36        = -1.0f;

  float maxHFEnergy_35_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_37        = -1.0f;

  float maxHFEnergy_35_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_38        = -1.0f;

  float maxHFEnergy_35_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_39        = -1.0f;

  float maxHFEnergy_35_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_40        = -1.0f;

  float maxHFEnergy_35_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_41        = -1.0f;

  float maxHFEnergy_35_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_42        = -1.0f;

  float maxHFEnergy_35_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_43        = -1.0f;

  float maxHFEnergy_35_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_44        = -1.0f;

  float maxHFEnergy_35_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_45        = -1.0f;

  float maxHFEnergy_35_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_46        = -1.0f;

  float maxHFEnergy_35_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_47        = -1.0f;

  float maxHFEnergy_35_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_48        = -1.0f;

  float maxHFEnergy_35_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_49        = -1.0f;

  float maxHFEnergy_35_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_50        = -1.0f;

  float maxHFEnergy_35_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_51        = -1.0f;

  float maxHFEnergy_35_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_35_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_35_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_35_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_35_52        = -1.0f;



  float maxHFEnergy_36_37       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_37       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_37   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_37  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_37        = -1.0f;

  float maxHFEnergy_36_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_38        = -1.0f;

  float maxHFEnergy_36_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_39        = -1.0f;

  float maxHFEnergy_36_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_40        = -1.0f;

  float maxHFEnergy_36_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_41        = -1.0f;

  float maxHFEnergy_36_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_42        = -1.0f;

  float maxHFEnergy_36_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_43        = -1.0f;

  float maxHFEnergy_36_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_44        = -1.0f;

  float maxHFEnergy_36_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_45        = -1.0f;

  float maxHFEnergy_36_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_46        = -1.0f;

  float maxHFEnergy_36_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_47        = -1.0f;

  float maxHFEnergy_36_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_48        = -1.0f;

  float maxHFEnergy_36_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_49        = -1.0f;

  float maxHFEnergy_36_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_50        = -1.0f;

  float maxHFEnergy_36_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_51        = -1.0f;

  float maxHFEnergy_36_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_36_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_36_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_36_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_36_52        = -1.0f;



  float maxHFEnergy_37_38       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_38       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_38   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_38  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_38        = -1.0f;

  float maxHFEnergy_37_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_39        = -1.0f;

  float maxHFEnergy_37_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_40        = -1.0f;

  float maxHFEnergy_37_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_41        = -1.0f;

  float maxHFEnergy_37_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_42        = -1.0f;

  float maxHFEnergy_37_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_43        = -1.0f;

  float maxHFEnergy_37_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_44        = -1.0f;

  float maxHFEnergy_37_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_45        = -1.0f;

  float maxHFEnergy_37_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_46        = -1.0f;

  float maxHFEnergy_37_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_47        = -1.0f;

  float maxHFEnergy_37_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_48        = -1.0f;

  float maxHFEnergy_37_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_49        = -1.0f;

  float maxHFEnergy_37_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_50        = -1.0f;

  float maxHFEnergy_37_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_51        = -1.0f;

  float maxHFEnergy_37_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_37_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_37_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_37_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_37_52        = -1.0f;



  float maxHFEnergy_38_39       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_39       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_39   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_39  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_39        = -1.0f;

  float maxHFEnergy_38_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_40        = -1.0f;

  float maxHFEnergy_38_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_41        = -1.0f;

  float maxHFEnergy_38_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_42        = -1.0f;

  float maxHFEnergy_38_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_43        = -1.0f;

  float maxHFEnergy_38_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_44        = -1.0f;

  float maxHFEnergy_38_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_45        = -1.0f;

  float maxHFEnergy_38_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_46        = -1.0f;

  float maxHFEnergy_38_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_47        = -1.0f;

  float maxHFEnergy_38_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_48        = -1.0f;

  float maxHFEnergy_38_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_49        = -1.0f;

  float maxHFEnergy_38_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_50        = -1.0f;

  float maxHFEnergy_38_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_51        = -1.0f;

  float maxHFEnergy_38_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_38_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_38_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_38_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_38_52        = -1.0f;



  float maxHFEnergy_39_40       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_40       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_40   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_40  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_40        = -1.0f;

  float maxHFEnergy_39_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_41        = -1.0f;

  float maxHFEnergy_39_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_42        = -1.0f;

  float maxHFEnergy_39_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_43        = -1.0f;

  float maxHFEnergy_39_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_44        = -1.0f;

  float maxHFEnergy_39_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_45        = -1.0f;

  float maxHFEnergy_39_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_46        = -1.0f;

  float maxHFEnergy_39_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_47        = -1.0f;

  float maxHFEnergy_39_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_48        = -1.0f;

  float maxHFEnergy_39_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_49        = -1.0f;

  float maxHFEnergy_39_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_50        = -1.0f;

  float maxHFEnergy_39_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_51        = -1.0f;

  float maxHFEnergy_39_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_39_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_39_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_39_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_39_52        = -1.0f;



  float maxHFEnergy_40_41       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_41       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_41   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_41  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_41        = -1.0f;

  float maxHFEnergy_40_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_42        = -1.0f;

  float maxHFEnergy_40_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_43        = -1.0f;

  float maxHFEnergy_40_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_44        = -1.0f;

  float maxHFEnergy_40_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_45        = -1.0f;

  float maxHFEnergy_40_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_46        = -1.0f;

  float maxHFEnergy_40_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_47        = -1.0f;

  float maxHFEnergy_40_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_48        = -1.0f;

  float maxHFEnergy_40_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_49        = -1.0f;

  float maxHFEnergy_40_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_50        = -1.0f;

  float maxHFEnergy_40_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_51        = -1.0f;

  float maxHFEnergy_40_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_40_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_40_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_40_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_40_52        = -1.0f;



  float maxHFEnergy_41_42       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_42       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_42   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_42  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_42        = -1.0f;

  float maxHFEnergy_41_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_43        = -1.0f;

  float maxHFEnergy_41_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_44        = -1.0f;

  float maxHFEnergy_41_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_45        = -1.0f;

  float maxHFEnergy_41_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_46        = -1.0f;

  float maxHFEnergy_41_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_47        = -1.0f;

  float maxHFEnergy_41_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_48        = -1.0f;

  float maxHFEnergy_41_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_49        = -1.0f;

  float maxHFEnergy_41_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_50        = -1.0f;

  float maxHFEnergy_41_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_51        = -1.0f;

  float maxHFEnergy_41_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_41_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_41_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_41_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_41_52        = -1.0f;



  float maxHFEnergy_42_43       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_43       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_43   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_43  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_43        = -1.0f;

  float maxHFEnergy_42_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_44        = -1.0f;

  float maxHFEnergy_42_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_45        = -1.0f;

  float maxHFEnergy_42_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_46        = -1.0f;

  float maxHFEnergy_42_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_47        = -1.0f;

  float maxHFEnergy_42_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_48        = -1.0f;

  float maxHFEnergy_42_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_49        = -1.0f;

  float maxHFEnergy_42_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_50        = -1.0f;

  float maxHFEnergy_42_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_51        = -1.0f;

  float maxHFEnergy_42_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_42_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_42_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_42_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_42_52        = -1.0f;



  float maxHFEnergy_43_44       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_44       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_44   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_44  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_44        = -1.0f;

  float maxHFEnergy_43_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_45        = -1.0f;

  float maxHFEnergy_43_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_46        = -1.0f;

  float maxHFEnergy_43_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_47        = -1.0f;

  float maxHFEnergy_43_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_48        = -1.0f;

  float maxHFEnergy_43_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_49        = -1.0f;

  float maxHFEnergy_43_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_50        = -1.0f;

  float maxHFEnergy_43_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_51        = -1.0f;

  float maxHFEnergy_43_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_43_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_43_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_43_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_43_52        = -1.0f;



  float maxHFEnergy_44_45       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_45       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_45   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_45  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_45        = -1.0f;

  float maxHFEnergy_44_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_46        = -1.0f;

  float maxHFEnergy_44_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_47        = -1.0f;

  float maxHFEnergy_44_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_48        = -1.0f;

  float maxHFEnergy_44_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_49        = -1.0f;

  float maxHFEnergy_44_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_50        = -1.0f;

  float maxHFEnergy_44_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_51        = -1.0f;

  float maxHFEnergy_44_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_44_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_44_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_44_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_44_52        = -1.0f;



  float maxHFEnergy_45_46       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_46       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_46   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_46  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_46        = -1.0f;

  float maxHFEnergy_45_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_47        = -1.0f;

  float maxHFEnergy_45_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_48        = -1.0f;

  float maxHFEnergy_45_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_49        = -1.0f;

  float maxHFEnergy_45_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_50        = -1.0f;

  float maxHFEnergy_45_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_51        = -1.0f;

  float maxHFEnergy_45_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_45_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_45_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_45_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_45_52        = -1.0f;



  float maxHFEnergy_46_47       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_46_47       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_46_47   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_46_47  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_46_47        = -1.0f;

  float maxHFEnergy_46_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_46_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_46_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_46_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_46_48        = -1.0f;

  float maxHFEnergy_46_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_46_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_46_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_46_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_46_49        = -1.0f;

  float maxHFEnergy_46_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_46_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_46_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_46_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_46_50        = -1.0f;

  float maxHFEnergy_46_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_46_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_46_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_46_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_46_51        = -1.0f;

  float maxHFEnergy_46_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_46_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_46_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_46_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_46_52        = -1.0f;



  float maxHFEnergy_47_48       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_47_48       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_47_48   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_47_48  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_47_48        = -1.0f;

  float maxHFEnergy_47_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_47_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_47_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_47_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_47_49        = -1.0f;

  float maxHFEnergy_47_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_47_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_47_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_47_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_47_50        = -1.0f;

  float maxHFEnergy_47_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_47_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_47_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_47_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_47_51        = -1.0f;

  float maxHFEnergy_47_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_47_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_47_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_47_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_47_52        = -1.0f;



  float maxHFEnergy_48_49       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_48_49       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_48_49   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_48_49  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_48_49        = -1.0f;

  float maxHFEnergy_48_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_48_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_48_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_48_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_48_50        = -1.0f;

  float maxHFEnergy_48_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_48_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_48_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_48_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_48_51        = -1.0f;

  float maxHFEnergy_48_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_48_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_48_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_48_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_48_52        = -1.0f;



  float maxHFEnergy_49_50       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_49_50       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_49_50   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_49_50  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_49_50        = -1.0f;

  float maxHFEnergy_49_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_49_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_49_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_49_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_49_51        = -1.0f;

  float maxHFEnergy_49_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_49_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_49_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_49_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_49_52        = -1.0f;



  float maxHFEnergy_50_51       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_50_51       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_50_51   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_50_51  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_50_51        = -1.0f;

  float maxHFEnergy_50_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_50_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_50_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_50_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_50_52        = -1.0f;



  float maxHFEnergy_51_52       =  std::numeric_limits<float>::lowest();
  float minHFEnergy_51_52       =  std::numeric_limits<float>::max();
  float maxHFEnergyPlus_51_52   = std::numeric_limits<float>::lowest();
  float maxHFEnergyMinus_51_52  = std::numeric_limits<float>::lowest();
  float minOfMaxHF_51_52        = -1.0f;
  // newnewfloats

  // New floats
  float hfSumPlus_30_31   = 0, hfSumMinus_30_31  = 0;
  float hfSumPlus_30_32   = 0, hfSumMinus_30_32  = 0;
  float hfSumPlus_30_33   = 0, hfSumMinus_30_33  = 0;
  float hfSumPlus_30_34   = 0, hfSumMinus_30_34  = 0;
  float hfSumPlus_30_35   = 0, hfSumMinus_30_35  = 0;
  float hfSumPlus_30_36   = 0, hfSumMinus_30_36  = 0;
  float hfSumPlus_30_37   = 0, hfSumMinus_30_37  = 0;
  float hfSumPlus_30_38   = 0, hfSumMinus_30_38  = 0;
  float hfSumPlus_30_39   = 0, hfSumMinus_30_39  = 0;
  float hfSumPlus_30_40   = 0, hfSumMinus_30_40  = 0;
  float hfSumPlus_30_41   = 0, hfSumMinus_30_41  = 0;
  float hfSumPlus_30_42   = 0, hfSumMinus_30_42  = 0;
  float hfSumPlus_30_43   = 0, hfSumMinus_30_43  = 0;
  float hfSumPlus_30_44   = 0, hfSumMinus_30_44  = 0;
  float hfSumPlus_30_45   = 0, hfSumMinus_30_45  = 0;
  float hfSumPlus_30_46   = 0, hfSumMinus_30_46  = 0;
  float hfSumPlus_30_47   = 0, hfSumMinus_30_47  = 0;
  float hfSumPlus_30_48   = 0, hfSumMinus_30_48  = 0;
  float hfSumPlus_30_49   = 0, hfSumMinus_30_49  = 0;
  float hfSumPlus_30_50   = 0, hfSumMinus_30_50  = 0;
  float hfSumPlus_30_51   = 0, hfSumMinus_30_51  = 0;
  float hfSumPlus_30_52   = 0, hfSumMinus_30_52  = 0;

  float hfSumPlus_31_32   = 0, hfSumMinus_31_32  = 0;
  float hfSumPlus_31_33   = 0, hfSumMinus_31_33  = 0;
  float hfSumPlus_31_34   = 0, hfSumMinus_31_34  = 0;
  float hfSumPlus_31_35   = 0, hfSumMinus_31_35  = 0;
  float hfSumPlus_31_36   = 0, hfSumMinus_31_36  = 0;
  float hfSumPlus_31_37   = 0, hfSumMinus_31_37  = 0;
  float hfSumPlus_31_38   = 0, hfSumMinus_31_38  = 0;
  float hfSumPlus_31_39   = 0, hfSumMinus_31_39  = 0;
  float hfSumPlus_31_40   = 0, hfSumMinus_31_40  = 0;
  float hfSumPlus_31_41   = 0, hfSumMinus_31_41  = 0;
  float hfSumPlus_31_42   = 0, hfSumMinus_31_42  = 0;
  float hfSumPlus_31_43   = 0, hfSumMinus_31_43  = 0;
  float hfSumPlus_31_44   = 0, hfSumMinus_31_44  = 0;
  float hfSumPlus_31_45   = 0, hfSumMinus_31_45  = 0;
  float hfSumPlus_31_46   = 0, hfSumMinus_31_46  = 0;
  float hfSumPlus_31_47   = 0, hfSumMinus_31_47  = 0;
  float hfSumPlus_31_48   = 0, hfSumMinus_31_48  = 0;
  float hfSumPlus_31_49   = 0, hfSumMinus_31_49  = 0;
  float hfSumPlus_31_50   = 0, hfSumMinus_31_50  = 0;
  float hfSumPlus_31_51   = 0, hfSumMinus_31_51  = 0;
  float hfSumPlus_31_52   = 0, hfSumMinus_31_52  = 0;

  float hfSumPlus_32_33   = 0, hfSumMinus_32_33  = 0;
  float hfSumPlus_32_34   = 0, hfSumMinus_32_34  = 0;
  float hfSumPlus_32_35   = 0, hfSumMinus_32_35  = 0;
  float hfSumPlus_32_36   = 0, hfSumMinus_32_36  = 0;
  float hfSumPlus_32_37   = 0, hfSumMinus_32_37  = 0;
  float hfSumPlus_32_38   = 0, hfSumMinus_32_38  = 0;
  float hfSumPlus_32_39   = 0, hfSumMinus_32_39  = 0;
  float hfSumPlus_32_40   = 0, hfSumMinus_32_40  = 0;
  float hfSumPlus_32_41   = 0, hfSumMinus_32_41  = 0;
  float hfSumPlus_32_42   = 0, hfSumMinus_32_42  = 0;
  float hfSumPlus_32_43   = 0, hfSumMinus_32_43  = 0;
  float hfSumPlus_32_44   = 0, hfSumMinus_32_44  = 0;
  float hfSumPlus_32_45   = 0, hfSumMinus_32_45  = 0;
  float hfSumPlus_32_46   = 0, hfSumMinus_32_46  = 0;
  float hfSumPlus_32_47   = 0, hfSumMinus_32_47  = 0;
  float hfSumPlus_32_48   = 0, hfSumMinus_32_48  = 0;
  float hfSumPlus_32_49   = 0, hfSumMinus_32_49  = 0;
  float hfSumPlus_32_50   = 0, hfSumMinus_32_50  = 0;
  float hfSumPlus_32_51   = 0, hfSumMinus_32_51  = 0;
  float hfSumPlus_32_52   = 0, hfSumMinus_32_52  = 0;

  float hfSumPlus_33_34   = 0, hfSumMinus_33_34  = 0;
  float hfSumPlus_33_35   = 0, hfSumMinus_33_35  = 0;
  float hfSumPlus_33_36   = 0, hfSumMinus_33_36  = 0;
  float hfSumPlus_33_37   = 0, hfSumMinus_33_37  = 0;
  float hfSumPlus_33_38   = 0, hfSumMinus_33_38  = 0;
  float hfSumPlus_33_39   = 0, hfSumMinus_33_39  = 0;
  float hfSumPlus_33_40   = 0, hfSumMinus_33_40  = 0;
  float hfSumPlus_33_41   = 0, hfSumMinus_33_41  = 0;
  float hfSumPlus_33_42   = 0, hfSumMinus_33_42  = 0;
  float hfSumPlus_33_43   = 0, hfSumMinus_33_43  = 0;
  float hfSumPlus_33_44   = 0, hfSumMinus_33_44  = 0;
  float hfSumPlus_33_45   = 0, hfSumMinus_33_45  = 0;
  float hfSumPlus_33_46   = 0, hfSumMinus_33_46  = 0;
  float hfSumPlus_33_47   = 0, hfSumMinus_33_47  = 0;
  float hfSumPlus_33_48   = 0, hfSumMinus_33_48  = 0;
  float hfSumPlus_33_49   = 0, hfSumMinus_33_49  = 0;
  float hfSumPlus_33_50   = 0, hfSumMinus_33_50  = 0;
  float hfSumPlus_33_51   = 0, hfSumMinus_33_51  = 0;
  float hfSumPlus_33_52   = 0, hfSumMinus_33_52  = 0;

  float hfSumPlus_34_35   = 0, hfSumMinus_34_35  = 0;
  float hfSumPlus_34_36   = 0, hfSumMinus_34_36  = 0;
  float hfSumPlus_34_37   = 0, hfSumMinus_34_37  = 0;
  float hfSumPlus_34_38   = 0, hfSumMinus_34_38  = 0;
  float hfSumPlus_34_39   = 0, hfSumMinus_34_39  = 0;
  float hfSumPlus_34_40   = 0, hfSumMinus_34_40  = 0;
  float hfSumPlus_34_41   = 0, hfSumMinus_34_41  = 0;
  float hfSumPlus_34_42   = 0, hfSumMinus_34_42  = 0;
  float hfSumPlus_34_43   = 0, hfSumMinus_34_43  = 0;
  float hfSumPlus_34_44   = 0, hfSumMinus_34_44  = 0;
  float hfSumPlus_34_45   = 0, hfSumMinus_34_45  = 0;
  float hfSumPlus_34_46   = 0, hfSumMinus_34_46  = 0;
  float hfSumPlus_34_47   = 0, hfSumMinus_34_47  = 0;
  float hfSumPlus_34_48   = 0, hfSumMinus_34_48  = 0;
  float hfSumPlus_34_49   = 0, hfSumMinus_34_49  = 0;
  float hfSumPlus_34_50   = 0, hfSumMinus_34_50  = 0;
  float hfSumPlus_34_51   = 0, hfSumMinus_34_51  = 0;
  float hfSumPlus_34_52   = 0, hfSumMinus_34_52  = 0;
  
  float hfSumPlus_35_36   = 0, hfSumMinus_35_36  = 0;
  float hfSumPlus_35_37   = 0, hfSumMinus_35_37  = 0;
  float hfSumPlus_35_38   = 0, hfSumMinus_35_38  = 0;
  float hfSumPlus_35_39   = 0, hfSumMinus_35_39  = 0;
  float hfSumPlus_35_40   = 0, hfSumMinus_35_40  = 0;
  float hfSumPlus_35_41   = 0, hfSumMinus_35_41  = 0;
  float hfSumPlus_35_42   = 0, hfSumMinus_35_42  = 0;
  float hfSumPlus_35_43   = 0, hfSumMinus_35_43  = 0;
  float hfSumPlus_35_44   = 0, hfSumMinus_35_44  = 0;
  float hfSumPlus_35_45   = 0, hfSumMinus_35_45  = 0;
  float hfSumPlus_35_46   = 0, hfSumMinus_35_46  = 0;
  float hfSumPlus_35_47   = 0, hfSumMinus_35_47  = 0;
  float hfSumPlus_35_48   = 0, hfSumMinus_35_48  = 0;
  float hfSumPlus_35_49   = 0, hfSumMinus_35_49  = 0;
  float hfSumPlus_35_50   = 0, hfSumMinus_35_50  = 0;
  float hfSumPlus_35_51   = 0, hfSumMinus_35_51  = 0;
  float hfSumPlus_35_52   = 0, hfSumMinus_35_52  = 0;

  float hfSumPlus_36_37   = 0, hfSumMinus_36_37  = 0;
  float hfSumPlus_36_38   = 0, hfSumMinus_36_38  = 0;
  float hfSumPlus_36_39   = 0, hfSumMinus_36_39  = 0;
  float hfSumPlus_36_40   = 0, hfSumMinus_36_40  = 0;
  float hfSumPlus_36_41   = 0, hfSumMinus_36_41  = 0;
  float hfSumPlus_36_42   = 0, hfSumMinus_36_42  = 0;
  float hfSumPlus_36_43   = 0, hfSumMinus_36_43  = 0;
  float hfSumPlus_36_44   = 0, hfSumMinus_36_44  = 0;
  float hfSumPlus_36_45   = 0, hfSumMinus_36_45  = 0;
  float hfSumPlus_36_46   = 0, hfSumMinus_36_46  = 0;
  float hfSumPlus_36_47   = 0, hfSumMinus_36_47  = 0;
  float hfSumPlus_36_48   = 0, hfSumMinus_36_48  = 0;
  float hfSumPlus_36_49   = 0, hfSumMinus_36_49  = 0;
  float hfSumPlus_36_50   = 0, hfSumMinus_36_50  = 0;
  float hfSumPlus_36_51   = 0, hfSumMinus_36_51  = 0;
  float hfSumPlus_36_52   = 0, hfSumMinus_36_52  = 0;

  float hfSumPlus_37_38   = 0, hfSumMinus_37_38  = 0;
  float hfSumPlus_37_39   = 0, hfSumMinus_37_39  = 0;
  float hfSumPlus_37_40   = 0, hfSumMinus_37_40  = 0;
  float hfSumPlus_37_41   = 0, hfSumMinus_37_41  = 0;
  float hfSumPlus_37_42   = 0, hfSumMinus_37_42  = 0;
  float hfSumPlus_37_43   = 0, hfSumMinus_37_43  = 0;
  float hfSumPlus_37_44   = 0, hfSumMinus_37_44  = 0;
  float hfSumPlus_37_45   = 0, hfSumMinus_37_45  = 0;
  float hfSumPlus_37_46   = 0, hfSumMinus_37_46  = 0;
  float hfSumPlus_37_47   = 0, hfSumMinus_37_47  = 0;
  float hfSumPlus_37_48   = 0, hfSumMinus_37_48  = 0;
  float hfSumPlus_37_49   = 0, hfSumMinus_37_49  = 0;
  float hfSumPlus_37_50   = 0, hfSumMinus_37_50  = 0;
  float hfSumPlus_37_51   = 0, hfSumMinus_37_51  = 0;
  float hfSumPlus_37_52   = 0, hfSumMinus_37_52  = 0;

  float hfSumPlus_38_39   = 0, hfSumMinus_38_39  = 0;
  float hfSumPlus_38_40   = 0, hfSumMinus_38_40  = 0;
  float hfSumPlus_38_41   = 0, hfSumMinus_38_41  = 0;
  float hfSumPlus_38_42   = 0, hfSumMinus_38_42  = 0;
  float hfSumPlus_38_43   = 0, hfSumMinus_38_43  = 0;
  float hfSumPlus_38_44   = 0, hfSumMinus_38_44  = 0;
  float hfSumPlus_38_45   = 0, hfSumMinus_38_45  = 0;
  float hfSumPlus_38_46   = 0, hfSumMinus_38_46  = 0;
  float hfSumPlus_38_47   = 0, hfSumMinus_38_47  = 0;
  float hfSumPlus_38_48   = 0, hfSumMinus_38_48  = 0;
  float hfSumPlus_38_49   = 0, hfSumMinus_38_49  = 0;
  float hfSumPlus_38_50   = 0, hfSumMinus_38_50  = 0;
  float hfSumPlus_38_51   = 0, hfSumMinus_38_51  = 0;
  float hfSumPlus_38_52   = 0, hfSumMinus_38_52  = 0;

  float hfSumPlus_39_40   = 0, hfSumMinus_39_40  = 0;
  float hfSumPlus_39_41   = 0, hfSumMinus_39_41  = 0;
  float hfSumPlus_39_42   = 0, hfSumMinus_39_42  = 0;
  float hfSumPlus_39_43   = 0, hfSumMinus_39_43  = 0;
  float hfSumPlus_39_44   = 0, hfSumMinus_39_44  = 0;
  float hfSumPlus_39_45   = 0, hfSumMinus_39_45  = 0;
  float hfSumPlus_39_46   = 0, hfSumMinus_39_46  = 0;
  float hfSumPlus_39_47   = 0, hfSumMinus_39_47  = 0;
  float hfSumPlus_39_48   = 0, hfSumMinus_39_48  = 0;
  float hfSumPlus_39_49   = 0, hfSumMinus_39_49  = 0;
  float hfSumPlus_39_50   = 0, hfSumMinus_39_50  = 0;
  float hfSumPlus_39_51   = 0, hfSumMinus_39_51  = 0;
  float hfSumPlus_39_52   = 0, hfSumMinus_39_52  = 0;

  float hfSumPlus_40_41   = 0, hfSumMinus_40_41  = 0;
  float hfSumPlus_40_42   = 0, hfSumMinus_40_42  = 0;
  float hfSumPlus_40_43   = 0, hfSumMinus_40_43  = 0;
  float hfSumPlus_40_44   = 0, hfSumMinus_40_44  = 0;
  float hfSumPlus_40_45   = 0, hfSumMinus_40_45  = 0;
  float hfSumPlus_40_46   = 0, hfSumMinus_40_46  = 0;
  float hfSumPlus_40_47   = 0, hfSumMinus_40_47  = 0;
  float hfSumPlus_40_48   = 0, hfSumMinus_40_48  = 0;
  float hfSumPlus_40_49   = 0, hfSumMinus_40_49  = 0;
  float hfSumPlus_40_50   = 0, hfSumMinus_40_50  = 0;
  float hfSumPlus_40_51   = 0, hfSumMinus_40_51  = 0;
  float hfSumPlus_40_52   = 0, hfSumMinus_40_52  = 0;

  float hfSumPlus_41_42   = 0, hfSumMinus_41_42  = 0;
  float hfSumPlus_41_43   = 0, hfSumMinus_41_43  = 0;
  float hfSumPlus_41_44   = 0, hfSumMinus_41_44  = 0;
  float hfSumPlus_41_45   = 0, hfSumMinus_41_45  = 0;
  float hfSumPlus_41_46   = 0, hfSumMinus_41_46  = 0;
  float hfSumPlus_41_47   = 0, hfSumMinus_41_47  = 0;
  float hfSumPlus_41_48   = 0, hfSumMinus_41_48  = 0;
  float hfSumPlus_41_49   = 0, hfSumMinus_41_49  = 0;
  float hfSumPlus_41_50   = 0, hfSumMinus_41_50  = 0;
  float hfSumPlus_41_51   = 0, hfSumMinus_41_51  = 0;
  float hfSumPlus_41_52   = 0, hfSumMinus_41_52  = 0;

  float hfSumPlus_42_43   = 0, hfSumMinus_42_43  = 0;
  float hfSumPlus_42_44   = 0, hfSumMinus_42_44  = 0;
  float hfSumPlus_42_45   = 0, hfSumMinus_42_45  = 0;
  float hfSumPlus_42_46   = 0, hfSumMinus_42_46  = 0;
  float hfSumPlus_42_47   = 0, hfSumMinus_42_47  = 0;
  float hfSumPlus_42_48   = 0, hfSumMinus_42_48  = 0;
  float hfSumPlus_42_49   = 0, hfSumMinus_42_49  = 0;
  float hfSumPlus_42_50   = 0, hfSumMinus_42_50  = 0;
  float hfSumPlus_42_51   = 0, hfSumMinus_42_51  = 0;
  float hfSumPlus_42_52   = 0, hfSumMinus_42_52  = 0;

  float hfSumPlus_43_44   = 0, hfSumMinus_43_44  = 0;
  float hfSumPlus_43_45   = 0, hfSumMinus_43_45  = 0;
  float hfSumPlus_43_46   = 0, hfSumMinus_43_46  = 0;
  float hfSumPlus_43_47   = 0, hfSumMinus_43_47  = 0;
  float hfSumPlus_43_48   = 0, hfSumMinus_43_48  = 0;
  float hfSumPlus_43_49   = 0, hfSumMinus_43_49  = 0;
  float hfSumPlus_43_50   = 0, hfSumMinus_43_50  = 0;
  float hfSumPlus_43_51   = 0, hfSumMinus_43_51  = 0;
  float hfSumPlus_43_52   = 0, hfSumMinus_43_52  = 0;

  float hfSumPlus_44_45   = 0, hfSumMinus_44_45  = 0;
  float hfSumPlus_44_46   = 0, hfSumMinus_44_46  = 0;
  float hfSumPlus_44_47   = 0, hfSumMinus_44_47  = 0;
  float hfSumPlus_44_48   = 0, hfSumMinus_44_48  = 0;
  float hfSumPlus_44_49   = 0, hfSumMinus_44_49  = 0;
  float hfSumPlus_44_50   = 0, hfSumMinus_44_50  = 0;
  float hfSumPlus_44_51   = 0, hfSumMinus_44_51  = 0;
  float hfSumPlus_44_52   = 0, hfSumMinus_44_52  = 0;

  float hfSumPlus_45_46   = 0, hfSumMinus_45_46  = 0;
  float hfSumPlus_45_47   = 0, hfSumMinus_45_47  = 0;
  float hfSumPlus_45_48   = 0, hfSumMinus_45_48  = 0;
  float hfSumPlus_45_49   = 0, hfSumMinus_45_49  = 0;
  float hfSumPlus_45_50   = 0, hfSumMinus_45_50  = 0;
  float hfSumPlus_45_51   = 0, hfSumMinus_45_51  = 0;
  float hfSumPlus_45_52   = 0, hfSumMinus_45_52  = 0;

  float hfSumPlus_46_47   = 0, hfSumMinus_46_47  = 0;
  float hfSumPlus_46_48   = 0, hfSumMinus_46_48  = 0;
  float hfSumPlus_46_49   = 0, hfSumMinus_46_49  = 0;
  float hfSumPlus_46_50   = 0, hfSumMinus_46_50  = 0;
  float hfSumPlus_46_51   = 0, hfSumMinus_46_51  = 0;
  float hfSumPlus_46_52   = 0, hfSumMinus_46_52  = 0;

  float hfSumPlus_47_48   = 0, hfSumMinus_47_48  = 0;
  float hfSumPlus_47_49   = 0, hfSumMinus_47_49  = 0;
  float hfSumPlus_47_50   = 0, hfSumMinus_47_50  = 0;
  float hfSumPlus_47_51   = 0, hfSumMinus_47_51  = 0;
  float hfSumPlus_47_52   = 0, hfSumMinus_47_52  = 0;

  float hfSumPlus_48_49   = 0, hfSumMinus_48_49  = 0;
  float hfSumPlus_48_50   = 0, hfSumMinus_48_50  = 0;
  float hfSumPlus_48_51   = 0, hfSumMinus_48_51  = 0;
  float hfSumPlus_48_52   = 0, hfSumMinus_48_52  = 0;

  float hfSumPlus_49_50   = 0, hfSumMinus_49_50  = 0;
  float hfSumPlus_49_51   = 0, hfSumMinus_49_51  = 0;
  float hfSumPlus_49_52   = 0, hfSumMinus_49_52  = 0;

  float hfSumPlus_50_51   = 0, hfSumMinus_50_51  = 0;
  float hfSumPlus_50_52   = 0, hfSumMinus_50_52  = 0;

  float hfSumPlus_51_52   = 0, hfSumMinus_51_52  = 0;


  float minOfHFSums_30_31 = -1, minOfHFSums_30_32 = -1, minOfHFSums_30_33 = -1;
  float minOfHFSums_30_34 = -1, minOfHFSums_30_35 = -1, minOfHFSums_30_36 = -1;
  float minOfHFSums_30_37 = -1, minOfHFSums_30_38 = -1, minOfHFSums_30_39 = -1;
  float minOfHFSums_30_40 = -1, minOfHFSums_30_41 = -1, minOfHFSums_30_42 = -1;
  float minOfHFSums_30_43 = -1, minOfHFSums_30_44 = -1, minOfHFSums_30_45 = -1;
  float minOfHFSums_30_46 = -1, minOfHFSums_30_47 = -1, minOfHFSums_30_48 = -1;
  float minOfHFSums_30_49 = -1, minOfHFSums_30_50 = -1, minOfHFSums_30_51 = -1;
  float minOfHFSums_30_52 = -1;

  float minOfHFSums_31_32 = -1, minOfHFSums_31_33 = -1;
  float minOfHFSums_31_34 = -1, minOfHFSums_31_35 = -1, minOfHFSums_31_36 = -1;
  float minOfHFSums_31_37 = -1, minOfHFSums_31_38 = -1, minOfHFSums_31_39 = -1;
  float minOfHFSums_31_40 = -1, minOfHFSums_31_41 = -1, minOfHFSums_31_42 = -1;
  float minOfHFSums_31_43 = -1, minOfHFSums_31_44 = -1, minOfHFSums_31_45 = -1;
  float minOfHFSums_31_46 = -1, minOfHFSums_31_47 = -1, minOfHFSums_31_48 = -1;
  float minOfHFSums_31_49 = -1, minOfHFSums_31_50 = -1, minOfHFSums_31_51 = -1;
  float minOfHFSums_31_52 = -1;

  float minOfHFSums_32_33 = -1, minOfHFSums_32_34 = -1, minOfHFSums_32_35 = -1;
  float minOfHFSums_32_36 = -1, minOfHFSums_32_37 = -1, minOfHFSums_32_38 = -1;
  float minOfHFSums_32_39 = -1, minOfHFSums_32_40 = -1, minOfHFSums_32_41 = -1;
  float minOfHFSums_32_42 = -1, minOfHFSums_32_43 = -1, minOfHFSums_32_44 = -1;
  float minOfHFSums_32_45 = -1, minOfHFSums_32_46 = -1, minOfHFSums_32_47 = -1;
  float minOfHFSums_32_48 = -1, minOfHFSums_32_49 = -1, minOfHFSums_32_50 = -1;
  float minOfHFSums_32_51 = -1, minOfHFSums_32_52 = -1;

  float minOfHFSums_33_34 = -1, minOfHFSums_33_35 = -1, minOfHFSums_33_36 = -1;
  float minOfHFSums_33_37 = -1, minOfHFSums_33_38 = -1, minOfHFSums_33_39 = -1;
  float minOfHFSums_33_40 = -1, minOfHFSums_33_41 = -1, minOfHFSums_33_42 = -1;
  float minOfHFSums_33_43 = -1, minOfHFSums_33_44 = -1, minOfHFSums_33_45 = -1;
  float minOfHFSums_33_46 = -1, minOfHFSums_33_47 = -1, minOfHFSums_33_48 = -1;
  float minOfHFSums_33_49 = -1, minOfHFSums_33_50 = -1, minOfHFSums_33_51 = -1;
  float minOfHFSums_33_52 = -1;

  float minOfHFSums_34_35 = -1, minOfHFSums_34_36 = -1, minOfHFSums_34_37 = -1;
  float minOfHFSums_34_38 = -1, minOfHFSums_34_39 = -1, minOfHFSums_34_40 = -1;
  float minOfHFSums_34_41 = -1, minOfHFSums_34_42 = -1, minOfHFSums_34_43 = -1;
  float minOfHFSums_34_44 = -1, minOfHFSums_34_45 = -1, minOfHFSums_34_46 = -1;
  float minOfHFSums_34_47 = -1, minOfHFSums_34_48 = -1, minOfHFSums_34_49 = -1;
  float minOfHFSums_34_50 = -1, minOfHFSums_34_51 = -1, minOfHFSums_34_52 = -1;

  float minOfHFSums_35_36 = -1, minOfHFSums_35_37 = -1, minOfHFSums_35_38 = -1;
  float minOfHFSums_35_39 = -1, minOfHFSums_35_40 = -1, minOfHFSums_35_41 = -1;
  float minOfHFSums_35_42 = -1, minOfHFSums_35_43 = -1, minOfHFSums_35_44 = -1;
  float minOfHFSums_35_45 = -1, minOfHFSums_35_46 = -1, minOfHFSums_35_47 = -1;
  float minOfHFSums_35_48 = -1, minOfHFSums_35_49 = -1, minOfHFSums_35_50 = -1;
  float minOfHFSums_35_51 = -1, minOfHFSums_35_52 = -1;

  float minOfHFSums_36_37 = -1, minOfHFSums_36_38 = -1, minOfHFSums_36_39 = -1;
  float minOfHFSums_36_40 = -1, minOfHFSums_36_41 = -1, minOfHFSums_36_42 = -1;
  float minOfHFSums_36_43 = -1, minOfHFSums_36_44 = -1, minOfHFSums_36_45 = -1;
  float minOfHFSums_36_46 = -1, minOfHFSums_36_47 = -1, minOfHFSums_36_48 = -1;
  float minOfHFSums_36_49 = -1, minOfHFSums_36_50 = -1, minOfHFSums_36_51 = -1;
  float minOfHFSums_36_52 = -1;

  float minOfHFSums_37_38 = -1, minOfHFSums_37_39 = -1, minOfHFSums_37_40 = -1;
  float minOfHFSums_37_41 = -1, minOfHFSums_37_42 = -1, minOfHFSums_37_43 = -1;
  float minOfHFSums_37_44 = -1, minOfHFSums_37_45 = -1, minOfHFSums_37_46 = -1;
  float minOfHFSums_37_47 = -1, minOfHFSums_37_48 = -1, minOfHFSums_37_49 = -1;
  float minOfHFSums_37_50 = -1, minOfHFSums_37_51 = -1, minOfHFSums_37_52 = -1;

  float minOfHFSums_38_39 = -1, minOfHFSums_38_40 = -1, minOfHFSums_38_41 = -1;
  float minOfHFSums_38_42 = -1, minOfHFSums_38_43 = -1, minOfHFSums_38_44 = -1;
  float minOfHFSums_38_45 = -1, minOfHFSums_38_46 = -1, minOfHFSums_38_47 = -1;
  float minOfHFSums_38_48 = -1, minOfHFSums_38_49 = -1, minOfHFSums_38_50 = -1;
  float minOfHFSums_38_51 = -1, minOfHFSums_38_52 = -1;

  float minOfHFSums_39_40 = -1, minOfHFSums_39_41 = -1, minOfHFSums_39_42 = -1;
  float minOfHFSums_39_43 = -1, minOfHFSums_39_44 = -1, minOfHFSums_39_45 = -1;
  float minOfHFSums_39_46 = -1, minOfHFSums_39_47 = -1, minOfHFSums_39_48 = -1;
  float minOfHFSums_39_49 = -1, minOfHFSums_39_50 = -1, minOfHFSums_39_51 = -1;
  float minOfHFSums_39_52 = -1;

  float minOfHFSums_40_41 = -1, minOfHFSums_40_42 = -1, minOfHFSums_40_43 = -1;
  float minOfHFSums_40_44 = -1, minOfHFSums_40_45 = -1, minOfHFSums_40_46 = -1;
  float minOfHFSums_40_47 = -1, minOfHFSums_40_48 = -1, minOfHFSums_40_49 = -1;
  float minOfHFSums_40_50 = -1, minOfHFSums_40_51 = -1, minOfHFSums_40_52 = -1;

  float minOfHFSums_41_42 = -1, minOfHFSums_41_43 = -1, minOfHFSums_41_44 = -1;
  float minOfHFSums_41_45 = -1, minOfHFSums_41_46 = -1, minOfHFSums_41_47 = -1;
  float minOfHFSums_41_48 = -1, minOfHFSums_41_49 = -1, minOfHFSums_41_50 = -1;
  float minOfHFSums_41_51 = -1, minOfHFSums_41_52 = -1;

  float minOfHFSums_42_43 = -1, minOfHFSums_42_44 = -1, minOfHFSums_42_45 = -1;
  float minOfHFSums_42_46 = -1, minOfHFSums_42_47 = -1, minOfHFSums_42_48 = -1;
  float minOfHFSums_42_49 = -1, minOfHFSums_42_50 = -1, minOfHFSums_42_51 = -1;
  float minOfHFSums_42_52 = -1;

  float minOfHFSums_43_44 = -1, minOfHFSums_43_45 = -1, minOfHFSums_43_46 = -1;
  float minOfHFSums_43_47 = -1, minOfHFSums_43_48 = -1, minOfHFSums_43_49 = -1;
  float minOfHFSums_43_50 = -1, minOfHFSums_43_51 = -1, minOfHFSums_43_52 = -1;

  float minOfHFSums_44_45 = -1, minOfHFSums_44_46 = -1, minOfHFSums_44_47 = -1;
  float minOfHFSums_44_48 = -1, minOfHFSums_44_49 = -1, minOfHFSums_44_50 = -1;
  float minOfHFSums_44_51 = -1, minOfHFSums_44_52 = -1;

  float minOfHFSums_45_46 = -1, minOfHFSums_45_47 = -1, minOfHFSums_45_48 = -1;
  float minOfHFSums_45_49 = -1, minOfHFSums_45_50 = -1, minOfHFSums_45_51 = -1;
  float minOfHFSums_45_52 = -1;

  float minOfHFSums_46_47 = -1, minOfHFSums_46_48 = -1, minOfHFSums_46_49 = -1;
  float minOfHFSums_46_50 = -1, minOfHFSums_46_51 = -1, minOfHFSums_46_52 = -1;

  float minOfHFSums_47_48 = -1, minOfHFSums_47_49 = -1, minOfHFSums_47_50 = -1;
  float minOfHFSums_47_51 = -1, minOfHFSums_47_52 = -1;

  float minOfHFSums_48_49 = -1, minOfHFSums_48_50 = -1, minOfHFSums_48_51 = -1;
  float minOfHFSums_48_52 = -1;

  float minOfHFSums_49_50 = -1, minOfHFSums_49_51 = -1, minOfHFSums_49_52 = -1;

  float minOfHFSums_50_51 = -1, minOfHFSums_50_52 = -1;

  float minOfHFSums_51_52 = -1;
  // New floats --end

  tree->Branch("rapidityGapPlus",  &rapidityGapPlus,  "rapidityGapPlus/F");
  tree->Branch("rapidityGapMinus", &rapidityGapMinus, "rapidityGapMinus/F");
  tree->Branch("hfSumPlus",        &hfSumPlus,        "hfSumPlus/F");
  tree->Branch("hfSumMinus",       &hfSumMinus,       "hfSumMinus/F");
  tree->Branch("nchPV",            &nchPV,            "nchPV/F");
  tree->Branch("nForwardJets",     &nForwardJets,     "nForwardJets/F");
  tree->Branch("nTotalJets",       &nTotalJets,       "nTotalJets/F");
  tree->Branch("leadLightJetPt",   &leadLightJetPt,   "leadLightJetPt/F");
  tree->Branch("maxHFEnergy",      &maxHFEnergy,      "maxHFEnergy/F");
  tree->Branch("minHFEnergy",      &minHFEnergy,      "minHFEnergy/F");
  tree->Branch("MET",              &METval,           "MET/F");

  tree->Branch("maxHFEnergyPlus",   &maxHFEnergyPlus,   "maxHFEnergyPlus/F");
  tree->Branch("maxHFEnergyMinus",  &maxHFEnergyMinus,  "maxHFEnergyMinus/F");
  tree->Branch("minOfMaxHF",        &minOfMaxHF,        "minOfMaxHF/F");
  tree->Branch("minOfHFSums",       &minOfHFSums,       "minOfHFSums/F");

  // minofmaxhf
  tree->Branch("minOfMaxHF_30_31", &minOfMaxHF_30_31, "minOfMaxHF_30_31/F");
  tree->Branch("minOfMaxHF_30_32", &minOfMaxHF_30_32, "minOfMaxHF_30_32/F");
  tree->Branch("minOfMaxHF_30_33", &minOfMaxHF_30_33, "minOfMaxHF_30_33/F");
  tree->Branch("minOfMaxHF_30_34", &minOfMaxHF_30_34, "minOfMaxHF_30_34/F");
  tree->Branch("minOfMaxHF_30_35", &minOfMaxHF_30_35, "minOfMaxHF_30_35/F");
  tree->Branch("minOfMaxHF_30_36", &minOfMaxHF_30_36, "minOfMaxHF_30_36/F");
  tree->Branch("minOfMaxHF_30_37", &minOfMaxHF_30_37, "minOfMaxHF_30_37/F");
  tree->Branch("minOfMaxHF_30_38", &minOfMaxHF_30_38, "minOfMaxHF_30_38/F");
  tree->Branch("minOfMaxHF_30_39", &minOfMaxHF_30_39, "minOfMaxHF_30_39/F");
  tree->Branch("minOfMaxHF_30_40", &minOfMaxHF_30_40, "minOfMaxHF_30_40/F");
  tree->Branch("minOfMaxHF_30_41", &minOfMaxHF_30_41, "minOfMaxHF_30_41/F");
  tree->Branch("minOfMaxHF_30_42", &minOfMaxHF_30_42, "minOfMaxHF_30_42/F");
  tree->Branch("minOfMaxHF_30_43", &minOfMaxHF_30_43, "minOfMaxHF_30_43/F");
  tree->Branch("minOfMaxHF_30_44", &minOfMaxHF_30_44, "minOfMaxHF_30_44/F");
  tree->Branch("minOfMaxHF_30_45", &minOfMaxHF_30_45, "minOfMaxHF_30_45/F");
  tree->Branch("minOfMaxHF_30_46", &minOfMaxHF_30_46, "minOfMaxHF_30_46/F");
  tree->Branch("minOfMaxHF_30_47", &minOfMaxHF_30_47, "minOfMaxHF_30_47/F");
  tree->Branch("minOfMaxHF_30_48", &minOfMaxHF_30_48, "minOfMaxHF_30_48/F");
  tree->Branch("minOfMaxHF_30_49", &minOfMaxHF_30_49, "minOfMaxHF_30_49/F");
  tree->Branch("minOfMaxHF_30_50", &minOfMaxHF_30_50, "minOfMaxHF_30_50/F");
  tree->Branch("minOfMaxHF_30_51", &minOfMaxHF_30_51, "minOfMaxHF_30_51/F");
  tree->Branch("minOfMaxHF_30_52", &minOfMaxHF_30_52, "minOfMaxHF_30_52/F");

  tree->Branch("minOfMaxHF_31_32", &minOfMaxHF_31_32, "minOfMaxHF_31_32/F");
  tree->Branch("minOfMaxHF_31_33", &minOfMaxHF_31_33, "minOfMaxHF_31_33/F");
  tree->Branch("minOfMaxHF_31_34", &minOfMaxHF_31_34, "minOfMaxHF_31_34/F");
  tree->Branch("minOfMaxHF_31_35", &minOfMaxHF_31_35, "minOfMaxHF_31_35/F");
  tree->Branch("minOfMaxHF_31_36", &minOfMaxHF_31_36, "minOfMaxHF_31_36/F");
  tree->Branch("minOfMaxHF_31_37", &minOfMaxHF_31_37, "minOfMaxHF_31_37/F");
  tree->Branch("minOfMaxHF_31_38", &minOfMaxHF_31_38, "minOfMaxHF_31_38/F");
  tree->Branch("minOfMaxHF_31_39", &minOfMaxHF_31_39, "minOfMaxHF_31_39/F");
  tree->Branch("minOfMaxHF_31_40", &minOfMaxHF_31_40, "minOfMaxHF_31_40/F");
  tree->Branch("minOfMaxHF_31_41", &minOfMaxHF_31_41, "minOfMaxHF_31_41/F");
  tree->Branch("minOfMaxHF_31_42", &minOfMaxHF_31_42, "minOfMaxHF_31_42/F");
  tree->Branch("minOfMaxHF_31_43", &minOfMaxHF_31_43, "minOfMaxHF_31_43/F");
  tree->Branch("minOfMaxHF_31_44", &minOfMaxHF_31_44, "minOfMaxHF_31_44/F");
  tree->Branch("minOfMaxHF_31_45", &minOfMaxHF_31_45, "minOfMaxHF_31_45/F");
  tree->Branch("minOfMaxHF_31_46", &minOfMaxHF_31_46, "minOfMaxHF_31_46/F");
  tree->Branch("minOfMaxHF_31_47", &minOfMaxHF_31_47, "minOfMaxHF_31_47/F");
  tree->Branch("minOfMaxHF_31_48", &minOfMaxHF_31_48, "minOfMaxHF_31_48/F");
  tree->Branch("minOfMaxHF_31_49", &minOfMaxHF_31_49, "minOfMaxHF_31_49/F");
  tree->Branch("minOfMaxHF_31_50", &minOfMaxHF_31_50, "minOfMaxHF_31_50/F");
  tree->Branch("minOfMaxHF_31_51", &minOfMaxHF_31_51, "minOfMaxHF_31_51/F");
  tree->Branch("minOfMaxHF_31_52", &minOfMaxHF_31_52, "minOfMaxHF_31_52/F");

  tree->Branch("minOfMaxHF_32_33", &minOfMaxHF_32_33, "minOfMaxHF_32_33/F");
  tree->Branch("minOfMaxHF_32_34", &minOfMaxHF_32_34, "minOfMaxHF_32_34/F");
  tree->Branch("minOfMaxHF_32_35", &minOfMaxHF_32_35, "minOfMaxHF_32_35/F");
  tree->Branch("minOfMaxHF_32_36", &minOfMaxHF_32_36, "minOfMaxHF_32_36/F");
  tree->Branch("minOfMaxHF_32_37", &minOfMaxHF_32_37, "minOfMaxHF_32_37/F");
  tree->Branch("minOfMaxHF_32_38", &minOfMaxHF_32_38, "minOfMaxHF_32_38/F");
  tree->Branch("minOfMaxHF_32_39", &minOfMaxHF_32_39, "minOfMaxHF_32_39/F");
  tree->Branch("minOfMaxHF_32_40", &minOfMaxHF_32_40, "minOfMaxHF_32_40/F");
  tree->Branch("minOfMaxHF_32_41", &minOfMaxHF_32_41, "minOfMaxHF_32_41/F");
  tree->Branch("minOfMaxHF_32_42", &minOfMaxHF_32_42, "minOfMaxHF_32_42/F");
  tree->Branch("minOfMaxHF_32_43", &minOfMaxHF_32_43, "minOfMaxHF_32_43/F");
  tree->Branch("minOfMaxHF_32_44", &minOfMaxHF_32_44, "minOfMaxHF_32_44/F");
  tree->Branch("minOfMaxHF_32_45", &minOfMaxHF_32_45, "minOfMaxHF_32_45/F");
  tree->Branch("minOfMaxHF_32_46", &minOfMaxHF_32_46, "minOfMaxHF_32_46/F");
  tree->Branch("minOfMaxHF_32_47", &minOfMaxHF_32_47, "minOfMaxHF_32_47/F");
  tree->Branch("minOfMaxHF_32_48", &minOfMaxHF_32_48, "minOfMaxHF_32_48/F");
  tree->Branch("minOfMaxHF_32_49", &minOfMaxHF_32_49, "minOfMaxHF_32_49/F");
  tree->Branch("minOfMaxHF_32_50", &minOfMaxHF_32_50, "minOfMaxHF_32_50/F");
  tree->Branch("minOfMaxHF_32_51", &minOfMaxHF_32_51, "minOfMaxHF_32_51/F");
  tree->Branch("minOfMaxHF_32_52", &minOfMaxHF_32_52, "minOfMaxHF_32_52/F");

  tree->Branch("minOfMaxHF_33_34", &minOfMaxHF_33_34, "minOfMaxHF_33_34/F");
  tree->Branch("minOfMaxHF_33_35", &minOfMaxHF_33_35, "minOfMaxHF_33_35/F");
  tree->Branch("minOfMaxHF_33_36", &minOfMaxHF_33_36, "minOfMaxHF_33_36/F");
  tree->Branch("minOfMaxHF_33_37", &minOfMaxHF_33_37, "minOfMaxHF_33_37/F");
  tree->Branch("minOfMaxHF_33_38", &minOfMaxHF_33_38, "minOfMaxHF_33_38/F");
  tree->Branch("minOfMaxHF_33_39", &minOfMaxHF_33_39, "minOfMaxHF_33_39/F");
  tree->Branch("minOfMaxHF_33_40", &minOfMaxHF_33_40, "minOfMaxHF_33_40/F");
  tree->Branch("minOfMaxHF_33_41", &minOfMaxHF_33_41, "minOfMaxHF_33_41/F");
  tree->Branch("minOfMaxHF_33_42", &minOfMaxHF_33_42, "minOfMaxHF_33_42/F");
  tree->Branch("minOfMaxHF_33_43", &minOfMaxHF_33_43, "minOfMaxHF_33_43/F");
  tree->Branch("minOfMaxHF_33_44", &minOfMaxHF_33_44, "minOfMaxHF_33_44/F");
  tree->Branch("minOfMaxHF_33_45", &minOfMaxHF_33_45, "minOfMaxHF_33_45/F");
  tree->Branch("minOfMaxHF_33_46", &minOfMaxHF_33_46, "minOfMaxHF_33_46/F");
  tree->Branch("minOfMaxHF_33_47", &minOfMaxHF_33_47, "minOfMaxHF_33_47/F");
  tree->Branch("minOfMaxHF_33_48", &minOfMaxHF_33_48, "minOfMaxHF_33_48/F");
  tree->Branch("minOfMaxHF_33_49", &minOfMaxHF_33_49, "minOfMaxHF_33_49/F");
  tree->Branch("minOfMaxHF_33_50", &minOfMaxHF_33_50, "minOfMaxHF_33_50/F");
  tree->Branch("minOfMaxHF_33_51", &minOfMaxHF_33_51, "minOfMaxHF_33_51/F");
  tree->Branch("minOfMaxHF_33_52", &minOfMaxHF_33_52, "minOfMaxHF_33_52/F");

  tree->Branch("minOfMaxHF_34_35", &minOfMaxHF_34_35, "minOfMaxHF_34_35/F");
  tree->Branch("minOfMaxHF_34_36", &minOfMaxHF_34_36, "minOfMaxHF_34_36/F");
  tree->Branch("minOfMaxHF_34_37", &minOfMaxHF_34_37, "minOfMaxHF_34_37/F");
  tree->Branch("minOfMaxHF_34_38", &minOfMaxHF_34_38, "minOfMaxHF_34_38/F");
  tree->Branch("minOfMaxHF_34_39", &minOfMaxHF_34_39, "minOfMaxHF_34_39/F");
  tree->Branch("minOfMaxHF_34_40", &minOfMaxHF_34_40, "minOfMaxHF_34_40/F");
  tree->Branch("minOfMaxHF_34_41", &minOfMaxHF_34_41, "minOfMaxHF_34_41/F");
  tree->Branch("minOfMaxHF_34_42", &minOfMaxHF_34_42, "minOfMaxHF_34_42/F");
  tree->Branch("minOfMaxHF_34_43", &minOfMaxHF_34_43, "minOfMaxHF_34_43/F");
  tree->Branch("minOfMaxHF_34_44", &minOfMaxHF_34_44, "minOfMaxHF_34_44/F");
  tree->Branch("minOfMaxHF_34_45", &minOfMaxHF_34_45, "minOfMaxHF_34_45/F");
  tree->Branch("minOfMaxHF_34_46", &minOfMaxHF_34_46, "minOfMaxHF_34_46/F");
  tree->Branch("minOfMaxHF_34_47", &minOfMaxHF_34_47, "minOfMaxHF_34_47/F");
  tree->Branch("minOfMaxHF_34_48", &minOfMaxHF_34_48, "minOfMaxHF_34_48/F");
  tree->Branch("minOfMaxHF_34_49", &minOfMaxHF_34_49, "minOfMaxHF_34_49/F");
  tree->Branch("minOfMaxHF_34_50", &minOfMaxHF_34_50, "minOfMaxHF_34_50/F");
  tree->Branch("minOfMaxHF_34_51", &minOfMaxHF_34_51, "minOfMaxHF_34_51/F");
  tree->Branch("minOfMaxHF_34_52", &minOfMaxHF_34_52, "minOfMaxHF_34_52/F");

  tree->Branch("minOfMaxHF_35_36", &minOfMaxHF_35_36, "minOfMaxHF_35_36/F");
  tree->Branch("minOfMaxHF_35_37", &minOfMaxHF_35_37, "minOfMaxHF_35_37/F");
  tree->Branch("minOfMaxHF_35_38", &minOfMaxHF_35_38, "minOfMaxHF_35_38/F");
  tree->Branch("minOfMaxHF_35_39", &minOfMaxHF_35_39, "minOfMaxHF_35_39/F");
  tree->Branch("minOfMaxHF_35_40", &minOfMaxHF_35_40, "minOfMaxHF_35_40/F");
  tree->Branch("minOfMaxHF_35_41", &minOfMaxHF_35_41, "minOfMaxHF_35_41/F");
  tree->Branch("minOfMaxHF_35_42", &minOfMaxHF_35_42, "minOfMaxHF_35_42/F");
  tree->Branch("minOfMaxHF_35_43", &minOfMaxHF_35_43, "minOfMaxHF_35_43/F");
  tree->Branch("minOfMaxHF_35_44", &minOfMaxHF_35_44, "minOfMaxHF_35_44/F");
  tree->Branch("minOfMaxHF_35_45", &minOfMaxHF_35_45, "minOfMaxHF_35_45/F");
  tree->Branch("minOfMaxHF_35_46", &minOfMaxHF_35_46, "minOfMaxHF_35_46/F");
  tree->Branch("minOfMaxHF_35_47", &minOfMaxHF_35_47, "minOfMaxHF_35_47/F");
  tree->Branch("minOfMaxHF_35_48", &minOfMaxHF_35_48, "minOfMaxHF_35_48/F");
  tree->Branch("minOfMaxHF_35_49", &minOfMaxHF_35_49, "minOfMaxHF_35_49/F");
  tree->Branch("minOfMaxHF_35_50", &minOfMaxHF_35_50, "minOfMaxHF_35_50/F");
  tree->Branch("minOfMaxHF_35_51", &minOfMaxHF_35_51, "minOfMaxHF_35_51/F");
  tree->Branch("minOfMaxHF_35_52", &minOfMaxHF_35_52, "minOfMaxHF_35_52/F");

  tree->Branch("minOfMaxHF_36_37", &minOfMaxHF_36_37, "minOfMaxHF_36_37/F");
  tree->Branch("minOfMaxHF_36_38", &minOfMaxHF_36_38, "minOfMaxHF_36_38/F");
  tree->Branch("minOfMaxHF_36_39", &minOfMaxHF_36_39, "minOfMaxHF_36_39/F");
  tree->Branch("minOfMaxHF_36_40", &minOfMaxHF_36_40, "minOfMaxHF_36_40/F");
  tree->Branch("minOfMaxHF_36_41", &minOfMaxHF_36_41, "minOfMaxHF_36_41/F");
  tree->Branch("minOfMaxHF_36_42", &minOfMaxHF_36_42, "minOfMaxHF_36_42/F");
  tree->Branch("minOfMaxHF_36_43", &minOfMaxHF_36_43, "minOfMaxHF_36_43/F");
  tree->Branch("minOfMaxHF_36_44", &minOfMaxHF_36_44, "minOfMaxHF_36_44/F");
  tree->Branch("minOfMaxHF_36_45", &minOfMaxHF_36_45, "minOfMaxHF_36_45/F");
  tree->Branch("minOfMaxHF_36_46", &minOfMaxHF_36_46, "minOfMaxHF_36_46/F");
  tree->Branch("minOfMaxHF_36_47", &minOfMaxHF_36_47, "minOfMaxHF_36_47/F");
  tree->Branch("minOfMaxHF_36_48", &minOfMaxHF_36_48, "minOfMaxHF_36_48/F");
  tree->Branch("minOfMaxHF_36_49", &minOfMaxHF_36_49, "minOfMaxHF_36_49/F");
  tree->Branch("minOfMaxHF_36_50", &minOfMaxHF_36_50, "minOfMaxHF_36_50/F");
  tree->Branch("minOfMaxHF_36_51", &minOfMaxHF_36_51, "minOfMaxHF_36_51/F");
  tree->Branch("minOfMaxHF_36_52", &minOfMaxHF_36_52, "minOfMaxHF_36_52/F");

  tree->Branch("minOfMaxHF_37_38", &minOfMaxHF_37_38, "minOfMaxHF_37_38/F");
  tree->Branch("minOfMaxHF_37_39", &minOfMaxHF_37_39, "minOfMaxHF_37_39/F");
  tree->Branch("minOfMaxHF_37_40", &minOfMaxHF_37_40, "minOfMaxHF_37_40/F");
  tree->Branch("minOfMaxHF_37_41", &minOfMaxHF_37_41, "minOfMaxHF_37_41/F");
  tree->Branch("minOfMaxHF_37_42", &minOfMaxHF_37_42, "minOfMaxHF_37_42/F");
  tree->Branch("minOfMaxHF_37_43", &minOfMaxHF_37_43, "minOfMaxHF_37_43/F");
  tree->Branch("minOfMaxHF_37_44", &minOfMaxHF_37_44, "minOfMaxHF_37_44/F");
  tree->Branch("minOfMaxHF_37_45", &minOfMaxHF_37_45, "minOfMaxHF_37_45/F");
  tree->Branch("minOfMaxHF_37_46", &minOfMaxHF_37_46, "minOfMaxHF_37_46/F");
  tree->Branch("minOfMaxHF_37_47", &minOfMaxHF_37_47, "minOfMaxHF_37_47/F");
  tree->Branch("minOfMaxHF_37_48", &minOfMaxHF_37_48, "minOfMaxHF_37_48/F");
  tree->Branch("minOfMaxHF_37_49", &minOfMaxHF_37_49, "minOfMaxHF_37_49/F");
  tree->Branch("minOfMaxHF_37_50", &minOfMaxHF_37_50, "minOfMaxHF_37_50/F");
  tree->Branch("minOfMaxHF_37_51", &minOfMaxHF_37_51, "minOfMaxHF_37_51/F");
  tree->Branch("minOfMaxHF_37_52", &minOfMaxHF_37_52, "minOfMaxHF_37_52/F");

  tree->Branch("minOfMaxHF_38_39", &minOfMaxHF_38_39, "minOfMaxHF_38_39/F");
  tree->Branch("minOfMaxHF_38_40", &minOfMaxHF_38_40, "minOfMaxHF_38_40/F");
  tree->Branch("minOfMaxHF_38_41", &minOfMaxHF_38_41, "minOfMaxHF_38_41/F");
  tree->Branch("minOfMaxHF_38_42", &minOfMaxHF_38_42, "minOfMaxHF_38_42/F");
  tree->Branch("minOfMaxHF_38_43", &minOfMaxHF_38_43, "minOfMaxHF_38_43/F");
  tree->Branch("minOfMaxHF_38_44", &minOfMaxHF_38_44, "minOfMaxHF_38_44/F");
  tree->Branch("minOfMaxHF_38_45", &minOfMaxHF_38_45, "minOfMaxHF_38_45/F");
  tree->Branch("minOfMaxHF_38_46", &minOfMaxHF_38_46, "minOfMaxHF_38_46/F");
  tree->Branch("minOfMaxHF_38_47", &minOfMaxHF_38_47, "minOfMaxHF_38_47/F");
  tree->Branch("minOfMaxHF_38_48", &minOfMaxHF_38_48, "minOfMaxHF_38_48/F");
  tree->Branch("minOfMaxHF_38_49", &minOfMaxHF_38_49, "minOfMaxHF_38_49/F");
  tree->Branch("minOfMaxHF_38_50", &minOfMaxHF_38_50, "minOfMaxHF_38_50/F");
  tree->Branch("minOfMaxHF_38_51", &minOfMaxHF_38_51, "minOfMaxHF_38_51/F");
  tree->Branch("minOfMaxHF_38_52", &minOfMaxHF_38_52, "minOfMaxHF_38_52/F");

  tree->Branch("minOfMaxHF_39_40", &minOfMaxHF_39_40, "minOfMaxHF_39_40/F");
  tree->Branch("minOfMaxHF_39_41", &minOfMaxHF_39_41, "minOfMaxHF_39_41/F");
  tree->Branch("minOfMaxHF_39_42", &minOfMaxHF_39_42, "minOfMaxHF_39_42/F");
  tree->Branch("minOfMaxHF_39_43", &minOfMaxHF_39_43, "minOfMaxHF_39_43/F");
  tree->Branch("minOfMaxHF_39_44", &minOfMaxHF_39_44, "minOfMaxHF_39_44/F");
  tree->Branch("minOfMaxHF_39_45", &minOfMaxHF_39_45, "minOfMaxHF_39_45/F");
  tree->Branch("minOfMaxHF_39_46", &minOfMaxHF_39_46, "minOfMaxHF_39_46/F");
  tree->Branch("minOfMaxHF_39_47", &minOfMaxHF_39_47, "minOfMaxHF_39_47/F");
  tree->Branch("minOfMaxHF_39_48", &minOfMaxHF_39_48, "minOfMaxHF_39_48/F");
  tree->Branch("minOfMaxHF_39_49", &minOfMaxHF_39_49, "minOfMaxHF_39_49/F");
  tree->Branch("minOfMaxHF_39_50", &minOfMaxHF_39_50, "minOfMaxHF_39_50/F");
  tree->Branch("minOfMaxHF_39_51", &minOfMaxHF_39_51, "minOfMaxHF_39_51/F");
  tree->Branch("minOfMaxHF_39_52", &minOfMaxHF_39_52, "minOfMaxHF_39_52/F");

  tree->Branch("minOfMaxHF_40_41", &minOfMaxHF_40_41, "minOfMaxHF_40_41/F");
  tree->Branch("minOfMaxHF_40_42", &minOfMaxHF_40_42, "minOfMaxHF_40_42/F");
  tree->Branch("minOfMaxHF_40_43", &minOfMaxHF_40_43, "minOfMaxHF_40_43/F");
  tree->Branch("minOfMaxHF_40_44", &minOfMaxHF_40_44, "minOfMaxHF_40_44/F");
  tree->Branch("minOfMaxHF_40_45", &minOfMaxHF_40_45, "minOfMaxHF_40_45/F");
  tree->Branch("minOfMaxHF_40_46", &minOfMaxHF_40_46, "minOfMaxHF_40_46/F");
  tree->Branch("minOfMaxHF_40_47", &minOfMaxHF_40_47, "minOfMaxHF_40_47/F");
  tree->Branch("minOfMaxHF_40_48", &minOfMaxHF_40_48, "minOfMaxHF_40_48/F");
  tree->Branch("minOfMaxHF_40_49", &minOfMaxHF_40_49, "minOfMaxHF_40_49/F");
  tree->Branch("minOfMaxHF_40_50", &minOfMaxHF_40_50, "minOfMaxHF_40_50/F");
  tree->Branch("minOfMaxHF_40_51", &minOfMaxHF_40_51, "minOfMaxHF_40_51/F");
  tree->Branch("minOfMaxHF_40_52", &minOfMaxHF_40_52, "minOfMaxHF_40_52/F");

  tree->Branch("minOfMaxHF_41_42", &minOfMaxHF_41_42, "minOfMaxHF_41_42/F");
  tree->Branch("minOfMaxHF_41_43", &minOfMaxHF_41_43, "minOfMaxHF_41_43/F");
  tree->Branch("minOfMaxHF_41_44", &minOfMaxHF_41_44, "minOfMaxHF_41_44/F");
  tree->Branch("minOfMaxHF_41_45", &minOfMaxHF_41_45, "minOfMaxHF_41_45/F");
  tree->Branch("minOfMaxHF_41_46", &minOfMaxHF_41_46, "minOfMaxHF_41_46/F");
  tree->Branch("minOfMaxHF_41_47", &minOfMaxHF_41_47, "minOfMaxHF_41_47/F");
  tree->Branch("minOfMaxHF_41_48", &minOfMaxHF_41_48, "minOfMaxHF_41_48/F");
  tree->Branch("minOfMaxHF_41_49", &minOfMaxHF_41_49, "minOfMaxHF_41_49/F");
  tree->Branch("minOfMaxHF_41_50", &minOfMaxHF_41_50, "minOfMaxHF_41_50/F");
  tree->Branch("minOfMaxHF_41_51", &minOfMaxHF_41_51, "minOfMaxHF_41_51/F");
  tree->Branch("minOfMaxHF_41_52", &minOfMaxHF_41_52, "minOfMaxHF_41_52/F");

  tree->Branch("minOfMaxHF_42_43", &minOfMaxHF_42_43, "minOfMaxHF_42_43/F");
  tree->Branch("minOfMaxHF_42_44", &minOfMaxHF_42_44, "minOfMaxHF_42_44/F");
  tree->Branch("minOfMaxHF_42_45", &minOfMaxHF_42_45, "minOfMaxHF_42_45/F");
  tree->Branch("minOfMaxHF_42_46", &minOfMaxHF_42_46, "minOfMaxHF_42_46/F");
  tree->Branch("minOfMaxHF_42_47", &minOfMaxHF_42_47, "minOfMaxHF_42_47/F");
  tree->Branch("minOfMaxHF_42_48", &minOfMaxHF_42_48, "minOfMaxHF_42_48/F");
  tree->Branch("minOfMaxHF_42_49", &minOfMaxHF_42_49, "minOfMaxHF_42_49/F");
  tree->Branch("minOfMaxHF_42_50", &minOfMaxHF_42_50, "minOfMaxHF_42_50/F");
  tree->Branch("minOfMaxHF_42_51", &minOfMaxHF_42_51, "minOfMaxHF_42_51/F");
  tree->Branch("minOfMaxHF_42_52", &minOfMaxHF_42_52, "minOfMaxHF_42_52/F");

  tree->Branch("minOfMaxHF_43_44", &minOfMaxHF_43_44, "minOfMaxHF_43_44/F");
  tree->Branch("minOfMaxHF_43_45", &minOfMaxHF_43_45, "minOfMaxHF_43_45/F");
  tree->Branch("minOfMaxHF_43_46", &minOfMaxHF_43_46, "minOfMaxHF_43_46/F");
  tree->Branch("minOfMaxHF_43_47", &minOfMaxHF_43_47, "minOfMaxHF_43_47/F");
  tree->Branch("minOfMaxHF_43_48", &minOfMaxHF_43_48, "minOfMaxHF_43_48/F");
  tree->Branch("minOfMaxHF_43_49", &minOfMaxHF_43_49, "minOfMaxHF_43_49/F");
  tree->Branch("minOfMaxHF_43_50", &minOfMaxHF_43_50, "minOfMaxHF_43_50/F");
  tree->Branch("minOfMaxHF_43_51", &minOfMaxHF_43_51, "minOfMaxHF_43_51/F");
  tree->Branch("minOfMaxHF_43_52", &minOfMaxHF_43_52, "minOfMaxHF_43_52/F");

  tree->Branch("minOfMaxHF_44_45", &minOfMaxHF_44_45, "minOfMaxHF_44_45/F");
  tree->Branch("minOfMaxHF_44_46", &minOfMaxHF_44_46, "minOfMaxHF_44_46/F");
  tree->Branch("minOfMaxHF_44_47", &minOfMaxHF_44_47, "minOfMaxHF_44_47/F");
  tree->Branch("minOfMaxHF_44_48", &minOfMaxHF_44_48, "minOfMaxHF_44_48/F");
  tree->Branch("minOfMaxHF_44_49", &minOfMaxHF_44_49, "minOfMaxHF_44_49/F");
  tree->Branch("minOfMaxHF_44_50", &minOfMaxHF_44_50, "minOfMaxHF_44_50/F");
  tree->Branch("minOfMaxHF_44_51", &minOfMaxHF_44_51, "minOfMaxHF_44_51/F");
  tree->Branch("minOfMaxHF_44_52", &minOfMaxHF_44_52, "minOfMaxHF_44_52/F");

  tree->Branch("minOfMaxHF_45_46", &minOfMaxHF_45_46, "minOfMaxHF_45_46/F");
  tree->Branch("minOfMaxHF_45_47", &minOfMaxHF_45_47, "minOfMaxHF_45_47/F");
  tree->Branch("minOfMaxHF_45_48", &minOfMaxHF_45_48, "minOfMaxHF_45_48/F");
  tree->Branch("minOfMaxHF_45_49", &minOfMaxHF_45_49, "minOfMaxHF_45_49/F");
  tree->Branch("minOfMaxHF_45_50", &minOfMaxHF_45_50, "minOfMaxHF_45_50/F");
  tree->Branch("minOfMaxHF_45_51", &minOfMaxHF_45_51, "minOfMaxHF_45_51/F");
  tree->Branch("minOfMaxHF_45_52", &minOfMaxHF_45_52, "minOfMaxHF_45_52/F");

  tree->Branch("minOfMaxHF_46_47", &minOfMaxHF_46_47, "minOfMaxHF_46_47/F");
  tree->Branch("minOfMaxHF_46_48", &minOfMaxHF_46_48, "minOfMaxHF_46_48/F");
  tree->Branch("minOfMaxHF_46_49", &minOfMaxHF_46_49, "minOfMaxHF_46_49/F");
  tree->Branch("minOfMaxHF_46_50", &minOfMaxHF_46_50, "minOfMaxHF_46_50/F");
  tree->Branch("minOfMaxHF_46_51", &minOfMaxHF_46_51, "minOfMaxHF_46_51/F");
  tree->Branch("minOfMaxHF_46_52", &minOfMaxHF_46_52, "minOfMaxHF_46_52/F");

  tree->Branch("minOfMaxHF_47_48", &minOfMaxHF_47_48, "minOfMaxHF_47_48/F");
  tree->Branch("minOfMaxHF_47_49", &minOfMaxHF_47_49, "minOfMaxHF_47_49/F");
  tree->Branch("minOfMaxHF_47_50", &minOfMaxHF_47_50, "minOfMaxHF_47_50/F");
  tree->Branch("minOfMaxHF_47_51", &minOfMaxHF_47_51, "minOfMaxHF_47_51/F");
  tree->Branch("minOfMaxHF_47_52", &minOfMaxHF_47_52, "minOfMaxHF_47_52/F");

  tree->Branch("minOfMaxHF_48_49", &minOfMaxHF_48_49, "minOfMaxHF_48_49/F");
  tree->Branch("minOfMaxHF_48_50", &minOfMaxHF_48_50, "minOfMaxHF_48_50/F");
  tree->Branch("minOfMaxHF_48_51", &minOfMaxHF_48_51, "minOfMaxHF_48_51/F");
  tree->Branch("minOfMaxHF_48_52", &minOfMaxHF_48_52, "minOfMaxHF_48_52/F");

  tree->Branch("minOfMaxHF_49_50", &minOfMaxHF_49_50, "minOfMaxHF_49_50/F");
  tree->Branch("minOfMaxHF_49_51", &minOfMaxHF_49_51, "minOfMaxHF_49_51/F");
  tree->Branch("minOfMaxHF_49_52", &minOfMaxHF_49_52, "minOfMaxHF_49_52/F");

  tree->Branch("minOfMaxHF_50_51", &minOfMaxHF_50_51, "minOfMaxHF_50_51/F");
  tree->Branch("minOfMaxHF_50_52", &minOfMaxHF_50_52, "minOfMaxHF_50_52/F");

  tree->Branch("minOfMaxHF_51_52", &minOfMaxHF_51_52, "minOfMaxHF_51_52/F");
  // minofhfmax --end

  // maxhfenergyminus
  tree->Branch("maxHFEnergyMinus_30_31", &maxHFEnergyMinus_30_31, "maxHFEnergyMinus_30_31/F");
  tree->Branch("maxHFEnergyMinus_30_32", &maxHFEnergyMinus_30_32, "maxHFEnergyMinus_30_32/F");
  tree->Branch("maxHFEnergyMinus_30_33", &maxHFEnergyMinus_30_33, "maxHFEnergyMinus_30_33/F");
  tree->Branch("maxHFEnergyMinus_30_34", &maxHFEnergyMinus_30_34, "maxHFEnergyMinus_30_34/F");
  tree->Branch("maxHFEnergyMinus_30_35", &maxHFEnergyMinus_30_35, "maxHFEnergyMinus_30_35/F");
  tree->Branch("maxHFEnergyMinus_30_36", &maxHFEnergyMinus_30_36, "maxHFEnergyMinus_30_36/F");
  tree->Branch("maxHFEnergyMinus_30_37", &maxHFEnergyMinus_30_37, "maxHFEnergyMinus_30_37/F");
  tree->Branch("maxHFEnergyMinus_30_38", &maxHFEnergyMinus_30_38, "maxHFEnergyMinus_30_38/F");
  tree->Branch("maxHFEnergyMinus_30_39", &maxHFEnergyMinus_30_39, "maxHFEnergyMinus_30_39/F");
  tree->Branch("maxHFEnergyMinus_30_40", &maxHFEnergyMinus_30_40, "maxHFEnergyMinus_30_40/F");
  tree->Branch("maxHFEnergyMinus_30_41", &maxHFEnergyMinus_30_41, "maxHFEnergyMinus_30_41/F");
  tree->Branch("maxHFEnergyMinus_30_42", &maxHFEnergyMinus_30_42, "maxHFEnergyMinus_30_42/F");
  tree->Branch("maxHFEnergyMinus_30_43", &maxHFEnergyMinus_30_43, "maxHFEnergyMinus_30_43/F");
  tree->Branch("maxHFEnergyMinus_30_44", &maxHFEnergyMinus_30_44, "maxHFEnergyMinus_30_44/F");
  tree->Branch("maxHFEnergyMinus_30_45", &maxHFEnergyMinus_30_45, "maxHFEnergyMinus_30_45/F");
  tree->Branch("maxHFEnergyMinus_30_46", &maxHFEnergyMinus_30_46, "maxHFEnergyMinus_30_46/F");
  tree->Branch("maxHFEnergyMinus_30_47", &maxHFEnergyMinus_30_47, "maxHFEnergyMinus_30_47/F");
  tree->Branch("maxHFEnergyMinus_30_48", &maxHFEnergyMinus_30_48, "maxHFEnergyMinus_30_48/F");
  tree->Branch("maxHFEnergyMinus_30_49", &maxHFEnergyMinus_30_49, "maxHFEnergyMinus_30_49/F");
  tree->Branch("maxHFEnergyMinus_30_50", &maxHFEnergyMinus_30_50, "maxHFEnergyMinus_30_50/F");
  tree->Branch("maxHFEnergyMinus_30_51", &maxHFEnergyMinus_30_51, "maxHFEnergyMinus_30_51/F");
  tree->Branch("maxHFEnergyMinus_30_52", &maxHFEnergyMinus_30_52, "maxHFEnergyMinus_30_52/F");

  tree->Branch("maxHFEnergyMinus_31_32", &maxHFEnergyMinus_31_32, "maxHFEnergyMinus_31_32/F");
  tree->Branch("maxHFEnergyMinus_31_33", &maxHFEnergyMinus_31_33, "maxHFEnergyMinus_31_33/F");
  tree->Branch("maxHFEnergyMinus_31_34", &maxHFEnergyMinus_31_34, "maxHFEnergyMinus_31_34/F");
  tree->Branch("maxHFEnergyMinus_31_35", &maxHFEnergyMinus_31_35, "maxHFEnergyMinus_31_35/F");
  tree->Branch("maxHFEnergyMinus_31_36", &maxHFEnergyMinus_31_36, "maxHFEnergyMinus_31_36/F");
  tree->Branch("maxHFEnergyMinus_31_37", &maxHFEnergyMinus_31_37, "maxHFEnergyMinus_31_37/F");
  tree->Branch("maxHFEnergyMinus_31_38", &maxHFEnergyMinus_31_38, "maxHFEnergyMinus_31_38/F");
  tree->Branch("maxHFEnergyMinus_31_39", &maxHFEnergyMinus_31_39, "maxHFEnergyMinus_31_39/F");
  tree->Branch("maxHFEnergyMinus_31_40", &maxHFEnergyMinus_31_40, "maxHFEnergyMinus_31_40/F");
  tree->Branch("maxHFEnergyMinus_31_41", &maxHFEnergyMinus_31_41, "maxHFEnergyMinus_31_41/F");
  tree->Branch("maxHFEnergyMinus_31_42", &maxHFEnergyMinus_31_42, "maxHFEnergyMinus_31_42/F");
  tree->Branch("maxHFEnergyMinus_31_43", &maxHFEnergyMinus_31_43, "maxHFEnergyMinus_31_43/F");
  tree->Branch("maxHFEnergyMinus_31_44", &maxHFEnergyMinus_31_44, "maxHFEnergyMinus_31_44/F");
  tree->Branch("maxHFEnergyMinus_31_45", &maxHFEnergyMinus_31_45, "maxHFEnergyMinus_31_45/F");
  tree->Branch("maxHFEnergyMinus_31_46", &maxHFEnergyMinus_31_46, "maxHFEnergyMinus_31_46/F");
  tree->Branch("maxHFEnergyMinus_31_47", &maxHFEnergyMinus_31_47, "maxHFEnergyMinus_31_47/F");
  tree->Branch("maxHFEnergyMinus_31_48", &maxHFEnergyMinus_31_48, "maxHFEnergyMinus_31_48/F");
  tree->Branch("maxHFEnergyMinus_31_49", &maxHFEnergyMinus_31_49, "maxHFEnergyMinus_31_49/F");
  tree->Branch("maxHFEnergyMinus_31_50", &maxHFEnergyMinus_31_50, "maxHFEnergyMinus_31_50/F");
  tree->Branch("maxHFEnergyMinus_31_51", &maxHFEnergyMinus_31_51, "maxHFEnergyMinus_31_51/F");
  tree->Branch("maxHFEnergyMinus_31_52", &maxHFEnergyMinus_31_52, "maxHFEnergyMinus_31_52/F");

  tree->Branch("maxHFEnergyMinus_32_33", &maxHFEnergyMinus_32_33, "maxHFEnergyMinus_32_33/F");
  tree->Branch("maxHFEnergyMinus_32_34", &maxHFEnergyMinus_32_34, "maxHFEnergyMinus_32_34/F");
  tree->Branch("maxHFEnergyMinus_32_35", &maxHFEnergyMinus_32_35, "maxHFEnergyMinus_32_35/F");
  tree->Branch("maxHFEnergyMinus_32_36", &maxHFEnergyMinus_32_36, "maxHFEnergyMinus_32_36/F");
  tree->Branch("maxHFEnergyMinus_32_37", &maxHFEnergyMinus_32_37, "maxHFEnergyMinus_32_37/F");
  tree->Branch("maxHFEnergyMinus_32_38", &maxHFEnergyMinus_32_38, "maxHFEnergyMinus_32_38/F");
  tree->Branch("maxHFEnergyMinus_32_39", &maxHFEnergyMinus_32_39, "maxHFEnergyMinus_32_39/F");
  tree->Branch("maxHFEnergyMinus_32_40", &maxHFEnergyMinus_32_40, "maxHFEnergyMinus_32_40/F");
  tree->Branch("maxHFEnergyMinus_32_41", &maxHFEnergyMinus_32_41, "maxHFEnergyMinus_32_41/F");
  tree->Branch("maxHFEnergyMinus_32_42", &maxHFEnergyMinus_32_42, "maxHFEnergyMinus_32_42/F");
  tree->Branch("maxHFEnergyMinus_32_43", &maxHFEnergyMinus_32_43, "maxHFEnergyMinus_32_43/F");
  tree->Branch("maxHFEnergyMinus_32_44", &maxHFEnergyMinus_32_44, "maxHFEnergyMinus_32_44/F");
  tree->Branch("maxHFEnergyMinus_32_45", &maxHFEnergyMinus_32_45, "maxHFEnergyMinus_32_45/F");
  tree->Branch("maxHFEnergyMinus_32_46", &maxHFEnergyMinus_32_46, "maxHFEnergyMinus_32_46/F");
  tree->Branch("maxHFEnergyMinus_32_47", &maxHFEnergyMinus_32_47, "maxHFEnergyMinus_32_47/F");
  tree->Branch("maxHFEnergyMinus_32_48", &maxHFEnergyMinus_32_48, "maxHFEnergyMinus_32_48/F");
  tree->Branch("maxHFEnergyMinus_32_49", &maxHFEnergyMinus_32_49, "maxHFEnergyMinus_32_49/F");
  tree->Branch("maxHFEnergyMinus_32_50", &maxHFEnergyMinus_32_50, "maxHFEnergyMinus_32_50/F");
  tree->Branch("maxHFEnergyMinus_32_51", &maxHFEnergyMinus_32_51, "maxHFEnergyMinus_32_51/F");
  tree->Branch("maxHFEnergyMinus_32_52", &maxHFEnergyMinus_32_52, "maxHFEnergyMinus_32_52/F");

  tree->Branch("maxHFEnergyMinus_33_34", &maxHFEnergyMinus_33_34, "maxHFEnergyMinus_33_34/F");
  tree->Branch("maxHFEnergyMinus_33_35", &maxHFEnergyMinus_33_35, "maxHFEnergyMinus_33_35/F");
  tree->Branch("maxHFEnergyMinus_33_36", &maxHFEnergyMinus_33_36, "maxHFEnergyMinus_33_36/F");
  tree->Branch("maxHFEnergyMinus_33_37", &maxHFEnergyMinus_33_37, "maxHFEnergyMinus_33_37/F");
  tree->Branch("maxHFEnergyMinus_33_38", &maxHFEnergyMinus_33_38, "maxHFEnergyMinus_33_38/F");
  tree->Branch("maxHFEnergyMinus_33_39", &maxHFEnergyMinus_33_39, "maxHFEnergyMinus_33_39/F");
  tree->Branch("maxHFEnergyMinus_33_40", &maxHFEnergyMinus_33_40, "maxHFEnergyMinus_33_40/F");
  tree->Branch("maxHFEnergyMinus_33_41", &maxHFEnergyMinus_33_41, "maxHFEnergyMinus_33_41/F");
  tree->Branch("maxHFEnergyMinus_33_42", &maxHFEnergyMinus_33_42, "maxHFEnergyMinus_33_42/F");
  tree->Branch("maxHFEnergyMinus_33_43", &maxHFEnergyMinus_33_43, "maxHFEnergyMinus_33_43/F");
  tree->Branch("maxHFEnergyMinus_33_44", &maxHFEnergyMinus_33_44, "maxHFEnergyMinus_33_44/F");
  tree->Branch("maxHFEnergyMinus_33_45", &maxHFEnergyMinus_33_45, "maxHFEnergyMinus_33_45/F");
  tree->Branch("maxHFEnergyMinus_33_46", &maxHFEnergyMinus_33_46, "maxHFEnergyMinus_33_46/F");
  tree->Branch("maxHFEnergyMinus_33_47", &maxHFEnergyMinus_33_47, "maxHFEnergyMinus_33_47/F");
  tree->Branch("maxHFEnergyMinus_33_48", &maxHFEnergyMinus_33_48, "maxHFEnergyMinus_33_48/F");
  tree->Branch("maxHFEnergyMinus_33_49", &maxHFEnergyMinus_33_49, "maxHFEnergyMinus_33_49/F");
  tree->Branch("maxHFEnergyMinus_33_50", &maxHFEnergyMinus_33_50, "maxHFEnergyMinus_33_50/F");
  tree->Branch("maxHFEnergyMinus_33_51", &maxHFEnergyMinus_33_51, "maxHFEnergyMinus_33_51/F");
  tree->Branch("maxHFEnergyMinus_33_52", &maxHFEnergyMinus_33_52, "maxHFEnergyMinus_33_52/F");

  tree->Branch("maxHFEnergyMinus_34_35", &maxHFEnergyMinus_34_35, "maxHFEnergyMinus_34_35/F");
  tree->Branch("maxHFEnergyMinus_34_36", &maxHFEnergyMinus_34_36, "maxHFEnergyMinus_34_36/F");
  tree->Branch("maxHFEnergyMinus_34_37", &maxHFEnergyMinus_34_37, "maxHFEnergyMinus_34_37/F");
  tree->Branch("maxHFEnergyMinus_34_38", &maxHFEnergyMinus_34_38, "maxHFEnergyMinus_34_38/F");
  tree->Branch("maxHFEnergyMinus_34_39", &maxHFEnergyMinus_34_39, "maxHFEnergyMinus_34_39/F");
  tree->Branch("maxHFEnergyMinus_34_40", &maxHFEnergyMinus_34_40, "maxHFEnergyMinus_34_40/F");
  tree->Branch("maxHFEnergyMinus_34_41", &maxHFEnergyMinus_34_41, "maxHFEnergyMinus_34_41/F");
  tree->Branch("maxHFEnergyMinus_34_42", &maxHFEnergyMinus_34_42, "maxHFEnergyMinus_34_42/F");
  tree->Branch("maxHFEnergyMinus_34_43", &maxHFEnergyMinus_34_43, "maxHFEnergyMinus_34_43/F");
  tree->Branch("maxHFEnergyMinus_34_44", &maxHFEnergyMinus_34_44, "maxHFEnergyMinus_34_44/F");
  tree->Branch("maxHFEnergyMinus_34_45", &maxHFEnergyMinus_34_45, "maxHFEnergyMinus_34_45/F");
  tree->Branch("maxHFEnergyMinus_34_46", &maxHFEnergyMinus_34_46, "maxHFEnergyMinus_34_46/F");
  tree->Branch("maxHFEnergyMinus_34_47", &maxHFEnergyMinus_34_47, "maxHFEnergyMinus_34_47/F");
  tree->Branch("maxHFEnergyMinus_34_48", &maxHFEnergyMinus_34_48, "maxHFEnergyMinus_34_48/F");
  tree->Branch("maxHFEnergyMinus_34_49", &maxHFEnergyMinus_34_49, "maxHFEnergyMinus_34_49/F");
  tree->Branch("maxHFEnergyMinus_34_50", &maxHFEnergyMinus_34_50, "maxHFEnergyMinus_34_50/F");
  tree->Branch("maxHFEnergyMinus_34_51", &maxHFEnergyMinus_34_51, "maxHFEnergyMinus_34_51/F");
  tree->Branch("maxHFEnergyMinus_34_52", &maxHFEnergyMinus_34_52, "maxHFEnergyMinus_34_52/F");

  tree->Branch("maxHFEnergyMinus_35_36", &maxHFEnergyMinus_35_36, "maxHFEnergyMinus_35_36/F");
  tree->Branch("maxHFEnergyMinus_35_37", &maxHFEnergyMinus_35_37, "maxHFEnergyMinus_35_37/F");
  tree->Branch("maxHFEnergyMinus_35_38", &maxHFEnergyMinus_35_38, "maxHFEnergyMinus_35_38/F");
  tree->Branch("maxHFEnergyMinus_35_39", &maxHFEnergyMinus_35_39, "maxHFEnergyMinus_35_39/F");
  tree->Branch("maxHFEnergyMinus_35_40", &maxHFEnergyMinus_35_40, "maxHFEnergyMinus_35_40/F");
  tree->Branch("maxHFEnergyMinus_35_41", &maxHFEnergyMinus_35_41, "maxHFEnergyMinus_35_41/F");
  tree->Branch("maxHFEnergyMinus_35_42", &maxHFEnergyMinus_35_42, "maxHFEnergyMinus_35_42/F");
  tree->Branch("maxHFEnergyMinus_35_43", &maxHFEnergyMinus_35_43, "maxHFEnergyMinus_35_43/F");
  tree->Branch("maxHFEnergyMinus_35_44", &maxHFEnergyMinus_35_44, "maxHFEnergyMinus_35_44/F");
  tree->Branch("maxHFEnergyMinus_35_45", &maxHFEnergyMinus_35_45, "maxHFEnergyMinus_35_45/F");
  tree->Branch("maxHFEnergyMinus_35_46", &maxHFEnergyMinus_35_46, "maxHFEnergyMinus_35_46/F");
  tree->Branch("maxHFEnergyMinus_35_47", &maxHFEnergyMinus_35_47, "maxHFEnergyMinus_35_47/F");
  tree->Branch("maxHFEnergyMinus_35_48", &maxHFEnergyMinus_35_48, "maxHFEnergyMinus_35_48/F");
  tree->Branch("maxHFEnergyMinus_35_49", &maxHFEnergyMinus_35_49, "maxHFEnergyMinus_35_49/F");
  tree->Branch("maxHFEnergyMinus_35_50", &maxHFEnergyMinus_35_50, "maxHFEnergyMinus_35_50/F");
  tree->Branch("maxHFEnergyMinus_35_51", &maxHFEnergyMinus_35_51, "maxHFEnergyMinus_35_51/F");
  tree->Branch("maxHFEnergyMinus_35_52", &maxHFEnergyMinus_35_52, "maxHFEnergyMinus_35_52/F");

  tree->Branch("maxHFEnergyMinus_36_37", &maxHFEnergyMinus_36_37, "maxHFEnergyMinus_36_37/F");
  tree->Branch("maxHFEnergyMinus_36_38", &maxHFEnergyMinus_36_38, "maxHFEnergyMinus_36_38/F");
  tree->Branch("maxHFEnergyMinus_36_39", &maxHFEnergyMinus_36_39, "maxHFEnergyMinus_36_39/F");
  tree->Branch("maxHFEnergyMinus_36_40", &maxHFEnergyMinus_36_40, "maxHFEnergyMinus_36_40/F");
  tree->Branch("maxHFEnergyMinus_36_41", &maxHFEnergyMinus_36_41, "maxHFEnergyMinus_36_41/F");
  tree->Branch("maxHFEnergyMinus_36_42", &maxHFEnergyMinus_36_42, "maxHFEnergyMinus_36_42/F");
  tree->Branch("maxHFEnergyMinus_36_43", &maxHFEnergyMinus_36_43, "maxHFEnergyMinus_36_43/F");
  tree->Branch("maxHFEnergyMinus_36_44", &maxHFEnergyMinus_36_44, "maxHFEnergyMinus_36_44/F");
  tree->Branch("maxHFEnergyMinus_36_45", &maxHFEnergyMinus_36_45, "maxHFEnergyMinus_36_45/F");
  tree->Branch("maxHFEnergyMinus_36_46", &maxHFEnergyMinus_36_46, "maxHFEnergyMinus_36_46/F");
  tree->Branch("maxHFEnergyMinus_36_47", &maxHFEnergyMinus_36_47, "maxHFEnergyMinus_36_47/F");
  tree->Branch("maxHFEnergyMinus_36_48", &maxHFEnergyMinus_36_48, "maxHFEnergyMinus_36_48/F");
  tree->Branch("maxHFEnergyMinus_36_49", &maxHFEnergyMinus_36_49, "maxHFEnergyMinus_36_49/F");
  tree->Branch("maxHFEnergyMinus_36_50", &maxHFEnergyMinus_36_50, "maxHFEnergyMinus_36_50/F");
  tree->Branch("maxHFEnergyMinus_36_51", &maxHFEnergyMinus_36_51, "maxHFEnergyMinus_36_51/F");
  tree->Branch("maxHFEnergyMinus_36_52", &maxHFEnergyMinus_36_52, "maxHFEnergyMinus_36_52/F");

  tree->Branch("maxHFEnergyMinus_37_38", &maxHFEnergyMinus_37_38, "maxHFEnergyMinus_37_38/F");
  tree->Branch("maxHFEnergyMinus_37_39", &maxHFEnergyMinus_37_39, "maxHFEnergyMinus_37_39/F");
  tree->Branch("maxHFEnergyMinus_37_40", &maxHFEnergyMinus_37_40, "maxHFEnergyMinus_37_40/F");
  tree->Branch("maxHFEnergyMinus_37_41", &maxHFEnergyMinus_37_41, "maxHFEnergyMinus_37_41/F");
  tree->Branch("maxHFEnergyMinus_37_42", &maxHFEnergyMinus_37_42, "maxHFEnergyMinus_37_42/F");
  tree->Branch("maxHFEnergyMinus_37_43", &maxHFEnergyMinus_37_43, "maxHFEnergyMinus_37_43/F");
  tree->Branch("maxHFEnergyMinus_37_44", &maxHFEnergyMinus_37_44, "maxHFEnergyMinus_37_44/F");
  tree->Branch("maxHFEnergyMinus_37_45", &maxHFEnergyMinus_37_45, "maxHFEnergyMinus_37_45/F");
  tree->Branch("maxHFEnergyMinus_37_46", &maxHFEnergyMinus_37_46, "maxHFEnergyMinus_37_46/F");
  tree->Branch("maxHFEnergyMinus_37_47", &maxHFEnergyMinus_37_47, "maxHFEnergyMinus_37_47/F");
  tree->Branch("maxHFEnergyMinus_37_48", &maxHFEnergyMinus_37_48, "maxHFEnergyMinus_37_48/F");
  tree->Branch("maxHFEnergyMinus_37_49", &maxHFEnergyMinus_37_49, "maxHFEnergyMinus_37_49/F");
  tree->Branch("maxHFEnergyMinus_37_50", &maxHFEnergyMinus_37_50, "maxHFEnergyMinus_37_50/F");
  tree->Branch("maxHFEnergyMinus_37_51", &maxHFEnergyMinus_37_51, "maxHFEnergyMinus_37_51/F");
  tree->Branch("maxHFEnergyMinus_37_52", &maxHFEnergyMinus_37_52, "maxHFEnergyMinus_37_52/F");

  tree->Branch("maxHFEnergyMinus_38_39", &maxHFEnergyMinus_38_39, "maxHFEnergyMinus_38_39/F");
  tree->Branch("maxHFEnergyMinus_38_40", &maxHFEnergyMinus_38_40, "maxHFEnergyMinus_38_40/F");
  tree->Branch("maxHFEnergyMinus_38_41", &maxHFEnergyMinus_38_41, "maxHFEnergyMinus_38_41/F");
  tree->Branch("maxHFEnergyMinus_38_42", &maxHFEnergyMinus_38_42, "maxHFEnergyMinus_38_42/F");
  tree->Branch("maxHFEnergyMinus_38_43", &maxHFEnergyMinus_38_43, "maxHFEnergyMinus_38_43/F");
  tree->Branch("maxHFEnergyMinus_38_44", &maxHFEnergyMinus_38_44, "maxHFEnergyMinus_38_44/F");
  tree->Branch("maxHFEnergyMinus_38_45", &maxHFEnergyMinus_38_45, "maxHFEnergyMinus_38_45/F");
  tree->Branch("maxHFEnergyMinus_38_46", &maxHFEnergyMinus_38_46, "maxHFEnergyMinus_38_46/F");
  tree->Branch("maxHFEnergyMinus_38_47", &maxHFEnergyMinus_38_47, "maxHFEnergyMinus_38_47/F");
  tree->Branch("maxHFEnergyMinus_38_48", &maxHFEnergyMinus_38_48, "maxHFEnergyMinus_38_48/F");
  tree->Branch("maxHFEnergyMinus_38_49", &maxHFEnergyMinus_38_49, "maxHFEnergyMinus_38_49/F");
  tree->Branch("maxHFEnergyMinus_38_50", &maxHFEnergyMinus_38_50, "maxHFEnergyMinus_38_50/F");
  tree->Branch("maxHFEnergyMinus_38_51", &maxHFEnergyMinus_38_51, "maxHFEnergyMinus_38_51/F");
  tree->Branch("maxHFEnergyMinus_38_52", &maxHFEnergyMinus_38_52, "maxHFEnergyMinus_38_52/F");

  tree->Branch("maxHFEnergyMinus_39_40", &maxHFEnergyMinus_39_40, "maxHFEnergyMinus_39_40/F");
  tree->Branch("maxHFEnergyMinus_39_41", &maxHFEnergyMinus_39_41, "maxHFEnergyMinus_39_41/F");
  tree->Branch("maxHFEnergyMinus_39_42", &maxHFEnergyMinus_39_42, "maxHFEnergyMinus_39_42/F");
  tree->Branch("maxHFEnergyMinus_39_43", &maxHFEnergyMinus_39_43, "maxHFEnergyMinus_39_43/F");
  tree->Branch("maxHFEnergyMinus_39_44", &maxHFEnergyMinus_39_44, "maxHFEnergyMinus_39_44/F");
  tree->Branch("maxHFEnergyMinus_39_45", &maxHFEnergyMinus_39_45, "maxHFEnergyMinus_39_45/F");
  tree->Branch("maxHFEnergyMinus_39_46", &maxHFEnergyMinus_39_46, "maxHFEnergyMinus_39_46/F");
  tree->Branch("maxHFEnergyMinus_39_47", &maxHFEnergyMinus_39_47, "maxHFEnergyMinus_39_47/F");
  tree->Branch("maxHFEnergyMinus_39_48", &maxHFEnergyMinus_39_48, "maxHFEnergyMinus_39_48/F");
  tree->Branch("maxHFEnergyMinus_39_49", &maxHFEnergyMinus_39_49, "maxHFEnergyMinus_39_49/F");
  tree->Branch("maxHFEnergyMinus_39_50", &maxHFEnergyMinus_39_50, "maxHFEnergyMinus_39_50/F");
  tree->Branch("maxHFEnergyMinus_39_51", &maxHFEnergyMinus_39_51, "maxHFEnergyMinus_39_51/F");
  tree->Branch("maxHFEnergyMinus_39_52", &maxHFEnergyMinus_39_52, "maxHFEnergyMinus_39_52/F");

  tree->Branch("maxHFEnergyMinus_40_41", &maxHFEnergyMinus_40_41, "maxHFEnergyMinus_40_41/F");
  tree->Branch("maxHFEnergyMinus_40_42", &maxHFEnergyMinus_40_42, "maxHFEnergyMinus_40_42/F");
  tree->Branch("maxHFEnergyMinus_40_43", &maxHFEnergyMinus_40_43, "maxHFEnergyMinus_40_43/F");
  tree->Branch("maxHFEnergyMinus_40_44", &maxHFEnergyMinus_40_44, "maxHFEnergyMinus_40_44/F");
  tree->Branch("maxHFEnergyMinus_40_45", &maxHFEnergyMinus_40_45, "maxHFEnergyMinus_40_45/F");
  tree->Branch("maxHFEnergyMinus_40_46", &maxHFEnergyMinus_40_46, "maxHFEnergyMinus_40_46/F");
  tree->Branch("maxHFEnergyMinus_40_47", &maxHFEnergyMinus_40_47, "maxHFEnergyMinus_40_47/F");
  tree->Branch("maxHFEnergyMinus_40_48", &maxHFEnergyMinus_40_48, "maxHFEnergyMinus_40_48/F");
  tree->Branch("maxHFEnergyMinus_40_49", &maxHFEnergyMinus_40_49, "maxHFEnergyMinus_40_49/F");
  tree->Branch("maxHFEnergyMinus_40_50", &maxHFEnergyMinus_40_50, "maxHFEnergyMinus_40_50/F");
  tree->Branch("maxHFEnergyMinus_40_51", &maxHFEnergyMinus_40_51, "maxHFEnergyMinus_40_51/F");
  tree->Branch("maxHFEnergyMinus_40_52", &maxHFEnergyMinus_40_52, "maxHFEnergyMinus_40_52/F");

  tree->Branch("maxHFEnergyMinus_41_42", &maxHFEnergyMinus_41_42, "maxHFEnergyMinus_41_42/F");
  tree->Branch("maxHFEnergyMinus_41_43", &maxHFEnergyMinus_41_43, "maxHFEnergyMinus_41_43/F");
  tree->Branch("maxHFEnergyMinus_41_44", &maxHFEnergyMinus_41_44, "maxHFEnergyMinus_41_44/F");
  tree->Branch("maxHFEnergyMinus_41_45", &maxHFEnergyMinus_41_45, "maxHFEnergyMinus_41_45/F");
  tree->Branch("maxHFEnergyMinus_41_46", &maxHFEnergyMinus_41_46, "maxHFEnergyMinus_41_46/F");
  tree->Branch("maxHFEnergyMinus_41_47", &maxHFEnergyMinus_41_47, "maxHFEnergyMinus_41_47/F");
  tree->Branch("maxHFEnergyMinus_41_48", &maxHFEnergyMinus_41_48, "maxHFEnergyMinus_41_48/F");
  tree->Branch("maxHFEnergyMinus_41_49", &maxHFEnergyMinus_41_49, "maxHFEnergyMinus_41_49/F");
  tree->Branch("maxHFEnergyMinus_41_50", &maxHFEnergyMinus_41_50, "maxHFEnergyMinus_41_50/F");
  tree->Branch("maxHFEnergyMinus_41_51", &maxHFEnergyMinus_41_51, "maxHFEnergyMinus_41_51/F");
  tree->Branch("maxHFEnergyMinus_41_52", &maxHFEnergyMinus_41_52, "maxHFEnergyMinus_41_52/F");

  tree->Branch("maxHFEnergyMinus_42_43", &maxHFEnergyMinus_42_43, "maxHFEnergyMinus_42_43/F");
  tree->Branch("maxHFEnergyMinus_42_44", &maxHFEnergyMinus_42_44, "maxHFEnergyMinus_42_44/F");
  tree->Branch("maxHFEnergyMinus_42_45", &maxHFEnergyMinus_42_45, "maxHFEnergyMinus_42_45/F");
  tree->Branch("maxHFEnergyMinus_42_46", &maxHFEnergyMinus_42_46, "maxHFEnergyMinus_42_46/F");
  tree->Branch("maxHFEnergyMinus_42_47", &maxHFEnergyMinus_42_47, "maxHFEnergyMinus_42_47/F");
  tree->Branch("maxHFEnergyMinus_42_48", &maxHFEnergyMinus_42_48, "maxHFEnergyMinus_42_48/F");
  tree->Branch("maxHFEnergyMinus_42_49", &maxHFEnergyMinus_42_49, "maxHFEnergyMinus_42_49/F");
  tree->Branch("maxHFEnergyMinus_42_50", &maxHFEnergyMinus_42_50, "maxHFEnergyMinus_42_50/F");
  tree->Branch("maxHFEnergyMinus_42_51", &maxHFEnergyMinus_42_51, "maxHFEnergyMinus_42_51/F");
  tree->Branch("maxHFEnergyMinus_42_52", &maxHFEnergyMinus_42_52, "maxHFEnergyMinus_42_52/F");

  tree->Branch("maxHFEnergyMinus_43_44", &maxHFEnergyMinus_43_44, "maxHFEnergyMinus_43_44/F");
  tree->Branch("maxHFEnergyMinus_43_45", &maxHFEnergyMinus_43_45, "maxHFEnergyMinus_43_45/F");
  tree->Branch("maxHFEnergyMinus_43_46", &maxHFEnergyMinus_43_46, "maxHFEnergyMinus_43_46/F");
  tree->Branch("maxHFEnergyMinus_43_47", &maxHFEnergyMinus_43_47, "maxHFEnergyMinus_43_47/F");
  tree->Branch("maxHFEnergyMinus_43_48", &maxHFEnergyMinus_43_48, "maxHFEnergyMinus_43_48/F");
  tree->Branch("maxHFEnergyMinus_43_49", &maxHFEnergyMinus_43_49, "maxHFEnergyMinus_43_49/F");
  tree->Branch("maxHFEnergyMinus_43_50", &maxHFEnergyMinus_43_50, "maxHFEnergyMinus_43_50/F");
  tree->Branch("maxHFEnergyMinus_43_51", &maxHFEnergyMinus_43_51, "maxHFEnergyMinus_43_51/F");
  tree->Branch("maxHFEnergyMinus_43_52", &maxHFEnergyMinus_43_52, "maxHFEnergyMinus_43_52/F");

  tree->Branch("maxHFEnergyMinus_44_45", &maxHFEnergyMinus_44_45, "maxHFEnergyMinus_44_45/F");
  tree->Branch("maxHFEnergyMinus_44_46", &maxHFEnergyMinus_44_46, "maxHFEnergyMinus_44_46/F");
  tree->Branch("maxHFEnergyMinus_44_47", &maxHFEnergyMinus_44_47, "maxHFEnergyMinus_44_47/F");
  tree->Branch("maxHFEnergyMinus_44_48", &maxHFEnergyMinus_44_48, "maxHFEnergyMinus_44_48/F");
  tree->Branch("maxHFEnergyMinus_44_49", &maxHFEnergyMinus_44_49, "maxHFEnergyMinus_44_49/F");
  tree->Branch("maxHFEnergyMinus_44_50", &maxHFEnergyMinus_44_50, "maxHFEnergyMinus_44_50/F");
  tree->Branch("maxHFEnergyMinus_44_51", &maxHFEnergyMinus_44_51, "maxHFEnergyMinus_44_51/F");
  tree->Branch("maxHFEnergyMinus_44_52", &maxHFEnergyMinus_44_52, "maxHFEnergyMinus_44_52/F");

  tree->Branch("maxHFEnergyMinus_45_46", &maxHFEnergyMinus_45_46, "maxHFEnergyMinus_45_46/F");
  tree->Branch("maxHFEnergyMinus_45_47", &maxHFEnergyMinus_45_47, "maxHFEnergyMinus_45_47/F");
  tree->Branch("maxHFEnergyMinus_45_48", &maxHFEnergyMinus_45_48, "maxHFEnergyMinus_45_48/F");
  tree->Branch("maxHFEnergyMinus_45_49", &maxHFEnergyMinus_45_49, "maxHFEnergyMinus_45_49/F");
  tree->Branch("maxHFEnergyMinus_45_50", &maxHFEnergyMinus_45_50, "maxHFEnergyMinus_45_50/F");
  tree->Branch("maxHFEnergyMinus_45_51", &maxHFEnergyMinus_45_51, "maxHFEnergyMinus_45_51/F");
  tree->Branch("maxHFEnergyMinus_45_52", &maxHFEnergyMinus_45_52, "maxHFEnergyMinus_45_52/F");

  tree->Branch("maxHFEnergyMinus_46_47", &maxHFEnergyMinus_46_47, "maxHFEnergyMinus_46_47/F");
  tree->Branch("maxHFEnergyMinus_46_48", &maxHFEnergyMinus_46_48, "maxHFEnergyMinus_46_48/F");
  tree->Branch("maxHFEnergyMinus_46_49", &maxHFEnergyMinus_46_49, "maxHFEnergyMinus_46_49/F");
  tree->Branch("maxHFEnergyMinus_46_50", &maxHFEnergyMinus_46_50, "maxHFEnergyMinus_46_50/F");
  tree->Branch("maxHFEnergyMinus_46_51", &maxHFEnergyMinus_46_51, "maxHFEnergyMinus_46_51/F");
  tree->Branch("maxHFEnergyMinus_46_52", &maxHFEnergyMinus_46_52, "maxHFEnergyMinus_46_52/F");

  tree->Branch("maxHFEnergyMinus_47_48", &maxHFEnergyMinus_47_48, "maxHFEnergyMinus_47_48/F");
  tree->Branch("maxHFEnergyMinus_47_49", &maxHFEnergyMinus_47_49, "maxHFEnergyMinus_47_49/F");
  tree->Branch("maxHFEnergyMinus_47_50", &maxHFEnergyMinus_47_50, "maxHFEnergyMinus_47_50/F");
  tree->Branch("maxHFEnergyMinus_47_51", &maxHFEnergyMinus_47_51, "maxHFEnergyMinus_47_51/F");
  tree->Branch("maxHFEnergyMinus_47_52", &maxHFEnergyMinus_47_52, "maxHFEnergyMinus_47_52/F");

  tree->Branch("maxHFEnergyMinus_48_49", &maxHFEnergyMinus_48_49, "maxHFEnergyMinus_48_49/F");
  tree->Branch("maxHFEnergyMinus_48_50", &maxHFEnergyMinus_48_50, "maxHFEnergyMinus_48_50/F");
  tree->Branch("maxHFEnergyMinus_48_51", &maxHFEnergyMinus_48_51, "maxHFEnergyMinus_48_51/F");
  tree->Branch("maxHFEnergyMinus_48_52", &maxHFEnergyMinus_48_52, "maxHFEnergyMinus_48_52/F");

  tree->Branch("maxHFEnergyMinus_49_50", &maxHFEnergyMinus_49_50, "maxHFEnergyMinus_49_50/F");
  tree->Branch("maxHFEnergyMinus_49_51", &maxHFEnergyMinus_49_51, "maxHFEnergyMinus_49_51/F");
  tree->Branch("maxHFEnergyMinus_49_52", &maxHFEnergyMinus_49_52, "maxHFEnergyMinus_49_52/F");

  tree->Branch("maxHFEnergyMinus_50_51", &maxHFEnergyMinus_50_51, "maxHFEnergyMinus_50_51/F");
  tree->Branch("maxHFEnergyMinus_50_52", &maxHFEnergyMinus_50_52, "maxHFEnergyMinus_50_52/F");

  tree->Branch("maxHFEnergyMinus_51_52", &maxHFEnergyMinus_51_52, "maxHFEnergyMinus_51_52/F");
  // maxhfenergyminus --end

  // maxhfenergyplus
  tree->Branch("maxHFEnergyPlus_30_31", &maxHFEnergyPlus_30_31, "maxHFEnergyPlus_30_31/F");
  tree->Branch("maxHFEnergyPlus_30_32", &maxHFEnergyPlus_30_32, "maxHFEnergyPlus_30_32/F");
  tree->Branch("maxHFEnergyPlus_30_33", &maxHFEnergyPlus_30_33, "maxHFEnergyPlus_30_33/F");
  tree->Branch("maxHFEnergyPlus_30_34", &maxHFEnergyPlus_30_34, "maxHFEnergyPlus_30_34/F");
  tree->Branch("maxHFEnergyPlus_30_35", &maxHFEnergyPlus_30_35, "maxHFEnergyPlus_30_35/F");
  tree->Branch("maxHFEnergyPlus_30_36", &maxHFEnergyPlus_30_36, "maxHFEnergyPlus_30_36/F");
  tree->Branch("maxHFEnergyPlus_30_37", &maxHFEnergyPlus_30_37, "maxHFEnergyPlus_30_37/F");
  tree->Branch("maxHFEnergyPlus_30_38", &maxHFEnergyPlus_30_38, "maxHFEnergyPlus_30_38/F");
  tree->Branch("maxHFEnergyPlus_30_39", &maxHFEnergyPlus_30_39, "maxHFEnergyPlus_30_39/F");
  tree->Branch("maxHFEnergyPlus_30_40", &maxHFEnergyPlus_30_40, "maxHFEnergyPlus_30_40/F");
  tree->Branch("maxHFEnergyPlus_30_41", &maxHFEnergyPlus_30_41, "maxHFEnergyPlus_30_41/F");
  tree->Branch("maxHFEnergyPlus_30_42", &maxHFEnergyPlus_30_42, "maxHFEnergyPlus_30_42/F");
  tree->Branch("maxHFEnergyPlus_30_43", &maxHFEnergyPlus_30_43, "maxHFEnergyPlus_30_43/F");
  tree->Branch("maxHFEnergyPlus_30_44", &maxHFEnergyPlus_30_44, "maxHFEnergyPlus_30_44/F");
  tree->Branch("maxHFEnergyPlus_30_45", &maxHFEnergyPlus_30_45, "maxHFEnergyPlus_30_45/F");
  tree->Branch("maxHFEnergyPlus_30_46", &maxHFEnergyPlus_30_46, "maxHFEnergyPlus_30_46/F");
  tree->Branch("maxHFEnergyPlus_30_47", &maxHFEnergyPlus_30_47, "maxHFEnergyPlus_30_47/F");
  tree->Branch("maxHFEnergyPlus_30_48", &maxHFEnergyPlus_30_48, "maxHFEnergyPlus_30_48/F");
  tree->Branch("maxHFEnergyPlus_30_49", &maxHFEnergyPlus_30_49, "maxHFEnergyPlus_30_49/F");
  tree->Branch("maxHFEnergyPlus_30_50", &maxHFEnergyPlus_30_50, "maxHFEnergyPlus_30_50/F");
  tree->Branch("maxHFEnergyPlus_30_51", &maxHFEnergyPlus_30_51, "maxHFEnergyPlus_30_51/F");
  tree->Branch("maxHFEnergyPlus_30_52", &maxHFEnergyPlus_30_52, "maxHFEnergyPlus_30_52/F");

  tree->Branch("maxHFEnergyPlus_31_32", &maxHFEnergyPlus_31_32, "maxHFEnergyPlus_31_32/F");
  tree->Branch("maxHFEnergyPlus_31_33", &maxHFEnergyPlus_31_33, "maxHFEnergyPlus_31_33/F");
  tree->Branch("maxHFEnergyPlus_31_34", &maxHFEnergyPlus_31_34, "maxHFEnergyPlus_31_34/F");
  tree->Branch("maxHFEnergyPlus_31_35", &maxHFEnergyPlus_31_35, "maxHFEnergyPlus_31_35/F");
  tree->Branch("maxHFEnergyPlus_31_36", &maxHFEnergyPlus_31_36, "maxHFEnergyPlus_31_36/F");
  tree->Branch("maxHFEnergyPlus_31_37", &maxHFEnergyPlus_31_37, "maxHFEnergyPlus_31_37/F");
  tree->Branch("maxHFEnergyPlus_31_38", &maxHFEnergyPlus_31_38, "maxHFEnergyPlus_31_38/F");
  tree->Branch("maxHFEnergyPlus_31_39", &maxHFEnergyPlus_31_39, "maxHFEnergyPlus_31_39/F");
  tree->Branch("maxHFEnergyPlus_31_40", &maxHFEnergyPlus_31_40, "maxHFEnergyPlus_31_40/F");
  tree->Branch("maxHFEnergyPlus_31_41", &maxHFEnergyPlus_31_41, "maxHFEnergyPlus_31_41/F");
  tree->Branch("maxHFEnergyPlus_31_42", &maxHFEnergyPlus_31_42, "maxHFEnergyPlus_31_42/F");
  tree->Branch("maxHFEnergyPlus_31_43", &maxHFEnergyPlus_31_43, "maxHFEnergyPlus_31_43/F");
  tree->Branch("maxHFEnergyPlus_31_44", &maxHFEnergyPlus_31_44, "maxHFEnergyPlus_31_44/F");
  tree->Branch("maxHFEnergyPlus_31_45", &maxHFEnergyPlus_31_45, "maxHFEnergyPlus_31_45/F");
  tree->Branch("maxHFEnergyPlus_31_46", &maxHFEnergyPlus_31_46, "maxHFEnergyPlus_31_46/F");
  tree->Branch("maxHFEnergyPlus_31_47", &maxHFEnergyPlus_31_47, "maxHFEnergyPlus_31_47/F");
  tree->Branch("maxHFEnergyPlus_31_48", &maxHFEnergyPlus_31_48, "maxHFEnergyPlus_31_48/F");
  tree->Branch("maxHFEnergyPlus_31_49", &maxHFEnergyPlus_31_49, "maxHFEnergyPlus_31_49/F");
  tree->Branch("maxHFEnergyPlus_31_50", &maxHFEnergyPlus_31_50, "maxHFEnergyPlus_31_50/F");
  tree->Branch("maxHFEnergyPlus_31_51", &maxHFEnergyPlus_31_51, "maxHFEnergyPlus_31_51/F");
  tree->Branch("maxHFEnergyPlus_31_52", &maxHFEnergyPlus_31_52, "maxHFEnergyPlus_31_52/F");

  tree->Branch("maxHFEnergyPlus_32_33", &maxHFEnergyPlus_32_33, "maxHFEnergyPlus_32_33/F");
  tree->Branch("maxHFEnergyPlus_32_34", &maxHFEnergyPlus_32_34, "maxHFEnergyPlus_32_34/F");
  tree->Branch("maxHFEnergyPlus_32_35", &maxHFEnergyPlus_32_35, "maxHFEnergyPlus_32_35/F");
  tree->Branch("maxHFEnergyPlus_32_36", &maxHFEnergyPlus_32_36, "maxHFEnergyPlus_32_36/F");
  tree->Branch("maxHFEnergyPlus_32_37", &maxHFEnergyPlus_32_37, "maxHFEnergyPlus_32_37/F");
  tree->Branch("maxHFEnergyPlus_32_38", &maxHFEnergyPlus_32_38, "maxHFEnergyPlus_32_38/F");
  tree->Branch("maxHFEnergyPlus_32_39", &maxHFEnergyPlus_32_39, "maxHFEnergyPlus_32_39/F");
  tree->Branch("maxHFEnergyPlus_32_40", &maxHFEnergyPlus_32_40, "maxHFEnergyPlus_32_40/F");
  tree->Branch("maxHFEnergyPlus_32_41", &maxHFEnergyPlus_32_41, "maxHFEnergyPlus_32_41/F");
  tree->Branch("maxHFEnergyPlus_32_42", &maxHFEnergyPlus_32_42, "maxHFEnergyPlus_32_42/F");
  tree->Branch("maxHFEnergyPlus_32_43", &maxHFEnergyPlus_32_43, "maxHFEnergyPlus_32_43/F");
  tree->Branch("maxHFEnergyPlus_32_44", &maxHFEnergyPlus_32_44, "maxHFEnergyPlus_32_44/F");
  tree->Branch("maxHFEnergyPlus_32_45", &maxHFEnergyPlus_32_45, "maxHFEnergyPlus_32_45/F");
  tree->Branch("maxHFEnergyPlus_32_46", &maxHFEnergyPlus_32_46, "maxHFEnergyPlus_32_46/F");
  tree->Branch("maxHFEnergyPlus_32_47", &maxHFEnergyPlus_32_47, "maxHFEnergyPlus_32_47/F");
  tree->Branch("maxHFEnergyPlus_32_48", &maxHFEnergyPlus_32_48, "maxHFEnergyPlus_32_48/F");
  tree->Branch("maxHFEnergyPlus_32_49", &maxHFEnergyPlus_32_49, "maxHFEnergyPlus_32_49/F");
  tree->Branch("maxHFEnergyPlus_32_50", &maxHFEnergyPlus_32_50, "maxHFEnergyPlus_32_50/F");
  tree->Branch("maxHFEnergyPlus_32_51", &maxHFEnergyPlus_32_51, "maxHFEnergyPlus_32_51/F");
  tree->Branch("maxHFEnergyPlus_32_52", &maxHFEnergyPlus_32_52, "maxHFEnergyPlus_32_52/F");

  tree->Branch("maxHFEnergyPlus_33_34", &maxHFEnergyPlus_33_34, "maxHFEnergyPlus_33_34/F");
  tree->Branch("maxHFEnergyPlus_33_35", &maxHFEnergyPlus_33_35, "maxHFEnergyPlus_33_35/F");
  tree->Branch("maxHFEnergyPlus_33_36", &maxHFEnergyPlus_33_36, "maxHFEnergyPlus_33_36/F");
  tree->Branch("maxHFEnergyPlus_33_37", &maxHFEnergyPlus_33_37, "maxHFEnergyPlus_33_37/F");
  tree->Branch("maxHFEnergyPlus_33_38", &maxHFEnergyPlus_33_38, "maxHFEnergyPlus_33_38/F");
  tree->Branch("maxHFEnergyPlus_33_39", &maxHFEnergyPlus_33_39, "maxHFEnergyPlus_33_39/F");
  tree->Branch("maxHFEnergyPlus_33_40", &maxHFEnergyPlus_33_40, "maxHFEnergyPlus_33_40/F");
  tree->Branch("maxHFEnergyPlus_33_41", &maxHFEnergyPlus_33_41, "maxHFEnergyPlus_33_41/F");
  tree->Branch("maxHFEnergyPlus_33_42", &maxHFEnergyPlus_33_42, "maxHFEnergyPlus_33_42/F");
  tree->Branch("maxHFEnergyPlus_33_43", &maxHFEnergyPlus_33_43, "maxHFEnergyPlus_33_43/F");
  tree->Branch("maxHFEnergyPlus_33_44", &maxHFEnergyPlus_33_44, "maxHFEnergyPlus_33_44/F");
  tree->Branch("maxHFEnergyPlus_33_45", &maxHFEnergyPlus_33_45, "maxHFEnergyPlus_33_45/F");
  tree->Branch("maxHFEnergyPlus_33_46", &maxHFEnergyPlus_33_46, "maxHFEnergyPlus_33_46/F");
  tree->Branch("maxHFEnergyPlus_33_47", &maxHFEnergyPlus_33_47, "maxHFEnergyPlus_33_47/F");
  tree->Branch("maxHFEnergyPlus_33_48", &maxHFEnergyPlus_33_48, "maxHFEnergyPlus_33_48/F");
  tree->Branch("maxHFEnergyPlus_33_49", &maxHFEnergyPlus_33_49, "maxHFEnergyPlus_33_49/F");
  tree->Branch("maxHFEnergyPlus_33_50", &maxHFEnergyPlus_33_50, "maxHFEnergyPlus_33_50/F");
  tree->Branch("maxHFEnergyPlus_33_51", &maxHFEnergyPlus_33_51, "maxHFEnergyPlus_33_51/F");
  tree->Branch("maxHFEnergyPlus_33_52", &maxHFEnergyPlus_33_52, "maxHFEnergyPlus_33_52/F");

  tree->Branch("maxHFEnergyPlus_34_35", &maxHFEnergyPlus_34_35, "maxHFEnergyPlus_34_35/F");
  tree->Branch("maxHFEnergyPlus_34_36", &maxHFEnergyPlus_34_36, "maxHFEnergyPlus_34_36/F");
  tree->Branch("maxHFEnergyPlus_34_37", &maxHFEnergyPlus_34_37, "maxHFEnergyPlus_34_37/F");
  tree->Branch("maxHFEnergyPlus_34_38", &maxHFEnergyPlus_34_38, "maxHFEnergyPlus_34_38/F");
  tree->Branch("maxHFEnergyPlus_34_39", &maxHFEnergyPlus_34_39, "maxHFEnergyPlus_34_39/F");
  tree->Branch("maxHFEnergyPlus_34_40", &maxHFEnergyPlus_34_40, "maxHFEnergyPlus_34_40/F");
  tree->Branch("maxHFEnergyPlus_34_41", &maxHFEnergyPlus_34_41, "maxHFEnergyPlus_34_41/F");
  tree->Branch("maxHFEnergyPlus_34_42", &maxHFEnergyPlus_34_42, "maxHFEnergyPlus_34_42/F");
  tree->Branch("maxHFEnergyPlus_34_43", &maxHFEnergyPlus_34_43, "maxHFEnergyPlus_34_43/F");
  tree->Branch("maxHFEnergyPlus_34_44", &maxHFEnergyPlus_34_44, "maxHFEnergyPlus_34_44/F");
  tree->Branch("maxHFEnergyPlus_34_45", &maxHFEnergyPlus_34_45, "maxHFEnergyPlus_34_45/F");
  tree->Branch("maxHFEnergyPlus_34_46", &maxHFEnergyPlus_34_46, "maxHFEnergyPlus_34_46/F");
  tree->Branch("maxHFEnergyPlus_34_47", &maxHFEnergyPlus_34_47, "maxHFEnergyPlus_34_47/F");
  tree->Branch("maxHFEnergyPlus_34_48", &maxHFEnergyPlus_34_48, "maxHFEnergyPlus_34_48/F");
  tree->Branch("maxHFEnergyPlus_34_49", &maxHFEnergyPlus_34_49, "maxHFEnergyPlus_34_49/F");
  tree->Branch("maxHFEnergyPlus_34_50", &maxHFEnergyPlus_34_50, "maxHFEnergyPlus_34_50/F");
  tree->Branch("maxHFEnergyPlus_34_51", &maxHFEnergyPlus_34_51, "maxHFEnergyPlus_34_51/F");
  tree->Branch("maxHFEnergyPlus_34_52", &maxHFEnergyPlus_34_52, "maxHFEnergyPlus_34_52/F");

  tree->Branch("maxHFEnergyPlus_35_36", &maxHFEnergyPlus_35_36, "maxHFEnergyPlus_35_36/F");
  tree->Branch("maxHFEnergyPlus_35_37", &maxHFEnergyPlus_35_37, "maxHFEnergyPlus_35_37/F");
  tree->Branch("maxHFEnergyPlus_35_38", &maxHFEnergyPlus_35_38, "maxHFEnergyPlus_35_38/F");
  tree->Branch("maxHFEnergyPlus_35_39", &maxHFEnergyPlus_35_39, "maxHFEnergyPlus_35_39/F");
  tree->Branch("maxHFEnergyPlus_35_40", &maxHFEnergyPlus_35_40, "maxHFEnergyPlus_35_40/F");
  tree->Branch("maxHFEnergyPlus_35_41", &maxHFEnergyPlus_35_41, "maxHFEnergyPlus_35_41/F");
  tree->Branch("maxHFEnergyPlus_35_42", &maxHFEnergyPlus_35_42, "maxHFEnergyPlus_35_42/F");
  tree->Branch("maxHFEnergyPlus_35_43", &maxHFEnergyPlus_35_43, "maxHFEnergyPlus_35_43/F");
  tree->Branch("maxHFEnergyPlus_35_44", &maxHFEnergyPlus_35_44, "maxHFEnergyPlus_35_44/F");
  tree->Branch("maxHFEnergyPlus_35_45", &maxHFEnergyPlus_35_45, "maxHFEnergyPlus_35_45/F");
  tree->Branch("maxHFEnergyPlus_35_46", &maxHFEnergyPlus_35_46, "maxHFEnergyPlus_35_46/F");
  tree->Branch("maxHFEnergyPlus_35_47", &maxHFEnergyPlus_35_47, "maxHFEnergyPlus_35_47/F");
  tree->Branch("maxHFEnergyPlus_35_48", &maxHFEnergyPlus_35_48, "maxHFEnergyPlus_35_48/F");
  tree->Branch("maxHFEnergyPlus_35_49", &maxHFEnergyPlus_35_49, "maxHFEnergyPlus_35_49/F");
  tree->Branch("maxHFEnergyPlus_35_50", &maxHFEnergyPlus_35_50, "maxHFEnergyPlus_35_50/F");
  tree->Branch("maxHFEnergyPlus_35_51", &maxHFEnergyPlus_35_51, "maxHFEnergyPlus_35_51/F");
  tree->Branch("maxHFEnergyPlus_35_52", &maxHFEnergyPlus_35_52, "maxHFEnergyPlus_35_52/F");

  tree->Branch("maxHFEnergyPlus_36_37", &maxHFEnergyPlus_36_37, "maxHFEnergyPlus_36_37/F");
  tree->Branch("maxHFEnergyPlus_36_38", &maxHFEnergyPlus_36_38, "maxHFEnergyPlus_36_38/F");
  tree->Branch("maxHFEnergyPlus_36_39", &maxHFEnergyPlus_36_39, "maxHFEnergyPlus_36_39/F");
  tree->Branch("maxHFEnergyPlus_36_40", &maxHFEnergyPlus_36_40, "maxHFEnergyPlus_36_40/F");
  tree->Branch("maxHFEnergyPlus_36_41", &maxHFEnergyPlus_36_41, "maxHFEnergyPlus_36_41/F");
  tree->Branch("maxHFEnergyPlus_36_42", &maxHFEnergyPlus_36_42, "maxHFEnergyPlus_36_42/F");
  tree->Branch("maxHFEnergyPlus_36_43", &maxHFEnergyPlus_36_43, "maxHFEnergyPlus_36_43/F");
  tree->Branch("maxHFEnergyPlus_36_44", &maxHFEnergyPlus_36_44, "maxHFEnergyPlus_36_44/F");
  tree->Branch("maxHFEnergyPlus_36_45", &maxHFEnergyPlus_36_45, "maxHFEnergyPlus_36_45/F");
  tree->Branch("maxHFEnergyPlus_36_46", &maxHFEnergyPlus_36_46, "maxHFEnergyPlus_36_46/F");
  tree->Branch("maxHFEnergyPlus_36_47", &maxHFEnergyPlus_36_47, "maxHFEnergyPlus_36_47/F");
  tree->Branch("maxHFEnergyPlus_36_48", &maxHFEnergyPlus_36_48, "maxHFEnergyPlus_36_48/F");
  tree->Branch("maxHFEnergyPlus_36_49", &maxHFEnergyPlus_36_49, "maxHFEnergyPlus_36_49/F");
  tree->Branch("maxHFEnergyPlus_36_50", &maxHFEnergyPlus_36_50, "maxHFEnergyPlus_36_50/F");
  tree->Branch("maxHFEnergyPlus_36_51", &maxHFEnergyPlus_36_51, "maxHFEnergyPlus_36_51/F");
  tree->Branch("maxHFEnergyPlus_36_52", &maxHFEnergyPlus_36_52, "maxHFEnergyPlus_36_52/F");

  tree->Branch("maxHFEnergyPlus_37_38", &maxHFEnergyPlus_37_38, "maxHFEnergyPlus_37_38/F");
  tree->Branch("maxHFEnergyPlus_37_39", &maxHFEnergyPlus_37_39, "maxHFEnergyPlus_37_39/F");
  tree->Branch("maxHFEnergyPlus_37_40", &maxHFEnergyPlus_37_40, "maxHFEnergyPlus_37_40/F");
  tree->Branch("maxHFEnergyPlus_37_41", &maxHFEnergyPlus_37_41, "maxHFEnergyPlus_37_41/F");
  tree->Branch("maxHFEnergyPlus_37_42", &maxHFEnergyPlus_37_42, "maxHFEnergyPlus_37_42/F");
  tree->Branch("maxHFEnergyPlus_37_43", &maxHFEnergyPlus_37_43, "maxHFEnergyPlus_37_43/F");
  tree->Branch("maxHFEnergyPlus_37_44", &maxHFEnergyPlus_37_44, "maxHFEnergyPlus_37_44/F");
  tree->Branch("maxHFEnergyPlus_37_45", &maxHFEnergyPlus_37_45, "maxHFEnergyPlus_37_45/F");
  tree->Branch("maxHFEnergyPlus_37_46", &maxHFEnergyPlus_37_46, "maxHFEnergyPlus_37_46/F");
  tree->Branch("maxHFEnergyPlus_37_47", &maxHFEnergyPlus_37_47, "maxHFEnergyPlus_37_47/F");
  tree->Branch("maxHFEnergyPlus_37_48", &maxHFEnergyPlus_37_48, "maxHFEnergyPlus_37_48/F");
  tree->Branch("maxHFEnergyPlus_37_49", &maxHFEnergyPlus_37_49, "maxHFEnergyPlus_37_49/F");
  tree->Branch("maxHFEnergyPlus_37_50", &maxHFEnergyPlus_37_50, "maxHFEnergyPlus_37_50/F");
  tree->Branch("maxHFEnergyPlus_37_51", &maxHFEnergyPlus_37_51, "maxHFEnergyPlus_37_51/F");
  tree->Branch("maxHFEnergyPlus_37_52", &maxHFEnergyPlus_37_52, "maxHFEnergyPlus_37_52/F");

  tree->Branch("maxHFEnergyPlus_38_39", &maxHFEnergyPlus_38_39, "maxHFEnergyPlus_38_39/F");
  tree->Branch("maxHFEnergyPlus_38_40", &maxHFEnergyPlus_38_40, "maxHFEnergyPlus_38_40/F");
  tree->Branch("maxHFEnergyPlus_38_41", &maxHFEnergyPlus_38_41, "maxHFEnergyPlus_38_41/F");
  tree->Branch("maxHFEnergyPlus_38_42", &maxHFEnergyPlus_38_42, "maxHFEnergyPlus_38_42/F");
  tree->Branch("maxHFEnergyPlus_38_43", &maxHFEnergyPlus_38_43, "maxHFEnergyPlus_38_43/F");
  tree->Branch("maxHFEnergyPlus_38_44", &maxHFEnergyPlus_38_44, "maxHFEnergyPlus_38_44/F");
  tree->Branch("maxHFEnergyPlus_38_45", &maxHFEnergyPlus_38_45, "maxHFEnergyPlus_38_45/F");
  tree->Branch("maxHFEnergyPlus_38_46", &maxHFEnergyPlus_38_46, "maxHFEnergyPlus_38_46/F");
  tree->Branch("maxHFEnergyPlus_38_47", &maxHFEnergyPlus_38_47, "maxHFEnergyPlus_38_47/F");
  tree->Branch("maxHFEnergyPlus_38_48", &maxHFEnergyPlus_38_48, "maxHFEnergyPlus_38_48/F");
  tree->Branch("maxHFEnergyPlus_38_49", &maxHFEnergyPlus_38_49, "maxHFEnergyPlus_38_49/F");
  tree->Branch("maxHFEnergyPlus_38_50", &maxHFEnergyPlus_38_50, "maxHFEnergyPlus_38_50/F");
  tree->Branch("maxHFEnergyPlus_38_51", &maxHFEnergyPlus_38_51, "maxHFEnergyPlus_38_51/F");
  tree->Branch("maxHFEnergyPlus_38_52", &maxHFEnergyPlus_38_52, "maxHFEnergyPlus_38_52/F");

  tree->Branch("maxHFEnergyPlus_39_40", &maxHFEnergyPlus_39_40, "maxHFEnergyPlus_39_40/F");
  tree->Branch("maxHFEnergyPlus_39_41", &maxHFEnergyPlus_39_41, "maxHFEnergyPlus_39_41/F");
  tree->Branch("maxHFEnergyPlus_39_42", &maxHFEnergyPlus_39_42, "maxHFEnergyPlus_39_42/F");
  tree->Branch("maxHFEnergyPlus_39_43", &maxHFEnergyPlus_39_43, "maxHFEnergyPlus_39_43/F");
  tree->Branch("maxHFEnergyPlus_39_44", &maxHFEnergyPlus_39_44, "maxHFEnergyPlus_39_44/F");
  tree->Branch("maxHFEnergyPlus_39_45", &maxHFEnergyPlus_39_45, "maxHFEnergyPlus_39_45/F");
  tree->Branch("maxHFEnergyPlus_39_46", &maxHFEnergyPlus_39_46, "maxHFEnergyPlus_39_46/F");
  tree->Branch("maxHFEnergyPlus_39_47", &maxHFEnergyPlus_39_47, "maxHFEnergyPlus_39_47/F");
  tree->Branch("maxHFEnergyPlus_39_48", &maxHFEnergyPlus_39_48, "maxHFEnergyPlus_39_48/F");
  tree->Branch("maxHFEnergyPlus_39_49", &maxHFEnergyPlus_39_49, "maxHFEnergyPlus_39_49/F");
  tree->Branch("maxHFEnergyPlus_39_50", &maxHFEnergyPlus_39_50, "maxHFEnergyPlus_39_50/F");
  tree->Branch("maxHFEnergyPlus_39_51", &maxHFEnergyPlus_39_51, "maxHFEnergyPlus_39_51/F");
  tree->Branch("maxHFEnergyPlus_39_52", &maxHFEnergyPlus_39_52, "maxHFEnergyPlus_39_52/F");

  tree->Branch("maxHFEnergyPlus_40_41", &maxHFEnergyPlus_40_41, "maxHFEnergyPlus_40_41/F");
  tree->Branch("maxHFEnergyPlus_40_42", &maxHFEnergyPlus_40_42, "maxHFEnergyPlus_40_42/F");
  tree->Branch("maxHFEnergyPlus_40_43", &maxHFEnergyPlus_40_43, "maxHFEnergyPlus_40_43/F");
  tree->Branch("maxHFEnergyPlus_40_44", &maxHFEnergyPlus_40_44, "maxHFEnergyPlus_40_44/F");
  tree->Branch("maxHFEnergyPlus_40_45", &maxHFEnergyPlus_40_45, "maxHFEnergyPlus_40_45/F");
  tree->Branch("maxHFEnergyPlus_40_46", &maxHFEnergyPlus_40_46, "maxHFEnergyPlus_40_46/F");
  tree->Branch("maxHFEnergyPlus_40_47", &maxHFEnergyPlus_40_47, "maxHFEnergyPlus_40_47/F");
  tree->Branch("maxHFEnergyPlus_40_48", &maxHFEnergyPlus_40_48, "maxHFEnergyPlus_40_48/F");
  tree->Branch("maxHFEnergyPlus_40_49", &maxHFEnergyPlus_40_49, "maxHFEnergyPlus_40_49/F");
  tree->Branch("maxHFEnergyPlus_40_50", &maxHFEnergyPlus_40_50, "maxHFEnergyPlus_40_50/F");
  tree->Branch("maxHFEnergyPlus_40_51", &maxHFEnergyPlus_40_51, "maxHFEnergyPlus_40_51/F");
  tree->Branch("maxHFEnergyPlus_40_52", &maxHFEnergyPlus_40_52, "maxHFEnergyPlus_40_52/F");

  tree->Branch("maxHFEnergyPlus_41_42", &maxHFEnergyPlus_41_42, "maxHFEnergyPlus_41_42/F");
  tree->Branch("maxHFEnergyPlus_41_43", &maxHFEnergyPlus_41_43, "maxHFEnergyPlus_41_43/F");
  tree->Branch("maxHFEnergyPlus_41_44", &maxHFEnergyPlus_41_44, "maxHFEnergyPlus_41_44/F");
  tree->Branch("maxHFEnergyPlus_41_45", &maxHFEnergyPlus_41_45, "maxHFEnergyPlus_41_45/F");
  tree->Branch("maxHFEnergyPlus_41_46", &maxHFEnergyPlus_41_46, "maxHFEnergyPlus_41_46/F");
  tree->Branch("maxHFEnergyPlus_41_47", &maxHFEnergyPlus_41_47, "maxHFEnergyPlus_41_47/F");
  tree->Branch("maxHFEnergyPlus_41_48", &maxHFEnergyPlus_41_48, "maxHFEnergyPlus_41_48/F");
  tree->Branch("maxHFEnergyPlus_41_49", &maxHFEnergyPlus_41_49, "maxHFEnergyPlus_41_49/F");
  tree->Branch("maxHFEnergyPlus_41_50", &maxHFEnergyPlus_41_50, "maxHFEnergyPlus_41_50/F");
  tree->Branch("maxHFEnergyPlus_41_51", &maxHFEnergyPlus_41_51, "maxHFEnergyPlus_41_51/F");
  tree->Branch("maxHFEnergyPlus_41_52", &maxHFEnergyPlus_41_52, "maxHFEnergyPlus_41_52/F");

  tree->Branch("maxHFEnergyPlus_42_43", &maxHFEnergyPlus_42_43, "maxHFEnergyPlus_42_43/F");
  tree->Branch("maxHFEnergyPlus_42_44", &maxHFEnergyPlus_42_44, "maxHFEnergyPlus_42_44/F");
  tree->Branch("maxHFEnergyPlus_42_45", &maxHFEnergyPlus_42_45, "maxHFEnergyPlus_42_45/F");
  tree->Branch("maxHFEnergyPlus_42_46", &maxHFEnergyPlus_42_46, "maxHFEnergyPlus_42_46/F");
  tree->Branch("maxHFEnergyPlus_42_47", &maxHFEnergyPlus_42_47, "maxHFEnergyPlus_42_47/F");
  tree->Branch("maxHFEnergyPlus_42_48", &maxHFEnergyPlus_42_48, "maxHFEnergyPlus_42_48/F");
  tree->Branch("maxHFEnergyPlus_42_49", &maxHFEnergyPlus_42_49, "maxHFEnergyPlus_42_49/F");
  tree->Branch("maxHFEnergyPlus_42_50", &maxHFEnergyPlus_42_50, "maxHFEnergyPlus_42_50/F");
  tree->Branch("maxHFEnergyPlus_42_51", &maxHFEnergyPlus_42_51, "maxHFEnergyPlus_42_51/F");
  tree->Branch("maxHFEnergyPlus_42_52", &maxHFEnergyPlus_42_52, "maxHFEnergyPlus_42_52/F");

  tree->Branch("maxHFEnergyPlus_43_44", &maxHFEnergyPlus_43_44, "maxHFEnergyPlus_43_44/F");
  tree->Branch("maxHFEnergyPlus_43_45", &maxHFEnergyPlus_43_45, "maxHFEnergyPlus_43_45/F");
  tree->Branch("maxHFEnergyPlus_43_46", &maxHFEnergyPlus_43_46, "maxHFEnergyPlus_43_46/F");
  tree->Branch("maxHFEnergyPlus_43_47", &maxHFEnergyPlus_43_47, "maxHFEnergyPlus_43_47/F");
  tree->Branch("maxHFEnergyPlus_43_48", &maxHFEnergyPlus_43_48, "maxHFEnergyPlus_43_48/F");
  tree->Branch("maxHFEnergyPlus_43_49", &maxHFEnergyPlus_43_49, "maxHFEnergyPlus_43_49/F");
  tree->Branch("maxHFEnergyPlus_43_50", &maxHFEnergyPlus_43_50, "maxHFEnergyPlus_43_50/F");
  tree->Branch("maxHFEnergyPlus_43_51", &maxHFEnergyPlus_43_51, "maxHFEnergyPlus_43_51/F");
  tree->Branch("maxHFEnergyPlus_43_52", &maxHFEnergyPlus_43_52, "maxHFEnergyPlus_43_52/F");

  tree->Branch("maxHFEnergyPlus_44_45", &maxHFEnergyPlus_44_45, "maxHFEnergyPlus_44_45/F");
  tree->Branch("maxHFEnergyPlus_44_46", &maxHFEnergyPlus_44_46, "maxHFEnergyPlus_44_46/F");
  tree->Branch("maxHFEnergyPlus_44_47", &maxHFEnergyPlus_44_47, "maxHFEnergyPlus_44_47/F");
  tree->Branch("maxHFEnergyPlus_44_48", &maxHFEnergyPlus_44_48, "maxHFEnergyPlus_44_48/F");
  tree->Branch("maxHFEnergyPlus_44_49", &maxHFEnergyPlus_44_49, "maxHFEnergyPlus_44_49/F");
  tree->Branch("maxHFEnergyPlus_44_50", &maxHFEnergyPlus_44_50, "maxHFEnergyPlus_44_50/F");
  tree->Branch("maxHFEnergyPlus_44_51", &maxHFEnergyPlus_44_51, "maxHFEnergyPlus_44_51/F");
  tree->Branch("maxHFEnergyPlus_44_52", &maxHFEnergyPlus_44_52, "maxHFEnergyPlus_44_52/F");

  tree->Branch("maxHFEnergyPlus_45_46", &maxHFEnergyPlus_45_46, "maxHFEnergyPlus_45_46/F");
  tree->Branch("maxHFEnergyPlus_45_47", &maxHFEnergyPlus_45_47, "maxHFEnergyPlus_45_47/F");
  tree->Branch("maxHFEnergyPlus_45_48", &maxHFEnergyPlus_45_48, "maxHFEnergyPlus_45_48/F");
  tree->Branch("maxHFEnergyPlus_45_49", &maxHFEnergyPlus_45_49, "maxHFEnergyPlus_45_49/F");
  tree->Branch("maxHFEnergyPlus_45_50", &maxHFEnergyPlus_45_50, "maxHFEnergyPlus_45_50/F");
  tree->Branch("maxHFEnergyPlus_45_51", &maxHFEnergyPlus_45_51, "maxHFEnergyPlus_45_51/F");
  tree->Branch("maxHFEnergyPlus_45_52", &maxHFEnergyPlus_45_52, "maxHFEnergyPlus_45_52/F");

  tree->Branch("maxHFEnergyPlus_46_47", &maxHFEnergyPlus_46_47, "maxHFEnergyPlus_46_47/F");
  tree->Branch("maxHFEnergyPlus_46_48", &maxHFEnergyPlus_46_48, "maxHFEnergyPlus_46_48/F");
  tree->Branch("maxHFEnergyPlus_46_49", &maxHFEnergyPlus_46_49, "maxHFEnergyPlus_46_49/F");
  tree->Branch("maxHFEnergyPlus_46_50", &maxHFEnergyPlus_46_50, "maxHFEnergyPlus_46_50/F");
  tree->Branch("maxHFEnergyPlus_46_51", &maxHFEnergyPlus_46_51, "maxHFEnergyPlus_46_51/F");
  tree->Branch("maxHFEnergyPlus_46_52", &maxHFEnergyPlus_46_52, "maxHFEnergyPlus_46_52/F");

  tree->Branch("maxHFEnergyPlus_47_48", &maxHFEnergyPlus_47_48, "maxHFEnergyPlus_47_48/F");
  tree->Branch("maxHFEnergyPlus_47_49", &maxHFEnergyPlus_47_49, "maxHFEnergyPlus_47_49/F");
  tree->Branch("maxHFEnergyPlus_47_50", &maxHFEnergyPlus_47_50, "maxHFEnergyPlus_47_50/F");
  tree->Branch("maxHFEnergyPlus_47_51", &maxHFEnergyPlus_47_51, "maxHFEnergyPlus_47_51/F");
  tree->Branch("maxHFEnergyPlus_47_52", &maxHFEnergyPlus_47_52, "maxHFEnergyPlus_47_52/F");

  tree->Branch("maxHFEnergyPlus_48_49", &maxHFEnergyPlus_48_49, "maxHFEnergyPlus_48_49/F");
  tree->Branch("maxHFEnergyPlus_48_50", &maxHFEnergyPlus_48_50, "maxHFEnergyPlus_48_50/F");
  tree->Branch("maxHFEnergyPlus_48_51", &maxHFEnergyPlus_48_51, "maxHFEnergyPlus_48_51/F");
  tree->Branch("maxHFEnergyPlus_48_52", &maxHFEnergyPlus_48_52, "maxHFEnergyPlus_48_52/F");

  tree->Branch("maxHFEnergyPlus_49_50", &maxHFEnergyPlus_49_50, "maxHFEnergyPlus_49_50/F");
  tree->Branch("maxHFEnergyPlus_49_51", &maxHFEnergyPlus_49_51, "maxHFEnergyPlus_49_51/F");
  tree->Branch("maxHFEnergyPlus_49_52", &maxHFEnergyPlus_49_52, "maxHFEnergyPlus_49_52/F");

  tree->Branch("maxHFEnergyPlus_50_51", &maxHFEnergyPlus_50_51, "maxHFEnergyPlus_50_51/F");
  tree->Branch("maxHFEnergyPlus_50_52", &maxHFEnergyPlus_50_52, "maxHFEnergyPlus_50_52/F");

  tree->Branch("maxHFEnergyPlus_51_52", &maxHFEnergyPlus_51_52, "maxHFEnergyPlus_51_52/F");
  // maxhfenergyplus

  // minhfenergy
  tree->Branch("minHFEnergy_30_31", &minHFEnergy_30_31, "minHFEnergy_30_31/F");
  tree->Branch("minHFEnergy_30_32", &minHFEnergy_30_32, "minHFEnergy_30_32/F");
  tree->Branch("minHFEnergy_30_33", &minHFEnergy_30_33, "minHFEnergy_30_33/F");
  tree->Branch("minHFEnergy_30_34", &minHFEnergy_30_34, "minHFEnergy_30_34/F");
  tree->Branch("minHFEnergy_30_35", &minHFEnergy_30_35, "minHFEnergy_30_35/F");
  tree->Branch("minHFEnergy_30_36", &minHFEnergy_30_36, "minHFEnergy_30_36/F");
  tree->Branch("minHFEnergy_30_37", &minHFEnergy_30_37, "minHFEnergy_30_37/F");
  tree->Branch("minHFEnergy_30_38", &minHFEnergy_30_38, "minHFEnergy_30_38/F");
  tree->Branch("minHFEnergy_30_39", &minHFEnergy_30_39, "minHFEnergy_30_39/F");
  tree->Branch("minHFEnergy_30_40", &minHFEnergy_30_40, "minHFEnergy_30_40/F");
  tree->Branch("minHFEnergy_30_41", &minHFEnergy_30_41, "minHFEnergy_30_41/F");
  tree->Branch("minHFEnergy_30_42", &minHFEnergy_30_42, "minHFEnergy_30_42/F");
  tree->Branch("minHFEnergy_30_43", &minHFEnergy_30_43, "minHFEnergy_30_43/F");
  tree->Branch("minHFEnergy_30_44", &minHFEnergy_30_44, "minHFEnergy_30_44/F");
  tree->Branch("minHFEnergy_30_45", &minHFEnergy_30_45, "minHFEnergy_30_45/F");
  tree->Branch("minHFEnergy_30_46", &minHFEnergy_30_46, "minHFEnergy_30_46/F");
  tree->Branch("minHFEnergy_30_47", &minHFEnergy_30_47, "minHFEnergy_30_47/F");
  tree->Branch("minHFEnergy_30_48", &minHFEnergy_30_48, "minHFEnergy_30_48/F");
  tree->Branch("minHFEnergy_30_49", &minHFEnergy_30_49, "minHFEnergy_30_49/F");
  tree->Branch("minHFEnergy_30_50", &minHFEnergy_30_50, "minHFEnergy_30_50/F");
  tree->Branch("minHFEnergy_30_51", &minHFEnergy_30_51, "minHFEnergy_30_51/F");
  tree->Branch("minHFEnergy_30_52", &minHFEnergy_30_52, "minHFEnergy_30_52/F");

  tree->Branch("minHFEnergy_31_32", &minHFEnergy_31_32, "minHFEnergy_31_32/F");
  tree->Branch("minHFEnergy_31_33", &minHFEnergy_31_33, "minHFEnergy_31_33/F");
  tree->Branch("minHFEnergy_31_34", &minHFEnergy_31_34, "minHFEnergy_31_34/F");
  tree->Branch("minHFEnergy_31_35", &minHFEnergy_31_35, "minHFEnergy_31_35/F");
  tree->Branch("minHFEnergy_31_36", &minHFEnergy_31_36, "minHFEnergy_31_36/F");
  tree->Branch("minHFEnergy_31_37", &minHFEnergy_31_37, "minHFEnergy_31_37/F");
  tree->Branch("minHFEnergy_31_38", &minHFEnergy_31_38, "minHFEnergy_31_38/F");
  tree->Branch("minHFEnergy_31_39", &minHFEnergy_31_39, "minHFEnergy_31_39/F");
  tree->Branch("minHFEnergy_31_40", &minHFEnergy_31_40, "minHFEnergy_31_40/F");
  tree->Branch("minHFEnergy_31_41", &minHFEnergy_31_41, "minHFEnergy_31_41/F");
  tree->Branch("minHFEnergy_31_42", &minHFEnergy_31_42, "minHFEnergy_31_42/F");
  tree->Branch("minHFEnergy_31_43", &minHFEnergy_31_43, "minHFEnergy_31_43/F");
  tree->Branch("minHFEnergy_31_44", &minHFEnergy_31_44, "minHFEnergy_31_44/F");
  tree->Branch("minHFEnergy_31_45", &minHFEnergy_31_45, "minHFEnergy_31_45/F");
  tree->Branch("minHFEnergy_31_46", &minHFEnergy_31_46, "minHFEnergy_31_46/F");
  tree->Branch("minHFEnergy_31_47", &minHFEnergy_31_47, "minHFEnergy_31_47/F");
  tree->Branch("minHFEnergy_31_48", &minHFEnergy_31_48, "minHFEnergy_31_48/F");
  tree->Branch("minHFEnergy_31_49", &minHFEnergy_31_49, "minHFEnergy_31_49/F");
  tree->Branch("minHFEnergy_31_50", &minHFEnergy_31_50, "minHFEnergy_31_50/F");
  tree->Branch("minHFEnergy_31_51", &minHFEnergy_31_51, "minHFEnergy_31_51/F");
  tree->Branch("minHFEnergy_31_52", &minHFEnergy_31_52, "minHFEnergy_31_52/F");

  tree->Branch("minHFEnergy_32_33", &minHFEnergy_32_33, "minHFEnergy_32_33/F");
  tree->Branch("minHFEnergy_32_34", &minHFEnergy_32_34, "minHFEnergy_32_34/F");
  tree->Branch("minHFEnergy_32_35", &minHFEnergy_32_35, "minHFEnergy_32_35/F");
  tree->Branch("minHFEnergy_32_36", &minHFEnergy_32_36, "minHFEnergy_32_36/F");
  tree->Branch("minHFEnergy_32_37", &minHFEnergy_32_37, "minHFEnergy_32_37/F");
  tree->Branch("minHFEnergy_32_38", &minHFEnergy_32_38, "minHFEnergy_32_38/F");
  tree->Branch("minHFEnergy_32_39", &minHFEnergy_32_39, "minHFEnergy_32_39/F");
  tree->Branch("minHFEnergy_32_40", &minHFEnergy_32_40, "minHFEnergy_32_40/F");
  tree->Branch("minHFEnergy_32_41", &minHFEnergy_32_41, "minHFEnergy_32_41/F");
  tree->Branch("minHFEnergy_32_42", &minHFEnergy_32_42, "minHFEnergy_32_42/F");
  tree->Branch("minHFEnergy_32_43", &minHFEnergy_32_43, "minHFEnergy_32_43/F");
  tree->Branch("minHFEnergy_32_44", &minHFEnergy_32_44, "minHFEnergy_32_44/F");
  tree->Branch("minHFEnergy_32_45", &minHFEnergy_32_45, "minHFEnergy_32_45/F");
  tree->Branch("minHFEnergy_32_46", &minHFEnergy_32_46, "minHFEnergy_32_46/F");
  tree->Branch("minHFEnergy_32_47", &minHFEnergy_32_47, "minHFEnergy_32_47/F");
  tree->Branch("minHFEnergy_32_48", &minHFEnergy_32_48, "minHFEnergy_32_48/F");
  tree->Branch("minHFEnergy_32_49", &minHFEnergy_32_49, "minHFEnergy_32_49/F");
  tree->Branch("minHFEnergy_32_50", &minHFEnergy_32_50, "minHFEnergy_32_50/F");
  tree->Branch("minHFEnergy_32_51", &minHFEnergy_32_51, "minHFEnergy_32_51/F");
  tree->Branch("minHFEnergy_32_52", &minHFEnergy_32_52, "minHFEnergy_32_52/F");

  tree->Branch("minHFEnergy_33_34", &minHFEnergy_33_34, "minHFEnergy_33_34/F");
  tree->Branch("minHFEnergy_33_35", &minHFEnergy_33_35, "minHFEnergy_33_35/F");
  tree->Branch("minHFEnergy_33_36", &minHFEnergy_33_36, "minHFEnergy_33_36/F");
  tree->Branch("minHFEnergy_33_37", &minHFEnergy_33_37, "minHFEnergy_33_37/F");
  tree->Branch("minHFEnergy_33_38", &minHFEnergy_33_38, "minHFEnergy_33_38/F");
  tree->Branch("minHFEnergy_33_39", &minHFEnergy_33_39, "minHFEnergy_33_39/F");
  tree->Branch("minHFEnergy_33_40", &minHFEnergy_33_40, "minHFEnergy_33_40/F");
  tree->Branch("minHFEnergy_33_41", &minHFEnergy_33_41, "minHFEnergy_33_41/F");
  tree->Branch("minHFEnergy_33_42", &minHFEnergy_33_42, "minHFEnergy_33_42/F");
  tree->Branch("minHFEnergy_33_43", &minHFEnergy_33_43, "minHFEnergy_33_43/F");
  tree->Branch("minHFEnergy_33_44", &minHFEnergy_33_44, "minHFEnergy_33_44/F");
  tree->Branch("minHFEnergy_33_45", &minHFEnergy_33_45, "minHFEnergy_33_45/F");
  tree->Branch("minHFEnergy_33_46", &minHFEnergy_33_46, "minHFEnergy_33_46/F");
  tree->Branch("minHFEnergy_33_47", &minHFEnergy_33_47, "minHFEnergy_33_47/F");
  tree->Branch("minHFEnergy_33_48", &minHFEnergy_33_48, "minHFEnergy_33_48/F");
  tree->Branch("minHFEnergy_33_49", &minHFEnergy_33_49, "minHFEnergy_33_49/F");
  tree->Branch("minHFEnergy_33_50", &minHFEnergy_33_50, "minHFEnergy_33_50/F");
  tree->Branch("minHFEnergy_33_51", &minHFEnergy_33_51, "minHFEnergy_33_51/F");
  tree->Branch("minHFEnergy_33_52", &minHFEnergy_33_52, "minHFEnergy_33_52/F");

  tree->Branch("minHFEnergy_34_35", &minHFEnergy_34_35, "minHFEnergy_34_35/F");
  tree->Branch("minHFEnergy_34_36", &minHFEnergy_34_36, "minHFEnergy_34_36/F");
  tree->Branch("minHFEnergy_34_37", &minHFEnergy_34_37, "minHFEnergy_34_37/F");
  tree->Branch("minHFEnergy_34_38", &minHFEnergy_34_38, "minHFEnergy_34_38/F");
  tree->Branch("minHFEnergy_34_39", &minHFEnergy_34_39, "minHFEnergy_34_39/F");
  tree->Branch("minHFEnergy_34_40", &minHFEnergy_34_40, "minHFEnergy_34_40/F");
  tree->Branch("minHFEnergy_34_41", &minHFEnergy_34_41, "minHFEnergy_34_41/F");
  tree->Branch("minHFEnergy_34_42", &minHFEnergy_34_42, "minHFEnergy_34_42/F");
  tree->Branch("minHFEnergy_34_43", &minHFEnergy_34_43, "minHFEnergy_34_43/F");
  tree->Branch("minHFEnergy_34_44", &minHFEnergy_34_44, "minHFEnergy_34_44/F");
  tree->Branch("minHFEnergy_34_45", &minHFEnergy_34_45, "minHFEnergy_34_45/F");
  tree->Branch("minHFEnergy_34_46", &minHFEnergy_34_46, "minHFEnergy_34_46/F");
  tree->Branch("minHFEnergy_34_47", &minHFEnergy_34_47, "minHFEnergy_34_47/F");
  tree->Branch("minHFEnergy_34_48", &minHFEnergy_34_48, "minHFEnergy_34_48/F");
  tree->Branch("minHFEnergy_34_49", &minHFEnergy_34_49, "minHFEnergy_34_49/F");
  tree->Branch("minHFEnergy_34_50", &minHFEnergy_34_50, "minHFEnergy_34_50/F");
  tree->Branch("minHFEnergy_34_51", &minHFEnergy_34_51, "minHFEnergy_34_51/F");
  tree->Branch("minHFEnergy_34_52", &minHFEnergy_34_52, "minHFEnergy_34_52/F");

  tree->Branch("minHFEnergy_35_36", &minHFEnergy_35_36, "minHFEnergy_35_36/F");
  tree->Branch("minHFEnergy_35_37", &minHFEnergy_35_37, "minHFEnergy_35_37/F");
  tree->Branch("minHFEnergy_35_38", &minHFEnergy_35_38, "minHFEnergy_35_38/F");
  tree->Branch("minHFEnergy_35_39", &minHFEnergy_35_39, "minHFEnergy_35_39/F");
  tree->Branch("minHFEnergy_35_40", &minHFEnergy_35_40, "minHFEnergy_35_40/F");
  tree->Branch("minHFEnergy_35_41", &minHFEnergy_35_41, "minHFEnergy_35_41/F");
  tree->Branch("minHFEnergy_35_42", &minHFEnergy_35_42, "minHFEnergy_35_42/F");
  tree->Branch("minHFEnergy_35_43", &minHFEnergy_35_43, "minHFEnergy_35_43/F");
  tree->Branch("minHFEnergy_35_44", &minHFEnergy_35_44, "minHFEnergy_35_44/F");
  tree->Branch("minHFEnergy_35_45", &minHFEnergy_35_45, "minHFEnergy_35_45/F");
  tree->Branch("minHFEnergy_35_46", &minHFEnergy_35_46, "minHFEnergy_35_46/F");
  tree->Branch("minHFEnergy_35_47", &minHFEnergy_35_47, "minHFEnergy_35_47/F");
  tree->Branch("minHFEnergy_35_48", &minHFEnergy_35_48, "minHFEnergy_35_48/F");
  tree->Branch("minHFEnergy_35_49", &minHFEnergy_35_49, "minHFEnergy_35_49/F");
  tree->Branch("minHFEnergy_35_50", &minHFEnergy_35_50, "minHFEnergy_35_50/F");
  tree->Branch("minHFEnergy_35_51", &minHFEnergy_35_51, "minHFEnergy_35_51/F");
  tree->Branch("minHFEnergy_35_52", &minHFEnergy_35_52, "minHFEnergy_35_52/F");

  tree->Branch("minHFEnergy_36_37", &minHFEnergy_36_37, "minHFEnergy_36_37/F");
  tree->Branch("minHFEnergy_36_38", &minHFEnergy_36_38, "minHFEnergy_36_38/F");
  tree->Branch("minHFEnergy_36_39", &minHFEnergy_36_39, "minHFEnergy_36_39/F");
  tree->Branch("minHFEnergy_36_40", &minHFEnergy_36_40, "minHFEnergy_36_40/F");
  tree->Branch("minHFEnergy_36_41", &minHFEnergy_36_41, "minHFEnergy_36_41/F");
  tree->Branch("minHFEnergy_36_42", &minHFEnergy_36_42, "minHFEnergy_36_42/F");
  tree->Branch("minHFEnergy_36_43", &minHFEnergy_36_43, "minHFEnergy_36_43/F");
  tree->Branch("minHFEnergy_36_44", &minHFEnergy_36_44, "minHFEnergy_36_44/F");
  tree->Branch("minHFEnergy_36_45", &minHFEnergy_36_45, "minHFEnergy_36_45/F");
  tree->Branch("minHFEnergy_36_46", &minHFEnergy_36_46, "minHFEnergy_36_46/F");
  tree->Branch("minHFEnergy_36_47", &minHFEnergy_36_47, "minHFEnergy_36_47/F");
  tree->Branch("minHFEnergy_36_48", &minHFEnergy_36_48, "minHFEnergy_36_48/F");
  tree->Branch("minHFEnergy_36_49", &minHFEnergy_36_49, "minHFEnergy_36_49/F");
  tree->Branch("minHFEnergy_36_50", &minHFEnergy_36_50, "minHFEnergy_36_50/F");
  tree->Branch("minHFEnergy_36_51", &minHFEnergy_36_51, "minHFEnergy_36_51/F");
  tree->Branch("minHFEnergy_36_52", &minHFEnergy_36_52, "minHFEnergy_36_52/F");

  tree->Branch("minHFEnergy_37_38", &minHFEnergy_37_38, "minHFEnergy_37_38/F");
  tree->Branch("minHFEnergy_37_39", &minHFEnergy_37_39, "minHFEnergy_37_39/F");
  tree->Branch("minHFEnergy_37_40", &minHFEnergy_37_40, "minHFEnergy_37_40/F");
  tree->Branch("minHFEnergy_37_41", &minHFEnergy_37_41, "minHFEnergy_37_41/F");
  tree->Branch("minHFEnergy_37_42", &minHFEnergy_37_42, "minHFEnergy_37_42/F");
  tree->Branch("minHFEnergy_37_43", &minHFEnergy_37_43, "minHFEnergy_37_43/F");
  tree->Branch("minHFEnergy_37_44", &minHFEnergy_37_44, "minHFEnergy_37_44/F");
  tree->Branch("minHFEnergy_37_45", &minHFEnergy_37_45, "minHFEnergy_37_45/F");
  tree->Branch("minHFEnergy_37_46", &minHFEnergy_37_46, "minHFEnergy_37_46/F");
  tree->Branch("minHFEnergy_37_47", &minHFEnergy_37_47, "minHFEnergy_37_47/F");
  tree->Branch("minHFEnergy_37_48", &minHFEnergy_37_48, "minHFEnergy_37_48/F");
  tree->Branch("minHFEnergy_37_49", &minHFEnergy_37_49, "minHFEnergy_37_49/F");
  tree->Branch("minHFEnergy_37_50", &minHFEnergy_37_50, "minHFEnergy_37_50/F");
  tree->Branch("minHFEnergy_37_51", &minHFEnergy_37_51, "minHFEnergy_37_51/F");
  tree->Branch("minHFEnergy_37_52", &minHFEnergy_37_52, "minHFEnergy_37_52/F");

  tree->Branch("minHFEnergy_38_39", &minHFEnergy_38_39, "minHFEnergy_38_39/F");
  tree->Branch("minHFEnergy_38_40", &minHFEnergy_38_40, "minHFEnergy_38_40/F");
  tree->Branch("minHFEnergy_38_41", &minHFEnergy_38_41, "minHFEnergy_38_41/F");
  tree->Branch("minHFEnergy_38_42", &minHFEnergy_38_42, "minHFEnergy_38_42/F");
  tree->Branch("minHFEnergy_38_43", &minHFEnergy_38_43, "minHFEnergy_38_43/F");
  tree->Branch("minHFEnergy_38_44", &minHFEnergy_38_44, "minHFEnergy_38_44/F");
  tree->Branch("minHFEnergy_38_45", &minHFEnergy_38_45, "minHFEnergy_38_45/F");
  tree->Branch("minHFEnergy_38_46", &minHFEnergy_38_46, "minHFEnergy_38_46/F");
  tree->Branch("minHFEnergy_38_47", &minHFEnergy_38_47, "minHFEnergy_38_47/F");
  tree->Branch("minHFEnergy_38_48", &minHFEnergy_38_48, "minHFEnergy_38_48/F");
  tree->Branch("minHFEnergy_38_49", &minHFEnergy_38_49, "minHFEnergy_38_49/F");
  tree->Branch("minHFEnergy_38_50", &minHFEnergy_38_50, "minHFEnergy_38_50/F");
  tree->Branch("minHFEnergy_38_51", &minHFEnergy_38_51, "minHFEnergy_38_51/F");
  tree->Branch("minHFEnergy_38_52", &minHFEnergy_38_52, "minHFEnergy_38_52/F");

  tree->Branch("minHFEnergy_39_40", &minHFEnergy_39_40, "minHFEnergy_39_40/F");
  tree->Branch("minHFEnergy_39_41", &minHFEnergy_39_41, "minHFEnergy_39_41/F");
  tree->Branch("minHFEnergy_39_42", &minHFEnergy_39_42, "minHFEnergy_39_42/F");
  tree->Branch("minHFEnergy_39_43", &minHFEnergy_39_43, "minHFEnergy_39_43/F");
  tree->Branch("minHFEnergy_39_44", &minHFEnergy_39_44, "minHFEnergy_39_44/F");
  tree->Branch("minHFEnergy_39_45", &minHFEnergy_39_45, "minHFEnergy_39_45/F");
  tree->Branch("minHFEnergy_39_46", &minHFEnergy_39_46, "minHFEnergy_39_46/F");
  tree->Branch("minHFEnergy_39_47", &minHFEnergy_39_47, "minHFEnergy_39_47/F");
  tree->Branch("minHFEnergy_39_48", &minHFEnergy_39_48, "minHFEnergy_39_48/F");
  tree->Branch("minHFEnergy_39_49", &minHFEnergy_39_49, "minHFEnergy_39_49/F");
  tree->Branch("minHFEnergy_39_50", &minHFEnergy_39_50, "minHFEnergy_39_50/F");
  tree->Branch("minHFEnergy_39_51", &minHFEnergy_39_51, "minHFEnergy_39_51/F");
  tree->Branch("minHFEnergy_39_52", &minHFEnergy_39_52, "minHFEnergy_39_52/F");

  tree->Branch("minHFEnergy_40_41", &minHFEnergy_40_41, "minHFEnergy_40_41/F");
  tree->Branch("minHFEnergy_40_42", &minHFEnergy_40_42, "minHFEnergy_40_42/F");
  tree->Branch("minHFEnergy_40_43", &minHFEnergy_40_43, "minHFEnergy_40_43/F");
  tree->Branch("minHFEnergy_40_44", &minHFEnergy_40_44, "minHFEnergy_40_44/F");
  tree->Branch("minHFEnergy_40_45", &minHFEnergy_40_45, "minHFEnergy_40_45/F");
  tree->Branch("minHFEnergy_40_46", &minHFEnergy_40_46, "minHFEnergy_40_46/F");
  tree->Branch("minHFEnergy_40_47", &minHFEnergy_40_47, "minHFEnergy_40_47/F");
  tree->Branch("minHFEnergy_40_48", &minHFEnergy_40_48, "minHFEnergy_40_48/F");
  tree->Branch("minHFEnergy_40_49", &minHFEnergy_40_49, "minHFEnergy_40_49/F");
  tree->Branch("minHFEnergy_40_50", &minHFEnergy_40_50, "minHFEnergy_40_50/F");
  tree->Branch("minHFEnergy_40_51", &minHFEnergy_40_51, "minHFEnergy_40_51/F");
  tree->Branch("minHFEnergy_40_52", &minHFEnergy_40_52, "minHFEnergy_40_52/F");

  tree->Branch("minHFEnergy_41_42", &minHFEnergy_41_42, "minHFEnergy_41_42/F");
  tree->Branch("minHFEnergy_41_43", &minHFEnergy_41_43, "minHFEnergy_41_43/F");
  tree->Branch("minHFEnergy_41_44", &minHFEnergy_41_44, "minHFEnergy_41_44/F");
  tree->Branch("minHFEnergy_41_45", &minHFEnergy_41_45, "minHFEnergy_41_45/F");
  tree->Branch("minHFEnergy_41_46", &minHFEnergy_41_46, "minHFEnergy_41_46/F");
  tree->Branch("minHFEnergy_41_47", &minHFEnergy_41_47, "minHFEnergy_41_47/F");
  tree->Branch("minHFEnergy_41_48", &minHFEnergy_41_48, "minHFEnergy_41_48/F");
  tree->Branch("minHFEnergy_41_49", &minHFEnergy_41_49, "minHFEnergy_41_49/F");
  tree->Branch("minHFEnergy_41_50", &minHFEnergy_41_50, "minHFEnergy_41_50/F");
  tree->Branch("minHFEnergy_41_51", &minHFEnergy_41_51, "minHFEnergy_41_51/F");
  tree->Branch("minHFEnergy_41_52", &minHFEnergy_41_52, "minHFEnergy_41_52/F");

  tree->Branch("minHFEnergy_42_43", &minHFEnergy_42_43, "minHFEnergy_42_43/F");
  tree->Branch("minHFEnergy_42_44", &minHFEnergy_42_44, "minHFEnergy_42_44/F");
  tree->Branch("minHFEnergy_42_45", &minHFEnergy_42_45, "minHFEnergy_42_45/F");
  tree->Branch("minHFEnergy_42_46", &minHFEnergy_42_46, "minHFEnergy_42_46/F");
  tree->Branch("minHFEnergy_42_47", &minHFEnergy_42_47, "minHFEnergy_42_47/F");
  tree->Branch("minHFEnergy_42_48", &minHFEnergy_42_48, "minHFEnergy_42_48/F");
  tree->Branch("minHFEnergy_42_49", &minHFEnergy_42_49, "minHFEnergy_42_49/F");
  tree->Branch("minHFEnergy_42_50", &minHFEnergy_42_50, "minHFEnergy_42_50/F");
  tree->Branch("minHFEnergy_42_51", &minHFEnergy_42_51, "minHFEnergy_42_51/F");
  tree->Branch("minHFEnergy_42_52", &minHFEnergy_42_52, "minHFEnergy_42_52/F");

  tree->Branch("minHFEnergy_43_44", &minHFEnergy_43_44, "minHFEnergy_43_44/F");
  tree->Branch("minHFEnergy_43_45", &minHFEnergy_43_45, "minHFEnergy_43_45/F");
  tree->Branch("minHFEnergy_43_46", &minHFEnergy_43_46, "minHFEnergy_43_46/F");
  tree->Branch("minHFEnergy_43_47", &minHFEnergy_43_47, "minHFEnergy_43_47/F");
  tree->Branch("minHFEnergy_43_48", &minHFEnergy_43_48, "minHFEnergy_43_48/F");
  tree->Branch("minHFEnergy_43_49", &minHFEnergy_43_49, "minHFEnergy_43_49/F");
  tree->Branch("minHFEnergy_43_50", &minHFEnergy_43_50, "minHFEnergy_43_50/F");
  tree->Branch("minHFEnergy_43_51", &minHFEnergy_43_51, "minHFEnergy_43_51/F");
  tree->Branch("minHFEnergy_43_52", &minHFEnergy_43_52, "minHFEnergy_43_52/F");

  tree->Branch("minHFEnergy_44_45", &minHFEnergy_44_45, "minHFEnergy_44_45/F");
  tree->Branch("minHFEnergy_44_46", &minHFEnergy_44_46, "minHFEnergy_44_46/F");
  tree->Branch("minHFEnergy_44_47", &minHFEnergy_44_47, "minHFEnergy_44_47/F");
  tree->Branch("minHFEnergy_44_48", &minHFEnergy_44_48, "minHFEnergy_44_48/F");
  tree->Branch("minHFEnergy_44_49", &minHFEnergy_44_49, "minHFEnergy_44_49/F");
  tree->Branch("minHFEnergy_44_50", &minHFEnergy_44_50, "minHFEnergy_44_50/F");
  tree->Branch("minHFEnergy_44_51", &minHFEnergy_44_51, "minHFEnergy_44_51/F");
  tree->Branch("minHFEnergy_44_52", &minHFEnergy_44_52, "minHFEnergy_44_52/F");

  tree->Branch("minHFEnergy_45_46", &minHFEnergy_45_46, "minHFEnergy_45_46/F");
  tree->Branch("minHFEnergy_45_47", &minHFEnergy_45_47, "minHFEnergy_45_47/F");
  tree->Branch("minHFEnergy_45_48", &minHFEnergy_45_48, "minHFEnergy_45_48/F");
  tree->Branch("minHFEnergy_45_49", &minHFEnergy_45_49, "minHFEnergy_45_49/F");
  tree->Branch("minHFEnergy_45_50", &minHFEnergy_45_50, "minHFEnergy_45_50/F");
  tree->Branch("minHFEnergy_45_51", &minHFEnergy_45_51, "minHFEnergy_45_51/F");
  tree->Branch("minHFEnergy_45_52", &minHFEnergy_45_52, "minHFEnergy_45_52/F");

  tree->Branch("minHFEnergy_46_47", &minHFEnergy_46_47, "minHFEnergy_46_47/F");
  tree->Branch("minHFEnergy_46_48", &minHFEnergy_46_48, "minHFEnergy_46_48/F");
  tree->Branch("minHFEnergy_46_49", &minHFEnergy_46_49, "minHFEnergy_46_49/F");
  tree->Branch("minHFEnergy_46_50", &minHFEnergy_46_50, "minHFEnergy_46_50/F");
  tree->Branch("minHFEnergy_46_51", &minHFEnergy_46_51, "minHFEnergy_46_51/F");
  tree->Branch("minHFEnergy_46_52", &minHFEnergy_46_52, "minHFEnergy_46_52/F");

  tree->Branch("minHFEnergy_47_48", &minHFEnergy_47_48, "minHFEnergy_47_48/F");
  tree->Branch("minHFEnergy_47_49", &minHFEnergy_47_49, "minHFEnergy_47_49/F");
  tree->Branch("minHFEnergy_47_50", &minHFEnergy_47_50, "minHFEnergy_47_50/F");
  tree->Branch("minHFEnergy_47_51", &minHFEnergy_47_51, "minHFEnergy_47_51/F");
  tree->Branch("minHFEnergy_47_52", &minHFEnergy_47_52, "minHFEnergy_47_52/F");

  tree->Branch("minHFEnergy_48_49", &minHFEnergy_48_49, "minHFEnergy_48_49/F");
  tree->Branch("minHFEnergy_48_50", &minHFEnergy_48_50, "minHFEnergy_48_50/F");
  tree->Branch("minHFEnergy_48_51", &minHFEnergy_48_51, "minHFEnergy_48_51/F");
  tree->Branch("minHFEnergy_48_52", &minHFEnergy_48_52, "minHFEnergy_48_52/F");

  tree->Branch("minHFEnergy_49_50", &minHFEnergy_49_50, "minHFEnergy_49_50/F");
  tree->Branch("minHFEnergy_49_51", &minHFEnergy_49_51, "minHFEnergy_49_51/F");
  tree->Branch("minHFEnergy_49_52", &minHFEnergy_49_52, "minHFEnergy_49_52/F");

  tree->Branch("minHFEnergy_50_51", &minHFEnergy_50_51, "minHFEnergy_50_51/F");
  tree->Branch("minHFEnergy_50_52", &minHFEnergy_50_52, "minHFEnergy_50_52/F");

  tree->Branch("minHFEnergy_51_52", &minHFEnergy_51_52, "minHFEnergy_51_52/F");
  // minhfenergy --end

  // maxhfenergy
  tree->Branch("maxHFEnergy_30_31", &maxHFEnergy_30_31, "maxHFEnergy_30_31/F");
  tree->Branch("maxHFEnergy_30_32", &maxHFEnergy_30_32, "maxHFEnergy_30_32/F");
  tree->Branch("maxHFEnergy_30_33", &maxHFEnergy_30_33, "maxHFEnergy_30_33/F");
  tree->Branch("maxHFEnergy_30_34", &maxHFEnergy_30_34, "maxHFEnergy_30_34/F");
  tree->Branch("maxHFEnergy_30_35", &maxHFEnergy_30_35, "maxHFEnergy_30_35/F");
  tree->Branch("maxHFEnergy_30_36", &maxHFEnergy_30_36, "maxHFEnergy_30_36/F");
  tree->Branch("maxHFEnergy_30_37", &maxHFEnergy_30_37, "maxHFEnergy_30_37/F");
  tree->Branch("maxHFEnergy_30_38", &maxHFEnergy_30_38, "maxHFEnergy_30_38/F");
  tree->Branch("maxHFEnergy_30_39", &maxHFEnergy_30_39, "maxHFEnergy_30_39/F");
  tree->Branch("maxHFEnergy_30_40", &maxHFEnergy_30_40, "maxHFEnergy_30_40/F");
  tree->Branch("maxHFEnergy_30_41", &maxHFEnergy_30_41, "maxHFEnergy_30_41/F");
  tree->Branch("maxHFEnergy_30_42", &maxHFEnergy_30_42, "maxHFEnergy_30_42/F");
  tree->Branch("maxHFEnergy_30_43", &maxHFEnergy_30_43, "maxHFEnergy_30_43/F");
  tree->Branch("maxHFEnergy_30_44", &maxHFEnergy_30_44, "maxHFEnergy_30_44/F");
  tree->Branch("maxHFEnergy_30_45", &maxHFEnergy_30_45, "maxHFEnergy_30_45/F");
  tree->Branch("maxHFEnergy_30_46", &maxHFEnergy_30_46, "maxHFEnergy_30_46/F");
  tree->Branch("maxHFEnergy_30_47", &maxHFEnergy_30_47, "maxHFEnergy_30_47/F");
  tree->Branch("maxHFEnergy_30_48", &maxHFEnergy_30_48, "maxHFEnergy_30_48/F");
  tree->Branch("maxHFEnergy_30_49", &maxHFEnergy_30_49, "maxHFEnergy_30_49/F");
  tree->Branch("maxHFEnergy_30_50", &maxHFEnergy_30_50, "maxHFEnergy_30_50/F");
  tree->Branch("maxHFEnergy_30_51", &maxHFEnergy_30_51, "maxHFEnergy_30_51/F");
  tree->Branch("maxHFEnergy_30_52", &maxHFEnergy_30_52, "maxHFEnergy_30_52/F");

  tree->Branch("maxHFEnergy_31_32", &maxHFEnergy_31_32, "maxHFEnergy_31_32/F");
  tree->Branch("maxHFEnergy_31_33", &maxHFEnergy_31_33, "maxHFEnergy_31_33/F");
  tree->Branch("maxHFEnergy_31_34", &maxHFEnergy_31_34, "maxHFEnergy_31_34/F");
  tree->Branch("maxHFEnergy_31_35", &maxHFEnergy_31_35, "maxHFEnergy_31_35/F");
  tree->Branch("maxHFEnergy_31_36", &maxHFEnergy_31_36, "maxHFEnergy_31_36/F");
  tree->Branch("maxHFEnergy_31_37", &maxHFEnergy_31_37, "maxHFEnergy_31_37/F");
  tree->Branch("maxHFEnergy_31_38", &maxHFEnergy_31_38, "maxHFEnergy_31_38/F");
  tree->Branch("maxHFEnergy_31_39", &maxHFEnergy_31_39, "maxHFEnergy_31_39/F");
  tree->Branch("maxHFEnergy_31_40", &maxHFEnergy_31_40, "maxHFEnergy_31_40/F");
  tree->Branch("maxHFEnergy_31_41", &maxHFEnergy_31_41, "maxHFEnergy_31_41/F");
  tree->Branch("maxHFEnergy_31_42", &maxHFEnergy_31_42, "maxHFEnergy_31_42/F");
  tree->Branch("maxHFEnergy_31_43", &maxHFEnergy_31_43, "maxHFEnergy_31_43/F");
  tree->Branch("maxHFEnergy_31_44", &maxHFEnergy_31_44, "maxHFEnergy_31_44/F");
  tree->Branch("maxHFEnergy_31_45", &maxHFEnergy_31_45, "maxHFEnergy_31_45/F");
  tree->Branch("maxHFEnergy_31_46", &maxHFEnergy_31_46, "maxHFEnergy_31_46/F");
  tree->Branch("maxHFEnergy_31_47", &maxHFEnergy_31_47, "maxHFEnergy_31_47/F");
  tree->Branch("maxHFEnergy_31_48", &maxHFEnergy_31_48, "maxHFEnergy_31_48/F");
  tree->Branch("maxHFEnergy_31_49", &maxHFEnergy_31_49, "maxHFEnergy_31_49/F");
  tree->Branch("maxHFEnergy_31_50", &maxHFEnergy_31_50, "maxHFEnergy_31_50/F");
  tree->Branch("maxHFEnergy_31_51", &maxHFEnergy_31_51, "maxHFEnergy_31_51/F");
  tree->Branch("maxHFEnergy_31_52", &maxHFEnergy_31_52, "maxHFEnergy_31_52/F");

  tree->Branch("maxHFEnergy_32_33", &maxHFEnergy_32_33, "maxHFEnergy_32_33/F");
  tree->Branch("maxHFEnergy_32_34", &maxHFEnergy_32_34, "maxHFEnergy_32_34/F");
  tree->Branch("maxHFEnergy_32_35", &maxHFEnergy_32_35, "maxHFEnergy_32_35/F");
  tree->Branch("maxHFEnergy_32_36", &maxHFEnergy_32_36, "maxHFEnergy_32_36/F");
  tree->Branch("maxHFEnergy_32_37", &maxHFEnergy_32_37, "maxHFEnergy_32_37/F");
  tree->Branch("maxHFEnergy_32_38", &maxHFEnergy_32_38, "maxHFEnergy_32_38/F");
  tree->Branch("maxHFEnergy_32_39", &maxHFEnergy_32_39, "maxHFEnergy_32_39/F");
  tree->Branch("maxHFEnergy_32_40", &maxHFEnergy_32_40, "maxHFEnergy_32_40/F");
  tree->Branch("maxHFEnergy_32_41", &maxHFEnergy_32_41, "maxHFEnergy_32_41/F");
  tree->Branch("maxHFEnergy_32_42", &maxHFEnergy_32_42, "maxHFEnergy_32_42/F");
  tree->Branch("maxHFEnergy_32_43", &maxHFEnergy_32_43, "maxHFEnergy_32_43/F");
  tree->Branch("maxHFEnergy_32_44", &maxHFEnergy_32_44, "maxHFEnergy_32_44/F");
  tree->Branch("maxHFEnergy_32_45", &maxHFEnergy_32_45, "maxHFEnergy_32_45/F");
  tree->Branch("maxHFEnergy_32_46", &maxHFEnergy_32_46, "maxHFEnergy_32_46/F");
  tree->Branch("maxHFEnergy_32_47", &maxHFEnergy_32_47, "maxHFEnergy_32_47/F");
  tree->Branch("maxHFEnergy_32_48", &maxHFEnergy_32_48, "maxHFEnergy_32_48/F");
  tree->Branch("maxHFEnergy_32_49", &maxHFEnergy_32_49, "maxHFEnergy_32_49/F");
  tree->Branch("maxHFEnergy_32_50", &maxHFEnergy_32_50, "maxHFEnergy_32_50/F");
  tree->Branch("maxHFEnergy_32_51", &maxHFEnergy_32_51, "maxHFEnergy_32_51/F");
  tree->Branch("maxHFEnergy_32_52", &maxHFEnergy_32_52, "maxHFEnergy_32_52/F");

  tree->Branch("maxHFEnergy_33_34", &maxHFEnergy_33_34, "maxHFEnergy_33_34/F");
  tree->Branch("maxHFEnergy_33_35", &maxHFEnergy_33_35, "maxHFEnergy_33_35/F");
  tree->Branch("maxHFEnergy_33_36", &maxHFEnergy_33_36, "maxHFEnergy_33_36/F");
  tree->Branch("maxHFEnergy_33_37", &maxHFEnergy_33_37, "maxHFEnergy_33_37/F");
  tree->Branch("maxHFEnergy_33_38", &maxHFEnergy_33_38, "maxHFEnergy_33_38/F");
  tree->Branch("maxHFEnergy_33_39", &maxHFEnergy_33_39, "maxHFEnergy_33_39/F");
  tree->Branch("maxHFEnergy_33_40", &maxHFEnergy_33_40, "maxHFEnergy_33_40/F");
  tree->Branch("maxHFEnergy_33_41", &maxHFEnergy_33_41, "maxHFEnergy_33_41/F");
  tree->Branch("maxHFEnergy_33_42", &maxHFEnergy_33_42, "maxHFEnergy_33_42/F");
  tree->Branch("maxHFEnergy_33_43", &maxHFEnergy_33_43, "maxHFEnergy_33_43/F");
  tree->Branch("maxHFEnergy_33_44", &maxHFEnergy_33_44, "maxHFEnergy_33_44/F");
  tree->Branch("maxHFEnergy_33_45", &maxHFEnergy_33_45, "maxHFEnergy_33_45/F");
  tree->Branch("maxHFEnergy_33_46", &maxHFEnergy_33_46, "maxHFEnergy_33_46/F");
  tree->Branch("maxHFEnergy_33_47", &maxHFEnergy_33_47, "maxHFEnergy_33_47/F");
  tree->Branch("maxHFEnergy_33_48", &maxHFEnergy_33_48, "maxHFEnergy_33_48/F");
  tree->Branch("maxHFEnergy_33_49", &maxHFEnergy_33_49, "maxHFEnergy_33_49/F");
  tree->Branch("maxHFEnergy_33_50", &maxHFEnergy_33_50, "maxHFEnergy_33_50/F");
  tree->Branch("maxHFEnergy_33_51", &maxHFEnergy_33_51, "maxHFEnergy_33_51/F");
  tree->Branch("maxHFEnergy_33_52", &maxHFEnergy_33_52, "maxHFEnergy_33_52/F");

  tree->Branch("maxHFEnergy_34_35", &maxHFEnergy_34_35, "maxHFEnergy_34_35/F");
  tree->Branch("maxHFEnergy_34_36", &maxHFEnergy_34_36, "maxHFEnergy_34_36/F");
  tree->Branch("maxHFEnergy_34_37", &maxHFEnergy_34_37, "maxHFEnergy_34_37/F");
  tree->Branch("maxHFEnergy_34_38", &maxHFEnergy_34_38, "maxHFEnergy_34_38/F");
  tree->Branch("maxHFEnergy_34_39", &maxHFEnergy_34_39, "maxHFEnergy_34_39/F");
  tree->Branch("maxHFEnergy_34_40", &maxHFEnergy_34_40, "maxHFEnergy_34_40/F");
  tree->Branch("maxHFEnergy_34_41", &maxHFEnergy_34_41, "maxHFEnergy_34_41/F");
  tree->Branch("maxHFEnergy_34_42", &maxHFEnergy_34_42, "maxHFEnergy_34_42/F");
  tree->Branch("maxHFEnergy_34_43", &maxHFEnergy_34_43, "maxHFEnergy_34_43/F");
  tree->Branch("maxHFEnergy_34_44", &maxHFEnergy_34_44, "maxHFEnergy_34_44/F");
  tree->Branch("maxHFEnergy_34_45", &maxHFEnergy_34_45, "maxHFEnergy_34_45/F");
  tree->Branch("maxHFEnergy_34_46", &maxHFEnergy_34_46, "maxHFEnergy_34_46/F");
  tree->Branch("maxHFEnergy_34_47", &maxHFEnergy_34_47, "maxHFEnergy_34_47/F");
  tree->Branch("maxHFEnergy_34_48", &maxHFEnergy_34_48, "maxHFEnergy_34_48/F");
  tree->Branch("maxHFEnergy_34_49", &maxHFEnergy_34_49, "maxHFEnergy_34_49/F");
  tree->Branch("maxHFEnergy_34_50", &maxHFEnergy_34_50, "maxHFEnergy_34_50/F");
  tree->Branch("maxHFEnergy_34_51", &maxHFEnergy_34_51, "maxHFEnergy_34_51/F");
  tree->Branch("maxHFEnergy_34_52", &maxHFEnergy_34_52, "maxHFEnergy_34_52/F");

  tree->Branch("maxHFEnergy_35_36", &maxHFEnergy_35_36, "maxHFEnergy_35_36/F");
  tree->Branch("maxHFEnergy_35_37", &maxHFEnergy_35_37, "maxHFEnergy_35_37/F");
  tree->Branch("maxHFEnergy_35_38", &maxHFEnergy_35_38, "maxHFEnergy_35_38/F");
  tree->Branch("maxHFEnergy_35_39", &maxHFEnergy_35_39, "maxHFEnergy_35_39/F");
  tree->Branch("maxHFEnergy_35_40", &maxHFEnergy_35_40, "maxHFEnergy_35_40/F");
  tree->Branch("maxHFEnergy_35_41", &maxHFEnergy_35_41, "maxHFEnergy_35_41/F");
  tree->Branch("maxHFEnergy_35_42", &maxHFEnergy_35_42, "maxHFEnergy_35_42/F");
  tree->Branch("maxHFEnergy_35_43", &maxHFEnergy_35_43, "maxHFEnergy_35_43/F");
  tree->Branch("maxHFEnergy_35_44", &maxHFEnergy_35_44, "maxHFEnergy_35_44/F");
  tree->Branch("maxHFEnergy_35_45", &maxHFEnergy_35_45, "maxHFEnergy_35_45/F");
  tree->Branch("maxHFEnergy_35_46", &maxHFEnergy_35_46, "maxHFEnergy_35_46/F");
  tree->Branch("maxHFEnergy_35_47", &maxHFEnergy_35_47, "maxHFEnergy_35_47/F");
  tree->Branch("maxHFEnergy_35_48", &maxHFEnergy_35_48, "maxHFEnergy_35_48/F");
  tree->Branch("maxHFEnergy_35_49", &maxHFEnergy_35_49, "maxHFEnergy_35_49/F");
  tree->Branch("maxHFEnergy_35_50", &maxHFEnergy_35_50, "maxHFEnergy_35_50/F");
  tree->Branch("maxHFEnergy_35_51", &maxHFEnergy_35_51, "maxHFEnergy_35_51/F");
  tree->Branch("maxHFEnergy_35_52", &maxHFEnergy_35_52, "maxHFEnergy_35_52/F");

  tree->Branch("maxHFEnergy_36_37", &maxHFEnergy_36_37, "maxHFEnergy_36_37/F");
  tree->Branch("maxHFEnergy_36_38", &maxHFEnergy_36_38, "maxHFEnergy_36_38/F");
  tree->Branch("maxHFEnergy_36_39", &maxHFEnergy_36_39, "maxHFEnergy_36_39/F");
  tree->Branch("maxHFEnergy_36_40", &maxHFEnergy_36_40, "maxHFEnergy_36_40/F");
  tree->Branch("maxHFEnergy_36_41", &maxHFEnergy_36_41, "maxHFEnergy_36_41/F");
  tree->Branch("maxHFEnergy_36_42", &maxHFEnergy_36_42, "maxHFEnergy_36_42/F");
  tree->Branch("maxHFEnergy_36_43", &maxHFEnergy_36_43, "maxHFEnergy_36_43/F");
  tree->Branch("maxHFEnergy_36_44", &maxHFEnergy_36_44, "maxHFEnergy_36_44/F");
  tree->Branch("maxHFEnergy_36_45", &maxHFEnergy_36_45, "maxHFEnergy_36_45/F");
  tree->Branch("maxHFEnergy_36_46", &maxHFEnergy_36_46, "maxHFEnergy_36_46/F");
  tree->Branch("maxHFEnergy_36_47", &maxHFEnergy_36_47, "maxHFEnergy_36_47/F");
  tree->Branch("maxHFEnergy_36_48", &maxHFEnergy_36_48, "maxHFEnergy_36_48/F");
  tree->Branch("maxHFEnergy_36_49", &maxHFEnergy_36_49, "maxHFEnergy_36_49/F");
  tree->Branch("maxHFEnergy_36_50", &maxHFEnergy_36_50, "maxHFEnergy_36_50/F");
  tree->Branch("maxHFEnergy_36_51", &maxHFEnergy_36_51, "maxHFEnergy_36_51/F");
  tree->Branch("maxHFEnergy_36_52", &maxHFEnergy_36_52, "maxHFEnergy_36_52/F");

  tree->Branch("maxHFEnergy_37_38", &maxHFEnergy_37_38, "maxHFEnergy_37_38/F");
  tree->Branch("maxHFEnergy_37_39", &maxHFEnergy_37_39, "maxHFEnergy_37_39/F");
  tree->Branch("maxHFEnergy_37_40", &maxHFEnergy_37_40, "maxHFEnergy_37_40/F");
  tree->Branch("maxHFEnergy_37_41", &maxHFEnergy_37_41, "maxHFEnergy_37_41/F");
  tree->Branch("maxHFEnergy_37_42", &maxHFEnergy_37_42, "maxHFEnergy_37_42/F");
  tree->Branch("maxHFEnergy_37_43", &maxHFEnergy_37_43, "maxHFEnergy_37_43/F");
  tree->Branch("maxHFEnergy_37_44", &maxHFEnergy_37_44, "maxHFEnergy_37_44/F");
  tree->Branch("maxHFEnergy_37_45", &maxHFEnergy_37_45, "maxHFEnergy_37_45/F");
  tree->Branch("maxHFEnergy_37_46", &maxHFEnergy_37_46, "maxHFEnergy_37_46/F");
  tree->Branch("maxHFEnergy_37_47", &maxHFEnergy_37_47, "maxHFEnergy_37_47/F");
  tree->Branch("maxHFEnergy_37_48", &maxHFEnergy_37_48, "maxHFEnergy_37_48/F");
  tree->Branch("maxHFEnergy_37_49", &maxHFEnergy_37_49, "maxHFEnergy_37_49/F");
  tree->Branch("maxHFEnergy_37_50", &maxHFEnergy_37_50, "maxHFEnergy_37_50/F");
  tree->Branch("maxHFEnergy_37_51", &maxHFEnergy_37_51, "maxHFEnergy_37_51/F");
  tree->Branch("maxHFEnergy_37_52", &maxHFEnergy_37_52, "maxHFEnergy_37_52/F");

  tree->Branch("maxHFEnergy_38_39", &maxHFEnergy_38_39, "maxHFEnergy_38_39/F");
  tree->Branch("maxHFEnergy_38_40", &maxHFEnergy_38_40, "maxHFEnergy_38_40/F");
  tree->Branch("maxHFEnergy_38_41", &maxHFEnergy_38_41, "maxHFEnergy_38_41/F");
  tree->Branch("maxHFEnergy_38_42", &maxHFEnergy_38_42, "maxHFEnergy_38_42/F");
  tree->Branch("maxHFEnergy_38_43", &maxHFEnergy_38_43, "maxHFEnergy_38_43/F");
  tree->Branch("maxHFEnergy_38_44", &maxHFEnergy_38_44, "maxHFEnergy_38_44/F");
  tree->Branch("maxHFEnergy_38_45", &maxHFEnergy_38_45, "maxHFEnergy_38_45/F");
  tree->Branch("maxHFEnergy_38_46", &maxHFEnergy_38_46, "maxHFEnergy_38_46/F");
  tree->Branch("maxHFEnergy_38_47", &maxHFEnergy_38_47, "maxHFEnergy_38_47/F");
  tree->Branch("maxHFEnergy_38_48", &maxHFEnergy_38_48, "maxHFEnergy_38_48/F");
  tree->Branch("maxHFEnergy_38_49", &maxHFEnergy_38_49, "maxHFEnergy_38_49/F");
  tree->Branch("maxHFEnergy_38_50", &maxHFEnergy_38_50, "maxHFEnergy_38_50/F");
  tree->Branch("maxHFEnergy_38_51", &maxHFEnergy_38_51, "maxHFEnergy_38_51/F");
  tree->Branch("maxHFEnergy_38_52", &maxHFEnergy_38_52, "maxHFEnergy_38_52/F");

  tree->Branch("maxHFEnergy_39_40", &maxHFEnergy_39_40, "maxHFEnergy_39_40/F");
  tree->Branch("maxHFEnergy_39_41", &maxHFEnergy_39_41, "maxHFEnergy_39_41/F");
  tree->Branch("maxHFEnergy_39_42", &maxHFEnergy_39_42, "maxHFEnergy_39_42/F");
  tree->Branch("maxHFEnergy_39_43", &maxHFEnergy_39_43, "maxHFEnergy_39_43/F");
  tree->Branch("maxHFEnergy_39_44", &maxHFEnergy_39_44, "maxHFEnergy_39_44/F");
  tree->Branch("maxHFEnergy_39_45", &maxHFEnergy_39_45, "maxHFEnergy_39_45/F");
  tree->Branch("maxHFEnergy_39_46", &maxHFEnergy_39_46, "maxHFEnergy_39_46/F");
  tree->Branch("maxHFEnergy_39_47", &maxHFEnergy_39_47, "maxHFEnergy_39_47/F");
  tree->Branch("maxHFEnergy_39_48", &maxHFEnergy_39_48, "maxHFEnergy_39_48/F");
  tree->Branch("maxHFEnergy_39_49", &maxHFEnergy_39_49, "maxHFEnergy_39_49/F");
  tree->Branch("maxHFEnergy_39_50", &maxHFEnergy_39_50, "maxHFEnergy_39_50/F");
  tree->Branch("maxHFEnergy_39_51", &maxHFEnergy_39_51, "maxHFEnergy_39_51/F");
  tree->Branch("maxHFEnergy_39_52", &maxHFEnergy_39_52, "maxHFEnergy_39_52/F");

  tree->Branch("maxHFEnergy_40_41", &maxHFEnergy_40_41, "maxHFEnergy_40_41/F");
  tree->Branch("maxHFEnergy_40_42", &maxHFEnergy_40_42, "maxHFEnergy_40_42/F");
  tree->Branch("maxHFEnergy_40_43", &maxHFEnergy_40_43, "maxHFEnergy_40_43/F");
  tree->Branch("maxHFEnergy_40_44", &maxHFEnergy_40_44, "maxHFEnergy_40_44/F");
  tree->Branch("maxHFEnergy_40_45", &maxHFEnergy_40_45, "maxHFEnergy_40_45/F");
  tree->Branch("maxHFEnergy_40_46", &maxHFEnergy_40_46, "maxHFEnergy_40_46/F");
  tree->Branch("maxHFEnergy_40_47", &maxHFEnergy_40_47, "maxHFEnergy_40_47/F");
  tree->Branch("maxHFEnergy_40_48", &maxHFEnergy_40_48, "maxHFEnergy_40_48/F");
  tree->Branch("maxHFEnergy_40_49", &maxHFEnergy_40_49, "maxHFEnergy_40_49/F");
  tree->Branch("maxHFEnergy_40_50", &maxHFEnergy_40_50, "maxHFEnergy_40_50/F");
  tree->Branch("maxHFEnergy_40_51", &maxHFEnergy_40_51, "maxHFEnergy_40_51/F");
  tree->Branch("maxHFEnergy_40_52", &maxHFEnergy_40_52, "maxHFEnergy_40_52/F");

  tree->Branch("maxHFEnergy_41_42", &maxHFEnergy_41_42, "maxHFEnergy_41_42/F");
  tree->Branch("maxHFEnergy_41_43", &maxHFEnergy_41_43, "maxHFEnergy_41_43/F");
  tree->Branch("maxHFEnergy_41_44", &maxHFEnergy_41_44, "maxHFEnergy_41_44/F");
  tree->Branch("maxHFEnergy_41_45", &maxHFEnergy_41_45, "maxHFEnergy_41_45/F");
  tree->Branch("maxHFEnergy_41_46", &maxHFEnergy_41_46, "maxHFEnergy_41_46/F");
  tree->Branch("maxHFEnergy_41_47", &maxHFEnergy_41_47, "maxHFEnergy_41_47/F");
  tree->Branch("maxHFEnergy_41_48", &maxHFEnergy_41_48, "maxHFEnergy_41_48/F");
  tree->Branch("maxHFEnergy_41_49", &maxHFEnergy_41_49, "maxHFEnergy_41_49/F");
  tree->Branch("maxHFEnergy_41_50", &maxHFEnergy_41_50, "maxHFEnergy_41_50/F");
  tree->Branch("maxHFEnergy_41_51", &maxHFEnergy_41_51, "maxHFEnergy_41_51/F");
  tree->Branch("maxHFEnergy_41_52", &maxHFEnergy_41_52, "maxHFEnergy_41_52/F");

  tree->Branch("maxHFEnergy_42_43", &maxHFEnergy_42_43, "maxHFEnergy_42_43/F");
  tree->Branch("maxHFEnergy_42_44", &maxHFEnergy_42_44, "maxHFEnergy_42_44/F");
  tree->Branch("maxHFEnergy_42_45", &maxHFEnergy_42_45, "maxHFEnergy_42_45/F");
  tree->Branch("maxHFEnergy_42_46", &maxHFEnergy_42_46, "maxHFEnergy_42_46/F");
  tree->Branch("maxHFEnergy_42_47", &maxHFEnergy_42_47, "maxHFEnergy_42_47/F");
  tree->Branch("maxHFEnergy_42_48", &maxHFEnergy_42_48, "maxHFEnergy_42_48/F");
  tree->Branch("maxHFEnergy_42_49", &maxHFEnergy_42_49, "maxHFEnergy_42_49/F");
  tree->Branch("maxHFEnergy_42_50", &maxHFEnergy_42_50, "maxHFEnergy_42_50/F");
  tree->Branch("maxHFEnergy_42_51", &maxHFEnergy_42_51, "maxHFEnergy_42_51/F");
  tree->Branch("maxHFEnergy_42_52", &maxHFEnergy_42_52, "maxHFEnergy_42_52/F");

  tree->Branch("maxHFEnergy_43_44", &maxHFEnergy_43_44, "maxHFEnergy_43_44/F");
  tree->Branch("maxHFEnergy_43_45", &maxHFEnergy_43_45, "maxHFEnergy_43_45/F");
  tree->Branch("maxHFEnergy_43_46", &maxHFEnergy_43_46, "maxHFEnergy_43_46/F");
  tree->Branch("maxHFEnergy_43_47", &maxHFEnergy_43_47, "maxHFEnergy_43_47/F");
  tree->Branch("maxHFEnergy_43_48", &maxHFEnergy_43_48, "maxHFEnergy_43_48/F");
  tree->Branch("maxHFEnergy_43_49", &maxHFEnergy_43_49, "maxHFEnergy_43_49/F");
  tree->Branch("maxHFEnergy_43_50", &maxHFEnergy_43_50, "maxHFEnergy_43_50/F");
  tree->Branch("maxHFEnergy_43_51", &maxHFEnergy_43_51, "maxHFEnergy_43_51/F");
  tree->Branch("maxHFEnergy_43_52", &maxHFEnergy_43_52, "maxHFEnergy_43_52/F");

  tree->Branch("maxHFEnergy_44_45", &maxHFEnergy_44_45, "maxHFEnergy_44_45/F");
  tree->Branch("maxHFEnergy_44_46", &maxHFEnergy_44_46, "maxHFEnergy_44_46/F");
  tree->Branch("maxHFEnergy_44_47", &maxHFEnergy_44_47, "maxHFEnergy_44_47/F");
  tree->Branch("maxHFEnergy_44_48", &maxHFEnergy_44_48, "maxHFEnergy_44_48/F");
  tree->Branch("maxHFEnergy_44_49", &maxHFEnergy_44_49, "maxHFEnergy_44_49/F");
  tree->Branch("maxHFEnergy_44_50", &maxHFEnergy_44_50, "maxHFEnergy_44_50/F");
  tree->Branch("maxHFEnergy_44_51", &maxHFEnergy_44_51, "maxHFEnergy_44_51/F");
  tree->Branch("maxHFEnergy_44_52", &maxHFEnergy_44_52, "maxHFEnergy_44_52/F");

  tree->Branch("maxHFEnergy_45_46", &maxHFEnergy_45_46, "maxHFEnergy_45_46/F");
  tree->Branch("maxHFEnergy_45_47", &maxHFEnergy_45_47, "maxHFEnergy_45_47/F");
  tree->Branch("maxHFEnergy_45_48", &maxHFEnergy_45_48, "maxHFEnergy_45_48/F");
  tree->Branch("maxHFEnergy_45_49", &maxHFEnergy_45_49, "maxHFEnergy_45_49/F");
  tree->Branch("maxHFEnergy_45_50", &maxHFEnergy_45_50, "maxHFEnergy_45_50/F");
  tree->Branch("maxHFEnergy_45_51", &maxHFEnergy_45_51, "maxHFEnergy_45_51/F");
  tree->Branch("maxHFEnergy_45_52", &maxHFEnergy_45_52, "maxHFEnergy_45_52/F");

  tree->Branch("maxHFEnergy_46_47", &maxHFEnergy_46_47, "maxHFEnergy_46_47/F");
  tree->Branch("maxHFEnergy_46_48", &maxHFEnergy_46_48, "maxHFEnergy_46_48/F");
  tree->Branch("maxHFEnergy_46_49", &maxHFEnergy_46_49, "maxHFEnergy_46_49/F");
  tree->Branch("maxHFEnergy_46_50", &maxHFEnergy_46_50, "maxHFEnergy_46_50/F");
  tree->Branch("maxHFEnergy_46_51", &maxHFEnergy_46_51, "maxHFEnergy_46_51/F");
  tree->Branch("maxHFEnergy_46_52", &maxHFEnergy_46_52, "maxHFEnergy_46_52/F");

  tree->Branch("maxHFEnergy_47_48", &maxHFEnergy_47_48, "maxHFEnergy_47_48/F");
  tree->Branch("maxHFEnergy_47_49", &maxHFEnergy_47_49, "maxHFEnergy_47_49/F");
  tree->Branch("maxHFEnergy_47_50", &maxHFEnergy_47_50, "maxHFEnergy_47_50/F");
  tree->Branch("maxHFEnergy_47_51", &maxHFEnergy_47_51, "maxHFEnergy_47_51/F");
  tree->Branch("maxHFEnergy_47_52", &maxHFEnergy_47_52, "maxHFEnergy_47_52/F");

  tree->Branch("maxHFEnergy_48_49", &maxHFEnergy_48_49, "maxHFEnergy_48_49/F");
  tree->Branch("maxHFEnergy_48_50", &maxHFEnergy_48_50, "maxHFEnergy_48_50/F");
  tree->Branch("maxHFEnergy_48_51", &maxHFEnergy_48_51, "maxHFEnergy_48_51/F");
  tree->Branch("maxHFEnergy_48_52", &maxHFEnergy_48_52, "maxHFEnergy_48_52/F");

  tree->Branch("maxHFEnergy_49_50", &maxHFEnergy_49_50, "maxHFEnergy_49_50/F");
  tree->Branch("maxHFEnergy_49_51", &maxHFEnergy_49_51, "maxHFEnergy_49_51/F");
  tree->Branch("maxHFEnergy_49_52", &maxHFEnergy_49_52, "maxHFEnergy_49_52/F");

  tree->Branch("maxHFEnergy_50_51", &maxHFEnergy_50_51, "maxHFEnergy_50_51/F");
  tree->Branch("maxHFEnergy_50_52", &maxHFEnergy_50_52, "maxHFEnergy_50_52/F");

  tree->Branch("maxHFEnergy_51_52", &maxHFEnergy_51_52, "maxHFEnergy_51_52/F");
  // maxhfenergy

  
  // hfsumplus 
  tree->Branch("hfSumPlus_30_31", &hfSumPlus_30_31, "hfSumPlus_30_31/F");
  tree->Branch("hfSumPlus_30_32", &hfSumPlus_30_32, "hfSumPlus_30_32/F");
  tree->Branch("hfSumPlus_30_33", &hfSumPlus_30_33, "hfSumPlus_30_33/F");
  tree->Branch("hfSumPlus_30_34", &hfSumPlus_30_34, "hfSumPlus_30_34/F");
  tree->Branch("hfSumPlus_30_35", &hfSumPlus_30_35, "hfSumPlus_30_35/F");
  tree->Branch("hfSumPlus_30_36", &hfSumPlus_30_36, "hfSumPlus_30_36/F");
  tree->Branch("hfSumPlus_30_37", &hfSumPlus_30_37, "hfSumPlus_30_37/F");
  tree->Branch("hfSumPlus_30_38", &hfSumPlus_30_38, "hfSumPlus_30_38/F");
  tree->Branch("hfSumPlus_30_39", &hfSumPlus_30_39, "hfSumPlus_30_39/F");
  tree->Branch("hfSumPlus_30_40", &hfSumPlus_30_40, "hfSumPlus_30_40/F");
  tree->Branch("hfSumPlus_30_41", &hfSumPlus_30_41, "hfSumPlus_30_41/F");
  tree->Branch("hfSumPlus_30_42", &hfSumPlus_30_42, "hfSumPlus_30_42/F");
  tree->Branch("hfSumPlus_30_43", &hfSumPlus_30_43, "hfSumPlus_30_43/F");
  tree->Branch("hfSumPlus_30_44", &hfSumPlus_30_44, "hfSumPlus_30_44/F");
  tree->Branch("hfSumPlus_30_45", &hfSumPlus_30_45, "hfSumPlus_30_45/F");
  tree->Branch("hfSumPlus_30_46", &hfSumPlus_30_46, "hfSumPlus_30_46/F");
  tree->Branch("hfSumPlus_30_47", &hfSumPlus_30_47, "hfSumPlus_30_47/F");
  tree->Branch("hfSumPlus_30_48", &hfSumPlus_30_48, "hfSumPlus_30_48/F");
  tree->Branch("hfSumPlus_30_49", &hfSumPlus_30_49, "hfSumPlus_30_49/F");
  tree->Branch("hfSumPlus_30_50", &hfSumPlus_30_50, "hfSumPlus_30_50/F");
  tree->Branch("hfSumPlus_30_51", &hfSumPlus_30_51, "hfSumPlus_30_51/F");
  tree->Branch("hfSumPlus_30_52", &hfSumPlus_30_52, "hfSumPlus_30_52/F");

  tree->Branch("hfSumPlus_31_32", &hfSumPlus_31_32, "hfSumPlus_31_32/F");
  tree->Branch("hfSumPlus_31_33", &hfSumPlus_31_33, "hfSumPlus_31_33/F");
  tree->Branch("hfSumPlus_31_34", &hfSumPlus_31_34, "hfSumPlus_31_34/F");
  tree->Branch("hfSumPlus_31_35", &hfSumPlus_31_35, "hfSumPlus_31_35/F");
  tree->Branch("hfSumPlus_31_36", &hfSumPlus_31_36, "hfSumPlus_31_36/F");
  tree->Branch("hfSumPlus_31_37", &hfSumPlus_31_37, "hfSumPlus_31_37/F");
  tree->Branch("hfSumPlus_31_38", &hfSumPlus_31_38, "hfSumPlus_31_38/F");
  tree->Branch("hfSumPlus_31_39", &hfSumPlus_31_39, "hfSumPlus_31_39/F");
  tree->Branch("hfSumPlus_31_40", &hfSumPlus_31_40, "hfSumPlus_31_40/F");
  tree->Branch("hfSumPlus_31_41", &hfSumPlus_31_41, "hfSumPlus_31_41/F");
  tree->Branch("hfSumPlus_31_42", &hfSumPlus_31_42, "hfSumPlus_31_42/F");
  tree->Branch("hfSumPlus_31_43", &hfSumPlus_31_43, "hfSumPlus_31_43/F");
  tree->Branch("hfSumPlus_31_44", &hfSumPlus_31_44, "hfSumPlus_31_44/F");
  tree->Branch("hfSumPlus_31_45", &hfSumPlus_31_45, "hfSumPlus_31_45/F");
  tree->Branch("hfSumPlus_31_46", &hfSumPlus_31_46, "hfSumPlus_31_46/F");
  tree->Branch("hfSumPlus_31_47", &hfSumPlus_31_47, "hfSumPlus_31_47/F");
  tree->Branch("hfSumPlus_31_48", &hfSumPlus_31_48, "hfSumPlus_31_48/F");
  tree->Branch("hfSumPlus_31_49", &hfSumPlus_31_49, "hfSumPlus_31_49/F");
  tree->Branch("hfSumPlus_31_50", &hfSumPlus_31_50, "hfSumPlus_31_50/F");
  tree->Branch("hfSumPlus_31_51", &hfSumPlus_31_51, "hfSumPlus_31_51/F");
  tree->Branch("hfSumPlus_31_52", &hfSumPlus_31_52, "hfSumPlus_31_52/F");

  tree->Branch("hfSumPlus_32_33", &hfSumPlus_32_33, "hfSumPlus_32_33/F");
  tree->Branch("hfSumPlus_32_34", &hfSumPlus_32_34, "hfSumPlus_32_34/F");
  tree->Branch("hfSumPlus_32_35", &hfSumPlus_32_35, "hfSumPlus_32_35/F");
  tree->Branch("hfSumPlus_32_36", &hfSumPlus_32_36, "hfSumPlus_32_36/F");
  tree->Branch("hfSumPlus_32_37", &hfSumPlus_32_37, "hfSumPlus_32_37/F");
  tree->Branch("hfSumPlus_32_38", &hfSumPlus_32_38, "hfSumPlus_32_38/F");
  tree->Branch("hfSumPlus_32_39", &hfSumPlus_32_39, "hfSumPlus_32_39/F");
  tree->Branch("hfSumPlus_32_40", &hfSumPlus_32_40, "hfSumPlus_32_40/F");
  tree->Branch("hfSumPlus_32_41", &hfSumPlus_32_41, "hfSumPlus_32_41/F");
  tree->Branch("hfSumPlus_32_42", &hfSumPlus_32_42, "hfSumPlus_32_42/F");
  tree->Branch("hfSumPlus_32_43", &hfSumPlus_32_43, "hfSumPlus_32_43/F");
  tree->Branch("hfSumPlus_32_44", &hfSumPlus_32_44, "hfSumPlus_32_44/F");
  tree->Branch("hfSumPlus_32_45", &hfSumPlus_32_45, "hfSumPlus_32_45/F");
  tree->Branch("hfSumPlus_32_46", &hfSumPlus_32_46, "hfSumPlus_32_46/F");
  tree->Branch("hfSumPlus_32_47", &hfSumPlus_32_47, "hfSumPlus_32_47/F");
  tree->Branch("hfSumPlus_32_48", &hfSumPlus_32_48, "hfSumPlus_32_48/F");
  tree->Branch("hfSumPlus_32_49", &hfSumPlus_32_49, "hfSumPlus_32_49/F");
  tree->Branch("hfSumPlus_32_50", &hfSumPlus_32_50, "hfSumPlus_32_50/F");
  tree->Branch("hfSumPlus_32_51", &hfSumPlus_32_51, "hfSumPlus_32_51/F");
  tree->Branch("hfSumPlus_32_52", &hfSumPlus_32_52, "hfSumPlus_32_52/F");

  tree->Branch("hfSumPlus_33_34", &hfSumPlus_33_34, "hfSumPlus_33_34/F");
  tree->Branch("hfSumPlus_33_35", &hfSumPlus_33_35, "hfSumPlus_33_35/F");
  tree->Branch("hfSumPlus_33_36", &hfSumPlus_33_36, "hfSumPlus_33_36/F");
  tree->Branch("hfSumPlus_33_37", &hfSumPlus_33_37, "hfSumPlus_33_37/F");
  tree->Branch("hfSumPlus_33_38", &hfSumPlus_33_38, "hfSumPlus_33_38/F");
  tree->Branch("hfSumPlus_33_39", &hfSumPlus_33_39, "hfSumPlus_33_39/F");
  tree->Branch("hfSumPlus_33_40", &hfSumPlus_33_40, "hfSumPlus_33_40/F");
  tree->Branch("hfSumPlus_33_41", &hfSumPlus_33_41, "hfSumPlus_33_41/F");
  tree->Branch("hfSumPlus_33_42", &hfSumPlus_33_42, "hfSumPlus_33_42/F");
  tree->Branch("hfSumPlus_33_43", &hfSumPlus_33_43, "hfSumPlus_33_43/F");
  tree->Branch("hfSumPlus_33_44", &hfSumPlus_33_44, "hfSumPlus_33_44/F");
  tree->Branch("hfSumPlus_33_45", &hfSumPlus_33_45, "hfSumPlus_33_45/F");
  tree->Branch("hfSumPlus_33_46", &hfSumPlus_33_46, "hfSumPlus_33_46/F");
  tree->Branch("hfSumPlus_33_47", &hfSumPlus_33_47, "hfSumPlus_33_47/F");
  tree->Branch("hfSumPlus_33_48", &hfSumPlus_33_48, "hfSumPlus_33_48/F");
  tree->Branch("hfSumPlus_33_49", &hfSumPlus_33_49, "hfSumPlus_33_49/F");
  tree->Branch("hfSumPlus_33_50", &hfSumPlus_33_50, "hfSumPlus_33_50/F");
  tree->Branch("hfSumPlus_33_51", &hfSumPlus_33_51, "hfSumPlus_33_51/F");
  tree->Branch("hfSumPlus_33_52", &hfSumPlus_33_52, "hfSumPlus_33_52/F");

  tree->Branch("hfSumPlus_34_35", &hfSumPlus_34_35, "hfSumPlus_34_35/F");
  tree->Branch("hfSumPlus_34_36", &hfSumPlus_34_36, "hfSumPlus_34_36/F");
  tree->Branch("hfSumPlus_34_37", &hfSumPlus_34_37, "hfSumPlus_34_37/F");
  tree->Branch("hfSumPlus_34_38", &hfSumPlus_34_38, "hfSumPlus_34_38/F");
  tree->Branch("hfSumPlus_34_39", &hfSumPlus_34_39, "hfSumPlus_34_39/F");
  tree->Branch("hfSumPlus_34_40", &hfSumPlus_34_40, "hfSumPlus_34_40/F");
  tree->Branch("hfSumPlus_34_41", &hfSumPlus_34_41, "hfSumPlus_34_41/F");
  tree->Branch("hfSumPlus_34_42", &hfSumPlus_34_42, "hfSumPlus_34_42/F");
  tree->Branch("hfSumPlus_34_43", &hfSumPlus_34_43, "hfSumPlus_34_43/F");
  tree->Branch("hfSumPlus_34_44", &hfSumPlus_34_44, "hfSumPlus_34_44/F");
  tree->Branch("hfSumPlus_34_45", &hfSumPlus_34_45, "hfSumPlus_34_45/F");
  tree->Branch("hfSumPlus_34_46", &hfSumPlus_34_46, "hfSumPlus_34_46/F");
  tree->Branch("hfSumPlus_34_47", &hfSumPlus_34_47, "hfSumPlus_34_47/F");
  tree->Branch("hfSumPlus_34_48", &hfSumPlus_34_48, "hfSumPlus_34_48/F");
  tree->Branch("hfSumPlus_34_49", &hfSumPlus_34_49, "hfSumPlus_34_49/F");
  tree->Branch("hfSumPlus_34_50", &hfSumPlus_34_50, "hfSumPlus_34_50/F");
  tree->Branch("hfSumPlus_34_51", &hfSumPlus_34_51, "hfSumPlus_34_51/F");
  tree->Branch("hfSumPlus_34_52", &hfSumPlus_34_52, "hfSumPlus_34_52/F");

  tree->Branch("hfSumPlus_35_36", &hfSumPlus_35_36, "hfSumPlus_35_36/F");
  tree->Branch("hfSumPlus_35_37", &hfSumPlus_35_37, "hfSumPlus_35_37/F");
  tree->Branch("hfSumPlus_35_38", &hfSumPlus_35_38, "hfSumPlus_35_38/F");
  tree->Branch("hfSumPlus_35_39", &hfSumPlus_35_39, "hfSumPlus_35_39/F");
  tree->Branch("hfSumPlus_35_40", &hfSumPlus_35_40, "hfSumPlus_35_40/F");
  tree->Branch("hfSumPlus_35_41", &hfSumPlus_35_41, "hfSumPlus_35_41/F");
  tree->Branch("hfSumPlus_35_42", &hfSumPlus_35_42, "hfSumPlus_35_42/F");
  tree->Branch("hfSumPlus_35_43", &hfSumPlus_35_43, "hfSumPlus_35_43/F");
  tree->Branch("hfSumPlus_35_44", &hfSumPlus_35_44, "hfSumPlus_35_44/F");
  tree->Branch("hfSumPlus_35_45", &hfSumPlus_35_45, "hfSumPlus_35_45/F");
  tree->Branch("hfSumPlus_35_46", &hfSumPlus_35_46, "hfSumPlus_35_46/F");
  tree->Branch("hfSumPlus_35_47", &hfSumPlus_35_47, "hfSumPlus_35_47/F");
  tree->Branch("hfSumPlus_35_48", &hfSumPlus_35_48, "hfSumPlus_35_48/F");
  tree->Branch("hfSumPlus_35_49", &hfSumPlus_35_49, "hfSumPlus_35_49/F");
  tree->Branch("hfSumPlus_35_50", &hfSumPlus_35_50, "hfSumPlus_35_50/F");
  tree->Branch("hfSumPlus_35_51", &hfSumPlus_35_51, "hfSumPlus_35_51/F");
  tree->Branch("hfSumPlus_35_52", &hfSumPlus_35_52, "hfSumPlus_35_52/F");

  tree->Branch("hfSumPlus_36_37", &hfSumPlus_36_37, "hfSumPlus_36_37/F");
  tree->Branch("hfSumPlus_36_38", &hfSumPlus_36_38, "hfSumPlus_36_38/F");
  tree->Branch("hfSumPlus_36_39", &hfSumPlus_36_39, "hfSumPlus_36_39/F");
  tree->Branch("hfSumPlus_36_40", &hfSumPlus_36_40, "hfSumPlus_36_40/F");
  tree->Branch("hfSumPlus_36_41", &hfSumPlus_36_41, "hfSumPlus_36_41/F");
  tree->Branch("hfSumPlus_36_42", &hfSumPlus_36_42, "hfSumPlus_36_42/F");
  tree->Branch("hfSumPlus_36_43", &hfSumPlus_36_43, "hfSumPlus_36_43/F");
  tree->Branch("hfSumPlus_36_44", &hfSumPlus_36_44, "hfSumPlus_36_44/F");
  tree->Branch("hfSumPlus_36_45", &hfSumPlus_36_45, "hfSumPlus_36_45/F");
  tree->Branch("hfSumPlus_36_46", &hfSumPlus_36_46, "hfSumPlus_36_46/F");
  tree->Branch("hfSumPlus_36_47", &hfSumPlus_36_47, "hfSumPlus_36_47/F");
  tree->Branch("hfSumPlus_36_48", &hfSumPlus_36_48, "hfSumPlus_36_48/F");
  tree->Branch("hfSumPlus_36_49", &hfSumPlus_36_49, "hfSumPlus_36_49/F");
  tree->Branch("hfSumPlus_36_50", &hfSumPlus_36_50, "hfSumPlus_36_50/F");
  tree->Branch("hfSumPlus_36_51", &hfSumPlus_36_51, "hfSumPlus_36_51/F");
  tree->Branch("hfSumPlus_36_52", &hfSumPlus_36_52, "hfSumPlus_36_52/F");

  tree->Branch("hfSumPlus_37_38", &hfSumPlus_37_38, "hfSumPlus_37_38/F");
  tree->Branch("hfSumPlus_37_39", &hfSumPlus_37_39, "hfSumPlus_37_39/F");
  tree->Branch("hfSumPlus_37_40", &hfSumPlus_37_40, "hfSumPlus_37_40/F");
  tree->Branch("hfSumPlus_37_41", &hfSumPlus_37_41, "hfSumPlus_37_41/F");
  tree->Branch("hfSumPlus_37_42", &hfSumPlus_37_42, "hfSumPlus_37_42/F");
  tree->Branch("hfSumPlus_37_43", &hfSumPlus_37_43, "hfSumPlus_37_43/F");
  tree->Branch("hfSumPlus_37_44", &hfSumPlus_37_44, "hfSumPlus_37_44/F");
  tree->Branch("hfSumPlus_37_45", &hfSumPlus_37_45, "hfSumPlus_37_45/F");
  tree->Branch("hfSumPlus_37_46", &hfSumPlus_37_46, "hfSumPlus_37_46/F");
  tree->Branch("hfSumPlus_37_47", &hfSumPlus_37_47, "hfSumPlus_37_47/F");
  tree->Branch("hfSumPlus_37_48", &hfSumPlus_37_48, "hfSumPlus_37_48/F");
  tree->Branch("hfSumPlus_37_49", &hfSumPlus_37_49, "hfSumPlus_37_49/F");
  tree->Branch("hfSumPlus_37_50", &hfSumPlus_37_50, "hfSumPlus_37_50/F");
  tree->Branch("hfSumPlus_37_51", &hfSumPlus_37_51, "hfSumPlus_37_51/F");
  tree->Branch("hfSumPlus_37_52", &hfSumPlus_37_52, "hfSumPlus_37_52/F");

  tree->Branch("hfSumPlus_38_39", &hfSumPlus_38_39, "hfSumPlus_38_39/F");
  tree->Branch("hfSumPlus_38_40", &hfSumPlus_38_40, "hfSumPlus_38_40/F");
  tree->Branch("hfSumPlus_38_41", &hfSumPlus_38_41, "hfSumPlus_38_41/F");
  tree->Branch("hfSumPlus_38_42", &hfSumPlus_38_42, "hfSumPlus_38_42/F");
  tree->Branch("hfSumPlus_38_43", &hfSumPlus_38_43, "hfSumPlus_38_43/F");
  tree->Branch("hfSumPlus_38_44", &hfSumPlus_38_44, "hfSumPlus_38_44/F");
  tree->Branch("hfSumPlus_38_45", &hfSumPlus_38_45, "hfSumPlus_38_45/F");
  tree->Branch("hfSumPlus_38_46", &hfSumPlus_38_46, "hfSumPlus_38_46/F");
  tree->Branch("hfSumPlus_38_47", &hfSumPlus_38_47, "hfSumPlus_38_47/F");
  tree->Branch("hfSumPlus_38_48", &hfSumPlus_38_48, "hfSumPlus_38_48/F");
  tree->Branch("hfSumPlus_38_49", &hfSumPlus_38_49, "hfSumPlus_38_49/F");
  tree->Branch("hfSumPlus_38_50", &hfSumPlus_38_50, "hfSumPlus_38_50/F");
  tree->Branch("hfSumPlus_38_51", &hfSumPlus_38_51, "hfSumPlus_38_51/F");
  tree->Branch("hfSumPlus_38_52", &hfSumPlus_38_52, "hfSumPlus_38_52/F");

  tree->Branch("hfSumPlus_39_40", &hfSumPlus_39_40, "hfSumPlus_39_40/F");
  tree->Branch("hfSumPlus_39_41", &hfSumPlus_39_41, "hfSumPlus_39_41/F");
  tree->Branch("hfSumPlus_39_42", &hfSumPlus_39_42, "hfSumPlus_39_42/F");
  tree->Branch("hfSumPlus_39_43", &hfSumPlus_39_43, "hfSumPlus_39_43/F");
  tree->Branch("hfSumPlus_39_44", &hfSumPlus_39_44, "hfSumPlus_39_44/F");
  tree->Branch("hfSumPlus_39_45", &hfSumPlus_39_45, "hfSumPlus_39_45/F");
  tree->Branch("hfSumPlus_39_46", &hfSumPlus_39_46, "hfSumPlus_39_46/F");
  tree->Branch("hfSumPlus_39_47", &hfSumPlus_39_47, "hfSumPlus_39_47/F");
  tree->Branch("hfSumPlus_39_48", &hfSumPlus_39_48, "hfSumPlus_39_48/F");
  tree->Branch("hfSumPlus_39_49", &hfSumPlus_39_49, "hfSumPlus_39_49/F");
  tree->Branch("hfSumPlus_39_50", &hfSumPlus_39_50, "hfSumPlus_39_50/F");
  tree->Branch("hfSumPlus_39_51", &hfSumPlus_39_51, "hfSumPlus_39_51/F");
  tree->Branch("hfSumPlus_39_52", &hfSumPlus_39_52, "hfSumPlus_39_52/F");

  tree->Branch("hfSumPlus_40_41", &hfSumPlus_40_41, "hfSumPlus_40_41/F");
  tree->Branch("hfSumPlus_40_42", &hfSumPlus_40_42, "hfSumPlus_40_42/F");
  tree->Branch("hfSumPlus_40_43", &hfSumPlus_40_43, "hfSumPlus_40_43/F");
  tree->Branch("hfSumPlus_40_44", &hfSumPlus_40_44, "hfSumPlus_40_44/F");
  tree->Branch("hfSumPlus_40_45", &hfSumPlus_40_45, "hfSumPlus_40_45/F");
  tree->Branch("hfSumPlus_40_46", &hfSumPlus_40_46, "hfSumPlus_40_46/F");
  tree->Branch("hfSumPlus_40_47", &hfSumPlus_40_47, "hfSumPlus_40_47/F");
  tree->Branch("hfSumPlus_40_48", &hfSumPlus_40_48, "hfSumPlus_40_48/F");
  tree->Branch("hfSumPlus_40_49", &hfSumPlus_40_49, "hfSumPlus_40_49/F");
  tree->Branch("hfSumPlus_40_50", &hfSumPlus_40_50, "hfSumPlus_40_50/F");
  tree->Branch("hfSumPlus_40_51", &hfSumPlus_40_51, "hfSumPlus_40_51/F");
  tree->Branch("hfSumPlus_40_52", &hfSumPlus_40_52, "hfSumPlus_40_52/F");

  tree->Branch("hfSumPlus_41_42", &hfSumPlus_41_42, "hfSumPlus_41_42/F");
  tree->Branch("hfSumPlus_41_43", &hfSumPlus_41_43, "hfSumPlus_41_43/F");
  tree->Branch("hfSumPlus_41_44", &hfSumPlus_41_44, "hfSumPlus_41_44/F");
  tree->Branch("hfSumPlus_41_45", &hfSumPlus_41_45, "hfSumPlus_41_45/F");
  tree->Branch("hfSumPlus_41_46", &hfSumPlus_41_46, "hfSumPlus_41_46/F");
  tree->Branch("hfSumPlus_41_47", &hfSumPlus_41_47, "hfSumPlus_41_47/F");
  tree->Branch("hfSumPlus_41_48", &hfSumPlus_41_48, "hfSumPlus_41_48/F");
  tree->Branch("hfSumPlus_41_49", &hfSumPlus_41_49, "hfSumPlus_41_49/F");
  tree->Branch("hfSumPlus_41_50", &hfSumPlus_41_50, "hfSumPlus_41_50/F");
  tree->Branch("hfSumPlus_41_51", &hfSumPlus_41_51, "hfSumPlus_41_51/F");
  tree->Branch("hfSumPlus_41_52", &hfSumPlus_41_52, "hfSumPlus_41_52/F");

  tree->Branch("hfSumPlus_42_43", &hfSumPlus_42_43, "hfSumPlus_42_43/F");
  tree->Branch("hfSumPlus_42_44", &hfSumPlus_42_44, "hfSumPlus_42_44/F");
  tree->Branch("hfSumPlus_42_45", &hfSumPlus_42_45, "hfSumPlus_42_45/F");
  tree->Branch("hfSumPlus_42_46", &hfSumPlus_42_46, "hfSumPlus_42_46/F");
  tree->Branch("hfSumPlus_42_47", &hfSumPlus_42_47, "hfSumPlus_42_47/F");
  tree->Branch("hfSumPlus_42_48", &hfSumPlus_42_48, "hfSumPlus_42_48/F");
  tree->Branch("hfSumPlus_42_49", &hfSumPlus_42_49, "hfSumPlus_42_49/F");
  tree->Branch("hfSumPlus_42_50", &hfSumPlus_42_50, "hfSumPlus_42_50/F");
  tree->Branch("hfSumPlus_42_51", &hfSumPlus_42_51, "hfSumPlus_42_51/F");
  tree->Branch("hfSumPlus_42_52", &hfSumPlus_42_52, "hfSumPlus_42_52/F");

  tree->Branch("hfSumPlus_43_44", &hfSumPlus_43_44, "hfSumPlus_43_44/F");
  tree->Branch("hfSumPlus_43_45", &hfSumPlus_43_45, "hfSumPlus_43_45/F");
  tree->Branch("hfSumPlus_43_46", &hfSumPlus_43_46, "hfSumPlus_43_46/F");
  tree->Branch("hfSumPlus_43_47", &hfSumPlus_43_47, "hfSumPlus_43_47/F");
  tree->Branch("hfSumPlus_43_48", &hfSumPlus_43_48, "hfSumPlus_43_48/F");
  tree->Branch("hfSumPlus_43_49", &hfSumPlus_43_49, "hfSumPlus_43_49/F");
  tree->Branch("hfSumPlus_43_50", &hfSumPlus_43_50, "hfSumPlus_43_50/F");
  tree->Branch("hfSumPlus_43_51", &hfSumPlus_43_51, "hfSumPlus_43_51/F");
  tree->Branch("hfSumPlus_43_52", &hfSumPlus_43_52, "hfSumPlus_43_52/F");

  tree->Branch("hfSumPlus_44_45", &hfSumPlus_44_45, "hfSumPlus_44_45/F");
  tree->Branch("hfSumPlus_44_46", &hfSumPlus_44_46, "hfSumPlus_44_46/F");
  tree->Branch("hfSumPlus_44_47", &hfSumPlus_44_47, "hfSumPlus_44_47/F");
  tree->Branch("hfSumPlus_44_48", &hfSumPlus_44_48, "hfSumPlus_44_48/F");
  tree->Branch("hfSumPlus_44_49", &hfSumPlus_44_49, "hfSumPlus_44_49/F");
  tree->Branch("hfSumPlus_44_50", &hfSumPlus_44_50, "hfSumPlus_44_50/F");
  tree->Branch("hfSumPlus_44_51", &hfSumPlus_44_51, "hfSumPlus_44_51/F");
  tree->Branch("hfSumPlus_44_52", &hfSumPlus_44_52, "hfSumPlus_44_52/F");

  tree->Branch("hfSumPlus_45_46", &hfSumPlus_45_46, "hfSumPlus_45_46/F");
  tree->Branch("hfSumPlus_45_47", &hfSumPlus_45_47, "hfSumPlus_45_47/F");
  tree->Branch("hfSumPlus_45_48", &hfSumPlus_45_48, "hfSumPlus_45_48/F");
  tree->Branch("hfSumPlus_45_49", &hfSumPlus_45_49, "hfSumPlus_45_49/F");
  tree->Branch("hfSumPlus_45_50", &hfSumPlus_45_50, "hfSumPlus_45_50/F");
  tree->Branch("hfSumPlus_45_51", &hfSumPlus_45_51, "hfSumPlus_45_51/F");
  tree->Branch("hfSumPlus_45_52", &hfSumPlus_45_52, "hfSumPlus_45_52/F");

  tree->Branch("hfSumPlus_46_47", &hfSumPlus_46_47, "hfSumPlus_46_47/F");
  tree->Branch("hfSumPlus_46_48", &hfSumPlus_46_48, "hfSumPlus_46_48/F");
  tree->Branch("hfSumPlus_46_49", &hfSumPlus_46_49, "hfSumPlus_46_49/F");
  tree->Branch("hfSumPlus_46_50", &hfSumPlus_46_50, "hfSumPlus_46_50/F");
  tree->Branch("hfSumPlus_46_51", &hfSumPlus_46_51, "hfSumPlus_46_51/F");
  tree->Branch("hfSumPlus_46_52", &hfSumPlus_46_52, "hfSumPlus_46_52/F");

  tree->Branch("hfSumPlus_47_48", &hfSumPlus_47_48, "hfSumPlus_47_48/F");
  tree->Branch("hfSumPlus_47_49", &hfSumPlus_47_49, "hfSumPlus_47_49/F");
  tree->Branch("hfSumPlus_47_50", &hfSumPlus_47_50, "hfSumPlus_47_50/F");
  tree->Branch("hfSumPlus_47_51", &hfSumPlus_47_51, "hfSumPlus_47_51/F");
  tree->Branch("hfSumPlus_47_52", &hfSumPlus_47_52, "hfSumPlus_47_52/F");

  tree->Branch("hfSumPlus_48_49", &hfSumPlus_48_49, "hfSumPlus_48_49/F");
  tree->Branch("hfSumPlus_48_50", &hfSumPlus_48_50, "hfSumPlus_48_50/F");
  tree->Branch("hfSumPlus_48_51", &hfSumPlus_48_51, "hfSumPlus_48_51/F");
  tree->Branch("hfSumPlus_48_52", &hfSumPlus_48_52, "hfSumPlus_48_52/F");

  tree->Branch("hfSumPlus_49_50", &hfSumPlus_49_50, "hfSumPlus_49_50/F");
  tree->Branch("hfSumPlus_49_51", &hfSumPlus_49_51, "hfSumPlus_49_51/F");
  tree->Branch("hfSumPlus_49_52", &hfSumPlus_49_52, "hfSumPlus_49_52/F");

  tree->Branch("hfSumPlus_50_51", &hfSumPlus_50_51, "hfSumPlus_50_51/F");
  tree->Branch("hfSumPlus_50_52", &hfSumPlus_50_52, "hfSumPlus_50_52/F");

  tree->Branch("hfSumPlus_51_52", &hfSumPlus_51_52, "hfSumPlus_51_52/F");

  // hfsumplus --end
  
  // hfsumminus
  tree->Branch("hfSumMinus_30_31", &hfSumMinus_30_31, "hfSumMinus_30_31/F");
  tree->Branch("hfSumMinus_30_32", &hfSumMinus_30_32, "hfSumMinus_30_32/F");
  tree->Branch("hfSumMinus_30_33", &hfSumMinus_30_33, "hfSumMinus_30_33/F");
  tree->Branch("hfSumMinus_30_34", &hfSumMinus_30_34, "hfSumMinus_30_34/F");
  tree->Branch("hfSumMinus_30_35", &hfSumMinus_30_35, "hfSumMinus_30_35/F");
  tree->Branch("hfSumMinus_30_36", &hfSumMinus_30_36, "hfSumMinus_30_36/F");
  tree->Branch("hfSumMinus_30_37", &hfSumMinus_30_37, "hfSumMinus_30_37/F");
  tree->Branch("hfSumMinus_30_38", &hfSumMinus_30_38, "hfSumMinus_30_38/F");
  tree->Branch("hfSumMinus_30_39", &hfSumMinus_30_39, "hfSumMinus_30_39/F");
  tree->Branch("hfSumMinus_30_40", &hfSumMinus_30_40, "hfSumMinus_30_40/F");
  tree->Branch("hfSumMinus_30_41", &hfSumMinus_30_41, "hfSumMinus_30_41/F");
  tree->Branch("hfSumMinus_30_42", &hfSumMinus_30_42, "hfSumMinus_30_42/F");
  tree->Branch("hfSumMinus_30_43", &hfSumMinus_30_43, "hfSumMinus_30_43/F");
  tree->Branch("hfSumMinus_30_44", &hfSumMinus_30_44, "hfSumMinus_30_44/F");
  tree->Branch("hfSumMinus_30_45", &hfSumMinus_30_45, "hfSumMinus_30_45/F");
  tree->Branch("hfSumMinus_30_46", &hfSumMinus_30_46, "hfSumMinus_30_46/F");
  tree->Branch("hfSumMinus_30_47", &hfSumMinus_30_47, "hfSumMinus_30_47/F");
  tree->Branch("hfSumMinus_30_48", &hfSumMinus_30_48, "hfSumMinus_30_48/F");
  tree->Branch("hfSumMinus_30_49", &hfSumMinus_30_49, "hfSumMinus_30_49/F");
  tree->Branch("hfSumMinus_30_50", &hfSumMinus_30_50, "hfSumMinus_30_50/F");
  tree->Branch("hfSumMinus_30_51", &hfSumMinus_30_51, "hfSumMinus_30_51/F");
  tree->Branch("hfSumMinus_30_52", &hfSumMinus_30_52, "hfSumMinus_30_52/F");

  tree->Branch("hfSumMinus_31_32", &hfSumMinus_31_32, "hfSumMinus_31_32/F");
  tree->Branch("hfSumMinus_31_33", &hfSumMinus_31_33, "hfSumMinus_31_33/F");
  tree->Branch("hfSumMinus_31_34", &hfSumMinus_31_34, "hfSumMinus_31_34/F");
  tree->Branch("hfSumMinus_31_35", &hfSumMinus_31_35, "hfSumMinus_31_35/F");
  tree->Branch("hfSumMinus_31_36", &hfSumMinus_31_36, "hfSumMinus_31_36/F");
  tree->Branch("hfSumMinus_31_37", &hfSumMinus_31_37, "hfSumMinus_31_37/F");
  tree->Branch("hfSumMinus_31_38", &hfSumMinus_31_38, "hfSumMinus_31_38/F");
  tree->Branch("hfSumMinus_31_39", &hfSumMinus_31_39, "hfSumMinus_31_39/F");
  tree->Branch("hfSumMinus_31_40", &hfSumMinus_31_40, "hfSumMinus_31_40/F");
  tree->Branch("hfSumMinus_31_41", &hfSumMinus_31_41, "hfSumMinus_31_41/F");
  tree->Branch("hfSumMinus_31_42", &hfSumMinus_31_42, "hfSumMinus_31_42/F");
  tree->Branch("hfSumMinus_31_43", &hfSumMinus_31_43, "hfSumMinus_31_43/F");
  tree->Branch("hfSumMinus_31_44", &hfSumMinus_31_44, "hfSumMinus_31_44/F");
  tree->Branch("hfSumMinus_31_45", &hfSumMinus_31_45, "hfSumMinus_31_45/F");
  tree->Branch("hfSumMinus_31_46", &hfSumMinus_31_46, "hfSumMinus_31_46/F");
  tree->Branch("hfSumMinus_31_47", &hfSumMinus_31_47, "hfSumMinus_31_47/F");
  tree->Branch("hfSumMinus_31_48", &hfSumMinus_31_48, "hfSumMinus_31_48/F");
  tree->Branch("hfSumMinus_31_49", &hfSumMinus_31_49, "hfSumMinus_31_49/F");
  tree->Branch("hfSumMinus_31_50", &hfSumMinus_31_50, "hfSumMinus_31_50/F");
  tree->Branch("hfSumMinus_31_51", &hfSumMinus_31_51, "hfSumMinus_31_51/F");
  tree->Branch("hfSumMinus_31_52", &hfSumMinus_31_52, "hfSumMinus_31_52/F");

  tree->Branch("hfSumMinus_32_33", &hfSumMinus_32_33, "hfSumMinus_32_33/F");
  tree->Branch("hfSumMinus_32_34", &hfSumMinus_32_34, "hfSumMinus_32_34/F");
  tree->Branch("hfSumMinus_32_35", &hfSumMinus_32_35, "hfSumMinus_32_35/F");
  tree->Branch("hfSumMinus_32_36", &hfSumMinus_32_36, "hfSumMinus_32_36/F");
  tree->Branch("hfSumMinus_32_37", &hfSumMinus_32_37, "hfSumMinus_32_37/F");
  tree->Branch("hfSumMinus_32_38", &hfSumMinus_32_38, "hfSumMinus_32_38/F");
  tree->Branch("hfSumMinus_32_39", &hfSumMinus_32_39, "hfSumMinus_32_39/F");
  tree->Branch("hfSumMinus_32_40", &hfSumMinus_32_40, "hfSumMinus_32_40/F");
  tree->Branch("hfSumMinus_32_41", &hfSumMinus_32_41, "hfSumMinus_32_41/F");
  tree->Branch("hfSumMinus_32_42", &hfSumMinus_32_42, "hfSumMinus_32_42/F");
  tree->Branch("hfSumMinus_32_43", &hfSumMinus_32_43, "hfSumMinus_32_43/F");
  tree->Branch("hfSumMinus_32_44", &hfSumMinus_32_44, "hfSumMinus_32_44/F");
  tree->Branch("hfSumMinus_32_45", &hfSumMinus_32_45, "hfSumMinus_32_45/F");
  tree->Branch("hfSumMinus_32_46", &hfSumMinus_32_46, "hfSumMinus_32_46/F");
  tree->Branch("hfSumMinus_32_47", &hfSumMinus_32_47, "hfSumMinus_32_47/F");
  tree->Branch("hfSumMinus_32_48", &hfSumMinus_32_48, "hfSumMinus_32_48/F");
  tree->Branch("hfSumMinus_32_49", &hfSumMinus_32_49, "hfSumMinus_32_49/F");
  tree->Branch("hfSumMinus_32_50", &hfSumMinus_32_50, "hfSumMinus_32_50/F");
  tree->Branch("hfSumMinus_32_51", &hfSumMinus_32_51, "hfSumMinus_32_51/F");
  tree->Branch("hfSumMinus_32_52", &hfSumMinus_32_52, "hfSumMinus_32_52/F");

  tree->Branch("hfSumMinus_33_34", &hfSumMinus_33_34, "hfSumMinus_33_34/F");
  tree->Branch("hfSumMinus_33_35", &hfSumMinus_33_35, "hfSumMinus_33_35/F");
  tree->Branch("hfSumMinus_33_36", &hfSumMinus_33_36, "hfSumMinus_33_36/F");
  tree->Branch("hfSumMinus_33_37", &hfSumMinus_33_37, "hfSumMinus_33_37/F");
  tree->Branch("hfSumMinus_33_38", &hfSumMinus_33_38, "hfSumMinus_33_38/F");
  tree->Branch("hfSumMinus_33_39", &hfSumMinus_33_39, "hfSumMinus_33_39/F");
  tree->Branch("hfSumMinus_33_40", &hfSumMinus_33_40, "hfSumMinus_33_40/F");
  tree->Branch("hfSumMinus_33_41", &hfSumMinus_33_41, "hfSumMinus_33_41/F");
  tree->Branch("hfSumMinus_33_42", &hfSumMinus_33_42, "hfSumMinus_33_42/F");
  tree->Branch("hfSumMinus_33_43", &hfSumMinus_33_43, "hfSumMinus_33_43/F");
  tree->Branch("hfSumMinus_33_44", &hfSumMinus_33_44, "hfSumMinus_33_44/F");
  tree->Branch("hfSumMinus_33_45", &hfSumMinus_33_45, "hfSumMinus_33_45/F");
  tree->Branch("hfSumMinus_33_46", &hfSumMinus_33_46, "hfSumMinus_33_46/F");
  tree->Branch("hfSumMinus_33_47", &hfSumMinus_33_47, "hfSumMinus_33_47/F");
  tree->Branch("hfSumMinus_33_48", &hfSumMinus_33_48, "hfSumMinus_33_48/F");
  tree->Branch("hfSumMinus_33_49", &hfSumMinus_33_49, "hfSumMinus_33_49/F");
  tree->Branch("hfSumMinus_33_50", &hfSumMinus_33_50, "hfSumMinus_33_50/F");
  tree->Branch("hfSumMinus_33_51", &hfSumMinus_33_51, "hfSumMinus_33_51/F");
  tree->Branch("hfSumMinus_33_52", &hfSumMinus_33_52, "hfSumMinus_33_52/F");

  tree->Branch("hfSumMinus_34_35", &hfSumMinus_34_35, "hfSumMinus_34_35/F");
  tree->Branch("hfSumMinus_34_36", &hfSumMinus_34_36, "hfSumMinus_34_36/F");
  tree->Branch("hfSumMinus_34_37", &hfSumMinus_34_37, "hfSumMinus_34_37/F");
  tree->Branch("hfSumMinus_34_38", &hfSumMinus_34_38, "hfSumMinus_34_38/F");
  tree->Branch("hfSumMinus_34_39", &hfSumMinus_34_39, "hfSumMinus_34_39/F");
  tree->Branch("hfSumMinus_34_40", &hfSumMinus_34_40, "hfSumMinus_34_40/F");
  tree->Branch("hfSumMinus_34_41", &hfSumMinus_34_41, "hfSumMinus_34_41/F");
  tree->Branch("hfSumMinus_34_42", &hfSumMinus_34_42, "hfSumMinus_34_42/F");
  tree->Branch("hfSumMinus_34_43", &hfSumMinus_34_43, "hfSumMinus_34_43/F");
  tree->Branch("hfSumMinus_34_44", &hfSumMinus_34_44, "hfSumMinus_34_44/F");
  tree->Branch("hfSumMinus_34_45", &hfSumMinus_34_45, "hfSumMinus_34_45/F");
  tree->Branch("hfSumMinus_34_46", &hfSumMinus_34_46, "hfSumMinus_34_46/F");
  tree->Branch("hfSumMinus_34_47", &hfSumMinus_34_47, "hfSumMinus_34_47/F");
  tree->Branch("hfSumMinus_34_48", &hfSumMinus_34_48, "hfSumMinus_34_48/F");
  tree->Branch("hfSumMinus_34_49", &hfSumMinus_34_49, "hfSumMinus_34_49/F");
  tree->Branch("hfSumMinus_34_50", &hfSumMinus_34_50, "hfSumMinus_34_50/F");
  tree->Branch("hfSumMinus_34_51", &hfSumMinus_34_51, "hfSumMinus_34_51/F");
  tree->Branch("hfSumMinus_34_52", &hfSumMinus_34_52, "hfSumMinus_34_52/F");

  tree->Branch("hfSumMinus_35_36", &hfSumMinus_35_36, "hfSumMinus_35_36/F");
  tree->Branch("hfSumMinus_35_37", &hfSumMinus_35_37, "hfSumMinus_35_37/F");
  tree->Branch("hfSumMinus_35_38", &hfSumMinus_35_38, "hfSumMinus_35_38/F");
  tree->Branch("hfSumMinus_35_39", &hfSumMinus_35_39, "hfSumMinus_35_39/F");
  tree->Branch("hfSumMinus_35_40", &hfSumMinus_35_40, "hfSumMinus_35_40/F");
  tree->Branch("hfSumMinus_35_41", &hfSumMinus_35_41, "hfSumMinus_35_41/F");
  tree->Branch("hfSumMinus_35_42", &hfSumMinus_35_42, "hfSumMinus_35_42/F");
  tree->Branch("hfSumMinus_35_43", &hfSumMinus_35_43, "hfSumMinus_35_43/F");
  tree->Branch("hfSumMinus_35_44", &hfSumMinus_35_44, "hfSumMinus_35_44/F");
  tree->Branch("hfSumMinus_35_45", &hfSumMinus_35_45, "hfSumMinus_35_45/F");
  tree->Branch("hfSumMinus_35_46", &hfSumMinus_35_46, "hfSumMinus_35_46/F");
  tree->Branch("hfSumMinus_35_47", &hfSumMinus_35_47, "hfSumMinus_35_47/F");
  tree->Branch("hfSumMinus_35_48", &hfSumMinus_35_48, "hfSumMinus_35_48/F");
  tree->Branch("hfSumMinus_35_49", &hfSumMinus_35_49, "hfSumMinus_35_49/F");
  tree->Branch("hfSumMinus_35_50", &hfSumMinus_35_50, "hfSumMinus_35_50/F");
  tree->Branch("hfSumMinus_35_51", &hfSumMinus_35_51, "hfSumMinus_35_51/F");
  tree->Branch("hfSumMinus_35_52", &hfSumMinus_35_52, "hfSumMinus_35_52/F");

  tree->Branch("hfSumMinus_36_37", &hfSumMinus_36_37, "hfSumMinus_36_37/F");
  tree->Branch("hfSumMinus_36_38", &hfSumMinus_36_38, "hfSumMinus_36_38/F");
  tree->Branch("hfSumMinus_36_39", &hfSumMinus_36_39, "hfSumMinus_36_39/F");
  tree->Branch("hfSumMinus_36_40", &hfSumMinus_36_40, "hfSumMinus_36_40/F");
  tree->Branch("hfSumMinus_36_41", &hfSumMinus_36_41, "hfSumMinus_36_41/F");
  tree->Branch("hfSumMinus_36_42", &hfSumMinus_36_42, "hfSumMinus_36_42/F");
  tree->Branch("hfSumMinus_36_43", &hfSumMinus_36_43, "hfSumMinus_36_43/F");
  tree->Branch("hfSumMinus_36_44", &hfSumMinus_36_44, "hfSumMinus_36_44/F");
  tree->Branch("hfSumMinus_36_45", &hfSumMinus_36_45, "hfSumMinus_36_45/F");
  tree->Branch("hfSumMinus_36_46", &hfSumMinus_36_46, "hfSumMinus_36_46/F");
  tree->Branch("hfSumMinus_36_47", &hfSumMinus_36_47, "hfSumMinus_36_47/F");
  tree->Branch("hfSumMinus_36_48", &hfSumMinus_36_48, "hfSumMinus_36_48/F");
  tree->Branch("hfSumMinus_36_49", &hfSumMinus_36_49, "hfSumMinus_36_49/F");
  tree->Branch("hfSumMinus_36_50", &hfSumMinus_36_50, "hfSumMinus_36_50/F");
  tree->Branch("hfSumMinus_36_51", &hfSumMinus_36_51, "hfSumMinus_36_51/F");
  tree->Branch("hfSumMinus_36_52", &hfSumMinus_36_52, "hfSumMinus_36_52/F");

  tree->Branch("hfSumMinus_37_38", &hfSumMinus_37_38, "hfSumMinus_37_38/F");
  tree->Branch("hfSumMinus_37_39", &hfSumMinus_37_39, "hfSumMinus_37_39/F");
  tree->Branch("hfSumMinus_37_40", &hfSumMinus_37_40, "hfSumMinus_37_40/F");
  tree->Branch("hfSumMinus_37_41", &hfSumMinus_37_41, "hfSumMinus_37_41/F");
  tree->Branch("hfSumMinus_37_42", &hfSumMinus_37_42, "hfSumMinus_37_42/F");
  tree->Branch("hfSumMinus_37_43", &hfSumMinus_37_43, "hfSumMinus_37_43/F");
  tree->Branch("hfSumMinus_37_44", &hfSumMinus_37_44, "hfSumMinus_37_44/F");
  tree->Branch("hfSumMinus_37_45", &hfSumMinus_37_45, "hfSumMinus_37_45/F");
  tree->Branch("hfSumMinus_37_46", &hfSumMinus_37_46, "hfSumMinus_37_46/F");
  tree->Branch("hfSumMinus_37_47", &hfSumMinus_37_47, "hfSumMinus_37_47/F");
  tree->Branch("hfSumMinus_37_48", &hfSumMinus_37_48, "hfSumMinus_37_48/F");
  tree->Branch("hfSumMinus_37_49", &hfSumMinus_37_49, "hfSumMinus_37_49/F");
  tree->Branch("hfSumMinus_37_50", &hfSumMinus_37_50, "hfSumMinus_37_50/F");
  tree->Branch("hfSumMinus_37_51", &hfSumMinus_37_51, "hfSumMinus_37_51/F");
  tree->Branch("hfSumMinus_37_52", &hfSumMinus_37_52, "hfSumMinus_37_52/F");

  tree->Branch("hfSumMinus_38_39", &hfSumMinus_38_39, "hfSumMinus_38_39/F");
  tree->Branch("hfSumMinus_38_40", &hfSumMinus_38_40, "hfSumMinus_38_40/F");
  tree->Branch("hfSumMinus_38_41", &hfSumMinus_38_41, "hfSumMinus_38_41/F");
  tree->Branch("hfSumMinus_38_42", &hfSumMinus_38_42, "hfSumMinus_38_42/F");
  tree->Branch("hfSumMinus_38_43", &hfSumMinus_38_43, "hfSumMinus_38_43/F");
  tree->Branch("hfSumMinus_38_44", &hfSumMinus_38_44, "hfSumMinus_38_44/F");
  tree->Branch("hfSumMinus_38_45", &hfSumMinus_38_45, "hfSumMinus_38_45/F");
  tree->Branch("hfSumMinus_38_46", &hfSumMinus_38_46, "hfSumMinus_38_46/F");
  tree->Branch("hfSumMinus_38_47", &hfSumMinus_38_47, "hfSumMinus_38_47/F");
  tree->Branch("hfSumMinus_38_48", &hfSumMinus_38_48, "hfSumMinus_38_48/F");
  tree->Branch("hfSumMinus_38_49", &hfSumMinus_38_49, "hfSumMinus_38_49/F");
  tree->Branch("hfSumMinus_38_50", &hfSumMinus_38_50, "hfSumMinus_38_50/F");
  tree->Branch("hfSumMinus_38_51", &hfSumMinus_38_51, "hfSumMinus_38_51/F");
  tree->Branch("hfSumMinus_38_52", &hfSumMinus_38_52, "hfSumMinus_38_52/F");

  tree->Branch("hfSumMinus_39_40", &hfSumMinus_39_40, "hfSumMinus_39_40/F");
  tree->Branch("hfSumMinus_39_41", &hfSumMinus_39_41, "hfSumMinus_39_41/F");
  tree->Branch("hfSumMinus_39_42", &hfSumMinus_39_42, "hfSumMinus_39_42/F");
  tree->Branch("hfSumMinus_39_43", &hfSumMinus_39_43, "hfSumMinus_39_43/F");
  tree->Branch("hfSumMinus_39_44", &hfSumMinus_39_44, "hfSumMinus_39_44/F");
  tree->Branch("hfSumMinus_39_45", &hfSumMinus_39_45, "hfSumMinus_39_45/F");
  tree->Branch("hfSumMinus_39_46", &hfSumMinus_39_46, "hfSumMinus_39_46/F");
  tree->Branch("hfSumMinus_39_47", &hfSumMinus_39_47, "hfSumMinus_39_47/F");
  tree->Branch("hfSumMinus_39_48", &hfSumMinus_39_48, "hfSumMinus_39_48/F");
  tree->Branch("hfSumMinus_39_49", &hfSumMinus_39_49, "hfSumMinus_39_49/F");
  tree->Branch("hfSumMinus_39_50", &hfSumMinus_39_50, "hfSumMinus_39_50/F");
  tree->Branch("hfSumMinus_39_51", &hfSumMinus_39_51, "hfSumMinus_39_51/F");
  tree->Branch("hfSumMinus_39_52", &hfSumMinus_39_52, "hfSumMinus_39_52/F");

  tree->Branch("hfSumMinus_40_41", &hfSumMinus_40_41, "hfSumMinus_40_41/F");
  tree->Branch("hfSumMinus_40_42", &hfSumMinus_40_42, "hfSumMinus_40_42/F");
  tree->Branch("hfSumMinus_40_43", &hfSumMinus_40_43, "hfSumMinus_40_43/F");
  tree->Branch("hfSumMinus_40_44", &hfSumMinus_40_44, "hfSumMinus_40_44/F");
  tree->Branch("hfSumMinus_40_45", &hfSumMinus_40_45, "hfSumMinus_40_45/F");
  tree->Branch("hfSumMinus_40_46", &hfSumMinus_40_46, "hfSumMinus_40_46/F");
  tree->Branch("hfSumMinus_40_47", &hfSumMinus_40_47, "hfSumMinus_40_47/F");
  tree->Branch("hfSumMinus_40_48", &hfSumMinus_40_48, "hfSumMinus_40_48/F");
  tree->Branch("hfSumMinus_40_49", &hfSumMinus_40_49, "hfSumMinus_40_49/F");
  tree->Branch("hfSumMinus_40_50", &hfSumMinus_40_50, "hfSumMinus_40_50/F");
  tree->Branch("hfSumMinus_40_51", &hfSumMinus_40_51, "hfSumMinus_40_51/F");
  tree->Branch("hfSumMinus_40_52", &hfSumMinus_40_52, "hfSumMinus_40_52/F");

  tree->Branch("hfSumMinus_41_42", &hfSumMinus_41_42, "hfSumMinus_41_42/F");
  tree->Branch("hfSumMinus_41_43", &hfSumMinus_41_43, "hfSumMinus_41_43/F");
  tree->Branch("hfSumMinus_41_44", &hfSumMinus_41_44, "hfSumMinus_41_44/F");
  tree->Branch("hfSumMinus_41_45", &hfSumMinus_41_45, "hfSumMinus_41_45/F");
  tree->Branch("hfSumMinus_41_46", &hfSumMinus_41_46, "hfSumMinus_41_46/F");
  tree->Branch("hfSumMinus_41_47", &hfSumMinus_41_47, "hfSumMinus_41_47/F");
  tree->Branch("hfSumMinus_41_48", &hfSumMinus_41_48, "hfSumMinus_41_48/F");
  tree->Branch("hfSumMinus_41_49", &hfSumMinus_41_49, "hfSumMinus_41_49/F");
  tree->Branch("hfSumMinus_41_50", &hfSumMinus_41_50, "hfSumMinus_41_50/F");
  tree->Branch("hfSumMinus_41_51", &hfSumMinus_41_51, "hfSumMinus_41_51/F");
  tree->Branch("hfSumMinus_41_52", &hfSumMinus_41_52, "hfSumMinus_41_52/F");

  tree->Branch("hfSumMinus_42_43", &hfSumMinus_42_43, "hfSumMinus_42_43/F");
  tree->Branch("hfSumMinus_42_44", &hfSumMinus_42_44, "hfSumMinus_42_44/F");
  tree->Branch("hfSumMinus_42_45", &hfSumMinus_42_45, "hfSumMinus_42_45/F");
  tree->Branch("hfSumMinus_42_46", &hfSumMinus_42_46, "hfSumMinus_42_46/F");
  tree->Branch("hfSumMinus_42_47", &hfSumMinus_42_47, "hfSumMinus_42_47/F");
  tree->Branch("hfSumMinus_42_48", &hfSumMinus_42_48, "hfSumMinus_42_48/F");
  tree->Branch("hfSumMinus_42_49", &hfSumMinus_42_49, "hfSumMinus_42_49/F");
  tree->Branch("hfSumMinus_42_50", &hfSumMinus_42_50, "hfSumMinus_42_50/F");
  tree->Branch("hfSumMinus_42_51", &hfSumMinus_42_51, "hfSumMinus_42_51/F");
  tree->Branch("hfSumMinus_42_52", &hfSumMinus_42_52, "hfSumMinus_42_52/F");

  tree->Branch("hfSumMinus_43_44", &hfSumMinus_43_44, "hfSumMinus_43_44/F");
  tree->Branch("hfSumMinus_43_45", &hfSumMinus_43_45, "hfSumMinus_43_45/F");
  tree->Branch("hfSumMinus_43_46", &hfSumMinus_43_46, "hfSumMinus_43_46/F");
  tree->Branch("hfSumMinus_43_47", &hfSumMinus_43_47, "hfSumMinus_43_47/F");
  tree->Branch("hfSumMinus_43_48", &hfSumMinus_43_48, "hfSumMinus_43_48/F");
  tree->Branch("hfSumMinus_43_49", &hfSumMinus_43_49, "hfSumMinus_43_49/F");
  tree->Branch("hfSumMinus_43_50", &hfSumMinus_43_50, "hfSumMinus_43_50/F");
  tree->Branch("hfSumMinus_43_51", &hfSumMinus_43_51, "hfSumMinus_43_51/F");
  tree->Branch("hfSumMinus_43_52", &hfSumMinus_43_52, "hfSumMinus_43_52/F");

  tree->Branch("hfSumMinus_44_45", &hfSumMinus_44_45, "hfSumMinus_44_45/F");
  tree->Branch("hfSumMinus_44_46", &hfSumMinus_44_46, "hfSumMinus_44_46/F");
  tree->Branch("hfSumMinus_44_47", &hfSumMinus_44_47, "hfSumMinus_44_47/F");
  tree->Branch("hfSumMinus_44_48", &hfSumMinus_44_48, "hfSumMinus_44_48/F");
  tree->Branch("hfSumMinus_44_49", &hfSumMinus_44_49, "hfSumMinus_44_49/F");
  tree->Branch("hfSumMinus_44_50", &hfSumMinus_44_50, "hfSumMinus_44_50/F");
  tree->Branch("hfSumMinus_44_51", &hfSumMinus_44_51, "hfSumMinus_44_51/F");
  tree->Branch("hfSumMinus_44_52", &hfSumMinus_44_52, "hfSumMinus_44_52/F");

  tree->Branch("hfSumMinus_45_46", &hfSumMinus_45_46, "hfSumMinus_45_46/F");
  tree->Branch("hfSumMinus_45_47", &hfSumMinus_45_47, "hfSumMinus_45_47/F");
  tree->Branch("hfSumMinus_45_48", &hfSumMinus_45_48, "hfSumMinus_45_48/F");
  tree->Branch("hfSumMinus_45_49", &hfSumMinus_45_49, "hfSumMinus_45_49/F");
  tree->Branch("hfSumMinus_45_50", &hfSumMinus_45_50, "hfSumMinus_45_50/F");
  tree->Branch("hfSumMinus_45_51", &hfSumMinus_45_51, "hfSumMinus_45_51/F");
  tree->Branch("hfSumMinus_45_52", &hfSumMinus_45_52, "hfSumMinus_45_52/F");

  tree->Branch("hfSumMinus_46_47", &hfSumMinus_46_47, "hfSumMinus_46_47/F");
  tree->Branch("hfSumMinus_46_48", &hfSumMinus_46_48, "hfSumMinus_46_48/F");
  tree->Branch("hfSumMinus_46_49", &hfSumMinus_46_49, "hfSumMinus_46_49/F");
  tree->Branch("hfSumMinus_46_50", &hfSumMinus_46_50, "hfSumMinus_46_50/F");
  tree->Branch("hfSumMinus_46_51", &hfSumMinus_46_51, "hfSumMinus_46_51/F");
  tree->Branch("hfSumMinus_46_52", &hfSumMinus_46_52, "hfSumMinus_46_52/F");

  tree->Branch("hfSumMinus_47_48", &hfSumMinus_47_48, "hfSumMinus_47_48/F");
  tree->Branch("hfSumMinus_47_49", &hfSumMinus_47_49, "hfSumMinus_47_49/F");
  tree->Branch("hfSumMinus_47_50", &hfSumMinus_47_50, "hfSumMinus_47_50/F");
  tree->Branch("hfSumMinus_47_51", &hfSumMinus_47_51, "hfSumMinus_47_51/F");
  tree->Branch("hfSumMinus_47_52", &hfSumMinus_47_52, "hfSumMinus_47_52/F");

  tree->Branch("hfSumMinus_48_49", &hfSumMinus_48_49, "hfSumMinus_48_49/F");
  tree->Branch("hfSumMinus_48_50", &hfSumMinus_48_50, "hfSumMinus_48_50/F");
  tree->Branch("hfSumMinus_48_51", &hfSumMinus_48_51, "hfSumMinus_48_51/F");
  tree->Branch("hfSumMinus_48_52", &hfSumMinus_48_52, "hfSumMinus_48_52/F");

  tree->Branch("hfSumMinus_49_50", &hfSumMinus_49_50, "hfSumMinus_49_50/F");
  tree->Branch("hfSumMinus_49_51", &hfSumMinus_49_51, "hfSumMinus_49_51/F");
  tree->Branch("hfSumMinus_49_52", &hfSumMinus_49_52, "hfSumMinus_49_52/F");

  tree->Branch("hfSumMinus_50_51", &hfSumMinus_50_51, "hfSumMinus_50_51/F");
  tree->Branch("hfSumMinus_50_52", &hfSumMinus_50_52, "hfSumMinus_50_52/F");

  tree->Branch("hfSumMinus_51_52", &hfSumMinus_51_52, "hfSumMinus_51_52/F");
  // hfsumminus --end



  // New trees
  tree->Branch("minOfHFSums_30_31", &minOfHFSums_30_31, "minOfHFSums_30_31/F");
  tree->Branch("minOfHFSums_30_32", &minOfHFSums_30_32, "minOfHFSums_30_32/F");
  tree->Branch("minOfHFSums_30_33", &minOfHFSums_30_33, "minOfHFSums_30_33/F");
  tree->Branch("minOfHFSums_30_34", &minOfHFSums_30_34, "minOfHFSums_30_34/F");
  tree->Branch("minOfHFSums_30_35", &minOfHFSums_30_35, "minOfHFSums_30_35/F");
  tree->Branch("minOfHFSums_30_36", &minOfHFSums_30_36, "minOfHFSums_30_36/F");
  tree->Branch("minOfHFSums_30_37", &minOfHFSums_30_37, "minOfHFSums_30_37/F");
  tree->Branch("minOfHFSums_30_38", &minOfHFSums_30_38, "minOfHFSums_30_38/F");
  tree->Branch("minOfHFSums_30_39", &minOfHFSums_30_39, "minOfHFSums_30_39/F");
  tree->Branch("minOfHFSums_30_40", &minOfHFSums_30_40, "minOfHFSums_30_40/F");
  tree->Branch("minOfHFSums_30_41", &minOfHFSums_30_41, "minOfHFSums_30_41/F");
  tree->Branch("minOfHFSums_30_42", &minOfHFSums_30_42, "minOfHFSums_30_42/F");
  tree->Branch("minOfHFSums_30_43", &minOfHFSums_30_43, "minOfHFSums_30_43/F");
  tree->Branch("minOfHFSums_30_44", &minOfHFSums_30_44, "minOfHFSums_30_44/F");
  tree->Branch("minOfHFSums_30_45", &minOfHFSums_30_45, "minOfHFSums_30_45/F");
  tree->Branch("minOfHFSums_30_46", &minOfHFSums_30_46, "minOfHFSums_30_46/F");
  tree->Branch("minOfHFSums_30_47", &minOfHFSums_30_47, "minOfHFSums_30_47/F");
  tree->Branch("minOfHFSums_30_48", &minOfHFSums_30_48, "minOfHFSums_30_48/F");
  tree->Branch("minOfHFSums_30_49", &minOfHFSums_30_49, "minOfHFSums_30_49/F");
  tree->Branch("minOfHFSums_30_50", &minOfHFSums_30_50, "minOfHFSums_30_50/F");
  tree->Branch("minOfHFSums_30_51", &minOfHFSums_30_51, "minOfHFSums_30_51/F");
  tree->Branch("minOfHFSums_30_52", &minOfHFSums_30_52, "minOfHFSums_30_52/F");

  tree->Branch("minOfHFSums_31_32", &minOfHFSums_31_32, "minOfHFSums_31_32/F");
  tree->Branch("minOfHFSums_31_33", &minOfHFSums_31_33, "minOfHFSums_31_33/F");
  tree->Branch("minOfHFSums_31_34", &minOfHFSums_31_34, "minOfHFSums_31_34/F");
  tree->Branch("minOfHFSums_31_35", &minOfHFSums_31_35, "minOfHFSums_31_35/F");
  tree->Branch("minOfHFSums_31_36", &minOfHFSums_31_36, "minOfHFSums_31_36/F");
  tree->Branch("minOfHFSums_31_37", &minOfHFSums_31_37, "minOfHFSums_31_37/F");
  tree->Branch("minOfHFSums_31_38", &minOfHFSums_31_38, "minOfHFSums_31_38/F");
  tree->Branch("minOfHFSums_31_39", &minOfHFSums_31_39, "minOfHFSums_31_39/F");
  tree->Branch("minOfHFSums_31_40", &minOfHFSums_31_40, "minOfHFSums_31_40/F");
  tree->Branch("minOfHFSums_31_41", &minOfHFSums_31_41, "minOfHFSums_31_41/F");
  tree->Branch("minOfHFSums_31_42", &minOfHFSums_31_42, "minOfHFSums_31_42/F");
  tree->Branch("minOfHFSums_31_43", &minOfHFSums_31_43, "minOfHFSums_31_43/F");
  tree->Branch("minOfHFSums_31_44", &minOfHFSums_31_44, "minOfHFSums_31_44/F");
  tree->Branch("minOfHFSums_31_45", &minOfHFSums_31_45, "minOfHFSums_31_45/F");
  tree->Branch("minOfHFSums_31_46", &minOfHFSums_31_46, "minOfHFSums_31_46/F");
  tree->Branch("minOfHFSums_31_47", &minOfHFSums_31_47, "minOfHFSums_31_47/F");
  tree->Branch("minOfHFSums_31_48", &minOfHFSums_31_48, "minOfHFSums_31_48/F");
  tree->Branch("minOfHFSums_31_49", &minOfHFSums_31_49, "minOfHFSums_31_49/F");
  tree->Branch("minOfHFSums_31_50", &minOfHFSums_31_50, "minOfHFSums_31_50/F");
  tree->Branch("minOfHFSums_31_51", &minOfHFSums_31_51, "minOfHFSums_31_51/F");
  tree->Branch("minOfHFSums_31_52", &minOfHFSums_31_52, "minOfHFSums_31_52/F");

  tree->Branch("minOfHFSums_32_33", &minOfHFSums_32_33, "minOfHFSums_32_33/F");
  tree->Branch("minOfHFSums_32_34", &minOfHFSums_32_34, "minOfHFSums_32_34/F");
  tree->Branch("minOfHFSums_32_35", &minOfHFSums_32_35, "minOfHFSums_32_35/F");
  tree->Branch("minOfHFSums_32_36", &minOfHFSums_32_36, "minOfHFSums_32_36/F");
  tree->Branch("minOfHFSums_32_37", &minOfHFSums_32_37, "minOfHFSums_32_37/F");
  tree->Branch("minOfHFSums_32_38", &minOfHFSums_32_38, "minOfHFSums_32_38/F");
  tree->Branch("minOfHFSums_32_39", &minOfHFSums_32_39, "minOfHFSums_32_39/F");
  tree->Branch("minOfHFSums_32_40", &minOfHFSums_32_40, "minOfHFSums_32_40/F");
  tree->Branch("minOfHFSums_32_41", &minOfHFSums_32_41, "minOfHFSums_32_41/F");
  tree->Branch("minOfHFSums_32_42", &minOfHFSums_32_42, "minOfHFSums_32_42/F");
  tree->Branch("minOfHFSums_32_43", &minOfHFSums_32_43, "minOfHFSums_32_43/F");
  tree->Branch("minOfHFSums_32_44", &minOfHFSums_32_44, "minOfHFSums_32_44/F");
  tree->Branch("minOfHFSums_32_45", &minOfHFSums_32_45, "minOfHFSums_32_45/F");
  tree->Branch("minOfHFSums_32_46", &minOfHFSums_32_46, "minOfHFSums_32_46/F");
  tree->Branch("minOfHFSums_32_47", &minOfHFSums_32_47, "minOfHFSums_32_47/F");
  tree->Branch("minOfHFSums_32_48", &minOfHFSums_32_48, "minOfHFSums_32_48/F");
  tree->Branch("minOfHFSums_32_49", &minOfHFSums_32_49, "minOfHFSums_32_49/F");
  tree->Branch("minOfHFSums_32_50", &minOfHFSums_32_50, "minOfHFSums_32_50/F");
  tree->Branch("minOfHFSums_32_51", &minOfHFSums_32_51, "minOfHFSums_32_51/F");
  tree->Branch("minOfHFSums_32_52", &minOfHFSums_32_52, "minOfHFSums_32_52/F");

  tree->Branch("minOfHFSums_33_34", &minOfHFSums_33_34, "minOfHFSums_33_34/F");
  tree->Branch("minOfHFSums_33_35", &minOfHFSums_33_35, "minOfHFSums_33_35/F");
  tree->Branch("minOfHFSums_33_36", &minOfHFSums_33_36, "minOfHFSums_33_36/F");
  tree->Branch("minOfHFSums_33_37", &minOfHFSums_33_37, "minOfHFSums_33_37/F");
  tree->Branch("minOfHFSums_33_38", &minOfHFSums_33_38, "minOfHFSums_33_38/F");
  tree->Branch("minOfHFSums_33_39", &minOfHFSums_33_39, "minOfHFSums_33_39/F");
  tree->Branch("minOfHFSums_33_40", &minOfHFSums_33_40, "minOfHFSums_33_40/F");
  tree->Branch("minOfHFSums_33_41", &minOfHFSums_33_41, "minOfHFSums_33_41/F");
  tree->Branch("minOfHFSums_33_42", &minOfHFSums_33_42, "minOfHFSums_33_42/F");
  tree->Branch("minOfHFSums_33_43", &minOfHFSums_33_43, "minOfHFSums_33_43/F");
  tree->Branch("minOfHFSums_33_44", &minOfHFSums_33_44, "minOfHFSums_33_44/F");
  tree->Branch("minOfHFSums_33_45", &minOfHFSums_33_45, "minOfHFSums_33_45/F");
  tree->Branch("minOfHFSums_33_46", &minOfHFSums_33_46, "minOfHFSums_33_46/F");
  tree->Branch("minOfHFSums_33_47", &minOfHFSums_33_47, "minOfHFSums_33_47/F");
  tree->Branch("minOfHFSums_33_48", &minOfHFSums_33_48, "minOfHFSums_33_48/F");
  tree->Branch("minOfHFSums_33_49", &minOfHFSums_33_49, "minOfHFSums_33_49/F");
  tree->Branch("minOfHFSums_33_50", &minOfHFSums_33_50, "minOfHFSums_33_50/F");
  tree->Branch("minOfHFSums_33_51", &minOfHFSums_33_51, "minOfHFSums_33_51/F");
  tree->Branch("minOfHFSums_33_52", &minOfHFSums_33_52, "minOfHFSums_33_52/F");

  tree->Branch("minOfHFSums_34_35", &minOfHFSums_34_35, "minOfHFSums_34_35/F");
  tree->Branch("minOfHFSums_34_36", &minOfHFSums_34_36, "minOfHFSums_34_36/F");
  tree->Branch("minOfHFSums_34_37", &minOfHFSums_34_37, "minOfHFSums_34_37/F");
  tree->Branch("minOfHFSums_34_38", &minOfHFSums_34_38, "minOfHFSums_34_38/F");
  tree->Branch("minOfHFSums_34_39", &minOfHFSums_34_39, "minOfHFSums_34_39/F");
  tree->Branch("minOfHFSums_34_40", &minOfHFSums_34_40, "minOfHFSums_34_40/F");
  tree->Branch("minOfHFSums_34_41", &minOfHFSums_34_41, "minOfHFSums_34_41/F");
  tree->Branch("minOfHFSums_34_42", &minOfHFSums_34_42, "minOfHFSums_34_42/F");
  tree->Branch("minOfHFSums_34_43", &minOfHFSums_34_43, "minOfHFSums_34_43/F");
  tree->Branch("minOfHFSums_34_44", &minOfHFSums_34_44, "minOfHFSums_34_44/F");
  tree->Branch("minOfHFSums_34_45", &minOfHFSums_34_45, "minOfHFSums_34_45/F");
  tree->Branch("minOfHFSums_34_46", &minOfHFSums_34_46, "minOfHFSums_34_46/F");
  tree->Branch("minOfHFSums_34_47", &minOfHFSums_34_47, "minOfHFSums_34_47/F");
  tree->Branch("minOfHFSums_34_48", &minOfHFSums_34_48, "minOfHFSums_34_48/F");
  tree->Branch("minOfHFSums_34_49", &minOfHFSums_34_49, "minOfHFSums_34_49/F");
  tree->Branch("minOfHFSums_34_50", &minOfHFSums_34_50, "minOfHFSums_34_50/F");
  tree->Branch("minOfHFSums_34_51", &minOfHFSums_34_51, "minOfHFSums_34_51/F");
  tree->Branch("minOfHFSums_34_52", &minOfHFSums_34_52, "minOfHFSums_34_52/F");

  tree->Branch("minOfHFSums_35_36", &minOfHFSums_35_36, "minOfHFSums_35_36/F");
  tree->Branch("minOfHFSums_35_37", &minOfHFSums_35_37, "minOfHFSums_35_37/F");
  tree->Branch("minOfHFSums_35_38", &minOfHFSums_35_38, "minOfHFSums_35_38/F");
  tree->Branch("minOfHFSums_35_39", &minOfHFSums_35_39, "minOfHFSums_35_39/F");
  tree->Branch("minOfHFSums_35_40", &minOfHFSums_35_40, "minOfHFSums_35_40/F");
  tree->Branch("minOfHFSums_35_41", &minOfHFSums_35_41, "minOfHFSums_35_41/F");
  tree->Branch("minOfHFSums_35_42", &minOfHFSums_35_42, "minOfHFSums_35_42/F");
  tree->Branch("minOfHFSums_35_43", &minOfHFSums_35_43, "minOfHFSums_35_43/F");
  tree->Branch("minOfHFSums_35_44", &minOfHFSums_35_44, "minOfHFSums_35_44/F");
  tree->Branch("minOfHFSums_35_45", &minOfHFSums_35_45, "minOfHFSums_35_45/F");
  tree->Branch("minOfHFSums_35_46", &minOfHFSums_35_46, "minOfHFSums_35_46/F");
  tree->Branch("minOfHFSums_35_47", &minOfHFSums_35_47, "minOfHFSums_35_47/F");
  tree->Branch("minOfHFSums_35_48", &minOfHFSums_35_48, "minOfHFSums_35_48/F");
  tree->Branch("minOfHFSums_35_49", &minOfHFSums_35_49, "minOfHFSums_35_49/F");
  tree->Branch("minOfHFSums_35_50", &minOfHFSums_35_50, "minOfHFSums_35_50/F");
  tree->Branch("minOfHFSums_35_51", &minOfHFSums_35_51, "minOfHFSums_35_51/F");
  tree->Branch("minOfHFSums_35_52", &minOfHFSums_35_52, "minOfHFSums_35_52/F");

  tree->Branch("minOfHFSums_36_37", &minOfHFSums_36_37, "minOfHFSums_36_37/F");
  tree->Branch("minOfHFSums_36_38", &minOfHFSums_36_38, "minOfHFSums_36_38/F");
  tree->Branch("minOfHFSums_36_39", &minOfHFSums_36_39, "minOfHFSums_36_39/F");
  tree->Branch("minOfHFSums_36_40", &minOfHFSums_36_40, "minOfHFSums_36_40/F");
  tree->Branch("minOfHFSums_36_41", &minOfHFSums_36_41, "minOfHFSums_36_41/F");
  tree->Branch("minOfHFSums_36_42", &minOfHFSums_36_42, "minOfHFSums_36_42/F");
  tree->Branch("minOfHFSums_36_43", &minOfHFSums_36_43, "minOfHFSums_36_43/F");
  tree->Branch("minOfHFSums_36_44", &minOfHFSums_36_44, "minOfHFSums_36_44/F");
  tree->Branch("minOfHFSums_36_45", &minOfHFSums_36_45, "minOfHFSums_36_45/F");
  tree->Branch("minOfHFSums_36_46", &minOfHFSums_36_46, "minOfHFSums_36_46/F");
  tree->Branch("minOfHFSums_36_47", &minOfHFSums_36_47, "minOfHFSums_36_47/F");
  tree->Branch("minOfHFSums_36_48", &minOfHFSums_36_48, "minOfHFSums_36_48/F");
  tree->Branch("minOfHFSums_36_49", &minOfHFSums_36_49, "minOfHFSums_36_49/F");
  tree->Branch("minOfHFSums_36_50", &minOfHFSums_36_50, "minOfHFSums_36_50/F");
  tree->Branch("minOfHFSums_36_51", &minOfHFSums_36_51, "minOfHFSums_36_51/F");
  tree->Branch("minOfHFSums_36_52", &minOfHFSums_36_52, "minOfHFSums_36_52/F");

  tree->Branch("minOfHFSums_37_38", &minOfHFSums_37_38, "minOfHFSums_37_38/F");
  tree->Branch("minOfHFSums_37_39", &minOfHFSums_37_39, "minOfHFSums_37_39/F");
  tree->Branch("minOfHFSums_37_40", &minOfHFSums_37_40, "minOfHFSums_37_40/F");
  tree->Branch("minOfHFSums_37_41", &minOfHFSums_37_41, "minOfHFSums_37_41/F");
  tree->Branch("minOfHFSums_37_42", &minOfHFSums_37_42, "minOfHFSums_37_42/F");
  tree->Branch("minOfHFSums_37_43", &minOfHFSums_37_43, "minOfHFSums_37_43/F");
  tree->Branch("minOfHFSums_37_44", &minOfHFSums_37_44, "minOfHFSums_37_44/F");
  tree->Branch("minOfHFSums_37_45", &minOfHFSums_37_45, "minOfHFSums_37_45/F");
  tree->Branch("minOfHFSums_37_46", &minOfHFSums_37_46, "minOfHFSums_37_46/F");
  tree->Branch("minOfHFSums_37_47", &minOfHFSums_37_47, "minOfHFSums_37_47/F");
  tree->Branch("minOfHFSums_37_48", &minOfHFSums_37_48, "minOfHFSums_37_48/F");
  tree->Branch("minOfHFSums_37_49", &minOfHFSums_37_49, "minOfHFSums_37_49/F");
  tree->Branch("minOfHFSums_37_50", &minOfHFSums_37_50, "minOfHFSums_37_50/F");
  tree->Branch("minOfHFSums_37_51", &minOfHFSums_37_51, "minOfHFSums_37_51/F");
  tree->Branch("minOfHFSums_37_52", &minOfHFSums_37_52, "minOfHFSums_37_52/F");

  tree->Branch("minOfHFSums_38_39", &minOfHFSums_38_39, "minOfHFSums_38_39/F");
  tree->Branch("minOfHFSums_38_40", &minOfHFSums_38_40, "minOfHFSums_38_40/F");
  tree->Branch("minOfHFSums_38_41", &minOfHFSums_38_41, "minOfHFSums_38_41/F");
  tree->Branch("minOfHFSums_38_42", &minOfHFSums_38_42, "minOfHFSums_38_42/F");
  tree->Branch("minOfHFSums_38_43", &minOfHFSums_38_43, "minOfHFSums_38_43/F");
  tree->Branch("minOfHFSums_38_44", &minOfHFSums_38_44, "minOfHFSums_38_44/F");
  tree->Branch("minOfHFSums_38_45", &minOfHFSums_38_45, "minOfHFSums_38_45/F");
  tree->Branch("minOfHFSums_38_46", &minOfHFSums_38_46, "minOfHFSums_38_46/F");
  tree->Branch("minOfHFSums_38_47", &minOfHFSums_38_47, "minOfHFSums_38_47/F");
  tree->Branch("minOfHFSums_38_48", &minOfHFSums_38_48, "minOfHFSums_38_48/F");
  tree->Branch("minOfHFSums_38_49", &minOfHFSums_38_49, "minOfHFSums_38_49/F");
  tree->Branch("minOfHFSums_38_50", &minOfHFSums_38_50, "minOfHFSums_38_50/F");
  tree->Branch("minOfHFSums_38_51", &minOfHFSums_38_51, "minOfHFSums_38_51/F");
  tree->Branch("minOfHFSums_38_52", &minOfHFSums_38_52, "minOfHFSums_38_52/F");

  tree->Branch("minOfHFSums_39_40", &minOfHFSums_39_40, "minOfHFSums_39_40/F");
  tree->Branch("minOfHFSums_39_41", &minOfHFSums_39_41, "minOfHFSums_39_41/F");
  tree->Branch("minOfHFSums_39_42", &minOfHFSums_39_42, "minOfHFSums_39_42/F");
  tree->Branch("minOfHFSums_39_43", &minOfHFSums_39_43, "minOfHFSums_39_43/F");
  tree->Branch("minOfHFSums_39_44", &minOfHFSums_39_44, "minOfHFSums_39_44/F");
  tree->Branch("minOfHFSums_39_45", &minOfHFSums_39_45, "minOfHFSums_39_45/F");
  tree->Branch("minOfHFSums_39_46", &minOfHFSums_39_46, "minOfHFSums_39_46/F");
  tree->Branch("minOfHFSums_39_47", &minOfHFSums_39_47, "minOfHFSums_39_47/F");
  tree->Branch("minOfHFSums_39_48", &minOfHFSums_39_48, "minOfHFSums_39_48/F");
  tree->Branch("minOfHFSums_39_49", &minOfHFSums_39_49, "minOfHFSums_39_49/F");
  tree->Branch("minOfHFSums_39_50", &minOfHFSums_39_50, "minOfHFSums_39_50/F");
  tree->Branch("minOfHFSums_39_51", &minOfHFSums_39_51, "minOfHFSums_39_51/F");
  tree->Branch("minOfHFSums_39_52", &minOfHFSums_39_52, "minOfHFSums_39_52/F");

  tree->Branch("minOfHFSums_40_41", &minOfHFSums_40_41, "minOfHFSums_40_41/F");
  tree->Branch("minOfHFSums_40_42", &minOfHFSums_40_42, "minOfHFSums_40_42/F");
  tree->Branch("minOfHFSums_40_43", &minOfHFSums_40_43, "minOfHFSums_40_43/F");
  tree->Branch("minOfHFSums_40_44", &minOfHFSums_40_44, "minOfHFSums_40_44/F");
  tree->Branch("minOfHFSums_40_45", &minOfHFSums_40_45, "minOfHFSums_40_45/F");
  tree->Branch("minOfHFSums_40_46", &minOfHFSums_40_46, "minOfHFSums_40_46/F");
  tree->Branch("minOfHFSums_40_47", &minOfHFSums_40_47, "minOfHFSums_40_47/F");
  tree->Branch("minOfHFSums_40_48", &minOfHFSums_40_48, "minOfHFSums_40_48/F");
  tree->Branch("minOfHFSums_40_49", &minOfHFSums_40_49, "minOfHFSums_40_49/F");
  tree->Branch("minOfHFSums_40_50", &minOfHFSums_40_50, "minOfHFSums_40_50/F");
  tree->Branch("minOfHFSums_40_51", &minOfHFSums_40_51, "minOfHFSums_40_51/F");
  tree->Branch("minOfHFSums_40_52", &minOfHFSums_40_52, "minOfHFSums_40_52/F");

  tree->Branch("minOfHFSums_41_42", &minOfHFSums_41_42, "minOfHFSums_41_42/F");
  tree->Branch("minOfHFSums_41_43", &minOfHFSums_41_43, "minOfHFSums_41_43/F");
  tree->Branch("minOfHFSums_41_44", &minOfHFSums_41_44, "minOfHFSums_41_44/F");
  tree->Branch("minOfHFSums_41_45", &minOfHFSums_41_45, "minOfHFSums_41_45/F");
  tree->Branch("minOfHFSums_41_46", &minOfHFSums_41_46, "minOfHFSums_41_46/F");
  tree->Branch("minOfHFSums_41_47", &minOfHFSums_41_47, "minOfHFSums_41_47/F");
  tree->Branch("minOfHFSums_41_48", &minOfHFSums_41_48, "minOfHFSums_41_48/F");
  tree->Branch("minOfHFSums_41_49", &minOfHFSums_41_49, "minOfHFSums_41_49/F");
  tree->Branch("minOfHFSums_41_50", &minOfHFSums_41_50, "minOfHFSums_41_50/F");
  tree->Branch("minOfHFSums_41_51", &minOfHFSums_41_51, "minOfHFSums_41_51/F");
  tree->Branch("minOfHFSums_41_52", &minOfHFSums_41_52, "minOfHFSums_41_52/F");

  tree->Branch("minOfHFSums_42_43", &minOfHFSums_42_43, "minOfHFSums_42_43/F");
  tree->Branch("minOfHFSums_42_44", &minOfHFSums_42_44, "minOfHFSums_42_44/F");
  tree->Branch("minOfHFSums_42_45", &minOfHFSums_42_45, "minOfHFSums_42_45/F");
  tree->Branch("minOfHFSums_42_46", &minOfHFSums_42_46, "minOfHFSums_42_46/F");
  tree->Branch("minOfHFSums_42_47", &minOfHFSums_42_47, "minOfHFSums_42_47/F");
  tree->Branch("minOfHFSums_42_48", &minOfHFSums_42_48, "minOfHFSums_42_48/F");
  tree->Branch("minOfHFSums_42_49", &minOfHFSums_42_49, "minOfHFSums_42_49/F");
  tree->Branch("minOfHFSums_42_50", &minOfHFSums_42_50, "minOfHFSums_42_50/F");
  tree->Branch("minOfHFSums_42_51", &minOfHFSums_42_51, "minOfHFSums_42_51/F");
  tree->Branch("minOfHFSums_42_52", &minOfHFSums_42_52, "minOfHFSums_42_52/F");

  tree->Branch("minOfHFSums_43_44", &minOfHFSums_43_44, "minOfHFSums_43_44/F");
  tree->Branch("minOfHFSums_43_45", &minOfHFSums_43_45, "minOfHFSums_43_45/F");
  tree->Branch("minOfHFSums_43_46", &minOfHFSums_43_46, "minOfHFSums_43_46/F");
  tree->Branch("minOfHFSums_43_47", &minOfHFSums_43_47, "minOfHFSums_43_47/F");
  tree->Branch("minOfHFSums_43_48", &minOfHFSums_43_48, "minOfHFSums_43_48/F");
  tree->Branch("minOfHFSums_43_49", &minOfHFSums_43_49, "minOfHFSums_43_49/F");
  tree->Branch("minOfHFSums_43_50", &minOfHFSums_43_50, "minOfHFSums_43_50/F");
  tree->Branch("minOfHFSums_43_51", &minOfHFSums_43_51, "minOfHFSums_43_51/F");
  tree->Branch("minOfHFSums_43_52", &minOfHFSums_43_52, "minOfHFSums_43_52/F");

  tree->Branch("minOfHFSums_44_45", &minOfHFSums_44_45, "minOfHFSums_44_45/F");
  tree->Branch("minOfHFSums_44_46", &minOfHFSums_44_46, "minOfHFSums_44_46/F");
  tree->Branch("minOfHFSums_44_47", &minOfHFSums_44_47, "minOfHFSums_44_47/F");
  tree->Branch("minOfHFSums_44_48", &minOfHFSums_44_48, "minOfHFSums_44_48/F");
  tree->Branch("minOfHFSums_44_49", &minOfHFSums_44_49, "minOfHFSums_44_49/F");
  tree->Branch("minOfHFSums_44_50", &minOfHFSums_44_50, "minOfHFSums_44_50/F");
  tree->Branch("minOfHFSums_44_51", &minOfHFSums_44_51, "minOfHFSums_44_51/F");
  tree->Branch("minOfHFSums_44_52", &minOfHFSums_44_52, "minOfHFSums_44_52/F");

  tree->Branch("minOfHFSums_45_46", &minOfHFSums_45_46, "minOfHFSums_45_46/F");
  tree->Branch("minOfHFSums_45_47", &minOfHFSums_45_47, "minOfHFSums_45_47/F");
  tree->Branch("minOfHFSums_45_48", &minOfHFSums_45_48, "minOfHFSums_45_48/F");
  tree->Branch("minOfHFSums_45_49", &minOfHFSums_45_49, "minOfHFSums_45_49/F");
  tree->Branch("minOfHFSums_45_50", &minOfHFSums_45_50, "minOfHFSums_45_50/F");
  tree->Branch("minOfHFSums_45_51", &minOfHFSums_45_51, "minOfHFSums_45_51/F");
  tree->Branch("minOfHFSums_45_52", &minOfHFSums_45_52, "minOfHFSums_45_52/F");

  tree->Branch("minOfHFSums_46_47", &minOfHFSums_46_47, "minOfHFSums_46_47/F");
  tree->Branch("minOfHFSums_46_48", &minOfHFSums_46_48, "minOfHFSums_46_48/F");
  tree->Branch("minOfHFSums_46_49", &minOfHFSums_46_49, "minOfHFSums_46_49/F");
  tree->Branch("minOfHFSums_46_50", &minOfHFSums_46_50, "minOfHFSums_46_50/F");
  tree->Branch("minOfHFSums_46_51", &minOfHFSums_46_51, "minOfHFSums_46_51/F");
  tree->Branch("minOfHFSums_46_52", &minOfHFSums_46_52, "minOfHFSums_46_52/F");

  tree->Branch("minOfHFSums_47_48", &minOfHFSums_47_48, "minOfHFSums_47_48/F");
  tree->Branch("minOfHFSums_47_49", &minOfHFSums_47_49, "minOfHFSums_47_49/F");
  tree->Branch("minOfHFSums_47_50", &minOfHFSums_47_50, "minOfHFSums_47_50/F");
  tree->Branch("minOfHFSums_47_51", &minOfHFSums_47_51, "minOfHFSums_47_51/F");
  tree->Branch("minOfHFSums_47_52", &minOfHFSums_47_52, "minOfHFSums_47_52/F");

  tree->Branch("minOfHFSums_48_49", &minOfHFSums_48_49, "minOfHFSums_48_49/F");
  tree->Branch("minOfHFSums_48_50", &minOfHFSums_48_50, "minOfHFSums_48_50/F");
  tree->Branch("minOfHFSums_48_51", &minOfHFSums_48_51, "minOfHFSums_48_51/F");
  tree->Branch("minOfHFSums_48_52", &minOfHFSums_48_52, "minOfHFSums_48_52/F");

  tree->Branch("minOfHFSums_49_50", &minOfHFSums_49_50, "minOfHFSums_49_50/F");
  tree->Branch("minOfHFSums_49_51", &minOfHFSums_49_51, "minOfHFSums_49_51/F");
  tree->Branch("minOfHFSums_49_52", &minOfHFSums_49_52, "minOfHFSums_49_52/F");

  tree->Branch("minOfHFSums_50_51", &minOfHFSums_50_51, "minOfHFSums_50_51/F");
  tree->Branch("minOfHFSums_50_52", &minOfHFSums_50_52, "minOfHFSums_50_52/F");

  tree->Branch("minOfHFSums_51_52", &minOfHFSums_51_52, "minOfHFSums_51_52/F");
  // New trees --end
    
  // Counters for cutflow
  const int Ncuts = 13;
  Int_t failCount[Ncuts] = {0};
  Int_t nPassed = 0;
  const char* cutNames[Ncuts] = {
    "rapidityGapPlus>0",
    "rapidityGapMinus>0",
    "hfSumPlus>=0",
    "hfSumMinus>=0",
    "0<=maxHFEnergy<100k",
    "0<=minHFEnergy<100k",
    ">=4 jets (pT>25,|eta|<2.4)",
    ">=2 b-tag jets",
    ">=2 light jets",
    "==1 lepton",
    "MET>20",
    "minOfMaxHF>=0",
    "minOfHFSums>=0"
  };

  // Event loop
  for (Long64_t i=0; i<nLoop; ++i) {
    reader->ReadEntry(offset + i);

    // --- 1) Tracks & nchPV (pT>0.4, |eta|<2.4) ---
    nchPV = 0;
    std::vector<double> trackEtas;
    Int_t nTrk = branchTrack->GetEntriesFast();
    trackEtas.reserve(nTrk);
    for (Int_t t=0; t<nTrk; ++t) {
      auto* tr = (Track*)branchTrack->At(t);
      if (tr->PT>0.4 && std::fabs(tr->Eta)<2.4) {
        nchPV++;
        trackEtas.push_back(tr->Eta);
      }
    }

    // --- 2) Jets ---
    nTotalJets    = branchJet->GetEntriesFast();
    nForwardJets  = 0;
    leadLightJetPt= 0;
    Int_t nJetsSel = 0, nBtagSel = 0, nLightSel = 0;
    std::vector<double> lightPts;
    lightPts.reserve(nTotalJets);
    for (Int_t j=0; j<nTotalJets; ++j) {
      auto* jet = (Jet*)branchJet->At(j);
      // count forward
      if (std::fabs(jet->Eta)>3.0) nForwardJets++;
      // leading light jet Pt
      if (jet->BTag==0) lightPts.push_back(jet->PT);
      // selection for multiplicity
      if (jet->PT>25.0 && std::fabs(jet->Eta)<2.4) {
        nJetsSel++;
        if (jet->BTag>0) nBtagSel++;
        else             nLightSel++;
      }
    }
    if (!lightPts.empty()) {
      leadLightJetPt = *std::max_element(lightPts.begin(), lightPts.end());
    }

    // --- 3) HF towers + fixed new variables ---
    hfSumPlus        = 0.0f;
    hfSumMinus       = 0.0f;
    // initialize per-hemisphere maxima to 0, not –1
    maxHFEnergyPlus  = 0.0f;
    maxHFEnergyMinus = 0.0f;
    // global HF extrema: sentinel for “none found”
    maxHFEnergy      = std::numeric_limits<float>::lowest();
    minHFEnergy      = std::numeric_limits<float>::max();

    // reset new
    maxHFEnergyPlus_30_31  = 0.0f;
    maxHFEnergyMinus_30_31 = 0.0f;
    maxHFEnergy_30_31      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_31      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_32  = 0.0f;
    maxHFEnergyMinus_30_32 = 0.0f;
    maxHFEnergy_30_32      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_32      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_33  = 0.0f;
    maxHFEnergyMinus_30_33 = 0.0f;
    maxHFEnergy_30_33      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_33      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_34  = 0.0f;
    maxHFEnergyMinus_30_34 = 0.0f;
    maxHFEnergy_30_34      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_34      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_35  = 0.0f;
    maxHFEnergyMinus_30_35 = 0.0f;
    maxHFEnergy_30_35      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_35      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_36  = 0.0f;
    maxHFEnergyMinus_30_36 = 0.0f;
    maxHFEnergy_30_36      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_36      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_37  = 0.0f;
    maxHFEnergyMinus_30_37 = 0.0f;
    maxHFEnergy_30_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_38  = 0.0f;
    maxHFEnergyMinus_30_38 = 0.0f;
    maxHFEnergy_30_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_39  = 0.0f;
    maxHFEnergyMinus_30_39 = 0.0f;
    maxHFEnergy_30_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_40  = 0.0f;
    maxHFEnergyMinus_30_40 = 0.0f;
    maxHFEnergy_30_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_41  = 0.0f;
    maxHFEnergyMinus_30_41 = 0.0f;
    maxHFEnergy_30_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_42  = 0.0f;
    maxHFEnergyMinus_30_42 = 0.0f;
    maxHFEnergy_30_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_43  = 0.0f;
    maxHFEnergyMinus_30_43 = 0.0f;
    maxHFEnergy_30_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_44  = 0.0f;
    maxHFEnergyMinus_30_44 = 0.0f;
    maxHFEnergy_30_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_45  = 0.0f;
    maxHFEnergyMinus_30_45 = 0.0f;
    maxHFEnergy_30_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_46  = 0.0f;
    maxHFEnergyMinus_30_46 = 0.0f;
    maxHFEnergy_30_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_47  = 0.0f;
    maxHFEnergyMinus_30_47 = 0.0f;
    maxHFEnergy_30_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_48  = 0.0f;
    maxHFEnergyMinus_30_48 = 0.0f;
    maxHFEnergy_30_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_49  = 0.0f;
    maxHFEnergyMinus_30_49 = 0.0f;
    maxHFEnergy_30_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_50  = 0.0f;
    maxHFEnergyMinus_30_50 = 0.0f;
    maxHFEnergy_30_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_51  = 0.0f;
    maxHFEnergyMinus_30_51 = 0.0f;
    maxHFEnergy_30_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_30_52  = 0.0f;
    maxHFEnergyMinus_30_52 = 0.0f;
    maxHFEnergy_30_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_30_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_31_32  = 0.0f;
    maxHFEnergyMinus_31_32 = 0.0f;
    maxHFEnergy_31_32      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_32      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_33  = 0.0f;
    maxHFEnergyMinus_31_33 = 0.0f;
    maxHFEnergy_31_33      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_33      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_34  = 0.0f;
    maxHFEnergyMinus_31_34 = 0.0f;
    maxHFEnergy_31_34      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_34      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_35  = 0.0f;
    maxHFEnergyMinus_31_35 = 0.0f;
    maxHFEnergy_31_35      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_35      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_36  = 0.0f;
    maxHFEnergyMinus_31_36 = 0.0f;
    maxHFEnergy_31_36      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_36      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_37  = 0.0f;
    maxHFEnergyMinus_31_37 = 0.0f;
    maxHFEnergy_31_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_38  = 0.0f;
    maxHFEnergyMinus_31_38 = 0.0f;
    maxHFEnergy_31_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_39  = 0.0f;
    maxHFEnergyMinus_31_39 = 0.0f;
    maxHFEnergy_31_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_40  = 0.0f;
    maxHFEnergyMinus_31_40 = 0.0f;
    maxHFEnergy_31_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_41  = 0.0f;
    maxHFEnergyMinus_31_41 = 0.0f;
    maxHFEnergy_31_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_42  = 0.0f;
    maxHFEnergyMinus_31_42 = 0.0f;
    maxHFEnergy_31_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_43  = 0.0f;
    maxHFEnergyMinus_31_43 = 0.0f;
    maxHFEnergy_31_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_44  = 0.0f;
    maxHFEnergyMinus_31_44 = 0.0f;
    maxHFEnergy_31_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_45  = 0.0f;
    maxHFEnergyMinus_31_45 = 0.0f;
    maxHFEnergy_31_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_46  = 0.0f;
    maxHFEnergyMinus_31_46 = 0.0f;
    maxHFEnergy_31_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_47  = 0.0f;
    maxHFEnergyMinus_31_47 = 0.0f;
    maxHFEnergy_31_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_48  = 0.0f;
    maxHFEnergyMinus_31_48 = 0.0f;
    maxHFEnergy_31_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_49  = 0.0f;
    maxHFEnergyMinus_31_49 = 0.0f;
    maxHFEnergy_31_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_50  = 0.0f;
    maxHFEnergyMinus_31_50 = 0.0f;
    maxHFEnergy_31_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_51  = 0.0f;
    maxHFEnergyMinus_31_51 = 0.0f;
    maxHFEnergy_31_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_31_52  = 0.0f;
    maxHFEnergyMinus_31_52 = 0.0f;
    maxHFEnergy_31_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_31_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_32_33  = 0.0f;
    maxHFEnergyMinus_32_33 = 0.0f;
    maxHFEnergy_32_33      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_33      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_34  = 0.0f;
    maxHFEnergyMinus_32_34 = 0.0f;
    maxHFEnergy_32_34      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_34      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_35  = 0.0f;
    maxHFEnergyMinus_32_35 = 0.0f;
    maxHFEnergy_32_35      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_35      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_36  = 0.0f;
    maxHFEnergyMinus_32_36 = 0.0f;
    maxHFEnergy_32_36      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_36      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_37  = 0.0f;
    maxHFEnergyMinus_32_37 = 0.0f;
    maxHFEnergy_32_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_38  = 0.0f;
    maxHFEnergyMinus_32_38 = 0.0f;
    maxHFEnergy_32_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_39  = 0.0f;
    maxHFEnergyMinus_32_39 = 0.0f;
    maxHFEnergy_32_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_40  = 0.0f;
    maxHFEnergyMinus_32_40 = 0.0f;
    maxHFEnergy_32_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_41  = 0.0f;
    maxHFEnergyMinus_32_41 = 0.0f;
    maxHFEnergy_32_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_42  = 0.0f;
    maxHFEnergyMinus_32_42 = 0.0f;
    maxHFEnergy_32_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_43  = 0.0f;
    maxHFEnergyMinus_32_43 = 0.0f;
    maxHFEnergy_32_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_44  = 0.0f;
    maxHFEnergyMinus_32_44 = 0.0f;
    maxHFEnergy_32_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_45  = 0.0f;
    maxHFEnergyMinus_32_45 = 0.0f;
    maxHFEnergy_32_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_46  = 0.0f;
    maxHFEnergyMinus_32_46 = 0.0f;
    maxHFEnergy_32_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_47  = 0.0f;
    maxHFEnergyMinus_32_47 = 0.0f;
    maxHFEnergy_32_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_48  = 0.0f;
    maxHFEnergyMinus_32_48 = 0.0f;
    maxHFEnergy_32_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_49  = 0.0f;
    maxHFEnergyMinus_32_49 = 0.0f;
    maxHFEnergy_32_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_50  = 0.0f;
    maxHFEnergyMinus_32_50 = 0.0f;
    maxHFEnergy_32_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_51  = 0.0f;
    maxHFEnergyMinus_32_51 = 0.0f;
    maxHFEnergy_32_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_32_52  = 0.0f;
    maxHFEnergyMinus_32_52 = 0.0f;
    maxHFEnergy_32_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_32_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_33_34  = 0.0f;
    maxHFEnergyMinus_33_34 = 0.0f;
    maxHFEnergy_33_34      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_34      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_35  = 0.0f;
    maxHFEnergyMinus_33_35 = 0.0f;
    maxHFEnergy_33_35      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_35      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_36  = 0.0f;
    maxHFEnergyMinus_33_36 = 0.0f;
    maxHFEnergy_33_36      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_36      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_37  = 0.0f;
    maxHFEnergyMinus_33_37 = 0.0f;
    maxHFEnergy_33_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_38  = 0.0f;
    maxHFEnergyMinus_33_38 = 0.0f;
    maxHFEnergy_33_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_39  = 0.0f;
    maxHFEnergyMinus_33_39 = 0.0f;
    maxHFEnergy_33_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_40  = 0.0f;
    maxHFEnergyMinus_33_40 = 0.0f;
    maxHFEnergy_33_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_41  = 0.0f;
    maxHFEnergyMinus_33_41 = 0.0f;
    maxHFEnergy_33_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_42  = 0.0f;
    maxHFEnergyMinus_33_42 = 0.0f;
    maxHFEnergy_33_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_43  = 0.0f;
    maxHFEnergyMinus_33_43 = 0.0f;
    maxHFEnergy_33_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_44  = 0.0f;
    maxHFEnergyMinus_33_44 = 0.0f;
    maxHFEnergy_33_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_45  = 0.0f;
    maxHFEnergyMinus_33_45 = 0.0f;
    maxHFEnergy_33_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_46  = 0.0f;
    maxHFEnergyMinus_33_46 = 0.0f;
    maxHFEnergy_33_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_47  = 0.0f;
    maxHFEnergyMinus_33_47 = 0.0f;
    maxHFEnergy_33_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_48  = 0.0f;
    maxHFEnergyMinus_33_48 = 0.0f;
    maxHFEnergy_33_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_49  = 0.0f;
    maxHFEnergyMinus_33_49 = 0.0f;
    maxHFEnergy_33_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_50  = 0.0f;
    maxHFEnergyMinus_33_50 = 0.0f;
    maxHFEnergy_33_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_51  = 0.0f;
    maxHFEnergyMinus_33_51 = 0.0f;
    maxHFEnergy_33_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_33_52  = 0.0f;
    maxHFEnergyMinus_33_52 = 0.0f;
    maxHFEnergy_33_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_33_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_34_35  = 0.0f;
    maxHFEnergyMinus_34_35 = 0.0f;
    maxHFEnergy_34_35      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_35      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_36  = 0.0f;
    maxHFEnergyMinus_34_36 = 0.0f;
    maxHFEnergy_34_36      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_36      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_37  = 0.0f;
    maxHFEnergyMinus_34_37 = 0.0f;
    maxHFEnergy_34_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_38  = 0.0f;
    maxHFEnergyMinus_34_38 = 0.0f;
    maxHFEnergy_34_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_39  = 0.0f;
    maxHFEnergyMinus_34_39 = 0.0f;
    maxHFEnergy_34_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_40  = 0.0f;
    maxHFEnergyMinus_34_40 = 0.0f;
    maxHFEnergy_34_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_41  = 0.0f;
    maxHFEnergyMinus_34_41 = 0.0f;
    maxHFEnergy_34_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_42  = 0.0f;
    maxHFEnergyMinus_34_42 = 0.0f;
    maxHFEnergy_34_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_43  = 0.0f;
    maxHFEnergyMinus_34_43 = 0.0f;
    maxHFEnergy_34_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_44  = 0.0f;
    maxHFEnergyMinus_34_44 = 0.0f;
    maxHFEnergy_34_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_45  = 0.0f;
    maxHFEnergyMinus_34_45 = 0.0f;
    maxHFEnergy_34_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_46  = 0.0f;
    maxHFEnergyMinus_34_46 = 0.0f;
    maxHFEnergy_34_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_47  = 0.0f;
    maxHFEnergyMinus_34_47 = 0.0f;
    maxHFEnergy_34_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_48  = 0.0f;
    maxHFEnergyMinus_34_48 = 0.0f;
    maxHFEnergy_34_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_49  = 0.0f;
    maxHFEnergyMinus_34_49 = 0.0f;
    maxHFEnergy_34_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_50  = 0.0f;
    maxHFEnergyMinus_34_50 = 0.0f;
    maxHFEnergy_34_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_51  = 0.0f;
    maxHFEnergyMinus_34_51 = 0.0f;
    maxHFEnergy_34_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_34_52  = 0.0f;
    maxHFEnergyMinus_34_52 = 0.0f;
    maxHFEnergy_34_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_34_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_35_36  = 0.0f;
    maxHFEnergyMinus_35_36 = 0.0f;
    maxHFEnergy_35_36      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_36      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_37  = 0.0f;
    maxHFEnergyMinus_35_37 = 0.0f;
    maxHFEnergy_35_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_38  = 0.0f;
    maxHFEnergyMinus_35_38 = 0.0f;
    maxHFEnergy_35_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_39  = 0.0f;
    maxHFEnergyMinus_35_39 = 0.0f;
    maxHFEnergy_35_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_40  = 0.0f;
    maxHFEnergyMinus_35_40 = 0.0f;
    maxHFEnergy_35_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_41  = 0.0f;
    maxHFEnergyMinus_35_41 = 0.0f;
    maxHFEnergy_35_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_42  = 0.0f;
    maxHFEnergyMinus_35_42 = 0.0f;
    maxHFEnergy_35_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_43  = 0.0f;
    maxHFEnergyMinus_35_43 = 0.0f;
    maxHFEnergy_35_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_44  = 0.0f;
    maxHFEnergyMinus_35_44 = 0.0f;
    maxHFEnergy_35_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_45  = 0.0f;
    maxHFEnergyMinus_35_45 = 0.0f;
    maxHFEnergy_35_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_46  = 0.0f;
    maxHFEnergyMinus_35_46 = 0.0f;
    maxHFEnergy_35_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_47  = 0.0f;
    maxHFEnergyMinus_35_47 = 0.0f;
    maxHFEnergy_35_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_48  = 0.0f;
    maxHFEnergyMinus_35_48 = 0.0f;
    maxHFEnergy_35_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_49  = 0.0f;
    maxHFEnergyMinus_35_49 = 0.0f;
    maxHFEnergy_35_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_50  = 0.0f;
    maxHFEnergyMinus_35_50 = 0.0f;
    maxHFEnergy_35_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_51  = 0.0f;
    maxHFEnergyMinus_35_51 = 0.0f;
    maxHFEnergy_35_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_35_52  = 0.0f;
    maxHFEnergyMinus_35_52 = 0.0f;
    maxHFEnergy_35_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_35_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_36_37  = 0.0f;
    maxHFEnergyMinus_36_37 = 0.0f;
    maxHFEnergy_36_37      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_37      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_38  = 0.0f;
    maxHFEnergyMinus_36_38 = 0.0f;
    maxHFEnergy_36_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_39  = 0.0f;
    maxHFEnergyMinus_36_39 = 0.0f;
    maxHFEnergy_36_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_40  = 0.0f;
    maxHFEnergyMinus_36_40 = 0.0f;
    maxHFEnergy_36_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_41  = 0.0f;
    maxHFEnergyMinus_36_41 = 0.0f;
    maxHFEnergy_36_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_42  = 0.0f;
    maxHFEnergyMinus_36_42 = 0.0f;
    maxHFEnergy_36_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_43  = 0.0f;
    maxHFEnergyMinus_36_43 = 0.0f;
    maxHFEnergy_36_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_44  = 0.0f;
    maxHFEnergyMinus_36_44 = 0.0f;
    maxHFEnergy_36_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_45  = 0.0f;
    maxHFEnergyMinus_36_45 = 0.0f;
    maxHFEnergy_36_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_46  = 0.0f;
    maxHFEnergyMinus_36_46 = 0.0f;
    maxHFEnergy_36_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_47  = 0.0f;
    maxHFEnergyMinus_36_47 = 0.0f;
    maxHFEnergy_36_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_48  = 0.0f;
    maxHFEnergyMinus_36_48 = 0.0f;
    maxHFEnergy_36_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_49  = 0.0f;
    maxHFEnergyMinus_36_49 = 0.0f;
    maxHFEnergy_36_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_50  = 0.0f;
    maxHFEnergyMinus_36_50 = 0.0f;
    maxHFEnergy_36_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_51  = 0.0f;
    maxHFEnergyMinus_36_51 = 0.0f;
    maxHFEnergy_36_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_36_52  = 0.0f;
    maxHFEnergyMinus_36_52 = 0.0f;
    maxHFEnergy_36_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_36_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_37_38  = 0.0f;
    maxHFEnergyMinus_37_38 = 0.0f;
    maxHFEnergy_37_38      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_38      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_39  = 0.0f;
    maxHFEnergyMinus_37_39 = 0.0f;
    maxHFEnergy_37_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_40  = 0.0f;
    maxHFEnergyMinus_37_40 = 0.0f;
    maxHFEnergy_37_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_41  = 0.0f;
    maxHFEnergyMinus_37_41 = 0.0f;
    maxHFEnergy_37_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_42  = 0.0f;
    maxHFEnergyMinus_37_42 = 0.0f;
    maxHFEnergy_37_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_43  = 0.0f;
    maxHFEnergyMinus_37_43 = 0.0f;
    maxHFEnergy_37_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_44  = 0.0f;
    maxHFEnergyMinus_37_44 = 0.0f;
    maxHFEnergy_37_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_45  = 0.0f;
    maxHFEnergyMinus_37_45 = 0.0f;
    maxHFEnergy_37_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_46  = 0.0f;
    maxHFEnergyMinus_37_46 = 0.0f;
    maxHFEnergy_37_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_47  = 0.0f;
    maxHFEnergyMinus_37_47 = 0.0f;
    maxHFEnergy_37_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_48  = 0.0f;
    maxHFEnergyMinus_37_48 = 0.0f;
    maxHFEnergy_37_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_49  = 0.0f;
    maxHFEnergyMinus_37_49 = 0.0f;
    maxHFEnergy_37_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_50  = 0.0f;
    maxHFEnergyMinus_37_50 = 0.0f;
    maxHFEnergy_37_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_51  = 0.0f;
    maxHFEnergyMinus_37_51 = 0.0f;
    maxHFEnergy_37_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_37_52  = 0.0f;
    maxHFEnergyMinus_37_52 = 0.0f;
    maxHFEnergy_37_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_37_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_38_39  = 0.0f;
    maxHFEnergyMinus_38_39 = 0.0f;
    maxHFEnergy_38_39      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_39      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_40  = 0.0f;
    maxHFEnergyMinus_38_40 = 0.0f;
    maxHFEnergy_38_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_41  = 0.0f;
    maxHFEnergyMinus_38_41 = 0.0f;
    maxHFEnergy_38_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_42  = 0.0f;
    maxHFEnergyMinus_38_42 = 0.0f;
    maxHFEnergy_38_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_43  = 0.0f;
    maxHFEnergyMinus_38_43 = 0.0f;
    maxHFEnergy_38_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_44  = 0.0f;
    maxHFEnergyMinus_38_44 = 0.0f;
    maxHFEnergy_38_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_45  = 0.0f;
    maxHFEnergyMinus_38_45 = 0.0f;
    maxHFEnergy_38_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_46  = 0.0f;
    maxHFEnergyMinus_38_46 = 0.0f;
    maxHFEnergy_38_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_47  = 0.0f;
    maxHFEnergyMinus_38_47 = 0.0f;
    maxHFEnergy_38_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_48  = 0.0f;
    maxHFEnergyMinus_38_48 = 0.0f;
    maxHFEnergy_38_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_49  = 0.0f;
    maxHFEnergyMinus_38_49 = 0.0f;
    maxHFEnergy_38_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_50  = 0.0f;
    maxHFEnergyMinus_38_50 = 0.0f;
    maxHFEnergy_38_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_51  = 0.0f;
    maxHFEnergyMinus_38_51 = 0.0f;
    maxHFEnergy_38_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_38_52  = 0.0f;
    maxHFEnergyMinus_38_52 = 0.0f;
    maxHFEnergy_38_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_38_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_39_40  = 0.0f;
    maxHFEnergyMinus_39_40 = 0.0f;
    maxHFEnergy_39_40      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_40      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_41  = 0.0f;
    maxHFEnergyMinus_39_41 = 0.0f;
    maxHFEnergy_39_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_42  = 0.0f;
    maxHFEnergyMinus_39_42 = 0.0f;
    maxHFEnergy_39_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_43  = 0.0f;
    maxHFEnergyMinus_39_43 = 0.0f;
    maxHFEnergy_39_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_44  = 0.0f;
    maxHFEnergyMinus_39_44 = 0.0f;
    maxHFEnergy_39_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_45  = 0.0f;
    maxHFEnergyMinus_39_45 = 0.0f;
    maxHFEnergy_39_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_46  = 0.0f;
    maxHFEnergyMinus_39_46 = 0.0f;
    maxHFEnergy_39_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_47  = 0.0f;
    maxHFEnergyMinus_39_47 = 0.0f;
    maxHFEnergy_39_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_48  = 0.0f;
    maxHFEnergyMinus_39_48 = 0.0f;
    maxHFEnergy_39_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_49  = 0.0f;
    maxHFEnergyMinus_39_49 = 0.0f;
    maxHFEnergy_39_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_50  = 0.0f;
    maxHFEnergyMinus_39_50 = 0.0f;
    maxHFEnergy_39_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_51  = 0.0f;
    maxHFEnergyMinus_39_51 = 0.0f;
    maxHFEnergy_39_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_39_52  = 0.0f;
    maxHFEnergyMinus_39_52 = 0.0f;
    maxHFEnergy_39_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_39_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_40_41  = 0.0f;
    maxHFEnergyMinus_40_41 = 0.0f;
    maxHFEnergy_40_41      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_41      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_42  = 0.0f;
    maxHFEnergyMinus_40_42 = 0.0f;
    maxHFEnergy_40_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_43  = 0.0f;
    maxHFEnergyMinus_40_43 = 0.0f;
    maxHFEnergy_40_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_44  = 0.0f;
    maxHFEnergyMinus_40_44 = 0.0f;
    maxHFEnergy_40_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_45  = 0.0f;
    maxHFEnergyMinus_40_45 = 0.0f;
    maxHFEnergy_40_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_46  = 0.0f;
    maxHFEnergyMinus_40_46 = 0.0f;
    maxHFEnergy_40_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_47  = 0.0f;
    maxHFEnergyMinus_40_47 = 0.0f;
    maxHFEnergy_40_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_48  = 0.0f;
    maxHFEnergyMinus_40_48 = 0.0f;
    maxHFEnergy_40_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_49  = 0.0f;
    maxHFEnergyMinus_40_49 = 0.0f;
    maxHFEnergy_40_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_50  = 0.0f;
    maxHFEnergyMinus_40_50 = 0.0f;
    maxHFEnergy_40_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_51  = 0.0f;
    maxHFEnergyMinus_40_51 = 0.0f;
    maxHFEnergy_40_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_40_52  = 0.0f;
    maxHFEnergyMinus_40_52 = 0.0f;
    maxHFEnergy_40_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_40_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_41_42  = 0.0f;
    maxHFEnergyMinus_41_42 = 0.0f;
    maxHFEnergy_41_42      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_42      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_43  = 0.0f;
    maxHFEnergyMinus_41_43 = 0.0f;
    maxHFEnergy_41_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_44  = 0.0f;
    maxHFEnergyMinus_41_44 = 0.0f;
    maxHFEnergy_41_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_45  = 0.0f;
    maxHFEnergyMinus_41_45 = 0.0f;
    maxHFEnergy_41_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_46  = 0.0f;
    maxHFEnergyMinus_41_46 = 0.0f;
    maxHFEnergy_41_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_47  = 0.0f;
    maxHFEnergyMinus_41_47 = 0.0f;
    maxHFEnergy_41_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_48  = 0.0f;
    maxHFEnergyMinus_41_48 = 0.0f;
    maxHFEnergy_41_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_49  = 0.0f;
    maxHFEnergyMinus_41_49 = 0.0f;
    maxHFEnergy_41_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_50  = 0.0f;
    maxHFEnergyMinus_41_50 = 0.0f;
    maxHFEnergy_41_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_51  = 0.0f;
    maxHFEnergyMinus_41_51 = 0.0f;
    maxHFEnergy_41_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_41_52  = 0.0f;
    maxHFEnergyMinus_41_52 = 0.0f;
    maxHFEnergy_41_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_41_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_42_43  = 0.0f;
    maxHFEnergyMinus_42_43 = 0.0f;
    maxHFEnergy_42_43      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_43      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_44  = 0.0f;
    maxHFEnergyMinus_42_44 = 0.0f;
    maxHFEnergy_42_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_45  = 0.0f;
    maxHFEnergyMinus_42_45 = 0.0f;
    maxHFEnergy_42_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_46  = 0.0f;
    maxHFEnergyMinus_42_46 = 0.0f;
    maxHFEnergy_42_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_47  = 0.0f;
    maxHFEnergyMinus_42_47 = 0.0f;
    maxHFEnergy_42_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_48  = 0.0f;
    maxHFEnergyMinus_42_48 = 0.0f;
    maxHFEnergy_42_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_49  = 0.0f;
    maxHFEnergyMinus_42_49 = 0.0f;
    maxHFEnergy_42_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_50  = 0.0f;
    maxHFEnergyMinus_42_50 = 0.0f;
    maxHFEnergy_42_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_51  = 0.0f;
    maxHFEnergyMinus_42_51 = 0.0f;
    maxHFEnergy_42_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_42_52  = 0.0f;
    maxHFEnergyMinus_42_52 = 0.0f;
    maxHFEnergy_42_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_42_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_43_44  = 0.0f;
    maxHFEnergyMinus_43_44 = 0.0f;
    maxHFEnergy_43_44      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_44      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_45  = 0.0f;
    maxHFEnergyMinus_43_45 = 0.0f;
    maxHFEnergy_43_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_46  = 0.0f;
    maxHFEnergyMinus_43_46 = 0.0f;
    maxHFEnergy_43_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_47  = 0.0f;
    maxHFEnergyMinus_43_47 = 0.0f;
    maxHFEnergy_43_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_48  = 0.0f;
    maxHFEnergyMinus_43_48 = 0.0f;
    maxHFEnergy_43_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_49  = 0.0f;
    maxHFEnergyMinus_43_49 = 0.0f;
    maxHFEnergy_43_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_50  = 0.0f;
    maxHFEnergyMinus_43_50 = 0.0f;
    maxHFEnergy_43_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_51  = 0.0f;
    maxHFEnergyMinus_43_51 = 0.0f;
    maxHFEnergy_43_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_43_52  = 0.0f;
    maxHFEnergyMinus_43_52 = 0.0f;
    maxHFEnergy_43_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_43_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_44_45  = 0.0f;
    maxHFEnergyMinus_44_45 = 0.0f;
    maxHFEnergy_44_45      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_45      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_46  = 0.0f;
    maxHFEnergyMinus_44_46 = 0.0f;
    maxHFEnergy_44_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_47  = 0.0f;
    maxHFEnergyMinus_44_47 = 0.0f;
    maxHFEnergy_44_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_48  = 0.0f;
    maxHFEnergyMinus_44_48 = 0.0f;
    maxHFEnergy_44_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_49  = 0.0f;
    maxHFEnergyMinus_44_49 = 0.0f;
    maxHFEnergy_44_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_50  = 0.0f;
    maxHFEnergyMinus_44_50 = 0.0f;
    maxHFEnergy_44_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_51  = 0.0f;
    maxHFEnergyMinus_44_51 = 0.0f;
    maxHFEnergy_44_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_44_52  = 0.0f;
    maxHFEnergyMinus_44_52 = 0.0f;
    maxHFEnergy_44_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_44_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_45_46  = 0.0f;
    maxHFEnergyMinus_45_46 = 0.0f;
    maxHFEnergy_45_46      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_46      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_45_47  = 0.0f;
    maxHFEnergyMinus_45_47 = 0.0f;
    maxHFEnergy_45_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_45_48  = 0.0f;
    maxHFEnergyMinus_45_48 = 0.0f;
    maxHFEnergy_45_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_45_49  = 0.0f;
    maxHFEnergyMinus_45_49 = 0.0f;
    maxHFEnergy_45_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_45_50  = 0.0f;
    maxHFEnergyMinus_45_50 = 0.0f;
    maxHFEnergy_45_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_45_51  = 0.0f;
    maxHFEnergyMinus_45_51 = 0.0f;
    maxHFEnergy_45_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_45_52  = 0.0f;
    maxHFEnergyMinus_45_52 = 0.0f;
    maxHFEnergy_45_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_45_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_46_47  = 0.0f;
    maxHFEnergyMinus_46_47 = 0.0f;
    maxHFEnergy_46_47      = std::numeric_limits<float>::lowest();
    minHFEnergy_46_47      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_46_48  = 0.0f;
    maxHFEnergyMinus_46_48 = 0.0f;
    maxHFEnergy_46_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_46_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_46_49  = 0.0f;
    maxHFEnergyMinus_46_49 = 0.0f;
    maxHFEnergy_46_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_46_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_46_50  = 0.0f;
    maxHFEnergyMinus_46_50 = 0.0f;
    maxHFEnergy_46_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_46_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_46_51  = 0.0f;
    maxHFEnergyMinus_46_51 = 0.0f;
    maxHFEnergy_46_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_46_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_46_52  = 0.0f;
    maxHFEnergyMinus_46_52 = 0.0f;
    maxHFEnergy_46_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_46_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_47_48  = 0.0f;
    maxHFEnergyMinus_47_48 = 0.0f;
    maxHFEnergy_47_48      = std::numeric_limits<float>::lowest();
    minHFEnergy_47_48      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_47_49  = 0.0f;
    maxHFEnergyMinus_47_49 = 0.0f;
    maxHFEnergy_47_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_47_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_47_50  = 0.0f;
    maxHFEnergyMinus_47_50 = 0.0f;
    maxHFEnergy_47_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_47_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_47_51  = 0.0f;
    maxHFEnergyMinus_47_51 = 0.0f;
    maxHFEnergy_47_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_47_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_47_52  = 0.0f;
    maxHFEnergyMinus_47_52 = 0.0f;
    maxHFEnergy_47_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_47_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_48_49  = 0.0f;
    maxHFEnergyMinus_48_49 = 0.0f;
    maxHFEnergy_48_49      = std::numeric_limits<float>::lowest();
    minHFEnergy_48_49      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_48_50  = 0.0f;
    maxHFEnergyMinus_48_50 = 0.0f;
    maxHFEnergy_48_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_48_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_48_51  = 0.0f;
    maxHFEnergyMinus_48_51 = 0.0f;
    maxHFEnergy_48_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_48_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_48_52  = 0.0f;
    maxHFEnergyMinus_48_52 = 0.0f;
    maxHFEnergy_48_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_48_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_49_50  = 0.0f;
    maxHFEnergyMinus_49_50 = 0.0f;
    maxHFEnergy_49_50      = std::numeric_limits<float>::lowest();
    minHFEnergy_49_50      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_49_51  = 0.0f;
    maxHFEnergyMinus_49_51 = 0.0f;
    maxHFEnergy_49_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_49_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_49_52  = 0.0f;
    maxHFEnergyMinus_49_52 = 0.0f;
    maxHFEnergy_49_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_49_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_50_51  = 0.0f;
    maxHFEnergyMinus_50_51 = 0.0f;
    maxHFEnergy_50_51      = std::numeric_limits<float>::lowest();
    minHFEnergy_50_51      = std::numeric_limits<float>::max();
    maxHFEnergyPlus_50_52  = 0.0f;
    maxHFEnergyMinus_50_52 = 0.0f;
    maxHFEnergy_50_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_50_52      = std::numeric_limits<float>::max();


    maxHFEnergyPlus_51_52  = 0.0f;
    maxHFEnergyMinus_51_52 = 0.0f;
    maxHFEnergy_51_52      = std::numeric_limits<float>::lowest();
    minHFEnergy_51_52      = std::numeric_limits<float>::max();
    // reset new

    // reset
    hfSumPlus_30_31 = hfSumMinus_30_31 = 0.0f;
    hfSumPlus_30_32 = hfSumMinus_30_32 = 0.0f;
    hfSumPlus_30_33 = hfSumMinus_30_33 = 0.0f;
    hfSumPlus_30_34 = hfSumMinus_30_34 = 0.0f;
    hfSumPlus_30_35 = hfSumMinus_30_35 = 0.0f;
    hfSumPlus_30_36 = hfSumMinus_30_36 = 0.0f;
    hfSumPlus_30_37 = hfSumMinus_30_37 = 0.0f;
    hfSumPlus_30_38 = hfSumMinus_30_38 = 0.0f;
    hfSumPlus_30_39 = hfSumMinus_30_39 = 0.0f;
    hfSumPlus_30_40 = hfSumMinus_30_40 = 0.0f;
    hfSumPlus_30_41 = hfSumMinus_30_41 = 0.0f;
    hfSumPlus_30_42 = hfSumMinus_30_42 = 0.0f;
    hfSumPlus_30_43 = hfSumMinus_30_43 = 0.0f;
    hfSumPlus_30_44 = hfSumMinus_30_44 = 0.0f;
    hfSumPlus_30_45 = hfSumMinus_30_45 = 0.0f;
    hfSumPlus_30_46 = hfSumMinus_30_46 = 0.0f;
    hfSumPlus_30_47 = hfSumMinus_30_47 = 0.0f;
    hfSumPlus_30_48 = hfSumMinus_30_48 = 0.0f;
    hfSumPlus_30_49 = hfSumMinus_30_49 = 0.0f;
    hfSumPlus_30_50 = hfSumMinus_30_50 = 0.0f;
    hfSumPlus_30_51 = hfSumMinus_30_51 = 0.0f;
    hfSumPlus_30_52 = hfSumMinus_30_52 = 0.0f;

    hfSumPlus_31_32 = hfSumMinus_31_32 = 0.0f;
    hfSumPlus_31_33 = hfSumMinus_31_33 = 0.0f;
    hfSumPlus_31_34 = hfSumMinus_31_34 = 0.0f;
    hfSumPlus_31_35 = hfSumMinus_31_35 = 0.0f;
    hfSumPlus_31_36 = hfSumMinus_31_36 = 0.0f;
    hfSumPlus_31_37 = hfSumMinus_31_37 = 0.0f;
    hfSumPlus_31_38 = hfSumMinus_31_38 = 0.0f;
    hfSumPlus_31_39 = hfSumMinus_31_39 = 0.0f;
    hfSumPlus_31_40 = hfSumMinus_31_40 = 0.0f;
    hfSumPlus_31_41 = hfSumMinus_31_41 = 0.0f;
    hfSumPlus_31_42 = hfSumMinus_31_42 = 0.0f;
    hfSumPlus_31_43 = hfSumMinus_31_43 = 0.0f;
    hfSumPlus_31_44 = hfSumMinus_31_44 = 0.0f;
    hfSumPlus_31_45 = hfSumMinus_31_45 = 0.0f;
    hfSumPlus_31_46 = hfSumMinus_31_46 = 0.0f;
    hfSumPlus_31_47 = hfSumMinus_31_47 = 0.0f;
    hfSumPlus_31_48 = hfSumMinus_31_48 = 0.0f;
    hfSumPlus_31_49 = hfSumMinus_31_49 = 0.0f;
    hfSumPlus_31_50 = hfSumMinus_31_50 = 0.0f;
    hfSumPlus_31_51 = hfSumMinus_31_51 = 0.0f;
    hfSumPlus_31_52 = hfSumMinus_31_52 = 0.0f;

    hfSumPlus_32_33 = hfSumMinus_32_33 = 0.0f;
    hfSumPlus_32_34 = hfSumMinus_32_34 = 0.0f;
    hfSumPlus_32_35 = hfSumMinus_32_35 = 0.0f;
    hfSumPlus_32_36 = hfSumMinus_32_36 = 0.0f;
    hfSumPlus_32_37 = hfSumMinus_32_37 = 0.0f;
    hfSumPlus_32_38 = hfSumMinus_32_38 = 0.0f;
    hfSumPlus_32_39 = hfSumMinus_32_39 = 0.0f;
    hfSumPlus_32_40 = hfSumMinus_32_40 = 0.0f;
    hfSumPlus_32_41 = hfSumMinus_32_41 = 0.0f;
    hfSumPlus_32_42 = hfSumMinus_32_42 = 0.0f;
    hfSumPlus_32_43 = hfSumMinus_32_43 = 0.0f;
    hfSumPlus_32_44 = hfSumMinus_32_44 = 0.0f;
    hfSumPlus_32_45 = hfSumMinus_32_45 = 0.0f;
    hfSumPlus_32_46 = hfSumMinus_32_46 = 0.0f;
    hfSumPlus_32_47 = hfSumMinus_32_47 = 0.0f;
    hfSumPlus_32_48 = hfSumMinus_32_48 = 0.0f;
    hfSumPlus_32_49 = hfSumMinus_32_49 = 0.0f;
    hfSumPlus_32_50 = hfSumMinus_32_50 = 0.0f;
    hfSumPlus_32_51 = hfSumMinus_32_51 = 0.0f;
    hfSumPlus_32_52 = hfSumMinus_32_52 = 0.0f;

    hfSumPlus_33_34 = hfSumMinus_33_34 = 0.0f;
    hfSumPlus_33_35 = hfSumMinus_33_35 = 0.0f;
    hfSumPlus_33_36 = hfSumMinus_33_36 = 0.0f;
    hfSumPlus_33_37 = hfSumMinus_33_37 = 0.0f;
    hfSumPlus_33_38 = hfSumMinus_33_38 = 0.0f;
    hfSumPlus_33_39 = hfSumMinus_33_39 = 0.0f;
    hfSumPlus_33_40 = hfSumMinus_33_40 = 0.0f;
    hfSumPlus_33_41 = hfSumMinus_33_41 = 0.0f;
    hfSumPlus_33_42 = hfSumMinus_33_42 = 0.0f;
    hfSumPlus_33_43 = hfSumMinus_33_43 = 0.0f;
    hfSumPlus_33_44 = hfSumMinus_33_44 = 0.0f;
    hfSumPlus_33_45 = hfSumMinus_33_45 = 0.0f;
    hfSumPlus_33_46 = hfSumMinus_33_46 = 0.0f;
    hfSumPlus_33_47 = hfSumMinus_33_47 = 0.0f;
    hfSumPlus_33_48 = hfSumMinus_33_48 = 0.0f;
    hfSumPlus_33_49 = hfSumMinus_33_49 = 0.0f;
    hfSumPlus_33_50 = hfSumMinus_33_50 = 0.0f;
    hfSumPlus_33_51 = hfSumMinus_33_51 = 0.0f;
    hfSumPlus_33_52 = hfSumMinus_33_52 = 0.0f;

    hfSumPlus_34_35 = hfSumMinus_34_35 = 0.0f;
    hfSumPlus_34_36 = hfSumMinus_34_36 = 0.0f;
    hfSumPlus_34_37 = hfSumMinus_34_37 = 0.0f;
    hfSumPlus_34_38 = hfSumMinus_34_38 = 0.0f;
    hfSumPlus_34_39 = hfSumMinus_34_39 = 0.0f;
    hfSumPlus_34_40 = hfSumMinus_34_40 = 0.0f;
    hfSumPlus_34_41 = hfSumMinus_34_41 = 0.0f;
    hfSumPlus_34_42 = hfSumMinus_34_42 = 0.0f;
    hfSumPlus_34_43 = hfSumMinus_34_43 = 0.0f;
    hfSumPlus_34_44 = hfSumMinus_34_44 = 0.0f;
    hfSumPlus_34_45 = hfSumMinus_34_45 = 0.0f;
    hfSumPlus_34_46 = hfSumMinus_34_46 = 0.0f;
    hfSumPlus_34_47 = hfSumMinus_34_47 = 0.0f;
    hfSumPlus_34_48 = hfSumMinus_34_48 = 0.0f;
    hfSumPlus_34_49 = hfSumMinus_34_49 = 0.0f;
    hfSumPlus_34_50 = hfSumMinus_34_50 = 0.0f;
    hfSumPlus_34_51 = hfSumMinus_34_51 = 0.0f;
    hfSumPlus_34_52 = hfSumMinus_34_52 = 0.0f;

    hfSumPlus_35_36 = hfSumMinus_35_36 = 0.0f;
    hfSumPlus_35_37 = hfSumMinus_35_37 = 0.0f;
    hfSumPlus_35_38 = hfSumMinus_35_38 = 0.0f;
    hfSumPlus_35_39 = hfSumMinus_35_39 = 0.0f;
    hfSumPlus_35_40 = hfSumMinus_35_40 = 0.0f;
    hfSumPlus_35_41 = hfSumMinus_35_41 = 0.0f;
    hfSumPlus_35_42 = hfSumMinus_35_42 = 0.0f;
    hfSumPlus_35_43 = hfSumMinus_35_43 = 0.0f;
    hfSumPlus_35_44 = hfSumMinus_35_44 = 0.0f;
    hfSumPlus_35_45 = hfSumMinus_35_45 = 0.0f;
    hfSumPlus_35_46 = hfSumMinus_35_46 = 0.0f;
    hfSumPlus_35_47 = hfSumMinus_35_47 = 0.0f;
    hfSumPlus_35_48 = hfSumMinus_35_48 = 0.0f;
    hfSumPlus_35_49 = hfSumMinus_35_49 = 0.0f;
    hfSumPlus_35_50 = hfSumMinus_35_50 = 0.0f;
    hfSumPlus_35_51 = hfSumMinus_35_51 = 0.0f;
    hfSumPlus_35_52 = hfSumMinus_35_52 = 0.0f;

    hfSumPlus_36_37 = hfSumMinus_36_37 = 0.0f;
    hfSumPlus_36_38 = hfSumMinus_36_38 = 0.0f;
    hfSumPlus_36_39 = hfSumMinus_36_39 = 0.0f;
    hfSumPlus_36_40 = hfSumMinus_36_40 = 0.0f;
    hfSumPlus_36_41 = hfSumMinus_36_41 = 0.0f;
    hfSumPlus_36_42 = hfSumMinus_36_42 = 0.0f;
    hfSumPlus_36_43 = hfSumMinus_36_43 = 0.0f;
    hfSumPlus_36_44 = hfSumMinus_36_44 = 0.0f;
    hfSumPlus_36_45 = hfSumMinus_36_45 = 0.0f;
    hfSumPlus_36_46 = hfSumMinus_36_46 = 0.0f;
    hfSumPlus_36_47 = hfSumMinus_36_47 = 0.0f;
    hfSumPlus_36_48 = hfSumMinus_36_48 = 0.0f;
    hfSumPlus_36_49 = hfSumMinus_36_49 = 0.0f;
    hfSumPlus_36_50 = hfSumMinus_36_50 = 0.0f;
    hfSumPlus_36_51 = hfSumMinus_36_51 = 0.0f;
    hfSumPlus_36_52 = hfSumMinus_36_52 = 0.0f;

    hfSumPlus_37_38 = hfSumMinus_37_38 = 0.0f;
    hfSumPlus_37_39 = hfSumMinus_37_39 = 0.0f;
    hfSumPlus_37_40 = hfSumMinus_37_40 = 0.0f;
    hfSumPlus_37_41 = hfSumMinus_37_41 = 0.0f;
    hfSumPlus_37_42 = hfSumMinus_37_42 = 0.0f;
    hfSumPlus_37_43 = hfSumMinus_37_43 = 0.0f;
    hfSumPlus_37_44 = hfSumMinus_37_44 = 0.0f;
    hfSumPlus_37_45 = hfSumMinus_37_45 = 0.0f;
    hfSumPlus_37_46 = hfSumMinus_37_46 = 0.0f;
    hfSumPlus_37_47 = hfSumMinus_37_47 = 0.0f;
    hfSumPlus_37_48 = hfSumMinus_37_48 = 0.0f;
    hfSumPlus_37_49 = hfSumMinus_37_49 = 0.0f;
    hfSumPlus_37_50 = hfSumMinus_37_50 = 0.0f;
    hfSumPlus_37_51 = hfSumMinus_37_51 = 0.0f;
    hfSumPlus_37_52 = hfSumMinus_37_52 = 0.0f;

    hfSumPlus_38_39 = hfSumMinus_38_39 = 0.0f;
    hfSumPlus_38_40 = hfSumMinus_38_40 = 0.0f;
    hfSumPlus_38_41 = hfSumMinus_38_41 = 0.0f;
    hfSumPlus_38_42 = hfSumMinus_38_42 = 0.0f;
    hfSumPlus_38_43 = hfSumMinus_38_43 = 0.0f;
    hfSumPlus_38_44 = hfSumMinus_38_44 = 0.0f;
    hfSumPlus_38_45 = hfSumMinus_38_45 = 0.0f;
    hfSumPlus_38_46 = hfSumMinus_38_46 = 0.0f;
    hfSumPlus_38_47 = hfSumMinus_38_47 = 0.0f;
    hfSumPlus_38_48 = hfSumMinus_38_48 = 0.0f;
    hfSumPlus_38_49 = hfSumMinus_38_49 = 0.0f;
    hfSumPlus_38_50 = hfSumMinus_38_50 = 0.0f;
    hfSumPlus_38_51 = hfSumMinus_38_51 = 0.0f;
    hfSumPlus_38_52 = hfSumMinus_38_52 = 0.0f;

    hfSumPlus_39_40 = hfSumMinus_39_40 = 0.0f;
    hfSumPlus_39_41 = hfSumMinus_39_41 = 0.0f;
    hfSumPlus_39_42 = hfSumMinus_39_42 = 0.0f;
    hfSumPlus_39_43 = hfSumMinus_39_43 = 0.0f;
    hfSumPlus_39_44 = hfSumMinus_39_44 = 0.0f;
    hfSumPlus_39_45 = hfSumMinus_39_45 = 0.0f;
    hfSumPlus_39_46 = hfSumMinus_39_46 = 0.0f;
    hfSumPlus_39_47 = hfSumMinus_39_47 = 0.0f;
    hfSumPlus_39_48 = hfSumMinus_39_48 = 0.0f;
    hfSumPlus_39_49 = hfSumMinus_39_49 = 0.0f;
    hfSumPlus_39_50 = hfSumMinus_39_50 = 0.0f;
    hfSumPlus_39_51 = hfSumMinus_39_51 = 0.0f;
    hfSumPlus_39_52 = hfSumMinus_39_52 = 0.0f;

    hfSumPlus_40_41 = hfSumMinus_40_41 = 0.0f;
    hfSumPlus_40_42 = hfSumMinus_40_42 = 0.0f;
    hfSumPlus_40_43 = hfSumMinus_40_43 = 0.0f;
    hfSumPlus_40_44 = hfSumMinus_40_44 = 0.0f;
    hfSumPlus_40_45 = hfSumMinus_40_45 = 0.0f;
    hfSumPlus_40_46 = hfSumMinus_40_46 = 0.0f;
    hfSumPlus_40_47 = hfSumMinus_40_47 = 0.0f;
    hfSumPlus_40_48 = hfSumMinus_40_48 = 0.0f;
    hfSumPlus_40_49 = hfSumMinus_40_49 = 0.0f;
    hfSumPlus_40_50 = hfSumMinus_40_50 = 0.0f;
    hfSumPlus_40_51 = hfSumMinus_40_51 = 0.0f;
    hfSumPlus_40_52 = hfSumMinus_40_52 = 0.0f;

    hfSumPlus_41_42 = hfSumMinus_41_42 = 0.0f;
    hfSumPlus_41_43 = hfSumMinus_41_43 = 0.0f;
    hfSumPlus_41_44 = hfSumMinus_41_44 = 0.0f;
    hfSumPlus_41_45 = hfSumMinus_41_45 = 0.0f;
    hfSumPlus_41_46 = hfSumMinus_41_46 = 0.0f;
    hfSumPlus_41_47 = hfSumMinus_41_47 = 0.0f;
    hfSumPlus_41_48 = hfSumMinus_41_48 = 0.0f;
    hfSumPlus_41_49 = hfSumMinus_41_49 = 0.0f;
    hfSumPlus_41_50 = hfSumMinus_41_50 = 0.0f;
    hfSumPlus_41_51 = hfSumMinus_41_51 = 0.0f;
    hfSumPlus_41_52 = hfSumMinus_41_52 = 0.0f;

    hfSumPlus_42_43 = hfSumMinus_42_43 = 0.0f;
    hfSumPlus_42_44 = hfSumMinus_42_44 = 0.0f;
    hfSumPlus_42_45 = hfSumMinus_42_45 = 0.0f;
    hfSumPlus_42_46 = hfSumMinus_42_46 = 0.0f;
    hfSumPlus_42_47 = hfSumMinus_42_47 = 0.0f;
    hfSumPlus_42_48 = hfSumMinus_42_48 = 0.0f;
    hfSumPlus_42_49 = hfSumMinus_42_49 = 0.0f;
    hfSumPlus_42_50 = hfSumMinus_42_50 = 0.0f;
    hfSumPlus_42_51 = hfSumMinus_42_51 = 0.0f;
    hfSumPlus_42_52 = hfSumMinus_42_52 = 0.0f;

    hfSumPlus_43_44 = hfSumMinus_43_44 = 0.0f;
    hfSumPlus_43_45 = hfSumMinus_43_45 = 0.0f;
    hfSumPlus_43_46 = hfSumMinus_43_46 = 0.0f;
    hfSumPlus_43_47 = hfSumMinus_43_47 = 0.0f;
    hfSumPlus_43_48 = hfSumMinus_43_48 = 0.0f;
    hfSumPlus_43_49 = hfSumMinus_43_49 = 0.0f;
    hfSumPlus_43_50 = hfSumMinus_43_50 = 0.0f;
    hfSumPlus_43_51 = hfSumMinus_43_51 = 0.0f;
    hfSumPlus_43_52 = hfSumMinus_43_52 = 0.0f;

    hfSumPlus_44_45 = hfSumMinus_44_45 = 0.0f;
    hfSumPlus_44_46 = hfSumMinus_44_46 = 0.0f;
    hfSumPlus_44_47 = hfSumMinus_44_47 = 0.0f;
    hfSumPlus_44_48 = hfSumMinus_44_48 = 0.0f;
    hfSumPlus_44_49 = hfSumMinus_44_49 = 0.0f;
    hfSumPlus_44_50 = hfSumMinus_44_50 = 0.0f;
    hfSumPlus_44_51 = hfSumMinus_44_51 = 0.0f;
    hfSumPlus_44_52 = hfSumMinus_44_52 = 0.0f;

    hfSumPlus_45_46 = hfSumMinus_45_46 = 0.0f;
    hfSumPlus_45_47 = hfSumMinus_45_47 = 0.0f;
    hfSumPlus_45_48 = hfSumMinus_45_48 = 0.0f;
    hfSumPlus_45_49 = hfSumMinus_45_49 = 0.0f;
    hfSumPlus_45_50 = hfSumMinus_45_50 = 0.0f;
    hfSumPlus_45_51 = hfSumMinus_45_51 = 0.0f;
    hfSumPlus_45_52 = hfSumMinus_45_52 = 0.0f;

    hfSumPlus_46_47 = hfSumMinus_46_47 = 0.0f;
    hfSumPlus_46_48 = hfSumMinus_46_48 = 0.0f;
    hfSumPlus_46_49 = hfSumMinus_46_49 = 0.0f;
    hfSumPlus_46_50 = hfSumMinus_46_50 = 0.0f;
    hfSumPlus_46_51 = hfSumMinus_46_51 = 0.0f;
    hfSumPlus_46_52 = hfSumMinus_46_52 = 0.0f;

    hfSumPlus_47_48 = hfSumMinus_47_48 = 0.0f;
    hfSumPlus_47_49 = hfSumMinus_47_49 = 0.0f;
    hfSumPlus_47_50 = hfSumMinus_47_50 = 0.0f;
    hfSumPlus_47_51 = hfSumMinus_47_51 = 0.0f;
    hfSumPlus_47_52 = hfSumMinus_47_52 = 0.0f;

    hfSumPlus_48_49 = hfSumMinus_48_49 = 0.0f;
    hfSumPlus_48_50 = hfSumMinus_48_50 = 0.0f;
    hfSumPlus_48_51 = hfSumMinus_48_51 = 0.0f;
    hfSumPlus_48_52 = hfSumMinus_48_52 = 0.0f;

    hfSumPlus_49_50 = hfSumMinus_49_50 = 0.0f;
    hfSumPlus_49_51 = hfSumMinus_49_51 = 0.0f;
    hfSumPlus_49_52 = hfSumMinus_49_52 = 0.0f;

    hfSumPlus_50_51 = hfSumMinus_50_51 = 0.0f;
    hfSumPlus_50_52 = hfSumMinus_50_52 = 0.0f;

    hfSumPlus_51_52 = hfSumMinus_51_52 = 0.0f;
    // reset --end
    

    Int_t nTow = branchTower->GetEntriesFast();
    for (Int_t t = 0; t < nTow; ++t) {
      auto* tow = (Tower*)branchTower->At(t);
      float η = tow->Eta, E = tow->E;
    
      // HF+ region (3.0 ≤ η ≤ 5.2)
      if (η >= 3.0f && η <= 5.2f && E >= 0.0f) {
        hfSumPlus          += E;
        maxHFEnergyPlus     = std::max(maxHFEnergyPlus, E);
      }
      // HF– region (–5.2 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -5.2f && E >= 0.0f) {
        hfSumMinus         += E;
        maxHFEnergyMinus    = std::max(maxHFEnergyMinus, E);
      }

      // RANGOS: +eta
      // HF+ region (3.0 ≤ η ≤ 3.1)
      if (η >= 3.0f && η <= 3.1f && E >= Emin) {
        hfSumPlus_30_31          += E;
        maxHFEnergyPlus_30_31     = std::max(maxHFEnergyPlus_30_31, E);
      }
      // HF– region (-3.1 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.1f && E >= Emin) {
        hfSumMinus_30_31         += E;
        maxHFEnergyMinus_30_31    = std::max(maxHFEnergyMinus_30_31, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.2)
      if (η >= 3.0f && η <= 3.2f && E >= Emin) {
        hfSumPlus_30_32          += E;
        maxHFEnergyPlus_30_32     = std::max(maxHFEnergyPlus_30_32, E);
      }
      // HF– region (-3.2 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.2f && E >= Emin) {
        hfSumMinus_30_32         += E;
        maxHFEnergyMinus_30_32    = std::max(maxHFEnergyMinus_30_32, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.3)
      if (η >= 3.0f && η <= 3.3f && E >= Emin) {
        hfSumPlus_30_33          += E;
        maxHFEnergyPlus_30_33     = std::max(maxHFEnergyPlus_30_33, E);
      }
      // HF– region (-3.3 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.3f && E >= Emin) {
        hfSumMinus_30_33         += E;
        maxHFEnergyMinus_30_33    = std::max(maxHFEnergyMinus_30_33, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.4)
      if (η >= 3.0f && η <= 3.4f && E >= Emin) {
        hfSumPlus_30_34          += E;
        maxHFEnergyPlus_30_34     = std::max(maxHFEnergyPlus_30_34, E);
      }
      // HF– region (-3.4 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.4f && E >= Emin) {
        hfSumMinus_30_34         += E;
        maxHFEnergyMinus_30_34    = std::max(maxHFEnergyMinus_30_34, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.5)
      if (η >= 3.0f && η <= 3.5f && E >= Emin) {
        hfSumPlus_30_35          += E;
        maxHFEnergyPlus_30_35     = std::max(maxHFEnergyPlus_30_35, E);
      }
      // HF– region (-3.5 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.5f && E >= Emin) {
        hfSumMinus_30_35         += E;
        maxHFEnergyMinus_30_35    = std::max(maxHFEnergyMinus_30_35, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.6)
      if (η >= 3.0f && η <= 3.6f && E >= Emin) {
        hfSumPlus_30_36          += E;
        maxHFEnergyPlus_30_36     = std::max(maxHFEnergyPlus_30_36, E);
      }
      // HF– region (-3.6 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.6f && E >= Emin) {
        hfSumMinus_30_36         += E;
        maxHFEnergyMinus_30_36    = std::max(maxHFEnergyMinus_30_36, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.7)
      if (η >= 3.0f && η <= 3.7f && E >= Emin) {
        hfSumPlus_30_37          += E;
        maxHFEnergyPlus_30_37     = std::max(maxHFEnergyPlus_30_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.7f && E >= Emin) {
        hfSumMinus_30_37         += E;
        maxHFEnergyMinus_30_37    = std::max(maxHFEnergyMinus_30_37, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.8)
      if (η >= 3.0f && η <= 3.8f && E >= Emin) {
        hfSumPlus_30_38          += E;
        maxHFEnergyPlus_30_38     = std::max(maxHFEnergyPlus_30_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.8f && E >= Emin) {
        hfSumMinus_30_38         += E;
        maxHFEnergyMinus_30_38    = std::max(maxHFEnergyMinus_30_38, E);
      }

      // HF+ region (3.0 ≤ η ≤ 3.9)
      if (η >= 3.0f && η <= 3.9f && E >= Emin) {
        hfSumPlus_30_39          += E;
        maxHFEnergyPlus_30_39     = std::max(maxHFEnergyPlus_30_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -3.9f && E >= Emin) {
        hfSumMinus_30_39         += E;
        maxHFEnergyMinus_30_39    = std::max(maxHFEnergyMinus_30_39, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.0)
      if (η >= 3.0f && η <= 4.0f && E >= Emin) {
        hfSumPlus_30_40          += E;
        maxHFEnergyPlus_30_40     = std::max(maxHFEnergyPlus_30_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.0f && E >= Emin) {
        hfSumMinus_30_40         += E;
        maxHFEnergyMinus_30_40    = std::max(maxHFEnergyMinus_30_40, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.1)
      if (η >= 3.0f && η <= 4.1f && E >= Emin) {
        hfSumPlus_30_41          += E;
        maxHFEnergyPlus_30_41     = std::max(maxHFEnergyPlus_30_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.1f && E >= Emin) {
        hfSumMinus_30_41         += E;
        maxHFEnergyMinus_30_41    = std::max(maxHFEnergyMinus_30_41, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.2)
      if (η >= 3.0f && η <= 4.2f && E >= Emin) {
        hfSumPlus_30_42          += E;
        maxHFEnergyPlus_30_42     = std::max(maxHFEnergyPlus_30_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.2f && E >= Emin) {
        hfSumMinus_30_42         += E;
        maxHFEnergyMinus_30_42    = std::max(maxHFEnergyMinus_30_42, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.3)
      if (η >= 3.0f && η <= 4.3f && E >= Emin) {
        hfSumPlus_30_43          += E;
        maxHFEnergyPlus_30_43     = std::max(maxHFEnergyPlus_30_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.3f && E >= Emin) {
        hfSumMinus_30_43         += E;
        maxHFEnergyMinus_30_43    = std::max(maxHFEnergyMinus_30_43, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.4)
      if (η >= 3.0f && η <= 4.4f && E >= Emin) {
        hfSumPlus_30_44          += E;
        maxHFEnergyPlus_30_44     = std::max(maxHFEnergyPlus_30_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.4f && E >= Emin) {
        hfSumMinus_30_44         += E;
        maxHFEnergyMinus_30_44    = std::max(maxHFEnergyMinus_30_44, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.5)
      if (η >= 3.0f && η <= 4.5f && E >= Emin) {
        hfSumPlus_30_45          += E;
        maxHFEnergyPlus_30_45     = std::max(maxHFEnergyPlus_30_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.5f && E >= Emin) {
        hfSumMinus_30_45         += E;
        maxHFEnergyMinus_30_45    = std::max(maxHFEnergyMinus_30_45, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.6)
      if (η >= 3.0f && η <= 4.6f && E >= Emin) {
        hfSumPlus_30_46          += E;
        maxHFEnergyPlus_30_46     = std::max(maxHFEnergyPlus_30_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.6f && E >= Emin) {
        hfSumMinus_30_46         += E;
        maxHFEnergyMinus_30_46    = std::max(maxHFEnergyMinus_30_46, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.7)
      if (η >= 3.0f && η <= 4.7f && E >= Emin) {
        hfSumPlus_30_47          += E;
        maxHFEnergyPlus_30_47     = std::max(maxHFEnergyPlus_30_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.7f && E >= Emin) {
        hfSumMinus_30_47         += E;
        maxHFEnergyMinus_30_47    = std::max(maxHFEnergyMinus_30_47, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.8)
      if (η >= 3.0f && η <= 4.8f && E >= Emin) {
        hfSumPlus_30_48          += E;
        maxHFEnergyPlus_30_48     = std::max(maxHFEnergyPlus_30_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.8f && E >= Emin) {
        hfSumMinus_30_48         += E;
        maxHFEnergyMinus_30_48    = std::max(maxHFEnergyMinus_30_48, E);
      }

      // HF+ region (3.0 ≤ η ≤ 4.9)
      if (η >= 3.0f && η <= 4.9f && E >= Emin) {
        hfSumPlus_30_49          += E;
        maxHFEnergyPlus_30_49     = std::max(maxHFEnergyPlus_30_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -4.9f && E >= Emin) {
        hfSumMinus_30_49         += E;
        maxHFEnergyMinus_30_49    = std::max(maxHFEnergyMinus_30_49, E);
      }

      // HF+ region (3.0 ≤ η ≤ 5.0)
      if (η >= 3.0f && η <= 5.0f && E >= Emin) {
        hfSumPlus_30_50          += E;
        maxHFEnergyPlus_30_50     = std::max(maxHFEnergyPlus_30_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -5.0f && E >= Emin) {
        hfSumMinus_30_50         += E;
        maxHFEnergyMinus_30_50    = std::max(maxHFEnergyMinus_30_50, E);
      }

      // HF+ region (3.0 ≤ η ≤ 5.1)
      if (η >= 3.0f && η <= 5.1f && E >= Emin) {
        hfSumPlus_30_51          += E;
        maxHFEnergyPlus_30_51     = std::max(maxHFEnergyPlus_30_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -5.1f && E >= Emin) {
        hfSumMinus_30_51         += E;
        maxHFEnergyMinus_30_51    = std::max(maxHFEnergyMinus_30_51, E);
      }

      // HF+ region (3.0 ≤ η ≤ 5.2)
      if (η >= 3.0f && η <= 5.2f && E >= Emin) {
        hfSumPlus_30_52          += E;
        maxHFEnergyPlus_30_52     = std::max(maxHFEnergyPlus_30_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.0)
      if (η <= -3.0f && η >= -5.2f && E >= Emin) {
        hfSumMinus_30_52         += E;
        maxHFEnergyMinus_30_52    = std::max(maxHFEnergyMinus_30_52, E);
      }



      // HF+ region (3.1 ≤ η ≤ 3.2)
      if (η >= 3.1f && η <= 3.2f && E >= Emin) {
        hfSumPlus_31_32          += E;
        maxHFEnergyPlus_31_32     = std::max(maxHFEnergyPlus_31_32, E);
      }
      // HF– region (-3.2 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.2f && E >= Emin) {
        hfSumMinus_31_32         += E;
        maxHFEnergyMinus_31_32    = std::max(maxHFEnergyMinus_31_32, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.3)
      if (η >= 3.1f && η <= 3.3f && E >= Emin) {
        hfSumPlus_31_33          += E;
        maxHFEnergyPlus_31_33     = std::max(maxHFEnergyPlus_31_33, E);
      }
      // HF– region (-3.3 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.3f && E >= Emin) {
        hfSumMinus_31_33         += E;
        maxHFEnergyMinus_31_33    = std::max(maxHFEnergyMinus_31_33, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.4)
      if (η >= 3.1f && η <= 3.4f && E >= Emin) {
        hfSumPlus_31_34          += E;
        maxHFEnergyPlus_31_34     = std::max(maxHFEnergyPlus_31_34, E);
      }
      // HF– region (-3.4 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.4f && E >= Emin) {
        hfSumMinus_31_34         += E;
        maxHFEnergyMinus_31_34    = std::max(maxHFEnergyMinus_31_34, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.5)
      if (η >= 3.1f && η <= 3.5f && E >= Emin) {
        hfSumPlus_31_35          += E;
        maxHFEnergyPlus_31_35     = std::max(maxHFEnergyPlus_31_35, E);
      }
      // HF– region (-3.5 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.5f && E >= Emin) {
        hfSumMinus_31_35         += E;
        maxHFEnergyMinus_31_35    = std::max(maxHFEnergyMinus_31_35, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.6)
      if (η >= 3.1f && η <= 3.6f && E >= Emin) {
        hfSumPlus_31_36          += E;
        maxHFEnergyPlus_31_36     = std::max(maxHFEnergyPlus_31_36, E);
      }
      // HF– region (-3.6 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.6f && E >= Emin) {
        hfSumMinus_31_36         += E;
        maxHFEnergyMinus_31_36    = std::max(maxHFEnergyMinus_31_36, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.7)
      if (η >= 3.1f && η <= 3.7f && E >= Emin) {
        hfSumPlus_31_37          += E;
        maxHFEnergyPlus_31_37     = std::max(maxHFEnergyPlus_31_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.7f && E >= Emin) {
        hfSumMinus_31_37         += E;
        maxHFEnergyMinus_31_37    = std::max(maxHFEnergyMinus_31_37, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.8)
      if (η >= 3.1f && η <= 3.8f && E >= Emin) {
        hfSumPlus_31_38          += E;
        maxHFEnergyPlus_31_38     = std::max(maxHFEnergyPlus_31_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.8f && E >= Emin) {
        hfSumMinus_31_38         += E;
        maxHFEnergyMinus_31_38    = std::max(maxHFEnergyMinus_31_38, E);
      }

      // HF+ region (3.1 ≤ η ≤ 3.9)
      if (η >= 3.1f && η <= 3.9f && E >= Emin) {
        hfSumPlus_31_39          += E;
        maxHFEnergyPlus_31_39     = std::max(maxHFEnergyPlus_31_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -3.9f && E >= Emin) {
        hfSumMinus_31_39         += E;
        maxHFEnergyMinus_31_39    = std::max(maxHFEnergyMinus_31_39, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.0)
      if (η >= 3.1f && η <= 4.0f && E >= Emin) {
        hfSumPlus_31_40          += E;
        maxHFEnergyPlus_31_40     = std::max(maxHFEnergyPlus_31_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.0f && E >= Emin) {
        hfSumMinus_31_40         += E;
        maxHFEnergyMinus_31_40    = std::max(maxHFEnergyMinus_31_40, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.1)
      if (η >= 3.1f && η <= 4.1f && E >= Emin) {
        hfSumPlus_31_41          += E;
        maxHFEnergyPlus_31_41     = std::max(maxHFEnergyPlus_31_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.1f && E >= Emin) {
        hfSumMinus_31_41         += E;
        maxHFEnergyMinus_31_41    = std::max(maxHFEnergyMinus_31_41, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.2)
      if (η >= 3.1f && η <= 4.2f && E >= Emin) {
        hfSumPlus_31_42          += E;
        maxHFEnergyPlus_31_42     = std::max(maxHFEnergyPlus_31_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.2f && E >= Emin) {
        hfSumMinus_31_42         += E;
        maxHFEnergyMinus_31_42    = std::max(maxHFEnergyMinus_31_42, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.3)
      if (η >= 3.1f && η <= 4.3f && E >= Emin) {
        hfSumPlus_31_43          += E;
        maxHFEnergyPlus_31_43     = std::max(maxHFEnergyPlus_31_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.3f && E >= Emin) {
        hfSumMinus_31_43         += E;
        maxHFEnergyMinus_31_43    = std::max(maxHFEnergyMinus_31_43, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.4)
      if (η >= 3.1f && η <= 4.4f && E >= Emin) {
        hfSumPlus_31_44          += E;
        maxHFEnergyPlus_31_44     = std::max(maxHFEnergyPlus_31_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.4f && E >= Emin) {
        hfSumMinus_31_44         += E;
        maxHFEnergyMinus_31_44    = std::max(maxHFEnergyMinus_31_44, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.5)
      if (η >= 3.1f && η <= 4.5f && E >= Emin) {
        hfSumPlus_31_45          += E;
        maxHFEnergyPlus_31_45     = std::max(maxHFEnergyPlus_31_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.5f && E >= Emin) {
        hfSumMinus_31_45         += E;
        maxHFEnergyMinus_31_45    = std::max(maxHFEnergyMinus_31_45, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.6)
      if (η >= 3.1f && η <= 4.6f && E >= Emin) {
        hfSumPlus_31_46          += E;
        maxHFEnergyPlus_31_46     = std::max(maxHFEnergyPlus_31_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.6f && E >= Emin) {
        hfSumMinus_31_46         += E;
        maxHFEnergyMinus_31_46    = std::max(maxHFEnergyMinus_31_46, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.7)
      if (η >= 3.1f && η <= 4.7f && E >= Emin) {
        hfSumPlus_31_47          += E;
        maxHFEnergyPlus_31_47     = std::max(maxHFEnergyPlus_31_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.7f && E >= Emin) {
        hfSumMinus_31_47         += E;
        maxHFEnergyMinus_31_47    = std::max(maxHFEnergyMinus_31_47, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.8)
      if (η >= 3.1f && η <= 4.8f && E >= Emin) {
        hfSumPlus_31_48          += E;
        maxHFEnergyPlus_31_48     = std::max(maxHFEnergyPlus_31_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.8f && E >= Emin) {
        hfSumMinus_31_48         += E;
        maxHFEnergyMinus_31_48    = std::max(maxHFEnergyMinus_31_48, E);
      }

      // HF+ region (3.1 ≤ η ≤ 4.9)
      if (η >= 3.1f && η <= 4.9f && E >= Emin) {
        hfSumPlus_31_49          += E;
        maxHFEnergyPlus_31_49     = std::max(maxHFEnergyPlus_31_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -4.9f && E >= Emin) {
        hfSumMinus_31_49         += E;
        maxHFEnergyMinus_31_49    = std::max(maxHFEnergyMinus_31_49, E);
      }

      // HF+ region (3.1 ≤ η ≤ 5.0)
      if (η >= 3.1f && η <= 5.0f && E >= Emin) {
        hfSumPlus_31_50          += E;
        maxHFEnergyPlus_31_50     = std::max(maxHFEnergyPlus_31_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -5.0f && E >= Emin) {
        hfSumMinus_31_50         += E;
        maxHFEnergyMinus_31_50    = std::max(maxHFEnergyMinus_31_50, E);
      }

      // HF+ region (3.1 ≤ η ≤ 5.1)
      if (η >= 3.1f && η <= 5.1f && E >= Emin) {
        hfSumPlus_31_51          += E;
        maxHFEnergyPlus_31_51     = std::max(maxHFEnergyPlus_31_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -5.1f && E >= Emin) {
        hfSumMinus_31_51         += E;
        maxHFEnergyMinus_31_51    = std::max(maxHFEnergyMinus_31_51, E);
      }

      // HF+ region (3.1 ≤ η ≤ 5.2)
      if (η >= 3.1f && η <= 5.2f && E >= Emin) {
        hfSumPlus_31_52          += E;
        maxHFEnergyPlus_31_52     = std::max(maxHFEnergyPlus_31_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.1)
      if (η <= -3.1f && η >= -5.2f && E >= Emin) {
        hfSumMinus_31_52         += E;
        maxHFEnergyMinus_31_52    = std::max(maxHFEnergyMinus_31_52, E);
      }



      // HF+ region (3.2 ≤ η ≤ 3.3)
      if (η >= 3.2f && η <= 3.3f && E >= Emin) {
        hfSumPlus_32_33          += E;
        maxHFEnergyPlus_32_33     = std::max(maxHFEnergyPlus_32_33, E);
      }
      // HF– region (-3.3 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.3f && E >= Emin) {
        hfSumMinus_32_33         += E;
        maxHFEnergyMinus_32_33    = std::max(maxHFEnergyMinus_32_33, E);
      }

      // HF+ region (3.2 ≤ η ≤ 3.4)
      if (η >= 3.2f && η <= 3.4f && E >= Emin) {
        hfSumPlus_32_34          += E;
        maxHFEnergyPlus_32_34     = std::max(maxHFEnergyPlus_32_34, E);
      }
      // HF– region (-3.4 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.4f && E >= Emin) {
        hfSumMinus_32_34         += E;
        maxHFEnergyMinus_32_34    = std::max(maxHFEnergyMinus_32_34, E);
      }

      // HF+ region (3.2 ≤ η ≤ 3.5)
      if (η >= 3.2f && η <= 3.5f && E >= Emin) {
        hfSumPlus_32_35          += E;
        maxHFEnergyPlus_32_35     = std::max(maxHFEnergyPlus_32_35, E);
      }
      // HF– region (-3.5 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.5f && E >= Emin) {
        hfSumMinus_32_35         += E;
        maxHFEnergyMinus_32_35    = std::max(maxHFEnergyMinus_32_35, E);
      }

      // HF+ region (3.2 ≤ η ≤ 3.6)
      if (η >= 3.2f && η <= 3.6f && E >= Emin) {
        hfSumPlus_32_36          += E;
        maxHFEnergyPlus_32_36     = std::max(maxHFEnergyPlus_32_36, E);
      }
      // HF– region (-3.6 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.6f && E >= Emin) {
        hfSumMinus_32_36         += E;
        maxHFEnergyMinus_32_36    = std::max(maxHFEnergyMinus_32_36, E);
      }

      // HF+ region (3.2 ≤ η ≤ 3.7)
      if (η >= 3.2f && η <= 3.7f && E >= Emin) {
        hfSumPlus_32_37          += E;
        maxHFEnergyPlus_32_37     = std::max(maxHFEnergyPlus_32_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.7f && E >= Emin) {
        hfSumMinus_32_37         += E;
        maxHFEnergyMinus_32_37    = std::max(maxHFEnergyMinus_32_37, E);
      }

      // HF+ region (3.2 ≤ η ≤ 3.8)
      if (η >= 3.2f && η <= 3.8f && E >= Emin) {
        hfSumPlus_32_38          += E;
        maxHFEnergyPlus_32_38     = std::max(maxHFEnergyPlus_32_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.8f && E >= Emin) {
        hfSumMinus_32_38         += E;
        maxHFEnergyMinus_32_38    = std::max(maxHFEnergyMinus_32_38, E);
      }

      // HF+ region (3.2 ≤ η ≤ 3.9)
      if (η >= 3.2f && η <= 3.9f && E >= Emin) {
        hfSumPlus_32_39          += E;
        maxHFEnergyPlus_32_39     = std::max(maxHFEnergyPlus_32_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -3.9f && E >= Emin) {
        hfSumMinus_32_39         += E;
        maxHFEnergyMinus_32_39    = std::max(maxHFEnergyMinus_32_39, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.0)
      if (η >= 3.2f && η <= 4.0f && E >= Emin) {
        hfSumPlus_32_40          += E;
        maxHFEnergyPlus_32_40     = std::max(maxHFEnergyPlus_32_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.0f && E >= Emin) {
        hfSumMinus_32_40         += E;
        maxHFEnergyMinus_32_40    = std::max(maxHFEnergyMinus_32_40, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.1)
      if (η >= 3.2f && η <= 4.1f && E >= Emin) {
        hfSumPlus_32_41          += E;
        maxHFEnergyPlus_32_41     = std::max(maxHFEnergyPlus_32_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.1f && E >= Emin) {
        hfSumMinus_32_41         += E;
        maxHFEnergyMinus_32_41    = std::max(maxHFEnergyMinus_32_41, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.2)
      if (η >= 3.2f && η <= 4.2f && E >= Emin) {
        hfSumPlus_32_42          += E;
        maxHFEnergyPlus_32_42     = std::max(maxHFEnergyPlus_32_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.2f && E >= Emin) {
        hfSumMinus_32_42         += E;
        maxHFEnergyMinus_32_42    = std::max(maxHFEnergyMinus_32_42, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.3)
      if (η >= 3.2f && η <= 4.3f && E >= Emin) {
        hfSumPlus_32_43          += E;
        maxHFEnergyPlus_32_43     = std::max(maxHFEnergyPlus_32_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.3f && E >= Emin) {
        hfSumMinus_32_43         += E;
        maxHFEnergyMinus_32_43    = std::max(maxHFEnergyMinus_32_43, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.4)
      if (η >= 3.2f && η <= 4.4f && E >= Emin) {
        hfSumPlus_32_44          += E;
        maxHFEnergyPlus_32_44     = std::max(maxHFEnergyPlus_32_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.4f && E >= Emin) {
        hfSumMinus_32_44         += E;
        maxHFEnergyMinus_32_44    = std::max(maxHFEnergyMinus_32_44, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.5)
      if (η >= 3.2f && η <= 4.5f && E >= Emin) {
        hfSumPlus_32_45          += E;
        maxHFEnergyPlus_32_45     = std::max(maxHFEnergyPlus_32_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.5f && E >= Emin) {
        hfSumMinus_32_45         += E;
        maxHFEnergyMinus_32_45    = std::max(maxHFEnergyMinus_32_45, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.6)
      if (η >= 3.2f && η <= 4.6f && E >= Emin) {
        hfSumPlus_32_46          += E;
        maxHFEnergyPlus_32_46     = std::max(maxHFEnergyPlus_32_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.6f && E >= Emin) {
        hfSumMinus_32_46         += E;
        maxHFEnergyMinus_32_46    = std::max(maxHFEnergyMinus_32_46, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.7)
      if (η >= 3.2f && η <= 4.7f && E >= Emin) {
        hfSumPlus_32_47          += E;
        maxHFEnergyPlus_32_47     = std::max(maxHFEnergyPlus_32_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.7f && E >= Emin) {
        hfSumMinus_32_47         += E;
        maxHFEnergyMinus_32_47    = std::max(maxHFEnergyMinus_32_47, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.8)
      if (η >= 3.2f && η <= 4.8f && E >= Emin) {
        hfSumPlus_32_48          += E;
        maxHFEnergyPlus_32_48     = std::max(maxHFEnergyPlus_32_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.8f && E >= Emin) {
        hfSumMinus_32_48         += E;
        maxHFEnergyMinus_32_48    = std::max(maxHFEnergyMinus_32_48, E);
      }

      // HF+ region (3.2 ≤ η ≤ 4.9)
      if (η >= 3.2f && η <= 4.9f && E >= Emin) {
        hfSumPlus_32_49          += E;
        maxHFEnergyPlus_32_49     = std::max(maxHFEnergyPlus_32_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -4.9f && E >= Emin) {
        hfSumMinus_32_49         += E;
        maxHFEnergyMinus_32_49    = std::max(maxHFEnergyMinus_32_49, E);
      }

      // HF+ region (3.2 ≤ η ≤ 5.0)
      if (η >= 3.2f && η <= 5.0f && E >= Emin) {
        hfSumPlus_32_50          += E;
        maxHFEnergyPlus_32_50     = std::max(maxHFEnergyPlus_32_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -5.0f && E >= Emin) {
        hfSumMinus_32_50         += E;
        maxHFEnergyMinus_32_50    = std::max(maxHFEnergyMinus_32_50, E);
      }

      // HF+ region (3.2 ≤ η ≤ 5.1)
      if (η >= 3.2f && η <= 5.1f && E >= Emin) {
        hfSumPlus_32_51          += E;
        maxHFEnergyPlus_32_51     = std::max(maxHFEnergyPlus_32_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -5.1f && E >= Emin) {
        hfSumMinus_32_51         += E;
        maxHFEnergyMinus_32_51    = std::max(maxHFEnergyMinus_32_51, E);
      }

      // HF+ region (3.2 ≤ η ≤ 5.2)
      if (η >= 3.2f && η <= 5.2f && E >= Emin) {
        hfSumPlus_32_52          += E;
        maxHFEnergyPlus_32_52     = std::max(maxHFEnergyPlus_32_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.2)
      if (η <= -3.2f && η >= -5.2f && E >= Emin) {
        hfSumMinus_32_52         += E;
        maxHFEnergyMinus_32_52    = std::max(maxHFEnergyMinus_32_52, E);
      }



      // HF+ region (3.3 ≤ η ≤ 3.4)
      if (η >= 3.3f && η <= 3.4f && E >= Emin) {
        hfSumPlus_33_34          += E;
        maxHFEnergyPlus_33_34     = std::max(maxHFEnergyPlus_33_34, E);
      }
      // HF– region (-3.4 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -3.4f && E >= Emin) {
        hfSumMinus_33_34         += E;
        maxHFEnergyMinus_33_34    = std::max(maxHFEnergyMinus_33_34, E);
      }

      // HF+ region (3.3 ≤ η ≤ 3.5)
      if (η >= 3.3f && η <= 3.5f && E >= Emin) {
        hfSumPlus_33_35          += E;
        maxHFEnergyPlus_33_35     = std::max(maxHFEnergyPlus_33_35, E);
      }
      // HF– region (-3.5 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -3.5f && E >= Emin) {
        hfSumMinus_33_35         += E;
        maxHFEnergyMinus_33_35    = std::max(maxHFEnergyMinus_33_35, E);
      }

      // HF+ region (3.3 ≤ η ≤ 3.6)
      if (η >= 3.3f && η <= 3.6f && E >= Emin) {
        hfSumPlus_33_36          += E;
        maxHFEnergyPlus_33_36     = std::max(maxHFEnergyPlus_33_36, E);
      }
      // HF– region (-3.6 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -3.6f && E >= Emin) {
        hfSumMinus_33_36         += E;
        maxHFEnergyMinus_33_36    = std::max(maxHFEnergyMinus_33_36, E);
      }

      // HF+ region (3.3 ≤ η ≤ 3.7)
      if (η >= 3.3f && η <= 3.7f && E >= Emin) {
        hfSumPlus_33_37          += E;
        maxHFEnergyPlus_33_37     = std::max(maxHFEnergyPlus_33_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -3.7f && E >= Emin) {
        hfSumMinus_33_37         += E;
        maxHFEnergyMinus_33_37    = std::max(maxHFEnergyMinus_33_37, E);
      }

      // HF+ region (3.3 ≤ η ≤ 3.8)
      if (η >= 3.3f && η <= 3.8f && E >= Emin) {
        hfSumPlus_33_38          += E;
        maxHFEnergyPlus_33_38     = std::max(maxHFEnergyPlus_33_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -3.8f && E >= Emin) {
        hfSumMinus_33_38         += E;
        maxHFEnergyMinus_33_38    = std::max(maxHFEnergyMinus_33_38, E);
      }

      // HF+ region (3.3 ≤ η ≤ 3.9)
      if (η >= 3.3f && η <= 3.9f && E >= Emin) {
        hfSumPlus_33_39          += E;
        maxHFEnergyPlus_33_39     = std::max(maxHFEnergyPlus_33_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -3.9f && E >= Emin) {
        hfSumMinus_33_39         += E;
        maxHFEnergyMinus_33_39    = std::max(maxHFEnergyMinus_33_39, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.0)
      if (η >= 3.3f && η <= 4.0f && E >= Emin) {
        hfSumPlus_33_40          += E;
        maxHFEnergyPlus_33_40     = std::max(maxHFEnergyPlus_33_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.0f && E >= Emin) {
        hfSumMinus_33_40         += E;
        maxHFEnergyMinus_33_40    = std::max(maxHFEnergyMinus_33_40, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.1)
      if (η >= 3.3f && η <= 4.1f && E >= Emin) {
        hfSumPlus_33_41          += E;
        maxHFEnergyPlus_33_41     = std::max(maxHFEnergyPlus_33_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.1f && E >= Emin) {
        hfSumMinus_33_41         += E;
        maxHFEnergyMinus_33_41    = std::max(maxHFEnergyMinus_33_41, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.2)
      if (η >= 3.3f && η <= 4.2f && E >= Emin) {
        hfSumPlus_33_42          += E;
        maxHFEnergyPlus_33_42     = std::max(maxHFEnergyPlus_33_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.2f && E >= Emin) {
        hfSumMinus_33_42         += E;
        maxHFEnergyMinus_33_42    = std::max(maxHFEnergyMinus_33_42, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.3)
      if (η >= 3.3f && η <= 4.3f && E >= Emin) {
        hfSumPlus_33_43          += E;
        maxHFEnergyPlus_33_43     = std::max(maxHFEnergyPlus_33_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.3f && E >= Emin) {
        hfSumMinus_33_43         += E;
        maxHFEnergyMinus_33_43    = std::max(maxHFEnergyMinus_33_43, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.4)
      if (η >= 3.3f && η <= 4.4f && E >= Emin) {
        hfSumPlus_33_44          += E;
        maxHFEnergyPlus_33_44     = std::max(maxHFEnergyPlus_33_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.4f && E >= Emin) {
        hfSumMinus_33_44         += E;
        maxHFEnergyMinus_33_44    = std::max(maxHFEnergyMinus_33_44, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.5)
      if (η >= 3.3f && η <= 4.5f && E >= Emin) {
        hfSumPlus_33_45          += E;
        maxHFEnergyPlus_33_45     = std::max(maxHFEnergyPlus_33_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.5f && E >= Emin) {
        hfSumMinus_33_45         += E;
        maxHFEnergyMinus_33_45    = std::max(maxHFEnergyMinus_33_45, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.6)
      if (η >= 3.3f && η <= 4.6f && E >= Emin) {
        hfSumPlus_33_46          += E;
        maxHFEnergyPlus_33_46     = std::max(maxHFEnergyPlus_33_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.6f && E >= Emin) {
        hfSumMinus_33_46         += E;
        maxHFEnergyMinus_33_46    = std::max(maxHFEnergyMinus_33_46, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.7)
      if (η >= 3.3f && η <= 4.7f && E >= Emin) {
        hfSumPlus_33_47          += E;
        maxHFEnergyPlus_33_47     = std::max(maxHFEnergyPlus_33_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.7f && E >= Emin) {
        hfSumMinus_33_47         += E;
        maxHFEnergyMinus_33_47    = std::max(maxHFEnergyMinus_33_47, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.8)
      if (η >= 3.3f && η <= 4.8f && E >= Emin) {
        hfSumPlus_33_48          += E;
        maxHFEnergyPlus_33_48     = std::max(maxHFEnergyPlus_33_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.8f && E >= Emin) {
        hfSumMinus_33_48         += E;
        maxHFEnergyMinus_33_48    = std::max(maxHFEnergyMinus_33_48, E);
      }

      // HF+ region (3.3 ≤ η ≤ 4.9)
      if (η >= 3.3f && η <= 4.9f && E >= Emin) {
        hfSumPlus_33_49          += E;
        maxHFEnergyPlus_33_49     = std::max(maxHFEnergyPlus_33_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -4.9f && E >= Emin) {
        hfSumMinus_33_49         += E;
        maxHFEnergyMinus_33_49    = std::max(maxHFEnergyMinus_33_49, E);
      }

      // HF+ region (3.3 ≤ η ≤ 5.0)
      if (η >= 3.3f && η <= 5.0f && E >= Emin) {
        hfSumPlus_33_50          += E;
        maxHFEnergyPlus_33_50     = std::max(maxHFEnergyPlus_33_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -5.0f && E >= Emin) {
        hfSumMinus_33_50         += E;
        maxHFEnergyMinus_33_50    = std::max(maxHFEnergyMinus_33_50, E);
      }

      // HF+ region (3.3 ≤ η ≤ 5.1)
      if (η >= 3.3f && η <= 5.1f && E >= Emin) {
        hfSumPlus_33_51          += E;
        maxHFEnergyPlus_33_51     = std::max(maxHFEnergyPlus_33_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -5.1f && E >= Emin) {
        hfSumMinus_33_51         += E;
        maxHFEnergyMinus_33_51    = std::max(maxHFEnergyMinus_33_51, E);
      }

      // HF+ region (3.3 ≤ η ≤ 5.2)
      if (η >= 3.3f && η <= 5.2f && E >= Emin) {
        hfSumPlus_33_52          += E;
        maxHFEnergyPlus_33_52     = std::max(maxHFEnergyPlus_33_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.3)
      if (η <= -3.3f && η >= -5.2f && E >= Emin) {
        hfSumMinus_33_52         += E;
        maxHFEnergyMinus_33_52    = std::max(maxHFEnergyMinus_33_52, E);
      }



      // HF+ region (3.4 ≤ η ≤ 3.5)
      if (η >= 3.4f && η <= 3.5f && E >= Emin) {
        hfSumPlus_34_35          += E;
        maxHFEnergyPlus_34_35     = std::max(maxHFEnergyPlus_34_35, E);
      }
      // HF– region (-3.5 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -3.5f && E >= Emin) {
        hfSumMinus_34_35         += E;
        maxHFEnergyMinus_34_35    = std::max(maxHFEnergyMinus_34_35, E);
      }

      // HF+ region (3.4 ≤ η ≤ 3.6)
      if (η >= 3.4f && η <= 3.6f && E >= Emin) {
        hfSumPlus_34_36          += E;
        maxHFEnergyPlus_34_36     = std::max(maxHFEnergyPlus_34_36, E);
      }
      // HF– region (-3.6 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -3.6f && E >= Emin) {
        hfSumMinus_34_36         += E;
        maxHFEnergyMinus_34_36    = std::max(maxHFEnergyMinus_34_36, E);
      }

      // HF+ region (3.4 ≤ η ≤ 3.7)
      if (η >= 3.4f && η <= 3.7f && E >= Emin) {
        hfSumPlus_34_37          += E;
        maxHFEnergyPlus_34_37     = std::max(maxHFEnergyPlus_34_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -3.7f && E >= Emin) {
        hfSumMinus_34_37         += E;
        maxHFEnergyMinus_34_37    = std::max(maxHFEnergyMinus_34_37, E);
      }

      // HF+ region (3.4 ≤ η ≤ 3.8)
      if (η >= 3.4f && η <= 3.8f && E >= Emin) {
        hfSumPlus_34_38          += E;
        maxHFEnergyPlus_34_38     = std::max(maxHFEnergyPlus_34_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -3.8f && E >= Emin) {
        hfSumMinus_34_38         += E;
        maxHFEnergyMinus_34_38    = std::max(maxHFEnergyMinus_34_38, E);
      }

      // HF+ region (3.4 ≤ η ≤ 3.9)
      if (η >= 3.4f && η <= 3.9f && E >= Emin) {
        hfSumPlus_34_39          += E;
        maxHFEnergyPlus_34_39     = std::max(maxHFEnergyPlus_34_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -3.9f && E >= Emin) {
        hfSumMinus_34_39         += E;
        maxHFEnergyMinus_34_39    = std::max(maxHFEnergyMinus_34_39, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.0)
      if (η >= 3.4f && η <= 4.0f && E >= Emin) {
        hfSumPlus_34_40          += E;
        maxHFEnergyPlus_34_40     = std::max(maxHFEnergyPlus_34_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.0f && E >= Emin) {
        hfSumMinus_34_40         += E;
        maxHFEnergyMinus_34_40    = std::max(maxHFEnergyMinus_34_40, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.1)
      if (η >= 3.4f && η <= 4.1f && E >= Emin) {
        hfSumPlus_34_41          += E;
        maxHFEnergyPlus_34_41     = std::max(maxHFEnergyPlus_34_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.1f && E >= Emin) {
        hfSumMinus_34_41         += E;
        maxHFEnergyMinus_34_41    = std::max(maxHFEnergyMinus_34_41, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.2)
      if (η >= 3.4f && η <= 4.2f && E >= Emin) {
        hfSumPlus_34_42          += E;
        maxHFEnergyPlus_34_42     = std::max(maxHFEnergyPlus_34_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.2f && E >= Emin) {
        hfSumMinus_34_42         += E;
        maxHFEnergyMinus_34_42    = std::max(maxHFEnergyMinus_34_42, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.3)
      if (η >= 3.4f && η <= 4.3f && E >= Emin) {
        hfSumPlus_34_43          += E;
        maxHFEnergyPlus_34_43     = std::max(maxHFEnergyPlus_34_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.3f && E >= Emin) {
        hfSumMinus_34_43         += E;
        maxHFEnergyMinus_34_43    = std::max(maxHFEnergyMinus_34_43, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.4)
      if (η >= 3.4f && η <= 4.4f && E >= Emin) {
        hfSumPlus_34_44          += E;
        maxHFEnergyPlus_34_44     = std::max(maxHFEnergyPlus_34_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.4f && E >= Emin) {
        hfSumMinus_34_44         += E;
        maxHFEnergyMinus_34_44    = std::max(maxHFEnergyMinus_34_44, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.5)
      if (η >= 3.4f && η <= 4.5f && E >= Emin) {
        hfSumPlus_34_45          += E;
        maxHFEnergyPlus_34_45     = std::max(maxHFEnergyPlus_34_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.5f && E >= Emin) {
        hfSumMinus_34_45         += E;
        maxHFEnergyMinus_34_45    = std::max(maxHFEnergyMinus_34_45, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.6)
      if (η >= 3.4f && η <= 4.6f && E >= Emin) {
        hfSumPlus_34_46          += E;
        maxHFEnergyPlus_34_46     = std::max(maxHFEnergyPlus_34_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.6f && E >= Emin) {
        hfSumMinus_34_46         += E;
        maxHFEnergyMinus_34_46    = std::max(maxHFEnergyMinus_34_46, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.7)
      if (η >= 3.4f && η <= 4.7f && E >= Emin) {
        hfSumPlus_34_47          += E;
        maxHFEnergyPlus_34_47     = std::max(maxHFEnergyPlus_34_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.7f && E >= Emin) {
        hfSumMinus_34_47         += E;
        maxHFEnergyMinus_34_47    = std::max(maxHFEnergyMinus_34_47, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.8)
      if (η >= 3.4f && η <= 4.8f && E >= Emin) {
        hfSumPlus_34_48          += E;
        maxHFEnergyPlus_34_48     = std::max(maxHFEnergyPlus_34_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.8f && E >= Emin) {
        hfSumMinus_34_48         += E;
        maxHFEnergyMinus_34_48    = std::max(maxHFEnergyMinus_34_48, E);
      }

      // HF+ region (3.4 ≤ η ≤ 4.9)
      if (η >= 3.4f && η <= 4.9f && E >= Emin) {
        hfSumPlus_34_49          += E;
        maxHFEnergyPlus_34_49     = std::max(maxHFEnergyPlus_34_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -4.9f && E >= Emin) {
        hfSumMinus_34_49         += E;
        maxHFEnergyMinus_34_49    = std::max(maxHFEnergyMinus_34_49, E);
      }

      // HF+ region (3.4 ≤ η ≤ 5.0)
      if (η >= 3.4f && η <= 5.0f && E >= Emin) {
        hfSumPlus_34_50          += E;
        maxHFEnergyPlus_34_50     = std::max(maxHFEnergyPlus_34_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -5.0f && E >= Emin) {
        hfSumMinus_34_50         += E;
        maxHFEnergyMinus_34_50    = std::max(maxHFEnergyMinus_34_50, E);
      }

      // HF+ region (3.4 ≤ η ≤ 5.1)
      if (η >= 3.4f && η <= 5.1f && E >= Emin) {
        hfSumPlus_34_51          += E;
        maxHFEnergyPlus_34_51     = std::max(maxHFEnergyPlus_34_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -5.1f && E >= Emin) {
        hfSumMinus_34_51         += E;
        maxHFEnergyMinus_34_51    = std::max(maxHFEnergyMinus_34_51, E);
      }

      // HF+ region (3.4 ≤ η ≤ 5.2)
      if (η >= 3.4f && η <= 5.2f && E >= Emin) {
        hfSumPlus_34_52          += E;
        maxHFEnergyPlus_34_52     = std::max(maxHFEnergyPlus_34_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.4)
      if (η <= -3.4f && η >= -5.2f && E >= Emin) {
        hfSumMinus_34_52         += E;
        maxHFEnergyMinus_34_52    = std::max(maxHFEnergyMinus_34_52, E);
      }



      // HF+ region (3.5 ≤ η ≤ 3.6)
      if (η >= 3.5f && η <= 3.6f && E >= Emin) {
        hfSumPlus_35_36          += E;
        maxHFEnergyPlus_35_36     = std::max(maxHFEnergyPlus_35_36, E);
      }
      // HF– region (-3.6 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -3.6f && E >= Emin) {
        hfSumMinus_35_36         += E;
        maxHFEnergyMinus_35_36    = std::max(maxHFEnergyMinus_35_36, E);
      }

      // HF+ region (3.5 ≤ η ≤ 3.7)
      if (η >= 3.5f && η <= 3.7f && E >= Emin) {
        hfSumPlus_35_37          += E;
        maxHFEnergyPlus_35_37     = std::max(maxHFEnergyPlus_35_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -3.7f && E >= Emin) {
        hfSumMinus_35_37         += E;
        maxHFEnergyMinus_35_37    = std::max(maxHFEnergyMinus_35_37, E);
      }

      // HF+ region (3.5 ≤ η ≤ 3.8)
      if (η >= 3.5f && η <= 3.8f && E >= Emin) {
        hfSumPlus_35_38          += E;
        maxHFEnergyPlus_35_38     = std::max(maxHFEnergyPlus_35_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -3.8f && E >= Emin) {
        hfSumMinus_35_38         += E;
        maxHFEnergyMinus_35_38    = std::max(maxHFEnergyMinus_35_38, E);
      }

      // HF+ region (3.5 ≤ η ≤ 3.9)
      if (η >= 3.5f && η <= 3.9f && E >= Emin) {
        hfSumPlus_35_39          += E;
        maxHFEnergyPlus_35_39     = std::max(maxHFEnergyPlus_35_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -3.9f && E >= Emin) {
        hfSumMinus_35_39         += E;
        maxHFEnergyMinus_35_39    = std::max(maxHFEnergyMinus_35_39, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.0)
      if (η >= 3.5f && η <= 4.0f && E >= Emin) {
        hfSumPlus_35_40          += E;
        maxHFEnergyPlus_35_40     = std::max(maxHFEnergyPlus_35_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.0f && E >= Emin) {
        hfSumMinus_35_40         += E;
        maxHFEnergyMinus_35_40    = std::max(maxHFEnergyMinus_35_40, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.1)
      if (η >= 3.5f && η <= 4.1f && E >= Emin) {
        hfSumPlus_35_41          += E;
        maxHFEnergyPlus_35_41     = std::max(maxHFEnergyPlus_35_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.1f && E >= Emin) {
        hfSumMinus_35_41         += E;
        maxHFEnergyMinus_35_41    = std::max(maxHFEnergyMinus_35_41, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.2)
      if (η >= 3.5f && η <= 4.2f && E >= Emin) {
        hfSumPlus_35_42          += E;
        maxHFEnergyPlus_35_42     = std::max(maxHFEnergyPlus_35_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.2f && E >= Emin) {
        hfSumMinus_35_42         += E;
        maxHFEnergyMinus_35_42    = std::max(maxHFEnergyMinus_35_42, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.3)
      if (η >= 3.5f && η <= 4.3f && E >= Emin) {
        hfSumPlus_35_43          += E;
        maxHFEnergyPlus_35_43     = std::max(maxHFEnergyPlus_35_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.3f && E >= Emin) {
        hfSumMinus_35_43         += E;
        maxHFEnergyMinus_35_43    = std::max(maxHFEnergyMinus_35_43, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.4)
      if (η >= 3.5f && η <= 4.4f && E >= Emin) {
        hfSumPlus_35_44          += E;
        maxHFEnergyPlus_35_44     = std::max(maxHFEnergyPlus_35_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.4f && E >= Emin) {
        hfSumMinus_35_44         += E;
        maxHFEnergyMinus_35_44    = std::max(maxHFEnergyMinus_35_44, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.5)
      if (η >= 3.5f && η <= 4.5f && E >= Emin) {
        hfSumPlus_35_45          += E;
        maxHFEnergyPlus_35_45     = std::max(maxHFEnergyPlus_35_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.5f && E >= Emin) {
        hfSumMinus_35_45         += E;
        maxHFEnergyMinus_35_45    = std::max(maxHFEnergyMinus_35_45, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.6)
      if (η >= 3.5f && η <= 4.6f && E >= Emin) {
        hfSumPlus_35_46          += E;
        maxHFEnergyPlus_35_46     = std::max(maxHFEnergyPlus_35_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.6f && E >= Emin) {
        hfSumMinus_35_46         += E;
        maxHFEnergyMinus_35_46    = std::max(maxHFEnergyMinus_35_46, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.7)
      if (η >= 3.5f && η <= 4.7f && E >= Emin) {
        hfSumPlus_35_47          += E;
        maxHFEnergyPlus_35_47     = std::max(maxHFEnergyPlus_35_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.7f && E >= Emin) {
        hfSumMinus_35_47         += E;
        maxHFEnergyMinus_35_47    = std::max(maxHFEnergyMinus_35_47, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.8)
      if (η >= 3.5f && η <= 4.8f && E >= Emin) {
        hfSumPlus_35_48          += E;
        maxHFEnergyPlus_35_48     = std::max(maxHFEnergyPlus_35_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.8f && E >= Emin) {
        hfSumMinus_35_48         += E;
        maxHFEnergyMinus_35_48    = std::max(maxHFEnergyMinus_35_48, E);
      }

      // HF+ region (3.5 ≤ η ≤ 4.9)
      if (η >= 3.5f && η <= 4.9f && E >= Emin) {
        hfSumPlus_35_49          += E;
        maxHFEnergyPlus_35_49     = std::max(maxHFEnergyPlus_35_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -4.9f && E >= Emin) {
        hfSumMinus_35_49         += E;
        maxHFEnergyMinus_35_49    = std::max(maxHFEnergyMinus_35_49, E);
      }

      // HF+ region (3.5 ≤ η ≤ 5.0)
      if (η >= 3.5f && η <= 5.0f && E >= Emin) {
        hfSumPlus_35_50          += E;
        maxHFEnergyPlus_35_50     = std::max(maxHFEnergyPlus_35_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -5.0f && E >= Emin) {
        hfSumMinus_35_50         += E;
        maxHFEnergyMinus_35_50    = std::max(maxHFEnergyMinus_35_50, E);
      }

      // HF+ region (3.5 ≤ η ≤ 5.1)
      if (η >= 3.5f && η <= 5.1f && E >= Emin) {
        hfSumPlus_35_51          += E;
        maxHFEnergyPlus_35_51     = std::max(maxHFEnergyPlus_35_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -5.1f && E >= Emin) {
        hfSumMinus_35_51         += E;
        maxHFEnergyMinus_35_51    = std::max(maxHFEnergyMinus_35_51, E);
      }

      // HF+ region (3.5 ≤ η ≤ 5.2)
      if (η >= 3.5f && η <= 5.2f && E >= Emin) {
        hfSumPlus_35_52          += E;
        maxHFEnergyPlus_35_52     = std::max(maxHFEnergyPlus_35_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.5)
      if (η <= -3.5f && η >= -5.2f && E >= Emin) {
        hfSumMinus_35_52         += E;
        maxHFEnergyMinus_35_52    = std::max(maxHFEnergyMinus_35_52, E);
      }



      // HF+ region (3.6 ≤ η ≤ 3.7)
      if (η >= 3.6f && η <= 3.7f && E >= Emin) {
        hfSumPlus_36_37          += E;
        maxHFEnergyPlus_36_37     = std::max(maxHFEnergyPlus_36_37, E);
      }
      // HF– region (-3.7 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -3.7f && E >= Emin) {
        hfSumMinus_36_37         += E;
        maxHFEnergyMinus_36_37    = std::max(maxHFEnergyMinus_36_37, E);
      }

      // HF+ region (3.6 ≤ η ≤ 3.8)
      if (η >= 3.6f && η <= 3.8f && E >= Emin) {
        hfSumPlus_36_38          += E;
        maxHFEnergyPlus_36_38     = std::max(maxHFEnergyPlus_36_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -3.8f && E >= Emin) {
        hfSumMinus_36_38         += E;
        maxHFEnergyMinus_36_38    = std::max(maxHFEnergyMinus_36_38, E);
      }

      // HF+ region (3.6 ≤ η ≤ 3.9)
      if (η >= 3.6f && η <= 3.9f && E >= Emin) {
        hfSumPlus_36_39          += E;
        maxHFEnergyPlus_36_39     = std::max(maxHFEnergyPlus_36_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -3.9f && E >= Emin) {
        hfSumMinus_36_39         += E;
        maxHFEnergyMinus_36_39    = std::max(maxHFEnergyMinus_36_39, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.0)
      if (η >= 3.6f && η <= 4.0f && E >= Emin) {
        hfSumPlus_36_40          += E;
        maxHFEnergyPlus_36_40     = std::max(maxHFEnergyPlus_36_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.0f && E >= Emin) {
        hfSumMinus_36_40         += E;
        maxHFEnergyMinus_36_40    = std::max(maxHFEnergyMinus_36_40, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.1)
      if (η >= 3.6f && η <= 4.1f && E >= Emin) {
        hfSumPlus_36_41          += E;
        maxHFEnergyPlus_36_41     = std::max(maxHFEnergyPlus_36_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.1f && E >= Emin) {
        hfSumMinus_36_41         += E;
        maxHFEnergyMinus_36_41    = std::max(maxHFEnergyMinus_36_41, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.2)
      if (η >= 3.6f && η <= 4.2f && E >= Emin) {
        hfSumPlus_36_42          += E;
        maxHFEnergyPlus_36_42     = std::max(maxHFEnergyPlus_36_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.2f && E >= Emin) {
        hfSumMinus_36_42         += E;
        maxHFEnergyMinus_36_42    = std::max(maxHFEnergyMinus_36_42, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.3)
      if (η >= 3.6f && η <= 4.3f && E >= Emin) {
        hfSumPlus_36_43          += E;
        maxHFEnergyPlus_36_43     = std::max(maxHFEnergyPlus_36_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.3f && E >= Emin) {
        hfSumMinus_36_43         += E;
        maxHFEnergyMinus_36_43    = std::max(maxHFEnergyMinus_36_43, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.4)
      if (η >= 3.6f && η <= 4.4f && E >= Emin) {
        hfSumPlus_36_44          += E;
        maxHFEnergyPlus_36_44     = std::max(maxHFEnergyPlus_36_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.4f && E >= Emin) {
        hfSumMinus_36_44         += E;
        maxHFEnergyMinus_36_44    = std::max(maxHFEnergyMinus_36_44, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.5)
      if (η >= 3.6f && η <= 4.5f && E >= Emin) {
        hfSumPlus_36_45          += E;
        maxHFEnergyPlus_36_45     = std::max(maxHFEnergyPlus_36_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.5f && E >= Emin) {
        hfSumMinus_36_45         += E;
        maxHFEnergyMinus_36_45    = std::max(maxHFEnergyMinus_36_45, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.6)
      if (η >= 3.6f && η <= 4.6f && E >= Emin) {
        hfSumPlus_36_46          += E;
        maxHFEnergyPlus_36_46     = std::max(maxHFEnergyPlus_36_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.6f && E >= Emin) {
        hfSumMinus_36_46         += E;
        maxHFEnergyMinus_36_46    = std::max(maxHFEnergyMinus_36_46, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.7)
      if (η >= 3.6f && η <= 4.7f && E >= Emin) {
        hfSumPlus_36_47          += E;
        maxHFEnergyPlus_36_47     = std::max(maxHFEnergyPlus_36_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.7f && E >= Emin) {
        hfSumMinus_36_47         += E;
        maxHFEnergyMinus_36_47    = std::max(maxHFEnergyMinus_36_47, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.8)
      if (η >= 3.6f && η <= 4.8f && E >= Emin) {
        hfSumPlus_36_48          += E;
        maxHFEnergyPlus_36_48     = std::max(maxHFEnergyPlus_36_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.8f && E >= Emin) {
        hfSumMinus_36_48         += E;
        maxHFEnergyMinus_36_48    = std::max(maxHFEnergyMinus_36_48, E);
      }

      // HF+ region (3.6 ≤ η ≤ 4.9)
      if (η >= 3.6f && η <= 4.9f && E >= Emin) {
        hfSumPlus_36_49          += E;
        maxHFEnergyPlus_36_49     = std::max(maxHFEnergyPlus_36_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -4.9f && E >= Emin) {
        hfSumMinus_36_49         += E;
        maxHFEnergyMinus_36_49    = std::max(maxHFEnergyMinus_36_49, E);
      }

      // HF+ region (3.6 ≤ η ≤ 5.0)
      if (η >= 3.6f && η <= 5.0f && E >= Emin) {
        hfSumPlus_36_50          += E;
        maxHFEnergyPlus_36_50     = std::max(maxHFEnergyPlus_36_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -5.0f && E >= Emin) {
        hfSumMinus_36_50         += E;
        maxHFEnergyMinus_36_50    = std::max(maxHFEnergyMinus_36_50, E);
      }

      // HF+ region (3.6 ≤ η ≤ 5.1)
      if (η >= 3.6f && η <= 5.1f && E >= Emin) {
        hfSumPlus_36_51          += E;
        maxHFEnergyPlus_36_51     = std::max(maxHFEnergyPlus_36_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -5.1f && E >= Emin) {
        hfSumMinus_36_51         += E;
        maxHFEnergyMinus_36_51    = std::max(maxHFEnergyMinus_36_51, E);
      }

      // HF+ region (3.6 ≤ η ≤ 5.2)
      if (η >= 3.6f && η <= 5.2f && E >= Emin) {
        hfSumPlus_36_52          += E;
        maxHFEnergyPlus_36_52     = std::max(maxHFEnergyPlus_36_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.6)
      if (η <= -3.6f && η >= -5.2f && E >= Emin) {
        hfSumMinus_36_52         += E;
        maxHFEnergyMinus_36_52    = std::max(maxHFEnergyMinus_36_52, E);
      }



      // HF+ region (3.7 ≤ η ≤ 3.8)
      if (η >= 3.7f && η <= 3.8f && E >= Emin) {
        hfSumPlus_37_38          += E;
        maxHFEnergyPlus_37_38     = std::max(maxHFEnergyPlus_37_38, E);
      }
      // HF– region (-3.8 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -3.8f && E >= Emin) {
        hfSumMinus_37_38         += E;
        maxHFEnergyMinus_37_38    = std::max(maxHFEnergyMinus_37_38, E);
      }

      // HF+ region (3.7 ≤ η ≤ 3.9)
      if (η >= 3.7f && η <= 3.9f && E >= Emin) {
        hfSumPlus_37_39          += E;
        maxHFEnergyPlus_37_39     = std::max(maxHFEnergyPlus_37_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -3.9f && E >= Emin) {
        hfSumMinus_37_39         += E;
        maxHFEnergyMinus_37_39    = std::max(maxHFEnergyMinus_37_39, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.0)
      if (η >= 3.7f && η <= 4.0f && E >= Emin) {
        hfSumPlus_37_40          += E;
        maxHFEnergyPlus_37_40     = std::max(maxHFEnergyPlus_37_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.0f && E >= Emin) {
        hfSumMinus_37_40         += E;
        maxHFEnergyMinus_37_40    = std::max(maxHFEnergyMinus_37_40, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.1)
      if (η >= 3.7f && η <= 4.1f && E >= Emin) {
        hfSumPlus_37_41          += E;
        maxHFEnergyPlus_37_41     = std::max(maxHFEnergyPlus_37_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.1f && E >= Emin) {
        hfSumMinus_37_41         += E;
        maxHFEnergyMinus_37_41    = std::max(maxHFEnergyMinus_37_41, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.2)
      if (η >= 3.7f && η <= 4.2f && E >= Emin) {
        hfSumPlus_37_42          += E;
        maxHFEnergyPlus_37_42     = std::max(maxHFEnergyPlus_37_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.2f && E >= Emin) {
        hfSumMinus_37_42         += E;
        maxHFEnergyMinus_37_42    = std::max(maxHFEnergyMinus_37_42, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.3)
      if (η >= 3.7f && η <= 4.3f && E >= Emin) {
        hfSumPlus_37_43          += E;
        maxHFEnergyPlus_37_43     = std::max(maxHFEnergyPlus_37_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.3f && E >= Emin) {
        hfSumMinus_37_43         += E;
        maxHFEnergyMinus_37_43    = std::max(maxHFEnergyMinus_37_43, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.4)
      if (η >= 3.7f && η <= 4.4f && E >= Emin) {
        hfSumPlus_37_44          += E;
        maxHFEnergyPlus_37_44     = std::max(maxHFEnergyPlus_37_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.4f && E >= Emin) {
        hfSumMinus_37_44         += E;
        maxHFEnergyMinus_37_44    = std::max(maxHFEnergyMinus_37_44, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.5)
      if (η >= 3.7f && η <= 4.5f && E >= Emin) {
        hfSumPlus_37_45          += E;
        maxHFEnergyPlus_37_45     = std::max(maxHFEnergyPlus_37_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.5f && E >= Emin) {
        hfSumMinus_37_45         += E;
        maxHFEnergyMinus_37_45    = std::max(maxHFEnergyMinus_37_45, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.6)
      if (η >= 3.7f && η <= 4.6f && E >= Emin) {
        hfSumPlus_37_46          += E;
        maxHFEnergyPlus_37_46     = std::max(maxHFEnergyPlus_37_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.6f && E >= Emin) {
        hfSumMinus_37_46         += E;
        maxHFEnergyMinus_37_46    = std::max(maxHFEnergyMinus_37_46, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.7)
      if (η >= 3.7f && η <= 4.7f && E >= Emin) {
        hfSumPlus_37_47          += E;
        maxHFEnergyPlus_37_47     = std::max(maxHFEnergyPlus_37_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.7f && E >= Emin) {
        hfSumMinus_37_47         += E;
        maxHFEnergyMinus_37_47    = std::max(maxHFEnergyMinus_37_47, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.8)
      if (η >= 3.7f && η <= 4.8f && E >= Emin) {
        hfSumPlus_37_48          += E;
        maxHFEnergyPlus_37_48     = std::max(maxHFEnergyPlus_37_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.8f && E >= Emin) {
        hfSumMinus_37_48         += E;
        maxHFEnergyMinus_37_48    = std::max(maxHFEnergyMinus_37_48, E);
      }

      // HF+ region (3.7 ≤ η ≤ 4.9)
      if (η >= 3.7f && η <= 4.9f && E >= Emin) {
        hfSumPlus_37_49          += E;
        maxHFEnergyPlus_37_49     = std::max(maxHFEnergyPlus_37_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -4.9f && E >= Emin) {
        hfSumMinus_37_49         += E;
        maxHFEnergyMinus_37_49    = std::max(maxHFEnergyMinus_37_49, E);
      }

      // HF+ region (3.7 ≤ η ≤ 5.0)
      if (η >= 3.7f && η <= 5.0f && E >= Emin) {
        hfSumPlus_37_50          += E;
        maxHFEnergyPlus_37_50     = std::max(maxHFEnergyPlus_37_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -5.0f && E >= Emin) {
        hfSumMinus_37_50         += E;
        maxHFEnergyMinus_37_50    = std::max(maxHFEnergyMinus_37_50, E);
      }

      // HF+ region (3.7 ≤ η ≤ 5.1)
      if (η >= 3.7f && η <= 5.1f && E >= Emin) {
        hfSumPlus_37_51          += E;
        maxHFEnergyPlus_37_51     = std::max(maxHFEnergyPlus_37_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -5.1f && E >= Emin) {
        hfSumMinus_37_51         += E;
        maxHFEnergyMinus_37_51    = std::max(maxHFEnergyMinus_37_51, E);
      }

      // HF+ region (3.7 ≤ η ≤ 5.2)
      if (η >= 3.7f && η <= 5.2f && E >= Emin) {
        hfSumPlus_37_52          += E;
        maxHFEnergyPlus_37_52     = std::max(maxHFEnergyPlus_37_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.7)
      if (η <= -3.7f && η >= -5.2f && E >= Emin) {
        hfSumMinus_37_52         += E;
        maxHFEnergyMinus_37_52    = std::max(maxHFEnergyMinus_37_52, E);
      }



      // HF+ region (3.8 ≤ η ≤ 3.9)
      if (η >= 3.8f && η <= 3.9f && E >= Emin) {
        hfSumPlus_38_39          += E;
        maxHFEnergyPlus_38_39     = std::max(maxHFEnergyPlus_38_39, E);
      }
      // HF– region (-3.9 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -3.9f && E >= Emin) {
        hfSumMinus_38_39         += E;
        maxHFEnergyMinus_38_39    = std::max(maxHFEnergyMinus_38_39, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.0)
      if (η >= 3.8f && η <= 4.0f && E >= Emin) {
        hfSumPlus_38_40          += E;
        maxHFEnergyPlus_38_40     = std::max(maxHFEnergyPlus_38_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.0f && E >= Emin) {
        hfSumMinus_38_40         += E;
        maxHFEnergyMinus_38_40    = std::max(maxHFEnergyMinus_38_40, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.1)
      if (η >= 3.8f && η <= 4.1f && E >= Emin) {
        hfSumPlus_38_41          += E;
        maxHFEnergyPlus_38_41     = std::max(maxHFEnergyPlus_38_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.1f && E >= Emin) {
        hfSumMinus_38_41         += E;
        maxHFEnergyMinus_38_41    = std::max(maxHFEnergyMinus_38_41, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.2)
      if (η >= 3.8f && η <= 4.2f && E >= Emin) {
        hfSumPlus_38_42          += E;
        maxHFEnergyPlus_38_42     = std::max(maxHFEnergyPlus_38_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.2f && E >= Emin) {
        hfSumMinus_38_42         += E;
        maxHFEnergyMinus_38_42    = std::max(maxHFEnergyMinus_38_42, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.3)
      if (η >= 3.8f && η <= 4.3f && E >= Emin) {
        hfSumPlus_38_43          += E;
        maxHFEnergyPlus_38_43     = std::max(maxHFEnergyPlus_38_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.3f && E >= Emin) {
        hfSumMinus_38_43         += E;
        maxHFEnergyMinus_38_43    = std::max(maxHFEnergyMinus_38_43, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.4)
      if (η >= 3.8f && η <= 4.4f && E >= Emin) {
        hfSumPlus_38_44          += E;
        maxHFEnergyPlus_38_44     = std::max(maxHFEnergyPlus_38_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.4f && E >= Emin) {
        hfSumMinus_38_44         += E;
        maxHFEnergyMinus_38_44    = std::max(maxHFEnergyMinus_38_44, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.5)
      if (η >= 3.8f && η <= 4.5f && E >= Emin) {
        hfSumPlus_38_45          += E;
        maxHFEnergyPlus_38_45     = std::max(maxHFEnergyPlus_38_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.5f && E >= Emin) {
        hfSumMinus_38_45         += E;
        maxHFEnergyMinus_38_45    = std::max(maxHFEnergyMinus_38_45, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.6)
      if (η >= 3.8f && η <= 4.6f && E >= Emin) {
        hfSumPlus_38_46          += E;
        maxHFEnergyPlus_38_46     = std::max(maxHFEnergyPlus_38_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.6f && E >= Emin) {
        hfSumMinus_38_46         += E;
        maxHFEnergyMinus_38_46    = std::max(maxHFEnergyMinus_38_46, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.7)
      if (η >= 3.8f && η <= 4.7f && E >= Emin) {
        hfSumPlus_38_47          += E;
        maxHFEnergyPlus_38_47     = std::max(maxHFEnergyPlus_38_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.7f && E >= Emin) {
        hfSumMinus_38_47         += E;
        maxHFEnergyMinus_38_47    = std::max(maxHFEnergyMinus_38_47, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.8)
      if (η >= 3.8f && η <= 4.8f && E >= Emin) {
        hfSumPlus_38_48          += E;
        maxHFEnergyPlus_38_48     = std::max(maxHFEnergyPlus_38_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.8f && E >= Emin) {
        hfSumMinus_38_48         += E;
        maxHFEnergyMinus_38_48    = std::max(maxHFEnergyMinus_38_48, E);
      }

      // HF+ region (3.8 ≤ η ≤ 4.9)
      if (η >= 3.8f && η <= 4.9f && E >= Emin) {
        hfSumPlus_38_49          += E;
        maxHFEnergyPlus_38_49     = std::max(maxHFEnergyPlus_38_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -4.9f && E >= Emin) {
        hfSumMinus_38_49         += E;
        maxHFEnergyMinus_38_49    = std::max(maxHFEnergyMinus_38_49, E);
      }

      // HF+ region (3.8 ≤ η ≤ 5.0)
      if (η >= 3.8f && η <= 5.0f && E >= Emin) {
        hfSumPlus_38_50          += E;
        maxHFEnergyPlus_38_50     = std::max(maxHFEnergyPlus_38_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -5.0f && E >= Emin) {
        hfSumMinus_38_50         += E;
        maxHFEnergyMinus_38_50    = std::max(maxHFEnergyMinus_38_50, E);
      }

      // HF+ region (3.8 ≤ η ≤ 5.1)
      if (η >= 3.8f && η <= 5.1f && E >= Emin) {
        hfSumPlus_38_51          += E;
        maxHFEnergyPlus_38_51     = std::max(maxHFEnergyPlus_38_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -5.1f && E >= Emin) {
        hfSumMinus_38_51         += E;
        maxHFEnergyMinus_38_51    = std::max(maxHFEnergyMinus_38_51, E);
      }

      // HF+ region (3.8 ≤ η ≤ 5.2)
      if (η >= 3.8f && η <= 5.2f && E >= Emin) {
        hfSumPlus_38_52          += E;
        maxHFEnergyPlus_38_52     = std::max(maxHFEnergyPlus_38_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.8)
      if (η <= -3.8f && η >= -5.2f && E >= Emin) {
        hfSumMinus_38_52         += E;
        maxHFEnergyMinus_38_52    = std::max(maxHFEnergyMinus_38_52, E);
      }



      // HF+ region (3.9 ≤ η ≤ 4.0)
      if (η >= 3.9f && η <= 4.0f && E >= Emin) {
        hfSumPlus_39_40          += E;
        maxHFEnergyPlus_39_40     = std::max(maxHFEnergyPlus_39_40, E);
      }
      // HF– region (-4.0 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.0f && E >= Emin) {
        hfSumMinus_39_40         += E;
        maxHFEnergyMinus_39_40    = std::max(maxHFEnergyMinus_39_40, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.1)
      if (η >= 3.9f && η <= 4.1f && E >= Emin) {
        hfSumPlus_39_41          += E;
        maxHFEnergyPlus_39_41     = std::max(maxHFEnergyPlus_39_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.1f && E >= Emin) {
        hfSumMinus_39_41         += E;
        maxHFEnergyMinus_39_41    = std::max(maxHFEnergyMinus_39_41, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.2)
      if (η >= 3.9f && η <= 4.2f && E >= Emin) {
        hfSumPlus_39_42          += E;
        maxHFEnergyPlus_39_42     = std::max(maxHFEnergyPlus_39_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.2f && E >= Emin) {
        hfSumMinus_39_42         += E;
        maxHFEnergyMinus_39_42    = std::max(maxHFEnergyMinus_39_42, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.3)
      if (η >= 3.9f && η <= 4.3f && E >= Emin) {
        hfSumPlus_39_43          += E;
        maxHFEnergyPlus_39_43     = std::max(maxHFEnergyPlus_39_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.3f && E >= Emin) {
        hfSumMinus_39_43         += E;
        maxHFEnergyMinus_39_43    = std::max(maxHFEnergyMinus_39_43, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.4)
      if (η >= 3.9f && η <= 4.4f && E >= Emin) {
        hfSumPlus_39_44          += E;
        maxHFEnergyPlus_39_44     = std::max(maxHFEnergyPlus_39_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.4f && E >= Emin) {
        hfSumMinus_39_44         += E;
        maxHFEnergyMinus_39_44    = std::max(maxHFEnergyMinus_39_44, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.5)
      if (η >= 3.9f && η <= 4.5f && E >= Emin) {
        hfSumPlus_39_45          += E;
        maxHFEnergyPlus_39_45     = std::max(maxHFEnergyPlus_39_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.5f && E >= Emin) {
        hfSumMinus_39_45         += E;
        maxHFEnergyMinus_39_45    = std::max(maxHFEnergyMinus_39_45, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.6)
      if (η >= 3.9f && η <= 4.6f && E >= Emin) {
        hfSumPlus_39_46          += E;
        maxHFEnergyPlus_39_46     = std::max(maxHFEnergyPlus_39_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.6f && E >= Emin) {
        hfSumMinus_39_46         += E;
        maxHFEnergyMinus_39_46    = std::max(maxHFEnergyMinus_39_46, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.7)
      if (η >= 3.9f && η <= 4.7f && E >= Emin) {
        hfSumPlus_39_47          += E;
        maxHFEnergyPlus_39_47     = std::max(maxHFEnergyPlus_39_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.7f && E >= Emin) {
        hfSumMinus_39_47         += E;
        maxHFEnergyMinus_39_47    = std::max(maxHFEnergyMinus_39_47, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.8)
      if (η >= 3.9f && η <= 4.8f && E >= Emin) {
        hfSumPlus_39_48          += E;
        maxHFEnergyPlus_39_48     = std::max(maxHFEnergyPlus_39_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.8f && E >= Emin) {
        hfSumMinus_39_48         += E;
        maxHFEnergyMinus_39_48    = std::max(maxHFEnergyMinus_39_48, E);
      }

      // HF+ region (3.9 ≤ η ≤ 4.9)
      if (η >= 3.9f && η <= 4.9f && E >= Emin) {
        hfSumPlus_39_49          += E;
        maxHFEnergyPlus_39_49     = std::max(maxHFEnergyPlus_39_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -4.9f && E >= Emin) {
        hfSumMinus_39_49         += E;
        maxHFEnergyMinus_39_49    = std::max(maxHFEnergyMinus_39_49, E);
      }

      // HF+ region (3.9 ≤ η ≤ 5.0)
      if (η >= 3.9f && η <= 5.0f && E >= Emin) {
        hfSumPlus_39_50          += E;
        maxHFEnergyPlus_39_50     = std::max(maxHFEnergyPlus_39_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -5.0f && E >= Emin) {
        hfSumMinus_39_50         += E;
        maxHFEnergyMinus_39_50    = std::max(maxHFEnergyMinus_39_50, E);
      }

      // HF+ region (3.9 ≤ η ≤ 5.1)
      if (η >= 3.9f && η <= 5.1f && E >= Emin) {
        hfSumPlus_39_51          += E;
        maxHFEnergyPlus_39_51     = std::max(maxHFEnergyPlus_39_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -5.1f && E >= Emin) {
        hfSumMinus_39_51         += E;
        maxHFEnergyMinus_39_51    = std::max(maxHFEnergyMinus_39_51, E);
      }

      // HF+ region (3.9 ≤ η ≤ 5.2)
      if (η >= 3.9f && η <= 5.2f && E >= Emin) {
        hfSumPlus_39_52          += E;
        maxHFEnergyPlus_39_52     = std::max(maxHFEnergyPlus_39_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –3.9)
      if (η <= -3.9f && η >= -5.2f && E >= Emin) {
        hfSumMinus_39_52         += E;
        maxHFEnergyMinus_39_52    = std::max(maxHFEnergyMinus_39_52, E);
      }



      // HF+ region (4.0 ≤ η ≤ 4.1)
      if (η >= 4.0f && η <= 4.1f && E >= Emin) {
        hfSumPlus_40_41          += E;
        maxHFEnergyPlus_40_41     = std::max(maxHFEnergyPlus_40_41, E);
      }
      // HF– region (-4.1 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.1f && E >= Emin) {
        hfSumMinus_40_41         += E;
        maxHFEnergyMinus_40_41    = std::max(maxHFEnergyMinus_40_41, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.2)
      if (η >= 4.0f && η <= 4.2f && E >= Emin) {
        hfSumPlus_40_42          += E;
        maxHFEnergyPlus_40_42     = std::max(maxHFEnergyPlus_40_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.2f && E >= Emin) {
        hfSumMinus_40_42         += E;
        maxHFEnergyMinus_40_42    = std::max(maxHFEnergyMinus_40_42, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.3)
      if (η >= 4.0f && η <= 4.3f && E >= Emin) {
        hfSumPlus_40_43          += E;
        maxHFEnergyPlus_40_43     = std::max(maxHFEnergyPlus_40_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.3f && E >= Emin) {
        hfSumMinus_40_43         += E;
        maxHFEnergyMinus_40_43    = std::max(maxHFEnergyMinus_40_43, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.4)
      if (η >= 4.0f && η <= 4.4f && E >= Emin) {
        hfSumPlus_40_44          += E;
        maxHFEnergyPlus_40_44     = std::max(maxHFEnergyPlus_40_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.4f && E >= Emin) {
        hfSumMinus_40_44         += E;
        maxHFEnergyMinus_40_44    = std::max(maxHFEnergyMinus_40_44, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.5)
      if (η >= 4.0f && η <= 4.5f && E >= Emin) {
        hfSumPlus_40_45          += E;
        maxHFEnergyPlus_40_45     = std::max(maxHFEnergyPlus_40_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.5f && E >= Emin) {
        hfSumMinus_40_45         += E;
        maxHFEnergyMinus_40_45    = std::max(maxHFEnergyMinus_40_45, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.6)
      if (η >= 4.0f && η <= 4.6f && E >= Emin) {
        hfSumPlus_40_46          += E;
        maxHFEnergyPlus_40_46     = std::max(maxHFEnergyPlus_40_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.6f && E >= Emin) {
        hfSumMinus_40_46         += E;
        maxHFEnergyMinus_40_46    = std::max(maxHFEnergyMinus_40_46, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.7)
      if (η >= 4.0f && η <= 4.7f && E >= Emin) {
        hfSumPlus_40_47          += E;
        maxHFEnergyPlus_40_47     = std::max(maxHFEnergyPlus_40_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.7f && E >= Emin) {
        hfSumMinus_40_47         += E;
        maxHFEnergyMinus_40_47    = std::max(maxHFEnergyMinus_40_47, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.8)
      if (η >= 4.0f && η <= 4.8f && E >= Emin) {
        hfSumPlus_40_48          += E;
        maxHFEnergyPlus_40_48     = std::max(maxHFEnergyPlus_40_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.8f && E >= Emin) {
        hfSumMinus_40_48         += E;
        maxHFEnergyMinus_40_48    = std::max(maxHFEnergyMinus_40_48, E);
      }

      // HF+ region (4.0 ≤ η ≤ 4.9)
      if (η >= 4.0f && η <= 4.9f && E >= Emin) {
        hfSumPlus_40_49          += E;
        maxHFEnergyPlus_40_49     = std::max(maxHFEnergyPlus_40_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -4.9f && E >= Emin) {
        hfSumMinus_40_49         += E;
        maxHFEnergyMinus_40_49    = std::max(maxHFEnergyMinus_40_49, E);
      }

      // HF+ region (4.0 ≤ η ≤ 5.0)
      if (η >= 4.0f && η <= 5.0f && E >= Emin) {
        hfSumPlus_40_50          += E;
        maxHFEnergyPlus_40_50     = std::max(maxHFEnergyPlus_40_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -5.0f && E >= Emin) {
        hfSumMinus_40_50         += E;
        maxHFEnergyMinus_40_50    = std::max(maxHFEnergyMinus_40_50, E);
      }

      // HF+ region (4.0 ≤ η ≤ 5.1)
      if (η >= 4.0f && η <= 5.1f && E >= Emin) {
        hfSumPlus_40_51          += E;
        maxHFEnergyPlus_40_51     = std::max(maxHFEnergyPlus_40_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -5.1f && E >= Emin) {
        hfSumMinus_40_51         += E;
        maxHFEnergyMinus_40_51    = std::max(maxHFEnergyMinus_40_51, E);
      }

      // HF+ region (4.0 ≤ η ≤ 5.2)
      if (η >= 4.0f && η <= 5.2f && E >= Emin) {
        hfSumPlus_40_52          += E;
        maxHFEnergyPlus_40_52     = std::max(maxHFEnergyPlus_40_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.0)
      if (η <= -4.0f && η >= -5.2f && E >= Emin) {
        hfSumMinus_40_52         += E;
        maxHFEnergyMinus_40_52    = std::max(maxHFEnergyMinus_40_52, E);
      }



      // HF+ region (4.1 ≤ η ≤ 4.2)
      if (η >= 4.1f && η <= 4.2f && E >= Emin) {
        hfSumPlus_41_42          += E;
        maxHFEnergyPlus_41_42     = std::max(maxHFEnergyPlus_41_42, E);
      }
      // HF– region (-4.2 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.2f && E >= Emin) {
        hfSumMinus_41_42         += E;
        maxHFEnergyMinus_41_42    = std::max(maxHFEnergyMinus_41_42, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.3)
      if (η >= 4.1f && η <= 4.3f && E >= Emin) {
        hfSumPlus_41_43          += E;
        maxHFEnergyPlus_41_43     = std::max(maxHFEnergyPlus_41_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.3f && E >= Emin) {
        hfSumMinus_41_43         += E;
        maxHFEnergyMinus_41_43    = std::max(maxHFEnergyMinus_41_43, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.4)
      if (η >= 4.1f && η <= 4.4f && E >= Emin) {
        hfSumPlus_41_44          += E;
        maxHFEnergyPlus_41_44     = std::max(maxHFEnergyPlus_41_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.4f && E >= Emin) {
        hfSumMinus_41_44         += E;
        maxHFEnergyMinus_41_44    = std::max(maxHFEnergyMinus_41_44, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.5)
      if (η >= 4.1f && η <= 4.5f && E >= Emin) {
        hfSumPlus_41_45          += E;
        maxHFEnergyPlus_41_45     = std::max(maxHFEnergyPlus_41_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.5f && E >= Emin) {
        hfSumMinus_41_45         += E;
        maxHFEnergyMinus_41_45    = std::max(maxHFEnergyMinus_41_45, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.6)
      if (η >= 4.1f && η <= 4.6f && E >= Emin) {
        hfSumPlus_41_46          += E;
        maxHFEnergyPlus_41_46     = std::max(maxHFEnergyPlus_41_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.6f && E >= Emin) {
        hfSumMinus_41_46         += E;
        maxHFEnergyMinus_41_46    = std::max(maxHFEnergyMinus_41_46, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.7)
      if (η >= 4.1f && η <= 4.7f && E >= Emin) {
        hfSumPlus_41_47          += E;
        maxHFEnergyPlus_41_47     = std::max(maxHFEnergyPlus_41_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.7f && E >= Emin) {
        hfSumMinus_41_47         += E;
        maxHFEnergyMinus_41_47    = std::max(maxHFEnergyMinus_41_47, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.8)
      if (η >= 4.1f && η <= 4.8f && E >= Emin) {
        hfSumPlus_41_48          += E;
        maxHFEnergyPlus_41_48     = std::max(maxHFEnergyPlus_41_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.8f && E >= Emin) {
        hfSumMinus_41_48         += E;
        maxHFEnergyMinus_41_48    = std::max(maxHFEnergyMinus_41_48, E);
      }

      // HF+ region (4.1 ≤ η ≤ 4.9)
      if (η >= 4.1f && η <= 4.9f && E >= Emin) {
        hfSumPlus_41_49          += E;
        maxHFEnergyPlus_41_49     = std::max(maxHFEnergyPlus_41_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -4.9f && E >= Emin) {
        hfSumMinus_41_49         += E;
        maxHFEnergyMinus_41_49    = std::max(maxHFEnergyMinus_41_49, E);
      }

      // HF+ region (4.1 ≤ η ≤ 5.0)
      if (η >= 4.1f && η <= 5.0f && E >= Emin) {
        hfSumPlus_41_50          += E;
        maxHFEnergyPlus_41_50     = std::max(maxHFEnergyPlus_41_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -5.0f && E >= Emin) {
        hfSumMinus_41_50         += E;
        maxHFEnergyMinus_41_50    = std::max(maxHFEnergyMinus_41_50, E);
      }

      // HF+ region (4.1 ≤ η ≤ 5.1)
      if (η >= 4.1f && η <= 5.1f && E >= Emin) {
        hfSumPlus_41_51          += E;
        maxHFEnergyPlus_41_51     = std::max(maxHFEnergyPlus_41_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -5.1f && E >= Emin) {
        hfSumMinus_41_51         += E;
        maxHFEnergyMinus_41_51    = std::max(maxHFEnergyMinus_41_51, E);
      }

      // HF+ region (4.1 ≤ η ≤ 5.2)
      if (η >= 4.1f && η <= 5.2f && E >= Emin) {
        hfSumPlus_41_52          += E;
        maxHFEnergyPlus_41_52     = std::max(maxHFEnergyPlus_41_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.1)
      if (η <= -4.1f && η >= -5.2f && E >= Emin) {
        hfSumMinus_41_52         += E;
        maxHFEnergyMinus_41_52    = std::max(maxHFEnergyMinus_41_52, E);
      }



      // HF+ region (4.2 ≤ η ≤ 4.3)
      if (η >= 4.2f && η <= 4.3f && E >= Emin) {
        hfSumPlus_42_43          += E;
        maxHFEnergyPlus_42_43     = std::max(maxHFEnergyPlus_42_43, E);
      }
      // HF– region (-4.3 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.3f && E >= Emin) {
        hfSumMinus_42_43         += E;
        maxHFEnergyMinus_42_43    = std::max(maxHFEnergyMinus_42_43, E);
      }

      // HF+ region (4.2 ≤ η ≤ 4.4)
      if (η >= 4.2f && η <= 4.4f && E >= Emin) {
        hfSumPlus_42_44          += E;
        maxHFEnergyPlus_42_44     = std::max(maxHFEnergyPlus_42_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.4f && E >= Emin) {
        hfSumMinus_42_44         += E;
        maxHFEnergyMinus_42_44    = std::max(maxHFEnergyMinus_42_44, E);
      }

      // HF+ region (4.2 ≤ η ≤ 4.5)
      if (η >= 4.2f && η <= 4.5f && E >= Emin) {
        hfSumPlus_42_45          += E;
        maxHFEnergyPlus_42_45     = std::max(maxHFEnergyPlus_42_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.5f && E >= Emin) {
        hfSumMinus_42_45         += E;
        maxHFEnergyMinus_42_45    = std::max(maxHFEnergyMinus_42_45, E);
      }

      // HF+ region (4.2 ≤ η ≤ 4.6)
      if (η >= 4.2f && η <= 4.6f && E >= Emin) {
        hfSumPlus_42_46          += E;
        maxHFEnergyPlus_42_46     = std::max(maxHFEnergyPlus_42_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.6f && E >= Emin) {
        hfSumMinus_42_46         += E;
        maxHFEnergyMinus_42_46    = std::max(maxHFEnergyMinus_42_46, E);
      }

      // HF+ region (4.2 ≤ η ≤ 4.7)
      if (η >= 4.2f && η <= 4.7f && E >= Emin) {
        hfSumPlus_42_47          += E;
        maxHFEnergyPlus_42_47     = std::max(maxHFEnergyPlus_42_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.7f && E >= Emin) {
        hfSumMinus_42_47         += E;
        maxHFEnergyMinus_42_47    = std::max(maxHFEnergyMinus_42_47, E);
      }

      // HF+ region (4.2 ≤ η ≤ 4.8)
      if (η >= 4.2f && η <= 4.8f && E >= Emin) {
        hfSumPlus_42_48          += E;
        maxHFEnergyPlus_42_48     = std::max(maxHFEnergyPlus_42_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.8f && E >= Emin) {
        hfSumMinus_42_48         += E;
        maxHFEnergyMinus_42_48    = std::max(maxHFEnergyMinus_42_48, E);
      }

      // HF+ region (4.2 ≤ η ≤ 4.9)
      if (η >= 4.2f && η <= 4.9f && E >= Emin) {
        hfSumPlus_42_49          += E;
        maxHFEnergyPlus_42_49     = std::max(maxHFEnergyPlus_42_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -4.9f && E >= Emin) {
        hfSumMinus_42_49         += E;
        maxHFEnergyMinus_42_49    = std::max(maxHFEnergyMinus_42_49, E);
      }

      // HF+ region (4.2 ≤ η ≤ 5.0)
      if (η >= 4.2f && η <= 5.0f && E >= Emin) {
        hfSumPlus_42_50          += E;
        maxHFEnergyPlus_42_50     = std::max(maxHFEnergyPlus_42_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -5.0f && E >= Emin) {
        hfSumMinus_42_50         += E;
        maxHFEnergyMinus_42_50    = std::max(maxHFEnergyMinus_42_50, E);
      }

      // HF+ region (4.2 ≤ η ≤ 5.1)
      if (η >= 4.2f && η <= 5.1f && E >= Emin) {
        hfSumPlus_42_51          += E;
        maxHFEnergyPlus_42_51     = std::max(maxHFEnergyPlus_42_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -5.1f && E >= Emin) {
        hfSumMinus_42_51         += E;
        maxHFEnergyMinus_42_51    = std::max(maxHFEnergyMinus_42_51, E);
      }

      // HF+ region (4.2 ≤ η ≤ 5.2)
      if (η >= 4.2f && η <= 5.2f && E >= Emin) {
        hfSumPlus_42_52          += E;
        maxHFEnergyPlus_42_52     = std::max(maxHFEnergyPlus_42_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.2)
      if (η <= -4.2f && η >= -5.2f && E >= Emin) {
        hfSumMinus_42_52         += E;
        maxHFEnergyMinus_42_52    = std::max(maxHFEnergyMinus_42_52, E);
      }



      // HF+ region (4.3 ≤ η ≤ 4.4)
      if (η >= 4.3f && η <= 4.4f && E >= Emin) {
        hfSumPlus_43_44          += E;
        maxHFEnergyPlus_43_44     = std::max(maxHFEnergyPlus_43_44, E);
      }
      // HF– region (-4.4 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -4.4f && E >= Emin) {
        hfSumMinus_43_44         += E;
        maxHFEnergyMinus_43_44    = std::max(maxHFEnergyMinus_43_44, E);
      }

      // HF+ region (4.3 ≤ η ≤ 4.5)
      if (η >= 4.3f && η <= 4.5f && E >= Emin) {
        hfSumPlus_43_45          += E;
        maxHFEnergyPlus_43_45     = std::max(maxHFEnergyPlus_43_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -4.5f && E >= Emin) {
        hfSumMinus_43_45         += E;
        maxHFEnergyMinus_43_45    = std::max(maxHFEnergyMinus_43_45, E);
      }

      // HF+ region (4.3 ≤ η ≤ 4.6)
      if (η >= 4.3f && η <= 4.6f && E >= Emin) {
        hfSumPlus_43_46          += E;
        maxHFEnergyPlus_43_46     = std::max(maxHFEnergyPlus_43_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -4.6f && E >= Emin) {
        hfSumMinus_43_46         += E;
        maxHFEnergyMinus_43_46    = std::max(maxHFEnergyMinus_43_46, E);
      }

      // HF+ region (4.3 ≤ η ≤ 4.7)
      if (η >= 4.3f && η <= 4.7f && E >= Emin) {
        hfSumPlus_43_47          += E;
        maxHFEnergyPlus_43_47     = std::max(maxHFEnergyPlus_43_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -4.7f && E >= Emin) {
        hfSumMinus_43_47         += E;
        maxHFEnergyMinus_43_47    = std::max(maxHFEnergyMinus_43_47, E);
      }

      // HF+ region (4.3 ≤ η ≤ 4.8)
      if (η >= 4.3f && η <= 4.8f && E >= Emin) {
        hfSumPlus_43_48          += E;
        maxHFEnergyPlus_43_48     = std::max(maxHFEnergyPlus_43_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -4.8f && E >= Emin) {
        hfSumMinus_43_48         += E;
        maxHFEnergyMinus_43_48    = std::max(maxHFEnergyMinus_43_48, E);
      }

      // HF+ region (4.3 ≤ η ≤ 4.9)
      if (η >= 4.3f && η <= 4.9f && E >= Emin) {
        hfSumPlus_43_49          += E;
        maxHFEnergyPlus_43_49     = std::max(maxHFEnergyPlus_43_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -4.9f && E >= Emin) {
        hfSumMinus_43_49         += E;
        maxHFEnergyMinus_43_49    = std::max(maxHFEnergyMinus_43_49, E);
      }

      // HF+ region (4.3 ≤ η ≤ 5.0)
      if (η >= 4.3f && η <= 5.0f && E >= Emin) {
        hfSumPlus_43_50          += E;
        maxHFEnergyPlus_43_50     = std::max(maxHFEnergyPlus_43_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -5.0f && E >= Emin) {
        hfSumMinus_43_50         += E;
        maxHFEnergyMinus_43_50    = std::max(maxHFEnergyMinus_43_50, E);
      }

      // HF+ region (4.3 ≤ η ≤ 5.1)
      if (η >= 4.3f && η <= 5.1f && E >= Emin) {
        hfSumPlus_43_51          += E;
        maxHFEnergyPlus_43_51     = std::max(maxHFEnergyPlus_43_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -5.1f && E >= Emin) {
        hfSumMinus_43_51         += E;
        maxHFEnergyMinus_43_51    = std::max(maxHFEnergyMinus_43_51, E);
      }

      // HF+ region (4.3 ≤ η ≤ 5.2)
      if (η >= 4.3f && η <= 5.2f && E >= Emin) {
        hfSumPlus_43_52          += E;
        maxHFEnergyPlus_43_52     = std::max(maxHFEnergyPlus_43_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.3)
      if (η <= -4.3f && η >= -5.2f && E >= Emin) {
        hfSumMinus_43_52         += E;
        maxHFEnergyMinus_43_52    = std::max(maxHFEnergyMinus_43_52, E);
      }



      // HF+ region (4.4 ≤ η ≤ 4.5)
      if (η >= 4.4f && η <= 4.5f && E >= Emin) {
        hfSumPlus_44_45          += E;
        maxHFEnergyPlus_44_45     = std::max(maxHFEnergyPlus_44_45, E);
      }
      // HF– region (-4.5 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -4.5f && E >= Emin) {
        hfSumMinus_44_45         += E;
        maxHFEnergyMinus_44_45    = std::max(maxHFEnergyMinus_44_45, E);
      }

      // HF+ region (4.4 ≤ η ≤ 4.6)
      if (η >= 4.4f && η <= 4.6f && E >= Emin) {
        hfSumPlus_44_46          += E;
        maxHFEnergyPlus_44_46     = std::max(maxHFEnergyPlus_44_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -4.6f && E >= Emin) {
        hfSumMinus_44_46         += E;
        maxHFEnergyMinus_44_46    = std::max(maxHFEnergyMinus_44_46, E);
      }

      // HF+ region (4.4 ≤ η ≤ 4.7)
      if (η >= 4.4f && η <= 4.7f && E >= Emin) {
        hfSumPlus_44_47          += E;
        maxHFEnergyPlus_44_47     = std::max(maxHFEnergyPlus_44_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -4.7f && E >= Emin) {
        hfSumMinus_44_47         += E;
        maxHFEnergyMinus_44_47    = std::max(maxHFEnergyMinus_44_47, E);
      }

      // HF+ region (4.4 ≤ η ≤ 4.8)
      if (η >= 4.4f && η <= 4.8f && E >= Emin) {
        hfSumPlus_44_48          += E;
        maxHFEnergyPlus_44_48     = std::max(maxHFEnergyPlus_44_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -4.8f && E >= Emin) {
        hfSumMinus_44_48         += E;
        maxHFEnergyMinus_44_48    = std::max(maxHFEnergyMinus_44_48, E);
      }

      // HF+ region (4.4 ≤ η ≤ 4.9)
      if (η >= 4.4f && η <= 4.9f && E >= Emin) {
        hfSumPlus_44_49          += E;
        maxHFEnergyPlus_44_49     = std::max(maxHFEnergyPlus_44_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -4.9f && E >= Emin) {
        hfSumMinus_44_49         += E;
        maxHFEnergyMinus_44_49    = std::max(maxHFEnergyMinus_44_49, E);
      }

      // HF+ region (4.4 ≤ η ≤ 5.0)
      if (η >= 4.4f && η <= 5.0f && E >= Emin) {
        hfSumPlus_44_50          += E;
        maxHFEnergyPlus_44_50     = std::max(maxHFEnergyPlus_44_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -5.0f && E >= Emin) {
        hfSumMinus_44_50         += E;
        maxHFEnergyMinus_44_50    = std::max(maxHFEnergyMinus_44_50, E);
      }

      // HF+ region (4.4 ≤ η ≤ 5.1)
      if (η >= 4.4f && η <= 5.1f && E >= Emin) {
        hfSumPlus_44_51          += E;
        maxHFEnergyPlus_44_51     = std::max(maxHFEnergyPlus_44_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -5.1f && E >= Emin) {
        hfSumMinus_44_51         += E;
        maxHFEnergyMinus_44_51    = std::max(maxHFEnergyMinus_44_51, E);
      }

      // HF+ region (4.4 ≤ η ≤ 5.2)
      if (η >= 4.4f && η <= 5.2f && E >= Emin) {
        hfSumPlus_44_52          += E;
        maxHFEnergyPlus_44_52     = std::max(maxHFEnergyPlus_44_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.4)
      if (η <= -4.4f && η >= -5.2f && E >= Emin) {
        hfSumMinus_44_52         += E;
        maxHFEnergyMinus_44_52    = std::max(maxHFEnergyMinus_44_52, E);
      }



      // HF+ region (4.5 ≤ η ≤ 4.6)
      if (η >= 4.5f && η <= 4.6f && E >= Emin) {
        hfSumPlus_45_46          += E;
        maxHFEnergyPlus_45_46     = std::max(maxHFEnergyPlus_45_46, E);
      }
      // HF– region (-4.6 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -4.6f && E >= Emin) {
        hfSumMinus_45_46         += E;
        maxHFEnergyMinus_45_46    = std::max(maxHFEnergyMinus_45_46, E);
      }

      // HF+ region (4.5 ≤ η ≤ 4.7)
      if (η >= 4.5f && η <= 4.7f && E >= Emin) {
        hfSumPlus_45_47          += E;
        maxHFEnergyPlus_45_47     = std::max(maxHFEnergyPlus_45_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -4.7f && E >= Emin) {
        hfSumMinus_45_47         += E;
        maxHFEnergyMinus_45_47    = std::max(maxHFEnergyMinus_45_47, E);
      }

      // HF+ region (4.5 ≤ η ≤ 4.8)
      if (η >= 4.5f && η <= 4.8f && E >= Emin) {
        hfSumPlus_45_48          += E;
        maxHFEnergyPlus_45_48     = std::max(maxHFEnergyPlus_45_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -4.8f && E >= Emin) {
        hfSumMinus_45_48         += E;
        maxHFEnergyMinus_45_48    = std::max(maxHFEnergyMinus_45_48, E);
      }

      // HF+ region (4.5 ≤ η ≤ 4.9)
      if (η >= 4.5f && η <= 4.9f && E >= Emin) {
        hfSumPlus_45_49          += E;
        maxHFEnergyPlus_45_49     = std::max(maxHFEnergyPlus_45_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -4.9f && E >= Emin) {
        hfSumMinus_45_49         += E;
        maxHFEnergyMinus_45_49    = std::max(maxHFEnergyMinus_45_49, E);
      }

      // HF+ region (4.5 ≤ η ≤ 5.0)
      if (η >= 4.5f && η <= 5.0f && E >= Emin) {
        hfSumPlus_45_50          += E;
        maxHFEnergyPlus_45_50     = std::max(maxHFEnergyPlus_45_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -5.0f && E >= Emin) {
        hfSumMinus_45_50         += E;
        maxHFEnergyMinus_45_50    = std::max(maxHFEnergyMinus_45_50, E);
      }

      // HF+ region (4.5 ≤ η ≤ 5.1)
      if (η >= 4.5f && η <= 5.1f && E >= Emin) {
        hfSumPlus_45_51          += E;
        maxHFEnergyPlus_45_51     = std::max(maxHFEnergyPlus_45_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -5.1f && E >= Emin) {
        hfSumMinus_45_51         += E;
        maxHFEnergyMinus_45_51    = std::max(maxHFEnergyMinus_45_51, E);
      }

      // HF+ region (4.5 ≤ η ≤ 5.2)
      if (η >= 4.5f && η <= 5.2f && E >= Emin) {
        hfSumPlus_45_52          += E;
        maxHFEnergyPlus_45_52     = std::max(maxHFEnergyPlus_45_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.5)
      if (η <= -4.5f && η >= -5.2f && E >= Emin) {
        hfSumMinus_45_52         += E;
        maxHFEnergyMinus_45_52    = std::max(maxHFEnergyMinus_45_52, E);
      }



      // HF+ region (4.6 ≤ η ≤ 4.7)
      if (η >= 4.6f && η <= 4.7f && E >= Emin) {
        hfSumPlus_46_47          += E;
        maxHFEnergyPlus_46_47     = std::max(maxHFEnergyPlus_46_47, E);
      }
      // HF– region (-4.7 ≤ η ≤ –4.6)
      if (η <= -4.6f && η >= -4.7f && E >= Emin) {
        hfSumMinus_46_47         += E;
        maxHFEnergyMinus_46_47    = std::max(maxHFEnergyMinus_46_47, E);
      }

      // HF+ region (4.6 ≤ η ≤ 4.8)
      if (η >= 4.6f && η <= 4.8f && E >= Emin) {
        hfSumPlus_46_48          += E;
        maxHFEnergyPlus_46_48     = std::max(maxHFEnergyPlus_46_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.6)
      if (η <= -4.6f && η >= -4.8f && E >= Emin) {
        hfSumMinus_46_48         += E;
        maxHFEnergyMinus_46_48    = std::max(maxHFEnergyMinus_46_48, E);
      }

      // HF+ region (4.6 ≤ η ≤ 4.9)
      if (η >= 4.6f && η <= 4.9f && E >= Emin) {
        hfSumPlus_46_49          += E;
        maxHFEnergyPlus_46_49     = std::max(maxHFEnergyPlus_46_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.6)
      if (η <= -4.6f && η >= -4.9f && E >= Emin) {
        hfSumMinus_46_49         += E;
        maxHFEnergyMinus_46_49    = std::max(maxHFEnergyMinus_46_49, E);
      }

      // HF+ region (4.6 ≤ η ≤ 5.0)
      if (η >= 4.6f && η <= 5.0f && E >= Emin) {
        hfSumPlus_46_50          += E;
        maxHFEnergyPlus_46_50     = std::max(maxHFEnergyPlus_46_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.6)
      if (η <= -4.6f && η >= -5.0f && E >= Emin) {
        hfSumMinus_46_50         += E;
        maxHFEnergyMinus_46_50    = std::max(maxHFEnergyMinus_46_50, E);
      }

      // HF+ region (4.6 ≤ η ≤ 5.1)
      if (η >= 4.6f && η <= 5.1f && E >= Emin) {
        hfSumPlus_46_51          += E;
        maxHFEnergyPlus_46_51     = std::max(maxHFEnergyPlus_46_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.6)
      if (η <= -4.6f && η >= -5.1f && E >= Emin) {
        hfSumMinus_46_51         += E;
        maxHFEnergyMinus_46_51    = std::max(maxHFEnergyMinus_46_51, E);
      }

      // HF+ region (4.6 ≤ η ≤ 5.2)
      if (η >= 4.6f && η <= 5.2f && E >= Emin) {
        hfSumPlus_46_52          += E;
        maxHFEnergyPlus_46_52     = std::max(maxHFEnergyPlus_46_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.6)
      if (η <= -4.6f && η >= -5.2f && E >= Emin) {
        hfSumMinus_46_52         += E;
        maxHFEnergyMinus_46_52    = std::max(maxHFEnergyMinus_46_52, E);
      }



      // HF+ region (4.7 ≤ η ≤ 4.8)
      if (η >= 4.7f && η <= 4.8f && E >= Emin) {
        hfSumPlus_47_48          += E;
        maxHFEnergyPlus_47_48     = std::max(maxHFEnergyPlus_47_48, E);
      }
      // HF– region (-4.8 ≤ η ≤ –4.7)
      if (η <= -4.7f && η >= -4.8f && E >= Emin) {
        hfSumMinus_47_48         += E;
        maxHFEnergyMinus_47_48    = std::max(maxHFEnergyMinus_47_48, E);
      }

      // HF+ region (4.7 ≤ η ≤ 4.9)
      if (η >= 4.7f && η <= 4.9f && E >= Emin) {
        hfSumPlus_47_49          += E;
        maxHFEnergyPlus_47_49     = std::max(maxHFEnergyPlus_47_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.7)
      if (η <= -4.7f && η >= -4.9f && E >= Emin) {
        hfSumMinus_47_49         += E;
        maxHFEnergyMinus_47_49    = std::max(maxHFEnergyMinus_47_49, E);
      }

      // HF+ region (4.7 ≤ η ≤ 5.0)
      if (η >= 4.7f && η <= 5.0f && E >= Emin) {
        hfSumPlus_47_50          += E;
        maxHFEnergyPlus_47_50     = std::max(maxHFEnergyPlus_47_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.7)
      if (η <= -4.7f && η >= -5.0f && E >= Emin) {
        hfSumMinus_47_50         += E;
        maxHFEnergyMinus_47_50    = std::max(maxHFEnergyMinus_47_50, E);
      }

      // HF+ region (4.7 ≤ η ≤ 5.1)
      if (η >= 4.7f && η <= 5.1f && E >= Emin) {
        hfSumPlus_47_51          += E;
        maxHFEnergyPlus_47_51     = std::max(maxHFEnergyPlus_47_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.7)
      if (η <= -4.7f && η >= -5.1f && E >= Emin) {
        hfSumMinus_47_51         += E;
        maxHFEnergyMinus_47_51    = std::max(maxHFEnergyMinus_47_51, E);
      }

      // HF+ region (4.7 ≤ η ≤ 5.2)
      if (η >= 4.7f && η <= 5.2f && E >= Emin) {
        hfSumPlus_47_52          += E;
        maxHFEnergyPlus_47_52     = std::max(maxHFEnergyPlus_47_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.7)
      if (η <= -4.7f && η >= -5.2f && E >= Emin) {
        hfSumMinus_47_52         += E;
        maxHFEnergyMinus_47_52    = std::max(maxHFEnergyMinus_47_52, E);
      }



      // HF+ region (4.8 ≤ η ≤ 4.9)
      if (η >= 4.8f && η <= 4.9f && E >= Emin) {
        hfSumPlus_48_49          += E;
        maxHFEnergyPlus_48_49     = std::max(maxHFEnergyPlus_48_49, E);
      }
      // HF– region (-4.9 ≤ η ≤ –4.8)
      if (η <= -4.8f && η >= -4.9f && E >= Emin) {
        hfSumMinus_48_49         += E;
        maxHFEnergyMinus_48_49    = std::max(maxHFEnergyMinus_48_49, E);
      }

      // HF+ region (4.8 ≤ η ≤ 5.0)
      if (η >= 4.8f && η <= 5.0f && E >= Emin) {
        hfSumPlus_48_50          += E;
        maxHFEnergyPlus_48_50     = std::max(maxHFEnergyPlus_48_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.8)
      if (η <= -4.8f && η >= -5.0f && E >= Emin) {
        hfSumMinus_48_50         += E;
        maxHFEnergyMinus_48_50    = std::max(maxHFEnergyMinus_48_50, E);
      }

      // HF+ region (4.8 ≤ η ≤ 5.1)
      if (η >= 4.8f && η <= 5.1f && E >= Emin) {
        hfSumPlus_48_51          += E;
        maxHFEnergyPlus_48_51     = std::max(maxHFEnergyPlus_48_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.8)
      if (η <= -4.8f && η >= -5.1f && E >= Emin) {
        hfSumMinus_48_51         += E;
        maxHFEnergyMinus_48_51    = std::max(maxHFEnergyMinus_48_51, E);
      }

      // HF+ region (4.8 ≤ η ≤ 5.2)
      if (η >= 4.8f && η <= 5.2f && E >= Emin) {
        hfSumPlus_48_52          += E;
        maxHFEnergyPlus_48_52     = std::max(maxHFEnergyPlus_48_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.8)
      if (η <= -4.8f && η >= -5.2f && E >= Emin) {
        hfSumMinus_48_52         += E;
        maxHFEnergyMinus_48_52    = std::max(maxHFEnergyMinus_48_52, E);
      }



      // HF+ region (4.9 ≤ η ≤ 5.0)
      if (η >= 4.9f && η <= 5.0f && E >= Emin) {
        hfSumPlus_49_50          += E;
        maxHFEnergyPlus_49_50     = std::max(maxHFEnergyPlus_49_50, E);
      }
      // HF– region (-5.0 ≤ η ≤ –4.9)
      if (η <= -4.9f && η >= -5.0f && E >= Emin) {
        hfSumMinus_49_50         += E;
        maxHFEnergyMinus_49_50    = std::max(maxHFEnergyMinus_49_50, E);
      }

      // HF+ region (4.9 ≤ η ≤ 5.1)
      if (η >= 4.9f && η <= 5.1f && E >= Emin) {
        hfSumPlus_49_51          += E;
        maxHFEnergyPlus_49_51     = std::max(maxHFEnergyPlus_49_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –4.9)
      if (η <= -4.9f && η >= -5.1f && E >= Emin) {
        hfSumMinus_49_51         += E;
        maxHFEnergyMinus_49_51    = std::max(maxHFEnergyMinus_49_51, E);
      }

      // HF+ region (4.9 ≤ η ≤ 5.2)
      if (η >= 4.9f && η <= 5.2f && E >= Emin) {
        hfSumPlus_49_52          += E;
        maxHFEnergyPlus_49_52     = std::max(maxHFEnergyPlus_49_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –4.9)
      if (η <= -4.9f && η >= -5.2f && E >= Emin) {
        hfSumMinus_49_52         += E;
        maxHFEnergyMinus_49_52    = std::max(maxHFEnergyMinus_49_52, E);
      }



      // HF+ region (5.0 ≤ η ≤ 5.1)
      if (η >= 5.0f && η <= 5.1f && E >= Emin) {
        hfSumPlus_50_51          += E;
        maxHFEnergyPlus_50_51     = std::max(maxHFEnergyPlus_50_51, E);
      }
      // HF– region (-5.1 ≤ η ≤ –5.0)
      if (η <= -5.0f && η >= -5.1f && E >= Emin) {
        hfSumMinus_50_51         += E;
        maxHFEnergyMinus_50_51    = std::max(maxHFEnergyMinus_50_51, E);
      }

      // HF+ region (5.0 ≤ η ≤ 5.2)
      if (η >= 5.0f && η <= 5.2f && E >= Emin) {
        hfSumPlus_50_52          += E;
        maxHFEnergyPlus_50_52     = std::max(maxHFEnergyPlus_50_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –5.0)
      if (η <= -5.0f && η >= -5.2f && E >= Emin) {
        hfSumMinus_50_52         += E;
        maxHFEnergyMinus_50_52    = std::max(maxHFEnergyMinus_50_52, E);
      }



      // HF+ region (5.1 ≤ η ≤ 5.2)
      if (η >= 5.1f && η <= 5.2f && E >= Emin) {
        hfSumPlus_51_52          += E;
        maxHFEnergyPlus_51_52     = std::max(maxHFEnergyPlus_51_52, E);
      }
      // HF– region (-5.2 ≤ η ≤ –5.1)
      if (η <= -5.1f && η >= -5.2f && E >= Emin) {
        hfSumMinus_51_52         += E;
        maxHFEnergyMinus_51_52    = std::max(maxHFEnergyMinus_51_52, E);
      }
      // RANGOS: eta --end

    
      // global HF extrema (any tower in |η|∈[3.0,5.2])
      if ((η >=  3.0f && η <=  5.2f) ||
          (η <= -3.0f && η >= -5.2f)) {
        maxHFEnergy = std::max(maxHFEnergy, E);
        minHFEnergy = std::min(minHFEnergy, E);
      }

      



       // global HF extrema (any tower in |η|∈[3.0,3.1])
      if ((η >=  3.0f && η <=  3.1f) ||
          (η <= -3.0f && η >= -3.1f)) {
        maxHFEnergy_30_31 = std::max(maxHFEnergy_30_31, E);
        minHFEnergy_30_31 = std::min(minHFEnergy_30_31, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.2])
      if ((η >=  3.0f && η <=  3.2f) ||
          (η <= -3.0f && η >= -3.2f)) {
        maxHFEnergy_30_32 = std::max(maxHFEnergy_30_32, E);
        minHFEnergy_30_32 = std::min(minHFEnergy_30_32, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.3])
      if ((η >=  3.0f && η <=  3.3f) ||
          (η <= -3.0f && η >= -3.3f)) {
        maxHFEnergy_30_33 = std::max(maxHFEnergy_30_33, E);
        minHFEnergy_30_33 = std::min(minHFEnergy_30_33, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.4])
      if ((η >=  3.0f && η <=  3.4f) ||
          (η <= -3.0f && η >= -3.4f)) {
        maxHFEnergy_30_34 = std::max(maxHFEnergy_30_34, E);
        minHFEnergy_30_34 = std::min(minHFEnergy_30_34, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.5])
      if ((η >=  3.0f && η <=  3.5f) ||
          (η <= -3.0f && η >= -3.5f)) {
        maxHFEnergy_30_35 = std::max(maxHFEnergy_30_35, E);
        minHFEnergy_30_35 = std::min(minHFEnergy_30_35, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.6])
      if ((η >=  3.0f && η <=  3.6f) ||
          (η <= -3.0f && η >= -3.6f)) {
        maxHFEnergy_30_36 = std::max(maxHFEnergy_30_36, E);
        minHFEnergy_30_36 = std::min(minHFEnergy_30_36, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.7])
      if ((η >=  3.0f && η <=  3.7f) ||
          (η <= -3.0f && η >= -3.7f)) {
        maxHFEnergy_30_37 = std::max(maxHFEnergy_30_37, E);
        minHFEnergy_30_37 = std::min(minHFEnergy_30_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.8])
      if ((η >=  3.0f && η <=  3.8f) ||
          (η <= -3.0f && η >= -3.8f)) {
        maxHFEnergy_30_38 = std::max(maxHFEnergy_30_38, E);
        minHFEnergy_30_38 = std::min(minHFEnergy_30_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,3.9])
      if ((η >=  3.0f && η <=  3.9f) ||
          (η <= -3.0f && η >= -3.9f)) {
        maxHFEnergy_30_39 = std::max(maxHFEnergy_30_39, E);
        minHFEnergy_30_39 = std::min(minHFEnergy_30_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.0])
      if ((η >=  3.0f && η <=  4.0f) ||
          (η <= -3.0f && η >= -4.0f)) {
        maxHFEnergy_30_40 = std::max(maxHFEnergy_30_40, E);
        minHFEnergy_30_40 = std::min(minHFEnergy_30_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.1])
      if ((η >=  3.0f && η <=  4.1f) ||
          (η <= -3.0f && η >= -4.1f)) {
        maxHFEnergy_30_41 = std::max(maxHFEnergy_30_41, E);
        minHFEnergy_30_41 = std::min(minHFEnergy_30_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.2])
      if ((η >=  3.0f && η <=  4.2f) ||
          (η <= -3.0f && η >= -4.2f)) {
        maxHFEnergy_30_42 = std::max(maxHFEnergy_30_42, E);
        minHFEnergy_30_42 = std::min(minHFEnergy_30_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.3])
      if ((η >=  3.0f && η <=  4.3f) ||
          (η <= -3.0f && η >= -4.3f)) {
        maxHFEnergy_30_43 = std::max(maxHFEnergy_30_43, E);
        minHFEnergy_30_43 = std::min(minHFEnergy_30_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.4])
      if ((η >=  3.0f && η <=  4.4f) ||
          (η <= -3.0f && η >= -4.4f)) {
        maxHFEnergy_30_44 = std::max(maxHFEnergy_30_44, E);
        minHFEnergy_30_44 = std::min(minHFEnergy_30_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.5])
      if ((η >=  3.0f && η <=  4.5f) ||
          (η <= -3.0f && η >= -4.5f)) {
        maxHFEnergy_30_45 = std::max(maxHFEnergy_30_45, E);
        minHFEnergy_30_45 = std::min(minHFEnergy_30_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.6])
      if ((η >=  3.0f && η <=  4.6f) ||
          (η <= -3.0f && η >= -4.6f)) {
        maxHFEnergy_30_46 = std::max(maxHFEnergy_30_46, E);
        minHFEnergy_30_46 = std::min(minHFEnergy_30_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.7])
      if ((η >=  3.0f && η <=  4.7f) ||
          (η <= -3.0f && η >= -4.7f)) {
        maxHFEnergy_30_47 = std::max(maxHFEnergy_30_47, E);
        minHFEnergy_30_47 = std::min(minHFEnergy_30_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.8])
      if ((η >=  3.0f && η <=  4.8f) ||
          (η <= -3.0f && η >= -4.8f)) {
        maxHFEnergy_30_48 = std::max(maxHFEnergy_30_48, E);
        minHFEnergy_30_48 = std::min(minHFEnergy_30_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,4.9])
      if ((η >=  3.0f && η <=  4.9f) ||
          (η <= -3.0f && η >= -4.9f)) {
        maxHFEnergy_30_49 = std::max(maxHFEnergy_30_49, E);
        minHFEnergy_30_49 = std::min(minHFEnergy_30_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,5.0])
      if ((η >=  3.0f && η <=  5.0f) ||
          (η <= -3.0f && η >= -5.0f)) {
        maxHFEnergy_30_50 = std::max(maxHFEnergy_30_50, E);
        minHFEnergy_30_50 = std::min(minHFEnergy_30_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,5.1])
      if ((η >=  3.0f && η <=  5.1f) ||
          (η <= -3.0f && η >= -5.1f)) {
        maxHFEnergy_30_51 = std::max(maxHFEnergy_30_51, E);
        minHFEnergy_30_51 = std::min(minHFEnergy_30_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.0,5.2])
      if ((η >=  3.0f && η <=  5.2f) ||
          (η <= -3.0f && η >= -5.2f)) {
        maxHFEnergy_30_52 = std::max(maxHFEnergy_30_52, E);
        minHFEnergy_30_52 = std::min(minHFEnergy_30_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.1,3.2])
      if ((η >=  3.1f && η <=  3.2f) ||
          (η <= -3.1f && η >= -3.2f)) {
        maxHFEnergy_31_32 = std::max(maxHFEnergy_31_32, E);
        minHFEnergy_31_32 = std::min(minHFEnergy_31_32, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.3])
      if ((η >=  3.1f && η <=  3.3f) ||
          (η <= -3.1f && η >= -3.3f)) {
        maxHFEnergy_31_33 = std::max(maxHFEnergy_31_33, E);
        minHFEnergy_31_33 = std::min(minHFEnergy_31_33, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.4])
      if ((η >=  3.1f && η <=  3.4f) ||
          (η <= -3.1f && η >= -3.4f)) {
        maxHFEnergy_31_34 = std::max(maxHFEnergy_31_34, E);
        minHFEnergy_31_34 = std::min(minHFEnergy_31_34, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.5])
      if ((η >=  3.1f && η <=  3.5f) ||
          (η <= -3.1f && η >= -3.5f)) {
        maxHFEnergy_31_35 = std::max(maxHFEnergy_31_35, E);
        minHFEnergy_31_35 = std::min(minHFEnergy_31_35, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.6])
      if ((η >=  3.1f && η <=  3.6f) ||
          (η <= -3.1f && η >= -3.6f)) {
        maxHFEnergy_31_36 = std::max(maxHFEnergy_31_36, E);
        minHFEnergy_31_36 = std::min(minHFEnergy_31_36, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.7])
      if ((η >=  3.1f && η <=  3.7f) ||
          (η <= -3.1f && η >= -3.7f)) {
        maxHFEnergy_31_37 = std::max(maxHFEnergy_31_37, E);
        minHFEnergy_31_37 = std::min(minHFEnergy_31_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.8])
      if ((η >=  3.1f && η <=  3.8f) ||
          (η <= -3.1f && η >= -3.8f)) {
        maxHFEnergy_31_38 = std::max(maxHFEnergy_31_38, E);
        minHFEnergy_31_38 = std::min(minHFEnergy_31_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,3.9])
      if ((η >=  3.1f && η <=  3.9f) ||
          (η <= -3.1f && η >= -3.9f)) {
        maxHFEnergy_31_39 = std::max(maxHFEnergy_31_39, E);
        minHFEnergy_31_39 = std::min(minHFEnergy_31_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.0])
      if ((η >=  3.1f && η <=  4.0f) ||
          (η <= -3.1f && η >= -4.0f)) {
        maxHFEnergy_31_40 = std::max(maxHFEnergy_31_40, E);
        minHFEnergy_31_40 = std::min(minHFEnergy_31_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.1])
      if ((η >=  3.1f && η <=  4.1f) ||
          (η <= -3.1f && η >= -4.1f)) {
        maxHFEnergy_31_41 = std::max(maxHFEnergy_31_41, E);
        minHFEnergy_31_41 = std::min(minHFEnergy_31_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.2])
      if ((η >=  3.1f && η <=  4.2f) ||
          (η <= -3.1f && η >= -4.2f)) {
        maxHFEnergy_31_42 = std::max(maxHFEnergy_31_42, E);
        minHFEnergy_31_42 = std::min(minHFEnergy_31_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.3])
      if ((η >=  3.1f && η <=  4.3f) ||
          (η <= -3.1f && η >= -4.3f)) {
        maxHFEnergy_31_43 = std::max(maxHFEnergy_31_43, E);
        minHFEnergy_31_43 = std::min(minHFEnergy_31_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.4])
      if ((η >=  3.1f && η <=  4.4f) ||
          (η <= -3.1f && η >= -4.4f)) {
        maxHFEnergy_31_44 = std::max(maxHFEnergy_31_44, E);
        minHFEnergy_31_44 = std::min(minHFEnergy_31_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.5])
      if ((η >=  3.1f && η <=  4.5f) ||
          (η <= -3.1f && η >= -4.5f)) {
        maxHFEnergy_31_45 = std::max(maxHFEnergy_31_45, E);
        minHFEnergy_31_45 = std::min(minHFEnergy_31_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.6])
      if ((η >=  3.1f && η <=  4.6f) ||
          (η <= -3.1f && η >= -4.6f)) {
        maxHFEnergy_31_46 = std::max(maxHFEnergy_31_46, E);
        minHFEnergy_31_46 = std::min(minHFEnergy_31_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.7])
      if ((η >=  3.1f && η <=  4.7f) ||
          (η <= -3.1f && η >= -4.7f)) {
        maxHFEnergy_31_47 = std::max(maxHFEnergy_31_47, E);
        minHFEnergy_31_47 = std::min(minHFEnergy_31_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.8])
      if ((η >=  3.1f && η <=  4.8f) ||
          (η <= -3.1f && η >= -4.8f)) {
        maxHFEnergy_31_48 = std::max(maxHFEnergy_31_48, E);
        minHFEnergy_31_48 = std::min(minHFEnergy_31_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,4.9])
      if ((η >=  3.1f && η <=  4.9f) ||
          (η <= -3.1f && η >= -4.9f)) {
        maxHFEnergy_31_49 = std::max(maxHFEnergy_31_49, E);
        minHFEnergy_31_49 = std::min(minHFEnergy_31_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,5.0])
      if ((η >=  3.1f && η <=  5.0f) ||
          (η <= -3.1f && η >= -5.0f)) {
        maxHFEnergy_31_50 = std::max(maxHFEnergy_31_50, E);
        minHFEnergy_31_50 = std::min(minHFEnergy_31_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,5.1])
      if ((η >=  3.1f && η <=  5.1f) ||
          (η <= -3.1f && η >= -5.1f)) {
        maxHFEnergy_31_51 = std::max(maxHFEnergy_31_51, E);
        minHFEnergy_31_51 = std::min(minHFEnergy_31_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.1,5.2])
      if ((η >=  3.1f && η <=  5.2f) ||
          (η <= -3.1f && η >= -5.2f)) {
        maxHFEnergy_31_52 = std::max(maxHFEnergy_31_52, E);
        minHFEnergy_31_52 = std::min(minHFEnergy_31_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.2,3.3])
      if ((η >=  3.2f && η <=  3.3f) ||
          (η <= -3.2f && η >= -3.3f)) {
        maxHFEnergy_32_33 = std::max(maxHFEnergy_32_33, E);
        minHFEnergy_32_33 = std::min(minHFEnergy_32_33, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,3.4])
      if ((η >=  3.2f && η <=  3.4f) ||
          (η <= -3.2f && η >= -3.4f)) {
        maxHFEnergy_32_34 = std::max(maxHFEnergy_32_34, E);
        minHFEnergy_32_34 = std::min(minHFEnergy_32_34, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,3.5])
      if ((η >=  3.2f && η <=  3.5f) ||
          (η <= -3.2f && η >= -3.5f)) {
        maxHFEnergy_32_35 = std::max(maxHFEnergy_32_35, E);
        minHFEnergy_32_35 = std::min(minHFEnergy_32_35, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,3.6])
      if ((η >=  3.2f && η <=  3.6f) ||
          (η <= -3.2f && η >= -3.6f)) {
        maxHFEnergy_32_36 = std::max(maxHFEnergy_32_36, E);
        minHFEnergy_32_36 = std::min(minHFEnergy_32_36, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,3.7])
      if ((η >=  3.2f && η <=  3.7f) ||
          (η <= -3.2f && η >= -3.7f)) {
        maxHFEnergy_32_37 = std::max(maxHFEnergy_32_37, E);
        minHFEnergy_32_37 = std::min(minHFEnergy_32_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,3.8])
      if ((η >=  3.2f && η <=  3.8f) ||
          (η <= -3.2f && η >= -3.8f)) {
        maxHFEnergy_32_38 = std::max(maxHFEnergy_32_38, E);
        minHFEnergy_32_38 = std::min(minHFEnergy_32_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,3.9])
      if ((η >=  3.2f && η <=  3.9f) ||
          (η <= -3.2f && η >= -3.9f)) {
        maxHFEnergy_32_39 = std::max(maxHFEnergy_32_39, E);
        minHFEnergy_32_39 = std::min(minHFEnergy_32_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.0])
      if ((η >=  3.2f && η <=  4.0f) ||
          (η <= -3.2f && η >= -4.0f)) {
        maxHFEnergy_32_40 = std::max(maxHFEnergy_32_40, E);
        minHFEnergy_32_40 = std::min(minHFEnergy_32_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.1])
      if ((η >=  3.2f && η <=  4.1f) ||
          (η <= -3.2f && η >= -4.1f)) {
        maxHFEnergy_32_41 = std::max(maxHFEnergy_32_41, E);
        minHFEnergy_32_41 = std::min(minHFEnergy_32_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.2])
      if ((η >=  3.2f && η <=  4.2f) ||
          (η <= -3.2f && η >= -4.2f)) {
        maxHFEnergy_32_42 = std::max(maxHFEnergy_32_42, E);
        minHFEnergy_32_42 = std::min(minHFEnergy_32_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.3])
      if ((η >=  3.2f && η <=  4.3f) ||
          (η <= -3.2f && η >= -4.3f)) {
        maxHFEnergy_32_43 = std::max(maxHFEnergy_32_43, E);
        minHFEnergy_32_43 = std::min(minHFEnergy_32_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.4])
      if ((η >=  3.2f && η <=  4.4f) ||
          (η <= -3.2f && η >= -4.4f)) {
        maxHFEnergy_32_44 = std::max(maxHFEnergy_32_44, E);
        minHFEnergy_32_44 = std::min(minHFEnergy_32_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.5])
      if ((η >=  3.2f && η <=  4.5f) ||
          (η <= -3.2f && η >= -4.5f)) {
        maxHFEnergy_32_45 = std::max(maxHFEnergy_32_45, E);
        minHFEnergy_32_45 = std::min(minHFEnergy_32_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.6])
      if ((η >=  3.2f && η <=  4.6f) ||
          (η <= -3.2f && η >= -4.6f)) {
        maxHFEnergy_32_46 = std::max(maxHFEnergy_32_46, E);
        minHFEnergy_32_46 = std::min(minHFEnergy_32_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.7])
      if ((η >=  3.2f && η <=  4.7f) ||
          (η <= -3.2f && η >= -4.7f)) {
        maxHFEnergy_32_47 = std::max(maxHFEnergy_32_47, E);
        minHFEnergy_32_47 = std::min(minHFEnergy_32_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.8])
      if ((η >=  3.2f && η <=  4.8f) ||
          (η <= -3.2f && η >= -4.8f)) {
        maxHFEnergy_32_48 = std::max(maxHFEnergy_32_48, E);
        minHFEnergy_32_48 = std::min(minHFEnergy_32_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,4.9])
      if ((η >=  3.2f && η <=  4.9f) ||
          (η <= -3.2f && η >= -4.9f)) {
        maxHFEnergy_32_49 = std::max(maxHFEnergy_32_49, E);
        minHFEnergy_32_49 = std::min(minHFEnergy_32_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,5.0])
      if ((η >=  3.2f && η <=  5.0f) ||
          (η <= -3.2f && η >= -5.0f)) {
        maxHFEnergy_32_50 = std::max(maxHFEnergy_32_50, E);
        minHFEnergy_32_50 = std::min(minHFEnergy_32_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,5.1])
      if ((η >=  3.2f && η <=  5.1f) ||
          (η <= -3.2f && η >= -5.1f)) {
        maxHFEnergy_32_51 = std::max(maxHFEnergy_32_51, E);
        minHFEnergy_32_51 = std::min(minHFEnergy_32_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.2,5.2])
      if ((η >=  3.2f && η <=  5.2f) ||
          (η <= -3.2f && η >= -5.2f)) {
        maxHFEnergy_32_52 = std::max(maxHFEnergy_32_52, E);
        minHFEnergy_32_52 = std::min(minHFEnergy_32_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.3,3.4])
      if ((η >=  3.3f && η <=  3.4f) ||
          (η <= -3.3f && η >= -3.4f)) {
        maxHFEnergy_33_34 = std::max(maxHFEnergy_33_34, E);
        minHFEnergy_33_34 = std::min(minHFEnergy_33_34, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,3.5])
      if ((η >=  3.3f && η <=  3.5f) ||
          (η <= -3.3f && η >= -3.5f)) {
        maxHFEnergy_33_35 = std::max(maxHFEnergy_33_35, E);
        minHFEnergy_33_35 = std::min(minHFEnergy_33_35, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,3.6])
      if ((η >=  3.3f && η <=  3.6f) ||
          (η <= -3.3f && η >= -3.6f)) {
        maxHFEnergy_33_36 = std::max(maxHFEnergy_33_36, E);
        minHFEnergy_33_36 = std::min(minHFEnergy_33_36, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,3.7])
      if ((η >=  3.3f && η <=  3.7f) ||
          (η <= -3.3f && η >= -3.7f)) {
        maxHFEnergy_33_37 = std::max(maxHFEnergy_33_37, E);
        minHFEnergy_33_37 = std::min(minHFEnergy_33_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,3.8])
      if ((η >=  3.3f && η <=  3.8f) ||
          (η <= -3.3f && η >= -3.8f)) {
        maxHFEnergy_33_38 = std::max(maxHFEnergy_33_38, E);
        minHFEnergy_33_38 = std::min(minHFEnergy_33_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,3.9])
      if ((η >=  3.3f && η <=  3.9f) ||
          (η <= -3.3f && η >= -3.9f)) {
        maxHFEnergy_33_39 = std::max(maxHFEnergy_33_39, E);
        minHFEnergy_33_39 = std::min(minHFEnergy_33_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.0])
      if ((η >=  3.3f && η <=  4.0f) ||
          (η <= -3.3f && η >= -4.0f)) {
        maxHFEnergy_33_40 = std::max(maxHFEnergy_33_40, E);
        minHFEnergy_33_40 = std::min(minHFEnergy_33_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.1])
      if ((η >=  3.3f && η <=  4.1f) ||
          (η <= -3.3f && η >= -4.1f)) {
        maxHFEnergy_33_41 = std::max(maxHFEnergy_33_41, E);
        minHFEnergy_33_41 = std::min(minHFEnergy_33_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.2])
      if ((η >=  3.3f && η <=  4.2f) ||
          (η <= -3.3f && η >= -4.2f)) {
        maxHFEnergy_33_42 = std::max(maxHFEnergy_33_42, E);
        minHFEnergy_33_42 = std::min(minHFEnergy_33_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.3])
      if ((η >=  3.3f && η <=  4.3f) ||
          (η <= -3.3f && η >= -4.3f)) {
        maxHFEnergy_33_43 = std::max(maxHFEnergy_33_43, E);
        minHFEnergy_33_43 = std::min(minHFEnergy_33_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.4])
      if ((η >=  3.3f && η <=  4.4f) ||
          (η <= -3.3f && η >= -4.4f)) {
        maxHFEnergy_33_44 = std::max(maxHFEnergy_33_44, E);
        minHFEnergy_33_44 = std::min(minHFEnergy_33_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.5])
      if ((η >=  3.3f && η <=  4.5f) ||
          (η <= -3.3f && η >= -4.5f)) {
        maxHFEnergy_33_45 = std::max(maxHFEnergy_33_45, E);
        minHFEnergy_33_45 = std::min(minHFEnergy_33_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.6])
      if ((η >=  3.3f && η <=  4.6f) ||
          (η <= -3.3f && η >= -4.6f)) {
        maxHFEnergy_33_46 = std::max(maxHFEnergy_33_46, E);
        minHFEnergy_33_46 = std::min(minHFEnergy_33_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.7])
      if ((η >=  3.3f && η <=  4.7f) ||
          (η <= -3.3f && η >= -4.7f)) {
        maxHFEnergy_33_47 = std::max(maxHFEnergy_33_47, E);
        minHFEnergy_33_47 = std::min(minHFEnergy_33_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.8])
      if ((η >=  3.3f && η <=  4.8f) ||
          (η <= -3.3f && η >= -4.8f)) {
        maxHFEnergy_33_48 = std::max(maxHFEnergy_33_48, E);
        minHFEnergy_33_48 = std::min(minHFEnergy_33_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,4.9])
      if ((η >=  3.3f && η <=  4.9f) ||
          (η <= -3.3f && η >= -4.9f)) {
        maxHFEnergy_33_49 = std::max(maxHFEnergy_33_49, E);
        minHFEnergy_33_49 = std::min(minHFEnergy_33_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,5.0])
      if ((η >=  3.3f && η <=  5.0f) ||
          (η <= -3.3f && η >= -5.0f)) {
        maxHFEnergy_33_50 = std::max(maxHFEnergy_33_50, E);
        minHFEnergy_33_50 = std::min(minHFEnergy_33_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,5.1])
      if ((η >=  3.3f && η <=  5.1f) ||
          (η <= -3.3f && η >= -5.1f)) {
        maxHFEnergy_33_51 = std::max(maxHFEnergy_33_51, E);
        minHFEnergy_33_51 = std::min(minHFEnergy_33_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.3,5.2])
      if ((η >=  3.3f && η <=  5.2f) ||
          (η <= -3.3f && η >= -5.2f)) {
        maxHFEnergy_33_52 = std::max(maxHFEnergy_33_52, E);
        minHFEnergy_33_52 = std::min(minHFEnergy_33_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.4,3.5])
      if ((η >=  3.4f && η <=  3.5f) ||
          (η <= -3.4f && η >= -3.5f)) {
        maxHFEnergy_34_35 = std::max(maxHFEnergy_34_35, E);
        minHFEnergy_34_35 = std::min(minHFEnergy_34_35, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,3.6])
      if ((η >=  3.4f && η <=  3.6f) ||
          (η <= -3.4f && η >= -3.6f)) {
        maxHFEnergy_34_36 = std::max(maxHFEnergy_34_36, E);
        minHFEnergy_34_36 = std::min(minHFEnergy_34_36, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,3.7])
      if ((η >=  3.4f && η <=  3.7f) ||
          (η <= -3.4f && η >= -3.7f)) {
        maxHFEnergy_34_37 = std::max(maxHFEnergy_34_37, E);
        minHFEnergy_34_37 = std::min(minHFEnergy_34_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,3.8])
      if ((η >=  3.4f && η <=  3.8f) ||
          (η <= -3.4f && η >= -3.8f)) {
        maxHFEnergy_34_38 = std::max(maxHFEnergy_34_38, E);
        minHFEnergy_34_38 = std::min(minHFEnergy_34_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,3.9])
      if ((η >=  3.4f && η <=  3.9f) ||
          (η <= -3.4f && η >= -3.9f)) {
        maxHFEnergy_34_39 = std::max(maxHFEnergy_34_39, E);
        minHFEnergy_34_39 = std::min(minHFEnergy_34_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.0])
      if ((η >=  3.4f && η <=  4.0f) ||
          (η <= -3.4f && η >= -4.0f)) {
        maxHFEnergy_34_40 = std::max(maxHFEnergy_34_40, E);
        minHFEnergy_34_40 = std::min(minHFEnergy_34_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.1])
      if ((η >=  3.4f && η <=  4.1f) ||
          (η <= -3.4f && η >= -4.1f)) {
        maxHFEnergy_34_41 = std::max(maxHFEnergy_34_41, E);
        minHFEnergy_34_41 = std::min(minHFEnergy_34_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.2])
      if ((η >=  3.4f && η <=  4.2f) ||
          (η <= -3.4f && η >= -4.2f)) {
        maxHFEnergy_34_42 = std::max(maxHFEnergy_34_42, E);
        minHFEnergy_34_42 = std::min(minHFEnergy_34_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.3])
      if ((η >=  3.4f && η <=  4.3f) ||
          (η <= -3.4f && η >= -4.3f)) {
        maxHFEnergy_34_43 = std::max(maxHFEnergy_34_43, E);
        minHFEnergy_34_43 = std::min(minHFEnergy_34_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.4])
      if ((η >=  3.4f && η <=  4.4f) ||
          (η <= -3.4f && η >= -4.4f)) {
        maxHFEnergy_34_44 = std::max(maxHFEnergy_34_44, E);
        minHFEnergy_34_44 = std::min(minHFEnergy_34_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.5])
      if ((η >=  3.4f && η <=  4.5f) ||
          (η <= -3.4f && η >= -4.5f)) {
        maxHFEnergy_34_45 = std::max(maxHFEnergy_34_45, E);
        minHFEnergy_34_45 = std::min(minHFEnergy_34_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.6])
      if ((η >=  3.4f && η <=  4.6f) ||
          (η <= -3.4f && η >= -4.6f)) {
        maxHFEnergy_34_46 = std::max(maxHFEnergy_34_46, E);
        minHFEnergy_34_46 = std::min(minHFEnergy_34_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.7])
      if ((η >=  3.4f && η <=  4.7f) ||
          (η <= -3.4f && η >= -4.7f)) {
        maxHFEnergy_34_47 = std::max(maxHFEnergy_34_47, E);
        minHFEnergy_34_47 = std::min(minHFEnergy_34_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.8])
      if ((η >=  3.4f && η <=  4.8f) ||
          (η <= -3.4f && η >= -4.8f)) {
        maxHFEnergy_34_48 = std::max(maxHFEnergy_34_48, E);
        minHFEnergy_34_48 = std::min(minHFEnergy_34_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,4.9])
      if ((η >=  3.4f && η <=  4.9f) ||
          (η <= -3.4f && η >= -4.9f)) {
        maxHFEnergy_34_49 = std::max(maxHFEnergy_34_49, E);
        minHFEnergy_34_49 = std::min(minHFEnergy_34_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,5.0])
      if ((η >=  3.4f && η <=  5.0f) ||
          (η <= -3.4f && η >= -5.0f)) {
        maxHFEnergy_34_50 = std::max(maxHFEnergy_34_50, E);
        minHFEnergy_34_50 = std::min(minHFEnergy_34_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,5.1])
      if ((η >=  3.4f && η <=  5.1f) ||
          (η <= -3.4f && η >= -5.1f)) {
        maxHFEnergy_34_51 = std::max(maxHFEnergy_34_51, E);
        minHFEnergy_34_51 = std::min(minHFEnergy_34_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.4,5.2])
      if ((η >=  3.4f && η <=  5.2f) ||
          (η <= -3.4f && η >= -5.2f)) {
        maxHFEnergy_34_52 = std::max(maxHFEnergy_34_52, E);
        minHFEnergy_34_52 = std::min(minHFEnergy_34_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.5,3.6])
      if ((η >=  3.5f && η <=  3.6f) ||
          (η <= -3.5f && η >= -3.6f)) {
        maxHFEnergy_35_36 = std::max(maxHFEnergy_35_36, E);
        minHFEnergy_35_36 = std::min(minHFEnergy_35_36, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,3.7])
      if ((η >=  3.5f && η <=  3.7f) ||
          (η <= -3.5f && η >= -3.7f)) {
        maxHFEnergy_35_37 = std::max(maxHFEnergy_35_37, E);
        minHFEnergy_35_37 = std::min(minHFEnergy_35_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,3.8])
      if ((η >=  3.5f && η <=  3.8f) ||
          (η <= -3.5f && η >= -3.8f)) {
        maxHFEnergy_35_38 = std::max(maxHFEnergy_35_38, E);
        minHFEnergy_35_38 = std::min(minHFEnergy_35_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,3.9])
      if ((η >=  3.5f && η <=  3.9f) ||
          (η <= -3.5f && η >= -3.9f)) {
        maxHFEnergy_35_39 = std::max(maxHFEnergy_35_39, E);
        minHFEnergy_35_39 = std::min(minHFEnergy_35_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.0])
      if ((η >=  3.5f && η <=  4.0f) ||
          (η <= -3.5f && η >= -4.0f)) {
        maxHFEnergy_35_40 = std::max(maxHFEnergy_35_40, E);
        minHFEnergy_35_40 = std::min(minHFEnergy_35_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.1])
      if ((η >=  3.5f && η <=  4.1f) ||
          (η <= -3.5f && η >= -4.1f)) {
        maxHFEnergy_35_41 = std::max(maxHFEnergy_35_41, E);
        minHFEnergy_35_41 = std::min(minHFEnergy_35_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.2])
      if ((η >=  3.5f && η <=  4.2f) ||
          (η <= -3.5f && η >= -4.2f)) {
        maxHFEnergy_35_42 = std::max(maxHFEnergy_35_42, E);
        minHFEnergy_35_42 = std::min(minHFEnergy_35_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.3])
      if ((η >=  3.5f && η <=  4.3f) ||
          (η <= -3.5f && η >= -4.3f)) {
        maxHFEnergy_35_43 = std::max(maxHFEnergy_35_43, E);
        minHFEnergy_35_43 = std::min(minHFEnergy_35_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.4])
      if ((η >=  3.5f && η <=  4.4f) ||
          (η <= -3.5f && η >= -4.4f)) {
        maxHFEnergy_35_44 = std::max(maxHFEnergy_35_44, E);
        minHFEnergy_35_44 = std::min(minHFEnergy_35_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.5])
      if ((η >=  3.5f && η <=  4.5f) ||
          (η <= -3.5f && η >= -4.5f)) {
        maxHFEnergy_35_45 = std::max(maxHFEnergy_35_45, E);
        minHFEnergy_35_45 = std::min(minHFEnergy_35_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.6])
      if ((η >=  3.5f && η <=  4.6f) ||
          (η <= -3.5f && η >= -4.6f)) {
        maxHFEnergy_35_46 = std::max(maxHFEnergy_35_46, E);
        minHFEnergy_35_46 = std::min(minHFEnergy_35_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.7])
      if ((η >=  3.5f && η <=  4.7f) ||
          (η <= -3.5f && η >= -4.7f)) {
        maxHFEnergy_35_47 = std::max(maxHFEnergy_35_47, E);
        minHFEnergy_35_47 = std::min(minHFEnergy_35_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.8])
      if ((η >=  3.5f && η <=  4.8f) ||
          (η <= -3.5f && η >= -4.8f)) {
        maxHFEnergy_35_48 = std::max(maxHFEnergy_35_48, E);
        minHFEnergy_35_48 = std::min(minHFEnergy_35_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,4.9])
      if ((η >=  3.5f && η <=  4.9f) ||
          (η <= -3.5f && η >= -4.9f)) {
        maxHFEnergy_35_49 = std::max(maxHFEnergy_35_49, E);
        minHFEnergy_35_49 = std::min(minHFEnergy_35_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,5.0])
      if ((η >=  3.5f && η <=  5.0f) ||
          (η <= -3.5f && η >= -5.0f)) {
        maxHFEnergy_35_50 = std::max(maxHFEnergy_35_50, E);
        minHFEnergy_35_50 = std::min(minHFEnergy_35_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,5.1])
      if ((η >=  3.5f && η <=  5.1f) ||
          (η <= -3.5f && η >= -5.1f)) {
        maxHFEnergy_35_51 = std::max(maxHFEnergy_35_51, E);
        minHFEnergy_35_51 = std::min(minHFEnergy_35_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.5,5.2])
      if ((η >=  3.5f && η <=  5.2f) ||
          (η <= -3.5f && η >= -5.2f)) {
        maxHFEnergy_35_52 = std::max(maxHFEnergy_35_52, E);
        minHFEnergy_35_52 = std::min(minHFEnergy_35_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.6,3.7])
      if ((η >=  3.6f && η <=  3.7f) ||
          (η <= -3.6f && η >= -3.7f)) {
        maxHFEnergy_36_37 = std::max(maxHFEnergy_36_37, E);
        minHFEnergy_36_37 = std::min(minHFEnergy_36_37, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,3.8])
      if ((η >=  3.6f && η <=  3.8f) ||
          (η <= -3.6f && η >= -3.8f)) {
        maxHFEnergy_36_38 = std::max(maxHFEnergy_36_38, E);
        minHFEnergy_36_38 = std::min(minHFEnergy_36_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,3.9])
      if ((η >=  3.6f && η <=  3.9f) ||
          (η <= -3.6f && η >= -3.9f)) {
        maxHFEnergy_36_39 = std::max(maxHFEnergy_36_39, E);
        minHFEnergy_36_39 = std::min(minHFEnergy_36_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.0])
      if ((η >=  3.6f && η <=  4.0f) ||
          (η <= -3.6f && η >= -4.0f)) {
        maxHFEnergy_36_40 = std::max(maxHFEnergy_36_40, E);
        minHFEnergy_36_40 = std::min(minHFEnergy_36_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.1])
      if ((η >=  3.6f && η <=  4.1f) ||
          (η <= -3.6f && η >= -4.1f)) {
        maxHFEnergy_36_41 = std::max(maxHFEnergy_36_41, E);
        minHFEnergy_36_41 = std::min(minHFEnergy_36_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.2])
      if ((η >=  3.6f && η <=  4.2f) ||
          (η <= -3.6f && η >= -4.2f)) {
        maxHFEnergy_36_42 = std::max(maxHFEnergy_36_42, E);
        minHFEnergy_36_42 = std::min(minHFEnergy_36_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.3])
      if ((η >=  3.6f && η <=  4.3f) ||
          (η <= -3.6f && η >= -4.3f)) {
        maxHFEnergy_36_43 = std::max(maxHFEnergy_36_43, E);
        minHFEnergy_36_43 = std::min(minHFEnergy_36_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.4])
      if ((η >=  3.6f && η <=  4.4f) ||
          (η <= -3.6f && η >= -4.4f)) {
        maxHFEnergy_36_44 = std::max(maxHFEnergy_36_44, E);
        minHFEnergy_36_44 = std::min(minHFEnergy_36_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.5])
      if ((η >=  3.6f && η <=  4.5f) ||
          (η <= -3.6f && η >= -4.5f)) {
        maxHFEnergy_36_45 = std::max(maxHFEnergy_36_45, E);
        minHFEnergy_36_45 = std::min(minHFEnergy_36_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.6])
      if ((η >=  3.6f && η <=  4.6f) ||
          (η <= -3.6f && η >= -4.6f)) {
        maxHFEnergy_36_46 = std::max(maxHFEnergy_36_46, E);
        minHFEnergy_36_46 = std::min(minHFEnergy_36_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.7])
      if ((η >=  3.6f && η <=  4.7f) ||
          (η <= -3.6f && η >= -4.7f)) {
        maxHFEnergy_36_47 = std::max(maxHFEnergy_36_47, E);
        minHFEnergy_36_47 = std::min(minHFEnergy_36_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.8])
      if ((η >=  3.6f && η <=  4.8f) ||
          (η <= -3.6f && η >= -4.8f)) {
        maxHFEnergy_36_48 = std::max(maxHFEnergy_36_48, E);
        minHFEnergy_36_48 = std::min(minHFEnergy_36_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,4.9])
      if ((η >=  3.6f && η <=  4.9f) ||
          (η <= -3.6f && η >= -4.9f)) {
        maxHFEnergy_36_49 = std::max(maxHFEnergy_36_49, E);
        minHFEnergy_36_49 = std::min(minHFEnergy_36_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,5.0])
      if ((η >=  3.6f && η <=  5.0f) ||
          (η <= -3.6f && η >= -5.0f)) {
        maxHFEnergy_36_50 = std::max(maxHFEnergy_36_50, E);
        minHFEnergy_36_50 = std::min(minHFEnergy_36_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,5.1])
      if ((η >=  3.6f && η <=  5.1f) ||
          (η <= -3.6f && η >= -5.1f)) {
        maxHFEnergy_36_51 = std::max(maxHFEnergy_36_51, E);
        minHFEnergy_36_51 = std::min(minHFEnergy_36_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.6,5.2])
      if ((η >=  3.6f && η <=  5.2f) ||
          (η <= -3.6f && η >= -5.2f)) {
        maxHFEnergy_36_52 = std::max(maxHFEnergy_36_52, E);
        minHFEnergy_36_52 = std::min(minHFEnergy_36_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.7,3.8])
      if ((η >=  3.7f && η <=  3.8f) ||
          (η <= -3.7f && η >= -3.8f)) {
        maxHFEnergy_37_38 = std::max(maxHFEnergy_37_38, E);
        minHFEnergy_37_38 = std::min(minHFEnergy_37_38, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,3.9])
      if ((η >=  3.7f && η <=  3.9f) ||
          (η <= -3.7f && η >= -3.9f)) {
        maxHFEnergy_37_39 = std::max(maxHFEnergy_37_39, E);
        minHFEnergy_37_39 = std::min(minHFEnergy_37_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.0])
      if ((η >=  3.7f && η <=  4.0f) ||
          (η <= -3.7f && η >= -4.0f)) {
        maxHFEnergy_37_40 = std::max(maxHFEnergy_37_40, E);
        minHFEnergy_37_40 = std::min(minHFEnergy_37_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.1])
      if ((η >=  3.7f && η <=  4.1f) ||
          (η <= -3.7f && η >= -4.1f)) {
        maxHFEnergy_37_41 = std::max(maxHFEnergy_37_41, E);
        minHFEnergy_37_41 = std::min(minHFEnergy_37_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.2])
      if ((η >=  3.7f && η <=  4.2f) ||
          (η <= -3.7f && η >= -4.2f)) {
        maxHFEnergy_37_42 = std::max(maxHFEnergy_37_42, E);
        minHFEnergy_37_42 = std::min(minHFEnergy_37_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.3])
      if ((η >=  3.7f && η <=  4.3f) ||
          (η <= -3.7f && η >= -4.3f)) {
        maxHFEnergy_37_43 = std::max(maxHFEnergy_37_43, E);
        minHFEnergy_37_43 = std::min(minHFEnergy_37_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.4])
      if ((η >=  3.7f && η <=  4.4f) ||
          (η <= -3.7f && η >= -4.4f)) {
        maxHFEnergy_37_44 = std::max(maxHFEnergy_37_44, E);
        minHFEnergy_37_44 = std::min(minHFEnergy_37_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.5])
      if ((η >=  3.7f && η <=  4.5f) ||
          (η <= -3.7f && η >= -4.5f)) {
        maxHFEnergy_37_45 = std::max(maxHFEnergy_37_45, E);
        minHFEnergy_37_45 = std::min(minHFEnergy_37_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.6])
      if ((η >=  3.7f && η <=  4.6f) ||
          (η <= -3.7f && η >= -4.6f)) {
        maxHFEnergy_37_46 = std::max(maxHFEnergy_37_46, E);
        minHFEnergy_37_46 = std::min(minHFEnergy_37_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.7])
      if ((η >=  3.7f && η <=  4.7f) ||
          (η <= -3.7f && η >= -4.7f)) {
        maxHFEnergy_37_47 = std::max(maxHFEnergy_37_47, E);
        minHFEnergy_37_47 = std::min(minHFEnergy_37_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.8])
      if ((η >=  3.7f && η <=  4.8f) ||
          (η <= -3.7f && η >= -4.8f)) {
        maxHFEnergy_37_48 = std::max(maxHFEnergy_37_48, E);
        minHFEnergy_37_48 = std::min(minHFEnergy_37_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,4.9])
      if ((η >=  3.7f && η <=  4.9f) ||
          (η <= -3.7f && η >= -4.9f)) {
        maxHFEnergy_37_49 = std::max(maxHFEnergy_37_49, E);
        minHFEnergy_37_49 = std::min(minHFEnergy_37_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,5.0])
      if ((η >=  3.7f && η <=  5.0f) ||
          (η <= -3.7f && η >= -5.0f)) {
        maxHFEnergy_37_50 = std::max(maxHFEnergy_37_50, E);
        minHFEnergy_37_50 = std::min(minHFEnergy_37_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,5.1])
      if ((η >=  3.7f && η <=  5.1f) ||
          (η <= -3.7f && η >= -5.1f)) {
        maxHFEnergy_37_51 = std::max(maxHFEnergy_37_51, E);
        minHFEnergy_37_51 = std::min(minHFEnergy_37_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.7,5.2])
      if ((η >=  3.7f && η <=  5.2f) ||
          (η <= -3.7f && η >= -5.2f)) {
        maxHFEnergy_37_52 = std::max(maxHFEnergy_37_52, E);
        minHFEnergy_37_52 = std::min(minHFEnergy_37_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.8,3.9])
      if ((η >=  3.8f && η <=  3.9f) ||
          (η <= -3.8f && η >= -3.9f)) {
        maxHFEnergy_38_39 = std::max(maxHFEnergy_38_39, E);
        minHFEnergy_38_39 = std::min(minHFEnergy_38_39, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.0])
      if ((η >=  3.8f && η <=  4.0f) ||
          (η <= -3.8f && η >= -4.0f)) {
        maxHFEnergy_38_40 = std::max(maxHFEnergy_38_40, E);
        minHFEnergy_38_40 = std::min(minHFEnergy_38_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.1])
      if ((η >=  3.8f && η <=  4.1f) ||
          (η <= -3.8f && η >= -4.1f)) {
        maxHFEnergy_38_41 = std::max(maxHFEnergy_38_41, E);
        minHFEnergy_38_41 = std::min(minHFEnergy_38_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.2])
      if ((η >=  3.8f && η <=  4.2f) ||
          (η <= -3.8f && η >= -4.2f)) {
        maxHFEnergy_38_42 = std::max(maxHFEnergy_38_42, E);
        minHFEnergy_38_42 = std::min(minHFEnergy_38_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.3])
      if ((η >=  3.8f && η <=  4.3f) ||
          (η <= -3.8f && η >= -4.3f)) {
        maxHFEnergy_38_43 = std::max(maxHFEnergy_38_43, E);
        minHFEnergy_38_43 = std::min(minHFEnergy_38_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.4])
      if ((η >=  3.8f && η <=  4.4f) ||
          (η <= -3.8f && η >= -4.4f)) {
        maxHFEnergy_38_44 = std::max(maxHFEnergy_38_44, E);
        minHFEnergy_38_44 = std::min(minHFEnergy_38_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.5])
      if ((η >=  3.8f && η <=  4.5f) ||
          (η <= -3.8f && η >= -4.5f)) {
        maxHFEnergy_38_45 = std::max(maxHFEnergy_38_45, E);
        minHFEnergy_38_45 = std::min(minHFEnergy_38_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.6])
      if ((η >=  3.8f && η <=  4.6f) ||
          (η <= -3.8f && η >= -4.6f)) {
        maxHFEnergy_38_46 = std::max(maxHFEnergy_38_46, E);
        minHFEnergy_38_46 = std::min(minHFEnergy_38_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.7])
      if ((η >=  3.8f && η <=  4.7f) ||
          (η <= -3.8f && η >= -4.7f)) {
        maxHFEnergy_38_47 = std::max(maxHFEnergy_38_47, E);
        minHFEnergy_38_47 = std::min(minHFEnergy_38_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.8])
      if ((η >=  3.8f && η <=  4.8f) ||
          (η <= -3.8f && η >= -4.8f)) {
        maxHFEnergy_38_48 = std::max(maxHFEnergy_38_48, E);
        minHFEnergy_38_48 = std::min(minHFEnergy_38_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,4.9])
      if ((η >=  3.8f && η <=  4.9f) ||
          (η <= -3.8f && η >= -4.9f)) {
        maxHFEnergy_38_49 = std::max(maxHFEnergy_38_49, E);
        minHFEnergy_38_49 = std::min(minHFEnergy_38_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,5.0])
      if ((η >=  3.8f && η <=  5.0f) ||
          (η <= -3.8f && η >= -5.0f)) {
        maxHFEnergy_38_50 = std::max(maxHFEnergy_38_50, E);
        minHFEnergy_38_50 = std::min(minHFEnergy_38_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,5.1])
      if ((η >=  3.8f && η <=  5.1f) ||
          (η <= -3.8f && η >= -5.1f)) {
        maxHFEnergy_38_51 = std::max(maxHFEnergy_38_51, E);
        minHFEnergy_38_51 = std::min(minHFEnergy_38_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.8,5.2])
      if ((η >=  3.8f && η <=  5.2f) ||
          (η <= -3.8f && η >= -5.2f)) {
        maxHFEnergy_38_52 = std::max(maxHFEnergy_38_52, E);
        minHFEnergy_38_52 = std::min(minHFEnergy_38_52, E);
        }


      // global HF extrema (any tower in |η|∈[3.9,4.0])
      if ((η >=  3.9f && η <=  4.0f) ||
          (η <= -3.9f && η >= -4.0f)) {
        maxHFEnergy_39_40 = std::max(maxHFEnergy_39_40, E);
        minHFEnergy_39_40 = std::min(minHFEnergy_39_40, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.1])
      if ((η >=  3.9f && η <=  4.1f) ||
          (η <= -3.9f && η >= -4.1f)) {
        maxHFEnergy_39_41 = std::max(maxHFEnergy_39_41, E);
        minHFEnergy_39_41 = std::min(minHFEnergy_39_41, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.2])
      if ((η >=  3.9f && η <=  4.2f) ||
          (η <= -3.9f && η >= -4.2f)) {
        maxHFEnergy_39_42 = std::max(maxHFEnergy_39_42, E);
        minHFEnergy_39_42 = std::min(minHFEnergy_39_42, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.3])
      if ((η >=  3.9f && η <=  4.3f) ||
          (η <= -3.9f && η >= -4.3f)) {
        maxHFEnergy_39_43 = std::max(maxHFEnergy_39_43, E);
        minHFEnergy_39_43 = std::min(minHFEnergy_39_43, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.4])
      if ((η >=  3.9f && η <=  4.4f) ||
          (η <= -3.9f && η >= -4.4f)) {
        maxHFEnergy_39_44 = std::max(maxHFEnergy_39_44, E);
        minHFEnergy_39_44 = std::min(minHFEnergy_39_44, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.5])
      if ((η >=  3.9f && η <=  4.5f) ||
          (η <= -3.9f && η >= -4.5f)) {
        maxHFEnergy_39_45 = std::max(maxHFEnergy_39_45, E);
        minHFEnergy_39_45 = std::min(minHFEnergy_39_45, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.6])
      if ((η >=  3.9f && η <=  4.6f) ||
          (η <= -3.9f && η >= -4.6f)) {
        maxHFEnergy_39_46 = std::max(maxHFEnergy_39_46, E);
        minHFEnergy_39_46 = std::min(minHFEnergy_39_46, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.7])
      if ((η >=  3.9f && η <=  4.7f) ||
          (η <= -3.9f && η >= -4.7f)) {
        maxHFEnergy_39_47 = std::max(maxHFEnergy_39_47, E);
        minHFEnergy_39_47 = std::min(minHFEnergy_39_47, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.8])
      if ((η >=  3.9f && η <=  4.8f) ||
          (η <= -3.9f && η >= -4.8f)) {
        maxHFEnergy_39_48 = std::max(maxHFEnergy_39_48, E);
        minHFEnergy_39_48 = std::min(minHFEnergy_39_48, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,4.9])
      if ((η >=  3.9f && η <=  4.9f) ||
          (η <= -3.9f && η >= -4.9f)) {
        maxHFEnergy_39_49 = std::max(maxHFEnergy_39_49, E);
        minHFEnergy_39_49 = std::min(minHFEnergy_39_49, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,5.0])
      if ((η >=  3.9f && η <=  5.0f) ||
          (η <= -3.9f && η >= -5.0f)) {
        maxHFEnergy_39_50 = std::max(maxHFEnergy_39_50, E);
        minHFEnergy_39_50 = std::min(minHFEnergy_39_50, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,5.1])
      if ((η >=  3.9f && η <=  5.1f) ||
          (η <= -3.9f && η >= -5.1f)) {
        maxHFEnergy_39_51 = std::max(maxHFEnergy_39_51, E);
        minHFEnergy_39_51 = std::min(minHFEnergy_39_51, E);
        }
      // global HF extrema (any tower in |η|∈[3.9,5.2])
      if ((η >=  3.9f && η <=  5.2f) ||
          (η <= -3.9f && η >= -5.2f)) {
        maxHFEnergy_39_52 = std::max(maxHFEnergy_39_52, E);
        minHFEnergy_39_52 = std::min(minHFEnergy_39_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.0,4.1])
      if ((η >=  4.0f && η <=  4.1f) ||
          (η <= -4.0f && η >= -4.1f)) {
        maxHFEnergy_40_41 = std::max(maxHFEnergy_40_41, E);
        minHFEnergy_40_41 = std::min(minHFEnergy_40_41, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.2])
      if ((η >=  4.0f && η <=  4.2f) ||
          (η <= -4.0f && η >= -4.2f)) {
        maxHFEnergy_40_42 = std::max(maxHFEnergy_40_42, E);
        minHFEnergy_40_42 = std::min(minHFEnergy_40_42, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.3])
      if ((η >=  4.0f && η <=  4.3f) ||
          (η <= -4.0f && η >= -4.3f)) {
        maxHFEnergy_40_43 = std::max(maxHFEnergy_40_43, E);
        minHFEnergy_40_43 = std::min(minHFEnergy_40_43, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.4])
      if ((η >=  4.0f && η <=  4.4f) ||
          (η <= -4.0f && η >= -4.4f)) {
        maxHFEnergy_40_44 = std::max(maxHFEnergy_40_44, E);
        minHFEnergy_40_44 = std::min(minHFEnergy_40_44, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.5])
      if ((η >=  4.0f && η <=  4.5f) ||
          (η <= -4.0f && η >= -4.5f)) {
        maxHFEnergy_40_45 = std::max(maxHFEnergy_40_45, E);
        minHFEnergy_40_45 = std::min(minHFEnergy_40_45, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.6])
      if ((η >=  4.0f && η <=  4.6f) ||
          (η <= -4.0f && η >= -4.6f)) {
        maxHFEnergy_40_46 = std::max(maxHFEnergy_40_46, E);
        minHFEnergy_40_46 = std::min(minHFEnergy_40_46, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.7])
      if ((η >=  4.0f && η <=  4.7f) ||
          (η <= -4.0f && η >= -4.7f)) {
        maxHFEnergy_40_47 = std::max(maxHFEnergy_40_47, E);
        minHFEnergy_40_47 = std::min(minHFEnergy_40_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.8])
      if ((η >=  4.0f && η <=  4.8f) ||
          (η <= -4.0f && η >= -4.8f)) {
        maxHFEnergy_40_48 = std::max(maxHFEnergy_40_48, E);
        minHFEnergy_40_48 = std::min(minHFEnergy_40_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,4.9])
      if ((η >=  4.0f && η <=  4.9f) ||
          (η <= -4.0f && η >= -4.9f)) {
        maxHFEnergy_40_49 = std::max(maxHFEnergy_40_49, E);
        minHFEnergy_40_49 = std::min(minHFEnergy_40_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,5.0])
      if ((η >=  4.0f && η <=  5.0f) ||
          (η <= -4.0f && η >= -5.0f)) {
        maxHFEnergy_40_50 = std::max(maxHFEnergy_40_50, E);
        minHFEnergy_40_50 = std::min(minHFEnergy_40_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,5.1])
      if ((η >=  4.0f && η <=  5.1f) ||
          (η <= -4.0f && η >= -5.1f)) {
        maxHFEnergy_40_51 = std::max(maxHFEnergy_40_51, E);
        minHFEnergy_40_51 = std::min(minHFEnergy_40_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.0,5.2])
      if ((η >=  4.0f && η <=  5.2f) ||
          (η <= -4.0f && η >= -5.2f)) {
        maxHFEnergy_40_52 = std::max(maxHFEnergy_40_52, E);
        minHFEnergy_40_52 = std::min(minHFEnergy_40_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.1,4.2])
      if ((η >=  4.1f && η <=  4.2f) ||
          (η <= -4.1f && η >= -4.2f)) {
        maxHFEnergy_41_42 = std::max(maxHFEnergy_41_42, E);
        minHFEnergy_41_42 = std::min(minHFEnergy_41_42, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.3])
      if ((η >=  4.1f && η <=  4.3f) ||
          (η <= -4.1f && η >= -4.3f)) {
        maxHFEnergy_41_43 = std::max(maxHFEnergy_41_43, E);
        minHFEnergy_41_43 = std::min(minHFEnergy_41_43, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.4])
      if ((η >=  4.1f && η <=  4.4f) ||
          (η <= -4.1f && η >= -4.4f)) {
        maxHFEnergy_41_44 = std::max(maxHFEnergy_41_44, E);
        minHFEnergy_41_44 = std::min(minHFEnergy_41_44, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.5])
      if ((η >=  4.1f && η <=  4.5f) ||
          (η <= -4.1f && η >= -4.5f)) {
        maxHFEnergy_41_45 = std::max(maxHFEnergy_41_45, E);
        minHFEnergy_41_45 = std::min(minHFEnergy_41_45, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.6])
      if ((η >=  4.1f && η <=  4.6f) ||
          (η <= -4.1f && η >= -4.6f)) {
        maxHFEnergy_41_46 = std::max(maxHFEnergy_41_46, E);
        minHFEnergy_41_46 = std::min(minHFEnergy_41_46, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.7])
      if ((η >=  4.1f && η <=  4.7f) ||
          (η <= -4.1f && η >= -4.7f)) {
        maxHFEnergy_41_47 = std::max(maxHFEnergy_41_47, E);
        minHFEnergy_41_47 = std::min(minHFEnergy_41_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.8])
      if ((η >=  4.1f && η <=  4.8f) ||
          (η <= -4.1f && η >= -4.8f)) {
        maxHFEnergy_41_48 = std::max(maxHFEnergy_41_48, E);
        minHFEnergy_41_48 = std::min(minHFEnergy_41_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,4.9])
      if ((η >=  4.1f && η <=  4.9f) ||
          (η <= -4.1f && η >= -4.9f)) {
        maxHFEnergy_41_49 = std::max(maxHFEnergy_41_49, E);
        minHFEnergy_41_49 = std::min(minHFEnergy_41_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,5.0])
      if ((η >=  4.1f && η <=  5.0f) ||
          (η <= -4.1f && η >= -5.0f)) {
        maxHFEnergy_41_50 = std::max(maxHFEnergy_41_50, E);
        minHFEnergy_41_50 = std::min(minHFEnergy_41_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,5.1])
      if ((η >=  4.1f && η <=  5.1f) ||
          (η <= -4.1f && η >= -5.1f)) {
        maxHFEnergy_41_51 = std::max(maxHFEnergy_41_51, E);
        minHFEnergy_41_51 = std::min(minHFEnergy_41_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.1,5.2])
      if ((η >=  4.1f && η <=  5.2f) ||
          (η <= -4.1f && η >= -5.2f)) {
        maxHFEnergy_41_52 = std::max(maxHFEnergy_41_52, E);
        minHFEnergy_41_52 = std::min(minHFEnergy_41_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.2,4.3])
      if ((η >=  4.2f && η <=  4.3f) ||
          (η <= -4.2f && η >= -4.3f)) {
        maxHFEnergy_42_43 = std::max(maxHFEnergy_42_43, E);
        minHFEnergy_42_43 = std::min(minHFEnergy_42_43, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,4.4])
      if ((η >=  4.2f && η <=  4.4f) ||
          (η <= -4.2f && η >= -4.4f)) {
        maxHFEnergy_42_44 = std::max(maxHFEnergy_42_44, E);
        minHFEnergy_42_44 = std::min(minHFEnergy_42_44, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,4.5])
      if ((η >=  4.2f && η <=  4.5f) ||
          (η <= -4.2f && η >= -4.5f)) {
        maxHFEnergy_42_45 = std::max(maxHFEnergy_42_45, E);
        minHFEnergy_42_45 = std::min(minHFEnergy_42_45, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,4.6])
      if ((η >=  4.2f && η <=  4.6f) ||
          (η <= -4.2f && η >= -4.6f)) {
        maxHFEnergy_42_46 = std::max(maxHFEnergy_42_46, E);
        minHFEnergy_42_46 = std::min(minHFEnergy_42_46, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,4.7])
      if ((η >=  4.2f && η <=  4.7f) ||
          (η <= -4.2f && η >= -4.7f)) {
        maxHFEnergy_42_47 = std::max(maxHFEnergy_42_47, E);
        minHFEnergy_42_47 = std::min(minHFEnergy_42_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,4.8])
      if ((η >=  4.2f && η <=  4.8f) ||
          (η <= -4.2f && η >= -4.8f)) {
        maxHFEnergy_42_48 = std::max(maxHFEnergy_42_48, E);
        minHFEnergy_42_48 = std::min(minHFEnergy_42_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,4.9])
      if ((η >=  4.2f && η <=  4.9f) ||
          (η <= -4.2f && η >= -4.9f)) {
        maxHFEnergy_42_49 = std::max(maxHFEnergy_42_49, E);
        minHFEnergy_42_49 = std::min(minHFEnergy_42_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,5.0])
      if ((η >=  4.2f && η <=  5.0f) ||
          (η <= -4.2f && η >= -5.0f)) {
        maxHFEnergy_42_50 = std::max(maxHFEnergy_42_50, E);
        minHFEnergy_42_50 = std::min(minHFEnergy_42_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,5.1])
      if ((η >=  4.2f && η <=  5.1f) ||
          (η <= -4.2f && η >= -5.1f)) {
        maxHFEnergy_42_51 = std::max(maxHFEnergy_42_51, E);
        minHFEnergy_42_51 = std::min(minHFEnergy_42_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.2,5.2])
      if ((η >=  4.2f && η <=  5.2f) ||
          (η <= -4.2f && η >= -5.2f)) {
        maxHFEnergy_42_52 = std::max(maxHFEnergy_42_52, E);
        minHFEnergy_42_52 = std::min(minHFEnergy_42_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.3,4.4])
      if ((η >=  4.3f && η <=  4.4f) ||
          (η <= -4.3f && η >= -4.4f)) {
        maxHFEnergy_43_44 = std::max(maxHFEnergy_43_44, E);
        minHFEnergy_43_44 = std::min(minHFEnergy_43_44, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,4.5])
      if ((η >=  4.3f && η <=  4.5f) ||
          (η <= -4.3f && η >= -4.5f)) {
        maxHFEnergy_43_45 = std::max(maxHFEnergy_43_45, E);
        minHFEnergy_43_45 = std::min(minHFEnergy_43_45, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,4.6])
      if ((η >=  4.3f && η <=  4.6f) ||
          (η <= -4.3f && η >= -4.6f)) {
        maxHFEnergy_43_46 = std::max(maxHFEnergy_43_46, E);
        minHFEnergy_43_46 = std::min(minHFEnergy_43_46, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,4.7])
      if ((η >=  4.3f && η <=  4.7f) ||
          (η <= -4.3f && η >= -4.7f)) {
        maxHFEnergy_43_47 = std::max(maxHFEnergy_43_47, E);
        minHFEnergy_43_47 = std::min(minHFEnergy_43_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,4.8])
      if ((η >=  4.3f && η <=  4.8f) ||
          (η <= -4.3f && η >= -4.8f)) {
        maxHFEnergy_43_48 = std::max(maxHFEnergy_43_48, E);
        minHFEnergy_43_48 = std::min(minHFEnergy_43_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,4.9])
      if ((η >=  4.3f && η <=  4.9f) ||
          (η <= -4.3f && η >= -4.9f)) {
        maxHFEnergy_43_49 = std::max(maxHFEnergy_43_49, E);
        minHFEnergy_43_49 = std::min(minHFEnergy_43_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,5.0])
      if ((η >=  4.3f && η <=  5.0f) ||
          (η <= -4.3f && η >= -5.0f)) {
        maxHFEnergy_43_50 = std::max(maxHFEnergy_43_50, E);
        minHFEnergy_43_50 = std::min(minHFEnergy_43_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,5.1])
      if ((η >=  4.3f && η <=  5.1f) ||
          (η <= -4.3f && η >= -5.1f)) {
        maxHFEnergy_43_51 = std::max(maxHFEnergy_43_51, E);
        minHFEnergy_43_51 = std::min(minHFEnergy_43_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.3,5.2])
      if ((η >=  4.3f && η <=  5.2f) ||
          (η <= -4.3f && η >= -5.2f)) {
        maxHFEnergy_43_52 = std::max(maxHFEnergy_43_52, E);
        minHFEnergy_43_52 = std::min(minHFEnergy_43_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.4,4.5])
      if ((η >=  4.4f && η <=  4.5f) ||
          (η <= -4.4f && η >= -4.5f)) {
        maxHFEnergy_44_45 = std::max(maxHFEnergy_44_45, E);
        minHFEnergy_44_45 = std::min(minHFEnergy_44_45, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,4.6])
      if ((η >=  4.4f && η <=  4.6f) ||
          (η <= -4.4f && η >= -4.6f)) {
        maxHFEnergy_44_46 = std::max(maxHFEnergy_44_46, E);
        minHFEnergy_44_46 = std::min(minHFEnergy_44_46, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,4.7])
      if ((η >=  4.4f && η <=  4.7f) ||
          (η <= -4.4f && η >= -4.7f)) {
        maxHFEnergy_44_47 = std::max(maxHFEnergy_44_47, E);
        minHFEnergy_44_47 = std::min(minHFEnergy_44_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,4.8])
      if ((η >=  4.4f && η <=  4.8f) ||
          (η <= -4.4f && η >= -4.8f)) {
        maxHFEnergy_44_48 = std::max(maxHFEnergy_44_48, E);
        minHFEnergy_44_48 = std::min(minHFEnergy_44_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,4.9])
      if ((η >=  4.4f && η <=  4.9f) ||
          (η <= -4.4f && η >= -4.9f)) {
        maxHFEnergy_44_49 = std::max(maxHFEnergy_44_49, E);
        minHFEnergy_44_49 = std::min(minHFEnergy_44_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,5.0])
      if ((η >=  4.4f && η <=  5.0f) ||
          (η <= -4.4f && η >= -5.0f)) {
        maxHFEnergy_44_50 = std::max(maxHFEnergy_44_50, E);
        minHFEnergy_44_50 = std::min(minHFEnergy_44_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,5.1])
      if ((η >=  4.4f && η <=  5.1f) ||
          (η <= -4.4f && η >= -5.1f)) {
        maxHFEnergy_44_51 = std::max(maxHFEnergy_44_51, E);
        minHFEnergy_44_51 = std::min(minHFEnergy_44_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.4,5.2])
      if ((η >=  4.4f && η <=  5.2f) ||
          (η <= -4.4f && η >= -5.2f)) {
        maxHFEnergy_44_52 = std::max(maxHFEnergy_44_52, E);
        minHFEnergy_44_52 = std::min(minHFEnergy_44_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.5,4.6])
      if ((η >=  4.5f && η <=  4.6f) ||
          (η <= -4.5f && η >= -4.6f)) {
        maxHFEnergy_45_46 = std::max(maxHFEnergy_45_46, E);
        minHFEnergy_45_46 = std::min(minHFEnergy_45_46, E);
        }
      // global HF extrema (any tower in |η|∈[4.5,4.7])
      if ((η >=  4.5f && η <=  4.7f) ||
          (η <= -4.5f && η >= -4.7f)) {
        maxHFEnergy_45_47 = std::max(maxHFEnergy_45_47, E);
        minHFEnergy_45_47 = std::min(minHFEnergy_45_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.5,4.8])
      if ((η >=  4.5f && η <=  4.8f) ||
          (η <= -4.5f && η >= -4.8f)) {
        maxHFEnergy_45_48 = std::max(maxHFEnergy_45_48, E);
        minHFEnergy_45_48 = std::min(minHFEnergy_45_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.5,4.9])
      if ((η >=  4.5f && η <=  4.9f) ||
          (η <= -4.5f && η >= -4.9f)) {
        maxHFEnergy_45_49 = std::max(maxHFEnergy_45_49, E);
        minHFEnergy_45_49 = std::min(minHFEnergy_45_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.5,5.0])
      if ((η >=  4.5f && η <=  5.0f) ||
          (η <= -4.5f && η >= -5.0f)) {
        maxHFEnergy_45_50 = std::max(maxHFEnergy_45_50, E);
        minHFEnergy_45_50 = std::min(minHFEnergy_45_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.5,5.1])
      if ((η >=  4.5f && η <=  5.1f) ||
          (η <= -4.5f && η >= -5.1f)) {
        maxHFEnergy_45_51 = std::max(maxHFEnergy_45_51, E);
        minHFEnergy_45_51 = std::min(minHFEnergy_45_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.5,5.2])
      if ((η >=  4.5f && η <=  5.2f) ||
          (η <= -4.5f && η >= -5.2f)) {
        maxHFEnergy_45_52 = std::max(maxHFEnergy_45_52, E);
        minHFEnergy_45_52 = std::min(minHFEnergy_45_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.6,4.7])
      if ((η >=  4.6f && η <=  4.7f) ||
          (η <= -4.6f && η >= -4.7f)) {
        maxHFEnergy_46_47 = std::max(maxHFEnergy_46_47, E);
        minHFEnergy_46_47 = std::min(minHFEnergy_46_47, E);
        }
      // global HF extrema (any tower in |η|∈[4.6,4.8])
      if ((η >=  4.6f && η <=  4.8f) ||
          (η <= -4.6f && η >= -4.8f)) {
        maxHFEnergy_46_48 = std::max(maxHFEnergy_46_48, E);
        minHFEnergy_46_48 = std::min(minHFEnergy_46_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.6,4.9])
      if ((η >=  4.6f && η <=  4.9f) ||
          (η <= -4.6f && η >= -4.9f)) {
        maxHFEnergy_46_49 = std::max(maxHFEnergy_46_49, E);
        minHFEnergy_46_49 = std::min(minHFEnergy_46_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.6,5.0])
      if ((η >=  4.6f && η <=  5.0f) ||
          (η <= -4.6f && η >= -5.0f)) {
        maxHFEnergy_46_50 = std::max(maxHFEnergy_46_50, E);
        minHFEnergy_46_50 = std::min(minHFEnergy_46_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.6,5.1])
      if ((η >=  4.6f && η <=  5.1f) ||
          (η <= -4.6f && η >= -5.1f)) {
        maxHFEnergy_46_51 = std::max(maxHFEnergy_46_51, E);
        minHFEnergy_46_51 = std::min(minHFEnergy_46_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.6,5.2])
      if ((η >=  4.6f && η <=  5.2f) ||
          (η <= -4.6f && η >= -5.2f)) {
        maxHFEnergy_46_52 = std::max(maxHFEnergy_46_52, E);
        minHFEnergy_46_52 = std::min(minHFEnergy_46_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.7,4.8])
      if ((η >=  4.7f && η <=  4.8f) ||
          (η <= -4.7f && η >= -4.8f)) {
        maxHFEnergy_47_48 = std::max(maxHFEnergy_47_48, E);
        minHFEnergy_47_48 = std::min(minHFEnergy_47_48, E);
        }
      // global HF extrema (any tower in |η|∈[4.7,4.9])
      if ((η >=  4.7f && η <=  4.9f) ||
          (η <= -4.7f && η >= -4.9f)) {
        maxHFEnergy_47_49 = std::max(maxHFEnergy_47_49, E);
        minHFEnergy_47_49 = std::min(minHFEnergy_47_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.7,5.0])
      if ((η >=  4.7f && η <=  5.0f) ||
          (η <= -4.7f && η >= -5.0f)) {
        maxHFEnergy_47_50 = std::max(maxHFEnergy_47_50, E);
        minHFEnergy_47_50 = std::min(minHFEnergy_47_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.7,5.1])
      if ((η >=  4.7f && η <=  5.1f) ||
          (η <= -4.7f && η >= -5.1f)) {
        maxHFEnergy_47_51 = std::max(maxHFEnergy_47_51, E);
        minHFEnergy_47_51 = std::min(minHFEnergy_47_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.7,5.2])
      if ((η >=  4.7f && η <=  5.2f) ||
          (η <= -4.7f && η >= -5.2f)) {
        maxHFEnergy_47_52 = std::max(maxHFEnergy_47_52, E);
        minHFEnergy_47_52 = std::min(minHFEnergy_47_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.8,4.9])
      if ((η >=  4.8f && η <=  4.9f) ||
          (η <= -4.8f && η >= -4.9f)) {
        maxHFEnergy_48_49 = std::max(maxHFEnergy_48_49, E);
        minHFEnergy_48_49 = std::min(minHFEnergy_48_49, E);
        }
      // global HF extrema (any tower in |η|∈[4.8,5.0])
      if ((η >=  4.8f && η <=  5.0f) ||
          (η <= -4.8f && η >= -5.0f)) {
        maxHFEnergy_48_50 = std::max(maxHFEnergy_48_50, E);
        minHFEnergy_48_50 = std::min(minHFEnergy_48_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.8,5.1])
      if ((η >=  4.8f && η <=  5.1f) ||
          (η <= -4.8f && η >= -5.1f)) {
        maxHFEnergy_48_51 = std::max(maxHFEnergy_48_51, E);
        minHFEnergy_48_51 = std::min(minHFEnergy_48_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.8,5.2])
      if ((η >=  4.8f && η <=  5.2f) ||
          (η <= -4.8f && η >= -5.2f)) {
        maxHFEnergy_48_52 = std::max(maxHFEnergy_48_52, E);
        minHFEnergy_48_52 = std::min(minHFEnergy_48_52, E);
        }


      // global HF extrema (any tower in |η|∈[4.9,5.0])
      if ((η >=  4.9f && η <=  5.0f) ||
          (η <= -4.9f && η >= -5.0f)) {
        maxHFEnergy_49_50 = std::max(maxHFEnergy_49_50, E);
        minHFEnergy_49_50 = std::min(minHFEnergy_49_50, E);
        }
      // global HF extrema (any tower in |η|∈[4.9,5.1])
      if ((η >=  4.9f && η <=  5.1f) ||
          (η <= -4.9f && η >= -5.1f)) {
        maxHFEnergy_49_51 = std::max(maxHFEnergy_49_51, E);
        minHFEnergy_49_51 = std::min(minHFEnergy_49_51, E);
        }
      // global HF extrema (any tower in |η|∈[4.9,5.2])
      if ((η >=  4.9f && η <=  5.2f) ||
          (η <= -4.9f && η >= -5.2f)) {
        maxHFEnergy_49_52 = std::max(maxHFEnergy_49_52, E);
        minHFEnergy_49_52 = std::min(minHFEnergy_49_52, E);
        }


      // global HF extrema (any tower in |η|∈[5.0,5.1])
      if ((η >=  5.0f && η <=  5.1f) ||
          (η <= -5.0f && η >= -5.1f)) {
        maxHFEnergy_50_51 = std::max(maxHFEnergy_50_51, E);
        minHFEnergy_50_51 = std::min(minHFEnergy_50_51, E);
        }
      // global HF extrema (any tower in |η|∈[5.0,5.2])
      if ((η >=  5.0f && η <=  5.2f) ||
          (η <= -5.0f && η >= -5.2f)) {
        maxHFEnergy_50_52 = std::max(maxHFEnergy_50_52, E);
        minHFEnergy_50_52 = std::min(minHFEnergy_50_52, E);
        }


      // global HF extrema (any tower in |η|∈[5.1,5.2])
      if ((η >=  5.1f && η <=  5.2f) ||
          (η <= -5.1f && η >= -5.2f)) {
        maxHFEnergy_51_52 = std::max(maxHFEnergy_51_52, E);
        minHFEnergy_51_52 = std::min(minHFEnergy_51_52, E);
        }
    }




    
    // if truly no HF at all, mark global extrema invalid
    if (maxHFEnergy == std::numeric_limits<float>::lowest()) {
      maxHFEnergy = 0.0f;
      minHFEnergy = 0.0f;
    }

    if (maxHFEnergy_30_31 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_31 = 0.0f;
      minHFEnergy_30_31 = 0.0f;
    }
    if (maxHFEnergy_30_32 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_32 = 0.0f;
      minHFEnergy_30_32 = 0.0f;
    }
    if (maxHFEnergy_30_33 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_33 = 0.0f;
      minHFEnergy_30_33 = 0.0f;
    }
    if (maxHFEnergy_30_34 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_34 = 0.0f;
      minHFEnergy_30_34 = 0.0f;
    }
    if (maxHFEnergy_30_35 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_35 = 0.0f;
      minHFEnergy_30_35 = 0.0f;
    }
    if (maxHFEnergy_30_36 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_36 = 0.0f;
      minHFEnergy_30_36 = 0.0f;
    }
    if (maxHFEnergy_30_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_37 = 0.0f;
      minHFEnergy_30_37 = 0.0f;
    }
    if (maxHFEnergy_30_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_38 = 0.0f;
      minHFEnergy_30_38 = 0.0f;
    }
    if (maxHFEnergy_30_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_39 = 0.0f;
      minHFEnergy_30_39 = 0.0f;
    }
    if (maxHFEnergy_30_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_40 = 0.0f;
      minHFEnergy_30_40 = 0.0f;
    }
    if (maxHFEnergy_30_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_41 = 0.0f;
      minHFEnergy_30_41 = 0.0f;
    }
    if (maxHFEnergy_30_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_42 = 0.0f;
      minHFEnergy_30_42 = 0.0f;
    }
    if (maxHFEnergy_30_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_43 = 0.0f;
      minHFEnergy_30_43 = 0.0f;
    }
    if (maxHFEnergy_30_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_44 = 0.0f;
      minHFEnergy_30_44 = 0.0f;
    }
    if (maxHFEnergy_30_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_45 = 0.0f;
      minHFEnergy_30_45 = 0.0f;
    }
    if (maxHFEnergy_30_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_46 = 0.0f;
      minHFEnergy_30_46 = 0.0f;
    }
    if (maxHFEnergy_30_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_47 = 0.0f;
      minHFEnergy_30_47 = 0.0f;
    }
    if (maxHFEnergy_30_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_48 = 0.0f;
      minHFEnergy_30_48 = 0.0f;
    }
    if (maxHFEnergy_30_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_49 = 0.0f;
      minHFEnergy_30_49 = 0.0f;
    }
    if (maxHFEnergy_30_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_50 = 0.0f;
      minHFEnergy_30_50 = 0.0f;
    }
    if (maxHFEnergy_30_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_51 = 0.0f;
      minHFEnergy_30_51 = 0.0f;
    }
    if (maxHFEnergy_30_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_30_52 = 0.0f;
      minHFEnergy_30_52 = 0.0f;
    }


    if (maxHFEnergy_31_32 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_32 = 0.0f;
      minHFEnergy_31_32 = 0.0f;
    }
    if (maxHFEnergy_31_33 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_33 = 0.0f;
      minHFEnergy_31_33 = 0.0f;
    }
    if (maxHFEnergy_31_34 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_34 = 0.0f;
      minHFEnergy_31_34 = 0.0f;
    }
    if (maxHFEnergy_31_35 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_35 = 0.0f;
      minHFEnergy_31_35 = 0.0f;
    }
    if (maxHFEnergy_31_36 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_36 = 0.0f;
      minHFEnergy_31_36 = 0.0f;
    }
    if (maxHFEnergy_31_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_37 = 0.0f;
      minHFEnergy_31_37 = 0.0f;
    }
    if (maxHFEnergy_31_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_38 = 0.0f;
      minHFEnergy_31_38 = 0.0f;
    }
    if (maxHFEnergy_31_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_39 = 0.0f;
      minHFEnergy_31_39 = 0.0f;
    }
    if (maxHFEnergy_31_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_40 = 0.0f;
      minHFEnergy_31_40 = 0.0f;
    }
    if (maxHFEnergy_31_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_41 = 0.0f;
      minHFEnergy_31_41 = 0.0f;
    }
    if (maxHFEnergy_31_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_42 = 0.0f;
      minHFEnergy_31_42 = 0.0f;
    }
    if (maxHFEnergy_31_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_43 = 0.0f;
      minHFEnergy_31_43 = 0.0f;
    }
    if (maxHFEnergy_31_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_44 = 0.0f;
      minHFEnergy_31_44 = 0.0f;
    }
    if (maxHFEnergy_31_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_45 = 0.0f;
      minHFEnergy_31_45 = 0.0f;
    }
    if (maxHFEnergy_31_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_46 = 0.0f;
      minHFEnergy_31_46 = 0.0f;
    }
    if (maxHFEnergy_31_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_47 = 0.0f;
      minHFEnergy_31_47 = 0.0f;
    }
    if (maxHFEnergy_31_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_48 = 0.0f;
      minHFEnergy_31_48 = 0.0f;
    }
    if (maxHFEnergy_31_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_49 = 0.0f;
      minHFEnergy_31_49 = 0.0f;
    }
    if (maxHFEnergy_31_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_50 = 0.0f;
      minHFEnergy_31_50 = 0.0f;
    }
    if (maxHFEnergy_31_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_51 = 0.0f;
      minHFEnergy_31_51 = 0.0f;
    }
    if (maxHFEnergy_31_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_31_52 = 0.0f;
      minHFEnergy_31_52 = 0.0f;
    }


    if (maxHFEnergy_32_33 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_33 = 0.0f;
      minHFEnergy_32_33 = 0.0f;
    }
    if (maxHFEnergy_32_34 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_34 = 0.0f;
      minHFEnergy_32_34 = 0.0f;
    }
    if (maxHFEnergy_32_35 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_35 = 0.0f;
      minHFEnergy_32_35 = 0.0f;
    }
    if (maxHFEnergy_32_36 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_36 = 0.0f;
      minHFEnergy_32_36 = 0.0f;
    }
    if (maxHFEnergy_32_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_37 = 0.0f;
      minHFEnergy_32_37 = 0.0f;
    }
    if (maxHFEnergy_32_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_38 = 0.0f;
      minHFEnergy_32_38 = 0.0f;
    }
    if (maxHFEnergy_32_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_39 = 0.0f;
      minHFEnergy_32_39 = 0.0f;
    }
    if (maxHFEnergy_32_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_40 = 0.0f;
      minHFEnergy_32_40 = 0.0f;
    }
    if (maxHFEnergy_32_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_41 = 0.0f;
      minHFEnergy_32_41 = 0.0f;
    }
    if (maxHFEnergy_32_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_42 = 0.0f;
      minHFEnergy_32_42 = 0.0f;
    }
    if (maxHFEnergy_32_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_43 = 0.0f;
      minHFEnergy_32_43 = 0.0f;
    }
    if (maxHFEnergy_32_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_44 = 0.0f;
      minHFEnergy_32_44 = 0.0f;
    }
    if (maxHFEnergy_32_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_45 = 0.0f;
      minHFEnergy_32_45 = 0.0f;
    }
    if (maxHFEnergy_32_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_46 = 0.0f;
      minHFEnergy_32_46 = 0.0f;
    }
    if (maxHFEnergy_32_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_47 = 0.0f;
      minHFEnergy_32_47 = 0.0f;
    }
    if (maxHFEnergy_32_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_48 = 0.0f;
      minHFEnergy_32_48 = 0.0f;
    }
    if (maxHFEnergy_32_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_49 = 0.0f;
      minHFEnergy_32_49 = 0.0f;
    }
    if (maxHFEnergy_32_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_50 = 0.0f;
      minHFEnergy_32_50 = 0.0f;
    }
    if (maxHFEnergy_32_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_51 = 0.0f;
      minHFEnergy_32_51 = 0.0f;
    }
    if (maxHFEnergy_32_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_32_52 = 0.0f;
      minHFEnergy_32_52 = 0.0f;
    }


    if (maxHFEnergy_33_34 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_34 = 0.0f;
      minHFEnergy_33_34 = 0.0f;
    }
    if (maxHFEnergy_33_35 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_35 = 0.0f;
      minHFEnergy_33_35 = 0.0f;
    }
    if (maxHFEnergy_33_36 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_36 = 0.0f;
      minHFEnergy_33_36 = 0.0f;
    }
    if (maxHFEnergy_33_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_37 = 0.0f;
      minHFEnergy_33_37 = 0.0f;
    }
    if (maxHFEnergy_33_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_38 = 0.0f;
      minHFEnergy_33_38 = 0.0f;
    }
    if (maxHFEnergy_33_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_39 = 0.0f;
      minHFEnergy_33_39 = 0.0f;
    }
    if (maxHFEnergy_33_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_40 = 0.0f;
      minHFEnergy_33_40 = 0.0f;
    }
    if (maxHFEnergy_33_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_41 = 0.0f;
      minHFEnergy_33_41 = 0.0f;
    }
    if (maxHFEnergy_33_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_42 = 0.0f;
      minHFEnergy_33_42 = 0.0f;
    }
    if (maxHFEnergy_33_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_43 = 0.0f;
      minHFEnergy_33_43 = 0.0f;
    }
    if (maxHFEnergy_33_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_44 = 0.0f;
      minHFEnergy_33_44 = 0.0f;
    }
    if (maxHFEnergy_33_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_45 = 0.0f;
      minHFEnergy_33_45 = 0.0f;
    }
    if (maxHFEnergy_33_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_46 = 0.0f;
      minHFEnergy_33_46 = 0.0f;
    }
    if (maxHFEnergy_33_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_47 = 0.0f;
      minHFEnergy_33_47 = 0.0f;
    }
    if (maxHFEnergy_33_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_48 = 0.0f;
      minHFEnergy_33_48 = 0.0f;
    }
    if (maxHFEnergy_33_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_49 = 0.0f;
      minHFEnergy_33_49 = 0.0f;
    }
    if (maxHFEnergy_33_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_50 = 0.0f;
      minHFEnergy_33_50 = 0.0f;
    }
    if (maxHFEnergy_33_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_51 = 0.0f;
      minHFEnergy_33_51 = 0.0f;
    }
    if (maxHFEnergy_33_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_33_52 = 0.0f;
      minHFEnergy_33_52 = 0.0f;
    }


    if (maxHFEnergy_34_35 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_35 = 0.0f;
      minHFEnergy_34_35 = 0.0f;
    }
    if (maxHFEnergy_34_36 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_36 = 0.0f;
      minHFEnergy_34_36 = 0.0f;
    }
    if (maxHFEnergy_34_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_37 = 0.0f;
      minHFEnergy_34_37 = 0.0f;
    }
    if (maxHFEnergy_34_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_38 = 0.0f;
      minHFEnergy_34_38 = 0.0f;
    }
    if (maxHFEnergy_34_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_39 = 0.0f;
      minHFEnergy_34_39 = 0.0f;
    }
    if (maxHFEnergy_34_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_40 = 0.0f;
      minHFEnergy_34_40 = 0.0f;
    }
    if (maxHFEnergy_34_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_41 = 0.0f;
      minHFEnergy_34_41 = 0.0f;
    }
    if (maxHFEnergy_34_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_42 = 0.0f;
      minHFEnergy_34_42 = 0.0f;
    }
    if (maxHFEnergy_34_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_43 = 0.0f;
      minHFEnergy_34_43 = 0.0f;
    }
    if (maxHFEnergy_34_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_44 = 0.0f;
      minHFEnergy_34_44 = 0.0f;
    }
    if (maxHFEnergy_34_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_45 = 0.0f;
      minHFEnergy_34_45 = 0.0f;
    }
    if (maxHFEnergy_34_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_46 = 0.0f;
      minHFEnergy_34_46 = 0.0f;
    }
    if (maxHFEnergy_34_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_47 = 0.0f;
      minHFEnergy_34_47 = 0.0f;
    }
    if (maxHFEnergy_34_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_48 = 0.0f;
      minHFEnergy_34_48 = 0.0f;
    }
    if (maxHFEnergy_34_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_49 = 0.0f;
      minHFEnergy_34_49 = 0.0f;
    }
    if (maxHFEnergy_34_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_50 = 0.0f;
      minHFEnergy_34_50 = 0.0f;
    }
    if (maxHFEnergy_34_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_51 = 0.0f;
      minHFEnergy_34_51 = 0.0f;
    }
    if (maxHFEnergy_34_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_34_52 = 0.0f;
      minHFEnergy_34_52 = 0.0f;
    }


    if (maxHFEnergy_35_36 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_36 = 0.0f;
      minHFEnergy_35_36 = 0.0f;
    }
    if (maxHFEnergy_35_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_37 = 0.0f;
      minHFEnergy_35_37 = 0.0f;
    }
    if (maxHFEnergy_35_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_38 = 0.0f;
      minHFEnergy_35_38 = 0.0f;
    }
    if (maxHFEnergy_35_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_39 = 0.0f;
      minHFEnergy_35_39 = 0.0f;
    }
    if (maxHFEnergy_35_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_40 = 0.0f;
      minHFEnergy_35_40 = 0.0f;
    }
    if (maxHFEnergy_35_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_41 = 0.0f;
      minHFEnergy_35_41 = 0.0f;
    }
    if (maxHFEnergy_35_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_42 = 0.0f;
      minHFEnergy_35_42 = 0.0f;
    }
    if (maxHFEnergy_35_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_43 = 0.0f;
      minHFEnergy_35_43 = 0.0f;
    }
    if (maxHFEnergy_35_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_44 = 0.0f;
      minHFEnergy_35_44 = 0.0f;
    }
    if (maxHFEnergy_35_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_45 = 0.0f;
      minHFEnergy_35_45 = 0.0f;
    }
    if (maxHFEnergy_35_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_46 = 0.0f;
      minHFEnergy_35_46 = 0.0f;
    }
    if (maxHFEnergy_35_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_47 = 0.0f;
      minHFEnergy_35_47 = 0.0f;
    }
    if (maxHFEnergy_35_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_48 = 0.0f;
      minHFEnergy_35_48 = 0.0f;
    }
    if (maxHFEnergy_35_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_49 = 0.0f;
      minHFEnergy_35_49 = 0.0f;
    }
    if (maxHFEnergy_35_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_50 = 0.0f;
      minHFEnergy_35_50 = 0.0f;
    }
    if (maxHFEnergy_35_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_51 = 0.0f;
      minHFEnergy_35_51 = 0.0f;
    }
    if (maxHFEnergy_35_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_35_52 = 0.0f;
      minHFEnergy_35_52 = 0.0f;
    }


    if (maxHFEnergy_36_37 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_37 = 0.0f;
      minHFEnergy_36_37 = 0.0f;
    }
    if (maxHFEnergy_36_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_38 = 0.0f;
      minHFEnergy_36_38 = 0.0f;
    }
    if (maxHFEnergy_36_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_39 = 0.0f;
      minHFEnergy_36_39 = 0.0f;
    }
    if (maxHFEnergy_36_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_40 = 0.0f;
      minHFEnergy_36_40 = 0.0f;
    }
    if (maxHFEnergy_36_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_41 = 0.0f;
      minHFEnergy_36_41 = 0.0f;
    }
    if (maxHFEnergy_36_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_42 = 0.0f;
      minHFEnergy_36_42 = 0.0f;
    }
    if (maxHFEnergy_36_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_43 = 0.0f;
      minHFEnergy_36_43 = 0.0f;
    }
    if (maxHFEnergy_36_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_44 = 0.0f;
      minHFEnergy_36_44 = 0.0f;
    }
    if (maxHFEnergy_36_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_45 = 0.0f;
      minHFEnergy_36_45 = 0.0f;
    }
    if (maxHFEnergy_36_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_46 = 0.0f;
      minHFEnergy_36_46 = 0.0f;
    }
    if (maxHFEnergy_36_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_47 = 0.0f;
      minHFEnergy_36_47 = 0.0f;
    }
    if (maxHFEnergy_36_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_48 = 0.0f;
      minHFEnergy_36_48 = 0.0f;
    }
    if (maxHFEnergy_36_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_49 = 0.0f;
      minHFEnergy_36_49 = 0.0f;
    }
    if (maxHFEnergy_36_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_50 = 0.0f;
      minHFEnergy_36_50 = 0.0f;
    }
    if (maxHFEnergy_36_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_51 = 0.0f;
      minHFEnergy_36_51 = 0.0f;
    }
    if (maxHFEnergy_36_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_36_52 = 0.0f;
      minHFEnergy_36_52 = 0.0f;
    }


    if (maxHFEnergy_37_38 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_38 = 0.0f;
      minHFEnergy_37_38 = 0.0f;
    }
    if (maxHFEnergy_37_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_39 = 0.0f;
      minHFEnergy_37_39 = 0.0f;
    }
    if (maxHFEnergy_37_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_40 = 0.0f;
      minHFEnergy_37_40 = 0.0f;
    }
    if (maxHFEnergy_37_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_41 = 0.0f;
      minHFEnergy_37_41 = 0.0f;
    }
    if (maxHFEnergy_37_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_42 = 0.0f;
      minHFEnergy_37_42 = 0.0f;
    }
    if (maxHFEnergy_37_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_43 = 0.0f;
      minHFEnergy_37_43 = 0.0f;
    }
    if (maxHFEnergy_37_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_44 = 0.0f;
      minHFEnergy_37_44 = 0.0f;
    }
    if (maxHFEnergy_37_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_45 = 0.0f;
      minHFEnergy_37_45 = 0.0f;
    }
    if (maxHFEnergy_37_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_46 = 0.0f;
      minHFEnergy_37_46 = 0.0f;
    }
    if (maxHFEnergy_37_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_47 = 0.0f;
      minHFEnergy_37_47 = 0.0f;
    }
    if (maxHFEnergy_37_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_48 = 0.0f;
      minHFEnergy_37_48 = 0.0f;
    }
    if (maxHFEnergy_37_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_49 = 0.0f;
      minHFEnergy_37_49 = 0.0f;
    }
    if (maxHFEnergy_37_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_50 = 0.0f;
      minHFEnergy_37_50 = 0.0f;
    }
    if (maxHFEnergy_37_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_51 = 0.0f;
      minHFEnergy_37_51 = 0.0f;
    }
    if (maxHFEnergy_37_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_37_52 = 0.0f;
      minHFEnergy_37_52 = 0.0f;
    }


    if (maxHFEnergy_38_39 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_39 = 0.0f;
      minHFEnergy_38_39 = 0.0f;
    }
    if (maxHFEnergy_38_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_40 = 0.0f;
      minHFEnergy_38_40 = 0.0f;
    }
    if (maxHFEnergy_38_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_41 = 0.0f;
      minHFEnergy_38_41 = 0.0f;
    }
    if (maxHFEnergy_38_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_42 = 0.0f;
      minHFEnergy_38_42 = 0.0f;
    }
    if (maxHFEnergy_38_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_43 = 0.0f;
      minHFEnergy_38_43 = 0.0f;
    }
    if (maxHFEnergy_38_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_44 = 0.0f;
      minHFEnergy_38_44 = 0.0f;
    }
    if (maxHFEnergy_38_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_45 = 0.0f;
      minHFEnergy_38_45 = 0.0f;
    }
    if (maxHFEnergy_38_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_46 = 0.0f;
      minHFEnergy_38_46 = 0.0f;
    }
    if (maxHFEnergy_38_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_47 = 0.0f;
      minHFEnergy_38_47 = 0.0f;
    }
    if (maxHFEnergy_38_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_48 = 0.0f;
      minHFEnergy_38_48 = 0.0f;
    }
    if (maxHFEnergy_38_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_49 = 0.0f;
      minHFEnergy_38_49 = 0.0f;
    }
    if (maxHFEnergy_38_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_50 = 0.0f;
      minHFEnergy_38_50 = 0.0f;
    }
    if (maxHFEnergy_38_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_51 = 0.0f;
      minHFEnergy_38_51 = 0.0f;
    }
    if (maxHFEnergy_38_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_38_52 = 0.0f;
      minHFEnergy_38_52 = 0.0f;
    }


    if (maxHFEnergy_39_40 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_40 = 0.0f;
      minHFEnergy_39_40 = 0.0f;
    }
    if (maxHFEnergy_39_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_41 = 0.0f;
      minHFEnergy_39_41 = 0.0f;
    }
    if (maxHFEnergy_39_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_42 = 0.0f;
      minHFEnergy_39_42 = 0.0f;
    }
    if (maxHFEnergy_39_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_43 = 0.0f;
      minHFEnergy_39_43 = 0.0f;
    }
    if (maxHFEnergy_39_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_44 = 0.0f;
      minHFEnergy_39_44 = 0.0f;
    }
    if (maxHFEnergy_39_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_45 = 0.0f;
      minHFEnergy_39_45 = 0.0f;
    }
    if (maxHFEnergy_39_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_46 = 0.0f;
      minHFEnergy_39_46 = 0.0f;
    }
    if (maxHFEnergy_39_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_47 = 0.0f;
      minHFEnergy_39_47 = 0.0f;
    }
    if (maxHFEnergy_39_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_48 = 0.0f;
      minHFEnergy_39_48 = 0.0f;
    }
    if (maxHFEnergy_39_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_49 = 0.0f;
      minHFEnergy_39_49 = 0.0f;
    }
    if (maxHFEnergy_39_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_50 = 0.0f;
      minHFEnergy_39_50 = 0.0f;
    }
    if (maxHFEnergy_39_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_51 = 0.0f;
      minHFEnergy_39_51 = 0.0f;
    }
    if (maxHFEnergy_39_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_39_52 = 0.0f;
      minHFEnergy_39_52 = 0.0f;
    }


    if (maxHFEnergy_40_41 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_41 = 0.0f;
      minHFEnergy_40_41 = 0.0f;
    }
    if (maxHFEnergy_40_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_42 = 0.0f;
      minHFEnergy_40_42 = 0.0f;
    }
    if (maxHFEnergy_40_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_43 = 0.0f;
      minHFEnergy_40_43 = 0.0f;
    }
    if (maxHFEnergy_40_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_44 = 0.0f;
      minHFEnergy_40_44 = 0.0f;
    }
    if (maxHFEnergy_40_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_45 = 0.0f;
      minHFEnergy_40_45 = 0.0f;
    }
    if (maxHFEnergy_40_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_46 = 0.0f;
      minHFEnergy_40_46 = 0.0f;
    }
    if (maxHFEnergy_40_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_47 = 0.0f;
      minHFEnergy_40_47 = 0.0f;
    }
    if (maxHFEnergy_40_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_48 = 0.0f;
      minHFEnergy_40_48 = 0.0f;
    }
    if (maxHFEnergy_40_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_49 = 0.0f;
      minHFEnergy_40_49 = 0.0f;
    }
    if (maxHFEnergy_40_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_50 = 0.0f;
      minHFEnergy_40_50 = 0.0f;
    }
    if (maxHFEnergy_40_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_51 = 0.0f;
      minHFEnergy_40_51 = 0.0f;
    }
    if (maxHFEnergy_40_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_40_52 = 0.0f;
      minHFEnergy_40_52 = 0.0f;
    }


    if (maxHFEnergy_41_42 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_42 = 0.0f;
      minHFEnergy_41_42 = 0.0f;
    }
    if (maxHFEnergy_41_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_43 = 0.0f;
      minHFEnergy_41_43 = 0.0f;
    }
    if (maxHFEnergy_41_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_44 = 0.0f;
      minHFEnergy_41_44 = 0.0f;
    }
    if (maxHFEnergy_41_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_45 = 0.0f;
      minHFEnergy_41_45 = 0.0f;
    }
    if (maxHFEnergy_41_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_46 = 0.0f;
      minHFEnergy_41_46 = 0.0f;
    }
    if (maxHFEnergy_41_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_47 = 0.0f;
      minHFEnergy_41_47 = 0.0f;
    }
    if (maxHFEnergy_41_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_48 = 0.0f;
      minHFEnergy_41_48 = 0.0f;
    }
    if (maxHFEnergy_41_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_49 = 0.0f;
      minHFEnergy_41_49 = 0.0f;
    }
    if (maxHFEnergy_41_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_50 = 0.0f;
      minHFEnergy_41_50 = 0.0f;
    }
    if (maxHFEnergy_41_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_51 = 0.0f;
      minHFEnergy_41_51 = 0.0f;
    }
    if (maxHFEnergy_41_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_41_52 = 0.0f;
      minHFEnergy_41_52 = 0.0f;
    }


    if (maxHFEnergy_42_43 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_43 = 0.0f;
      minHFEnergy_42_43 = 0.0f;
    }
    if (maxHFEnergy_42_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_44 = 0.0f;
      minHFEnergy_42_44 = 0.0f;
    }
    if (maxHFEnergy_42_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_45 = 0.0f;
      minHFEnergy_42_45 = 0.0f;
    }
    if (maxHFEnergy_42_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_46 = 0.0f;
      minHFEnergy_42_46 = 0.0f;
    }
    if (maxHFEnergy_42_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_47 = 0.0f;
      minHFEnergy_42_47 = 0.0f;
    }
    if (maxHFEnergy_42_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_48 = 0.0f;
      minHFEnergy_42_48 = 0.0f;
    }
    if (maxHFEnergy_42_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_49 = 0.0f;
      minHFEnergy_42_49 = 0.0f;
    }
    if (maxHFEnergy_42_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_50 = 0.0f;
      minHFEnergy_42_50 = 0.0f;
    }
    if (maxHFEnergy_42_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_51 = 0.0f;
      minHFEnergy_42_51 = 0.0f;
    }
    if (maxHFEnergy_42_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_42_52 = 0.0f;
      minHFEnergy_42_52 = 0.0f;
    }


    if (maxHFEnergy_43_44 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_44 = 0.0f;
      minHFEnergy_43_44 = 0.0f;
    }
    if (maxHFEnergy_43_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_45 = 0.0f;
      minHFEnergy_43_45 = 0.0f;
    }
    if (maxHFEnergy_43_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_46 = 0.0f;
      minHFEnergy_43_46 = 0.0f;
    }
    if (maxHFEnergy_43_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_47 = 0.0f;
      minHFEnergy_43_47 = 0.0f;
    }
    if (maxHFEnergy_43_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_48 = 0.0f;
      minHFEnergy_43_48 = 0.0f;
    }
    if (maxHFEnergy_43_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_49 = 0.0f;
      minHFEnergy_43_49 = 0.0f;
    }
    if (maxHFEnergy_43_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_50 = 0.0f;
      minHFEnergy_43_50 = 0.0f;
    }
    if (maxHFEnergy_43_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_51 = 0.0f;
      minHFEnergy_43_51 = 0.0f;
    }
    if (maxHFEnergy_43_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_43_52 = 0.0f;
      minHFEnergy_43_52 = 0.0f;
    }


    if (maxHFEnergy_44_45 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_45 = 0.0f;
      minHFEnergy_44_45 = 0.0f;
    }
    if (maxHFEnergy_44_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_46 = 0.0f;
      minHFEnergy_44_46 = 0.0f;
    }
    if (maxHFEnergy_44_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_47 = 0.0f;
      minHFEnergy_44_47 = 0.0f;
    }
    if (maxHFEnergy_44_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_48 = 0.0f;
      minHFEnergy_44_48 = 0.0f;
    }
    if (maxHFEnergy_44_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_49 = 0.0f;
      minHFEnergy_44_49 = 0.0f;
    }
    if (maxHFEnergy_44_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_50 = 0.0f;
      minHFEnergy_44_50 = 0.0f;
    }
    if (maxHFEnergy_44_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_51 = 0.0f;
      minHFEnergy_44_51 = 0.0f;
    }
    if (maxHFEnergy_44_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_44_52 = 0.0f;
      minHFEnergy_44_52 = 0.0f;
    }


    if (maxHFEnergy_45_46 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_46 = 0.0f;
      minHFEnergy_45_46 = 0.0f;
    }
    if (maxHFEnergy_45_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_47 = 0.0f;
      minHFEnergy_45_47 = 0.0f;
    }
    if (maxHFEnergy_45_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_48 = 0.0f;
      minHFEnergy_45_48 = 0.0f;
    }
    if (maxHFEnergy_45_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_49 = 0.0f;
      minHFEnergy_45_49 = 0.0f;
    }
    if (maxHFEnergy_45_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_50 = 0.0f;
      minHFEnergy_45_50 = 0.0f;
    }
    if (maxHFEnergy_45_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_51 = 0.0f;
      minHFEnergy_45_51 = 0.0f;
    }
    if (maxHFEnergy_45_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_45_52 = 0.0f;
      minHFEnergy_45_52 = 0.0f;
    }


    if (maxHFEnergy_46_47 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_46_47 = 0.0f;
      minHFEnergy_46_47 = 0.0f;
    }
    if (maxHFEnergy_46_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_46_48 = 0.0f;
      minHFEnergy_46_48 = 0.0f;
    }
    if (maxHFEnergy_46_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_46_49 = 0.0f;
      minHFEnergy_46_49 = 0.0f;
    }
    if (maxHFEnergy_46_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_46_50 = 0.0f;
      minHFEnergy_46_50 = 0.0f;
    }
    if (maxHFEnergy_46_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_46_51 = 0.0f;
      minHFEnergy_46_51 = 0.0f;
    }
    if (maxHFEnergy_46_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_46_52 = 0.0f;
      minHFEnergy_46_52 = 0.0f;
    }


    if (maxHFEnergy_47_48 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_47_48 = 0.0f;
      minHFEnergy_47_48 = 0.0f;
    }
    if (maxHFEnergy_47_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_47_49 = 0.0f;
      minHFEnergy_47_49 = 0.0f;
    }
    if (maxHFEnergy_47_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_47_50 = 0.0f;
      minHFEnergy_47_50 = 0.0f;
    }
    if (maxHFEnergy_47_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_47_51 = 0.0f;
      minHFEnergy_47_51 = 0.0f;
    }
    if (maxHFEnergy_47_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_47_52 = 0.0f;
      minHFEnergy_47_52 = 0.0f;
    }


    if (maxHFEnergy_48_49 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_48_49 = 0.0f;
      minHFEnergy_48_49 = 0.0f;
    }
    if (maxHFEnergy_48_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_48_50 = 0.0f;
      minHFEnergy_48_50 = 0.0f;
    }
    if (maxHFEnergy_48_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_48_51 = 0.0f;
      minHFEnergy_48_51 = 0.0f;
    }
    if (maxHFEnergy_48_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_48_52 = 0.0f;
      minHFEnergy_48_52 = 0.0f;
    }


    if (maxHFEnergy_49_50 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_49_50 = 0.0f;
      minHFEnergy_49_50 = 0.0f;
    }
    if (maxHFEnergy_49_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_49_51 = 0.0f;
      minHFEnergy_49_51 = 0.0f;
    }
    if (maxHFEnergy_49_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_49_52 = 0.0f;
      minHFEnergy_49_52 = 0.0f;
    }


    if (maxHFEnergy_50_51 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_50_51 = 0.0f;
      minHFEnergy_50_51 = 0.0f;
    }
    if (maxHFEnergy_50_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_50_52 = 0.0f;
      minHFEnergy_50_52 = 0.0f;
    }


    if (maxHFEnergy_51_52 == std::numeric_limits<float>::lowest()) {
      maxHFEnergy_51_52 = 0.0f;
      minHFEnergy_51_52 = 0.0f;
    }




    
    // derived variables (now always ≥ 0)
    minOfMaxHF  = std::min(maxHFEnergyPlus,  maxHFEnergyMinus);
    minOfHFSums = std::min(hfSumPlus,        hfSumMinus);



    minOfMaxHF_30_31  = std::min(maxHFEnergyPlus_30_31,  maxHFEnergyMinus_30_31);
    minOfMaxHF_30_32  = std::min(maxHFEnergyPlus_30_32,  maxHFEnergyMinus_30_32);
    minOfMaxHF_30_33  = std::min(maxHFEnergyPlus_30_33,  maxHFEnergyMinus_30_33);
    minOfMaxHF_30_34  = std::min(maxHFEnergyPlus_30_34,  maxHFEnergyMinus_30_34);
    minOfMaxHF_30_35  = std::min(maxHFEnergyPlus_30_35,  maxHFEnergyMinus_30_35);
    minOfMaxHF_30_36  = std::min(maxHFEnergyPlus_30_36,  maxHFEnergyMinus_30_36);
    minOfMaxHF_30_37  = std::min(maxHFEnergyPlus_30_37,  maxHFEnergyMinus_30_37);
    minOfMaxHF_30_38  = std::min(maxHFEnergyPlus_30_38,  maxHFEnergyMinus_30_38);
    minOfMaxHF_30_39  = std::min(maxHFEnergyPlus_30_39,  maxHFEnergyMinus_30_39);
    minOfMaxHF_30_40  = std::min(maxHFEnergyPlus_30_40,  maxHFEnergyMinus_30_40);
    minOfMaxHF_30_41  = std::min(maxHFEnergyPlus_30_41,  maxHFEnergyMinus_30_41);
    minOfMaxHF_30_42  = std::min(maxHFEnergyPlus_30_42,  maxHFEnergyMinus_30_42);
    minOfMaxHF_30_43  = std::min(maxHFEnergyPlus_30_43,  maxHFEnergyMinus_30_43);
    minOfMaxHF_30_44  = std::min(maxHFEnergyPlus_30_44,  maxHFEnergyMinus_30_44);
    minOfMaxHF_30_45  = std::min(maxHFEnergyPlus_30_45,  maxHFEnergyMinus_30_45);
    minOfMaxHF_30_46  = std::min(maxHFEnergyPlus_30_46,  maxHFEnergyMinus_30_46);
    minOfMaxHF_30_47  = std::min(maxHFEnergyPlus_30_47,  maxHFEnergyMinus_30_47);
    minOfMaxHF_30_48  = std::min(maxHFEnergyPlus_30_48,  maxHFEnergyMinus_30_48);
    minOfMaxHF_30_49  = std::min(maxHFEnergyPlus_30_49,  maxHFEnergyMinus_30_49);
    minOfMaxHF_30_50  = std::min(maxHFEnergyPlus_30_50,  maxHFEnergyMinus_30_50);
    minOfMaxHF_30_51  = std::min(maxHFEnergyPlus_30_51,  maxHFEnergyMinus_30_51);
    minOfMaxHF_30_52  = std::min(maxHFEnergyPlus_30_52,  maxHFEnergyMinus_30_52);


    minOfMaxHF_31_32  = std::min(maxHFEnergyPlus_31_32,  maxHFEnergyMinus_31_32);
    minOfMaxHF_31_33  = std::min(maxHFEnergyPlus_31_33,  maxHFEnergyMinus_31_33);
    minOfMaxHF_31_34  = std::min(maxHFEnergyPlus_31_34,  maxHFEnergyMinus_31_34);
    minOfMaxHF_31_35  = std::min(maxHFEnergyPlus_31_35,  maxHFEnergyMinus_31_35);
    minOfMaxHF_31_36  = std::min(maxHFEnergyPlus_31_36,  maxHFEnergyMinus_31_36);
    minOfMaxHF_31_37  = std::min(maxHFEnergyPlus_31_37,  maxHFEnergyMinus_31_37);
    minOfMaxHF_31_38  = std::min(maxHFEnergyPlus_31_38,  maxHFEnergyMinus_31_38);
    minOfMaxHF_31_39  = std::min(maxHFEnergyPlus_31_39,  maxHFEnergyMinus_31_39);
    minOfMaxHF_31_40  = std::min(maxHFEnergyPlus_31_40,  maxHFEnergyMinus_31_40);
    minOfMaxHF_31_41  = std::min(maxHFEnergyPlus_31_41,  maxHFEnergyMinus_31_41);
    minOfMaxHF_31_42  = std::min(maxHFEnergyPlus_31_42,  maxHFEnergyMinus_31_42);
    minOfMaxHF_31_43  = std::min(maxHFEnergyPlus_31_43,  maxHFEnergyMinus_31_43);
    minOfMaxHF_31_44  = std::min(maxHFEnergyPlus_31_44,  maxHFEnergyMinus_31_44);
    minOfMaxHF_31_45  = std::min(maxHFEnergyPlus_31_45,  maxHFEnergyMinus_31_45);
    minOfMaxHF_31_46  = std::min(maxHFEnergyPlus_31_46,  maxHFEnergyMinus_31_46);
    minOfMaxHF_31_47  = std::min(maxHFEnergyPlus_31_47,  maxHFEnergyMinus_31_47);
    minOfMaxHF_31_48  = std::min(maxHFEnergyPlus_31_48,  maxHFEnergyMinus_31_48);
    minOfMaxHF_31_49  = std::min(maxHFEnergyPlus_31_49,  maxHFEnergyMinus_31_49);
    minOfMaxHF_31_50  = std::min(maxHFEnergyPlus_31_50,  maxHFEnergyMinus_31_50);
    minOfMaxHF_31_51  = std::min(maxHFEnergyPlus_31_51,  maxHFEnergyMinus_31_51);
    minOfMaxHF_31_52  = std::min(maxHFEnergyPlus_31_52,  maxHFEnergyMinus_31_52);


    minOfMaxHF_32_33  = std::min(maxHFEnergyPlus_32_33,  maxHFEnergyMinus_32_33);
    minOfMaxHF_32_34  = std::min(maxHFEnergyPlus_32_34,  maxHFEnergyMinus_32_34);
    minOfMaxHF_32_35  = std::min(maxHFEnergyPlus_32_35,  maxHFEnergyMinus_32_35);
    minOfMaxHF_32_36  = std::min(maxHFEnergyPlus_32_36,  maxHFEnergyMinus_32_36);
    minOfMaxHF_32_37  = std::min(maxHFEnergyPlus_32_37,  maxHFEnergyMinus_32_37);
    minOfMaxHF_32_38  = std::min(maxHFEnergyPlus_32_38,  maxHFEnergyMinus_32_38);
    minOfMaxHF_32_39  = std::min(maxHFEnergyPlus_32_39,  maxHFEnergyMinus_32_39);
    minOfMaxHF_32_40  = std::min(maxHFEnergyPlus_32_40,  maxHFEnergyMinus_32_40);
    minOfMaxHF_32_41  = std::min(maxHFEnergyPlus_32_41,  maxHFEnergyMinus_32_41);
    minOfMaxHF_32_42  = std::min(maxHFEnergyPlus_32_42,  maxHFEnergyMinus_32_42);
    minOfMaxHF_32_43  = std::min(maxHFEnergyPlus_32_43,  maxHFEnergyMinus_32_43);
    minOfMaxHF_32_44  = std::min(maxHFEnergyPlus_32_44,  maxHFEnergyMinus_32_44);
    minOfMaxHF_32_45  = std::min(maxHFEnergyPlus_32_45,  maxHFEnergyMinus_32_45);
    minOfMaxHF_32_46  = std::min(maxHFEnergyPlus_32_46,  maxHFEnergyMinus_32_46);
    minOfMaxHF_32_47  = std::min(maxHFEnergyPlus_32_47,  maxHFEnergyMinus_32_47);
    minOfMaxHF_32_48  = std::min(maxHFEnergyPlus_32_48,  maxHFEnergyMinus_32_48);
    minOfMaxHF_32_49  = std::min(maxHFEnergyPlus_32_49,  maxHFEnergyMinus_32_49);
    minOfMaxHF_32_50  = std::min(maxHFEnergyPlus_32_50,  maxHFEnergyMinus_32_50);
    minOfMaxHF_32_51  = std::min(maxHFEnergyPlus_32_51,  maxHFEnergyMinus_32_51);
    minOfMaxHF_32_52  = std::min(maxHFEnergyPlus_32_52,  maxHFEnergyMinus_32_52);


    minOfMaxHF_33_34  = std::min(maxHFEnergyPlus_33_34,  maxHFEnergyMinus_33_34);
    minOfMaxHF_33_35  = std::min(maxHFEnergyPlus_33_35,  maxHFEnergyMinus_33_35);
    minOfMaxHF_33_36  = std::min(maxHFEnergyPlus_33_36,  maxHFEnergyMinus_33_36);
    minOfMaxHF_33_37  = std::min(maxHFEnergyPlus_33_37,  maxHFEnergyMinus_33_37);
    minOfMaxHF_33_38  = std::min(maxHFEnergyPlus_33_38,  maxHFEnergyMinus_33_38);
    minOfMaxHF_33_39  = std::min(maxHFEnergyPlus_33_39,  maxHFEnergyMinus_33_39);
    minOfMaxHF_33_40  = std::min(maxHFEnergyPlus_33_40,  maxHFEnergyMinus_33_40);
    minOfMaxHF_33_41  = std::min(maxHFEnergyPlus_33_41,  maxHFEnergyMinus_33_41);
    minOfMaxHF_33_42  = std::min(maxHFEnergyPlus_33_42,  maxHFEnergyMinus_33_42);
    minOfMaxHF_33_43  = std::min(maxHFEnergyPlus_33_43,  maxHFEnergyMinus_33_43);
    minOfMaxHF_33_44  = std::min(maxHFEnergyPlus_33_44,  maxHFEnergyMinus_33_44);
    minOfMaxHF_33_45  = std::min(maxHFEnergyPlus_33_45,  maxHFEnergyMinus_33_45);
    minOfMaxHF_33_46  = std::min(maxHFEnergyPlus_33_46,  maxHFEnergyMinus_33_46);
    minOfMaxHF_33_47  = std::min(maxHFEnergyPlus_33_47,  maxHFEnergyMinus_33_47);
    minOfMaxHF_33_48  = std::min(maxHFEnergyPlus_33_48,  maxHFEnergyMinus_33_48);
    minOfMaxHF_33_49  = std::min(maxHFEnergyPlus_33_49,  maxHFEnergyMinus_33_49);
    minOfMaxHF_33_50  = std::min(maxHFEnergyPlus_33_50,  maxHFEnergyMinus_33_50);
    minOfMaxHF_33_51  = std::min(maxHFEnergyPlus_33_51,  maxHFEnergyMinus_33_51);
    minOfMaxHF_33_52  = std::min(maxHFEnergyPlus_33_52,  maxHFEnergyMinus_33_52);


    minOfMaxHF_34_35  = std::min(maxHFEnergyPlus_34_35,  maxHFEnergyMinus_34_35);
    minOfMaxHF_34_36  = std::min(maxHFEnergyPlus_34_36,  maxHFEnergyMinus_34_36);
    minOfMaxHF_34_37  = std::min(maxHFEnergyPlus_34_37,  maxHFEnergyMinus_34_37);
    minOfMaxHF_34_38  = std::min(maxHFEnergyPlus_34_38,  maxHFEnergyMinus_34_38);
    minOfMaxHF_34_39  = std::min(maxHFEnergyPlus_34_39,  maxHFEnergyMinus_34_39);
    minOfMaxHF_34_40  = std::min(maxHFEnergyPlus_34_40,  maxHFEnergyMinus_34_40);
    minOfMaxHF_34_41  = std::min(maxHFEnergyPlus_34_41,  maxHFEnergyMinus_34_41);
    minOfMaxHF_34_42  = std::min(maxHFEnergyPlus_34_42,  maxHFEnergyMinus_34_42);
    minOfMaxHF_34_43  = std::min(maxHFEnergyPlus_34_43,  maxHFEnergyMinus_34_43);
    minOfMaxHF_34_44  = std::min(maxHFEnergyPlus_34_44,  maxHFEnergyMinus_34_44);
    minOfMaxHF_34_45  = std::min(maxHFEnergyPlus_34_45,  maxHFEnergyMinus_34_45);
    minOfMaxHF_34_46  = std::min(maxHFEnergyPlus_34_46,  maxHFEnergyMinus_34_46);
    minOfMaxHF_34_47  = std::min(maxHFEnergyPlus_34_47,  maxHFEnergyMinus_34_47);
    minOfMaxHF_34_48  = std::min(maxHFEnergyPlus_34_48,  maxHFEnergyMinus_34_48);
    minOfMaxHF_34_49  = std::min(maxHFEnergyPlus_34_49,  maxHFEnergyMinus_34_49);
    minOfMaxHF_34_50  = std::min(maxHFEnergyPlus_34_50,  maxHFEnergyMinus_34_50);
    minOfMaxHF_34_51  = std::min(maxHFEnergyPlus_34_51,  maxHFEnergyMinus_34_51);
    minOfMaxHF_34_52  = std::min(maxHFEnergyPlus_34_52,  maxHFEnergyMinus_34_52);


    minOfMaxHF_35_36  = std::min(maxHFEnergyPlus_35_36,  maxHFEnergyMinus_35_36);
    minOfMaxHF_35_37  = std::min(maxHFEnergyPlus_35_37,  maxHFEnergyMinus_35_37);
    minOfMaxHF_35_38  = std::min(maxHFEnergyPlus_35_38,  maxHFEnergyMinus_35_38);
    minOfMaxHF_35_39  = std::min(maxHFEnergyPlus_35_39,  maxHFEnergyMinus_35_39);
    minOfMaxHF_35_40  = std::min(maxHFEnergyPlus_35_40,  maxHFEnergyMinus_35_40);
    minOfMaxHF_35_41  = std::min(maxHFEnergyPlus_35_41,  maxHFEnergyMinus_35_41);
    minOfMaxHF_35_42  = std::min(maxHFEnergyPlus_35_42,  maxHFEnergyMinus_35_42);
    minOfMaxHF_35_43  = std::min(maxHFEnergyPlus_35_43,  maxHFEnergyMinus_35_43);
    minOfMaxHF_35_44  = std::min(maxHFEnergyPlus_35_44,  maxHFEnergyMinus_35_44);
    minOfMaxHF_35_45  = std::min(maxHFEnergyPlus_35_45,  maxHFEnergyMinus_35_45);
    minOfMaxHF_35_46  = std::min(maxHFEnergyPlus_35_46,  maxHFEnergyMinus_35_46);
    minOfMaxHF_35_47  = std::min(maxHFEnergyPlus_35_47,  maxHFEnergyMinus_35_47);
    minOfMaxHF_35_48  = std::min(maxHFEnergyPlus_35_48,  maxHFEnergyMinus_35_48);
    minOfMaxHF_35_49  = std::min(maxHFEnergyPlus_35_49,  maxHFEnergyMinus_35_49);
    minOfMaxHF_35_50  = std::min(maxHFEnergyPlus_35_50,  maxHFEnergyMinus_35_50);
    minOfMaxHF_35_51  = std::min(maxHFEnergyPlus_35_51,  maxHFEnergyMinus_35_51);
    minOfMaxHF_35_52  = std::min(maxHFEnergyPlus_35_52,  maxHFEnergyMinus_35_52);


    minOfMaxHF_36_37  = std::min(maxHFEnergyPlus_36_37,  maxHFEnergyMinus_36_37);
    minOfMaxHF_36_38  = std::min(maxHFEnergyPlus_36_38,  maxHFEnergyMinus_36_38);
    minOfMaxHF_36_39  = std::min(maxHFEnergyPlus_36_39,  maxHFEnergyMinus_36_39);
    minOfMaxHF_36_40  = std::min(maxHFEnergyPlus_36_40,  maxHFEnergyMinus_36_40);
    minOfMaxHF_36_41  = std::min(maxHFEnergyPlus_36_41,  maxHFEnergyMinus_36_41);
    minOfMaxHF_36_42  = std::min(maxHFEnergyPlus_36_42,  maxHFEnergyMinus_36_42);
    minOfMaxHF_36_43  = std::min(maxHFEnergyPlus_36_43,  maxHFEnergyMinus_36_43);
    minOfMaxHF_36_44  = std::min(maxHFEnergyPlus_36_44,  maxHFEnergyMinus_36_44);
    minOfMaxHF_36_45  = std::min(maxHFEnergyPlus_36_45,  maxHFEnergyMinus_36_45);
    minOfMaxHF_36_46  = std::min(maxHFEnergyPlus_36_46,  maxHFEnergyMinus_36_46);
    minOfMaxHF_36_47  = std::min(maxHFEnergyPlus_36_47,  maxHFEnergyMinus_36_47);
    minOfMaxHF_36_48  = std::min(maxHFEnergyPlus_36_48,  maxHFEnergyMinus_36_48);
    minOfMaxHF_36_49  = std::min(maxHFEnergyPlus_36_49,  maxHFEnergyMinus_36_49);
    minOfMaxHF_36_50  = std::min(maxHFEnergyPlus_36_50,  maxHFEnergyMinus_36_50);
    minOfMaxHF_36_51  = std::min(maxHFEnergyPlus_36_51,  maxHFEnergyMinus_36_51);
    minOfMaxHF_36_52  = std::min(maxHFEnergyPlus_36_52,  maxHFEnergyMinus_36_52);


    minOfMaxHF_37_38  = std::min(maxHFEnergyPlus_37_38,  maxHFEnergyMinus_37_38);
    minOfMaxHF_37_39  = std::min(maxHFEnergyPlus_37_39,  maxHFEnergyMinus_37_39);
    minOfMaxHF_37_40  = std::min(maxHFEnergyPlus_37_40,  maxHFEnergyMinus_37_40);
    minOfMaxHF_37_41  = std::min(maxHFEnergyPlus_37_41,  maxHFEnergyMinus_37_41);
    minOfMaxHF_37_42  = std::min(maxHFEnergyPlus_37_42,  maxHFEnergyMinus_37_42);
    minOfMaxHF_37_43  = std::min(maxHFEnergyPlus_37_43,  maxHFEnergyMinus_37_43);
    minOfMaxHF_37_44  = std::min(maxHFEnergyPlus_37_44,  maxHFEnergyMinus_37_44);
    minOfMaxHF_37_45  = std::min(maxHFEnergyPlus_37_45,  maxHFEnergyMinus_37_45);
    minOfMaxHF_37_46  = std::min(maxHFEnergyPlus_37_46,  maxHFEnergyMinus_37_46);
    minOfMaxHF_37_47  = std::min(maxHFEnergyPlus_37_47,  maxHFEnergyMinus_37_47);
    minOfMaxHF_37_48  = std::min(maxHFEnergyPlus_37_48,  maxHFEnergyMinus_37_48);
    minOfMaxHF_37_49  = std::min(maxHFEnergyPlus_37_49,  maxHFEnergyMinus_37_49);
    minOfMaxHF_37_50  = std::min(maxHFEnergyPlus_37_50,  maxHFEnergyMinus_37_50);
    minOfMaxHF_37_51  = std::min(maxHFEnergyPlus_37_51,  maxHFEnergyMinus_37_51);
    minOfMaxHF_37_52  = std::min(maxHFEnergyPlus_37_52,  maxHFEnergyMinus_37_52);


    minOfMaxHF_38_39  = std::min(maxHFEnergyPlus_38_39,  maxHFEnergyMinus_38_39);
    minOfMaxHF_38_40  = std::min(maxHFEnergyPlus_38_40,  maxHFEnergyMinus_38_40);
    minOfMaxHF_38_41  = std::min(maxHFEnergyPlus_38_41,  maxHFEnergyMinus_38_41);
    minOfMaxHF_38_42  = std::min(maxHFEnergyPlus_38_42,  maxHFEnergyMinus_38_42);
    minOfMaxHF_38_43  = std::min(maxHFEnergyPlus_38_43,  maxHFEnergyMinus_38_43);
    minOfMaxHF_38_44  = std::min(maxHFEnergyPlus_38_44,  maxHFEnergyMinus_38_44);
    minOfMaxHF_38_45  = std::min(maxHFEnergyPlus_38_45,  maxHFEnergyMinus_38_45);
    minOfMaxHF_38_46  = std::min(maxHFEnergyPlus_38_46,  maxHFEnergyMinus_38_46);
    minOfMaxHF_38_47  = std::min(maxHFEnergyPlus_38_47,  maxHFEnergyMinus_38_47);
    minOfMaxHF_38_48  = std::min(maxHFEnergyPlus_38_48,  maxHFEnergyMinus_38_48);
    minOfMaxHF_38_49  = std::min(maxHFEnergyPlus_38_49,  maxHFEnergyMinus_38_49);
    minOfMaxHF_38_50  = std::min(maxHFEnergyPlus_38_50,  maxHFEnergyMinus_38_50);
    minOfMaxHF_38_51  = std::min(maxHFEnergyPlus_38_51,  maxHFEnergyMinus_38_51);
    minOfMaxHF_38_52  = std::min(maxHFEnergyPlus_38_52,  maxHFEnergyMinus_38_52);


    minOfMaxHF_39_40  = std::min(maxHFEnergyPlus_39_40,  maxHFEnergyMinus_39_40);
    minOfMaxHF_39_41  = std::min(maxHFEnergyPlus_39_41,  maxHFEnergyMinus_39_41);
    minOfMaxHF_39_42  = std::min(maxHFEnergyPlus_39_42,  maxHFEnergyMinus_39_42);
    minOfMaxHF_39_43  = std::min(maxHFEnergyPlus_39_43,  maxHFEnergyMinus_39_43);
    minOfMaxHF_39_44  = std::min(maxHFEnergyPlus_39_44,  maxHFEnergyMinus_39_44);
    minOfMaxHF_39_45  = std::min(maxHFEnergyPlus_39_45,  maxHFEnergyMinus_39_45);
    minOfMaxHF_39_46  = std::min(maxHFEnergyPlus_39_46,  maxHFEnergyMinus_39_46);
    minOfMaxHF_39_47  = std::min(maxHFEnergyPlus_39_47,  maxHFEnergyMinus_39_47);
    minOfMaxHF_39_48  = std::min(maxHFEnergyPlus_39_48,  maxHFEnergyMinus_39_48);
    minOfMaxHF_39_49  = std::min(maxHFEnergyPlus_39_49,  maxHFEnergyMinus_39_49);
    minOfMaxHF_39_50  = std::min(maxHFEnergyPlus_39_50,  maxHFEnergyMinus_39_50);
    minOfMaxHF_39_51  = std::min(maxHFEnergyPlus_39_51,  maxHFEnergyMinus_39_51);
    minOfMaxHF_39_52  = std::min(maxHFEnergyPlus_39_52,  maxHFEnergyMinus_39_52);


    minOfMaxHF_40_41  = std::min(maxHFEnergyPlus_40_41,  maxHFEnergyMinus_40_41);
    minOfMaxHF_40_42  = std::min(maxHFEnergyPlus_40_42,  maxHFEnergyMinus_40_42);
    minOfMaxHF_40_43  = std::min(maxHFEnergyPlus_40_43,  maxHFEnergyMinus_40_43);
    minOfMaxHF_40_44  = std::min(maxHFEnergyPlus_40_44,  maxHFEnergyMinus_40_44);
    minOfMaxHF_40_45  = std::min(maxHFEnergyPlus_40_45,  maxHFEnergyMinus_40_45);
    minOfMaxHF_40_46  = std::min(maxHFEnergyPlus_40_46,  maxHFEnergyMinus_40_46);
    minOfMaxHF_40_47  = std::min(maxHFEnergyPlus_40_47,  maxHFEnergyMinus_40_47);
    minOfMaxHF_40_48  = std::min(maxHFEnergyPlus_40_48,  maxHFEnergyMinus_40_48);
    minOfMaxHF_40_49  = std::min(maxHFEnergyPlus_40_49,  maxHFEnergyMinus_40_49);
    minOfMaxHF_40_50  = std::min(maxHFEnergyPlus_40_50,  maxHFEnergyMinus_40_50);
    minOfMaxHF_40_51  = std::min(maxHFEnergyPlus_40_51,  maxHFEnergyMinus_40_51);
    minOfMaxHF_40_52  = std::min(maxHFEnergyPlus_40_52,  maxHFEnergyMinus_40_52);


    minOfMaxHF_41_42  = std::min(maxHFEnergyPlus_41_42,  maxHFEnergyMinus_41_42);
    minOfMaxHF_41_43  = std::min(maxHFEnergyPlus_41_43,  maxHFEnergyMinus_41_43);
    minOfMaxHF_41_44  = std::min(maxHFEnergyPlus_41_44,  maxHFEnergyMinus_41_44);
    minOfMaxHF_41_45  = std::min(maxHFEnergyPlus_41_45,  maxHFEnergyMinus_41_45);
    minOfMaxHF_41_46  = std::min(maxHFEnergyPlus_41_46,  maxHFEnergyMinus_41_46);
    minOfMaxHF_41_47  = std::min(maxHFEnergyPlus_41_47,  maxHFEnergyMinus_41_47);
    minOfMaxHF_41_48  = std::min(maxHFEnergyPlus_41_48,  maxHFEnergyMinus_41_48);
    minOfMaxHF_41_49  = std::min(maxHFEnergyPlus_41_49,  maxHFEnergyMinus_41_49);
    minOfMaxHF_41_50  = std::min(maxHFEnergyPlus_41_50,  maxHFEnergyMinus_41_50);
    minOfMaxHF_41_51  = std::min(maxHFEnergyPlus_41_51,  maxHFEnergyMinus_41_51);
    minOfMaxHF_41_52  = std::min(maxHFEnergyPlus_41_52,  maxHFEnergyMinus_41_52);


    minOfMaxHF_42_43  = std::min(maxHFEnergyPlus_42_43,  maxHFEnergyMinus_42_43);
    minOfMaxHF_42_44  = std::min(maxHFEnergyPlus_42_44,  maxHFEnergyMinus_42_44);
    minOfMaxHF_42_45  = std::min(maxHFEnergyPlus_42_45,  maxHFEnergyMinus_42_45);
    minOfMaxHF_42_46  = std::min(maxHFEnergyPlus_42_46,  maxHFEnergyMinus_42_46);
    minOfMaxHF_42_47  = std::min(maxHFEnergyPlus_42_47,  maxHFEnergyMinus_42_47);
    minOfMaxHF_42_48  = std::min(maxHFEnergyPlus_42_48,  maxHFEnergyMinus_42_48);
    minOfMaxHF_42_49  = std::min(maxHFEnergyPlus_42_49,  maxHFEnergyMinus_42_49);
    minOfMaxHF_42_50  = std::min(maxHFEnergyPlus_42_50,  maxHFEnergyMinus_42_50);
    minOfMaxHF_42_51  = std::min(maxHFEnergyPlus_42_51,  maxHFEnergyMinus_42_51);
    minOfMaxHF_42_52  = std::min(maxHFEnergyPlus_42_52,  maxHFEnergyMinus_42_52);


    minOfMaxHF_43_44  = std::min(maxHFEnergyPlus_43_44,  maxHFEnergyMinus_43_44);
    minOfMaxHF_43_45  = std::min(maxHFEnergyPlus_43_45,  maxHFEnergyMinus_43_45);
    minOfMaxHF_43_46  = std::min(maxHFEnergyPlus_43_46,  maxHFEnergyMinus_43_46);
    minOfMaxHF_43_47  = std::min(maxHFEnergyPlus_43_47,  maxHFEnergyMinus_43_47);
    minOfMaxHF_43_48  = std::min(maxHFEnergyPlus_43_48,  maxHFEnergyMinus_43_48);
    minOfMaxHF_43_49  = std::min(maxHFEnergyPlus_43_49,  maxHFEnergyMinus_43_49);
    minOfMaxHF_43_50  = std::min(maxHFEnergyPlus_43_50,  maxHFEnergyMinus_43_50);
    minOfMaxHF_43_51  = std::min(maxHFEnergyPlus_43_51,  maxHFEnergyMinus_43_51);
    minOfMaxHF_43_52  = std::min(maxHFEnergyPlus_43_52,  maxHFEnergyMinus_43_52);


    minOfMaxHF_44_45  = std::min(maxHFEnergyPlus_44_45,  maxHFEnergyMinus_44_45);
    minOfMaxHF_44_46  = std::min(maxHFEnergyPlus_44_46,  maxHFEnergyMinus_44_46);
    minOfMaxHF_44_47  = std::min(maxHFEnergyPlus_44_47,  maxHFEnergyMinus_44_47);
    minOfMaxHF_44_48  = std::min(maxHFEnergyPlus_44_48,  maxHFEnergyMinus_44_48);
    minOfMaxHF_44_49  = std::min(maxHFEnergyPlus_44_49,  maxHFEnergyMinus_44_49);
    minOfMaxHF_44_50  = std::min(maxHFEnergyPlus_44_50,  maxHFEnergyMinus_44_50);
    minOfMaxHF_44_51  = std::min(maxHFEnergyPlus_44_51,  maxHFEnergyMinus_44_51);
    minOfMaxHF_44_52  = std::min(maxHFEnergyPlus_44_52,  maxHFEnergyMinus_44_52);


    minOfMaxHF_45_46  = std::min(maxHFEnergyPlus_45_46,  maxHFEnergyMinus_45_46);
    minOfMaxHF_45_47  = std::min(maxHFEnergyPlus_45_47,  maxHFEnergyMinus_45_47);
    minOfMaxHF_45_48  = std::min(maxHFEnergyPlus_45_48,  maxHFEnergyMinus_45_48);
    minOfMaxHF_45_49  = std::min(maxHFEnergyPlus_45_49,  maxHFEnergyMinus_45_49);
    minOfMaxHF_45_50  = std::min(maxHFEnergyPlus_45_50,  maxHFEnergyMinus_45_50);
    minOfMaxHF_45_51  = std::min(maxHFEnergyPlus_45_51,  maxHFEnergyMinus_45_51);
    minOfMaxHF_45_52  = std::min(maxHFEnergyPlus_45_52,  maxHFEnergyMinus_45_52);


    minOfMaxHF_46_47  = std::min(maxHFEnergyPlus_46_47,  maxHFEnergyMinus_46_47);
    minOfMaxHF_46_48  = std::min(maxHFEnergyPlus_46_48,  maxHFEnergyMinus_46_48);
    minOfMaxHF_46_49  = std::min(maxHFEnergyPlus_46_49,  maxHFEnergyMinus_46_49);
    minOfMaxHF_46_50  = std::min(maxHFEnergyPlus_46_50,  maxHFEnergyMinus_46_50);
    minOfMaxHF_46_51  = std::min(maxHFEnergyPlus_46_51,  maxHFEnergyMinus_46_51);
    minOfMaxHF_46_52  = std::min(maxHFEnergyPlus_46_52,  maxHFEnergyMinus_46_52);


    minOfMaxHF_47_48  = std::min(maxHFEnergyPlus_47_48,  maxHFEnergyMinus_47_48);
    minOfMaxHF_47_49  = std::min(maxHFEnergyPlus_47_49,  maxHFEnergyMinus_47_49);
    minOfMaxHF_47_50  = std::min(maxHFEnergyPlus_47_50,  maxHFEnergyMinus_47_50);
    minOfMaxHF_47_51  = std::min(maxHFEnergyPlus_47_51,  maxHFEnergyMinus_47_51);
    minOfMaxHF_47_52  = std::min(maxHFEnergyPlus_47_52,  maxHFEnergyMinus_47_52);


    minOfMaxHF_48_49  = std::min(maxHFEnergyPlus_48_49,  maxHFEnergyMinus_48_49);
    minOfMaxHF_48_50  = std::min(maxHFEnergyPlus_48_50,  maxHFEnergyMinus_48_50);
    minOfMaxHF_48_51  = std::min(maxHFEnergyPlus_48_51,  maxHFEnergyMinus_48_51);
    minOfMaxHF_48_52  = std::min(maxHFEnergyPlus_48_52,  maxHFEnergyMinus_48_52);


    minOfMaxHF_49_50  = std::min(maxHFEnergyPlus_49_50,  maxHFEnergyMinus_49_50);
    minOfMaxHF_49_51  = std::min(maxHFEnergyPlus_49_51,  maxHFEnergyMinus_49_51);
    minOfMaxHF_49_52  = std::min(maxHFEnergyPlus_49_52,  maxHFEnergyMinus_49_52);


    minOfMaxHF_50_51  = std::min(maxHFEnergyPlus_50_51,  maxHFEnergyMinus_50_51);
    minOfMaxHF_50_52  = std::min(maxHFEnergyPlus_50_52,  maxHFEnergyMinus_50_52);


    minOfMaxHF_51_52  = std::min(maxHFEnergyPlus_51_52,  maxHFEnergyMinus_51_52);




    
    // min of 
    minOfHFSums_30_31 = std::min(hfSumPlus_30_31, hfSumMinus_30_31);
    minOfHFSums_30_32 = std::min(hfSumPlus_30_32, hfSumMinus_30_32);
    minOfHFSums_30_33 = std::min(hfSumPlus_30_33, hfSumMinus_30_33);
    minOfHFSums_30_34 = std::min(hfSumPlus_30_34, hfSumMinus_30_34);
    minOfHFSums_30_35 = std::min(hfSumPlus_30_35, hfSumMinus_30_35);
    minOfHFSums_30_36 = std::min(hfSumPlus_30_36, hfSumMinus_30_36);
    minOfHFSums_30_37 = std::min(hfSumPlus_30_37, hfSumMinus_30_37);
    minOfHFSums_30_38 = std::min(hfSumPlus_30_38, hfSumMinus_30_38);
    minOfHFSums_30_39 = std::min(hfSumPlus_30_39, hfSumMinus_30_39);
    minOfHFSums_30_40 = std::min(hfSumPlus_30_40, hfSumMinus_30_40);
    minOfHFSums_30_41 = std::min(hfSumPlus_30_41, hfSumMinus_30_41);
    minOfHFSums_30_42 = std::min(hfSumPlus_30_42, hfSumMinus_30_42);
    minOfHFSums_30_43 = std::min(hfSumPlus_30_43, hfSumMinus_30_43);
    minOfHFSums_30_44 = std::min(hfSumPlus_30_44, hfSumMinus_30_44);
    minOfHFSums_30_45 = std::min(hfSumPlus_30_45, hfSumMinus_30_45);
    minOfHFSums_30_46 = std::min(hfSumPlus_30_46, hfSumMinus_30_46);
    minOfHFSums_30_47 = std::min(hfSumPlus_30_47, hfSumMinus_30_47);
    minOfHFSums_30_48 = std::min(hfSumPlus_30_48, hfSumMinus_30_48);
    minOfHFSums_30_49 = std::min(hfSumPlus_30_49, hfSumMinus_30_49);
    minOfHFSums_30_50 = std::min(hfSumPlus_30_50, hfSumMinus_30_50);
    minOfHFSums_30_51 = std::min(hfSumPlus_30_51, hfSumMinus_30_51);
    minOfHFSums_30_52 = std::min(hfSumPlus_30_52, hfSumMinus_30_52);

    minOfHFSums_31_32 = std::min(hfSumPlus_31_32, hfSumMinus_31_32);
    minOfHFSums_31_33 = std::min(hfSumPlus_31_33, hfSumMinus_31_33);
    minOfHFSums_31_34 = std::min(hfSumPlus_31_34, hfSumMinus_31_34);
    minOfHFSums_31_35 = std::min(hfSumPlus_31_35, hfSumMinus_31_35);
    minOfHFSums_31_36 = std::min(hfSumPlus_31_36, hfSumMinus_31_36);
    minOfHFSums_31_37 = std::min(hfSumPlus_31_37, hfSumMinus_31_37);
    minOfHFSums_31_38 = std::min(hfSumPlus_31_38, hfSumMinus_31_38);
    minOfHFSums_31_39 = std::min(hfSumPlus_31_39, hfSumMinus_31_39);
    minOfHFSums_31_40 = std::min(hfSumPlus_31_40, hfSumMinus_31_40);
    minOfHFSums_31_41 = std::min(hfSumPlus_31_41, hfSumMinus_31_41);
    minOfHFSums_31_42 = std::min(hfSumPlus_31_42, hfSumMinus_31_42);
    minOfHFSums_31_43 = std::min(hfSumPlus_31_43, hfSumMinus_31_43);
    minOfHFSums_31_44 = std::min(hfSumPlus_31_44, hfSumMinus_31_44);
    minOfHFSums_31_45 = std::min(hfSumPlus_31_45, hfSumMinus_31_45);
    minOfHFSums_31_46 = std::min(hfSumPlus_31_46, hfSumMinus_31_46);
    minOfHFSums_31_47 = std::min(hfSumPlus_31_47, hfSumMinus_31_47);
    minOfHFSums_31_48 = std::min(hfSumPlus_31_48, hfSumMinus_31_48);
    minOfHFSums_31_49 = std::min(hfSumPlus_31_49, hfSumMinus_31_49);
    minOfHFSums_31_50 = std::min(hfSumPlus_31_50, hfSumMinus_31_50);
    minOfHFSums_31_51 = std::min(hfSumPlus_31_51, hfSumMinus_31_51);
    minOfHFSums_31_52 = std::min(hfSumPlus_31_52, hfSumMinus_31_52);

    minOfHFSums_32_33 = std::min(hfSumPlus_32_33, hfSumMinus_32_33);
    minOfHFSums_32_34 = std::min(hfSumPlus_32_34, hfSumMinus_32_34);
    minOfHFSums_32_35 = std::min(hfSumPlus_32_35, hfSumMinus_32_35);
    minOfHFSums_32_36 = std::min(hfSumPlus_32_36, hfSumMinus_32_36);
    minOfHFSums_32_37 = std::min(hfSumPlus_32_37, hfSumMinus_32_37);
    minOfHFSums_32_38 = std::min(hfSumPlus_32_38, hfSumMinus_32_38);
    minOfHFSums_32_39 = std::min(hfSumPlus_32_39, hfSumMinus_32_39);
    minOfHFSums_32_40 = std::min(hfSumPlus_32_40, hfSumMinus_32_40);
    minOfHFSums_32_41 = std::min(hfSumPlus_32_41, hfSumMinus_32_41);
    minOfHFSums_32_42 = std::min(hfSumPlus_32_42, hfSumMinus_32_42);
    minOfHFSums_32_43 = std::min(hfSumPlus_32_43, hfSumMinus_32_43);
    minOfHFSums_32_44 = std::min(hfSumPlus_32_44, hfSumMinus_32_44);
    minOfHFSums_32_45 = std::min(hfSumPlus_32_45, hfSumMinus_32_45);
    minOfHFSums_32_46 = std::min(hfSumPlus_32_46, hfSumMinus_32_46);
    minOfHFSums_32_47 = std::min(hfSumPlus_32_47, hfSumMinus_32_47);
    minOfHFSums_32_48 = std::min(hfSumPlus_32_48, hfSumMinus_32_48);
    minOfHFSums_32_49 = std::min(hfSumPlus_32_49, hfSumMinus_32_49);
    minOfHFSums_32_50 = std::min(hfSumPlus_32_50, hfSumMinus_32_50);
    minOfHFSums_32_51 = std::min(hfSumPlus_32_51, hfSumMinus_32_51);
    minOfHFSums_32_52 = std::min(hfSumPlus_32_52, hfSumMinus_32_52);

    minOfHFSums_33_34 = std::min(hfSumPlus_33_34, hfSumMinus_33_34);
    minOfHFSums_33_35 = std::min(hfSumPlus_33_35, hfSumMinus_33_35);
    minOfHFSums_33_36 = std::min(hfSumPlus_33_36, hfSumMinus_33_36);
    minOfHFSums_33_37 = std::min(hfSumPlus_33_37, hfSumMinus_33_37);
    minOfHFSums_33_38 = std::min(hfSumPlus_33_38, hfSumMinus_33_38);
    minOfHFSums_33_39 = std::min(hfSumPlus_33_39, hfSumMinus_33_39);
    minOfHFSums_33_40 = std::min(hfSumPlus_33_40, hfSumMinus_33_40);
    minOfHFSums_33_41 = std::min(hfSumPlus_33_41, hfSumMinus_33_41);
    minOfHFSums_33_42 = std::min(hfSumPlus_33_42, hfSumMinus_33_42);
    minOfHFSums_33_43 = std::min(hfSumPlus_33_43, hfSumMinus_33_43);
    minOfHFSums_33_44 = std::min(hfSumPlus_33_44, hfSumMinus_33_44);
    minOfHFSums_33_45 = std::min(hfSumPlus_33_45, hfSumMinus_33_45);
    minOfHFSums_33_46 = std::min(hfSumPlus_33_46, hfSumMinus_33_46);
    minOfHFSums_33_47 = std::min(hfSumPlus_33_47, hfSumMinus_33_47);
    minOfHFSums_33_48 = std::min(hfSumPlus_33_48, hfSumMinus_33_48);
    minOfHFSums_33_49 = std::min(hfSumPlus_33_49, hfSumMinus_33_49);
    minOfHFSums_33_50 = std::min(hfSumPlus_33_50, hfSumMinus_33_50);
    minOfHFSums_33_51 = std::min(hfSumPlus_33_51, hfSumMinus_33_51);
    minOfHFSums_33_52 = std::min(hfSumPlus_33_52, hfSumMinus_33_52);

    minOfHFSums_34_35 = std::min(hfSumPlus_34_35, hfSumMinus_34_35);
    minOfHFSums_34_36 = std::min(hfSumPlus_34_36, hfSumMinus_34_36);
    minOfHFSums_34_37 = std::min(hfSumPlus_34_37, hfSumMinus_34_37);
    minOfHFSums_34_38 = std::min(hfSumPlus_34_38, hfSumMinus_34_38);
    minOfHFSums_34_39 = std::min(hfSumPlus_34_39, hfSumMinus_34_39);
    minOfHFSums_34_40 = std::min(hfSumPlus_34_40, hfSumMinus_34_40);
    minOfHFSums_34_41 = std::min(hfSumPlus_34_41, hfSumMinus_34_41);
    minOfHFSums_34_42 = std::min(hfSumPlus_34_42, hfSumMinus_34_42);
    minOfHFSums_34_43 = std::min(hfSumPlus_34_43, hfSumMinus_34_43);
    minOfHFSums_34_44 = std::min(hfSumPlus_34_44, hfSumMinus_34_44);
    minOfHFSums_34_45 = std::min(hfSumPlus_34_45, hfSumMinus_34_45);
    minOfHFSums_34_46 = std::min(hfSumPlus_34_46, hfSumMinus_34_46);
    minOfHFSums_34_47 = std::min(hfSumPlus_34_47, hfSumMinus_34_47);
    minOfHFSums_34_48 = std::min(hfSumPlus_34_48, hfSumMinus_34_48);
    minOfHFSums_34_49 = std::min(hfSumPlus_34_49, hfSumMinus_34_49);
    minOfHFSums_34_50 = std::min(hfSumPlus_34_50, hfSumMinus_34_50);
    minOfHFSums_34_51 = std::min(hfSumPlus_34_51, hfSumMinus_34_51);
    minOfHFSums_34_52 = std::min(hfSumPlus_34_52, hfSumMinus_34_52);

    minOfHFSums_35_36 = std::min(hfSumPlus_35_36, hfSumMinus_35_36);
    minOfHFSums_35_37 = std::min(hfSumPlus_35_37, hfSumMinus_35_37);
    minOfHFSums_35_38 = std::min(hfSumPlus_35_38, hfSumMinus_35_38);
    minOfHFSums_35_39 = std::min(hfSumPlus_35_39, hfSumMinus_35_39);
    minOfHFSums_35_40 = std::min(hfSumPlus_35_40, hfSumMinus_35_40);
    minOfHFSums_35_41 = std::min(hfSumPlus_35_41, hfSumMinus_35_41);
    minOfHFSums_35_42 = std::min(hfSumPlus_35_42, hfSumMinus_35_42);
    minOfHFSums_35_43 = std::min(hfSumPlus_35_43, hfSumMinus_35_43);
    minOfHFSums_35_44 = std::min(hfSumPlus_35_44, hfSumMinus_35_44);
    minOfHFSums_35_45 = std::min(hfSumPlus_35_45, hfSumMinus_35_45);
    minOfHFSums_35_46 = std::min(hfSumPlus_35_46, hfSumMinus_35_46);
    minOfHFSums_35_47 = std::min(hfSumPlus_35_47, hfSumMinus_35_47);
    minOfHFSums_35_48 = std::min(hfSumPlus_35_48, hfSumMinus_35_48);
    minOfHFSums_35_49 = std::min(hfSumPlus_35_49, hfSumMinus_35_49);
    minOfHFSums_35_50 = std::min(hfSumPlus_35_50, hfSumMinus_35_50);
    minOfHFSums_35_51 = std::min(hfSumPlus_35_51, hfSumMinus_35_51);
    minOfHFSums_35_52 = std::min(hfSumPlus_35_52, hfSumMinus_35_52);

    minOfHFSums_36_37 = std::min(hfSumPlus_36_37, hfSumMinus_36_37);
    minOfHFSums_36_38 = std::min(hfSumPlus_36_38, hfSumMinus_36_38);
    minOfHFSums_36_39 = std::min(hfSumPlus_36_39, hfSumMinus_36_39);
    minOfHFSums_36_40 = std::min(hfSumPlus_36_40, hfSumMinus_36_40);
    minOfHFSums_36_41 = std::min(hfSumPlus_36_41, hfSumMinus_36_41);
    minOfHFSums_36_42 = std::min(hfSumPlus_36_42, hfSumMinus_36_42);
    minOfHFSums_36_43 = std::min(hfSumPlus_36_43, hfSumMinus_36_43);
    minOfHFSums_36_44 = std::min(hfSumPlus_36_44, hfSumMinus_36_44);
    minOfHFSums_36_45 = std::min(hfSumPlus_36_45, hfSumMinus_36_45);
    minOfHFSums_36_46 = std::min(hfSumPlus_36_46, hfSumMinus_36_46);
    minOfHFSums_36_47 = std::min(hfSumPlus_36_47, hfSumMinus_36_47);
    minOfHFSums_36_48 = std::min(hfSumPlus_36_48, hfSumMinus_36_48);
    minOfHFSums_36_49 = std::min(hfSumPlus_36_49, hfSumMinus_36_49);
    minOfHFSums_36_50 = std::min(hfSumPlus_36_50, hfSumMinus_36_50);
    minOfHFSums_36_51 = std::min(hfSumPlus_36_51, hfSumMinus_36_51);
    minOfHFSums_36_52 = std::min(hfSumPlus_36_52, hfSumMinus_36_52);

    minOfHFSums_37_38 = std::min(hfSumPlus_37_38, hfSumMinus_37_38);
    minOfHFSums_37_39 = std::min(hfSumPlus_37_39, hfSumMinus_37_39);
    minOfHFSums_37_40 = std::min(hfSumPlus_37_40, hfSumMinus_37_40);
    minOfHFSums_37_41 = std::min(hfSumPlus_37_41, hfSumMinus_37_41);
    minOfHFSums_37_42 = std::min(hfSumPlus_37_42, hfSumMinus_37_42);
    minOfHFSums_37_43 = std::min(hfSumPlus_37_43, hfSumMinus_37_43);
    minOfHFSums_37_44 = std::min(hfSumPlus_37_44, hfSumMinus_37_44);
    minOfHFSums_37_45 = std::min(hfSumPlus_37_45, hfSumMinus_37_45);
    minOfHFSums_37_46 = std::min(hfSumPlus_37_46, hfSumMinus_37_46);
    minOfHFSums_37_47 = std::min(hfSumPlus_37_47, hfSumMinus_37_47);
    minOfHFSums_37_48 = std::min(hfSumPlus_37_48, hfSumMinus_37_48);
    minOfHFSums_37_49 = std::min(hfSumPlus_37_49, hfSumMinus_37_49);
    minOfHFSums_37_50 = std::min(hfSumPlus_37_50, hfSumMinus_37_50);
    minOfHFSums_37_51 = std::min(hfSumPlus_37_51, hfSumMinus_37_51);
    minOfHFSums_37_52 = std::min(hfSumPlus_37_52, hfSumMinus_37_52);

    minOfHFSums_38_39 = std::min(hfSumPlus_38_39, hfSumMinus_38_39);
    minOfHFSums_38_40 = std::min(hfSumPlus_38_40, hfSumMinus_38_40);
    minOfHFSums_38_41 = std::min(hfSumPlus_38_41, hfSumMinus_38_41);
    minOfHFSums_38_42 = std::min(hfSumPlus_38_42, hfSumMinus_38_42);
    minOfHFSums_38_43 = std::min(hfSumPlus_38_43, hfSumMinus_38_43);
    minOfHFSums_38_44 = std::min(hfSumPlus_38_44, hfSumMinus_38_44);
    minOfHFSums_38_45 = std::min(hfSumPlus_38_45, hfSumMinus_38_45);
    minOfHFSums_38_46 = std::min(hfSumPlus_38_46, hfSumMinus_38_46);
    minOfHFSums_38_47 = std::min(hfSumPlus_38_47, hfSumMinus_38_47);
    minOfHFSums_38_48 = std::min(hfSumPlus_38_48, hfSumMinus_38_48);
    minOfHFSums_38_49 = std::min(hfSumPlus_38_49, hfSumMinus_38_49);
    minOfHFSums_38_50 = std::min(hfSumPlus_38_50, hfSumMinus_38_50);
    minOfHFSums_38_51 = std::min(hfSumPlus_38_51, hfSumMinus_38_51);
    minOfHFSums_38_52 = std::min(hfSumPlus_38_52, hfSumMinus_38_52);

    minOfHFSums_39_40 = std::min(hfSumPlus_39_40, hfSumMinus_39_40);
    minOfHFSums_39_41 = std::min(hfSumPlus_39_41, hfSumMinus_39_41);
    minOfHFSums_39_42 = std::min(hfSumPlus_39_42, hfSumMinus_39_42);
    minOfHFSums_39_43 = std::min(hfSumPlus_39_43, hfSumMinus_39_43);
    minOfHFSums_39_44 = std::min(hfSumPlus_39_44, hfSumMinus_39_44);
    minOfHFSums_39_45 = std::min(hfSumPlus_39_45, hfSumMinus_39_45);
    minOfHFSums_39_46 = std::min(hfSumPlus_39_46, hfSumMinus_39_46);
    minOfHFSums_39_47 = std::min(hfSumPlus_39_47, hfSumMinus_39_47);
    minOfHFSums_39_48 = std::min(hfSumPlus_39_48, hfSumMinus_39_48);
    minOfHFSums_39_49 = std::min(hfSumPlus_39_49, hfSumMinus_39_49);
    minOfHFSums_39_50 = std::min(hfSumPlus_39_50, hfSumMinus_39_50);
    minOfHFSums_39_51 = std::min(hfSumPlus_39_51, hfSumMinus_39_51);
    minOfHFSums_39_52 = std::min(hfSumPlus_39_52, hfSumMinus_39_52);

    minOfHFSums_40_41 = std::min(hfSumPlus_40_41, hfSumMinus_40_41);
    minOfHFSums_40_42 = std::min(hfSumPlus_40_42, hfSumMinus_40_42);
    minOfHFSums_40_43 = std::min(hfSumPlus_40_43, hfSumMinus_40_43);
    minOfHFSums_40_44 = std::min(hfSumPlus_40_44, hfSumMinus_40_44);
    minOfHFSums_40_45 = std::min(hfSumPlus_40_45, hfSumMinus_40_45);
    minOfHFSums_40_46 = std::min(hfSumPlus_40_46, hfSumMinus_40_46);
    minOfHFSums_40_47 = std::min(hfSumPlus_40_47, hfSumMinus_40_47);
    minOfHFSums_40_48 = std::min(hfSumPlus_40_48, hfSumMinus_40_48);
    minOfHFSums_40_49 = std::min(hfSumPlus_40_49, hfSumMinus_40_49);
    minOfHFSums_40_50 = std::min(hfSumPlus_40_50, hfSumMinus_40_50);
    minOfHFSums_40_51 = std::min(hfSumPlus_40_51, hfSumMinus_40_51);
    minOfHFSums_40_52 = std::min(hfSumPlus_40_52, hfSumMinus_40_52);

    minOfHFSums_41_42 = std::min(hfSumPlus_41_42, hfSumMinus_41_42);
    minOfHFSums_41_43 = std::min(hfSumPlus_41_43, hfSumMinus_41_43);
    minOfHFSums_41_44 = std::min(hfSumPlus_41_44, hfSumMinus_41_44);
    minOfHFSums_41_45 = std::min(hfSumPlus_41_45, hfSumMinus_41_45);
    minOfHFSums_41_46 = std::min(hfSumPlus_41_46, hfSumMinus_41_46);
    minOfHFSums_41_47 = std::min(hfSumPlus_41_47, hfSumMinus_41_47);
    minOfHFSums_41_48 = std::min(hfSumPlus_41_48, hfSumMinus_41_48);
    minOfHFSums_41_49 = std::min(hfSumPlus_41_49, hfSumMinus_41_49);
    minOfHFSums_41_50 = std::min(hfSumPlus_41_50, hfSumMinus_41_50);
    minOfHFSums_41_51 = std::min(hfSumPlus_41_51, hfSumMinus_41_51);
    minOfHFSums_41_52 = std::min(hfSumPlus_41_52, hfSumMinus_41_52);

    minOfHFSums_42_43 = std::min(hfSumPlus_42_43, hfSumMinus_42_43);
    minOfHFSums_42_44 = std::min(hfSumPlus_42_44, hfSumMinus_42_44);
    minOfHFSums_42_45 = std::min(hfSumPlus_42_45, hfSumMinus_42_45);
    minOfHFSums_42_46 = std::min(hfSumPlus_42_46, hfSumMinus_42_46);
    minOfHFSums_42_47 = std::min(hfSumPlus_42_47, hfSumMinus_42_47);
    minOfHFSums_42_48 = std::min(hfSumPlus_42_48, hfSumMinus_42_48);
    minOfHFSums_42_49 = std::min(hfSumPlus_42_49, hfSumMinus_42_49);
    minOfHFSums_42_50 = std::min(hfSumPlus_42_50, hfSumMinus_42_50);
    minOfHFSums_42_51 = std::min(hfSumPlus_42_51, hfSumMinus_42_51);
    minOfHFSums_42_52 = std::min(hfSumPlus_42_52, hfSumMinus_42_52);

    minOfHFSums_43_44 = std::min(hfSumPlus_43_44, hfSumMinus_43_44);
    minOfHFSums_43_45 = std::min(hfSumPlus_43_45, hfSumMinus_43_45);
    minOfHFSums_43_46 = std::min(hfSumPlus_43_46, hfSumMinus_43_46);
    minOfHFSums_43_47 = std::min(hfSumPlus_43_47, hfSumMinus_43_47);
    minOfHFSums_43_48 = std::min(hfSumPlus_43_48, hfSumMinus_43_48);
    minOfHFSums_43_49 = std::min(hfSumPlus_43_49, hfSumMinus_43_49);
    minOfHFSums_43_50 = std::min(hfSumPlus_43_50, hfSumMinus_43_50);
    minOfHFSums_43_51 = std::min(hfSumPlus_43_51, hfSumMinus_43_51);
    minOfHFSums_43_52 = std::min(hfSumPlus_43_52, hfSumMinus_43_52);

    minOfHFSums_44_45 = std::min(hfSumPlus_44_45, hfSumMinus_44_45);
    minOfHFSums_44_46 = std::min(hfSumPlus_44_46, hfSumMinus_44_46);
    minOfHFSums_44_47 = std::min(hfSumPlus_44_47, hfSumMinus_44_47);
    minOfHFSums_44_48 = std::min(hfSumPlus_44_48, hfSumMinus_44_48);
    minOfHFSums_44_49 = std::min(hfSumPlus_44_49, hfSumMinus_44_49);
    minOfHFSums_44_50 = std::min(hfSumPlus_44_50, hfSumMinus_44_50);
    minOfHFSums_44_51 = std::min(hfSumPlus_44_51, hfSumMinus_44_51);
    minOfHFSums_44_52 = std::min(hfSumPlus_44_52, hfSumMinus_44_52);

    minOfHFSums_45_46 = std::min(hfSumPlus_45_46, hfSumMinus_45_46);
    minOfHFSums_45_47 = std::min(hfSumPlus_45_47, hfSumMinus_45_47);
    minOfHFSums_45_48 = std::min(hfSumPlus_45_48, hfSumMinus_45_48);
    minOfHFSums_45_49 = std::min(hfSumPlus_45_49, hfSumMinus_45_49);
    minOfHFSums_45_50 = std::min(hfSumPlus_45_50, hfSumMinus_45_50);
    minOfHFSums_45_51 = std::min(hfSumPlus_45_51, hfSumMinus_45_51);
    minOfHFSums_45_52 = std::min(hfSumPlus_45_52, hfSumMinus_45_52);

    minOfHFSums_46_47 = std::min(hfSumPlus_46_47, hfSumMinus_46_47);
    minOfHFSums_46_48 = std::min(hfSumPlus_46_48, hfSumMinus_46_48);
    minOfHFSums_46_49 = std::min(hfSumPlus_46_49, hfSumMinus_46_49);
    minOfHFSums_46_50 = std::min(hfSumPlus_46_50, hfSumMinus_46_50);
    minOfHFSums_46_51 = std::min(hfSumPlus_46_51, hfSumMinus_46_51);
    minOfHFSums_46_52 = std::min(hfSumPlus_46_52, hfSumMinus_46_52);

    minOfHFSums_47_48 = std::min(hfSumPlus_47_48, hfSumMinus_47_48);
    minOfHFSums_47_49 = std::min(hfSumPlus_47_49, hfSumMinus_47_49);
    minOfHFSums_47_50 = std::min(hfSumPlus_47_50, hfSumMinus_47_50);
    minOfHFSums_47_51 = std::min(hfSumPlus_47_51, hfSumMinus_47_51);
    minOfHFSums_47_52 = std::min(hfSumPlus_47_52, hfSumMinus_47_52);

    minOfHFSums_48_49 = std::min(hfSumPlus_48_49, hfSumMinus_48_49);
    minOfHFSums_48_50 = std::min(hfSumPlus_48_50, hfSumMinus_48_50);
    minOfHFSums_48_51 = std::min(hfSumPlus_48_51, hfSumMinus_48_51);
    minOfHFSums_48_52 = std::min(hfSumPlus_48_52, hfSumMinus_48_52);

    minOfHFSums_49_50 = std::min(hfSumPlus_49_50, hfSumMinus_49_50);
    minOfHFSums_49_51 = std::min(hfSumPlus_49_51, hfSumMinus_49_51);
    minOfHFSums_49_52 = std::min(hfSumPlus_49_52, hfSumMinus_49_52);

    minOfHFSums_50_51 = std::min(hfSumPlus_50_51, hfSumMinus_50_51);
    minOfHFSums_50_52 = std::min(hfSumPlus_50_52, hfSumMinus_50_52);

    minOfHFSums_51_52 = std::min(hfSumPlus_51_52, hfSumMinus_51_52);
   // min of --end

      

    // --- 4) Leptons ---
    Int_t nLeptons=0;
    Int_t nEle = branchElectron->GetEntriesFast();
    for (Int_t e=0; e<nEle; ++e) {
      auto* el = (Electron*)branchElectron->At(e);
      if (el->PT>15.0 && std::fabs(el->Eta)<2.4) nLeptons++;
    }
    Int_t nMu = branchMuon->GetEntriesFast();
    for (Int_t m=0; m<nMu; ++m) {
      auto* mu = (Muon*)branchMuon->At(m);
      if (mu->PT>10.0 && std::fabs(mu->Eta)<2.4) nLeptons++;
    }

    // --- 5) Missing ET ---
    auto* metobj = (MissingET*)branchMissingET->At(0);
    METval = metobj->MET;

    // --- 6) Rapidity gaps ---
    rapidityGapPlus  = 1e6;
    rapidityGapMinus = 1e6;
    for (double eta: trackEtas) {
      double dP = std::fabs(eta - 2.4);
      double dM = std::fabs(eta + 2.4);
      if (dP < rapidityGapPlus)  rapidityGapPlus  = dP;
      if (dM < rapidityGapMinus) rapidityGapMinus = dM;
    }
    if (trackEtas.empty()) {
      rapidityGapPlus=-1.0; rapidityGapMinus=-1.0;
    }

    // --- 7) Define and apply cuts ---
    bool c0  = (rapidityGapPlus >  0.0f);
    bool c1  = (rapidityGapMinus>  0.0f);
    bool c2  = (hfSumPlus      >= 0.0f);
    bool c3  = (hfSumMinus     >= 0.0f);
    bool c4  = (maxHFEnergy    >= 0.0f && maxHFEnergy < 100000.0f);
    bool c5  = (minHFEnergy    >= 0.0f && minHFEnergy < 100000.0f);
    bool c6  = (nJetsSel       >= 2);
    bool c7  = (nBtagSel       >= 1);
    bool c8  = (nLightSel      >= 1);
    bool c9  = (nLeptons       == 1);
    bool c10 = (METval         > 20.0f);
    bool c11 = (minOfMaxHF        >= 0.0f);
    bool c12 = (minOfHFSums       >= 0.0f);

    if (!c0 ) { ++failCount[0];  continue; }
    if (!c1 ) { ++failCount[1];  continue; }
    if (!c2 ) { ++failCount[2];  continue; }
    if (!c3 ) { ++failCount[3];  continue; }
    if (!c4 ) { ++failCount[4];  continue; }
    if (!c5 ) { ++failCount[5];  continue; }
    if (!c6 ) { ++failCount[6];  continue; }
    if (!c7 ) { ++failCount[7];  continue; }
    if (!c8 ) { ++failCount[8];  continue; }
    if (!c9 ) { ++failCount[9];  continue; }
    if (!c10) { ++failCount[10]; continue; }
    if (!c11)  { ++failCount[11]; continue; }
    if (!c12)  { ++failCount[12]; continue; }

    ++nPassed;
    tree->Fill();

    if (((i+1)%10000)==0 || (i+1)==nLoop)
      std::cout<<"\r  filled "<<(i+1)<<"/"<<nLoop<<" events"<<std::flush;
  }

  // Cutflow report
  std::cout<<"\n=== Cutflow ===\n";
  Int_t survived = nLoop;
  for (int i=0; i<Ncuts; ++i) {
    Int_t lost = failCount[i];
    float frac = 100.0f * lost / nLoop;
    survived -= lost;
    std::cout<<"After cut "<<(i+1)<<" ("<<cutNames[i]<<"): "
             <<survived<<" events (lost "<<lost<<" = "<<frac<<"%)\n";
  }
  std::cout<<"Final passed              : "<<nPassed
           <<" ("<<100.0f*nPassed/nLoop<<"% of initial)\n";

  // Write & cleanup
  std::cout<<"\nWriting "<<outputFile<<"\n";
  tree->Write();
  outFile->Close();
  delete reader;
  std::cout<<"Done.\n";
}

//ready 3
