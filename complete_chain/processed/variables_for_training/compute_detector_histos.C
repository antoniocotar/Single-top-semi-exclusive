// compute_detector_histos.C
#ifdef __CLING__
  R__LOAD_LIBRARY(libDelphes)
  #include "classes/DelphesClasses.h"
  #include "external/ExRootAnalysis/ExRootTreeReader.h"
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

void compute_detector_histos(const char* inputFile,
                             const char* outputFile,
                             Long64_t    maxEvents = -1,
                             Long64_t    offset    = 0)
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
    
      // global HF extrema (any tower in |η|∈[3.0,5.2])
      if ((η >=  3.0f && η <=  5.2f) ||
          (η <= -3.0f && η >= -5.2f)) {
        maxHFEnergy = std::max(maxHFEnergy, E);
        minHFEnergy = std::min(minHFEnergy, E);
      }
    }
    
    // if truly no HF at all, mark global extrema invalid
    if (maxHFEnergy == std::numeric_limits<float>::lowest()) {
      maxHFEnergy = 0.0f;
      minHFEnergy = 0.0f;
    }
    
    // derived variables (now always ≥ 0)
    minOfMaxHF  = std::min(maxHFEnergyPlus,  maxHFEnergyMinus);
    minOfHFSums = std::min(hfSumPlus,        hfSumMinus);

      

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
