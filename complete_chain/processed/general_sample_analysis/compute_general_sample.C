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
#include <TClonesArray.h>
#include <TFile.h>
#include <TH1D.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

const char* sample_label(Int_t sampleType) {
  if (sampleType == 0) return "photon";
  if (sampleType == 1) return "pomeron";
  return "background";
}

bool event_matches_sample(TClonesArray* parts, Int_t sampleType) {
  const Int_t nPart = parts->GetEntriesFast();
  for (Int_t i = 0; i < nPart; ++i) {
    auto* p = (GenParticle*)parts->At(i);
    if (sampleType == 0) {
      if (p->PID == 22 && p->Eta > 0.0 && p->M1 >= 0 && p->M1 < nPart) {
        auto* mom = (GenParticle*)parts->At(p->M1);
        if (std::abs(mom->PID) == 2212) return true;
      }
    } else if (sampleType == 1) {
      if (p->PID == 990 && p->Eta > 0.0) return true;
    } else {
      return true;
    }
  }
  return false;
}

}  // namespace

void compute_general_sample(const char* inputFile,
                            const char* outputFile,
                            Int_t       sampleType,
                            Long64_t    maxEvents = -1,
                            Long64_t    offset    = 0) {
  gSystem->Load("libDelphes");

  TChain chain("Delphes");
  chain.Add(inputFile);

  ExRootTreeReader* reader = new ExRootTreeReader(&chain);
  const Long64_t nEntries = reader->GetEntries();

  if (offset < 0) offset = 0;
  if (offset >= nEntries) {
    std::cerr << "Offset (" << offset << ") >= entries (" << nEntries << "). Aborting.\n";
    return;
  }

  const Long64_t available = nEntries - offset;
  const Long64_t nLoop = (maxEvents > 0 && maxEvents < available) ? maxEvents : available;

  std::cout << "INPUT_FILE   = " << inputFile << "\n"
            << "OUTPUT_FILE  = " << outputFile << "\n"
            << "SAMPLE_TYPE  = " << sampleType << " (" << sample_label(sampleType) << ")\n"
            << "MAX_EVENTS   = " << maxEvents << "\n"
            << "OFFSET       = " << offset << "\n"
            << "EVENTS_TO_RUN= " << nLoop << "\n";

  TClonesArray* branchParticle  = reader->UseBranch("Particle");
  TClonesArray* branchJet       = reader->UseBranch("Jet");
  TClonesArray* branchTrack     = reader->UseBranch("Track");
  TClonesArray* branchTower     = reader->UseBranch("Tower");
  TClonesArray* branchElectron  = reader->UseBranch("Electron");
  TClonesArray* branchMuon      = reader->UseBranch("Muon");
  TClonesArray* branchMissingET = reader->UseBranch("MissingET");

  TH1D* hNTrk      = new TH1D("hNTrk", "Charged multiplicity;N_{trk};Events", 80, 0, 400);
  TH1D* hFRG       = new TH1D("hFRG", "Forward gap;#Delta#eta_{min}(+3);Events", 50, 0.0, 5.0);
  TH1D* hBRG       = new TH1D("hBRG", "Backward gap;#Delta#eta_{min}(-3);Events", 50, 0.0, 5.0);
  TH1D* hPhotonX   = new TH1D("hPhotonX", "Photon Bjorken x;x;Events", 50, 0.0, 1.0);
  TH1D* hGluonX    = new TH1D("hGluonX", "Gluon Bjorken x;x;Events", 50, 0.0, 1.0);
  TH1D* hPomeronX  = new TH1D("hPomeronX", "Pomeron Bjorken x;x_{P};Events", 50, 0.0, 1.0);
  TH1D* hPomQuark  = new TH1D("hPomQuark", "Pomeron light-quark x;x;Events", 50, 0.0, 1.0);
  TH1D* hPomGluon  = new TH1D("hPomGluon", "Pomeron gluon x;x;Events", 50, 0.0, 1.0);
  TH1D* hProtQuark = new TH1D("hProtQuark", "Proton light-quark x;x;Events", 50, 0.0, 1.0);
  TH1D* hProtGluon = new TH1D("hProtGluon", "Proton gluon x;x;Events", 50, 0.0, 1.0);
  TH1D* hU_Pom     = new TH1D("hU_Pom", "u from Pomeron;x;Events", 50, 0.0, 1.0);
  TH1D* hD_Pom     = new TH1D("hD_Pom", "d from Pomeron;x;Events", 50, 0.0, 1.0);
  TH1D* hS_Pom     = new TH1D("hS_Pom", "s from Pomeron;x;Events", 50, 0.0, 1.0);
  TH1D* hU_Pro     = new TH1D("hU_Pro", "u from Proton;x;Events", 50, 0.0, 1.0);
  TH1D* hD_Pro     = new TH1D("hD_Pro", "d from Proton;x;Events", 50, 0.0, 1.0);
  TH1D* hS_Pro     = new TH1D("hS_Pro", "s from Proton;x;Events", 50, 0.0, 1.0);

  TFile* outFile = TFile::Open(outputFile, "RECREATE");
  TTree* tree = new TTree("Events", "General observables per selected event");

  Int_t sampleTypeOut = sampleType;
  Float_t rapidityGapPlus = -1.0f;
  Float_t rapidityGapMinus = -1.0f;
  Float_t hfSumPlus = 0.0f;
  Float_t hfSumMinus = 0.0f;
  Float_t nchPV = 0.0f;
  Float_t nForwardJets = 0.0f;
  Float_t nTotalJets = 0.0f;
  Float_t nJetsSel = 0.0f;
  Float_t nBtagSel = 0.0f;
  Float_t nLightSel = 0.0f;
  Float_t leadLightJetPt = 0.0f;
  Float_t METval = 0.0f;
  Float_t maxHFEnergyPlus = 0.0f;
  Float_t maxHFEnergyMinus = 0.0f;
  Float_t minOfMaxHF = 0.0f;
  Float_t minOfHFSums = 0.0f;
  Float_t nLeptons = 0.0f;

  tree->Branch("sampleType", &sampleTypeOut, "sampleType/I");
  tree->Branch("rapidityGapPlus", &rapidityGapPlus, "rapidityGapPlus/F");
  tree->Branch("rapidityGapMinus", &rapidityGapMinus, "rapidityGapMinus/F");
  tree->Branch("hfSumPlus", &hfSumPlus, "hfSumPlus/F");
  tree->Branch("hfSumMinus", &hfSumMinus, "hfSumMinus/F");
  tree->Branch("nchPV", &nchPV, "nchPV/F");
  tree->Branch("nForwardJets", &nForwardJets, "nForwardJets/F");
  tree->Branch("nTotalJets", &nTotalJets, "nTotalJets/F");
  tree->Branch("nJetsSel", &nJetsSel, "nJetsSel/F");
  tree->Branch("nBtagSel", &nBtagSel, "nBtagSel/F");
  tree->Branch("nLightSel", &nLightSel, "nLightSel/F");
  tree->Branch("leadLightJetPt", &leadLightJetPt, "leadLightJetPt/F");
  tree->Branch("MET", &METval, "MET/F");
  tree->Branch("maxHFEnergyPlus", &maxHFEnergyPlus, "maxHFEnergyPlus/F");
  tree->Branch("maxHFEnergyMinus", &maxHFEnergyMinus, "maxHFEnergyMinus/F");
  tree->Branch("minOfMaxHF", &minOfMaxHF, "minOfMaxHF/F");
  tree->Branch("minOfHFSums", &minOfHFSums, "minOfHFSums/F");
  tree->Branch("nLeptons", &nLeptons, "nLeptons/F");

  Long64_t nSelected = 0;

  for (Long64_t i = 0; i < nLoop; ++i) {
    reader->ReadEntry(offset + i);

    if (!event_matches_sample(branchParticle, sampleType)) continue;
    ++nSelected;

    const Int_t nPart = branchParticle->GetEntriesFast();
    std::vector<double> gapEtas;
    gapEtas.reserve(nPart);

    for (Int_t j = 0; j < nPart; ++j) {
      auto* p = (GenParticle*)branchParticle->At(j);
      const Int_t pid = p->PID;

      if ((pid == 22 || pid == 21) && p->M1 >= 0 && p->M1 < nPart) {
        auto* mom = (GenParticle*)branchParticle->At(p->M1);
        if (mom->PID == 2212) {
          const Double_t x = p->E / 6500.0;
          if (pid == 22) hPhotonX->Fill(x);
          if (pid == 21) hGluonX->Fill(x);
        }
      }

      if (pid == 990) {
        hPomeronX->Fill(p->E / 6500.0);
      }

      const Int_t d1 = p->D1;
      if (d1 >= 0 && d1 < nPart) {
        auto* child = (GenParticle*)branchParticle->At(d1);
        const Int_t cpid = child->PID;
        const Double_t xChild = child->E / 6500.0;

        if (pid == 990) {
          if (cpid == 21) hPomGluon->Fill(xChild);
          else if (std::abs(cpid) <= 3) hPomQuark->Fill(xChild);

          if (std::abs(cpid) == 2) hU_Pom->Fill(xChild);
          else if (std::abs(cpid) == 1) hD_Pom->Fill(xChild);
          else if (std::abs(cpid) == 3) hS_Pom->Fill(xChild);
        } else if (pid == 2212) {
          if (cpid == 21) hProtGluon->Fill(xChild);
          else if (std::abs(cpid) <= 3) hProtQuark->Fill(xChild);

          if (std::abs(cpid) == 2) hU_Pro->Fill(xChild);
          else if (std::abs(cpid) == 1) hD_Pro->Fill(xChild);
          else if (std::abs(cpid) == 3) hS_Pro->Fill(xChild);
        }
      }

      if (p->Status != 1) continue;
      if (std::abs(pid) == 2212) continue;
      if (p->Charge == 0) continue;
      if (p->PT <= 0.2) continue;
      gapEtas.push_back(p->Eta);
    }

    Int_t nChargedTrk = 0;
    nchPV = 0.0f;
    std::vector<double> recoTrackEtas;
    const Int_t nTrk = branchTrack->GetEntriesFast();
    recoTrackEtas.reserve(nTrk);
    for (Int_t t = 0; t < nTrk; ++t) {
      auto* tr = (Track*)branchTrack->At(t);
      if (tr->PT > 0.4 && std::fabs(tr->Eta) < 2.4) {
        ++nChargedTrk;
        nchPV += 1.0f;
        recoTrackEtas.push_back(tr->Eta);
      }
    }
    hNTrk->Fill(nChargedTrk);

    rapidityGapPlus = -1.0f;
    rapidityGapMinus = -1.0f;
    if (!gapEtas.empty()) {
      double minPlus = 1e6;
      double minMinus = 1e6;
      for (double eta : gapEtas) {
        minPlus = std::min(minPlus, std::abs(eta + 3.0));
        minMinus = std::min(minMinus, std::abs(eta - 3.0));
      }
      rapidityGapPlus = static_cast<Float_t>(minPlus);
      rapidityGapMinus = static_cast<Float_t>(minMinus);
      hFRG->Fill(rapidityGapPlus);
      hBRG->Fill(rapidityGapMinus);
    }

    nTotalJets = static_cast<Float_t>(branchJet->GetEntriesFast());
    nForwardJets = 0.0f;
    nJetsSel = 0.0f;
    nBtagSel = 0.0f;
    nLightSel = 0.0f;
    leadLightJetPt = 0.0f;

    for (Int_t j = 0; j < branchJet->GetEntriesFast(); ++j) {
      auto* jet = (Jet*)branchJet->At(j);
      if (std::fabs(jet->Eta) > 3.0) nForwardJets += 1.0f;
      if (jet->BTag == 0 && jet->PT > leadLightJetPt) leadLightJetPt = jet->PT;
      if (jet->PT > 25.0 && std::fabs(jet->Eta) < 2.4) {
        nJetsSel += 1.0f;
        if (jet->BTag > 0) nBtagSel += 1.0f;
        else nLightSel += 1.0f;
      }
    }

    hfSumPlus = 0.0f;
    hfSumMinus = 0.0f;
    maxHFEnergyPlus = 0.0f;
    maxHFEnergyMinus = 0.0f;
    for (Int_t t = 0; t < branchTower->GetEntriesFast(); ++t) {
      auto* tow = (Tower*)branchTower->At(t);
      const float eta = tow->Eta;
      const float e = tow->E;
      if (eta >= 3.0f && eta <= 5.2f && e >= 0.0f) {
        hfSumPlus += e;
        maxHFEnergyPlus = std::max(maxHFEnergyPlus, e);
      }
      if (eta <= -3.0f && eta >= -5.2f && e >= 0.0f) {
        hfSumMinus += e;
        maxHFEnergyMinus = std::max(maxHFEnergyMinus, e);
      }
    }
    minOfMaxHF = std::min(maxHFEnergyPlus, maxHFEnergyMinus);
    minOfHFSums = std::min(hfSumPlus, hfSumMinus);

    nLeptons = 0.0f;
    for (Int_t e = 0; e < branchElectron->GetEntriesFast(); ++e) {
      auto* el = (Electron*)branchElectron->At(e);
      if (el->PT > 15.0 && std::fabs(el->Eta) < 2.4) nLeptons += 1.0f;
    }
    for (Int_t m = 0; m < branchMuon->GetEntriesFast(); ++m) {
      auto* mu = (Muon*)branchMuon->At(m);
      if (mu->PT > 10.0 && std::fabs(mu->Eta) < 2.4) nLeptons += 1.0f;
    }

    METval = 0.0f;
    if (branchMissingET->GetEntriesFast() > 0) {
      auto* metobj = (MissingET*)branchMissingET->At(0);
      METval = metobj->MET;
    }

    tree->Fill();

    if (((i + 1) % 10000) == 0 || (i + 1) == nLoop) {
      std::cout << "Processed " << (i + 1) << " / " << nLoop
                << " events, selected " << nSelected << "\r" << std::flush;
    }
  }

  std::cout << "\nSelected events = " << nSelected << " / " << nLoop << "\n";

  tree->Write();
  hNTrk->Write();
  hFRG->Write();
  hBRG->Write();
  hPhotonX->Write();
  hGluonX->Write();
  hPomeronX->Write();
  hPomQuark->Write();
  hPomGluon->Write();
  hProtQuark->Write();
  hProtGluon->Write();
  hU_Pom->Write();
  hD_Pom->Write();
  hS_Pom->Write();
  hU_Pro->Write();
  hD_Pro->Write();
  hS_Pro->Write();
  outFile->Close();

  delete reader;
  std::cout << "Saved output to " << outputFile << "\n";
}
