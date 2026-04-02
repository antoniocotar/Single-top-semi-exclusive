// extract_bjorken_pom_prot.C
// Batch-mode usage (no interactive prompts):
//   root -l -b -q 'extract_bjorken_pom_prot.C("in.root","out.root",10000)'

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#else
class ExRootTreeReader;
#endif

#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TMath.h>
#include <iostream>

void extract_bjorken_pom_prot(const char* inputFile,
                             const char* outputFile,
                             Long64_t    maxEvents) {
    // Load Delphes library
    gSystem->Load("libDelphes");

    // Open chain
    TChain chain("Delphes");
    chain.Add(inputFile);

    // Reader
    ExRootTreeReader* reader = new ExRootTreeReader(&chain);
    Long64_t nEntries = reader->GetEntries();
    Long64_t nLoop    = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;

    std::cout
        << "INPUT_FILE  = " << inputFile  << "\n"
        << "OUTPUT_FILE = " << outputFile << "\n"
        << "MAX_EVENTS  = " << maxEvents  << "\n"
        << "Looping over " << nLoop << " / " << nEntries << " events\n";

    // Grab Particle branch
    TClonesArray* parts = reader->UseBranch("Particle");

    // Booking histos (50 bins [0,1])
    const Int_t nBins = 50;
    TH1D* hPomQuark  = new TH1D("hPomQuark",  "Pomeron Light Quark x; x; Events", nBins, 0., 1.);
    TH1D* hPomGluon  = new TH1D("hPomGluon",  "Pomeron Gluon x; x; Events",        nBins, 0., 1.);
    TH1D* hProtQuark = new TH1D("hProtQuark", "Proton Light Quark x; x; Events",   nBins, 0., 1.);
    TH1D* hProtGluon = new TH1D("hProtGluon", "Proton Gluon x; x; Events",          nBins, 0., 1.);

    // Event loop
    for (Long64_t i = 0; i < nLoop; ++i) {
        reader->ReadEntry(i);
        Int_t nP = parts->GetEntriesFast();
        for (Int_t j = 0; j < nP; ++j) {
            GenParticle* p = (GenParticle*)parts->At(j);
            Int_t pid = p->PID, d1 = p->D1;
            if (d1 < 0 || d1 >= nP) continue;
            GenParticle* c = (GenParticle*)parts->At(d1);
            Int_t cpid = c->PID;
            Double_t x = c->E / 6500.0;

            if (pid == 990) {                // Pomeron parent
                if      (cpid == 21)           hPomGluon->Fill(x);
                else if (TMath::Abs(cpid)<=3) hPomQuark->Fill(x);
            }
            else if (pid == 2212) {         // Proton parent
                if      (cpid == 21)           hProtGluon->Fill(x);
                else if (TMath::Abs(cpid)<=3) hProtQuark->Fill(x);
            }
        }
        if (((i+1) % 10000) == 0 || (i+1)==nLoop)
            std::cout << "Processed " << (i+1) << " events\r" << std::flush;
    }

    std::cout << "\nWriting histograms to " << outputFile << "...\n";
    TFile* out = TFile::Open(outputFile, "RECREATE");
    hPomQuark->Write();
    hPomGluon->Write();
    hProtQuark->Write();
    hProtGluon->Write();
    out->Close();
    std::cout << "Done.\n";
}

