// compute_bjorken_pomeron.C
// Batch-mode usage (no interactive prompts):
//   root -l -b -q 'compute_bjorken_pomeron.C("in.root","out.root",10000,0)'

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
#include <TTree.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <Riostream.h>
#include <iostream>

void compute_bjorken_pomeron(const char* inputFile,
                             const char* outputFile,
                             Long64_t    maxEvents = -1,
                             Long64_t    offset    = 0) {
    // 1) Load Delphes library
    gSystem->Load("libDelphes");

    // 2) Chain up the Delphes tree
    TChain chain("Delphes");
    chain.Add(inputFile);
    ExRootTreeReader* reader = new ExRootTreeReader(&chain);
    TClonesArray* parts = reader->UseBranch("Particle");
    Long64_t nEntries = reader->GetEntries();

    // 3) Determine how many entries to process
    if (offset < 0) offset = 0;
    if (offset >= nEntries) {
        std::cerr << "Offset >= total entries, exiting.\n";
        return;
    }
    Long64_t avail = nEntries - offset;
    Long64_t nLoop = (maxEvents > 0 && maxEvents < avail ? maxEvents : avail);

    std::cout << "Processing " << nLoop 
              << " entries (from " << offset 
              << " to " << offset + nLoop - 1 << ")\n";

    // 4) Book histogram
    const Int_t nBins = 50;
    TH1D* hPomeronX = new TH1D(
        "hPomeronX",
        "Pomeron Bjorken x; x_{P}; Events",
        nBins, 0.0, 1.0
    );

    // 5) Prepare output file and TTree
    TFile outF(outputFile, "RECREATE");
    TTree* tree = new TTree("PomeronTree", "Pomeron x values");
    Double_t x_pomeron = 0.0;
    tree->Branch("x_pomeron", &x_pomeron, "x_pomeron/D");

    Long64_t fillCount = 0;

    // 6) Event loop: select all pomerons, ignore mother check
    for (Long64_t i = 0; i < nLoop; ++i) {
        reader->ReadEntry(offset + i);
        Int_t nPart = parts->GetEntriesFast();

        for (Int_t j = 0; j < nPart; ++j) {
            GenParticle* p = (GenParticle*) parts->At(j);
            if (p->PID != 990) continue;               // select pomeron only

            // compute Bjorken x_pomeron
            x_pomeron = p->E / 6500.0;
            hPomeronX->Fill(x_pomeron);
            tree->Fill();
            ++fillCount;
        }

        if (((i+1) % 10000) == 0 || (i+1) == nLoop) {
            std::cout << "\rProcessed " << (i+1)
                      << " / " << nLoop << std::flush;
        }
    }
    std::cout << "\nFilled tree entries: " << fillCount << "\n";

    // 7) Write outputs
    hPomeronX->Write();
    tree->Write();
    outF.Close();

    std::cout << "Saved histogram and tree to " << outputFile << "\n";

    // 8) Validation
    TFile inF(outputFile, "READ");
    TH1D* hCheck = (TH1D*) inF.Get("hPomeronX");
    if (!hCheck) {
        std::cerr << "ERROR: Histogram 'hPomeronX' not found in " << outputFile << "\n";
    } else {
        std::cout << "Validation: Histogram 'hPomeronX' has "
                  << hCheck->GetNbinsX() << " bins and "
                  << hCheck->GetEntries() << " entries.\n";
    }
    TTree* tCheck = (TTree*) inF.Get("PomeronTree");
    if (!tCheck) {
        std::cerr << "ERROR: Tree 'PomeronTree' not found in " << outputFile << "\n";
    } else {
        std::cout << "Validation: Tree 'PomeronTree' has "
                  << tCheck->GetEntries() << " entries.\n";
    }
    inF.Close();
}
