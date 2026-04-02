// compute_particle_size.C
// Usage (batch mode):
// root -l -b -q 'compute_particle_size.C("in.root","out.root",10000,0)'

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
#include <Riostream.h>
#include <iostream>

void compute_particle_size(const char* inputFile,
                           const char* outputFile,
                           Long64_t    maxEvents = -1,
                           Long64_t    offset    = 0) {

    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader* reader = new ExRootTreeReader(&chain);
    Long64_t nEntries = reader->GetEntries();

    if (offset < 0) offset = 0;
    if (offset >= nEntries) {
        std::cerr << "Offset (" << offset << ") >= total entries (" << nEntries << "), nothing to do.\n";
        return;
    }
    Long64_t available = nEntries - offset;
    Long64_t nLoop = (maxEvents > 0 && maxEvents < available) ? maxEvents : available;

    std::cout
        << "INPUT_FILE  = " << inputFile  << "\n"
        << "OUTPUT_FILE = " << outputFile << "\n"
        << "MAX_EVENTS  = " << maxEvents  << "\n"
        << "OFFSET      = " << offset     << "\n"
        << "Looping over " << nLoop
        << " entries (from " << offset
        << " to " << offset + nLoop - 1 << ")\n";

    TClonesArray* particles = reader->UseBranch("Particle");
    TH1D* hNTrk = new TH1D("hNTrk", "Number of charged particles with cuts;N_{trk};Events", 80, 0, 400);

    for (Long64_t i = 0; i < nLoop; ++i) {
        Long64_t entry = offset + i;
        reader->ReadEntry(entry);

        Int_t nPart = particles->GetEntriesFast();
        int nChargedTrk = 0;

        for (Int_t j = 0; j < nPart; ++j) {
            GenParticle* p = (GenParticle*) particles->At(j);

            if (p->Status != 1) continue;            // only stable
            if (abs(p->PID) == 2212) continue;        // exclude protons
            if (p->Charge == 0) continue;             // only charged
            if (p->PT <= 0.4) continue;               // minimum pT
            if (p->Eta < -2.4 || p->Eta > 2.4) continue; // central region

            ++nChargedTrk;
        }

        hNTrk->Fill(nChargedTrk);
        if ((i+1) % 10000 == 0 || (i+1) == nLoop) {
            std::cout << "Processed " << (i+1) << " / " << nLoop << " entries\r" << std::flush;
        }
    }

    std::cout << "\nDone loop, writing histogram...\n";
    TFile* outFile = TFile::Open(outputFile, "UPDATE");
    hNTrk->Write();
    outFile->Close();
    std::cout << "Saved to: " << outputFile << "\n";
}

