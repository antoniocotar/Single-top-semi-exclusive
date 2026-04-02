#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#else
class ExRootTreeReader;
#endif

#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TSystem.h"
#include "Riostream.h"

void compute_bjorken(const char* inputFile,
                     const char* outputFile,
                     Long64_t    maxEvents = -1) {

    // load Delphes
    gSystem->Load("libDelphes");

    // chain up the tree
    TChain chain("Delphes");
    chain.Add(inputFile);

    // reader
    ExRootTreeReader* reader = new ExRootTreeReader(&chain);
    Long64_t nEntries = reader->GetEntries();
    // decide how many to loop
    Long64_t nLoop = (maxEvents > 0 && maxEvents < nEntries) ? maxEvents : nEntries;

    std::cout << "DEBUG: Input = "   << inputFile
              << "\n       Output = " << outputFile
              << "\n       Entries = " << nEntries
              << "\n       Looping = " << nLoop
              << std::endl;

    // get the Particle branch
    TClonesArray* parts = reader->UseBranch("Particle");

    // histograms
    const Int_t nBins = 50;
    TH1D* hPhotonX = new TH1D("hPhotonX", "Photon Bjorken x;x;Events", nBins, 0., 1.);
    TH1D* hGluonX  = new TH1D("hGluonX",  "Gluon   Bjorken x;x;Events", nBins, 0., 1.);

    // loop
    for (Long64_t entry = 0; entry < nLoop; ++entry) {
        reader->ReadEntry(entry);
        Int_t nPart = parts->GetEntriesFast();
        for (Int_t i = 0; i < nPart; ++i) {
            GenParticle* p = (GenParticle*) parts->At(i);
            Int_t pid = p->PID;
            if (pid == 22 || pid == 21) {
                Int_t momIdx = p->M1;
                if (momIdx >= 0 && momIdx < nPart) {
                    GenParticle* mom = (GenParticle*) parts->At(momIdx);
                    if (mom->PID == 2212) {
                        Double_t x = p->E / 6500.0;
                        if (pid == 22)    hPhotonX->Fill(x);
                        else /* pid==21 */ hGluonX ->Fill(x);
                    }
                }
            }
        }
        if (entry % 10000 == 0) {
            std::cout << "DEBUG: processed " 
                      << entry << " / " << nLoop 
                      << " entries\r" << std::flush;
        }
    }
    std::cout << "\nDEBUG: Done loop, writing output...\n";

    // write out
    TFile* out = TFile::Open(outputFile, "RECREATE");
    hPhotonX->Write();
    hGluonX ->Write();
    out->Close();
    std::cout << "DEBUG: Histos written to " << outputFile << std::endl;
}

