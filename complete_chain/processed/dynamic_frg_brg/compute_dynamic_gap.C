// compute_dynamic_gap.C
// Usage: root -l -b -q 'compute_dynamic_gap.C("in.root","out.root",sampleType,10000,0)'
//   sampleType: 0 = photon, 1 = pomeron, 2 = qcd

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
#include <vector>
#include <cmath>
#include <iomanip>

void compute_dynamic_gap(const char* inputFile,
                         const char* outputFile,
                         Int_t       sampleType,
                         Long64_t    maxEvents = -1,
                         Long64_t    offset    = 0) {
    // --- setup
    gSystem->Load("libDelphes");
    TChain chain("Delphes");
    chain.Add(inputFile);

    ExRootTreeReader* reader = new ExRootTreeReader(&chain);
    Long64_t nEntries = reader->GetEntries();
    if (offset<0) offset=0;
    if (offset>=nEntries) {
        std::cerr<<"Offset >= total entries, nothing to do\n";
        return;
    }
    Long64_t avail = nEntries - offset;
    Long64_t nLoop = (maxEvents>0 && maxEvents<avail) ? maxEvents : avail;

    std::cout<<std::fixed<<std::setprecision(0)
             <<"INPUT_FILE    = "<<inputFile<<"\n"
             <<"OUTPUT_FILE   = "<<outputFile<<"\n"
             <<"SAMPLE_TYPE   = "<<sampleType<<"\n"
             <<"MAX_EVENTS    = "<<maxEvents<<"\n"
             <<"OFFSET        = "<<offset<<"\n"
             <<"EVENTS_TO_RUN = "<<nLoop<<"\n";

    TClonesArray* parts = reader->UseBranch("Particle");

    // --- book histograms
    const Int_t    nBins = 50;
    const Double_t xMin  = 0.0, xMax = 5.0;
    TH1D* hFRG = new TH1D("hFRG",
                          "Forward gap;#Delta#eta_{min}(+3);Events",
                          nBins, xMin, xMax);
    TH1D* hBRG = new TH1D("hBRG",
                          "Backward gap;#Delta#eta_{min}(−3);Events",
                          nBins, xMin, xMax);

    Long64_t passed=0;
    for(Long64_t ie=0; ie<nLoop; ++ie) {
        reader->ReadEntry(offset+ie);
        Int_t nPart = parts->GetEntriesFast();

        bool select=false;
        // --- seed loop (do NOT require status==1 here!)
        for(Int_t j=0; j<nPart; ++j) {
            GenParticle* p = (GenParticle*)parts->At(j);
            if (sampleType==0) {
                // photon: PID==22, eta>0, mother proton
                if (p->PID==22 && p->Eta>0.0 && p->M1>=0) {
                    GenParticle* mom = (GenParticle*)parts->At(p->M1);
                    if (std::abs(mom->PID)==2212) { select=true; break; }
                }
            }
            else if (sampleType==1) {
                // pomeron: PID==990, eta>0
                if (p->PID==990 && p->Eta>0.0) { select=true; break; }
            }
            else {
                // qcd: accept all
                select=true;
                break;
            }
        }
        if (!select) continue;
        ++passed;

        // --- final-state charged collection (status==1 here!)
        std::vector<Double_t> validEtas;
        validEtas.reserve(nPart);
        for(Int_t j=0; j<nPart; ++j) {
            GenParticle* q = (GenParticle*)parts->At(j);
            if (q->Status!=1)           continue; // stable
            if (std::abs(q->PID)==2212) continue; // no protons
            if (q->Charge==0)           continue; // charged only
            if (q->PT<=0.2)             continue; // pT cut
            validEtas.push_back(q->Eta);
        }
        if (validEtas.empty()) continue;

        // --- compute and fill gaps
        Double_t minPlus=1e6, minMinus=1e6;
        for (auto eta: validEtas) {
            minPlus  = std::min(minPlus,  std::abs(eta+3.0));
            minMinus = std::min(minMinus, std::abs(eta-3.0));
        }
        hFRG->Fill(minPlus);
        hBRG->Fill(minMinus);
    }

    std::cout<<"Selected events = "<<passed<<" / "<<nLoop<<"\n";

    // --- write out
    TFile* out = TFile::Open(outputFile,"RECREATE");
    hFRG->Write();
    hBRG->Write();
    out->Close();
    std::cout<<"Histograms saved to "<<outputFile<<"\n";
}

