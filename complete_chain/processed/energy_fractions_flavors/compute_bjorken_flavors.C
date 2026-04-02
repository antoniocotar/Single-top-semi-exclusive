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

void compute_bjorken_flavors(const char* inputFile,
                            const char* outputFile,
                            Long64_t    maxEvents = -1) {
    // load Delphes
    gSystem->Load("libDelphes");

    // chain up the tree
    TChain chain("Delphes");
    chain.Add(inputFile);

    // set up reader
    ExRootTreeReader* reader = new ExRootTreeReader(&chain);
    Long64_t nEntries = reader->GetEntries();
    Long64_t nLoop    = (maxEvents>0 && maxEvents<nEntries) ? maxEvents : nEntries;

    std::cout<<"DEBUG: input="<<inputFile
             <<"\n       output="<<outputFile
             <<"\n       entries="<<nEntries
             <<"\n       looping="<<nLoop<<std::endl;

    // get Particle branch
    TClonesArray* parts = reader->UseBranch("Particle");

    // book histograms
    const Int_t nBins=50;
    TH1D* hU_Pom = new TH1D("hU_Pom","u-quark from Pomeron;x;Events",nBins,0.,1.);
    TH1D* hD_Pom = new TH1D("hD_Pom","d-quark from Pomeron;x;Events",nBins,0.,1.);
    TH1D* hS_Pom = new TH1D("hS_Pom","s-quark from Pomeron;x;Events",nBins,0.,1.);
    TH1D* hU_Pro = new TH1D("hU_Pro","u-quark from Proton;x;Events",nBins,0.,1.);
    TH1D* hD_Pro = new TH1D("hD_Pro","d-quark from Proton;x;Events",nBins,0.,1.);
    TH1D* hS_Pro = new TH1D("hS_Pro","s-quark from Proton;x;Events",nBins,0.,1.);

    // loop events
    for(Long64_t entry=0; entry<nLoop; ++entry) {
        reader->ReadEntry(entry);
        Int_t nPart = parts->GetEntriesFast();
        for(Int_t i=0; i<nPart; ++i) {
            GenParticle* p = (GenParticle*)parts->At(i);
            Int_t pid = p->PID;
            Int_t d1  = p->D1;
            if(d1<0 || d1>=nPart) continue;
            GenParticle* child = (GenParticle*)parts->At(d1);
            Int_t cpid = child->PID;
            Double_t x = child->E/6500.0;

            if(pid==990) {
                // from Pomeron
                if(abs(cpid)==2)     hU_Pom->Fill(x);
                else if(abs(cpid)==1)hD_Pom->Fill(x);
                else if(abs(cpid)==3)hS_Pom->Fill(x);
            }
            else if(pid==2212) {
                // from Proton
                if(abs(cpid)==2)     hU_Pro->Fill(x);
                else if(abs(cpid)==1)hD_Pro->Fill(x);
                else if(abs(cpid)==3)hS_Pro->Fill(x);
            }
        }
        if(entry%10000==0) {
            std::cout<<"DEBUG: processed "<<entry<<" / "<<nLoop<<"\r"<<std::flush;
        }
    }
    std::cout<<"\nDEBUG: loop done, writing output...\n";

    // write out
    TFile* out = TFile::Open(outputFile,"RECREATE");
    hU_Pom->Write();
    hD_Pom->Write();
    hS_Pom->Write();
    hU_Pro->Write();
    hD_Pro->Write();
    hS_Pro->Write();
    out->Close();
    std::cout<<"DEBUG: histograms saved to "<<outputFile<<std::endl;
}

