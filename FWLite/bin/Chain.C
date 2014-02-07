#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

void Chain()
{
    gROOT->SetBatch(1);
    //gROOT->ProcessLine(".L ../interface/SimpleCandidate.h");
    gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc462/libJetMETTriggerAnalysisFWLite.so");

    TChain * chain  = new TChain("tree");
    TString outfilename = "compactified.L1ETM40.root";
    chain->Add("compactified_0.root");
    chain->Add("compactified_1.root");
    chain->Add("compactified_2.root");
    chain->Add("compactified_3.root");
    chain->Add("compactified_4.root");
    chain->Add("compactified_5.root");
    chain->Add("compactified_6.root");
    chain->Add("compactified_7.root");
    chain->Add("compactified_8.root");
    chain->Add("compactified_9.root");
    chain->Add("compactified_10.root");
    chain->Add("compactified_11.root");
    chain->Add("compactified_12.root");
    chain->Add("compactified_13.root");
    chain->Add("compactified_14.root");

    TFile* f1 = TFile::Open(outfilename, "RECREATE");
    TTree* t1 = (TTree*) chain->CopyTree("");
    std::clog << "From " << chain->GetEntriesFast() << " to " << t1->GetEntriesFast() << " entries." << std::endl;

    t1->Write();
    f1->Close();

    return;
}

// To run: root -l -b -q Chain.C+

