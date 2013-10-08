#include <string>
#include "TFile.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TTree.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// CMSSW: FWLite
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
//// CMSSW: DataFormats
//#include "DataFormats/Common/interface/TriggerResults.h"

#include "JetMETTriggerAnalysis/FWLite/interface/Handler.h"
#include "JetMETTriggerAnalysis/FWLite/interface/SimpleCandidate.h"
#endif

//______________________________________________________________________________
// PFJetID, refer PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
bool PFJetID(const pat::Jet& jet, std::string wp="L") {
    double eta = jet.eta();
    float chf = jet.chargedHadronEnergyFraction();
	float nhf = ( jet.neutralHadronEnergyFraction() + jet.HFHadronEnergyFraction() );
	float cef = jet.chargedEmEnergyFraction();
	float nef = jet.neutralEmEnergyFraction();
	int nch = jet.chargedMultiplicity();
	int nconstituents = jet.numberOfDaughters();
    bool ret = true;
    if (wp == "L") {
        ret = (nconstituents > 1) &&
              (nef < 0.99) &&
              (nhf < 0.99) &&
              ((cef < 0.99) || std::abs(eta) > 2.4) &&
              ((chf > 0.) || std::abs(eta) > 2.4) &&
              ((nch > 0) || std::abs(eta) > 2.4);
    } else if (wp == "M") {
        ret = (nconstituents > 1) &&
              (nef < 0.95) &&
              (nhf < 0.95) &&
              ((cef < 0.99) || std::abs(eta) > 2.4) &&
              ((chf > 0.) || std::abs(eta) > 2.4) &&
              ((nch > 0) || std::abs(eta) > 2.4);
    } else if (wp == "T") {
        ret = (nconstituents > 1) &&
              (nef < 0.90) &&
              (nhf < 0.90) &&
              ((cef < 0.99) || std::abs(eta) > 2.4) &&
              ((chf > 0.) || std::abs(eta) > 2.4) &&
              ((nch > 0) || std::abs(eta) > 2.4);
    }
    return ret;
}


//______________________________________________________________________________
template<typename T>
void simple_fill(const T& cand, SimpleLorentzVector& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
}

template<typename T>
void simple_fill(const T& cand, SimpleParticle& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
    simple.charge = cand.charge();
    simple.pdgId  = cand.pdgId();
    simple.isPU   = false;  // need to update manually!
}


template<typename T>
void simple_fill(const T& cand, SimpleJet& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
    simple.jetIdL = PFJetId(cand, "L");
    simple.jetIdT = PFJetId(cand, "T");
    return;
}

template<typename T>
void simple_fill(const T& cand, SimpleMET& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
    simple.sumEt  = cand.sumEt();
    return;
}


//______________________________________________________________________________
int main(int argc, char *argv[]) {
    // Load FWLite library
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
    //gSystem->Load("libDataFormatsFWLite");
    
    // Generate libraries
    //gInterpreter->GenerateDictionary("SimpleMET", "JetMETTriggerAnalysis/FWLite/interface/SimpleCandidate.h");

    // Load config file
    std::string cfgfilename = "compactify_cfg.py";
    if (argc >= 2) {
        cfgfilename = argv[1];
    }
    PythonProcessDesc builder(cfgfilename.c_str());
    const edm::ParameterSet& inputpset    = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("input");
    const edm::ParameterSet& outputpset   = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("output");
    const edm::ParameterSet& analyzerpset = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("analyzer");
    const edm::ParameterSet& handlerpset  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("handler");

    // Read from input pset
    const std::vector<std::string>& infilenames = inputpset.getParameter<std::vector<std::string> >("fileNames");
    const int maxEvents = inputpset.getParameter<int>("maxEvents");
    const int runMin = inputpset.getParameter<int>("runMin");
    const int runMax = inputpset.getParameter<int>("runMax");
    const int skipEvents = inputpset.getParameter<int>("skipEvents");
    //const int reportEvery = inputpset.getParameter<int>("reportEvery");
    std::vector<edm::LuminosityBlockRange> lumisToProcess;
    if (inputpset.exists("lumisToProcess") ) {
        const std::vector<edm::LuminosityBlockRange>& lumisTemp =
            inputpset.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess");
        lumisToProcess.resize( lumisTemp.size() );
        copy(lumisTemp.begin(), lumisTemp.end(), lumisToProcess.begin() );
    }
    
    // Read from output pset
    const std::string outfilename = outputpset.getParameter<std::string>("fileName");
    
    // Read from analyzer pset
    const std::vector<std::string>& triggers = analyzerpset.getParameter<std::vector<std::string> >("triggers");
    const std::vector<std::string>& metfilters = analyzerpset.getParameter<std::vector<std::string> >("metfilters");
    //const double pfjetPtMin = analyzerpset.getParameter<double>("pfjetPtMin");
    //const double pfjetEtaMax = analyzerpset.getParameter<double>("pfjetEtaMax");
    //const double pfjetEtaMaxCtr = analyzerpset.getParameter<double>("pfjetEtaMaxCtr");
    //const double calojetPtMin = analyzerpset.getParameter<double>("calojetPtMin");
    //const double calojetEtaMax = analyzerpset.getParameter<double>("calojetEtaMax");
    //const double calojetEtaMaxCtr = analyzerpset.getParameter<double>("calojetEtaMaxCtr");
    //const double calometcleanPtMin = analyzerpset.getParameter<double>("calometcleanPtMin");
    //const double calometjetidPtMin = analyzerpset.getParameter<double>("calometjetidPtMin");
    const bool isData = analyzerpset.getParameter<bool>("isData");
    const bool verbose = analyzerpset.getParameter<bool>("verbose");
    
    // Read from handler pset
    Handler handler(handlerpset, isData);
    
    
    //__________________________________________________________________________
    // Prepare output tree
    std::vector<SimpleLorentzVector> hltCaloJets;
    std::vector<SimpleLorentzVector> hltPFJets;
    SimpleMET hltCaloMET;
    SimpleMET hltPFMET;
    // event weights and event flags
    
    TFile* outfile = new TFile(outfilename.c_str(), "RECREATE");
    TTree* outtree = new TTree("Events", "Events");
    outtree->Branch("hltCaloJets", &hltCaloJets);
    outtree->Branch("hltCaloJets", &hltPFJets);
    outtree->Branch("hltCaloMET", &hltCaloMET);
    outtree->Branch("hltCaloMET", &hltPFMET);
    
    
    //__________________________________________________________________________
    // Loop over events
    fwlite::ChainEvent ev(infilenames);
    int ievent = 0;
    int jevent = 0;
    for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievent) {
        // Skip events
        if (ievent < skipEvents)  continue;
        if (runMin > 0 && (int)ev.id().run() < runMin)  continue; 
        if (runMax > 0 && (int)ev.id().run() > runMax)  continue;
        // Stop event loop
        if (maxEvents > 0 && (ievent >= (maxEvents + skipEvents)) )  break;

        const edm::EventBase& eventbase = ev;
        handler.get(eventbase);
        
        
        //______________________________________________________________________
        // hltCaloJets
        hltCaloJets.clear();
        for (unsigned int i = 0; i < handler.hltCaloJets->size(); ++i) {
            SimpleLorentzVector simple;
            simple_fill(handler.hltCaloJets->at(i), simple);
            hltCaloJets.push_back(simple);
        }
        
        // hltPFJets
        hltPFJets.clear();
        for (unsigned int i = 0; i < handler.hltPFJets->size(); ++i) {
            SimpleLorentzVector simple;
            simple_fill(handler.hltPFJets->at(i), simple);
            hltPFJets.push_back(simple);
        }
        
        // hltCaloMET
        simple_fill(handler.hltCaloMETs->at(0), hltCaloMET);
        
        // hltPFMET
        simple_fill(handler.hltPFMETs->at(0), hltPFMET);
        
        
        //______________________________________________________________________
        outtree->Fill();

        ++jevent;
    }
    
    outfile->cd();
    outtree->Write();
    outfile->Close();
    
    return 0;
}

