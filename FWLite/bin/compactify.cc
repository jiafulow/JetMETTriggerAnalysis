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
// CMSSW: Utilities
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
// CMSSW: JetMETTriggerAnalysis
#include "JetMETTriggerAnalysis/FWLite/interface/Handler.h"
#include "JetMETTriggerAnalysis/FWLite/interface/SimpleCandidate.h"
#include "JetMETTriggerAnalysis/FWLite/interface/JSONFilter.h"
#endif


//______________________________________________________________________________
template<typename T>
void simple_fill(const T& cand, simple::LorentzVector& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
    return;
}

template<typename T>
void simple_fill(const T& cand, bool isPU, simple::Particle& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
    simple.charge = cand.charge();
    simple.pdgId  = cand.pdgId();
    simple.isPU   = isPU;
    return;
}

template<typename T>
void simple_fill(const T& cand, simple::Particle& simple) {
    simple_fill(cand, false, simple);
    return;
}

template<typename T>
void simple_fill(const T& cand, bool jetID, simple::Jet& simple) {
    const pat::Jet* patJet = dynamic_cast<const pat::Jet*>(&cand);
    
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pz     = cand.p4().pz();
    simple.E      = cand.p4().E();
    simple.pt     = cand.p4().pt();
    simple.eta    = cand.p4().eta();
    simple.phi    = cand.p4().phi();
    simple.rawpt  = (patJet != 0 ? patJet->correctedP4(0).pt() : -999);
    simple.jetID  = jetID;
    return;
}

template<typename T>
void simple_fill(const T& cand, simple::Jet& simple) {
    simple_fill(cand, true, simple);
    return;
}

template<typename T>
void simple_fill(const T& cand, simple::MET& simple) {
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

void simple_fill(const edm::EventBase& iEvent, unsigned int nPV, bool json, simple::Event& simple) {
    simple.run    = iEvent.id().run();
    simple.lumi   = iEvent.id().luminosityBlock();
    simple.event  = iEvent.id().event();
    simple.nPV    = nPV;
    simple.json   = json;
    return;
}


//______________________________________________________________________________
int main(int argc, char *argv[]) {
    // Load FWLite library
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
    //gSystem->Load("libDataFormatsFWLite");
    
    // Generate libraries
    //gInterpreter->GenerateDictionary("simple::MET", "JetMETTriggerAnalysis/FWLite/interface/simple::Candidate.h");

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
    // #FIXME - PU weight, PDF weight
    simple::Event simpleEvent;
    bool triggerFlags[99]; int nTriggerFlags = triggers.size() + 1;  // last one is "OR" combination
    bool metfilterFlags[99]; int nMetfilterFlags = metfilters.size() + 1;  // last one is "AND" combination
    // HLT
    std::vector<simple::Particle> hltPFCandidates;
    std::vector<simple::CaloJet> hltCaloJets;
    std::vector<simple::PFJet> hltPFJets;
    simple::MET hltCaloMET;
    simple::MET hltCaloMETClean;
    simple::MET hltCaloMETJetIDClean;
    simple::MET hltPFMET;
    simple::MET hltTrackMET;
    double hltRho_kt6CaloJets, hltRho_kt6PFJets;
    // non-HLT
    std::vector<simple::Particle> recoPFCandidates;
    std::vector<simple::PFJet> patJets;
    simple::MET patMET;
    double recoRho_kt6CaloJets, recoRho_kt6PFJets;
    
    
    TFile* outfile = new TFile(outfilename.c_str(), "RECREATE");
    TTree* outtree = new TTree("Events", "Events");
    outtree->Branch("event", &simpleEvent);
    outtree->Branch("triggerFlags", triggerFlags, Form("triggerFlags[%i]/b", nTriggerFlags));
    outtree->Branch("metfilterFlags", metfilterFlags, Form("metfilterFlags[%i]/b", nMetfilterFlags));
    // HLT
    outtree->Branch("hltPFCandidates", &hltPFCandidates);
    outtree->Branch("hltCaloJets", &hltCaloJets);
    outtree->Branch("hltPFJets", &hltPFJets);
    outtree->Branch("hltCaloMET", &hltCaloMET);
    outtree->Branch("hltCaloMETClean", &hltCaloMETClean);
    outtree->Branch("hltCaloMETJetIDClean", &hltCaloMETJetIDClean);
    outtree->Branch("hltPFMET", &hltPFMET);
    outtree->Branch("hltTrackMET", &hltTrackMET);
    outtree->Branch("hltRho_kt6CaloJets", &hltRho_kt6CaloJets);
    outtree->Branch("hltRho_kt6PFJets", &hltRho_kt6PFJets);
    // non-HLT
    outtree->Branch("recoPFCandidates", &recoPFCandidates);
    outtree->Branch("patJets", &patJets);
    outtree->Branch("patMET", &patMET);
    outtree->Branch("recoRho_kt6CaloJets", &recoRho_kt6CaloJets);
    outtree->Branch("recoRho_kt6PFJets", &recoRho_kt6PFJets);
    
    //__________________________________________________________________________
    // JetID selector
    JetIDSelectionFunctor jetIDSelectionFunctor( analyzerpset.getParameter<edm::ParameterSet>("jetIDSelector") );
    pat::strbitset jetIDBitSet = jetIDSelectionFunctor.getBitTemplate();
    PFJetIDSelectionFunctor pfJetIDSelectionFunctor( analyzerpset.getParameter<edm::ParameterSet>("pfJetIDSelector") );
    pat::strbitset pfJetIDBitSet = pfJetIDSelectionFunctor.getBitTemplate();
    
    //__________________________________________________________________________
    // Loop over events
    fwlite::ChainEvent ev(infilenames);
    int ievent = 0;
    int jevent = 0;
    std::cout << "<<< number of input files: " << infilenames.size() << std::endl;
    for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievent) {
        // Skip events
        if (ievent < skipEvents)  continue;
        if (runMin > 0 && (int)ev.id().run() < runMin)  continue; 
        if (runMax > 0 && (int)ev.id().run() > runMax)  continue;
        // Stop event loop
        if (maxEvents > 0 && (ievent >= (maxEvents + skipEvents)) )  break;
        
        if (verbose)  std::cout << "=== iEvent " << ievent << ", jEvent " << jevent << " ===" << std::endl;
        const edm::EventBase& eventbase = ev;
        handler.get(eventbase);
        
        
        //______________________________________________________________________
        // Event info
        if (verbose)  std::cout << "compactify: Begin filling event info..." << std::endl;
        bool goodjson = jsonContainsEvent(lumisToProcess, ev);
        unsigned int nPV = handler.recoVertices->size();
        simple_fill(ev, nPV, goodjson, simpleEvent);
        
        //______________________________________________________________________
        // Trigger results
        if (verbose)  std::cout << "compactify: Begin filling trigger results..." << std::endl;
        edm::TriggerResultsByName resultsByName = ev.triggerResultsByName("HLT");  //< from original HLT
        //edm::TriggerResultsByName resultsByName = ev.triggerResultsByName("HLT3");  //< from production mode (always true)
        //edm::TriggerResultsByName resultsByName = ev.triggerResultsByName("HLT3PB");  //< from filtering mode
        bool triggerFlagsOR = false;
        for (unsigned int i = 0; i < triggers.size(); ++i) {
            bool accept = resultsByName.accept(triggers.at(i));
            triggerFlags[i] = accept;
            if (accept)  triggerFlagsOR = true;
        }
        triggerFlags[triggers.size()] = triggerFlagsOR;
        
        //______________________________________________________________________
        // MET filter results
        if (verbose)  std::cout << "compactify: Begin filling MET filter results..." << std::endl;
        edm::TriggerResultsByName resultsByNamePAT = ev.triggerResultsByName("PAT");
        bool metfilterFlagsAND = true;
        for (unsigned int i = 0; i < metfilters.size(); ++i) {
            bool accept = resultsByNamePAT.accept(metfilters.at(i));
            metfilterFlags[i] = accept;
            if (!accept)  metfilterFlagsAND = false;
        }
        metfilterFlags[metfilters.size()] = metfilterFlagsAND;
        
        //______________________________________________________________________
        // hltPFCandidates
        if (verbose)  std::cout << "compactify: Begin filling hltPFCandidates..." << std::endl;
        hltPFCandidates.clear();
        assert(handler.hltPFPileUpFlags->size() == handler.hltPFCandidates->size());
        for (unsigned int i = 0; i < handler.hltPFCandidates->size(); ++i) {
            simple::Particle simple;
            const reco::PFCandidate& cand = handler.hltPFCandidates->at(i);
            const bool& isPU = handler.hltPFPileUpFlags->at(i);
            simple_fill(cand, isPU, simple);
            hltPFCandidates.push_back(simple);
        }
        
        //______________________________________________________________________
        // hltCaloJets
        if (verbose)  std::cout << "compactify: Begin filling hltCaloJets..." << std::endl;
        hltCaloJets.clear();
        for (unsigned int i = 0; i < handler.hltCaloJets->size(); ++i) {
            const reco::CaloJet& jet = handler.hltCaloJets->at(i);
            simple::CaloJet simple;
            simple_fill(jet, simple);
            hltCaloJets.push_back(simple);
        }
        
        //______________________________________________________________________
        // hltPFJets
        if (verbose)  std::cout << "compactify: Begin filling hltPFJets..." << std::endl;
        hltPFJets.clear();
        for (unsigned int i = 0; i < handler.hltPFJets->size(); ++i) {
            const reco::PFJet& jet = handler.hltPFJets->at(i);
            simple::PFJet simple;
            //simple_fill(jet, pfJetIDSelectionFunctor(jet, pfJetIDBitSet), simple);
            simple_fill(jet, simple);
            hltPFJets.push_back(simple);
        }
        
        //______________________________________________________________________
        // hltCaloMET
        if (verbose)  std::cout << "compactify: Begin filling hltCaloMET..." << std::endl;
        simple_fill(handler.hltCaloMETs->at(0), hltCaloMET);
        
        //______________________________________________________________________
        // hltCaloMETClean
        if (verbose)  std::cout << "compactify: Begin filling hltCaloMETClean..." << std::endl;
        simple_fill(handler.hltCaloMETCleans->at(0), hltCaloMETClean);
        
        //______________________________________________________________________
        // hltCaloMETJetIDClean
        if (verbose)  std::cout << "compactify: Begin filling hltCaloMETJetIDClean..." << std::endl;
        simple_fill(handler.hltCaloMETJetIDCleans->at(0), hltCaloMETJetIDClean);
        
        //______________________________________________________________________
        // hltPFMET
        if (verbose)  std::cout << "compactify: Begin filling hltPFMET..." << std::endl;
        simple_fill(handler.hltPFMETs->at(0), hltPFMET);
        
        //______________________________________________________________________
        // hltTrackMET
        if (verbose)  std::cout << "compactify: Begin filling hltTrackMET..." << std::endl;
        simple_fill(handler.hltTrackMETs->at(0), hltTrackMET);
        
        //______________________________________________________________________
        // hltRhos
        if (verbose)  std::cout << "compactify: Begin filling hltRhos..." << std::endl;
        hltRho_kt6CaloJets = *(handler.hltRho_kt6CaloJets);
        hltRho_kt6PFJets = *(handler.hltRho_kt6PFJets);
        
        //______________________________________________________________________
        // recoPFCandidates
        if (verbose)  std::cout << "compactify: Begin filling recoPFCandidates..." << std::endl;
        recoPFCandidates.clear();
        assert(handler.patPFPileUpFlags->size() == handler.recoPFCandidates->size());
        for (unsigned int i = 0; i < handler.recoPFCandidates->size(); ++i) {
            simple::Particle simple;
            const reco::PFCandidate& cand = handler.recoPFCandidates->at(i);
            const bool& isPU = handler.patPFPileUpFlags->at(i);
            simple_fill(cand, isPU, simple);
            recoPFCandidates.push_back(simple);
        }
        
        //______________________________________________________________________
        // patJets
        if (verbose)  std::cout << "compactify: Begin filling patJets..." << std::endl;
        patJets.clear();
        for (unsigned int i = 0; i < handler.patJets->size(); ++i) {
            const pat::Jet& jet = handler.patJets->at(i);
            simple::PFJet simple;
            simple_fill(jet, pfJetIDSelectionFunctor(jet, pfJetIDBitSet), simple);
            patJets.push_back(simple);
        }
        
        //______________________________________________________________________
        // patMET
        if (verbose)  std::cout << "compactify: Begin filling patMET..." << std::endl;
        simple_fill(handler.patMETs->at(0), patMET);
        
        //______________________________________________________________________
        // patTrackMET
        // ???
        
        //______________________________________________________________________
        // recoRhos
        if (verbose)  std::cout << "compactify: Begin filling recoRhos..." << std::endl;
        recoRho_kt6CaloJets = *(handler.recoRho_kt6CaloJets);
        recoRho_kt6PFJets = *(handler.recoRho_kt6PFJets);
        
        //______________________________________________________________________
        // Fill
        outtree->Fill();
        ++jevent;
    }
    
    outfile->cd();
    outtree->Write();
    outfile->Close();
    
    std::cout << ">>> number of events: " << jevent << "/" << ievent << std::endl;
    
    return 0;
}

