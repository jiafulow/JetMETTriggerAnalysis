#include <string>
#include "TFile.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TTree.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
// CMSSW: DataFormats
#include "DataFormats/Math/interface/deltaR.h"
// CMSSW: FWLite
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
// CMSSW: JetMETTriggerAnalysis
#include "JetMETTriggerAnalysis/FWLite/interface/SimpleCandidate.h"
#include "JetMETTriggerAnalysis/FWLite/interface/Handler.h"
#include "JetMETTriggerAnalysis/FWLite/interface/JSONFilter.h"
#include "JetMETTriggerAnalysis/FWLite/interface/JetIDHelper.h"
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


//______________________________________________________________________________
// Fill
template<typename T>
void simple_fill(const T& cand, simple::XYVector& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pt     = cand.p4().pt();
    simple.phi    = cand.p4().phi();
    return;
}

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
void simple_fill(const T& cand, bool jetID, const reco::JetID& clsJetID, simple::CaloJet& simple) {
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
    simple.emf    = cand.emEnergyFraction();
    simple.fHPD   = clsJetID.fHPD;
    simple.n90Hits= clsJetID.n90Hits;
    return;
}

template<typename T>
void simple_fill(const T& cand, bool jetID, simple::PFJet& simple) {
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
    simple.chf    = cand.chargedHadronEnergyFraction();
    simple.nhf    = cand.neutralHadronEnergyFraction() + cand.HFHadronEnergyFraction();
    simple.cef    = cand.chargedEmEnergyFraction();
    simple.nef    = cand.neutralEmEnergyFraction();
    simple.nch    = cand.chargedMultiplicity();
    simple.ntot   = cand.numberOfDaughters();
    simple.csv    = (patJet != 0 ? patJet->bDiscriminator("combinedSecondaryVertexBJetTags") : -999);
    return;
}

template<typename T>
void simple_fill(const T& cand, simple::MET& simple) {
    simple.px     = cand.p4().px();
    simple.py     = cand.p4().py();
    simple.pt     = cand.p4().pt();
    simple.phi    = cand.p4().phi();
    simple.sumEt  = cand.sumEt();
    return;
}

void simple_fill(const edm::EventBase& iEvent, unsigned int nPV,
                 unsigned int nGoodPV, unsigned int nTruePV, bool json,
                 simple::Event& simple) {
    simple.run    = iEvent.id().run();
    simple.lumi   = iEvent.id().luminosityBlock();
    simple.event  = iEvent.id().event();
    simple.nPV    = nPV;
    simple.nGoodPV= nGoodPV;
    simple.json   = json;
    return;
}


//______________________________________________________________________________
// Functions
template<typename T>
unsigned int count_jets(const std::vector<T>& jets, double ptmin, double etamax) {
    unsigned int count = 0;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        if (jets.at(j).pt() > ptmin && fabs(jets.at(j).eta()) < etamax) {
            ++count;
        }
    }
    return count;
}

template<typename T>
double eval_ht(const std::vector<T>& jets, double ptmin, double etamax) {
    double ht = 0.0;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        const T& jet = jets.at(j);
        if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {
            ht += jets.at(j).pt();
        }
    }
    return ht;
}

template<typename T>
double eval_mindphi(double metphi, const std::vector<T>& jets, double ptmin, double etamax) {
    double mindphi = M_PI;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        if (jets.at(j).pt() > ptmin && fabs(jets.at(j).eta()) < etamax) {
            double dphi = fabs(reco::deltaPhi(jets.at(j).phi(), metphi));
            if (mindphi > dphi)
                mindphi = dphi;
        }
    }
    return mindphi;
}

template<typename T>
reco::Candidate::LorentzVector build_met(const std::vector<T>& particles, double ptmin, double etamax) {
    reco::Candidate::LorentzVector met(0,0,0,0);
    for (unsigned int j = 0; j < particles.size(); ++j) {
        const T& particle = particles.at(j);
        if (particle.pt() > ptmin && fabs(particle.eta()) < etamax) {
            met -= particle.p4();
        }
    }
    return met;
}

template<typename T>
bool isEqualEtaPhi(const T& j1, const T& j2) {
    double epsilon = 1e-3;
    bool equal = (fabs(j1.eta() - j2.eta()) < epsilon &&
                  fabs(reco::deltaPhi(j1.phi(), j2.phi())) < epsilon);
    return equal;
}


//______________________________________________________________________________
// Compactify
int main(int argc, char *argv[]) {
    // Load FWLite library
    gSystem->Load("libFWCoreFWLite");
    AutoLibraryLoader::enable();
    //gSystem->Load("libDataFormatsFWLite");

    // Generate libraries
    //gInterpreter->GenerateDictionary("simple::MET", "JetMETTriggerAnalysis/FWLite/interface/simple::Candidate.h");

    // Load config file
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " [parameters.py]" << std::endl;
        return 0;
    }
    PythonProcessDesc builder(argv[1], argc, argv);
    std::cout << "Read config: " << argv[1] << std::endl;

    const edm::ParameterSet cfg = *builder.parameterSet();
    const edm::ParameterSet& inputpset    = cfg.getParameter<edm::ParameterSet>("input");
    const edm::ParameterSet& outputpset   = cfg.getParameter<edm::ParameterSet>("output");
    const edm::ParameterSet& analyzerpset = cfg.getParameter<edm::ParameterSet>("analyzer");
    const edm::ParameterSet& handlerpset  = cfg.getParameter<edm::ParameterSet>("handler");
    const edm::ParameterSet& lumicalcpset  = cfg.getParameter<edm::ParameterSet>("lumicalc");

    // Read from input pset
    const std::vector<std::string>& infilenames_all = inputpset.getParameter<std::vector<std::string> >("fileNames");
    const int njobs = inputpset.getParameter<int>("njobs");
    const int jobid = inputpset.getParameter<int>("jobid");
    const int maxEvents = inputpset.getParameter<int>("maxEvents");
    const int runMin = inputpset.getParameter<int>("runMin");
    const int runMax = inputpset.getParameter<int>("runMax");
    const int skipEvents = inputpset.getParameter<int>("skipEvents");
    const int reportEvery = inputpset.getParameter<int>("reportEvery");
    std::vector<edm::LuminosityBlockRange> lumisToProcess;
    if (inputpset.exists("lumisToProcess") ) {
        const std::vector<edm::LuminosityBlockRange>& lumisTemp =
            inputpset.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> >("lumisToProcess");
        lumisToProcess.resize( lumisTemp.size() );
        copy(lumisTemp.begin(), lumisTemp.end(), lumisToProcess.begin() );
    }
    int ibeginjob = 0;
    int iendjob = infilenames_all.size();
    if (njobs != 1) {
        ibeginjob = jobid * (infilenames_all.size() / njobs);
        iendjob = (jobid+1) * (infilenames_all.size() / njobs);
        if (jobid == njobs - 1)  iendjob = infilenames_all.size();
    }
    const std::vector<std::string> infilenames(infilenames_all.begin()+ibeginjob,
                                               infilenames_all.begin()+iendjob);

    // Read from output pset
    TString outfilename = outputpset.getParameter<std::string>("fileName");
    if (njobs != 1) {
        outfilename.ReplaceAll(".root", Form("_%i.root", jobid));
    }

    // Read from analyzer pset
    const std::vector<std::string>& triggers = analyzerpset.getParameter<std::vector<std::string> >("triggers");
    const std::vector<std::string>& metfilters = analyzerpset.getParameter<std::vector<std::string> >("metfilters");
    const std::vector<std::string>& optmetfilters = analyzerpset.getParameter<std::vector<std::string> >("optmetfilters");
    const edm::ParameterSet& hltJetIDParams = analyzerpset.getParameter<edm::ParameterSet>("hltJetID");
    const edm::ParameterSet& jetIDParams = analyzerpset.getParameter<edm::ParameterSet>("jetID");
    const double hltCaloJetPtMin = analyzerpset.getParameter<double>("hltCaloJetPtMin");
    const double hltCaloJetEtaMax = analyzerpset.getParameter<double>("hltCaloJetEtaMax");
    const double hltPFJetPtMin = analyzerpset.getParameter<double>("hltPFJetPtMin");
    const double hltPFJetEtaMax = analyzerpset.getParameter<double>("hltPFJetEtaMax");
    const double patJetPtMin = analyzerpset.getParameter<double>("patJetPtMin");
    const double patJetEtaMax = analyzerpset.getParameter<double>("patJetEtaMax");
    const bool isData = analyzerpset.getParameter<bool>("isData");
    const bool verbose = analyzerpset.getParameter<bool>("verbose");
    HLTJetIDHelper hltJetIDHelper(hltJetIDParams);
    JetIDHelper jetIDHelper(jetIDParams);

    // Read from handler pset
    Handler handler(handlerpset, isData);

    // Read from lumicalc pset
    const std::vector<unsigned long long>& lumiA = lumicalcpset.getParameter<std::vector<unsigned long long> >("lumiA");
    const std::vector<unsigned long long>& lumiB = lumicalcpset.getParameter<std::vector<unsigned long long> >("lumiB");
    const std::vector<unsigned long long>& lumiC = lumicalcpset.getParameter<std::vector<unsigned long long> >("lumiC");


    //__________________________________________________________________________
    // Prepare output tree
    simple::Event simpleEvent;
    float weightGenEvent, weightPU, weightPdf;
    bool triggerFlags[99]; int nTriggerFlags = triggers.size() + 1;  // last one is "OR" combination
    bool metfilterFlags[99]; int nMetfilterFlags = metfilters.size() + 1;  // last one is "AND" combination
    bool optmetfilterFlags[99]; int nOptmetfilterFlags = optmetfilters.size() + 1;  // last one is "AND" combination
    // HLT
    std::vector<simple::Particle> hltPFCandidates;
    std::vector<simple::CaloJet> hltCaloJets;
    std::vector<simple::Jet> hltCaloJetsL1Fast;  // JetID info is not saved for this collection
    std::vector<simple::PFJet> hltPFJets;
    std::vector<simple::PFJet> hltPFJetsNoPU;
    std::vector<simple::PFJet> hltPFJetsL1FastL2L3;
    std::vector<simple::PFJet> hltPFJetsL1FastL2L3NoPU;
    simple::MET hltCaloMET;
    simple::MET hltCaloMETClean;
    simple::MET hltCaloMETCleanUsingJetID;
    simple::MET hltPFMET;
    simple::MET hltPFMETNoMu;
    simple::MET hltTrackMET;
    simple::MET hltHTMHT;
    simple::MET hltPFHTMHT;
    simple::MET hltPFHTMHTNoPU;
    double hltRho_kt6CaloJets, hltRho_kt6PFJets;
    unsigned int lumilevel;
    // non-HLT
    std::vector<simple::Particle> recoPFCandidates;
    std::vector<simple::PFJet> patJets;
    simple::MET recoPFMET;
    simple::MET recoPFMETT1;
    simple::MET recoPFMETT0T1;
    simple::MET recoPFMETMVA;
    simple::MET recoPFMETNoPU;
    simple::MET patMET;
    simple::MET patMPT;
    double recoRho_kt6CaloJets, recoRho_kt6PFJets;
    unsigned int topo_patJets;


    TFile* outfile = new TFile(outfilename, "RECREATE");
    TTree* outtree = new TTree("Events", "Events");
    outtree->Branch("event", &simpleEvent);
    outtree->Branch("weightGenEvent", &weightGenEvent);
    outtree->Branch("weightPU", &weightPU);  //FIXME not yet set
    outtree->Branch("weightPdf", &weightPdf);  //FIXME not yet set
    outtree->Branch("triggerFlags", triggerFlags, Form("triggerFlags[%i]/b", nTriggerFlags));
    outtree->Branch("metfilterFlags", metfilterFlags, Form("metfilterFlags[%i]/b", nMetfilterFlags));
    outtree->Branch("optmetfilterFlags", optmetfilterFlags, Form("optmetfilterFlags[%i]/b", nOptmetfilterFlags));
    // HLT
    outtree->Branch("hltPFCandidates", &hltPFCandidates);
    outtree->Branch("hltCaloJets", &hltCaloJets);
    outtree->Branch("hltCaloJetsL1Fast", & hltCaloJetsL1Fast);
    outtree->Branch("hltPFJets", &hltPFJets);
    outtree->Branch("hltPFJetsNoPU", &hltPFJetsNoPU);
    outtree->Branch("hltPFJetsL1FastL2L3", &hltPFJetsL1FastL2L3);
    outtree->Branch("hltPFJetsL1FastL2L3NoPU", &hltPFJetsL1FastL2L3NoPU);
    outtree->Branch("hltCaloMET", &hltCaloMET);
    outtree->Branch("hltCaloMETClean", &hltCaloMETClean);
    outtree->Branch("hltCaloMETCleanUsingJetID", &hltCaloMETCleanUsingJetID);
    outtree->Branch("hltPFMET", &hltPFMET);
    outtree->Branch("hltPFMETNoMu", &hltPFMETNoMu);
    outtree->Branch("hltTrackMET", &hltTrackMET);
    outtree->Branch("hltHTMHT", &hltHTMHT);
    outtree->Branch("hltPFHTMHT", &hltPFHTMHT);
    outtree->Branch("hltPFHTMHTNoPU", &hltPFHTMHTNoPU);
    outtree->Branch("hltRho_kt6CaloJets", &hltRho_kt6CaloJets);
    outtree->Branch("hltRho_kt6PFJets", &hltRho_kt6PFJets);
    outtree->Branch("lumilevel", &lumilevel);
    // non-HLT
    outtree->Branch("recoPFCandidates", &recoPFCandidates);
    outtree->Branch("patJets", &patJets);
    outtree->Branch("recoPFMET", &recoPFMET);
    outtree->Branch("recoPFMETT1", &recoPFMETT1);
    outtree->Branch("recoPFMETT0T1", &recoPFMETT0T1);
    outtree->Branch("recoPFMETMVA", &recoPFMETMVA);
    outtree->Branch("recoPFMETNoPU", &recoPFMETNoPU);
    outtree->Branch("patMET", &patMET);
    outtree->Branch("patMPT", &patMPT);
    outtree->Branch("recoRho_kt6CaloJets", &recoRho_kt6CaloJets);
    outtree->Branch("recoRho_kt6PFJets", &recoRho_kt6PFJets);
    outtree->Branch("topo_patJets", &topo_patJets);


    //__________________________________________________________________________
    // Loop over events
    std::cout << "<<< number of input files: " << infilenames.size() << std::endl;
    bool verboseInput = true;
    if (verbose || verboseInput) {
        for (unsigned int i = 0; i < infilenames.size(); ++i) {
            std::cout << "    --- " << infilenames.at(i) << std::endl;
        }
    }
    fwlite::ChainEvent ev(infilenames);
    int ievent = 0;
    int jevent = 0;
    unsigned long long runlumi_old = 0;
    unsigned int lumilevel_old = 0;
    for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievent) {
        // Skip events
        if (ievent < skipEvents)  continue;
        if (runMin > 0 && (int)ev.id().run() < runMin)  continue;
        if (runMax > 0 && (int)ev.id().run() > runMax)  continue;
        // Stop event loop
        if (maxEvents > 0 && (ievent >= (maxEvents + skipEvents)) )  break;
        if (ievent % reportEvery == 0)  std::cout << "--- ... Processing event: " << ievent << std::endl;

        const edm::EventBase& eventbase = ev;
        handler.get(eventbase);


        //______________________________________________________________________
        // Event info
        if (verbose)  std::cout << "compactify: Begin filling event info..." << std::endl;
        bool goodjson = jsonContainsEvent(lumisToProcess, ev);
        //unsigned int nPV = handler.recoVertices->size();
        unsigned int nPV = 999;
        unsigned int nGoodPV = handler.recoGoodVertices->size();
        unsigned int nTruePV = (isData ? 999 : handler.simPileupInfo->front().getTrueNumInteractions());
        simple_fill(ev, nPV, nGoodPV, nTruePV, goodjson, simpleEvent);
        weightGenEvent = (isData ? 1.0 : handler.genEventInfo->weight());


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

        bool optmetfilterFlagsAND = true;
        for (unsigned int i = 0; i < optmetfilters.size(); ++i) {
            bool accept = resultsByNamePAT.accept(optmetfilters.at(i));
            optmetfilterFlags[i] = accept;
            if (!accept)  optmetfilterFlagsAND = false;
        }
        optmetfilterFlags[optmetfilters.size()] = optmetfilterFlagsAND;

        //______________________________________________________________________
        // hltPFCandidates
        //if (verbose)  std::cout << "compactify: Begin filling hltPFCandidates..." << std::endl;
        //hltPFCandidates.clear();
        //assert(handler.hltPFPileUpFlags->size() == handler.hltPFCandidates->size());
        //for (unsigned int i = 0; i < handler.hltPFCandidates->size(); ++i) {
        //    simple::Particle simple;
        //    const reco::PFCandidate& cand = handler.hltPFCandidates->at(i);
        //    const bool& isPU = handler.hltPFPileUpFlags->at(i);
        //    simple_fill(cand, isPU, simple);
        //    hltPFCandidates.push_back(simple);
        //}

        //______________________________________________________________________
        // hltCaloJets
        if (verbose)  std::cout << "compactify: Begin filling hltCaloJets..." << std::endl;
        hltCaloJets.clear();
        assert(handler.hltCaloJets->size() == handler.hltCaloJetIDs->size());
        for (unsigned int i = 0; i < handler.hltCaloJets->size(); ++i) {
            const reco::CaloJet& jet = handler.hltCaloJets->at(i);
            if (jet.pt() < hltCaloJetPtMin || fabs(jet.eta()) > hltCaloJetEtaMax)  continue;
            // These 2 lines don't work in FWLite
            //const reco::CaloJetRef jetref(handler.hltCaloJets, i);
            //const reco::JetID jetID_ = (*handler.hltCaloJetIDs)[jetref];
            assert(handler.hltCaloJetIDs->idSize() == 1);
            edm::ProductID productid = handler.hltCaloJetIDs->ids().front().first;
            const reco::JetID clsJetID = handler.hltCaloJetIDs->get(productid, i);
            simple::CaloJet simple;
            bool jetID = hltJetIDHelper(jet, clsJetID);
            simple_fill(jet, jetID, clsJetID, simple);
            hltCaloJets.push_back(simple);
        }

        hltCaloJetsL1Fast.clear();
        for (unsigned int i = 0; i < handler.hltCaloJetsL1Fast->size(); ++i) {
            const reco::CaloJet& jet = handler.hltCaloJetsL1Fast->at(i);
            if (jet.pt() < hltCaloJetPtMin || fabs(jet.eta()) > hltCaloJetEtaMax)  continue;
            simple::Jet simple;
            simple_fill(jet, true, simple);  // already passed JetID
            hltCaloJetsL1Fast.push_back(simple);
        }

        //______________________________________________________________________
        // hltPFJets
        if (verbose)  std::cout << "compactify: Begin filling hltPFJets..." << std::endl;
        hltPFJets.clear();
        for (unsigned int i = 0; i < handler.hltPFJets->size(); ++i) {
            const reco::PFJet& jet = handler.hltPFJets->at(i);
            if (jet.pt() < hltPFJetPtMin || fabs(jet.eta()) > hltPFJetEtaMax)  continue;
            simple::PFJet simple;
            bool jetID = hltJetIDHelper(jet);
            simple_fill(jet, jetID, simple);
            hltPFJets.push_back(simple);
        }

        hltPFJetsNoPU.clear();
        for (unsigned int i = 0; i < handler.hltPFJetsNoPU->size(); ++i) {
            const reco::PFJet& jet = handler.hltPFJetsNoPU->at(i);
            if (jet.pt() < hltPFJetPtMin || fabs(jet.eta()) > hltPFJetEtaMax)  continue;
            simple::PFJet simple;
            bool jetID = hltJetIDHelper(jet);
            simple_fill(jet, jetID, simple);
            hltPFJetsNoPU.push_back(simple);
        }

        hltPFJetsL1FastL2L3.clear();
        assert(handler.hltPFJets->size() == handler.hltPFJetsL1FastL2L3->size());
        for (unsigned int i = 0; i < handler.hltPFJetsL1FastL2L3->size(); ++i) {
            const reco::PFJet& jet = handler.hltPFJetsL1FastL2L3->at(i);
            if (jet.pt() < hltPFJetPtMin || fabs(jet.eta()) > hltPFJetEtaMax)  continue;
            simple::PFJet simple;
            bool jetID = hltJetIDHelper(jet);
            simple_fill(jet, jetID, simple);
            hltPFJetsL1FastL2L3.push_back(simple);
        }

        hltPFJetsL1FastL2L3NoPU.clear();
        assert(handler.hltPFJetsNoPU->size() == handler.hltPFJetsL1FastL2L3NoPU->size());
        for (unsigned int i = 0; i < handler.hltPFJetsL1FastL2L3NoPU->size(); ++i) {
            const reco::PFJet& jet = handler.hltPFJetsL1FastL2L3NoPU->at(i);
            if (jet.pt() < hltPFJetPtMin || fabs(jet.eta()) > hltPFJetEtaMax)  continue;
            simple::PFJet simple;
            bool jetID = hltJetIDHelper(jet);
            simple_fill(jet, jetID, simple);
            hltPFJetsL1FastL2L3.push_back(simple);
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
        // hltCaloMETCleanUsingJetID
        if (verbose)  std::cout << "compactify: Begin filling hltCaloMETCleanUsingJetID..." << std::endl;
        simple_fill(handler.hltCaloMETCleansUsingJetID->at(0), hltCaloMETCleanUsingJetID);

        //______________________________________________________________________
        // hltPFMET
        if (verbose)  std::cout << "compactify: Begin filling hltPFMET..." << std::endl;
        simple_fill(handler.hltPFMETs->at(0), hltPFMET);

        //______________________________________________________________________
        // hltPFMETNoMu
        if (verbose)  std::cout << "compactify: Begin filling hltPFMETNoMu..." << std::endl;
        simple_fill(handler.hltPFMETsNoMu->at(0), hltPFMETNoMu);

        //______________________________________________________________________
        // hltTrackMET
        if (verbose)  std::cout << "compactify: Begin filling hltTrackMET..." << std::endl;
        simple_fill(handler.hltTrackMETs->at(0), hltTrackMET);

        //______________________________________________________________________
        // hltHTMHT
        if (verbose)  std::cout << "compactify: Begin filling hltHTMHT..." << std::endl;
        simple_fill(handler.hltHTMHTs->at(0), hltHTMHT);

        //______________________________________________________________________
        // hltPFHTMHT
        //if (verbose)  std::cout << "compactify: Begin filling hltPFHTMHT..." << std::endl;
        //simple_fill(handler.hltPFHTMHTs->at(0), hltPFHTMHT);

        //______________________________________________________________________
        // hltPFHTMHTNoPU
        if (verbose)  std::cout << "compactify: Begin filling hltPFHTMHTNoPU..." << std::endl;
        simple_fill(handler.hltPFHTMHTsNoPU->at(0), hltPFHTMHTNoPU);

        //______________________________________________________________________
        // hltRhos
        if (verbose)  std::cout << "compactify: Begin filling hltRhos..." << std::endl;
        hltRho_kt6CaloJets = *(handler.hltRho_kt6CaloJets);
        hltRho_kt6PFJets = *(handler.hltRho_kt6PFJets);

        //______________________________________________________________________
        // lumilevel
        if (verbose)  std::cout << "compactify: Begin filling lumilevel..." << std::endl;
        lumilevel = 0;
        unsigned long long runlumi = ev.id().run();
        runlumi *= 100000;
        runlumi += ev.id().luminosityBlock();
        if (runlumi == runlumi_old) {
            lumilevel = lumilevel_old;
        }

        if (ev.id().run() != 1) {  // not MC
            std::vector<unsigned long long>::const_iterator it;
            if (lumilevel == 0) {
                it = std::find(lumiC.begin(), lumiC.end(), runlumi);
                if (it != lumiC.end())  lumilevel = 3;
            }
            if (lumilevel == 0) {
                it = std::find(lumiB.begin(), lumiB.end(), runlumi);
                if (it != lumiB.end())  lumilevel = 2;
            }
            if (lumilevel == 0) {
                it = std::find(lumiA.begin(), lumiA.end(), runlumi);
                if (it != lumiA.end())  lumilevel = 1;
            }
            if (lumilevel == 0)
                std::cout << "ERROR: runlumi " << runlumi << " not found!" << std::endl;
            else if (lumilevel_old == 0 || runlumi_old != runlumi) {
                runlumi_old = runlumi;
                lumilevel_old = lumilevel;
            }
        }

        //______________________________________________________________________
        // recoPFCandidates
        //if (verbose)  std::cout << "compactify: Begin filling recoPFCandidates..." << std::endl;
        //recoPFCandidates.clear();
        //assert(handler.patPFPileUpFlags->size() == handler.recoPFCandidates->size());
        //for (unsigned int i = 0; i < handler.recoPFCandidates->size(); ++i) {
        //    simple::Particle simple;
        //    const reco::PFCandidate& cand = handler.recoPFCandidates->at(i);
        //    const bool& isPU = handler.patPFPileUpFlags->at(i);
        //    simple_fill(cand, isPU, simple);
        //    recoPFCandidates.push_back(simple);
        //}

        //______________________________________________________________________
        // patJets
        if (verbose)  std::cout << "compactify: Begin filling patJets..." << std::endl;
        patJets.clear();
        for (unsigned int i = 0; i < handler.patJets->size(); ++i) {
            const pat::Jet& jet = handler.patJets->at(i);
            if (jet.pt() < patJetPtMin || fabs(jet.eta()) > patJetEtaMax)  continue;
            simple::PFJet simple;
            bool jetID = jetIDHelper(jet);
            simple_fill(jet, jetID, simple);
            patJets.push_back(simple);
        }

        //______________________________________________________________________
        // recoPFMET
        if (verbose)  std::cout << "compactify: Begin filling recoPFMET..." << std::endl;
        simple_fill(handler.recoPFMETs->at(0), recoPFMET);
        simple_fill(handler.recoPFMETT1s->at(0), recoPFMETT1);
        simple_fill(handler.recoPFMETT0T1s->at(0), recoPFMETT0T1);
        simple_fill(handler.recoPFMETMVAs->at(0), recoPFMETMVA);
        simple_fill(handler.recoPFMETNoPUs->at(0), recoPFMETNoPU);

        //______________________________________________________________________
        // patMET
        if (verbose)  std::cout << "compactify: Begin filling patMET..." << std::endl;
        simple_fill(handler.patMETs->at(0), patMET);

        reco::PFCandidate::LorentzVector patMPT_p4(0,0,0,0);
        double patMPT_sumEt = 0.;
        for (unsigned int i = 0; i < handler.recoPFCandidates->size(); ++i) {
            const reco::PFCandidate& cand = handler.recoPFCandidates->at(i);
            const bool& isPU = handler.patPFPileUpFlags->at(i);
            if (!isPU && cand.charge() != 0) {
                patMPT_p4 -= cand.p4();
                patMPT_sumEt += cand.pt();
            }
        }
        patMPT.px     = patMPT_p4.px();
        patMPT.py     = patMPT_p4.py();
        patMPT.pt     = patMPT_p4.pt();
        patMPT.phi    = patMPT_p4.phi();
        patMPT.sumEt  = patMPT_sumEt;

        //______________________________________________________________________
        // recoRhos
        if (verbose)  std::cout << "compactify: Begin filling recoRhos..." << std::endl;
        recoRho_kt6CaloJets = *(handler.recoRho_kt6CaloJets);
        recoRho_kt6PFJets = *(handler.recoRho_kt6PFJets);

        //______________________________________________________________________
        // topo_patJets
        if (verbose)  std::cout << "compactify: Begin filling topo_patJets..." << std::endl;
        topo_patJets = 0;
        for (unsigned int i = 0; i < handler.patJets->size(); ++i) {
            const pat::Jet& jet = handler.patJets->at(i);
            if (jet.pt() < 30 || fabs(jet.eta()) > 2.5)  continue;
            topo_patJets += 1;
        }
        if (topo_patJets > 2)
            topo_patJets = 2;

        // Check for HT
        if (handler.patJets->size() >= 1) {
            double ht_patJets = 0.;
            for (unsigned int i = 0; i < handler.patJets->size(); ++i) {
                const pat::Jet& jet = handler.patJets->at(i);
                if (jet.pt() < 40 || fabs(jet.eta()) > 3)  continue;
                ht_patJets += jet.pt();
            }
            if (ht_patJets > 300)
                topo_patJets = 3;
        }

        // Check for VBF
        if (handler.patJets->size() >= 2) {
            for (unsigned int i = 0; i < handler.patJets->size()-1; ++i) {
                const pat::Jet& jet1 = handler.patJets->at(i);
                if (jet1.pt() < 30 || fabs(jet1.eta()) > 4.7)  continue;

                for (unsigned int j = i+1; j < handler.patJets->size(); ++j) {
                    const pat::Jet& jet2 = handler.patJets->at(j);
                    if (jet2.pt() < 30 || fabs(jet2.eta()) > 4.7)  continue;

                    if ((jet1.p4() + jet2.p4()).mass() > 600 && fabs(jet1.eta() - jet2.eta()) > 3.5)
                        topo_patJets = 4;

                    if (topo_patJets == 4)  break;
                }
                if (topo_patJets == 4)  break;
            }
        }


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

