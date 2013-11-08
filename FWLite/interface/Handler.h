#ifndef HANDLER_H_
#define HANDLER_H_

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
// L1
//#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
//#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
// HLT+RECO
//#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
//#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
//#include "DataFormats/EgammaCandidates/interface/Electron.h"
//#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
//#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
//#include "DataFormats/Candidate/interface/CompositeCandidate.h"
//#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
//#include "DataFormats/TauReco/interface/CaloTau.h"
//#include "DataFormats/TauReco/interface/CaloTauFwd.h"
//#include "DataFormats/TauReco/interface/HLTTau.h"
//#include "DataFormats/TauReco/interface/HLTTauFwd.h"
//#include "DataFormats/TauReco/interface/PFTau.h"
//#include "DataFormats/TauReco/interface/PFTauFwd.h"
//#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/JetReco/interface/PFJetCollection.h"
// PAT
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
// Trigger results
#include "DataFormats/Common/interface/TriggerResults.h"
// GEN
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


class Handler {
  public:
    Handler(const edm::ParameterSet& iConfig, bool isData)
          : iConfig_(iConfig), isData_(isData) {}
    ~Handler() {}

    template<typename PROD>
    void getByLabel(const edm::EventBase& iEvent,
                    const std::string& parname, edm::Handle<PROD>& result) {
        if (!iConfig_.exists(parname)) {
            std::cout << "Handler: " << parname + " not found!" << std::endl;
            return;
        }
        const std::string& parvalue = iConfig_.getParameter<std::string>(parname);
        if (parvalue == "") {
            return;
        }
        bool b = iEvent.getByLabel(parvalue, result);
        if (!b || !result.isValid()) {
            std::cout << "Handler: " << parname + " is not retrieved successfully!" << std::endl;
        }
        return;
    }

    void get(const edm::EventBase& iEvent) {
        getByLabel(iEvent, "l1METs", l1METs);
        getByLabel(iEvent, "l1MHTs", l1MHTs);
        getByLabel(iEvent, "l1CentralJets", l1CentralJets);
        getByLabel(iEvent, "l1ForwardJets", l1ForwardJets);
        getByLabel(iEvent, "l1TauJets", l1TauJets);
        getByLabel(iEvent, "hltCaloJets", hltCaloJets);
        getByLabel(iEvent, "hltCaloJetIDPasseds", hltCaloJetIDPasseds);
        getByLabel(iEvent, "hltCaloJetL1Fasts", hltCaloJetL1Fasts);
        getByLabel(iEvent, "hltCaloMETs", hltCaloMETs);
        getByLabel(iEvent, "hltCaloMETCleans", hltCaloMETCleans);
        getByLabel(iEvent, "hltCaloMETJetIDCleans", hltCaloMETJetIDCleans);
        getByLabel(iEvent, "hltPFMETs", hltPFMETs);
        getByLabel(iEvent, "hltPFMETNoMus", hltPFMETNoMus);
        getByLabel(iEvent, "hltTrackMETs", hltTrackMETs);
        getByLabel(iEvent, "hltMuons", hltMuons);
        getByLabel(iEvent, "hltPFCandidates", hltPFCandidates);
        getByLabel(iEvent, "hltPFJets", hltPFJets);
        getByLabel(iEvent, "hltPFJetsNoPU", hltPFJetsNoPU);
        getByLabel(iEvent, "hltPFJetsL1FastL2L3", hltPFJetsL1FastL2L3);
        getByLabel(iEvent, "hltPFJetsL1FastL2L3NoPU", hltPFJetsL1FastL2L3NoPU);
        getByLabel(iEvent, "hltPFRecTracks", hltPFRecTracks);
        getByLabel(iEvent, "hltVertices", hltVertices);
        getByLabel(iEvent, "hltRho_kt6CaloJets", hltRho_kt6CaloJets);
        getByLabel(iEvent, "hltRho_kt6PFJets", hltRho_kt6PFJets);
        getByLabel(iEvent, "recoTracks", recoTracks);
        getByLabel(iEvent, "recoCaloJets", recoCaloJets);
        getByLabel(iEvent, "recoCaloMETs", recoCaloMETs);
        getByLabel(iEvent, "recoPFCandidates", recoPFCandidates);
        getByLabel(iEvent, "recoPFJets", recoPFJets);
        getByLabel(iEvent, "recoPFMETs", recoPFMETs);
        getByLabel(iEvent, "recoVertices", recoVertices);
        getByLabel(iEvent, "recoGoodVertices", recoGoodVertices);
        getByLabel(iEvent, "recoRho_kt6CaloJets", recoRho_kt6CaloJets);
        getByLabel(iEvent, "recoRho_kt6PFJets", recoRho_kt6PFJets);
        getByLabel(iEvent, "patElectrons", patElectrons);
        getByLabel(iEvent, "patJets", patJets);
        getByLabel(iEvent, "patMETs", patMETs);
        getByLabel(iEvent, "patMETTypeIs", patMETTypeIs);
        getByLabel(iEvent, "patMuons", patMuons);
        getByLabel(iEvent, "patPhotons", patPhotons);
        getByLabel(iEvent, "patTaus", patTaus);
        getByLabel(iEvent, "triggerResults", triggerResults);
        if (!isData_) {
            getByLabel(iEvent, "genJets", genJets);
            getByLabel(iEvent, "genMETs", genMETs);
            getByLabel(iEvent, "genParticles", genParticles);
            getByLabel(iEvent, "genEventInfo", genEventInfo);
            getByLabel(iEvent, "simPileupInfo", simPileupInfo);
        }
        getByLabel(iEvent, "hltPFPileUpFlags", hltPFPileUpFlags);
        getByLabel(iEvent, "patPFPileUpFlags", patPFPileUpFlags);
        getByLabel(iEvent, "hltCaloJetIDs", hltCaloJetIDs);
    }

    // L1
    edm::Handle<std::vector<l1extra::L1EtMissParticle> > l1METs;
    edm::Handle<std::vector<l1extra::L1EtMissParticle> > l1MHTs;
    edm::Handle<std::vector<l1extra::L1JetParticle   > > l1CentralJets;
    edm::Handle<std::vector<l1extra::L1JetParticle   > > l1ForwardJets;
    edm::Handle<std::vector<l1extra::L1JetParticle   > > l1TauJets;

    // HLT
    edm::Handle<std::vector<reco::CaloJet      > > hltCaloJets;
    edm::Handle<std::vector<reco::CaloJet      > > hltCaloJetIDPasseds;
    edm::Handle<std::vector<reco::CaloJet      > > hltCaloJetL1Fasts;
    edm::Handle<std::vector<reco::CaloMET      > > hltCaloMETs;
    edm::Handle<std::vector<reco::CaloMET      > > hltCaloMETCleans;
    edm::Handle<std::vector<reco::CaloMET      > > hltCaloMETJetIDCleans;
    edm::Handle<std::vector<reco::MET          > > hltPFMETs;
    edm::Handle<std::vector<reco::MET          > > hltPFMETNoMus;
    edm::Handle<std::vector<reco::MET          > > hltTrackMETs;
    edm::Handle<std::vector<reco::Muon         > > hltMuons;
    edm::Handle<std::vector<reco::PFCandidate  > > hltPFCandidates;
    edm::Handle<std::vector<reco::PFJet        > > hltPFJets;
    edm::Handle<std::vector<reco::PFJet        > > hltPFJetsNoPU;
    edm::Handle<std::vector<reco::PFJet        > > hltPFJetsL1FastL2L3;
    edm::Handle<std::vector<reco::PFJet        > > hltPFJetsL1FastL2L3NoPU;
    edm::Handle<std::vector<reco::PFRecTrack   > > hltPFRecTracks;
    edm::Handle<std::vector<reco::Vertex       > > hltVertices;
    edm::Handle<double                           > hltRho_kt6CaloJets;
    edm::Handle<double                           > hltRho_kt6PFJets;

    // RECO
    edm::Handle<std::vector<reco::Track        > > recoTracks;
    edm::Handle<std::vector<reco::CaloJet      > > recoCaloJets;
    edm::Handle<std::vector<reco::CaloMET      > > recoCaloMETs;
    edm::Handle<std::vector<reco::PFCandidate  > > recoPFCandidates;
    edm::Handle<std::vector<reco::PFJet        > > recoPFJets;
    edm::Handle<std::vector<reco::PFMET        > > recoPFMETs;
    edm::Handle<std::vector<reco::Vertex       > > recoVertices;
    edm::Handle<std::vector<reco::Vertex       > > recoGoodVertices;
    edm::Handle<double                           > recoRho_kt6CaloJets;
    edm::Handle<double                           > recoRho_kt6PFJets;

    // PAT
    edm::Handle<std::vector<pat::Electron      > > patElectrons;
    edm::Handle<std::vector<pat::Jet           > > patJets;
    edm::Handle<std::vector<pat::MET           > > patMETs;
    edm::Handle<std::vector<pat::MET           > > patMETTypeIs;
    edm::Handle<std::vector<pat::Muon          > > patMuons;
    edm::Handle<std::vector<pat::Photon        > > patPhotons;
    edm::Handle<std::vector<pat::Tau           > > patTaus;

    // Trigger results
    edm::Handle<edm::TriggerResults              > triggerResults;

    // GEN/SIM
    edm::Handle<std::vector<reco::GenJet       > > genJets;
    edm::Handle<std::vector<reco::GenMET       > > genMETs;
    edm::Handle<std::vector<reco::GenParticle  > > genParticles;
    edm::Handle<GenEventInfoProduct              > genEventInfo;
    edm::Handle<std::vector<PileupSummaryInfo  > > simPileupInfo;

    // User
    edm::Handle<std::vector<bool               > > hltPFPileUpFlags;
    edm::Handle<std::vector<bool               > > patPFPileUpFlags;
    edm::Handle<edm::ValueMap<reco::JetID      > > hltCaloJetIDs;

  private:
    edm::ParameterSet iConfig_;
    bool isData_;
};

#endif  // HANDLER_H_
