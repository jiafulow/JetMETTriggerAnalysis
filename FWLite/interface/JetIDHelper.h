#ifndef JETIDHELPER_H_
#define JETIDHELPER_H_

#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/JetReco/interface/Jet.h"
//#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/JetID.h"


class JetIDHelper {
  public:
    JetIDHelper(const edm::ParameterSet& iConfig)
          : caloJetIDSelectionFunctor(iConfig.getParameter<edm::ParameterSet>("caloJetID") ),
            pfJetIDSelectionFunctor(iConfig.getParameter<edm::ParameterSet>("pfJetID") ) {
        caloJetIDBitSet = caloJetIDSelectionFunctor.getBitTemplate();
        pfJetIDBitSet = pfJetIDSelectionFunctor.getBitTemplate();
    }
    ~JetIDHelper() {}
    
    bool operator()(const pat::Jet& jet) {
        if (jet.isCaloJet()) {
            return caloJetIDSelectionFunctor(jet, caloJetIDBitSet);
        } else if (jet.isPFJet()) {
            return pfJetIDSelectionFunctor(jet, pfJetIDBitSet);
        } else {
            std::cout << "JetIDHelper: Unknown jet type!" << std::endl;
            return true;
        }
    }
    
    // JetID selector
    JetIDSelectionFunctor   caloJetIDSelectionFunctor;
    PFJetIDSelectionFunctor pfJetIDSelectionFunctor;
    pat::strbitset          caloJetIDBitSet;
    pat::strbitset          pfJetIDBitSet;
};


class HLTJetIDHelper {
  public:
    HLTJetIDHelper(const edm::ParameterSet& iConfig) {
        const edm::ParameterSet& caloJetIDConfig = iConfig.getParameter<edm::ParameterSet>("hltCaloJetID");
        min_EMF_     = caloJetIDConfig.getParameter<double>("min_EMF");
        max_EMF_     = caloJetIDConfig.getParameter<double>("max_EMF");
        min_N90_     = caloJetIDConfig.getParameter<int>("min_N90");
        min_N90hits_ = caloJetIDConfig.getParameter<int>("min_N90hits");
        
        const edm::ParameterSet& pfJetIDConfig = iConfig.getParameter<edm::ParameterSet>("hltPFJetID");
        min_CEEF_ = pfJetIDConfig.getParameter<double>("min_CEEF");
        max_CEEF_ = pfJetIDConfig.getParameter<double>("max_CEEF");
        min_NEEF_ = pfJetIDConfig.getParameter<double>("min_NEEF");
        max_NEEF_ = pfJetIDConfig.getParameter<double>("max_NEEF");
        min_CHEF_ = pfJetIDConfig.getParameter<double>("min_CHEF");
        max_CHEF_ = pfJetIDConfig.getParameter<double>("max_CHEF");
        min_NHEF_ = pfJetIDConfig.getParameter<double>("min_NHEF");
        max_NHEF_ = pfJetIDConfig.getParameter<double>("max_NHEF");
    }
    ~HLTJetIDHelper() {}
    
    bool operator()(const reco::CaloJet& caloJet, const reco::JetID& caloJetID);
    bool operator()(const reco::PFJet& pfJet);
    
    // caloJetID
    double min_EMF_;         // minimum EMF
    double max_EMF_;         // maximum EMF
    int min_N90_;            // mininum N90
    int min_N90hits_;        // mininum Nhit90
    
    // pfJetID
    double min_CEEF_;
    double max_CEEF_;
    double min_NEEF_;
    double max_NEEF_;
    double min_CHEF_;
    double max_CHEF_;
    double min_NHEF_;
    double max_NHEF_;
};

// refer HLTrigger/JetMET/src/HLTCaloJetIDProducer.cc
bool HLTJetIDHelper::operator()(const reco::CaloJet& caloJet, const reco::JetID& caloJetID) {
    const reco::CaloJet* calojetc = &caloJet;  // rename
    if (std::abs(calojetc->eta()) >= 2.6) {
      return true;
    } else {
        //if (min_N90hits_>0) jetID_.calculate(iEvent, *calojetc);
        //if ((calojetc->emEnergyFraction() >= min_EMF_) && ((min_N90hits_<=0) || (jetID_.n90Hits() >= min_N90hits_))  && (calojetc->n90() >= min_N90_) && (calojetc->emEnergyFraction() <= max_EMF_)) {
        if ((calojetc->emEnergyFraction() >= min_EMF_) && ((min_N90hits_<=0) || (caloJetID.n90Hits >= min_N90hits_))  && (calojetc->n90() >= min_N90_) && (calojetc->emEnergyFraction() <= max_EMF_)) {
            return true;
        }
    }
    return false;
}

// refer HLTrigger/JetMET/src/HLTPFEnergyFractionsFilter.cc
bool HLTJetIDHelper::operator()(const reco::PFJet& pfJet) {
    bool accept = true;
    const reco::PFJet* i = &pfJet;
    
    if(i->chargedEmEnergyFraction()<min_CEEF_) accept = false;
    if(i->chargedEmEnergyFraction()>max_CEEF_) accept = false;
    //
    
    if(i->neutralEmEnergyFraction()<min_NEEF_) accept = false;
    if(i->neutralEmEnergyFraction()>max_NEEF_) accept = false;
    //
    
    if(i->chargedHadronEnergyFraction()<min_CHEF_) accept = false;
    if(i->chargedHadronEnergyFraction()>max_CHEF_) accept = false;
    //
    
    if(i->neutralHadronEnergyFraction()<min_NHEF_) accept = false;
    if(i->neutralHadronEnergyFraction()>max_NHEF_) accept = false;
    //
    
    return accept;
}
    
#endif  // JETIDHELPER_H_

