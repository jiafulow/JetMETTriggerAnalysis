#ifndef JetMETTriggerAnalysis_Utilities_PFPileUpFlag_h_
#define JetMETTriggerAnalysis_Utilities_PFPileUpFlag_h_

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"


class PFPileUpFlag : public edm::EDProducer {
  public:
    typedef std::vector<edm::FwdPtr<reco::PFCandidate> >   PFCollection;
    typedef edm::View<reco::PFCandidate>                   PFView;
    typedef std::vector<reco::PFCandidate>                 PFCollectionByValue;
    
    explicit PFPileUpFlag(const edm::ParameterSet&);
    ~PFPileUpFlag();

  private:
    //virtual void beginJob() ;
    virtual void produce(edm::Event&, const edm::EventSetup&);
    //virtual void endJob() ;

    edm::InputTag pfCollection_;
    edm::InputTag pfPileUpCollection_;
};

#endif

