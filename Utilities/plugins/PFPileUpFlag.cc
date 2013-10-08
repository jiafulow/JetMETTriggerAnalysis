#include "JetMETTriggerAnalysis/Utilities/plugins/PFPileUpFlag.h"


PFPileUpFlag::PFPileUpFlag(const edm::ParameterSet& iConfig)
      : pfCollection_(iConfig.getParameter<edm::InputTag>("pfCollection")),
        pfPileUpCollection_(iConfig.getParameter<edm::InputTag>("pfPileUpCollection")) {
    produces<std::vector<bool> >();
}


PFPileUpFlag::~PFPileUpFlag() {}


void PFPileUpFlag::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    Handle<PFView> pfCandidates;
    iEvent.getByLabel(pfCollection_, pfCandidates);
    //assert(pfCandidates->size() > 0);
    
    Handle<PFView> pfPileUpCandidates;
    iEvent.getByLabel(pfPileUpCollection_, pfPileUpCandidates);
    
    std::auto_ptr<std::vector<bool> > pOutput(new std::vector<bool>);
    
    unsigned int count = 0;
    PFView::const_iterator it1 = pfCandidates->begin();
    PFView::const_iterator end1 = pfCandidates->end();
    for (; it1 != end1; ++it1) {
        bool found = false;
        
        const reco::Candidate::LorentzVector& p4_1 = (*it1).p4();
        
        PFView::const_iterator it2 = pfPileUpCandidates->begin();
        PFView::const_iterator end2 = pfPileUpCandidates->end();
        for (; it2 != end2; ++it2) {
            const reco::Candidate::LorentzVector& p4_2 = (*it2).p4();
            if (fabs(p4_1.px() - p4_2.px()) < 1e-6 &&
                fabs(p4_1.py() - p4_2.py()) < 1e-6 &&
                fabs(p4_1.pz() - p4_2.pz()) < 1e-6 &&
                fabs(p4_1.e()  - p4_2.e() ) < 1e-6) {
                found = true;
                ++count;
                break;
            }
        }
        pOutput->push_back(found);
    }
    // pfPileUpCandidates must be strictly a subset of pfCandidates
    assert(count == pfPileUpCandidates->size());
    // Sanity check
    assert(pOutput->size() == pfCandidates->size());
    
    iEvent.put(pOutput);
}


// Define this as a plug-in
DEFINE_FWK_MODULE(PFPileUpFlag);

