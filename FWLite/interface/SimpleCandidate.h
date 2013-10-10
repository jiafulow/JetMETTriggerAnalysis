#ifndef SIMPLECANDIDATE_H_
#define SIMPLECANDIDATE_H_

// NOTE: If you modify/add any class, do `scram b -r` to remake the dictionary

//#include "Math/LorentzVector.h"
//namespace math {
//    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  XYZTLorentzVector;
//}

namespace simple {
    class LorentzVector {
      public:
        float px, py, pz, E;
        float pt, eta, phi;  // trade memory for simplicity
    };
    
    // refer DataFormats/Candidate/interface/Particle.h
    class Particle : public LorentzVector {
      public:
        int charge;
        int pdgId;
        bool isPU;
    };
    
    // refer DataFormats/JetReco/interface/Jet.h
    class Jet : public LorentzVector {
      public:
        float rawpt;
        bool jetID;
    };
    
    // refer DataFormats/JetReco/interface/CaloJet.h
    class CaloJet : public Jet {
      public:
        float emf;
        float fhf;
        int nHit;
        //float etaWidth;
        //float phiWidth;
    };
    
    // refer DataFormats/JetReco/interface/PFJet.h
    class PFJet : public Jet {
      public:
        float chf;
        float nhf;
        float cef;
        float nef;
        int nch;
        int nconstituents;  
    };
    
    // refer DataFormats/METReco/interface/MET.h
    class MET : public LorentzVector {
      public:
        float sumEt;
    };

    class Event {
      public:
        unsigned int run;
        unsigned int lumi;
        unsigned int event;
        unsigned int nPV;
        bool json;
    };

}

#endif  // SIMPLECANDIDATE_H_
