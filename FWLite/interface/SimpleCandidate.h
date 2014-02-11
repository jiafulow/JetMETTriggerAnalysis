#ifndef SIMPLECANDIDATE_H_
#define SIMPLECANDIDATE_H_

// NOTE: If you modify/add any class, do `scram b -r` to remake the dictionary

//#include "Math/LorentzVector.h"
//namespace math {
//    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  XYZTLorentzVector;
//}

namespace simple {
    class XYVector {
      public:
        float px, py;
        float pt, phi;  // trade memory for simplicity
    };

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
        float fHPD;
        int n90Hits;
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
        int ntot;
        float csv;
    };

    // refer DataFormats/METReco/interface/MET.h
    class MET : public XYVector {
      public:
        float sumEt;
    };


    class Global {
      public:
        float monojet_maxpt_idx;
        float dijet_maxpt;
        float dijet_maxpt_mjj;
        float dijet_mindphi_j30;
        float dijet_mindphi_j40;
        float dijet_mindphi_cj30;
        float dijet_mindphi_cj40;
        float dijet_mindphi_2cj;
        float dijet_mindphi_3cj;
        float bjet_maxcsv;
        float bjet_maxcsv2;
        float vbf_maxmjj;
        float vbf_maxmjj_deta;
        float vbf_leadmjj;
        float vbf_leadmjj_deta;
        float rho_kt6;
        unsigned int njets_j30;
        unsigned int njets_cj30;
    };


    class Event {
      public:
        unsigned int run;
        unsigned int lumi;
        unsigned int event;
        unsigned int nPV;
        unsigned int nGoodPV;
        unsigned int nTruePV;
        bool json;
    };

}

#endif  // SIMPLECANDIDATE_H_

