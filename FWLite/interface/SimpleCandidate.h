#ifndef SIMPLECANDIDATE_H
#define SIMPLECANDIDATE_H

#include <cmath>

//#include "Math/LorentzVector.h"
//namespace math {
//    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >  XYZTLorentzVector;
//}

class SimpleLorentzVector {
  public:
    float px, py, pz, E;
    float pt, eta, phi;  // trade memory for simplicity
};

// refer DataFormats/Candidate/interface/Particle.h
class SimpleParticle : public SimpleLorentzVector {
  public:
    int charge;
    int pdgId;
    bool isPU;
};

// refer DataFormats/JetReco/interface/PFJet.h
class SimpleJet : public SimpleLorentzVector {
  public:
    //float chf;
    //float nhf;
    //float cef;
    //float nef;
    //int nch;
    //int nconstituents;
    bool jetIdL;
    bool jetIdT;
};

// refer DataFormats/METReco/interface/PFMET.h
class SimpleMET : public SimpleLorentzVector {
  public:
    float sumEt;
};



#endif
