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

typedef reco::Candidate::LorentzVector FourVector;
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
    //simple.nhf    = cand.neutralHadronEnergyFraction();  // used by HLTPFEnergyFractionsFilter.h
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
unsigned int eval_maxpt_idx(const std::vector<T>& jets, double ptmin, double etamax) {  // no jetID
    unsigned int idx = 0;
    double maxpt = -99.;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        const T& jet = jets.at(j);
        if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {
            if (maxpt < jet.pt()) {
                idx = j;
                maxpt = jet.pt();
            }
        }
    }
    return idx;
}

template<typename T>
unsigned int eval_njets(const std::vector<T>& jets, double ptmin, double etamax) {  // no jetID
    unsigned int count = 0;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        const T& jet = jets.at(j);
        if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {
            ++count;
        }
    }
    return count;
}

template<typename T>
double eval_ht(const std::vector<T>& jets, double ptmin, double etamax) {  // no jetID
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
FourVector eval_mht(const std::vector<T>& jets, double ptmin, double etamax) {  // no jetID
    FourVector mht(0., 0., 0., 0.);
    for (unsigned int j = 0; j < jets.size(); ++j) {
        const T& jet = jets.at(j);
        if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {
            mht -= jets.at(j).p4();
        }
    }
    return mht;
}

template<typename T>
double eval_mindphi(double metphi, const std::vector<T>& jets, double ptmin, double etamax, int njmax=999) {  // no jetID
    double mindphi = M_PI;
    int nj = 0;
    for (unsigned int j = 0; j < jets.size() && nj < njmax; ++j) {
        const T& jet = jets.at(j);
        if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {
            double dphi = fabs(reco::deltaPhi(jet.phi(), metphi));
            if (mindphi > dphi)
                mindphi = dphi;
            ++nj;
        }
    }
    return mindphi;
}

template<typename T>
std::vector<double> eval_maxcsv(const std::vector<T>& jets, double ptmin, double etamax) {
    std::vector<double> results(2, -99.);

    std::vector<double> values;
    for (unsigned int j = 0; j < jets.size(); ++j) {
        const T& jet = jets.at(j);
        if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {
            const pat::Jet* patJet = dynamic_cast<const pat::Jet*>(&jet);
            double csv = (patJet != 0 ? patJet->bDiscriminator("combinedSecondaryVertexBJetTags") : -999);
            values.push_back(csv);
        }
    }
    std::sort(values.begin(), values.end(), std::greater<double>());

    if (values.size() > 0)
        results[0] = values[0];
    if (values.size() > 1)
        results[1] = values[1];
    return results;
}

template<typename T>
std::vector<double> eval_maxcsv(const edm::AssociationVector<edm::RefToBaseProd<T>, std::vector<float> >& tags, double ptmin, double etamax) {
    std::vector<double> results(2, -99.);

    std::vector<double> values;
    //typename edm::AssociationVector<edm::RefToBaseProd<T>, std::vector<float> >::const_iterator tag;
    //for (tag = tags.begin(); tag != tags.end(); ++tag) {
    //    if (tag->first->pt() > ptmin && fabs(tag->first->eta()) < etamax) {
    //        double csv = tag->second;
    //        values.push_back(csv);
    //    }
    //}
    for (unsigned int i = 0; i < tags.size(); ++i) {  //FIXME: forgot to store the product
        double csv = tags.value(i);
        values.push_back(csv);
    }
    std::sort(values.begin(), values.end(), std::greater<double>());

    if (values.size() > 0)
        results[0] = values[0];
    if (values.size() > 1)
        results[1] = values[1];
    return results;
}

template<typename T>
std::vector<double> eval_maxmjj(const std::vector<T>& jets, double ptmin, double etamax, double minmjj, double mindeta, bool checkEtaOpposite, bool leadingJetOnly) {  // no jetID
    std::vector<double> results(2, -99.);  // mjj, deta
    double maxsumpt = -99;  // scalar sum

    if (jets.size() >= 2) {
        for (unsigned int j = 0; j < jets.size()-1; ++j) {
            if (leadingJetOnly && j >= 2)  break;
            const T& jet = jets.at(j);
            if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {

                for (unsigned int k = j+1; k < jets.size(); ++k) {
                    if (leadingJetOnly && k >= 3)  break;
                    const T& ket = jets.at(k);
                    if (ket.pt() > ptmin && fabs(ket.eta()) < etamax) {
                        double deta = fabs(jet.eta() - ket.eta());
                        double mjj = (jet.p4() + ket.p4()).mass();
                        double sumpt = jet.pt() + ket.pt();
                        bool etaOpposite = (jet.eta()*ket.eta() < 0);

                        if ( (maxsumpt < sumpt) &&
                             ((checkEtaOpposite && etaOpposite) || !checkEtaOpposite) &&
                             (deta > mindeta) &&
                             (mjj > minmjj) ) {

                            maxsumpt = sumpt;
                            results[0] = mjj;
                            results[1] = deta;
                        }
                    }
                }  // end for loop
            }
        }  // end for loop
    }
    return results;
}

template<typename T>
std::vector<double> eval_maxptjj(const std::vector<T>& jets, double ptmin, double etamax) {  // no jetID
    std::vector<double> results(2, -99.);  // ptjj, mjj
    double maxptjj = -99;  // vectorial sum

    if (jets.size() >= 2) {
        for (unsigned int j = 0; j < jets.size()-1; ++j) {
            const T& jet = jets.at(j);
            if (jet.pt() > ptmin && fabs(jet.eta()) < etamax) {

                for (unsigned int k = j+1; k < jets.size(); ++k) {
                    const T& ket = jets.at(k);
                    if (ket.pt() > ptmin && fabs(ket.eta()) < etamax) {
                        double mjj = (jet.p4() + ket.p4()).mass();
                        double ptjj = (jet.p4() + ket.p4()).pt();

                        if (maxptjj < ptjj) {
                            maxptjj = ptjj;
                            results[0] = ptjj;
                            results[1] = mjj;
                        }
                    }
                }  // end for loop

            }
        }  // end for loop
    }
    return results;
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
    const bool isGolden = analyzerpset.getParameter<bool>("isGolden");
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
    simple::MET l1MET;
    simple::MET l1MHT;
    simple::MET hltCaloMET;
    simple::MET hltCaloMETClean;
    simple::MET hltCaloMETCleanUsingJetID;
    simple::MET hltPFMET;
    simple::MET hltPFMETNoMu;
    simple::MET hltPFMETCleanUsingJetID;
    simple::MET hltTrackMET;
    simple::MET hltCaloHTMHT;
    simple::MET hltPFHTMHT;
    simple::MET hltPFHTMHTNoPU;
    simple::Global hltCaloGlobal;
    simple::Global hltPFGlobal;
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
    simple::MET patHTMHT;
    simple::Global patGlobal;
    unsigned int patJetTopo;


    TFile* outfile = new TFile(outfilename, "RECREATE");
    TTree* outtree = new TTree("tree", "tree");
    outtree->Branch("event", &simpleEvent);
    outtree->Branch("weightGenEvent", &weightGenEvent);
    outtree->Branch("weightPU", &weightPU);  //FIXME not yet set
    outtree->Branch("weightPdf", &weightPdf);  //FIXME not yet set
    outtree->Branch("triggerFlags", triggerFlags, Form("triggerFlags[%i]/b", nTriggerFlags));
    outtree->Branch("metfilterFlags", metfilterFlags, Form("metfilterFlags[%i]/b", nMetfilterFlags));
    outtree->Branch("optmetfilterFlags", optmetfilterFlags, Form("optmetfilterFlags[%i]/b", nOptmetfilterFlags));
    // L1
    outtree->Branch("l1MET", &l1MET);
    outtree->Branch("l1MHT", &l1MHT);
    // HLT
    outtree->Branch("hltPFCandidates", &hltPFCandidates);
    outtree->Branch("hltCaloJets", &hltCaloJets);
    outtree->Branch("hltCaloJetsL1Fast", &hltCaloJetsL1Fast);
    outtree->Branch("hltPFJets", &hltPFJets);
    outtree->Branch("hltPFJetsNoPU", &hltPFJetsNoPU);
    outtree->Branch("hltPFJetsL1FastL2L3", &hltPFJetsL1FastL2L3);
    outtree->Branch("hltPFJetsL1FastL2L3NoPU", &hltPFJetsL1FastL2L3NoPU);
    outtree->Branch("hltCaloMET", &hltCaloMET);
    outtree->Branch("hltCaloMETClean", &hltCaloMETClean);
    outtree->Branch("hltCaloMETCleanUsingJetID", &hltCaloMETCleanUsingJetID);
    outtree->Branch("hltPFMET", &hltPFMET);
    outtree->Branch("hltPFMETNoMu", &hltPFMETNoMu);
    outtree->Branch("hltPFMETCleanUsingJetID", &hltPFMETCleanUsingJetID);
    outtree->Branch("hltTrackMET", &hltTrackMET);
    outtree->Branch("hltCaloHTMHT", &hltCaloHTMHT);
    outtree->Branch("hltPFHTMHT", &hltPFHTMHT);
    outtree->Branch("hltPFHTMHTNoPU", &hltPFHTMHTNoPU);
    outtree->Branch("hltCaloGlobal", &hltCaloGlobal);
    outtree->Branch("hltPFGlobal", &hltPFGlobal);
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
    outtree->Branch("patHTMHT", &patHTMHT);
    outtree->Branch("patGlobal", &patGlobal);
    outtree->Branch("patJetTopo", &patJetTopo);


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
        if (isGolden && !goodjson)  continue;
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
            hltPFJetsL1FastL2L3NoPU.push_back(simple);
        }

        //______________________________________________________________________
        // l1MET
        if (verbose)  std::cout << "compactify: Begin filling l1MET..." << std::endl;
        l1MET.pt     = handler.l1METs->begin()->etMiss();
        l1MET.phi    = handler.l1METs->begin()->phi();
        l1MET.sumEt  = handler.l1METs->begin()->etTotal();
        l1MET.px     = l1MET.pt * cos(l1MET.phi);
        l1MET.py     = l1MET.pt * sin(l1MET.phi);

        //______________________________________________________________________
        // l1MHT
        if (verbose)  std::cout << "compactify: Begin filling l1MHT..." << std::endl;
        l1MHT.pt     = handler.l1MHTs->begin()->etMiss();
        l1MHT.phi    = handler.l1MHTs->begin()->phi();
        l1MHT.sumEt  = handler.l1MHTs->begin()->etTotal();
        l1MHT.px     = l1MHT.pt * cos(l1MHT.phi);
        l1MHT.py     = l1MHT.pt * sin(l1MHT.phi);

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
        //simple_fill(handler.hltCaloMETCleansUsingJetID->at(0), hltCaloMETCleanUsingJetID);  // FIXME

        FourVector hltCaloMETCleanUsingJetID_p4(hltCaloMET.px, hltCaloMET.py, 0, 0);
        double hltCaloMETCleanUsingJetID_sumEt = 0.;
        for (unsigned int i = 0; i < handler.hltCaloJets->size(); ++i) {
            const reco::CaloJet& jet = handler.hltCaloJets->at(i);
            if (jet.et() < 20 || fabs(jet.eta()) > 5.0)  continue;

            edm::ProductID productid = handler.hltCaloJetIDs->ids().front().first;
            const reco::JetID clsJetID = handler.hltCaloJetIDs->get(productid, i);
            bool jetID = hltJetIDHelper(jet, clsJetID);
            if (jetID)  continue;  // anti-jetID

            double pt  = jet.et();
            double phi = jet.phi();
            double px  = jet.et() * cos(phi);
            double py  = jet.et() * sin(phi);
            hltCaloMETCleanUsingJetID_p4 += FourVector(px, py, 0, 0);
            hltCaloMETCleanUsingJetID_sumEt -= pt;
        }

        hltCaloMETCleanUsingJetID.px     = hltCaloMETCleanUsingJetID_p4.px();
        hltCaloMETCleanUsingJetID.py     = hltCaloMETCleanUsingJetID_p4.py();
        hltCaloMETCleanUsingJetID.pt     = hltCaloMETCleanUsingJetID_p4.pt();
        hltCaloMETCleanUsingJetID.phi    = hltCaloMETCleanUsingJetID_p4.phi();
        hltCaloMETCleanUsingJetID.sumEt  = hltCaloMETCleanUsingJetID_sumEt;


        //______________________________________________________________________
        // hltPFMET
        if (verbose)  std::cout << "compactify: Begin filling hltPFMET..." << std::endl;
        //simple_fill(handler.hltPFMETs->at(0), hltPFMET);  // FIXME

        FourVector hltPFMET_p4(0,0,0,0);
        double hltPFMET_sumEt = 0.;
        for (unsigned int i = 0; i < handler.hltPFCandidates->size(); ++i) {
            const reco::PFCandidate& cand = handler.hltPFCandidates->at(i);
            hltPFMET_p4 -= cand.p4();
            hltPFMET_sumEt += cand.pt();
        }
        hltPFMET.px     = hltPFMET_p4.px();
        hltPFMET.py     = hltPFMET_p4.py();
        hltPFMET.pt     = hltPFMET_p4.pt();
        hltPFMET.phi    = hltPFMET_p4.phi();
        hltPFMET.sumEt  = hltPFMET_sumEt;

        //______________________________________________________________________
        // hltPFMETNoMu
        if (verbose)  std::cout << "compactify: Begin filling hltPFMETNoMu..." << std::endl;
        simple_fill(handler.hltPFMETsNoMu->at(0), hltPFMETNoMu);

        //______________________________________________________________________
        // hltPFMETCleanUsingJetID
        if (verbose)  std::cout << "compactify: Begin filling hltPFMETCleanUsingJetID..." << std::endl;
        simple_fill(handler.hltPFMETCleansUsingJetID->at(0), hltPFMETCleanUsingJetID);

        //______________________________________________________________________
        // hltTrackMET
        if (verbose)  std::cout << "compactify: Begin filling hltTrackMET..." << std::endl;
        //simple_fill(handler.hltTrackMETs->at(0), hltTrackMET);  // FIXME

        FourVector hltTrackMET_p4(0,0,0,0);
        double hltTrackMET_sumEt = 0.;
        for (unsigned int i = 0; i < handler.hltPFCandidates->size(); ++i) {
            const reco::PFCandidate& cand = handler.hltPFCandidates->at(i);
            const bool& isPU = handler.hltPFPileUpFlags->at(i);
            if (!isPU && cand.charge() != 0) {
                hltTrackMET_p4 -= cand.p4();
                hltTrackMET_sumEt += cand.pt();
            }
        }
        hltTrackMET.px     = hltTrackMET_p4.px();
        hltTrackMET.py     = hltTrackMET_p4.py();
        hltTrackMET.pt     = hltTrackMET_p4.pt();
        hltTrackMET.phi    = hltTrackMET_p4.phi();
        hltTrackMET.sumEt  = hltTrackMET_sumEt;

        //______________________________________________________________________
        // hltCaloHTMHT
        if (verbose)  std::cout << "compactify: Begin filling hltCaloHTMHT..." << std::endl;
        simple_fill(handler.hltCaloHTMHTs->at(0), hltCaloHTMHT);

        //______________________________________________________________________
        // hltPFHTMHT
        //if (verbose)  std::cout << "compactify: Begin filling hltPFHTMHT..." << std::endl;
        //simple_fill(handler.hltPFHTMHTs->at(0), hltPFHTMHT);

        //______________________________________________________________________
        // hltPFHTMHTNoPU
        if (verbose)  std::cout << "compactify: Begin filling hltPFHTMHTNoPU..." << std::endl;
        simple_fill(handler.hltPFHTMHTsNoPU->at(0), hltPFHTMHTNoPU);

        //______________________________________________________________________
        // hltCaloGlobal
        if (verbose)  std::cout << "compactify: Begin filling hltCaloGlobal..." << std::endl;
        std::vector<double> hltCaloGlobal_results;
        hltCaloGlobal.monojet_maxpt_idx  = eval_maxpt_idx(*handler.hltCaloJetsL1Fast, 80., 2.6);
        hltCaloGlobal_results            = eval_maxptjj(*handler.hltCaloJetsL1Fast, 15., 2.6);
        hltCaloGlobal.dijet_maxpt        = hltCaloGlobal_results[0];
        hltCaloGlobal.dijet_maxpt_mjj    = hltCaloGlobal_results[1];
        hltCaloGlobal.dijet_mindphi_j30  = eval_mindphi(hltCaloMET.phi, *handler.hltCaloJetsL1Fast, 30., 5.0);
        hltCaloGlobal.dijet_mindphi_j40  = eval_mindphi(hltCaloMET.phi, *handler.hltCaloJetsL1Fast, 40., 5.0);
        hltCaloGlobal.dijet_mindphi_cj30 = eval_mindphi(hltCaloMET.phi, *handler.hltCaloJetsL1Fast, 30., 2.6);
        hltCaloGlobal.dijet_mindphi_cj40 = eval_mindphi(hltCaloMET.phi, *handler.hltCaloJetsL1Fast, 40., 2.6);
        hltCaloGlobal.dijet_mindphi_2cj  = eval_mindphi(hltCaloMET.phi, *handler.hltCaloJetsL1Fast, 30., 2.6, 2);
        hltCaloGlobal.dijet_mindphi_3cj  = eval_mindphi(hltCaloMET.phi, *handler.hltCaloJetsL1Fast, 30., 2.6, 3);
        hltCaloGlobal_results            = eval_maxcsv(*handler.hltCSVBJetTags, 20., 2.6);
        hltCaloGlobal.bjet_maxcsv        = hltCaloGlobal_results[0];
        hltCaloGlobal.bjet_maxcsv2       = hltCaloGlobal_results[1];
        hltCaloGlobal_results            = eval_maxmjj(*handler.hltCaloJetsL1Fast, 30., 5.0, 500., 3.5, true, false);
        hltCaloGlobal.vbf_maxmjj         = hltCaloGlobal_results[0];
        hltCaloGlobal.vbf_maxmjj_deta    = hltCaloGlobal_results[1];
        hltCaloGlobal_results            = eval_maxmjj(*handler.hltCaloJetsL1Fast, 30., 5.0, 500., 3.5, true, true);
        hltCaloGlobal.vbf_leadmjj        = hltCaloGlobal_results[0];
        hltCaloGlobal.vbf_leadmjj_deta   = hltCaloGlobal_results[1];
        hltCaloGlobal.rho_kt6            = *handler.hltRho_kt6CaloJets;
        hltCaloGlobal.njets_j30          = eval_njets(*handler.hltCaloJetsL1Fast, 30., 5.0);
        hltCaloGlobal.njets_cj30         = eval_njets(*handler.hltCaloJetsL1Fast, 30., 2.5);

        //______________________________________________________________________
        // hltPFGlobal
        if (verbose)  std::cout << "compactify: Begin filling hltPFGlobal..." << std::endl;
        std::vector<double> hltPFGlobal_results;
        hltPFGlobal.monojet_maxpt_idx  = eval_maxpt_idx(*handler.hltPFJetsL1FastL2L3, 80., 2.6);
        hltPFGlobal_results            = eval_maxptjj(*handler.hltPFJetsL1FastL2L3, 25., 2.6);
        hltPFGlobal.dijet_maxpt        = hltPFGlobal_results[0];
        hltPFGlobal.dijet_maxpt_mjj    = hltPFGlobal_results[1];
        hltPFGlobal.dijet_mindphi_j30  = eval_mindphi(hltPFMET.phi, *handler.hltPFJetsL1FastL2L3, 30., 5.0);
        hltPFGlobal.dijet_mindphi_j40  = eval_mindphi(hltPFMET.phi, *handler.hltPFJetsL1FastL2L3, 40., 5.0);
        hltPFGlobal.dijet_mindphi_cj30 = eval_mindphi(hltPFMET.phi, *handler.hltPFJetsL1FastL2L3, 30., 2.6);
        hltPFGlobal.dijet_mindphi_cj40 = eval_mindphi(hltPFMET.phi, *handler.hltPFJetsL1FastL2L3, 40., 2.6);
        hltPFGlobal.dijet_mindphi_2cj  = eval_mindphi(hltPFMET.phi, *handler.hltPFJetsL1FastL2L3, 30., 2.6, 2);
        hltPFGlobal.dijet_mindphi_3cj  = eval_mindphi(hltPFMET.phi, *handler.hltPFJetsL1FastL2L3, 30., 2.6, 3);
        hltPFGlobal.bjet_maxcsv        = -99.;
        hltPFGlobal.bjet_maxcsv2       = -99.;
        hltPFGlobal_results            = eval_maxmjj(*handler.hltPFJetsL1FastL2L3, 40., 5.0, 500., 3.5, true, false);  // should be 600
        hltPFGlobal.vbf_maxmjj         = hltPFGlobal_results[0];
        hltPFGlobal.vbf_maxmjj_deta    = hltPFGlobal_results[1];
        hltPFGlobal_results            = eval_maxmjj(*handler.hltPFJetsL1FastL2L3, 40., 5.0, 500., 3.5, true, true);  // should be 600
        hltPFGlobal.vbf_leadmjj        = hltPFGlobal_results[0];
        hltPFGlobal.vbf_leadmjj_deta   = hltPFGlobal_results[1];
        hltPFGlobal.rho_kt6            = *handler.hltRho_kt6PFJets;
        hltPFGlobal.njets_j30          = eval_njets(*handler.hltPFJetsL1FastL2L3, 30., 5.0);
        hltPFGlobal.njets_cj30         = eval_njets(*handler.hltPFJetsL1FastL2L3, 30., 2.5);

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

        //______________________________________________________________________
        // patMPT
        FourVector patMPT_p4(0,0,0,0);
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
        // patHTMHT
        FourVector patHTMHT_p4(0,0,0,0);
        double patHTMHT_sumEt = 0.;
        for (unsigned int i = 0; i < handler.patJets->size(); ++i) {
            const pat::Jet& jet = handler.patJets->at(i);
            bool jetID = jetIDHelper(jet);
            if (jet.pt() > 30 && fabs(jet.eta()) < 4.5 && jetID) {  // RA2, RA2b definition
                patHTMHT_p4 -= jet.p4();
            }
            if (jet.pt() > 50 && fabs(jet.eta()) < 2.5 && jetID) {  // RA2, RA2b definition
                patHTMHT_sumEt += jet.pt();
            }
        }
        patHTMHT.px    = patHTMHT_p4.px();
        patHTMHT.py    = patHTMHT_p4.py();
        patHTMHT.pt    = patHTMHT_p4.pt();
        patHTMHT.phi   = patHTMHT_p4.phi();
        patHTMHT.sumEt = patHTMHT_sumEt;

        //______________________________________________________________________
        // patGlobal
        if (verbose)  std::cout << "compactify: Begin filling patGlobal..." << std::endl;
        std::vector<double> patGlobal_results;
        patGlobal.monojet_maxpt_idx  = eval_maxpt_idx(*handler.patJets, 30., 2.5);
        patGlobal_results            = eval_maxptjj(*handler.patJets, 30., 2.5);
        patGlobal.dijet_maxpt        = patGlobal_results[0];
        patGlobal.dijet_maxpt_mjj    = patGlobal_results[1];
        patGlobal.dijet_mindphi_j30  = eval_mindphi(recoPFMETT0T1.phi, *handler.patJets, 30., 4.7);
        patGlobal.dijet_mindphi_j40  = eval_mindphi(recoPFMETT0T1.phi, *handler.patJets, 40., 4.7);
        patGlobal.dijet_mindphi_cj30 = eval_mindphi(recoPFMETT0T1.phi, *handler.patJets, 30., 2.5);
        patGlobal.dijet_mindphi_cj40 = eval_mindphi(recoPFMETT0T1.phi, *handler.patJets, 40., 2.5);
        patGlobal.dijet_mindphi_2cj  = eval_mindphi(recoPFMETT0T1.phi, *handler.patJets, 30., 2.5, 2);
        patGlobal.dijet_mindphi_3cj  = eval_mindphi(recoPFMETT0T1.phi, *handler.patJets, 30., 2.5, 3);
        patGlobal_results            = eval_maxcsv(*handler.patJets, 30., 2.5);
        patGlobal.bjet_maxcsv        = patGlobal_results[0];
        patGlobal.bjet_maxcsv2       = patGlobal_results[1];
        patGlobal_results            = eval_maxmjj(*handler.patJets, 30., 4.7, 500., 3.5, true, false);
        patGlobal.vbf_maxmjj         = patGlobal_results[0];
        patGlobal.vbf_maxmjj_deta    = patGlobal_results[1];
        patGlobal_results            = eval_maxmjj(*handler.patJets, 30., 4.7, 500., 3.5, true, true);
        patGlobal.vbf_leadmjj        = patGlobal_results[0];
        patGlobal.vbf_leadmjj_deta   = patGlobal_results[1];
        patGlobal.rho_kt6            = *handler.recoRho_kt6PFJets;
        patGlobal.njets_j30          = eval_njets(*handler.patJets, 30., 4.7);
        patGlobal.njets_cj30         = eval_njets(*handler.patJets, 30., 2.5);


        //______________________________________________________________________
        // patJetTopo
        if (verbose)  std::cout << "compactify: Begin filling patJetTopo..." << std::endl;
        patJetTopo = 0;  // FIXME


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

