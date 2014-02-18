#!/usr/bin/env python

from ROOT import TFile, TH1F, gROOT, gSystem, gInterpreter

# from compactify_cfg.py
triggers = [
        # MET
        'HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5',
        'HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5',
        'HLT_DiCentralPFJet30_PFMET80_v6',
        'HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4',
        'HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9',
        'HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9',
        'HLT_L1ETM100_v2',
        'HLT_L1ETM30_v2',
        'HLT_L1ETM40_v2',
        'HLT_L1ETM70_v2',
        #'HLT_MET120_HBHENoiseCleaned_v6',
        #'HLT_MET120_v13',
        #'HLT_MET200_HBHENoiseCleaned_v5',
        #'HLT_MET200_v12',
        #'HLT_MET300_HBHENoiseCleaned_v5',
        #'HLT_MET300_v4',
        #'HLT_MET400_HBHENoiseCleaned_v5',
        #'HLT_MET400_v7',
        'HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4',
        'HLT_PFMET150_v7',
        'HLT_PFMET180_v7',
        # HTMHT
        'HLT_PFNoPUHT350_PFMET100_v4',
        'HLT_PFNoPUHT400_PFMET100_v4',
        # JetHT
        'HLT_MET80_Track50_dEdx3p6_v6',
        'HLT_MET80_Track60_dEdx3p7_v6',
        'HLT_MET80_v5',
        ]

reftrig = "HLT_L1ETM40_v2"
ireftrig = triggers.index(reftrig)

metfilters = [
    "p_HBHENoiseFilter",
    "p_CSCTightHaloFilter",
    "p_hcalLaserEventFilter",
    "p_EcalDeadCellTriggerPrimitiveFilter",
    "p_trackingFailureFilter",
    "p_eeBadScFilter",
    "p_ecalLaserCorrFilter",
    "p_trkPOGFilters",
    ]

optmetfilters = [
    "p_goodVerticesFilter",
    "p_noscraping",
    "p_hcallasereventfilter2012",
    "p_EcalDeadCellBoundaryEnergyFilter",
    "p_tobtecfakesFilters",
    ]


#_______________________________________________________________________________

#gSystem.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gInterpreter.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gROOT.LoadMacro("../src/SimpleCandidateLinkDef.h")
gROOT.LoadMacro("HelperFunctions.h")

#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
#tfile.tree = tfile.Events
tfile = TFile.Open("../bin/compactified.L1ETM40.3.root")
tree = tfile.tree

sections = {}
sections["puremet"]         = False
sections["puremet_clean"]   = False
sections["monojet"]         = False
sections["monojet_clean"]   = False
sections["higdijet"]        = False
sections["higdijet_clean"]  = False
sections["susdijet"]        = False
sections["susdijet_clean"]  = False
sections["multijet"]        = False
sections["multijet_clean"]  = False
sections["bjet"]            = False
sections["bjet_clean"]      = False
sections["vbf"]             = False
sections["vbf_clean"]       = False

sections["puremet_venn"]    = False

sections["future_triggers"] = True


# ______________________________________________________________________________
# Functions

def count(sel_trigs):
    #histos = []
    #counts = []
    #for i in xrange(len(sel_trigs)):
    #    histos.append(TH1F("h%i" % i, "h%i" %i, 5, 0, 5))
    #    tree.Project(histos[i].GetName(), "event.json", ("(triggerFlags[%i]) * " % ireftrig) + sel_trigs[i], "goff")
    #    count1 = float(histos[i].Integral())
    #    count0 = float(histos[0].Integral())
    #    delta  = (count1 - count0) / count0 * 100
    #    counts.append((count0, count1, delta))
    #return counts

    counts_ = []
    counts = []
    for i in xrange(len(sel_trigs)):
        count = tree.GetEntries(("(triggerFlags[%i]) * " % ireftrig) + sel_trigs[i])
        counts_.append(float(count))
        count1 = counts_[i]
        count0 = counts_[0]
        delta  = (count1 - count0) / count0 * 100
        counts.append((count0, count1, delta, 100.+delta))
    return counts

def pprint(counts, overlap=False):
    for c in counts:
        (count0, count1, delta, delta1) = c
        print "%6.0f, %6.0f, %6.1f, %6.1f" % (count0, count1, delta, delta1)
    if overlap:
        assert(len(counts)==4)
        print "  overlap:      A=%4.0f%%       B=%4.0f%%  A&&B=%4.0f%%" % (counts[1][3], counts[2][3], counts[3][3])
        print "  overlap: only A=%4.0f%%  only B=%4.0f%%  A&&B=%4.0f%%" % (counts[1][3]-counts[3][3], counts[2][3]-counts[3][3], counts[3][3])
    print

if sections["puremet"]:
    sel_trigs = [
        "(hltCaloMET.pt>80 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)",
        ]
    counts = count(sel_trigs); pprint(counts)

if sections["puremet_venn"]:
    sel_trigs = [
        "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<50 || hltCaloMETCleanUsingJetID.pt<50))",
        "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<50 && hltCaloMETCleanUsingJetID.pt<9999))",
        "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<9999 && hltCaloMETCleanUsingJetID.pt<50))",
        "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<50 && hltCaloMETCleanUsingJetID.pt<50))",
        ]
    counts = count(sel_trigs); pprint(counts, overlap=True)

if sections["puremet_clean"]:
    sel_trigs = [
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && abs(deltaPhi(hltTrackMET.phi,hltPFMET.phi))<0.7 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltTrackMET.pt>10 && hltPFMET.pt>150)",
        #"(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltTrackMET.pt>15 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.sumEt>0 && hltTrackMET.sumEt/hltPFMET.sumEt>0.02 && hltPFMET.pt>150)",
        ]
    counts = count(sel_trigs); pprint(counts)

if sections["monojet_clean"]:
    sel_trigs = [
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>0)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100 && hltPFJetsL1FastL2L3[0].nhf<0.95)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>0)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100 && hltPFJetsL1FastL2L3[0].nhf<0.90)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>0)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100 && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>0)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100 && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].chf>0.01)",
        ]
    counts = count(sel_trigs); pprint(counts)



if sections["future_triggers"]:
    ## PFMET150
    #sel_trig0 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    #sel_trig1 = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)"

    ## MonoCentralJet
    #sel_trig0 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    #sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)"
    ##sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>80 && hltCaloMETClean.pt>70 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>120||hltPFMET.pt>120) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)"

    # DiCentralJetHIG
    sel_trig0 = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    sel_trig0 = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFGlobal.dijet_mindphi_2cj>0.5)"

    # DiCentralJetSUS
    #sel_trig0 = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    sel_trig1 = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) && hltPFGlobal.dijet_mindphi_2cj>0.3)"


    # Count!
    sel_trigs = [
        sel_trig0 + "||" + sel_trig1,
        sel_trig0,
        sel_trig1,
        sel_trig0 + "&&" + sel_trig1,
        ]
    counts = count(sel_trigs); pprint(counts, overlap=True)





