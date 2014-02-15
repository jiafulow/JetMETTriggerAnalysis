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
sections["puremet"]         = True
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

sections["puremet_venn"]    = True


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
        print "  overlap:      A=%4.0f%%     B=%4.0f%%    A&&B=%4.0f%%" % (counts[1][3], counts[2][3], counts[3][3])
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

if sections["puremet_clean"]:
    sel_trigs = [
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && abs(deltaPhi(hltTrackMET.phi,hltPFMET.phi))<0.7 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && hltTrackMET.pt>15 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && hltTrackMET.sumEt/hltPFMET.pt>0.1 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && hltTrackMET.sumEt/hltPFMET.pt>0.02 && hltPFMET.pt>150)",
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

