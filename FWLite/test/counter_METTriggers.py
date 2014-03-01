#!/usr/bin/env python

from ROOT import TFile, TH1, TH1F, TLatex, TLine, TChain, TCanvas, TPad, TLegend, gROOT, gInterpreter, gStyle, gSystem, gPad, kGray, kWhite, kBlack
from math import sqrt
from ordereddict import OrderedDict

# For init
class CounterInit:
    def __init__(self):
        # ROOT
        gROOT.LoadMacro("tdrstyle.C")
        gROOT.LoadMacro("HelperFunctions.h")
        gROOT.ProcessLine("setTDRStyle()")

        #gSystem.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
        gInterpreter.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
        gROOT.LoadMacro("../src/SimpleCandidateLinkDef.h")

        gStyle.SetEndErrorSize(2)
        gStyle.SetLabelSize(0.04, "Y")

        TH1.SetDefaultSumw2()


# ______________________________________________________________________________
# Text
latex = TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.026)

line = TLine()
line.SetLineColor(kGray+2)
line.SetLineStyle(2)


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
reftrig2 = "HLT_PFMET150_v7"
ireftrig2 = triggers.index(reftrig2)

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
# Configurations
counterinit = CounterInit()
chain = TChain("tree", "tree")
infiles = [
    "../bin/compactified.L1ETM40.4.root",
    ]
for f in infiles:
    chain.Add(f)
tree = chain

tfile2 = TFile.Open("../bin/compactified.PFMETplusX.0.root")
tree2 = tfile2.tree


sections = {}
sections["puremet"]         = False
sections["puremet_clean"]   = False
sections["monojet"]         = False
sections["monojet_clean"]   = False

sections["future_triggers"] = False
sections["future_numbers"]  = False
sections["future_numbers2"] = True

#imgdir = "figures_20131130/"  # for Torino workshop
imgdir = "figures_20140206/"  # for first draft
if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)


# ______________________________________________________________________________
# Functions

def count(sel_trigs):
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

htemp = TH1F("htemp", "", 4, 0, 4)
def count_future(odict, sel, normalizeby, benchmarks={}, exclusive=False):
    if benchmarks:  assert(len(benchmarks) == len(odict))

    counts = []
    exclsel = "(0)"  # USE OR
    exclcount = 0
    for k, v in odict.iteritems():
        vv = v

        if benchmarks:
            benchkey = benchmarks.keys()[odict.keys().index(k)]
            vv = "(%s * %s)" % (vv, benchmarks[benchkey])

        if exclusive:
            exclsel += ("||" + vv)
            vv = "(%s)" % exclsel

        selection = "*".join([sel, vv, normalizeby])
        tree.Project("htemp", "event.json", selection, "goff")
        count = float(htemp.Integral())
        if exclusive:
            count -= exclcount
            exclcount += count
        counts.append((k, count))
    return counts

def pprint_future(counts, refcount, refrate):
    for c in counts:
        (name, count) = c
        print "%-20s\t%.2f" % (name, count / refcount * refrate)
    print


# ______________________________________________________________________________
if sections["puremet"]:
    sel_trigs = [
        "(hltCaloMET.pt>80 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)",
        ]
    counts = count(sel_trigs); pprint(counts)

    #sel_trigs = [
    #    "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<50 || hltCaloMETCleanUsingJetID.pt<50))",
    #    "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<50 && hltCaloMETCleanUsingJetID.pt<9999))",
    #    "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<9999 && hltCaloMETCleanUsingJetID.pt<50))",
    #    "(hltCaloMET.pt>80 && hltPFMET.pt>150 && (hltCaloMETClean.pt<50 && hltCaloMETCleanUsingJetID.pt<50))",
    #    ]
    #counts = count(sel_trigs); pprint(counts, overlap=True)

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


# ______________________________________________________________________________
sel = "(triggerFlags[%i])" % ireftrig
#sel_noNoise  = "(metfilterFlags[%i] && event.json)" % len(metfilters)
sel_noNoise  = "(metfilterFlags[%i] && event.json)" % len(metfilters)

sel2 = "(triggerFlags[%i] && event.run==207454 && (79<=event.lumi && event.lumi<=600))" % ireftrig2

triggers2012D = OrderedDict([
    ("PFMET150", "(hltCaloMET.pt>80 && hltPFMET.pt>150)"),
    ("MonoCentralJet", "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"),
    ("DiCentralJetHIG", "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"),
    ("DiCentralJetSUS", "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"),
    ("HT", "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"),
    ("btag", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"),
    #("btag_ctrl", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"),
    #("VBFAll", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"),
    #("VBFLead", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_leadmjj>600 && hltPFGlobal.vbf_leadmjj_deta>3.5 && hltPFMETNoMu.pt>65)"),
    ])

triggers2012D_half = OrderedDict([
    ("PFMET150", "(hltCaloMET.pt>80 && hltPFMET.pt>200)"),
    ("MonoCentralJet", "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>135)"),
    ("DiCentralJetHIG", "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>120)"),
    ("DiCentralJetSUS", "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) )"),
    ("HT", "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"),
    ("btag", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>90)"),
    #("btag_ctrl", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>90)"),
    #("VBFAll", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>75)"),
    #("VBFLead", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_leadmjj>600 && hltPFGlobal.vbf_leadmjj_deta>3.5 && hltPFMETNoMu.pt>75)"),
    ])

triggers2015a = OrderedDict([
    ("PFMET150", "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)"),
    ("PFMET100CJ100", "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)"),
    #("PFMET120CJ80", "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>80 && hltCaloMETClean.pt>70 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>120||hltPFMET.pt>120) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)"),
    ("PFMET100CJ60CJ30", "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFGlobal.dijet_mindphi_2cj>0.5)"),
    ("PFMET80CJ60x2", "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) && hltPFGlobal.dijet_mindphi_2cj>0.3)"),
    #("PFMET100HT400", "(hltCaloHTMHT.sumEt>350 && hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && hltPFHTMHTNoPU.sumEt>400 && hltPFMET.pt>100 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)"),
    ("PFMET80HT400", "(hltCaloHTMHT.sumEt>350 && hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && hltPFHTMHTNoPU.sumEt>400 && hltPFMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)"),
    ("PFMET80CJ30x2CSV07", "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && hltPFMET.pt>80)"),
    #("PFMET100CJ80CSV07", "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100)"),
    #("PFMET80CJ30x2ctrl", "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && hltPFMET.pt>80)"),
    #("PFMET100CJ80ctrl", "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100)"),
    #("PFMET65VBFMJJ800", "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"),
    ])

benchmarks = OrderedDict([
    ("PFMET150", sel_noNoise + "*(recoPFMETT0T1.pt>150)*" + "(Sum$(patJets.pt>20)>0 && patJets[0].jetID==1 && ((Sum$(patJets.pt>20)>1 && patJets[1].jetID==1) || Sum$(patJets.pt>20)==1))"),
    ("PFMET100CJ100", sel_noNoise + "*(recoPFMETT0T1.pt>100)*" + "((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && abs(deltaPhi(patJets[0].phi,patJets[1].phi))<2.5) || (Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1))"),
    #("PFMET120CJ80", sel_noNoise + "*(recoPFMETT0T1.pt>120)*" + "((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && abs(deltaPhi(patJets[0].phi,patJets[1].phi))<2.5) || (Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1))"),
    ("PFMET100CJ60CJ30", sel_noNoise + "*(recoPFMETT0T1.pt>100)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>60 && abs(patJets[0].eta)<2.5 && patJets[1].pt>30 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_maxpt>110 && patGlobal.dijet_maxpt_mjj<250 && patGlobal.dijet_mindphi_3cj>0.5)"),
    ("PFMET80CJ60x2", sel_noNoise + "*(recoPFMETT0T1.pt>80)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>70 && abs(patJets[0].eta)<2.5 && patJets[1].pt>70 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_mindphi_2cj>0.5 && ((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>2 && abs(patJets[2].eta)<2.5 && patJets[2].jetID==1 && patGlobal.dijet_mindphi_3cj>0.3)||Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==2))"),
    #("PFMET100HT400", sel_noNoise + "*(recoPFMETT0T1.pt>100)*" + "(Sum$(patJets.pt>50 && abs(patJets.eta)<2.5)>2 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[2].pt>50 && abs(patJets[2].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patJets[2].jetID==1 && patHTMHT.sumEt>450)"),
    ("PFMET80HT400", sel_noNoise + "*(recoPFMETT0T1.pt>80)*" + "(Sum$(patJets.pt>50 && abs(patJets.eta)<2.5)>2 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[2].pt>50 && abs(patJets[2].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patJets[2].jetID==1 && patHTMHT.sumEt>450)"),
    ("PFMET80CJ30x2CSV07", sel_noNoise + "*(recoPFMETT0T1.pt>80)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.bjet_maxcsv>0.898)"),
    #("PFMET100CJ80CSV07", sel_noNoise + "*(recoPFMETT0T1.pt>100)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>0 && patJets[0].pt>100 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1 && patGlobal.bjet_maxcsv>0.898)"),
    #("PFMET80CJ30x2ctrl", sel_noNoise + "*(recoPFMETT0T1.pt>80)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1)"),
    #("PFMET100CJ80ctrl", sel_noNoise + "*(recoPFMETT0T1.pt>100)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>0 && patJets[0].pt>100 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1)"),
    #("PFMET65VBFMJJ800", sel_noNoise + "*(recoPFMETT0T1.pt>65)*" + "(Sum$(patJets.pt>30 && abs(patJets.eta)<4.7)>1 && patJets[0].pt>50 && abs(patJets[0].eta)<4.7 && patJets[1].pt>50 && abs(patJets[1].eta)<4.7 && (patJets[0].eta*patJets[1].eta)<=0 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.vbf_maxmjj>900 && patGlobal.vbf_maxmjj_deta>3.5)"),
    ])


if sections["future_triggers"]:
    ## PFMET150
    #sel_trig0 = triggers2012D["PFMET150"]
    #sel_trig1 = triggers2015a["PFMET150"]

    ## MonoCentralJet
    #sel_trig0 = triggers2012D["MonoCentralJet"]
    #sel_trig1 = triggers2015a["PFMET100CJ100"]
    ##sel_trig1 = triggers2015a["PFMET120CJ80"]

    ## DiCentralJetHIG
    #sel_trig0 = triggers2012D["DiCentralJetHIG"]
    #sel_trig1 = triggers2015a["PFMET100CJ60CJ30"]

    ## DiCentralJetSUS
    #sel_trig0 = triggers2012D["DiCentralJetSUS"]
    #sel_trig1 = triggers2015a["PFMET80CJ60x2"]

    ## Multijet
    #sel_trig0 = triggers2012D["HT"]
    #sel_trig1 = triggers2015a["PFMET100HT400"]
    ##sel_trig1 = triggers2015a["PFMET80HT400"]

    ## btag
    sel_trig0 = triggers2012D["btag"]
    sel_trig1 = triggers2015a["PFMET80CJ30x2CSV07"]
    #sel_trig1 = triggers2015a["PFMET100CJ80CSV07"]

    # Count!
    sel_trigs = [
        sel_trig0 + "||" + sel_trig1,
        sel_trig0,
        sel_trig1,
        sel_trig0 + "&&" + sel_trig1,
        ]
    counts = count(sel_trigs); pprint(counts, overlap=True)


if sections["future_numbers"]:
    gROOT.LoadMacro("Helper_counter_METTriggers.h")

    htemp2 = TH1F("htemp2", "", 4, 0, 4)
    tree2.Project("htemp2", sel_noNoise, "*".join([sel2]), "goff")
    refcount = float(htemp2.Integral())

    #refscale = 101411.0 / 3120.0
    #refrate = 3.22
    refscale = 38104.0 / 3123.3 # or 3124.0?
    refrate = 3.13
    normalizeby = "(normalizeby_nGoodPV(event.nGoodPV) * normalizeby_noiseType2(metfilterFlags[0], metfilterFlags[8]&&optmetfilterFlags[5]) * %.3f)" % refscale

    # All triggers
    counts = count_future(triggers2012D, sel, normalizeby); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2012D_half, sel, normalizeby); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2015a, sel, normalizeby); pprint_future(counts, refcount, refrate)

    # All triggers + benchmarks
    counts = count_future(triggers2012D, sel, normalizeby, benchmarks=benchmarks); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2012D_half, sel, normalizeby, benchmarks=benchmarks); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2015a, sel, normalizeby, benchmarks=benchmarks); pprint_future(counts, refcount, refrate)

    # Exclusive triggers
    counts = count_future(triggers2012D, sel, normalizeby, exclusive=True); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2012D_half, sel, normalizeby, exclusive=True); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2015a, sel, normalizeby, exclusive=True); pprint_future(counts, refcount, refrate)

    # Exclusive triggers + benchmarks
    counts = count_future(triggers2012D, sel, normalizeby, benchmarks=benchmarks, exclusive=True); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2012D_half, sel, normalizeby, benchmarks=benchmarks, exclusive=True); pprint_future(counts, refcount, refrate)
    counts = count_future(triggers2015a, sel, normalizeby, benchmarks=benchmarks, exclusive=True); pprint_future(counts, refcount, refrate)




### Involve drawing below ######################################################

# ______________________________________________________________________________
# Functions

# Book
def prepare(variable, normalize=0):
    if normalize==0:
        params = [
            (variable[0]+"_all", kBlack, kWhite, variable[1], sel + "*" + triggers2012D["PFMET150"]),
            (variable[0]+"_sel",      2, kWhite, variable[1], sel2),
            ]
    elif normalize==1:
        params = [
            (variable[0]+"_all", kBlack, kWhite, variable[1], sel + "*" + triggers2012D["PFMET150"] + "* normalizeby_nGoodPV(event.nGoodPV)"),
            (variable[0]+"_sel",      2, kWhite, variable[1], sel2),
            ]
    elif normalize==2:
        params = [
            #(variable[0]+"_all", kBlack, kWhite, variable[1], sel + "*" + triggers2012D["PFMET150"] + "* normalizeby_nGoodPV(event.nGoodPV) * normalizeby_noiseType(metfilterFlags[0], metfilterFlags[1], (metfilterFlags[3]&&optmetfilterFlags[3]), (metfilterFlags[4]&&metfilterFlags[7]&&optmetfilterFlags[4]))" % (ireftrig, ireftrig2)),
            (variable[0]+"_all", kBlack, kWhite, variable[1], sel + "*" + triggers2012D["PFMET150"] + "* normalizeby_nGoodPV(event.nGoodPV) * normalizeby_noiseType2(metfilterFlags[0], metfilterFlags[8]&&optmetfilterFlags[5])" % (ireftrig, ireftrig2)),
            (variable[0]+"_sel",      2, kWhite, variable[1], sel2),
            ]
    return params

def book(params, binning):
    histos = []
    for i, p in enumerate(params):
        h = TH1F("h_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3])
        h.SetLineWidth(2)
        h.SetLineColor(p[1])
        h.SetMarkerColor(p[1])
        h.SetMarkerSize(0)
        h.SetFillColor(p[2])
        histos.append(h)
    return histos

# Project
def fixOverflow(h):
    nbins = h.GetNbinsX()
    if h.GetBinContent(nbins+1) > 0:
        h.SetBinContent(nbins, h.GetBinContent(nbins) + h.GetBinContent(nbins+1))
        h.SetBinError(nbins, sqrt(h.GetBinError(nbins)**2 + h.GetBinError(nbins+1)**2))
        h.SetBinContent(nbins+1, 0)
        h.SetBinError(nbins+1, 0)
        h.SetEntries(h.GetEntries() - 2)  # SetBinContent() automatically increases # entries by one

def project(params, histos, normalize=-1, drawOverflow=True):
    for i, p in enumerate(params):
        if i%2 == 0:
            tree.Project("h_"+p[0], p[3], p[4], "goff")
        else:
            tree2.Project("h_"+p[0], p[3], p[4], "goff")
        if normalize > 0:
            histos[i].Scale(normalize / histos[i].GetSumOfWeights())
        if drawOverflow:
            fixOverflow(histos[i])
    return

def prepare_canvas_withratio():
    c1 = TCanvas("c1", "c1", 700, 700)
    #c1 = TCanvas("c1", "c1", 600, 600)
    pad1 = TPad("pad1", "top pad"   , 0.0, 0.3, 1.0, 1.0)
    pad2 = TPad("pad2", "bottom pad", 0.0, 0.0, 1.0, 0.3)
    pad1.SetBottomMargin(0.0)
    pad1.Draw()
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(0.35)
    pad2.Draw()
    return (c1, pad1, pad2)

def style_ratio(h1, title):
    h1.Sumw2()
    h1.SetStats(0)
    h1.SetTitle(title)
    h1.SetMaximum(2.2)
    h1.SetMinimum(0)
    h1.SetMarkerSize(0)
    h1.GetXaxis().SetLabelSize(0.12)
    h1.GetXaxis().SetTitleSize(0.14)
    h1.GetXaxis().SetTitleOffset(1.10)
    h1.GetYaxis().CenterTitle()
    h1.GetYaxis().SetLabelSize(0.10)
    h1.GetYaxis().SetTitleSize(0.12)
    h1.GetYaxis().SetTitleOffset(0.6)
    h1.GetYaxis().SetNdivisions(505)

def CMS_label():
    old = (latex.GetTextFont(), latex.GetTextSize())
    latex.SetTextFont(42); latex.SetTextSize(0.026)
    latex.DrawLatex(0.665, 0.968, "Run2012D HLT_L1ETM40_v2")
    latex.SetTextFont(62); latex.SetTextSize(0.028)
    latex.DrawLatex(0.445, 0.968, "CMS Preliminary")
    latex.SetTextFont(old[0]); latex.SetTextSize(old[1])
    return

def draw_withratio(c1, pad1, pad2, histos, labels, rhistos, rlabels, ytitle="Events", logy=False, ymax=8200, rymax=8.2, legend=(0.52,0.78,0.92,0.94)):
    histos[0].SetMaximum(ymax)
    histos[0].GetYaxis().SetTitle(ytitle)
    rhistos[0].SetMaximum(rymax)

    leg1 = TLegend(*legend)
    leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
    leg1.AddEntry(histos[0], labels[0], "l")
    leg1.AddEntry(histos[1], labels[1], "l")
    leg2 = TLegend(legend[0], 0.84, legend[2], legend[3])
    leg2.SetFillStyle(0); leg2.SetLineColor(0); leg2.SetShadowColor(0); leg2.SetBorderSize(0)
    leg2.AddEntry(rhistos[0], rlabels[0], "lpf")

    pad1.cd()
    histos[0].Draw("hist")
    histos[1].Draw("hist same")
    leg1.Draw()
    #CMS_label()
    pad2.cd()
    rhistos[0].Draw()
    leg2.Draw()
    line.DrawLine(histos[0].GetBinLowEdge(1), 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    return (leg1, leg2)

def save_withratio(c1, pad1, pad2, imgdir, imgname):
    pad1.cd()
    gPad.RedrawAxis()
    gPad.Modified()
    gPad.Update()
    pad2.cd()
    gPad.RedrawAxis()
    gPad.Modified()
    gPad.Update()
    c1.cd()
    gPad.Print(imgdir+imgname+".pdf")
    gPad.Print(imgdir+imgname+".png")


#_______________________________________________________________________________
# Determine the normalization to R207454
if sections["future_numbers2"]:
    gROOT.LoadMacro("Helper_counter_METTriggers.h")

    # NPV
    variable = ("nGoodPV", "event.nGoodPV")
    binning = ("#scale[0.7]{RECO} # good PV", 40, 0, 40)
    params = prepare(variable, normalize=0)
    histos = book(params, binning)
    project(params, histos)

    h1 = histos[0]
    h2 = histos[1]
    print h1.Integral(), h2.Integral()
    scale = h2.Integral()/h1.Integral()
    h1.Scale(scale)
    print h1.Integral(), h2.Integral()

    h2by1 = h2.Clone("ratio")
    h2by1.Divide(h2, h1, 1., 1., "")
    title = ";%s ; Ratio" % binning[0]
    style_ratio(h2by1, title)
    #h2by1.SetLineColor(6)

    # Draw
    (c1, pad1, pad2) = prepare_canvas_withratio()
    #label1 = "PFMET150 #splitline{(All runs}{from L1 passthrough)}"
    label1 = "PFMET150 (passthru x%.1f)" % (scale)
    label2 = "PFMET150 (Run 207454)"
    label2by1 = "(Run 207454)/(passthru)"
    (leg1, leg2) = draw_withratio(c1, pad1, pad2, histos, [label1, label2], [h2by1], [label2by1], ymax=7000, rymax=4.7)
    save_withratio(c1, pad1, pad2, imgdir, "future_R207454_normalizeby_nGoodPV")
    for i in xrange(42):  print "%.4f," % h2by1.GetBinContent(i),
    print


    # IsNoise
    variable = ("noiseType", "1.0*(!metfilterFlags[0]) + 2.0*(metfilterFlags[0] && !metfilterFlags[1]) + 3.0*(metfilterFlags[0] && metfilterFlags[1] && !(metfilterFlags[3]&&optmetfilterFlags[3])) + 4.0*(metfilterFlags[0] && metfilterFlags[1] && metfilterFlags[3] && optmetfilterFlags[3] && !(metfilterFlags[4]&&metfilterFlags[7]&&optmetfilterFlags[4]))")
    #variable = ("noiseType", "1.0*(!metfilterFlags[0]) + 2.0*(metfilterFlags[0] && !(metfilterFlags[8]&&optmetfilterFlags[5]))")
    binning = ("#scale[0.7]{RECO} Type of noise", 7, 0, 7)
    params = prepare(variable, normalize=1)
    histos = book(params, binning)
    project(params, histos)

    h1 = histos[0]
    h2 = histos[1]
    scale = h2.Integral()/h1.Integral()
    h1.Scale(scale)
    print h1.Integral(), h2.Integral()

    h2by1 = h2.Clone("ratio")
    h2by1.Divide(h2, h1, 1., 1., "")
    title = ";%s ; Ratio" % binning[0]
    style_ratio(h2by1, title)
    #h2by1.SetLineColor(6)
    h2by1.GetXaxis().SetBinLabel(1, "no noise")
    h2by1.GetXaxis().SetBinLabel(2, "HBHE")
    h2by1.GetXaxis().SetBinLabel(3, "Halo")
    h2by1.GetXaxis().SetBinLabel(4, "ECAL")
    h2by1.GetXaxis().SetBinLabel(5, "TRK")
    h2by1.GetXaxis().SetBinLabel(6, "--")
    h2by1.GetXaxis().SetBinLabel(7, "--")


    # Draw
    (c1, pad1, pad2) = prepare_canvas_withratio()
    #label1 = "PFMET150 #splitline{(All runs}{from L1 passthrough)}"
    label1 = "PFMET150 (passthru x%.1f)" % (scale)
    label2 = "PFMET150 (Run 207454)"
    label2by1 = "(Run 207454)/(passthru)"
    (leg1, leg2) = draw_withratio(c1, pad1, pad2, histos, [label1, label2], [h2by1], [label2by1], ymax=60000, rymax=3.2)
    save_withratio(c1, pad1, pad2, imgdir, "future_R207454_normalizeby_noiseType")
    for i in xrange(7):  print "%.4f," % h2by1.GetBinContent(i),
    print

