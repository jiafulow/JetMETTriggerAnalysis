from ROOT import TH1F, TH2F, TProfile, TFile, TCanvas, TLegend, TColor, TLatex, gROOT, gInterpreter, gStyle, gPad, kRed, kBlue, kGreen
from math import sqrt

class DrawerInit:
    def __init__(self):
        # ROOT
        gROOT.LoadMacro("tdrstyle.C")
        gROOT.ProcessLine("setTDRStyle()")

        #gROOT.ProcessLine(".L ../interface/SimpleCandidate.h");
        gInterpreter.GenerateDictionary("simple::Event", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::Particle", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::XYVector", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::LorentzVector", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::Jet", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::CaloJet", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::PFJet", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::MET", "../interface/SimpleCandidate.h")

        gStyle.SetEndErrorSize(2)


# Colors
kRed2 = TColor.GetColor("#C02020")
kGreen2 = TColor.GetColor("#20C020")
kBlue2 = TColor.GetColor("#2020C0")
kCyan2 = TColor.GetColor("#20C0C0")
kMagenta2 = TColor.GetColor("#C020C0")
kYellow2 = TColor.GetColor("#C0C020")

# Text
latex = TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.03)

# Configurations
triggers = [
    "HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4",
    "HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9",
    "HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9",
    "HLT_PFMET150_v7",
    "HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5",
    "HLT_DiCentralPFJet30_PFMET80_v6",
    "HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5",
    "HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4",
    "HLT_PFNoPUHT350_PFMET100_v4",
    "HLT_PFNoPUHT400_PFMET100_v4",
    "HLT_MET120_v13",
    "HLT_MET120_HBHENoiseCleaned_v6",
    "HLT_MET200_v12",
    "HLT_MET200_HBHENoiseCleaned_v5",
    "HLT_MET300_HBHENoiseCleaned_v5",
    "HLT_MET400_v7",
    "HLT_MET400_HBHENoiseCleaned_v5",
    "HLT_L1ETM30_v2",
    "HLT_L1ETM40_v2",
    "HLT_L1ETM70_v2",
    "HLT_L1ETM100_v2",
    ]

reftrig = "HLT_L1ETM40_v2"
i_reftrig = triggers.index(reftrig)
ps_reftrig = 1500

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

runsel = "(event.run == 207454 && 139 <= event.lumi && event.lumi <= 1880)"
runsel1 = "(event.run == 207454 && 1200 <= event.lumi && event.lumi <= 1880)"
runsel2 = "(event.run == 207454 && 400 <= event.lumi && event.lumi < 1200)"
runsel3 = "(event.run == 207454 && 139 <= event.lumi && event.lumi < 400)"

imgdir = "figures/"
if not imgdir.endswith("/"):  imgdir += "/"
