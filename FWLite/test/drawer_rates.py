from ROOT import TH1F, TFile, TCanvas, TLegend, TColor, TLatex, gROOT, gInterpreter, gStyle, gPad, kRed, kBlue, kGreen
from math import sqrt


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

kRed2 = TColor.GetColor("#C02020")
kGreen2 = TColor.GetColor("#20C020")
kBlue2 = TColor.GetColor("#2020C0")
kCyan2 = TColor.GetColor("#20C0C0")
kMagenta2 = TColor.GetColor("#C020C0")
kYellow2 = TColor.GetColor("#C0C020")

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
runsel1 = "(event.run == 207454 && 139 <= event.lumi && event.lumi < 400)"
runsel2 = "(event.run == 207454 && 400 <= event.lumi && event.lumi < 1200)"
runsel3 = "(event.run == 207454 && 1200 <= event.lumi && event.lumi <= 1880)"

wait = False
if not wait:  gROOT.SetBatch(1)

#tfile = TFile.Open("../bin/compactified.root")
tfile = TFile.Open("../bin/compactified.0.root")
#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
tree = tfile.Events

imgdir = "figures/"
if not imgdir.endswith("/"):  imgdir += "/"
imgname = "rates_MET_R207454"


# Rates
h_rates_obs = TH1F("h_rates_obs", "; ; count", 30, 0, 30)
h_rates_exp = TH1F("h_rates_exp", "; ; count", 30, 0, 30)

htemp1 = TH1F("htemp1", ";", 10, 0, 10)
htemp2 = TH1F("htemp2", ";", 10, 0, 10)

ii = 1
for i, trig in enumerate(triggers):
    tree.Draw("triggerFlags[%i] >> htemp1" % i, runsel, "goff")
    tree.Draw("triggerFlags[%i] && triggerFlags[%i] >> htemp2" % (i, i_reftrig), runsel, "goff")

    print trig
    print "    %7i   %7i +/- %7i" % (htemp1.GetBinContent(2), htemp2.GetBinContent(2) * ps_reftrig, sqrt(htemp2.GetBinContent(2)) * ps_reftrig)
    print

    if not trig.startswith("HLT_L1ETM") and not trig.startswith("HLT_MET"):
        h_rates_obs.SetBinContent(ii, htemp1.GetBinContent(2))
        #h_rates_obs.SetBinError(ii, 0)
        h_rates_obs.GetXaxis().SetBinLabel(ii, trig)
        h_rates_exp.SetBinContent(ii, htemp2.GetBinContent(2) * ps_reftrig)
        h_rates_exp.SetBinError(ii, sqrt(htemp2.GetBinContent(2)) * ps_reftrig)
        h_rates_exp.GetXaxis().SetBinLabel(ii, trig)
        ii += 1

c1 = TCanvas("c1", "c1", 1100, 700)
#c1.SetLeftMargin(0.1)
c1.SetRightMargin(0.2)
c1.SetBottomMargin(0.4)
h_rates_exp.SetMaximum(h_rates_exp.GetMaximum() * 1.3)
h_rates_exp.SetMarkerSize(0)
h_rates_exp.SetLineColor(kRed)
h_rates_exp.SetLineWidth(2)
h_rates_exp.GetXaxis().SetLabelSize(0.04)
h_rates_exp.GetYaxis().SetLabelSize(0.04)
h_rates_exp.GetYaxis().SetTitleOffset(0.9)
h_rates_exp.GetYaxis().SetTitleSize(0.05)
h_rates_exp.Draw("hist")
h_rates_exp.GetXaxis().SetRangeUser(0,10)
h_rates_exp.Draw("same e1")
h_rates_obs.Draw("same p")

leg = TLegend(0.82, 0.74, 0.96, 0.96, "MET Run 207454")
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.SetShadowColor(0)
leg.AddEntry(h_rates_obs, "observed", "p")
leg.AddEntry(h_rates_exp, "predicted")
leg.Draw()

gPad.Print(imgdir + imgname + ".png")
gPad.Print(imgdir + imgname + ".pdf")
