# Import
from ROOT import TH1, TH1F, TH2F, TProfile, TFile, TCanvas, TLegend, TColor, TLatex, TLine, gROOT, gInterpreter, gStyle, gPad, kWhite, kGray, kBlack, kRed, kBlue, kGreen, kCyan, kMagenta, kYellow
from math import sqrt
from array import array

# For init
class DrawerInit:
    def __init__(self):
        # ROOT
        gROOT.LoadMacro("tdrstyle.C")
        gROOT.LoadMacro("HelperFunctions.h")
        gROOT.ProcessLine("setTDRStyle()")

        #gROOT.ProcessLine(".L ../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::Event", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::Particle", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::XYVector", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::LorentzVector", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::Jet", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::CaloJet", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::PFJet", "../interface/SimpleCandidate.h")
        gInterpreter.GenerateDictionary("simple::MET", "../interface/SimpleCandidate.h")

        gStyle.SetEndErrorSize(2)
        gStyle.SetLabelSize(0.04, "Y")

def fixOverflow(h):
    nbins = h.GetNbinsX()
    if h.GetBinContent(nbins+1) > 0:
        h.SetBinContent(nbins, h.GetBinContent(nbins) + h.GetBinContent(nbins+1))
        h.SetBinError(nbins, sqrt(h.GetBinError(nbins)**2 + h.GetBinError(nbins+1)**2))
        h.SetBinContent(nbins+1, 0)
        h.SetBinError(nbins+1, 0)

def getMaximum(histos):
    ymax = histos[0].GetMaximum()
    for h in histos[1:]:
        ymax = max(ymax, h.GetMaximum())
    return ymax


# Colors
kRed2 = TColor.GetColor("#C02020")
kGreen2 = TColor.GetColor("#20C020")
kBlue2 = TColor.GetColor("#2020C0")
kCyan2 = TColor.GetColor("#20C0C0")
kMagenta2 = TColor.GetColor("#C020C0")
kYellow2 = TColor.GetColor("#C0C020")
kOrange = TColor.GetColor("#FF9900")
kOrange2 = TColor.GetColor("#FFCC33")
kPurple2 = TColor.GetColor("#800080")
kOlive2 = TColor.GetColor("#808000")
kTeal2 = TColor.GetColor("#008080")
kGold2 = TColor.GetColor("#DAA520")
kSalmon2 = TColor.GetColor("#FA8072")

# Text
latex = TLatex()
latex.SetNDC()
latex.SetTextFont(42)
latex.SetTextSize(0.026)

line = TLine()
line.SetLineColor(kGray+2)
line.SetLineStyle(2)

# Configurations
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

imgdir = "figures_20131130/"  # for Torino workshop
if not imgdir.endswith("/"):  imgdir += "/"

sections = {}
sections["overview"] = True
sections["overview_prof"] = True
sections["topology"] = True
sections["topology_online"] = True
sections["puremet"] = True
sections["puremet_eff"] = True


# ______________________________________________________________________________
# Init
drawerInit = DrawerInit()
wait = True
if not wait:  gROOT.SetBatch(1)
TH1.SetDefaultSumw2()

tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
#tfile = TFile.Open("../bin/compactified.0.root")
tree = tfile.Events

# ______________________________________________________________________________
# Overview


if sections["overview"]:

    def prepare(variable):
        sel = "(triggerFlags[%i])" % ireftrig
        sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
        sel_lumiA = "(lumilevel<=1)"
        sel_lumiB = "(lumilevel==2)"
        sel_lumiC = "(lumilevel==3)"

        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel])),
            (variable[0]+""             , kBlack, kBlack , variable[1], "*".join([sel, sel_noNoise])),
            (variable[0]+"_lumiA_nofilt", kBlue , kWhite , variable[1], "*".join([sel, sel_lumiA])),
            (variable[0]+"_lumiA"       , kBlue , kBlue2 , variable[1], "*".join([sel, sel_lumiA, sel_noNoise])),
            (variable[0]+"_lumiB_nofilt", kGreen, kWhite , variable[1], "*".join([sel, sel_lumiB])),
            (variable[0]+"_lumiB"       , kGreen, kGreen2, variable[1], "*".join([sel, sel_lumiB, sel_noNoise])),
            (variable[0]+"_lumiC_nofilt", kRed  , kWhite , variable[1], "*".join([sel, sel_lumiC])),
            (variable[0]+"_lumiC"       , kRed  , kRed2  , variable[1], "*".join([sel, sel_lumiC, sel_noNoise])),
        ]
        return params

    def book(params, binning):
        histos = []
        for i, p in enumerate(params):
            h = TH1F("h_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3])
            #h = TH1F("h_"+p[0], "; "+binning[0]+"; Events", binning[1], binning[2], binning[3])
            h.SetLineWidth(2)
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            h.SetFillColor(p[2])
            histos.append(h)
        return histos

    def project(params, histos, normalize=-1, drawOverflow=True):
        for i, p in enumerate(params):
            tree.Project("h_"+p[0], p[3], p[4], "goff")
            if normalize > 0:
                histos[i].Scale(normalize / histos[i].GetSumOfWeights())
            if drawOverflow:
                fixOverflow(histos[i])
        return

    def draw(params, histos, text, logy=False):
        ymax = getMaximum(histos)
        histos[0].SetMaximum(ymax * 1.5)

        histos[1].SetFillStyle(3003)
        histos[3].SetFillStyle(3003)
        histos[5].SetFillStyle(3003)
        histos[7].SetFillStyle(3003)

        histos[0].Draw("hist")
        histos[1].Draw("hist same")

        histos[6].Draw("hist same")
        histos[7].Draw("hist same")

        histos[4].Draw("hist same")
        histos[5].Draw("hist same")

        histos[2].Draw("hist same")
        histos[3].Draw("hist same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label(histos, legend=(0.26,0.74,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[0]+0.5*(legend[2]-legend[0]), legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "all  #mu,#sigma = %.1f,%.1f" % (histos[0].GetMean(), histos[0].GetRMS()), "f")
        leg1.AddEntry(histos[2], "high #mu,#sigma = %.1f,%.1f" % (histos[2].GetMean(), histos[2].GetRMS()), "f")
        leg1.AddEntry(histos[4], "med  #mu,#sigma = %.1f,%.1f" % (histos[4].GetMean(), histos[4].GetRMS()), "f")
        leg1.AddEntry(histos[6], "low  #mu,#sigma = %.1f,%.1f" % (histos[6].GetMean(), histos[6].GetRMS()), "f")
        leg1.Draw()

        leg2 = TLegend(legend[0]+0.5*(legend[2]-legend[0]), legend[1], legend[2], legend[3])
        leg2.SetFillStyle(0)
        leg2.SetLineColor(0)
        leg2.SetShadowColor(0)
        leg2.SetBorderSize(0)
        leg2.AddEntry(histos[1], "(no noise) #mu,#sigma = %.1f,%.1f" % (histos[1].GetMean(), histos[1].GetRMS()), "f")
        leg2.AddEntry(histos[3], "(no noise) #mu,#sigma = %.1f,%.1f" % (histos[3].GetMean(), histos[3].GetRMS()), "f")
        leg2.AddEntry(histos[5], "(no noise) #mu,#sigma = %.1f,%.1f" % (histos[5].GetMean(), histos[5].GetRMS()), "f")
        leg2.AddEntry(histos[7], "(no noise) #mu,#sigma = %.1f,%.1f" % (histos[7].GetMean(), histos[7].GetRMS()), "f")
        leg2.Draw()

        return (leg1, leg2)

    def save(imgname):
        gPad.Print(imgdir+imgname+".pdf")
        gPad.Print(imgdir+imgname+".png")

    # nGoodPV
    variable = ("nGoodPV", "event.nGoodPV")
    binning = ("#scale[0.7]{RECO} # good PVs", 40, 0, 40)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    # nJets
    variable = ("nJets", "Sum$(patJets.pt > 30 && abs(patJets.eta) < 2.5)")
    binning = ("#scale[0.7]{RECO} # jets_{p_{T} > 30 GeV, |#eta| < 2.5}", 5, 0, 5)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    # HLT MET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    # RECO MET
    variable = ("recoPFMET", "recoPFMET.pt")
    binning = ("#scale[0.7]{RECO} PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    variable = ("recoPFMETT1", "recoPFMETT1.pt")
    binning = ("#scale[0.7]{RECO} Type-1 PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    # RECO MET variants
    variable = ("patMPT", "patMPT.pt")
    binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])

    variable = ("recoPFMETNoPU", "recoPFMETNoPU.pt")
    binning = ("#scale[0.7]{RECO} NoPU PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("overview_"+variable[0])


if sections["overview_prof"]:

    def book(params, binning):
        histos = []
        for i, p in enumerate(params):
            h = TProfile("p_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
            h.SetLineWidth(2)
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            #h.SetFillColor(p[2])
            histos.append(h)
        return histos

    def project(params, histos):
        for i, p in enumerate(params):
            tree.Project("p_"+p[0], p[3], p[4], "prof goff")
        return

    def draw(params, histos, text, logy=False):
        histos[0].SetMinimum(1e-6)
        histos[0].SetMinimum(1e2)

        histos[0].Draw()
        for h in histos[1:]:
            h.Draw("same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label_HLT(histos, legend=(0.26,0.79,0.79,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "CaloMET (no noise)")
        leg1.AddEntry(histos[1], "PFMET (no noise)")
        leg1.AddEntry(histos[2], "TrackMET (no noise)")
        leg1.Draw()
        return leg1

    def label_RECO(histos, legend=(0.26,0.74,0.79,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "Type-0+1 PFMET (no noise)")
        leg1.AddEntry(histos[1], "TrackMET (no noise)")
        leg1.AddEntry(histos[2], "MVA PFMET (no noise)")
        leg1.AddEntry(histos[3], "NoPU PFMET (no noise)")
        leg1.Draw()
        return leg1

    def save(imgname):
        gPad.Print(imgdir+imgname+".pdf")
        gPad.Print(imgdir+imgname+".png")

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)

    params = [
        ("hltCaloMET"    , kBlack, kBlack, "hltCaloMET.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("hltPFMET"      , kRed  , kRed  , "hltPFMET.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("hltTrackMET"   , kBlue , kBlue , "hltTrackMET.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PVs; #scale[0.7]{HLT} <MET> [GeV]", 20, 0, 40, 0, 200)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label_HLT(histos)
    save("overview_prof_hltMETs")

    params = [
        ("recoPFMETT0T1", kRed     , kRed     , "recoPFMETT0T1.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("patMPT"       , kBlue    , kBlue    , "patMPT.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("recoPFMETMVA" , kMagenta2, kMagenta2, "recoPFMETMVA.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
        ("recoPFMETNoPU", kCyan2   , kCyan2   , "recoPFMETNoPU.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PVs; #scale[0.7]{RECO} <MET> [GeV]", 20, 0, 40, 0, 200)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label_RECO(histos)
    save("overview_prof_recoMETs")


if sections["topology"]:

    def prepare(variable):
        sel = "(triggerFlags[%i])" % ireftrig
        sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
        sel_topo1 = "(topo_patJets!=0 && (topo_patJets<=1 || topo_patJets==4))"
        sel_topo2 = "(topo_patJets!=0 && (topo_patJets<=2 || topo_patJets==4))"
        sel_topo3 = "(topo_patJets!=0 && (topo_patJets<=3 || topo_patJets==4))"
        sel_topo4 = "(topo_patJets!=0 && topo_patJets==4)"

        params = [
            (variable[0]+"_nofilt", kBlack, kWhite , variable[1], "*".join([sel])),
            (variable[0]+"_topo0" , kBlack, kGray  , variable[1], "*".join([sel, sel_noNoise])),
            (variable[0]+"_topo1" , kBlack, kYellow, variable[1], "*".join([sel, sel_topo1, sel_noNoise])),
            (variable[0]+"_topo2" , kBlack, kOrange, variable[1], "*".join([sel, sel_topo2, sel_noNoise])),
            (variable[0]+"_topo3" , kBlack, kRed   , variable[1], "*".join([sel, sel_topo3, sel_noNoise])),
            (variable[0]+"_topo4" , kBlack, kGreen , variable[1], "*".join([sel, sel_topo4, sel_noNoise])),
        ]
        return params

    def book(params, binning):
        histos = []
        for i, p in enumerate(params):
            h = TH1F("h_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3])
            #h = TH1F("h_"+p[0], "; "+binning[0]+"; Events", binning[1], binning[2], binning[3])
            h.SetLineWidth(2)
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            h.SetFillColor(p[2])
            histos.append(h)
        return histos

    def project(params, histos, normalize=-1, drawOverflow=True):
        for i, p in enumerate(params):
            tree.Project("h_"+p[0], p[3], p[4], "goff")
            if normalize > 0:
                histos[i].Scale(normalize / histos[i].GetSumOfWeights())
            if drawOverflow:
                fixOverflow(histos[i])
        return

    def draw(params, histos, text, logy=False, zoom=False):
        ymax = getMaximum(histos)
        if zoom:
            histos[0].SetMaximum(ymax * 1.5 / 100)
        else:
            histos[0].SetMaximum(ymax * 1.5)

        histos[0].Draw("hist")
        histos[1].Draw("hist same")

        histos[4].Draw("hist same")
        histos[3].Draw("hist same")
        histos[2].Draw("hist same")

        histos[5].Draw("hist same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label(histos, legend=(0.52,0.70,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "Noise", "f")
        leg1.AddEntry(histos[1], "Uncategorized", "f")
        leg1.AddEntry(histos[5], "VBF", "f")
        leg1.AddEntry(histos[2], "MonoCentralJet", "f")
        leg1.AddEntry(histos[3], "DiCentralJets", "f")
        leg1.AddEntry(histos[4], "HT", "f")
        leg1.Draw()
        return (leg1)

    def save(imgname):
        gPad.Print(imgdir+imgname+".pdf")
        gPad.Print(imgdir+imgname+".png")

    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 50, 200)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save("topology_"+variable[0])

    draw(params, histos, "Run2012D HLT_L1ETM40_v2", zoom=True)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save("topology_"+variable[0]+"_zoom")

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 50, 200)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save("topology_"+variable[0])

    draw(params, histos, "Run2012D HLT_L1ETM40_v2", zoom=True)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save("topology_"+variable[0]+"_zoom")


if sections["topology_online"]:
    def prepare(variable):
        sel = "(triggerFlags[%i])" % ireftrig
        sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
        sel_trig1 = ("(triggerFlags[%i])"
                     %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
        sel_trig2 = ("(triggerFlags[%i] || triggerFlags[%i] || triggerFlags[%i] || triggerFlags[%i])"
                     %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5'),
                       triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5'),
                       triggers.index('HLT_DiCentralPFJet30_PFMET80_v6'),
                       triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
        sel_trig3 = ("(triggerFlags[%i] || triggerFlags[%i])"
                     %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'),
                       triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
        sel_trig4 = ("(triggerFlags[%i] || triggerFlags[%i])"
                     %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9'),
                       triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) )
        sel_trig5 = ("(triggerFlags[%i] || triggerFlags[%i])"
                     %(triggers.index('HLT_PFMET150_v7'),
                       triggers.index('HLT_PFMET180_v7')) )

        params = [
            (variable[0]+"_nofilt", kBlack, kWhite , variable[1], "*".join([sel])),
            (variable[0]+"_trig0" , kBlack, kGray  , variable[1], "*".join([sel, sel_noNoise])),
            (variable[0]+"_trig1" , kBlack, kYellow, variable[1], "*".join([sel, "("+("||".join([sel_trig1, sel_trig4, sel_trig5]))+")", sel_noNoise])),
            (variable[0]+"_trig2" , kBlack, kOrange, variable[1], "*".join([sel, "("+("||".join([sel_trig2, sel_trig1, sel_trig4, sel_trig5]))+")", sel_noNoise])),
            (variable[0]+"_trig3" , kBlack, kRed   , variable[1], "*".join([sel, "("+("||".join([sel_trig3, sel_trig2, sel_trig1, sel_trig4, sel_trig5]))+")", sel_noNoise])),
            (variable[0]+"_trig4" , kBlack, kGreen , variable[1], "*".join([sel, "("+("||".join([sel_trig4, sel_trig5]))+")", sel_noNoise])),
            (variable[0]+"_trig5" , kBlack, kBlue  , variable[1], "*".join([sel, "("+("||".join([sel_trig5]))+")", sel_noNoise])),
        ]
        return params

    def book(params, binning):
        histos = []
        for i, p in enumerate(params):
            h = TH1F("h_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3])
            #h = TH1F("h_"+p[0], "; "+binning[0]+"; Events", binning[1], binning[2], binning[3])
            h.SetLineWidth(2)
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            h.SetFillColor(p[2])
            histos.append(h)
        return histos

    def project(params, histos, normalize=-1, drawOverflow=True):
        for i, p in enumerate(params):
            tree.Project("h_"+p[0], p[3], p[4], "goff")
            if normalize > 0:
                histos[i].Scale(normalize / histos[i].GetSumOfWeights())
            if drawOverflow:
                fixOverflow(histos[i])
        return

    def draw(params, histos, text, logy=False, zoom=False):
        ymax = getMaximum(histos)
        if zoom:
            histos[0].SetMaximum(ymax * 1.5 / 100)
        else:
            histos[0].SetMaximum(ymax * 1.5)

        histos[0].Draw("hist")
        histos[1].Draw("hist same")

        histos[4].Draw("hist same")
        histos[3].Draw("hist same")
        histos[2].Draw("hist same")

        histos[5].Draw("hist same")
        histos[6].Draw("hist same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label(histos, legend=(0.52,0.66,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "Noise", "f")
        leg1.AddEntry(histos[1], "Uncollected", "f")
        leg1.AddEntry(histos[5], "VBF", "f")
        leg1.AddEntry(histos[2], "MonoCentralJet", "f")
        leg1.AddEntry(histos[3], "DiCentralJets", "f")
        leg1.AddEntry(histos[4], "HT", "f")
        leg1.AddEntry(histos[6], "Inclusive", "f")
        leg1.Draw()
        return (leg1)

    def save(imgname):
        gPad.Print(imgdir+imgname+".pdf")
        gPad.Print(imgdir+imgname+".png")

    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 50, 200)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.62, "categorized by MET trigger (online)")
    save("topology_trigger_"+variable[0])

    draw(params, histos, "Run2012D HLT_L1ETM40_v2", zoom=True)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.62, "categorized by MET trigger (online)")
    save("topology_trigger_"+variable[0]+"_zoom")

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 50, 200)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.62, "categorized by MET trigger (online)")
    save("topology_trigger_"+variable[0])

    draw(params, histos, "Run2012D HLT_L1ETM40_v2", zoom=True)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.62, "categorized by MET trigger (online)")
    save("topology_trigger_"+variable[0]+"_zoom")


if sections["puremet"]:

    def book(params, binning):
        histos = []
        for i, p in enumerate(params):
            h = TH1F("h_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3])
            #h = TH1F("h_"+p[0], "; "+binning[0]+"; Events", binning[1], binning[2], binning[3])
            h.SetLineWidth(2)
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            h.SetMarkerSize(0)
            h.SetFillColor(p[2])
            histos.append(h)
        return histos

    def project(params, histos, normalize=-1, drawOverflow=True):
        for i, p in enumerate(params):
            tree.Project("h_"+p[0], p[3], p[4], "goff")
            if normalize > 0:
                histos[i].Scale(normalize / histos[i].GetSumOfWeights())
            if drawOverflow:
                fixOverflow(histos[i])
        return

    def draw(params, histos, text, logy=False, zoom=False):
        ymax = histos[2].GetMaximum()
        if zoom:
            histos[0].SetMaximum(ymax * 2.5 / 100)
        else:
            histos[0].SetMaximum(ymax * 2.5)
        histos[0].SetMinimum(0)

        histos[0].Draw("hist")
        for h in histos[1:]:
            h.Draw("hist same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def draw2(params, histos, text, logy=False):
        histos[0].Draw("hist")
        for h in histos[1:]:
            h.Draw("hist same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label(histos, legend=(0.52,0.82,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "Noise", "f")
        leg1.AddEntry(histos[1], "Uncollected by Inclusive", "f")
        leg1.AddEntry(histos[2], "Inclusive", "f")
        leg1.Draw()
        return (leg1)

    def save(imgname):
        gPad.Print(imgdir+imgname+".pdf")
        gPad.Print(imgdir+imgname+".png")


    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig  = ("(triggerFlags[%i] || triggerFlags[%i])"
                 %(triggers.index('HLT_PFMET150_v7'),
                   triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    sel_trig2 = "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltPFMET.pt>150)"
    sel_trig3 = "(hltCaloMET.pt>80 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)"
    sel_trig4 = "(hltCaloMET.pt>80 && hltPFMET.pt>150 && hltTrackMET.pt>15)"
    sel_trig5 = "(hltCaloMET.pt>80 && hltPFMET.pt>150 && abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))<0.5)"
    sel_trig11 = "(hltCaloMET.pt>80 && hltTrackMET.pt>90)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 130, 300)
    addsel = "(1)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack, kBlue , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    save("puremet_"+variable[0])

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltPFMET.pt>150)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack, kBlue , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.78, "HLT PFMET > 150")
    save("puremet_"+variable[0])

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite , variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack, kGreen2, variable[1], "*".join([sel, addsel, sel_trig2, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.78, "HLT CaloMET > 80, PFMET > 150")
    save("puremet_"+variable[0])

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite , variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack, kGreen2, variable[1], "*".join([sel, addsel, sel_trig3, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.78, "HLT CaloMET > 80, PFMET > 150")
    save("puremet_"+variable[0])

    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>150 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite , variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack, kGreen2, variable[1], "*".join([sel, addsel, sel_trig4, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.78, "HLT CaloMET > 80, PFMET > 150,")
    latex.DrawLatex(0.56, 0.74, "HBHE+JetID cleaned")
    save("puremet_"+variable[0])

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>150 && hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite , variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack, kGreen2, variable[1], "*".join([sel, addsel, sel_trig5, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label(histos)
    latex.DrawLatex(0.56, 0.78, "HLT CaloMET > 80, PFMET > 150,")
    latex.DrawLatex(0.56, 0.74, "HBHE+JetID cleaned")
    save("puremet_"+variable[0])

    # hltTrackMET_standalone
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 34, 80, 250)
    addsel = "(1)"
    params = [
        (variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack, kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig11", kBlack, kMagenta, variable[1], "*".join([sel, addsel, "("+sel_trig1+"||"+sel_trig11+")", sel_noNoise])),
        (variable[0]+"_trig1", kBlack, kBlue, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1 = TLegend(0.52,0.78,0.96,0.94)
    leg1.SetFillStyle(0)
    leg1.SetLineColor(0)
    leg1.SetShadowColor(0)
    leg1.SetBorderSize(0)
    leg1.AddEntry(histos[0], "Noise", "f")
    leg1.AddEntry(histos[1], "Uncollected by Inclusive", "f")
    leg1.AddEntry(histos[3], "Collected by Inclusive", "f")
    leg1.AddEntry(histos[2], "TrackMET > 90", "f")
    latex.DrawLatex(0.56, 0.74, "HLT CaloMET > 80")
    leg1.Draw()
    save("puremet_trk_"+variable[0])

    #___________________________________________________________________________
    # recoPFMETT0T1
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 36, 0, 360)
    addsel = "(1)"
    params = [
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, "*".join([sel_trig1])])),
        (variable[0]+"_trig2" , kRed    , kRed2, variable[1], "*".join([sel, addsel, "*".join([sel_trig1, sel_trig2])])),
        (variable[0]+"_trig3" , kBlue   , kBlue2, variable[1], "*".join([sel, addsel, "*".join([sel_trig1, sel_trig2, sel_trig3])])),
        (variable[0]+"_trig4" , kMagenta, kMagenta2, variable[1], "*".join([sel, addsel, "*".join([sel_trig1, sel_trig2, sel_trig3, sel_trig4])])),
    ]
    histos = book(params, binning)
    for h in histos:
        h.SetFillStyle(3003)
    project(params, histos)
    draw2(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1 = TLegend(0.26,0.74,0.80,0.94)
    leg1.SetFillStyle(0)
    leg1.SetLineColor(0)
    leg1.SetShadowColor(0)
    leg1.SetBorderSize(0)
    leg1.AddEntry(histos[0], "As is")
    leg1.AddEntry(histos[1], "+ HBHENoiseCleaned (%.0f%% rate)" % ((histos[1].Integral() - histos[0].Integral())/histos[0].Integral() * 100))
    leg1.AddEntry(histos[2], "+ JetIDCleaned (%.0f%% rate)" % ((histos[2].Integral() - histos[1].Integral())/histos[0].Integral() * 100))
    leg1.AddEntry(histos[3], "+ TrackMET > 15 (%.0f%% rate)" % ((histos[3].Integral() - histos[2].Integral())/histos[0].Integral() * 100))
    leg1.Draw()
    save("puremet_"+variable[0])

    # recoPFMETMVA
    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", 36, 0, 360)
    addsel = "(1)"
    params = [
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, "*".join([sel_trig1])])),
        (variable[0]+"_trig2" , kRed    , kRed2, variable[1], "*".join([sel, addsel, "*".join([sel_trig1, sel_trig2])])),
        (variable[0]+"_trig3" , kBlue   , kBlue2, variable[1], "*".join([sel, addsel, "*".join([sel_trig1, sel_trig2, sel_trig3])])),
        (variable[0]+"_trig4" , kMagenta, kMagenta2, variable[1], "*".join([sel, addsel, "*".join([sel_trig1, sel_trig2, sel_trig3, sel_trig4])])),
    ]
    histos = book(params, binning)
    for h in histos:
        h.SetFillStyle(3003)
    project(params, histos)
    draw2(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1.Draw()
    save("puremet_"+variable[0])

    # recoPFMETTOT1 for hltTrackMET standalone
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 36, 0, 360)
    addsel = "(1)"
    params = [
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, "*".join([sel_trig1])])),
        (variable[0]+"_trig11", kMagenta, kWhite, variable[1], "*".join([sel, addsel, "("+sel_trig1+"||"+sel_trig11+")"])),
    ]
    histos = book(params, binning)
    project(params, histos)
    histos[0].SetMaximum(histos[1].GetMaximum())
    draw2(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1.Clear()
    leg1.AddEntry(histos[0], "As is")
    leg1.AddEntry(histos[1], "As is || TrackMET > 90 (+%.0f%% rate)" % ((histos[1].Integral() - histos[0].Integral())/histos[0].Integral() * 100))
    leg1.Draw()
    save("puremet_trk_"+variable[0])

    addsel = "(hltCaloMETClean.pt>50 && hltCaloMETCleanUsingJetID.pt>50)"
    params = [
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, "*".join([sel_trig1])])),
        (variable[0]+"_trig11", kMagenta, kWhite, variable[1], "*".join([sel, addsel, "("+sel_trig1+"||"+sel_trig11+")"])),
    ]
    histos = book(params, binning)
    project(params, histos)
    histos[0].SetMaximum(histos[1].GetMaximum())
    draw2(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1.Clear()
    leg1.AddEntry(histos[0], "As is (HBHE+JetID cleaned)")
    leg1.AddEntry(histos[1], "#splitline{As is (HBHE+JetID cleaned) ||}{TrackMET > 90 (+%.0f%% rate)}" % ((histos[1].Integral() - histos[0].Integral())/histos[0].Integral() * 100))
    leg1.Draw()
    save("puremet_trk_clean_"+variable[0])



if sections["puremet_eff"]:

    def book(params, binning):
        histos = []
        for i, p in enumerate(params):
            #h = TH1F("h_"+p[0], "; "+binning[0], binning[1], binning[2])
            h = TH1F("h_"+p[0], "; "+binning[0]+"; HLT Efficiency", binning[1], binning[2])
            h.SetLineWidth(2)
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            h.SetFillColor(p[2])
            histos.append(h)
        return histos

    def project(params, histos, normalize=-1, drawOverflow=False):
        for i, p in enumerate(params):
            tree.Project("h_"+p[0], p[3], p[4], "goff")
            if normalize > 0:
                histos[i].Scale(normalize / histos[i].GetSumOfWeights())
            if drawOverflow:
                fixOverflow(histos[i])
        return

    def draw(params, histos, text, logy=False, zoom=False):
        for h in histos[1:]:
            h.Divide(h, histos[0], 1, 1, "b")

        ymax = 1.3
        histos[1].SetMaximum(ymax)
        histos[1].SetMinimum(0)

        histos[1].Draw("e1")
        histos[1].Draw("lhist same")
        for h in histos[2:]:
            h.Draw("e1 same")
            h.Draw("lhist same")

        gPad.SetLogy(logy)
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label1(histos, legend=(0.26,0.78,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[1], "As is (covered by next graph)")
        leg1.AddEntry(histos[2], "HBHENoiseCleaned")
        leg1.AddEntry(histos[3], "JetIDCleaned (covered by next graph)")
        leg1.AddEntry(histos[4], "HBHE+JetIDCleaned")
        leg1.Draw()
        return (leg1)

    def label2(histos, legend=(0.26,0.78,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[1], "As is")
        leg1.AddEntry(histos[2], "TrackMET > 15")
        leg1.AddEntry(histos[3], "#Delta#phi(PFMET,TrackMET) < 0.5")
        leg1.AddEntry(histos[4], "TrackMET > 15 && #Delta#phi < 0.5")
        leg1.Draw()
        return (leg1)

    def save(imgname):
        gPad.Print(imgdir+imgname+".pdf")
        gPad.Print(imgdir+imgname+".png")


    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig  = ("(triggerFlags[%i] || triggerFlags[%i])"
                 %(triggers.index('HLT_PFMET150_v7'),
                   triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    sel_trig2 = "(hltCaloMET.pt>80 && hltCaloMETClean.pt>50 && hltPFMET.pt>150)"
    sel_trig3 = "(hltCaloMET.pt>80 && hltCaloMETCleanUsingJetID.pt>50 && hltPFMET.pt>150)"
    sel_trig4 = "(hltCaloMET.pt>80 && hltPFMET.pt>150 && hltTrackMET.pt>15)"
    sel_trig5 = "(hltCaloMET.pt>80 && hltPFMET.pt>150 && abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))<0.5)"
    sel_trig11 = "(hltCaloMET.pt>80 && hltTrackMET.pt>90)"

    # recoPFMETT0T1
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    bins = range(0,80,5) + range(80,120,10) + range(120, 180, 15) + range(180,260,20) + [260, 300, 360]
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", len(bins)-1, array('f',bins))
    addsel = "(1)"
    params = [
        #(variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack  , kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
        (variable[0]+"_trig2" , kPurple2, kWhite, variable[1], "*".join([sel, addsel, sel_trig2, sel_noNoise])),
        (variable[0]+"_trig3" , kOrange2, kWhite, variable[1], "*".join([sel, addsel, sel_trig3, sel_noNoise])),
        (variable[0]+"_trig23", kGreen2 , kWhite, variable[1], "*".join([sel, addsel, sel_trig2, sel_trig3, sel_noNoise])),
    ]

    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label1(histos)
    line.DrawLine(0, 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    save("puremet_eff_"+variable[0])

    # recoPFMETMVA
    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", len(bins)-1, array('f',bins))
    params = [
        #(variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack  , kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
        (variable[0]+"_trig2" , kPurple2, kWhite, variable[1], "*".join([sel, addsel, sel_trig2, sel_noNoise])),
        (variable[0]+"_trig3" , kOrange2, kWhite, variable[1], "*".join([sel, addsel, sel_trig3, sel_noNoise])),
        (variable[0]+"_trig23", kGreen2 , kWhite, variable[1], "*".join([sel, addsel, sel_trig2, sel_trig3, sel_noNoise])),
    ]

    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label1(histos)
    line.DrawLine(0, 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    save("puremet_eff_"+variable[0])

    # recoPFMETT0T1
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", len(bins)-1, array('f',bins))
    params = [
        #(variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack  , kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
        (variable[0]+"_trig4" , kTeal2  , kWhite, variable[1], "*".join([sel, addsel, sel_trig4, sel_noNoise])),
        (variable[0]+"_trig5" , kSalmon2, kWhite, variable[1], "*".join([sel, addsel, sel_trig5, sel_noNoise])),
        (variable[0]+"_trig45", kOlive2 , kWhite, variable[1], "*".join([sel, addsel, sel_trig4, sel_trig5, sel_noNoise])),
    ]

    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label2(histos)
    line.DrawLine(0, 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    save("puremet_eff2_"+variable[0])

    # recoPFMETMVA
    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", len(bins)-1, array('f',bins))
    params = [
        #(variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""       , kBlack  , kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1" , kBlack  , kWhite, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
        (variable[0]+"_trig4" , kTeal2  , kWhite, variable[1], "*".join([sel, addsel, sel_trig4, sel_noNoise])),
        (variable[0]+"_trig5" , kSalmon2, kWhite, variable[1], "*".join([sel, addsel, sel_trig5, sel_noNoise])),
        (variable[0]+"_trig45", kOlive2 , kWhite, variable[1], "*".join([sel, addsel, sel_trig4, sel_trig5, sel_noNoise])),
    ]

    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    legs = label2(histos)
    line.DrawLine(0, 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    save("puremet_eff2_"+variable[0])


    # recoPFMETT0T1 for hltTrackMET_standalone
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", len(bins)-1, array('f',bins))
    params = [
        #(variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""          , kBlack  , kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1"    , kBlack  , kWhite, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
        (variable[0]+"_trig11"   , kMagenta, kWhite, variable[1], "*".join([sel, addsel, sel_trig11, sel_noNoise])),
        (variable[0]+"_trig1or11", kGreen-2, kWhite, variable[1], "*".join([sel, addsel, "("+sel_trig1+"||"+sel_trig11+")", sel_noNoise])),
    ]

    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1 = TLegend(0.26,0.82,0.96,0.94)
    leg1.SetFillStyle(0)
    leg1.SetLineColor(0)
    leg1.SetShadowColor(0)
    leg1.SetBorderSize(0)
    leg1.AddEntry(histos[1], "As is")
    leg1.AddEntry(histos[2], "TrackMET > 90")
    leg1.AddEntry(histos[3], "As is || TrackMET > 90")
    leg1.Draw()
    line.DrawLine(0, 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    save("puremet_eff_trk_"+variable[0])

    # recoPFMETMVA for hltTrackMET_standalone
    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", len(bins)-1, array('f',bins))
    params = [
        #(variable[0]+"_nofilt", kBlack, kWhite, variable[1], "*".join([sel, addsel])),
        (variable[0]+""          , kBlack  , kGray , variable[1], "*".join([sel, addsel, sel_noNoise])),
        (variable[0]+"_trig1"    , kBlack  , kWhite, variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
        (variable[0]+"_trig11"   , kMagenta, kWhite, variable[1], "*".join([sel, addsel, sel_trig11, sel_noNoise])),
        (variable[0]+"_trig1or11", kGreen-2, kWhite, variable[1], "*".join([sel, addsel, "("+sel_trig1+"||"+sel_trig11+")", sel_noNoise])),
    ]

    histos = book(params, binning)
    project(params, histos)
    draw(params, histos, "Run2012D HLT_L1ETM40_v2")
    leg1.Draw()
    line.DrawLine(0, 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    save("puremet_eff_trk_"+variable[0])
