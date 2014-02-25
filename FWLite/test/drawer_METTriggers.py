#!/usr/bin/env python

from ROOT import TH1, TH1F, TH2F, TProfile, TFile, TChain, TCanvas, TLegend, TLatex, TLine, gROOT, gInterpreter, gStyle, gSystem, gPad
from rootcolors import *
from math import sqrt
from operator import itemgetter
import numpy

# For init
class DrawerInit:
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
        'HLT_L1ETM40_v2',  # 8
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
        'HLT_PFMET150_v7', # 11
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
    ]  # 8

optmetfilters = [
    "p_goodVerticesFilter",
    "p_noscraping",
    "p_hcallasereventfilter2012",
    "p_EcalDeadCellBoundaryEnergyFilter",  # should apply?
    "p_tobtecfakesFilters",  # should apply?
    ]

sel = "(triggerFlags[%i])" % ireftrig
sel_noNoise  = "(metfilterFlags[%i] && event.json)" % len(metfilters)
sel_noNoise0 = "(metfilterFlags[%i] && event.json && (Sum$(patJets.pt>20)>0 && patJets[0].jetID==1 && ((Sum$(patJets.pt>20)>1 && patJets[1].jetID==1) || Sum$(patJets.pt>20)==1)) )" % len(metfilters)
sel_noNoise1 = "(metfilterFlags[%i] && event.json && Sum$(patJets.pt>20)>0 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1)" % len(metfilters)
sel_noNoise2 = "(metfilterFlags[%i] && event.json && Sum$(patJets.pt>20 && abs(patJets.eta)<2.5)>1 && patJets[0].jetID==1 && ((abs(patJets[1].eta)<2.5 && patJets[1].jetID==1)||abs(patJets[1].eta)>=2.5) )" % len(metfilters)
sel_noQCD = "(patGlobal.dijet_mindphi_2cj>0.5)"
sel_noQCD1  = "(patGlobal.dijet_mindphi_j30>0.5)"
sel_noQCD2  = "(alphaT(patJets[0].pt, patJets[0].px, patJets[0].py, patJets[1].pt, patJets[1].px, patJets[1].py)>0.55)"


# ______________________________________________________________________________
# Configurations
drawerInit = DrawerInit()
chain = TChain("tree", "tree")
infiles = [
    "../bin/compactified.L1ETM40.4.root",
    ]
for f in infiles:
    chain.Add(f)
tree = chain


sections = {}
sections["overview"]          = True
sections["overview_prof"]     = False
sections["overview_scat"]     = False
sections["topology"]          = False
sections["topology_hlt"]      = False
sections["puremet"]           = False
sections["puremet_eff"]       = False
sections["puremet_clean"]     = False
sections["puremet_clean_eff"] = False
sections["monojet"]           = False
sections["monojet_eff"]       = False
sections["higdijet"]          = False
sections["higdijet_eff"]      = False
sections["susdijet"]          = False
sections["susdijet_eff"]      = False
sections["susdijet_eff2"]     = False
sections["multijet"]          = False
sections["multijet_eff"]      = False
sections["bjet"]              = False
sections["bjet_eff"]          = False
sections["vbf"]               = False
sections["vbf_eff"]           = False

sections["future_topology_hlt"] = False
sections["future_triggers"]     = False

#imgdir = "figures_20131130/"  # for Torino workshop
imgdir = "figures_20140206/"  # for first draft
if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)

# interactive: python -i drawer_METTriggers.py
# batch      : python drawer_METTriggers.py -b
#wait = True
#if not wait:  gROOT.SetBatch(1)

plotting = []
tones = [
    TColor.GetColor("#9CE0EB"),  #? use TColor.GetColor("#CCF5FF"),
    TColor.GetColor("#66D6E5"),
    TColor.GetColor("#20ABBF"),
    TColor.GetColor("#0D8799"),
    TColor.GetColor("#055A66"),
    TColor.GetColor("#002E33"),
    ]
toneC = TColor.GetColor("#D67F20")


# ______________________________________________________________________________
# Functions

# Book
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

def book_ratio(params, binning):
    histos = []
    xlow, xup = binning[1], binning[2]
    tmpbins = [0,10,20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,280,320,360,400,480,600,750,1000]
    tmpbins2 = [xlow] + [x for x in tmpbins if (x > xlow and x < xup)] + [xup]
    for i, p in enumerate(params):
        h = TH1F("h_"+p[0], "; "+binning[0], len(tmpbins2)-1, numpy.array(tmpbins2, dtype=float))
        h.SetLineWidth(2)
        h.SetLineColor(p[1])
        h.SetMarkerColor(p[1])
        #h.SetMarkerSize(0)
        h.SetFillColor(p[2])
        histos.append(h)
    return histos

def book_prof(params, binning, option=""):
    histos = []
    for i, p in enumerate(params):
        h = TProfile("p_"+p[0], "; "+binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
        h.SetLineWidth(2)
        h.SetLineColor(p[1])
        h.SetMarkerColor(p[1])
        histos.append(h)
    if "s" in option:
        for h in histos:
            h.SetErrorOption("s")
    return histos

def book_scat(params, ybinnings, binning):
    histos = []
    for i, (p, ybinning) in enumerate(zip(params, ybinnings)):
        h = TH2F("h2_"+p[0], "; "+binning[0]+"; "+ybinning[0], binning[1], binning[2], binning[3], ybinning[1], ybinning[2], ybinning[3])
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
        tree.Project("h_"+p[0], p[3], p[4], "goff")
        if normalize > 0:
            histos[i].Scale(normalize / histos[i].GetSumOfWeights())
        if drawOverflow:
            fixOverflow(histos[i])
    return

def project_prof(params, histos):
    for i, p in enumerate(params):
        # Error is standard error of the mean, e(J)  =  s(J)/sqrt(L(J))
        tree.Project("p_"+p[0], p[3], p[4], "prof goff")
    return

def project_scat(params, histos):
    for i, p in enumerate(params):
        tree.Project("h2_"+p[0], p[3], p[4], "goff")
    return

def CMS_label():
    old = (latex.GetTextFont(), latex.GetTextSize())
    latex.SetTextFont(42); latex.SetTextSize(0.026)
    latex.DrawLatex(0.665, 0.968, "Run2012D HLT_L1ETM40_v2")
    latex.SetTextFont(62); latex.SetTextSize(0.028)
    latex.DrawLatex(0.445, 0.968, "CMS Preliminary")
    latex.SetTextFont(old[0]); latex.SetTextSize(old[1])
    return

def draw_rate(params, histos, ytitle="Events", logy=False):
    ymax = histos[3].GetMaximum()
    histos[0].SetMaximum(ymax * 2.5)
    histos[0].SetMinimum(0)
    histos[0].GetYaxis().SetTitle(ytitle)
    histos[0].SetFillStyle(0)
    histos[2].SetFillStyle(0)
    histos[1].SetFillStyle(3004)

    histos[0].Draw("hist")
    histos[1].Draw("hist same")
    histos[3].Draw("hist same")
    histos[0].Draw("hist same")
    histos[2].Draw("hist same")
    gPad.SetLogy(logy)
    CMS_label()
    return

def label_rate(histos, trig="Trigger", legend=(0.52,0.82,0.94,0.94), benchmark=False):
    leg1 = TLegend(legend[0], legend[1], legend[0]+0.55*(legend[2]-legend[0]), legend[3])
    leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
    leg1.AddEntry(histos[0], "all", "f")
    leg1.AddEntry(histos[1], "%s" % trig, "f")
    leg1.Draw()
    leg2 = TLegend(legend[0]+0.55*(legend[2]-legend[0]), legend[1], legend[2], legend[3])
    leg2.SetFillStyle(0); leg2.SetLineColor(0); leg2.SetShadowColor(0); leg2.SetBorderSize(0)
    filt = "(no noise)" if not benchmark else "(benchmark)"
    leg2.AddEntry(histos[2], filt, "f")
    leg2.AddEntry(histos[3], filt, "f")
    leg2.Draw()
    latex.DrawLatex(0.54, 0.79, "All HLT filters except this one")
    return (leg1, leg2)

def draw_effnum(params, histos, k, trigs, ytitle="Events", ymin=0, ymax=-2.0, logy=False, legend=(0.46,0.625,0.96,0.94), normalize=False):
    if normalize:
        ytitle = "(normalized)"
        for i, h in enumerate(histos):
            h.Scale(1.0/h.GetSumOfWeights())
            if i%2 == 0:
                h.SetFillStyle(3004)
            else:
                h.SetFillStyle(3005)
    else:
        histos[0].SetFillStyle(3004)

    if ymax < 0:
        ymax1 = histos[k].GetMaximum()
        histos[k].SetMaximum(abs(ymax) * ymax1)
    else:
        histos[k].SetMaximum(ymax)
    histos[k].SetMinimum(ymin)
    histos[k].GetYaxis().SetTitle(ytitle)

    histos[k].Draw("hist")
    for h in histos:
        h.SetLineWidth(1)
        h.Draw("hist same")
    gPad.SetLogy(logy)

    leg1 = TLegend(*legend)
    leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
    for h, t in zip(histos, trigs):
        leg1.AddEntry(h, t, "f")
    leg1.Draw()
    CMS_label()

    tmphistos = []
    if normalize:
        for h in histos:
            hclone = h.Clone(h.GetName() + "_clone")
            hclone.SetFillStyle(0)
            hclone.SetLineWidth(2)
            hclone.Draw("hist same")
            tmphistos.append(hclone)  # persistent
    return (leg1, tmphistos)  # persistent

def draw_eff(params, histos, k, trigs, ytitle="HLT Efficiency", ymin=0, ymax=1.5, logy=False, legend=(0.48,0.695,0.96,0.94)):
    for h in histos[1:]:
        h.Divide(h, histos[0], 1, 1, "b")

    if ymax < 0:
        ymax1 = histos[k].GetMaximum()
        histos[k].SetMaximum(abs(ymax) * ymax1)
    else:
        histos[k].SetMaximum(ymax)
    histos[k].SetMinimum(ymin)
    histos[k].GetYaxis().SetTitle(ytitle)

    histos[k].Draw("hist")
    for h in histos[1:]:
        h.SetLineWidth(1)
        h.Draw("hist same")
    gPad.SetLogy(logy)

    leg1 = TLegend(*legend)
    leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
    for h, t in zip(histos[1:], trigs[1:]):
        leg1.AddEntry(h, t, "f")
    leg1.Draw()
    CMS_label()
    line.DrawLine(histos[0].GetBinLowEdge(1), 1, histos[0].GetBinLowEdge(histos[0].GetNbinsX()+1), 1)
    return (leg1)  # persistent


# Blah
def getMaximum(histos):
    ymax = histos[0].GetMaximum()
    for h in histos[1:]:
        ymax = max(ymax, h.GetMaximum())
    return ymax

def save(imgdir, imgname):
    gPad.RedrawAxis()
    gPad.Print(imgdir+imgname+".pdf")
    gPad.Print(imgdir+imgname+".png")


# ______________________________________________________________________________
# Overview
if sections["overview"]:

    def prepare(variable):
        sel_lumiA = "(lumilevel<=0)"
        sel_lumiB = "(0<lumilevel && lumilevel<=4)"
        sel_lumiC = "(4<lumilevel)"

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

    def draw(params, histos, ytitle="Events", logy=False):
        ymax = getMaximum(histos)
        histos[0].SetMaximum(ymax * 1.5)
        histos[0].GetYaxis().SetTitle(ytitle)

        histos[1].SetFillStyle(3004)
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
        CMS_label()
        return

    def draw_norm(params, histos, ytitle="(normalized)", logy=False):
        ymax = getMaximum(histos)
        histos[0].SetMaximum(ymax * 1.5)
        histos[0].GetYaxis().SetTitle(ytitle)

        histos[0].SetFillStyle(3004)
        histos[1].SetFillStyle(3003)
        histos[0].Draw("hist")
        histos[1].Draw("hist same")

        gPad.SetLogy(logy)
        CMS_label()
        return

    def label(histos, legend=(0.26,0.74,0.96,0.94)):
        leg1 = TLegend(legend[0], legend[1], legend[0]+0.45*(legend[2]-legend[0]), legend[3])
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "all  #mu,#sigma = %.1f,%.1f" % (histos[0].GetMean(), histos[0].GetRMS()), "f")
        leg1.AddEntry(histos[2], "low  #mu,#sigma = %.1f,%.1f" % (histos[2].GetMean(), histos[2].GetRMS()), "f")
        leg1.AddEntry(histos[4], "med  #mu,#sigma = %.1f,%.1f" % (histos[4].GetMean(), histos[4].GetRMS()), "f")
        leg1.AddEntry(histos[6], "high #mu,#sigma = %.1f,%.1f" % (histos[6].GetMean(), histos[6].GetRMS()), "f")
        leg1.Draw()

        leg2 = TLegend(legend[0]+0.45*(legend[2]-legend[0]), legend[1], legend[2], legend[3])
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

    def label_norm(histos, legend=(0.26,0.80,0.96,0.94)):
        leg1 = TLegend(0.26,0.84,0.96,0.94)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "1 #leq #PV #leq 10 (no noise) #mu,#sigma = %.1f,%.1f" % (histos[0].GetMean(), histos[0].GetRMS()), "f")
        leg1.AddEntry(histos[1], "25 #leq #PV #leq 35(no noise) #mu,#sigma = %.1f,%.1f" % (histos[1].GetMean(), histos[1].GetRMS()), "f")
        leg1.Draw()
        return (leg1)

    # nGoodPV
    variable = ("nGoodPV", "event.nGoodPV")
    binning = ("#scale[0.7]{RECO} # good PV", 40, 0, 40)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos); draw(params, histos); legs = label(histos)
    save(imgdir, "overview_"+variable[0])
    #
    ## nJets
    #variable = ("nJets", "Sum$(patJets.pt > 30 && abs(patJets.eta) < 2.5)")
    #binning = ("#scale[0.7]{RECO} # jets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## Jet 1 pT
    #variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    #binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## Jet 2 pT
    #variable = ("ptj2", "Alt$(patJets[1].pt, 0)")
    #binning = ("#scale[0.7]{RECO} jet 2 p_{T} [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## Jet 1 eta
    #variable = ("etaj1", "Alt$(patJets[0].eta, -99)")
    #binning = ("#scale[0.7]{RECO} jet 1 #eta", 20, -5, 5)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## Jet 2 eta
    #variable = ("etaj2", "Alt$(patJets[1].eta, -99)")
    #binning = ("#scale[0.7]{RECO} jet 2 #eta", 20, -5, 5)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## HLT MET
    #variable = ("hltCaloMET", "hltCaloMET.pt")
    #binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("hltPFMET", "hltPFMET.pt")
    #binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("hltTrackMET", "hltTrackMET.pt")
    #binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("hltPFMETphi", "hltPFMET.phi")
    #binning = ("#scale[0.7]{HLT} PFMET #phi", 32, -3.2, 3.2)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## RECO MET
    #variable = ("recoPFMET", "recoPFMET.pt")
    #binning = ("#scale[0.7]{RECO} PFMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("recoPFMETT1", "recoPFMETT1.pt")
    #binning = ("#scale[0.7]{RECO} T1 PFMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    #binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("recoPFMETT0T1phi", "recoPFMETT0T1.phi")
    #binning = ("#scale[0.7]{RECO} T0T1 PFMET #phi", 32, -3.2, 3.2)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## RECO MET variants
    #variable = ("patMPT", "patMPT.pt")
    #binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    #binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("recoPFMETNoPU", "recoPFMETNoPU.pt")
    #binning = ("#scale[0.7]{RECO} NoPU PFMET [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## RECO mindphi
    #variable = ("patmindphi_2cj", "patGlobal.dijet_mindphi_2cj")
    #binning = ("#scale[0.7]{RECO} min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("patmindphi_j30", "patGlobal.dijet_mindphi_j30")
    #binning = ("#scale[0.7]{RECO} min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    ## RECO PFHT, PFMHT
    #variable = ("patHT", "patHTMHT.sumEt")
    #binning = ("#scale[0.7]{RECO} PFHT [GeV]", 30, 0, 300)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])
    #
    #variable = ("patMHT", "patHTMHT.pt")
    #binning = ("#scale[0.7]{RECO} PFMHT [GeV]", 30, 0, 150)
    #params = prepare(variable)
    #histos = book(params, binning)
    #project(params, histos); draw(params, histos); legs = label(histos)
    #save(imgdir, "overview_"+variable[0])


    #___________________________________________________________________________
    ## Shape comparisons
    #sel_lowPU  = "(1<=event.nGoodPV && event.nGoodPV<=10)"
    #sel_highPU = "(25<=event.nGoodPV && event.nGoodPV<=35)"
    #
    ## HLT MET
    #variable = ("hltCaloMET", "hltCaloMET.pt")
    #binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kGray+2, kGray+2, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kBlack, kBlack, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])
    #
    #variable = ("hltPFMET", "hltPFMET.pt")
    #binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kRed, kRed, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kMaroon2, kMaroon2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])
    #
    #variable = ("hltTrackMET", "hltTrackMET.pt")
    #binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kBlue, kBlue, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kNavy2, kNavy2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])
    #
    ## RECO MET
    #variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    #binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kRed, kRed, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kMaroon2, kMaroon2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])
    #
    #variable = ("patMPT", "patMPT.pt")
    #binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kBlue, kBlue, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kNavy2, kNavy2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])
    #
    #variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    #binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kMagenta, kMagenta, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kPurple2, kPurple2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])
    #
    #variable = ("recoPFMETNoPU", "recoPFMETNoPU.pt")
    #binning = ("#scale[0.7]{RECO} NoPU PFMET [GeV]", 30, 0, 150)
    #params = [
    #    (variable[0]+"_lowPU", kCyan, kCyan, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
    #    (variable[0]+"_highPU", kTeal2, kTeal2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    #]
    #histos = book(params, binning)
    #project(params, histos, normalize=1); draw_norm(params, histos); legs = label_norm(histos)
    #save(imgdir, "overview_norm_"+variable[0])


if sections["overview_prof"]:

    def rmserror(binning, histos):
        hhistos = []
        for h in histos:
            hh = TH1F("h_"+h.GetName(), "; "+binning[0], binning[1], binning[2], binning[3])
            hh.SetLineWidth(h.GetLineWidth())
            hh.SetLineColor(h.GetLineColor())
            hh.SetMarkerColor(h.GetMarkerColor())
            #hh.SetFillColor(h.GetFillColor())
            for b in xrange(1, h.GetNbinsX()+1):
                err = h.GetBinError(b)
                n = h.GetBinEffectiveEntries(b)
                # This formula is only valid for normal distribution and unweighted entries
                errerr = sqrt(2.0 * (err ** 4) / (n - 1))
                hh.SetBinContent(b, err)
                hh.SetBinError(b, errerr)
            hhistos.append(hh)
        return hhistos

    def draw(params, histos, ytitle="Events", logy=False, ymin=1e-6, ymax=1e2):
        histos[0].SetMinimum(ymin)
        histos[0].SetMaximum(ymax)
        #histos[0].GetYaxis().SetTitle(ytitle)

        histos[0].Draw()
        for h in histos[1:]:
            h.Draw("same")

        gPad.SetLogy(logy)
        CMS_label()
        return

    def label_HLT(histos, legend=(0.26,0.79,0.79,0.94)):
        leg1 = TLegend(*legend)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "CaloMET (no noise)")
        leg1.AddEntry(histos[1], "PFMET (no noise)")
        leg1.AddEntry(histos[2], "TrackMET (no noise)")
        leg1.Draw()
        return leg1

    def label_RECO(histos, legend=(0.26,0.74,0.79,0.94)):
        leg1 = TLegend(*legend)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "T0T1 PFMET (no noise)")
        leg1.AddEntry(histos[1], "TrackMET (no noise)")
        leg1.AddEntry(histos[2], "MVA PFMET (no noise)")
        leg1.AddEntry(histos[3], "NoPU PFMET (no noise)")
        leg1.Draw()
        return leg1

    # Mean
    params = [
        ("hltCaloMET"    , kBlack, kBlack, "hltCaloMET.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("hltPFMET"      , kRed  , kRed  , "hltPFMET.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("hltTrackMET"   , kBlue , kBlue , "hltTrackMET.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PV; #scale[0.7]{HLT} #mu_{MET} [GeV]", 20, 0, 40, 0, 200)
    histos = book_prof(params, binning)
    project_prof(params, histos)
    draw(params, histos, ymax=80)
    legs = label_HLT(histos)
    save(imgdir, "overview_prof_hltMETs")

    params = [
        ("recoPFMETT0T1", kRed     , kRed     , "recoPFMETT0T1.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("patMPT"       , kBlue    , kBlue    , "patMPT.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("recoPFMETMVA" , kMagenta2, kMagenta2, "recoPFMETMVA.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
        ("recoPFMETNoPU", kCyan2   , kCyan2   , "recoPFMETNoPU.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PV; #scale[0.7]{RECO} #mu_{MET} [GeV]", 20, 0, 40, 0, 200)
    histos = book_prof(params, binning)
    project_prof(params, histos)
    draw(params, histos, ymax=80)
    legs = label_RECO(histos)
    save(imgdir, "overview_prof_recoMETs")

    # Stdev
    params = [
        ("hltCaloMET"    , kBlack, kBlack, "hltCaloMET.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("hltPFMET"      , kRed  , kRed  , "hltPFMET.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("hltTrackMET"   , kBlue , kBlue , "hltTrackMET.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PV; #scale[0.7]{HLT} #sigma_{MET} [GeV]", 20, 0, 40, 0, 200)
    histos = book_prof(params, binning, option="s")
    project_prof(params, histos)
    hhistos = rmserror(binning, histos)
    draw(params, hhistos, ymax=40)
    legs = label_HLT(histos)
    save(imgdir, "overview_prof_sigma_hltMETs")

    params = [
        ("recoPFMETT0T1", kRed     , kRed     , "recoPFMETT0T1.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("patMPT"       , kBlue    , kBlue    , "patMPT.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("recoPFMETMVA" , kMagenta2, kMagenta2, "recoPFMETMVA.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
        ("recoPFMETNoPU", kCyan2   , kCyan2   , "recoPFMETNoPU.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PV; #scale[0.7]{RECO} #sigma_{MET} [GeV]", 20, 0, 40, 0, 200)
    histos = book_prof(params, binning, option="s")
    project_prof(params, histos)
    hhistos = rmserror(binning, histos)
    draw(params, hhistos, ymax=40)
    legs = label_RECO(histos)
    save(imgdir, "overview_prof_sigma_recoMETs")


if sections["overview_scat"]:

    def draw(histo, legend=(0.52,0.88,0.88,0.94), pfmet50=False):
        histo.GetYaxis().SetTitleOffset(1)
        histo.Draw("COLZ")
        gPad.SetLeftMargin(0.13); gPad.SetRightMargin(0.10)
        gPad.Modified(); gPad.Update()
        palette = histo.FindObject("palette")
        xy = (0.91, 0.13, 0.95, 0.95)
        palette.SetX1NDC(xy[0]); palette.SetY1NDC(xy[1]); palette.SetX2NDC(xy[2]); palette.SetY2NDC(xy[3])
        palette.SetTitleSize(0.024); palette.SetLabelSize(0.024)
        gPad.Modified(); gPad.Update()

        leg1 = TLegend(*legend)
        leg1.SetFillColor(0) #leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histo, "(no noise) corr(X,Y) = %.2f" % histo.GetCorrelationFactor(), "f")
        if pfmet50:  latex.DrawLatex(0.58, 0.84, "HLT PFMET>50")
        leg1.Draw()
        CMS_label()
        return (palette, leg1)


    # 2D vs. RECO PFMET
    palettes = []
    params = [
        ("hltCaloMET"    , kBlack, kBlack, "hltCaloMET.pt:recoPFMETT0T1.pt" , "*".join([sel, sel_noNoise])),
        ("hltPFMET"      , kRed  , kRed  , "hltPFMET.pt:recoPFMETT0T1.pt"   , "*".join([sel, sel_noNoise])),
        ("hltTrackMET"   , kBlue , kBlue , "hltTrackMET.pt:recoPFMETT0T1.pt", "*".join([sel, sel_noNoise])),
        ]
    ybinnings = [
        ("#scale[0.7]{HLT} CaloMET [GeV]" , 20, 0, 200),
        ("#scale[0.7]{HLT} PFMET [GeV]"   , 20, 0, 200),
        ("#scale[0.7]{HLT} TrackMET [GeV]", 20, 0, 200),
        ]
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 20, 0, 200)
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    #for p, h in zip(params, histos):
    #    palettes.append(draw(h))
    #    save(imgdir, "overview_scat_recoPFMETT0T1_" + p[0])

    # 2D vs. RECO PFMET phi
    del palettes[:]
    params = [
        ("hltCaloMETphi"    , kBlack, kBlack, "hltCaloMET.phi:recoPFMETT0T1.phi" , "*".join([sel, sel_noNoise])),
        ("hltPFMETphi"      , kRed  , kRed  , "hltPFMET.phi:recoPFMETT0T1.phi"   , "*".join([sel, sel_noNoise])),
        ("hltTrackMETphi"   , kBlue , kBlue , "hltTrackMET.phi:recoPFMETT0T1.phi", "*".join([sel, sel_noNoise])),
        ]
    ybinnings = [
        ("#scale[0.7]{HLT} CaloMET #phi" , 32, -3.2, 3.2),
        ("#scale[0.7]{HLT} PFMET #phi"   , 32, -3.2, 3.2),
        ("#scale[0.7]{HLT} TrackMET #phi", 32, -3.2, 3.2),
        ]
    binning = ("#scale[0.7]{RECO} T0T1 PFMET #phi", 32, -3.2, 3.2)
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    #for p, h in zip(params, histos):
    #    palettes.append(draw(h))
    #    save(imgdir, "overview_scat_recoPFMETT0T1phi_" + p[0])

    # .. and PFMET > 50 GeV
    del palettes[:]
    addsel = "hltPFMET.pt>50"
    params = [
        ("hltCaloMETphi"    , kBlack, kBlack, "hltCaloMET.phi:recoPFMETT0T1.phi" , "*".join([sel, sel_noNoise, addsel])),
        ("hltPFMETphi"      , kRed  , kRed  , "hltPFMET.phi:recoPFMETT0T1.phi"   , "*".join([sel, sel_noNoise, addsel])),
        ("hltTrackMETphi"   , kBlue , kBlue , "hltTrackMET.phi:recoPFMETT0T1.phi", "*".join([sel, sel_noNoise, addsel])),
        ]
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    #for p, h in zip(params, histos):
    #    palettes.append(draw(h, pfmet50=True))
    #    save(imgdir, "overview_scat_pfmet50_recoPFMETT0T1phi_" + p[0])


    # 2D vs. RECO TrackMET
    palettes = []
    params = [
        ("hltTrackMET"   , kBlue , kBlue , "hltTrackMET.pt:patMPT.pt", "*".join([sel, sel_noNoise])),
        ]
    ybinnings = [
        ("#scale[0.7]{HLT} TrackMET [GeV]", 20, 0, 200),
        ]
    binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 20, 0, 200)
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    #for p, h in zip(params, histos):
    #    palettes.append(draw(h))
    #    save(imgdir, "overview_scat_patMPT_" + p[0])

    # 2D vs. RECO PFMET phi
    del palettes[:]
    params = [
        ("hltTrackMETphi"   , kBlue , kBlue , "hltTrackMET.phi:patMPT.phi", "*".join([sel, sel_noNoise])),
        ]
    ybinnings = [
        ("#scale[0.7]{HLT} TrackMET #phi", 32, -3.2, 3.2),
        ]
    binning = ("#scale[0.7]{RECO} TrackMET #phi", 32, -3.2, 3.2)
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    #for p, h in zip(params, histos):
    #    palettes.append(draw(h))
    #    save(imgdir, "overview_scat_patMPTphi_" + p[0])

    # .. and PFMET > 50 GeV
    del palettes[:]
    addsel = "hltPFMET.pt>50"
    params = [
        ("hltTrackMETphi"   , kBlue , kBlue , "hltTrackMET.phi:patMPT.phi", "*".join([sel, sel_noNoise, addsel])),
        ]
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    #for p, h in zip(params, histos):
    #    palettes.append(draw(h, pfmet50=True))
    #    save(imgdir, "overview_scat_pfmet50_patMPTphi_" + p[0])

    # CaloMET vs CaloMHT for recoPFMETT0T1>100
    palettes = []
    addsel = "recoPFMETT0T1.pt>100"
    params = [
        ("hltCaloMHT"    , kBlack, kBlack, "hltCaloHTMHT.pt:hltCaloMET.pt" , "*".join([sel, sel_noNoise])),
        ]
    ybinnings = [
        ("#scale[0.7]{HLT} CaloMHT [GeV]" , 20, 0, 200),
        ]
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 20, 0, 200)
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    for p, h in zip(params, histos):
        palettes.append(draw(h))
        latex.DrawLatex(0.58, 0.84, "RECO T0T1 PFMET>100")
        save(imgdir, "overview_scat_hltCaloMET_" + p[0])

    del palettes[:]
    params = [
        ("hltCaloMHT"    , kBlack, kBlack, "hltCaloHTMHT.pt:hltPFMET.pt" , "*".join([sel, sel_noNoise])),
        ]
    ybinnings = [
        ("#scale[0.7]{HLT} CaloMHT [GeV]" , 20, 0, 200),
        ]
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 20, 0, 200)
    histos = book_scat(params, ybinnings, binning)
    project_scat(params, histos)
    for p, h in zip(params, histos):
        palettes.append(draw(h))
        latex.DrawLatex(0.58, 0.84, "RECO T0T1 PFMET>100")
        save(imgdir, "overview_scat_hltPFMET_" + p[0])



# ______________________________________________________________________________
# Topology
if sections["topology"]:
    pass

if sections["topology_hlt"]:
    if plotting: del plotting[:]

    sel_trig0  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1  = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
    sel_trig2a = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) )
    sel_trig2b = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
    sel_trig2  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4'), triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')))
    sel_trig3  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
    sel_trig4a = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) )
    sel_trig4b = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) )
    sel_trig4  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5'), triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')))
    sel_trig5a = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) )
    sel_trig5b = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) )
    sel_trig5  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9'), triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')))

    def prepare_all(variable, plotting, addsel):
        params = []
        for i, p in enumerate(reversed(plotting)):
            (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
            params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])))
        return params

    def draw(params, histos, ytitle="Events", logy=False, zoom=False):
        ymax = histos[0].GetMaximum()
        if zoom:
            histos[0].SetMaximum(ymax * 1.5 / 100)
        else:
            histos[0].SetMaximum(ymax * 1.5)
        histos[0].SetMinimum(0)
        histos[0].GetYaxis().SetTitle(ytitle)

        histos[0].Draw("hist")
        histos[1].Draw("hist same")
        histos[2].Draw("hist same")

        gPad.SetLogy(logy)
        CMS_label()
        return

    def draw_all(params, histos, ytitle="Events", logy=False, zoom=False):
        ymax = histos[0].GetMaximum()
        if zoom:
            histos[0].SetMaximum(ymax * 1.5 / 100)
        else:
            histos[0].SetMaximum(ymax * 1.5)
        histos[0].SetMinimum(0)
        histos[0].GetYaxis().SetTitle(ytitle)

        histos[0].Draw("hist")
        for h in histos[1:]:
            h.Draw("hist same")

        gPad.SetLogy(logy)
        CMS_label()
        return

    def label(histos, trig="MET trigger", legend=(0.52,0.82,0.96,0.94), noNoise=False, noQCD=False):
        rate0 = (histos[0].Integral() - histos[1].Integral()) / histos[0].Integral() * 100.
        rate1 = (histos[1].Integral() - histos[2].Integral()) / histos[0].Integral() * 100.
        rate2 = (histos[2].Integral()) / histos[0].Integral() * 100.
        leg1 = TLegend(*legend)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "A && !B (%.0f%%)" % rate0, "f")
        leg1.AddEntry(histos[1], "A && B (%.0f%%)" % rate1, "f")
        leg1.AddEntry(histos[2], "!A && B (%.0f%%)" % rate2, "f")
        leg1.Draw()
        latex.DrawLatex(0.56, 0.78, "A=%s  B=existing" % trig)
        if noNoise:  latex.DrawLatex(0.56, 0.75, "after offline noise cleaning")
        if noQCD  :  latex.DrawLatex(0.56, 0.72, "after offline QCD rejection")
        return (leg1)

    def label_all(histos, trig="MET trigger", legend=(0.52,0.70,0.96,0.94), noNoise=False, noQCD=False):
        leg1 = TLegend(*legend)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[7], "Inclusive", "f")
        leg1.AddEntry(histos[6], "MonoCentralJet", "f")
        leg1.AddEntry(histos[5], "DiCentralJet", "f")
        #leg1.AddEntry(histos[4], "DiCentralJet", "f")
        leg1.AddEntry(histos[3], "HT", "f")
        leg1.AddEntry(histos[2], "btag", "f")
        leg1.AddEntry(histos[1], "VBF", "f")
        #leg1.AddEntry(histos[0], "VBF", "f")
        leg1.Draw()
        if noNoise:  latex.DrawLatex(0.56, 0.66, "after offline noise cleaning")
        if noQCD  :  latex.DrawLatex(0.56, 0.62, "after offline QCD rejection")
        return (leg1)

    kColor = kBlue2
    kTrig = "Inclusive"
    sel_trigC  = "(%s &&  %s)" % (sel_trig0, sel_trig0)
    sel_trigD  = "(%s ||  %s)" % (sel_trig0, sel_trig0)
    sel_trigB  = "(%s && !%s)" % (sel_trig0, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trig0, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kYellow3
    kTrig = "MonoCentralJet"
    sel_trigC  = "(%s &&  %s)" % (sel_trig0, sel_trig1)
    sel_trigD  = "(%s ||  %s)" % (sel_trig0, sel_trig1)
    sel_trigB  = "(%s && !%s)" % (sel_trig0, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trig0, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kOrange3
    kTrig = "DiCentralJetHIG"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig2a)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig2a)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kOrange3
    kTrig = "DiCentralJetSUS"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig2b)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig2b)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kCocoa2
    kTrig = "HT"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig3)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig3)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kBoson
    kTrig = "btag"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig4a)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig4a)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kQuark
    kTrig = "VBFAll"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig5a)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig5a)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kQuark
    kTrig = "VBFLead"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig5b)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig5b)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))


    for p in plotting[1:]:  # skip Inclusive
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p

        # HLT PFMET
        variable = ("hltPFMET", "hltPFMET.pt")
        binning = ("#scale[0.7]{HLT} PFMET [GeV]", 50, 0, 250)
        addsel = "(1)"
        params = [
            (variable[0]+"_nofilt_trigA" , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])),
            (variable[0]+"_nofilt_trigAB", kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB])),
            (variable[0]+"_nofilt_trigB" , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB])),
            ]
        histos = book(params, binning)
        project(params, histos); draw(params, histos); legs = label(histos, kTrig)
        save(imgdir, "topology_hlt_nofilt_"+kTrig+"_"+variable[0])

        params = [
            (variable[0]+"_filt_trigA"   , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA , sel_noNoise])),
            (variable[0]+"_filt_trigAB"  , kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB, sel_noNoise])),
            (variable[0]+"_filt_trigB"   , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB , sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos); draw(params, histos); legs = label(histos, kTrig, noNoise=True)
        save(imgdir, "topology_hlt_"+kTrig+"_"+variable[0])

        params = [
            (variable[0]+"_kill_trigA"   , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA , sel_noNoise, sel_noQCD])),
            (variable[0]+"_kill_trigAB"  , kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB, sel_noNoise, sel_noQCD])),
            (variable[0]+"_kill_trigB"   , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB , sel_noNoise, sel_noQCD])),
            ]
        histos = book(params, binning)
        project(params, histos); draw(params, histos); legs = label(histos, kTrig, noNoise=True, noQCD=True)
        save(imgdir, "topology_hlt_kill_"+kTrig+"_"+variable[0])

        # RECO PFMET
        variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
        binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 0, 250)
        addsel = "(1)"
        params = [
            (variable[0]+"_nofilt_trigA" , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])),
            (variable[0]+"_nofilt_trigAB", kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB])),
            (variable[0]+"_nofilt_trigB" , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB])),
            ]
        histos = book(params, binning)
        project(params, histos); draw(params, histos); legs = label(histos, kTrig)
        save(imgdir, "topology_hlt_nofilt_"+kTrig+"_"+variable[0])

        params = [
            (variable[0]+"_filt_trigA"   , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA , sel_noNoise])),
            (variable[0]+"_filt_trigAB"  , kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB, sel_noNoise])),
            (variable[0]+"_filt_trigB"   , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB , sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos); draw(params, histos); legs = label(histos, kTrig, noNoise=True)
        save(imgdir, "topology_hlt_"+kTrig+"_"+variable[0])

        params = [
            (variable[0]+"_kill_trigA"   , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA , sel_noNoise, sel_noQCD])),
            (variable[0]+"_kill_trigAB"  , kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB, sel_noNoise, sel_noQCD])),
            (variable[0]+"_kill_trigB"   , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB , sel_noNoise, sel_noQCD])),
            ]
        histos = book(params, binning)
        project(params, histos); draw(params, histos); legs = label(histos, kTrig, noNoise=True, noQCD=True)
        save(imgdir, "topology_hlt_kill_"+kTrig+"_"+variable[0])


    # All
    params = []

    # HLT PFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 50, 0, 250)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig)
    save(imgdir, "topology_hlt_nofilt_all_"+variable[0])

    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])

    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise, sel_noQCD])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True, noQCD=True)
    save(imgdir, "topology_hlt_kill_all_"+variable[0])

    # RECO PFMET
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 0, 250)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig)
    save(imgdir, "topology_hlt_nofilt_all_"+variable[0])

    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])

    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise, sel_noQCD])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True, noQCD=True)
    save(imgdir, "topology_hlt_kill_all_"+variable[0])


    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig)
    save(imgdir, "topology_hlt_nofilt_all_"+variable[0])

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig)
    save(imgdir, "topology_hlt_nofilt_all_"+variable[0])

    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig)
    save(imgdir, "topology_hlt_nofilt_all_"+variable[0])

    # HLT mindphi
    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{RECO} min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])

    # RECO mindphi
    variable = ("patmindphi_2cj", "patGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{RECO} min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])

    variable = ("patmindphi_j30", "patGlobal.dijet_mindphi_j30")
    binning = ("#scale[0.7]{RECO} min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])

    # RECO alphaT
    variable = ("patalphat", "alphaT(patJets[0].pt, patJets[0].px, patJets[0].py, patJets[1].pt, patJets[1].px, patJets[1].py)")
    binning = ("#scale[0.7]{RECO} #alpha_{T}", 32, 0, 3.2)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])


    # RECO HT, MHT
    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 30, 0, 300)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])

    variable = ("patMHT", "patHTMHT.pt")
    binning = ("#scale[0.7]{RECO} PFMHT [GeV]", 30, 0, 150)
    addsel = "(1)"
    del params[:]
    for i, p in enumerate(reversed(plotting)):
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        params.append((variable[0]+"_nofilt_trigA_%i" % i, kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA, sel_noNoise])))
    histos = book(params, binning)
    project(params, histos); draw_all(params, histos); legs = label_all(histos, kTrig, noNoise=True)
    save(imgdir, "topology_hlt_all_"+variable[0])


# ______________________________________________________________________________
if sections["puremet"]:
    if plotting: del plotting[:]

    kColor = kBlue2
    kTrig = "PFMET150"

    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    sel_trig2 = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 130, 300)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 34, 130, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} PFMET [GeV]", 34, 130, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    #addsel = "(hltPFMETNoMu.pt>105)"
    addsel = "(hltCaloMET.pt>-99 && hltPFMET.pt>150)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # hltCaloMETphi
    variable = ("hltCaloMETphi", "hltCaloMET.phi")
    binning = ("#scale[0.7]{HLT} CaloMET #phi", 32, -3.2, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMETphi
    variable = ("hltPFMETphi", "hltPFMET.phi")
    binning = ("#scale[0.7]{HLT} PFMET #phi", 32, -3.2, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptchf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["puremet_clean"]:
    if plotting: del plotting[:]

    kColor = kGreen2
    kTrig = "PFMET150_clean"

    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    sel_trig2 = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)"

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = sel_trig2
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = sel_trig2
    plotting.append((variable, binning, addsel))

    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig2
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig2
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptchf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")  # 0.02 or 0.1
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig2
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptchf1
    variable = ("hltPFMET_sumptchf1", "hltTrackMET.sumEt/hltPFMET.pt")  # 0.02 or 0.1
    binning = ("#scale[0.7]{HLT} sum track p_{T}/PFMET", 40, 0, 2)
    addsel = sel_trig2
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig2"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig2])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_noNoise_trig2", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig2, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


# ______________________________________________________________________________
if sections["monojet"]:
    if plotting: del plotting[:]

    kColor = kYellow3
    kTrig = "MonoCentralJet"

    sel_bench = "((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && abs(deltaPhi(patJets[0].phi,patJets[1].phi))<2.5) || (Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1))"
    sel_bench = sel_noNoise + "*" + sel_bench
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
    sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 80, 250)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 34, 80, 250)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} PFMET [GeV]", 34, 80, 250)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    #addsel = "(hltPFMETNoMu.pt>105)"
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>-99 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # nhf
    variable = ("hltPFJetsL1FastL2L3_nhf1", "hltPFJetsL1FastL2L3[0].nhf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral HAD)", 40, 0, 1)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<99 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nef, chf, cef, nch, ntot
    variable = ("hltPFJetsL1FastL2L3_nef1", "hltPFJetsL1FastL2L3[0].nef")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral EM)", 40, 0, 1)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>60)"  # add hltCaloMETClean.pt>60
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_chf1", "hltPFJetsL1FastL2L3[0].chf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(charged HAD)", 40, 0, 1)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_cef1", "hltPFJetsL1FastL2L3[0].cef")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(charged EM)", 40, 0, 1)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_nch1", "hltPFJetsL1FastL2L3[0].nch")
    binning = ("#scale[0.7]{HLT} Leading PFJet # charged constits", 40, 0, 40)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_ntot1", "hltPFJetsL1FastL2L3[0].ntot")
    binning = ("#scale[0.7]{HLT} Leading PFJet # constits", 40, 0, 80)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # dphijj
    variable = ("hltPFJetsL1FastL2L3_dphijj", "abs(deltaPhi(hltPFJetsL1FastL2L3[0].phi, hltPFJetsL1FastL2L3[1].phi))")
    binning = ("#scale[0.7]{HLT} Leading-two PFJets #Delta#phi", 32, 0, 3.2)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>60 && Sum$(hltPFJetsL1FastL2L3.pt>30)>1)"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_bench])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_bench])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig, benchmark=True)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["higdijet"]:
    if plotting: del plotting[:]

    kColor = kOrange3
    kTrig = "DiCentralJetHIG"

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>60 && abs(patJets[0].eta)<2.5 && patJets[1].pt>30 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_maxpt>110 && patGlobal.dijet_maxpt_mjj<250 && patGlobal.dijet_mindphi_3cj>0.5)"
    sel_bench = sel_noNoise + "*" + sel_bench

    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) )
    sel_trig1 = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 32, 60, 220)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 32, 60, 220)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT} PFMETCleanUsingJetID [GeV]", 32, 60, 220)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && hltCaloMETClean.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # Dijet
    variable = ("hltCaloJetsL1Fast_maxptjj", "hltCaloGlobal.dijet_maxpt")
    binning = ("#scale[0.7]{HLT} Leading Calo dijet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>0 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_maxptjj_mjj", "hltCaloGlobal.dijet_maxpt_mjj")
    binning = ("#scale[0.7]{HLT} Leading-p_{T} Calo dijet mass [GeV]", 25, 0, 250)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_dphijj", "abs(deltaPhi(hltCaloJetsL1Fast[0].phi,hltCaloJetsL1Fast[1].phi))")
    binning = ("#scale[0.7]{HLT} Calo #Delta#phi(j_{1},j_{2})", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_maxptjj", "hltPFGlobal.dijet_maxpt")
    binning = ("#scale[0.7]{HLT} Leading PF dijet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_maxptjj_mjj", "hltPFGlobal.dijet_maxpt_mjj")
    binning = ("#scale[0.7]{HLT} Leading-p_{T} PF dijet mass [GeV]", 25, 0, 250)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_dphijj", "abs(deltaPhi(hltPFJetsL1FastL2L3[0].phi,hltPFJetsL1FastL2L3[1].phi))")
    binning = ("#scale[0.7]{HLT} PF #Delta#phi(j_{1},j_{2})", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # mindphi
    variable = ("hltCaloJetsL1Fast_mindphi_j30", "hltCaloGlobal.dijet_mindphi_j30")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>-99 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_mindphi_j40", "hltCaloGlobal.dijet_mindphi_j40")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{j40}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_mindphi_cj30", "hltCaloGlobal.dijet_mindphi_cj30")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{cj30}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_mindphi_cj40", "hltCaloGlobal.dijet_mindphi_cj40")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{cj40}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_mindphi_2cj", "hltCaloGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_mindphi_3cj", "hltCaloGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # (PF)
    variable = ("hltPFJetsL1FastL2L3_mindphi_j30", "hltPFGlobal.dijet_mindphi_j30")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_j40", "hltPFGlobal.dijet_mindphi_j40")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{j40}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_cj30", "hltPFGlobal.dijet_mindphi_cj30")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{cj30}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_cj40", "hltPFGlobal.dijet_mindphi_cj40")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{cj40}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_3cj", "hltPFGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nhf
    variable = ("hltPFJetsL1FastL2L3_nhf1", "hltPFJetsL1FastL2L3[0].nhf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral HAD)", 40, 0, 1)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nch
    variable = ("hltPFJetsL1FastL2L3_nch1", "hltPFJetsL1FastL2L3[0].nch")
    binning = ("#scale[0.7]{HLT} Leading PFJet # charged constits", 40, 0, 40)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))


    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_bench])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_bench])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig, benchmark=True)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["susdijet"]:
    if plotting: del plotting[:]

    kColor = kOrange3
    kTrig = "DiCentralJetSUS"

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>70 && abs(patJets[0].eta)<2.5 && patJets[1].pt>70 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_mindphi_2cj>0.5 && ((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>2 && abs(patJets[2].eta)<2.5 && patJets[2].jetID==1 && patGlobal.dijet_mindphi_3cj>0.3)||Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==2))"
    sel_bench = sel_noNoise + "*" + sel_bench

    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
    sel_trig1 = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 32, 60, 220)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>-99||hltPFMETNoMu.pt>-99) )"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 32, 60, 220)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT} PFMETCleanUsingJetID [GeV]", 32, 60, 220)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>-99)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nhf
    variable = ("hltPFJetsL1FastL2L3_nhf1", "hltPFJetsL1FastL2L3[0].nhf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral HAD)", 40, 0, 1)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nch
    variable = ("hltPFJetsL1FastL2L3_nch1", "hltPFJetsL1FastL2L3[0].nch")
    binning = ("#scale[0.7]{HLT} Leading PFJet # charged constits", 40, 0, 40)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # mindphi
    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_3cj", "hltPFGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_j40", "hltPFGlobal.dijet_mindphi_j40")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{j40}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_mindphi_j40", "hltCaloGlobal.dijet_mindphi_j40")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{j40}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # dphijj
    variable = ("hltPFJetsL1FastL2L3_dphijj", "abs(deltaPhi(hltPFJetsL1FastL2L3[0].phi,hltPFJetsL1FastL2L3[1].phi))")
    binning = ("#scale[0.7]{HLT} PF #Delta#phi(j_{1},j_{2})", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))


    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_bench])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_bench])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig, benchmark=True)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["multijet"]:
    if plotting: del plotting[:]

    kColor = kCocoa2
    kTrig = "HT"

    sel_bench = "(Sum$(patJets.pt>50 && abs(patJets.eta)<2.5)>2 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[2].pt>50 && abs(patJets[2].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patJets[2].jetID==1 && patHTMHT.sumEt>450)"
    sel_bench = sel_noNoise + "*" + sel_bench

    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
    sel_trig1 = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 80, 250)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>-99))"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 34, 80, 250)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT} PFMETCleanUsingJetID [GeV]", 34, 80, 250)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltPFNoPUHT
    variable = ("hltPFNoPUHT", "hltPFHTMHTNoPU.sumEt")
    binning = ("#scale[0.7]{HLT} PFNoPUHT [GeV]", 28, 300, 1000)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>-99 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    # hltCaloHT
    variable = ("hltCaloHT", "hltCaloHTMHT.sumEt")
    binning = ("#scale[0.7]{HLT} CaloHT [GeV]", 28, 250, 950)
    addsel = "(hltCaloHTMHT.sumEt>-99 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    # hltCaloMHT
    variable = ("hltCaloMHT", "hltCaloHTMHT.pt")
    binning = ("#scale[0.7]{HLT} CaloMHT [GeV]", 27, 50, 320)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>-99 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>-99 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # l1MET
    variable = ("l1ETM", "l1MET.pt")
    binning = ("#scale[0.7]{L1T} ETM [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # l1HTT
    variable = ("l1HTT", "l1MHT.sumEt")
    binning = ("#scale[0.7]{L1T} HTT [GeV]", 30, 100, 400)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 40, 0, 400)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 40, 0, 200)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_eta1", "Sum$(hltPFJetsL1FastL2L3NoPU.eta * (hltPFJetsL1FastL2L3NoPU.pt==Max$(hltPFJetsL1FastL2L3NoPU.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 40, 0, 400)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nhf
    variable = ("hltPFJetsL1FastL2L3_nhf1", "hltPFJetsL1FastL2L3[0].nhf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral HAD)", 40, 0, 1)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nch
    variable = ("hltPFJetsL1FastL2L3_nch1", "hltPFJetsL1FastL2L3[0].nch")
    binning = ("#scale[0.7]{HLT} Leading PFJet # charged constits", 40, 0, 40)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # mindphi
    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_3cj", "hltPFGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # dphijj
    variable = ("hltPFJetsL1FastL2L3_dphijj", "abs(deltaPhi(hltPFJetsL1FastL2L3[0].phi,hltPFJetsL1FastL2L3[1].phi))")
    binning = ("#scale[0.7]{HLT} PF #Delta#phi(j_{1},j_{2})", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))


    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_bench])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_bench])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig, benchmark=True)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["bjet"]:
    if plotting: del plotting[:]

    kColor = kBoson
    kTrig = "btag"

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.bjet_maxcsv>0.898)"
    sel_bench = sel_noNoise + "*" + sel_bench

    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) )
    sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    #sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) )
    #sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 32, 60, 220)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 32, 60, 220)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT} PFMETCleanUsingJetID [GeV]", 32, 60, 220)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # maxcsv
    variable = ("hltCaloJetsL1Fast_maxcsv", "hltCaloGlobal.bjet_maxcsv")
    binning = ("#scale[0.7]{HLT} CaloJet max CSV", 22, 0, 1.1)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>-99 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_maxcsv2", "hltCaloGlobal.bjet_maxcsv2")
    binning = ("#scale[0.7]{HLT} CaloJet 2nd max CSV", 22, 0, 1.1)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nhf
    variable = ("hltPFJetsL1FastL2L3_nhf1", "hltPFJetsL1FastL2L3[0].nhf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral HAD)", 40, 0, 1)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # nch
    variable = ("hltPFJetsL1FastL2L3_nch1", "hltPFJetsL1FastL2L3[0].nch")
    binning = ("#scale[0.7]{HLT} Leading PFJet # charged constits", 40, 0, 40)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # mindphi
    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_3cj", "hltPFGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # dphijj
    variable = ("hltPFJetsL1FastL2L3_dphijj", "abs(deltaPhi(hltPFJetsL1FastL2L3[0].phi,hltPFJetsL1FastL2L3[1].phi))")
    binning = ("#scale[0.7]{HLT} PF #Delta#phi(j_{1},j_{2})", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))


    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_bench])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_bench])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig, benchmark=True)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["vbf"]:
    if plotting: del plotting[:]

    kColor = kQuark
    kTrig = "VBFAll"

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>50 && abs(patJets[0].eta)<4.7 && patJets[1].pt>50 && abs(patJets[1].eta)<4.7 && patJets[0].eta*patJets[1].eta<=0 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.vbf_maxmjj>900 && patGlobal.vbf_maxmjj_deta>3.5)"
    sel_bench = sel_noNoise + "*" + sel_bench

    #sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) )
    #sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_leadmjj>600 && hltPFGlobal.vbf_leadmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) )
    sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 50, 200)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 30, 50, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT} PFMETCleanUsingJetID [GeV]", 30, 50, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>-99)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_njets", "Sum$(hltPFJetsL1FastL2L3.pt>30 && abs(hltPFJetsL1FastL2L3.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_njets", "Sum$(hltPFJetsL1FastL2L3NoPU.pt>30 && abs(hltPFJetsL1FastL2L3NoPU.eta)<2.5)")
    binning = ("#scale[0.7]{HLT} # PFNoPUJets_{p_{T} > 30 GeV, |#eta| < 2.5}", 8, 0, 8)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>-99)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    # VBF jets
    variable = ("hltPFJetsL1FastL2L3_maxmjj", "hltPFGlobal.vbf_maxmjj")
    binning = ("#scale[0.7]{HLT} VBF max m_{jj} [GeV]", 34, 500, 2200)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>-99 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_leadmjj", "hltPFGlobal.vbf_leadmjj")
    binning = ("#scale[0.7]{HLT} VBF max m_{jj,leading} [GeV]", 34, 500, 2200)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_maxmjj_deta", "hltPFGlobal.vbf_maxmjj_deta")
    binning = ("#scale[0.7]{HLT} Leading-m_{jj} VBF #Delta#eta", 18, 2.5, 7)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>-99 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_leadmjj_deta", "hltPFGlobal.vbf_leadmjj_deta")
    binning = ("#scale[0.7]{HLT} Leading-m_{jj,leading} VBF #Delta#eta", 18, 2.5, 7)
    addsel = addsel
    plotting.append((variable, binning, addsel))

    #___________________________________________________________________________
    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltTrackMETDPhi
    variable = ("hltTrackMETDPhi", "abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))")
    binning = ("#scale[0.7]{HLT} #Delta#phi(PFMET,TrackMET)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_sumptf
    variable = ("hltPFMET_sumptchf", "hltTrackMET.sumEt/hltPFMET.sumEt")
    binning = ("#scale[0.7]{HLT} sum track p_{T}/sum E_{T}", 40, 0, 0.4)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # hltPFMET_signif
    variable = ("hltPFMET_signif", "hltPFMET.pt/sqrt(hltPFMET.sumEt)")
    binning = ("#scale[0.7]{HLT} PFMET/sqrt(sum E_{T})", 40, 0, 8)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # mindphi
    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_mindphi_3cj", "hltPFGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    # dphijj
    variable = ("hltPFJetsL1FastL2L3_dphijj", "abs(deltaPhi(hltPFJetsL1FastL2L3[0].phi,hltPFJetsL1FastL2L3[1].phi))")
    binning = ("#scale[0.7]{HLT} PF #Delta#phi(j_{1},j_{2})", 32, 0, 3.2)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))


    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_all"          , kBlack , kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_trig1"        , kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_noNoise_all"  , kGray+1, kWhite , variable[1], "*".join([sel, addsel, sel_bench])),
            (variable[0]+"_noNoise_trig1", kColor , kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_bench])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig, benchmark=True)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["puremet_eff"]:
    if plotting: del plotting[:]

    # use a slightly tighter offline noise cleaning
    sel_noNoise = sel_noNoise0
    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>%i && hltPFMET.pt>%i)"
    sel_trig2 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltPFMET.pt>%i)"
    sel_trig3 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltCaloMETCleanUsingJetID.pt>%i && hltPFMET.pt>%i)"
    sel_trig4 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltCaloMETCleanUsingJetID.pt>%i && hltPFMET.pt>%i && hltTrackMET.pt>%i)"

    addsel = "(1)"

    # recoPFMETT0T1: 200, 150, 120, 100, 80, 65
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 80, 440)
    thresholds = [0, 80, 90, 100, 110, 120]
    kHist = thresholds.index(90)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 200, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 80, 400)
    thresholds = [0, 65, 70, 80, 90, 100]
    kHist = thresholds.index(90)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 150, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 70, 400)
    thresholds = [0, 55, 65, 70, 80, 90]
    kHist = thresholds.index(80)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 120, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 360)
    thresholds = [0, 50, 55, 65, 70, 80]
    kHist = thresholds.index(70)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 100, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 40, 320)
    thresholds = [0, 50, 55, 60, 65, 70]
    kHist = thresholds.index(65)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 80, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 240)
    thresholds = [0, 45, 50, 55, 60, 65]
    kHist = thresholds.index(50)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 65, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    for p in plotting[:0]:
        (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig4 % (x, calomet1, calomet2, pfmet, trackmet), sel_noNoise])),
                ]
        kTrig = "PFMET%i" % pfmet
        trigs = ["HLT PFMET<%i" % pfmet]+[("HLT PFMET>%i && CaloMET>%i" % (pfmet,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # CaloMETClean or CaloMETCleanUsingJetID optimization
    if plotting: del plotting[:]

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 80, 440)
    thresholds = [0, 60, 70, 80, 90, 100]
    kHist = thresholds.index(70)
    calomet, calomet1, calomet2, pfmet, trackmet = 100, 0, 0, 220, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 80, 400)
    thresholds = [0, 50, 60, 70, 80, 90]
    kHist = thresholds.index(70)
    calomet, calomet1, calomet2, pfmet, trackmet = 90, 0, 0, 150, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 70, 400)
    thresholds = [0, 40, 50, 60, 70, 80]
    kHist = thresholds.index(60)
    calomet, calomet1, calomet2, pfmet, trackmet = 80, 0, 0, 120, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 360)
    thresholds = [0, 40, 45, 50, 60, 70]
    kHist = thresholds.index(60)
    calomet, calomet1, calomet2, pfmet, trackmet = 70, 0, 0, 100, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 40, 320)
    thresholds = [0, 30, 40, 45, 50, 60]
    kHist = thresholds.index(50)
    calomet, calomet1, calomet2, pfmet, trackmet = 65, 0, 0, 80, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 240)
    thresholds = [0, 30, 35, 40, 45, 50]
    kHist = thresholds.index(50)
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 0, 0, 65, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    for p in plotting:
        (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                #(variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, x, calomet2, pfmet, trackmet), sel_noNoise])),
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet1, x, pfmet, trackmet), sel_noNoise])),
                ]
        #kTrig = "CaloMET%i" % calomet
        kTrig = "CaloMET%i_2" % calomet
        trigs = ["HLT CaloMET<%i" % pfmet]+[("HLT CaloMET>%i && CaloMETClean>%i" % (calomet,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # Other stuff
    if plotting: del plotting[:]

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 0, 160)
    thresholds = [0, 65, 80, 100, 120, 150]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 0, 0
    kHist = thresholds.index(65)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 0, 160)
    thresholds = [0, 65, 80, 100, 120, 150]
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 0, 0, 0, 0
    kHist = thresholds.index(65)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 0, 160)
    thresholds = [0, 65, 80, 100, 120, 150]
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 0, 0, 0, 0
    kHist = thresholds.index(65)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    # hltTrackMET
    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 0, 160)
    thresholds = [0, 65, 80, 100, 120, 150]
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 0, 0, 0, 0
    kHist = thresholds.index(65)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    for p in plotting:
        (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p

        params = []
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet1, calomet2, x, trackmet), sel_noNoise])),
                ]
        kTrig = "PFMET%i" % pfmet
        trigs = [("HLT CaloMET>%i && PFMET>%i" % (calomet,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs, normalize=True)
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])


if sections["puremet_clean_eff"]:
    if plotting: del plotting[:]

    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>%i && hltPFMET.pt>%i)"
    sel_trig2 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltPFMET.pt>%i)"
    sel_trig3 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltCaloMETCleanUsingJetID.pt>%i && hltPFMET.pt>%i)"
    sel_trig4 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltCaloMETCleanUsingJetID.pt>%i && hltPFMET.pt>%i && hltTrackMET.pt>%i)"
    sel_trig5 = "(abs(deltaPhi(hltPFMET.phi,hltTrackMET.phi))<%.2f && hltTrackMET.pt>%i && hltPFMET.sumEt>0 && hltTrackMET.sumEt/hltPFMET.sumEt>%.2f)"

    addsel = "(1)"

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 80, 440)
    thresholds = [0, 10, 15, 20, 30, 40]
    kHist = thresholds.index(20)
    calomet, calomet1, calomet2, pfmet, trackmet = 100, 90, 90, 220, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 80, 400)
    thresholds = [0, 10, 15, 20, 25, 30]
    kHist = thresholds.index(15)
    calomet, calomet1, calomet2, pfmet, trackmet = 90, 80, 80, 150, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 70, 400)
    thresholds = [0, 10, 12, 15, 20, 25]
    kHist = thresholds.index(15)
    calomet, calomet1, calomet2, pfmet, trackmet = 80, 70, 70, 120, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 360)
    thresholds = [0, 5, 8, 10, 12, 15]
    kHist = thresholds.index(10)
    calomet, calomet1, calomet2, pfmet, trackmet = 70, 60, 60, 100, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 40, 320)
    thresholds = [0, 5, 8, 10, 12, 15]
    kHist = thresholds.index(10)
    calomet, calomet1, calomet2, pfmet, trackmet = 65, 55, 55, 80, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 240)
    thresholds = [0, 5, 8, 10, 12, 15]
    kHist = thresholds.index(10)
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 40, 40, 65, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    for p in plotting[:0]:
        (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise0])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet1, calomet2, pfmet, x), sel_noNoise0])),
                ]
        kTrig = "PFMET%i_2" % pfmet
        trigs = ["HLT (PFMET<%i || CaloMET<%i)" % (pfmet, calomet)]+[("HLT PFMET>%i && CaloMET>%i && TrackMET>%i" % (pfmet,calomet,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs, legend=(0.23,0.625,0.96,0.94))
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs, legend=(0.24,0.695,0.96,0.94))
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # TrackMET only
    if plotting: del plotting[:]

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 320)
    thresholds = [0, 80, 90, 100, 110, 120]
    kHist = thresholds.index(90)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 0, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    # TrackMET as x variable
    variable = ("patMPT", "patMPT.pt")
    binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 30, 320)
    thresholds = [0, 80, 90, 100, 110, 120]
    kHist = thresholds.index(90)
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 0, 0
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    for p in plotting[:0]:
        (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise0])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet1, calomet2, pfmet, x), sel_noNoise0])),
                ]
        kTrig = "TrackMET%i" % trackmet
        trigs = ["HLT all"]+[("HLT TrackMET>%i" % (x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs, legend=(0.23,0.625,0.96,0.94))
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs, legend=(0.24,0.695,0.96,0.94))
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # TrackMETdphi, TrackMET, sumpt TrackMET
    if plotting: del plotting[:]

    addsel = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150 && hltPFMET.sumEt>0)"

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 30, 320)
    thresholds = [99, 1.00, 0.70]
    kHist = thresholds.index(0.70)
    trackmetdphi, trackmet, sumptchf = 99., -99, -99.
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise0])),
        ]
    for i, x in enumerate(thresholds):
        trackmetdphi = x
        params += [
            (variable[0]+("_trig%i" % i), tones[1+i*2], tones[1+i*2], variable[1], "*".join([sel, addsel, sel_trig5 % (trackmetdphi, trackmet, sumptchf), sel_noNoise0])),
            ]
    kTrig = "trackmetdphi"
    trigs = ["HLT TrackMET #Delta#phi<%.2f" % thresholds[0]]+[("HLT TrackMET #Delta#phi<%.2f" % (x)) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])

    thresholds = [-99, 10, 20]
    kHist = thresholds.index(10)
    trackmetdphi, trackmet, sumptchf = 99., -99, -99.
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise0])),
        ]
    for i, x in enumerate(thresholds):
        trackmet = x
        params += [
            (variable[0]+("_trig%i" % i), tones[1+i*2], tones[1+i*2], variable[1], "*".join([sel, addsel, sel_trig5 % (trackmetdphi, trackmet, sumptchf), sel_noNoise0])),
            ]
    kTrig = "trackmet"
    trigs = ["HLT TrackMET>%i" % thresholds[0]]+[("HLT TrackMET>%i" % (x)) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    thresholds = [-99, 0.02, 0.04]
    kHist = thresholds.index(0.02)
    trackmetdphi, trackmet, sumptchf = 99., -99, -99.

    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise0])),
        ]
    for i, x in enumerate(thresholds):
        sumptchf = x
        params += [
            (variable[0]+("_trig%i" % i), tones[1+i*2], tones[1+i*2], variable[1], "*".join([sel, addsel, sel_trig5 % (trackmetdphi, trackmet, sumptchf), sel_noNoise0])),
            ]
    kTrig = "sumptchf"
    trigs = ["HLT sumptchf>%.2f" % thresholds[0]]+[("HLT sumptchf>%.2f" % (x)) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()



if sections["monojet_eff"]:
    if plotting: del plotting[:]

    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
    sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>%i)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>%i)>0 && hltPFMETNoMu.pt>105)"
    sel_trig2 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>%i)>0 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>%i)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>%i)>0)"
    sel_trig3 = "(hltPFJetsL1FastL2L3[0].nhf<%.2f && ((abs(hltPFJetsL1FastL2L3[0].eta)<2.4 && hltPFJetsL1FastL2L3[0].nch>%i && hltPFJetsL1FastL2L3[0].chf>%.3f) || abs(hltPFJetsL1FastL2L3[0].eta)>=2.4))"  # use |eta|<2.4

    addsel = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && hltPFMET.pt>100)"

    # ptj1
    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 320)
    thresholds = [0, 90, 100, 105, 110, 120]
    kHist = thresholds.index(100)
    calojet, pfjet, pfnopujet = 0, 120, 0
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 280)
    thresholds = [0, 70, 80, 85, 90, 100]
    kHist = thresholds.index(80)
    calojet, pfjet, pfnopujet = 0, 100, 0
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 240)
    thresholds = [0, 50, 60, 65, 70, 80]
    kHist = thresholds.index(60)
    calojet, pfjet, pfnopujet = 0, 80, 0
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 40, 50, 55, 60, 70]
    kHist = thresholds.index(50)
    calojet, pfjet, pfnopujet = 0, 70, 0
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 30, 40, 45, 50, 60]
    kHist = thresholds.index(40)
    calojet, pfjet, pfnopujet = 0, 60, 0
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 20, 30, 35, 40, 50]
    kHist = thresholds.index(30)
    calojet, pfjet, pfnopujet = 0, 50, 0
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    for p in plotting[:0]:
        (variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise1])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig2 % (x, pfjet, pfnopujet), sel_noNoise1])),
                ]
        kTrig = "PFJet%i" % pfjet
        trigs = ["HLT PFJet1 p_{T}<%i" % pfjet]+[("HLT PFJet1 p_{T}>%i && CaloJet1 p_{T}>%i" % (pfjet,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs, legend=(0.48,0.695,0.94,0.94))
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])

    #___________________________________________________________________________
    # PFnoPUJet
    if plotting: del plotting[:]

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 320)
    thresholds = [0, 90, 100, 105, 110, 120]
    kHist = thresholds.index(100)
    calojet, pfjet, pfnopujet = 0, 0, 120
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 280)
    thresholds = [0, 70, 80, 85, 90, 100]
    kHist = thresholds.index(80)
    calojet, pfjet, pfnopujet = 0, 0, 100
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 240)
    thresholds = [0, 50, 60, 65, 70, 80]
    kHist = thresholds.index(60)
    calojet, pfjet, pfnopujet = 0, 0, 80
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 40, 50, 55, 60, 70]
    kHist = thresholds.index(50)
    calojet, pfjet, pfnopujet = 0, 0, 70
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 30, 40, 45, 50, 60]
    kHist = thresholds.index(40)
    calojet, pfjet, pfnopujet = 0, 0, 60
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 20, 30, 35, 40, 50]
    kHist = thresholds.index(30)
    calojet, pfjet, pfnopujet = 0, 0, 50
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    for p in plotting[:0]:
        (variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise1])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig2 % (x, pfjet, pfnopujet), sel_noNoise1])),
                ]
        kTrig = "PFNoPUJet%i" % pfnopujet
        trigs = ["HLT PFNoPUJet1 p_{T}<%i" % pfnopujet]+[("HLT PFNoPUJet1 p_{T}>%i && CaloJet1 p_{T}>%i" % (pfnopujet,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs, legend=(0.48,0.695,0.94,0.94))
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # nhf & nch & chf
    if plotting: del plotting[:]

    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>0)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100)"
    sel_bench_met = "((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && abs(patJets[0].eta)<2.5 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && abs(deltaPhi(patJets[0].phi,patJets[1].phi))<2.5) || (Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1)) * (recoPFMETT0T1.pt>100)"
    sel_bench_met = sel_noNoise + "*" + sel_bench_met

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 50, 360)
    thresholds = [99., 0.95, 0.90]
    kHist = thresholds.index(0.95)
    nhf, nch, chf = 99., -99, -99.
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench_met])),
        ]
    for i, x in enumerate(thresholds):
        nhf = x
        params += [
            (variable[0]+("_trig%i" % i), tones[1+i*2], tones[1+i*2], variable[1], "*".join([sel, addsel, sel_trig3 % (nhf, nch, chf), sel_bench_met])),
            ]
    kTrig = "nhf"
    trigs = ["HLT PFJet1 nhf<%.2f" % thresholds[0]]+[("HLT PFJet1 nhf<%.2f" % x) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94), ymin=0.8, ymax=1.1)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()

    thresholds = [-99, 0, 1]
    kHist = thresholds.index(0)
    nhf, nch, chf = 0.95, -99, -99.
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench_met])),
        ]
    for i, x in enumerate(thresholds):
        nch = x
        params += [
            (variable[0]+("_trig%i" % i), tones[1+i*2], tones[1+i*2], variable[1], "*".join([sel, addsel, sel_trig3 % (nhf, nch, chf), sel_bench_met])),
            ]
    kTrig = "nch"
    trigs = ["HLT PFJet1 nhf<%.2f" % nhf]+[("HLT PFJet1 nhf<%.2f && nch>%i" % (nhf,x)) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94), ymin=0.8, ymax=1.1)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()

    thresholds = [-99., 0.001, 0.01]
    kHist = thresholds.index(0.01)
    nhf, nch, chf = 0.95, -99, -99.
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench_met])),
        ]
    for i, x in enumerate(thresholds):
        chf = x
        params += [
            (variable[0]+("_trig%i" % i), tones[1+i*2], tones[1+i*2], variable[1], "*".join([sel, addsel, sel_trig3 % (nhf, nch, chf), sel_bench_met])),
            ]
    kTrig = "chf"
    trigs = ["HLT PFJet1 nhf<%.2f" % nhf]+[("HLT PFJet1 nhf<%.2f && chf>%.3f" % (nhf,x)) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94), ymin=0.8, ymax=1.1)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()


if sections["higdijet_eff"]:
    if plotting: del plotting[:]

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>60 && abs(patJets[0].eta)<2.5 && patJets[1].pt>30 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_maxpt>110 && patGlobal.dijet_maxpt_mjj<250 && patGlobal.dijet_mindphi_3cj>0.5)"
    sel_bench = sel_noNoise + "*" + sel_bench
    sel_trig1 = "(hltPFGlobal.dijet_mindphi_2cj>%.1f && hltPFGlobal.dijet_mindphi_3cj>%.1f)"

    addsel = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && hltPFMET.pt>100)"

    # recoPFMETT0T1
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 360)
    thresholds = [-99., 0.2, 0.3, 0.4, 0.5, 0.7]
    kHist = thresholds.index(0.5)
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench])),
        ]
    for i, x in enumerate(thresholds):
        params += [
            (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig1 % (x, -99.), sel_bench])),
            ]
    kTrig = "mindphi_2cj"
    trigs = ["HLT PF min#Delta#phi_{2cj}>%.1f" % thresholds[0]]+[("HLT PF min#Delta#phi_{2cj}>%.1f" % x) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()

    thresholds = [-99., 0.2, 0.3, 0.4, 0.5, 0.7]
    kHist = thresholds.index(0.5)
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench])),
        ]
    for i, x in enumerate(thresholds):
        params += [
            (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig1 % (-99., x), sel_bench])),
            ]
    kTrig = "mindphi_3cj"
    trigs = ["HLT PF min#Delta#phi_{3cj}>%.1f" % thresholds[0]]+[("HLT PF min#Delta#phi_{3cj}>%.1f" % x) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()


if sections["susdijet_eff"]:
    if plotting: del plotting[:]

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>70 && abs(patJets[0].eta)<2.5 && patJets[1].pt>70 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_mindphi_2cj>0.5 && ((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>2 && abs(patJets[2].eta)<2.5 && patJets[2].jetID==1 && patGlobal.dijet_mindphi_3cj>0.3)||Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==2))"
    sel_bench = sel_noNoise + "*" + sel_bench
    sel_trig1 = "(hltPFGlobal.dijet_mindphi_2cj>%.1f && hltPFGlobal.dijet_mindphi_3cj>%.1f)"

    addsel = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1 && hltPFMET.pt>80)"

    # recoPFMETT0T1
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 360)
    thresholds = [-99., 0.2, 0.3, 0.4, 0.5, 0.7]
    kHist = thresholds.index(0.5)
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench])),
        ]
    for i, x in enumerate(thresholds):
        params += [
            (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig1 % (x, -99.), sel_bench])),
            ]
    kTrig = "mindphi_2cj_2"
    trigs = ["HLT PF min#Delta#phi_{2cj}>%.1f" % thresholds[0]]+[("HLT PF min#Delta#phi_{2cj}>%.1f" % x) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()

    thresholds = [-99., 0.2, 0.3, 0.4, 0.5, 0.7]
    kHist = thresholds.index(0.5)
    params = [
        (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_bench])),
        ]
    for i, x in enumerate(thresholds):
        params += [
            (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig1 % (-99., x), sel_bench])),
            ]
    kTrig = "mindphi_3cj_2"
    trigs = ["HLT PF min#Delta#phi_{3cj}>%.1f" % thresholds[0]]+[("HLT PF min#Delta#phi_{3cj}>%.1f" % x) for x in thresholds]

    histos = book_ratio(params, binning)
    project(params, histos, drawOverflow=False)
    legs = draw_effnum(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    legs = draw_eff(params, histos, kHist, trigs)
    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])
    print histos[0].GetEntries()


if sections["susdijet_eff2"]:
    if plotting: del plotting[:]

    def label_var(histos, signal="PFMET80+CJ60x2", legend=(0.52,0.84,0.96,0.94)):
        leg1 = TLegend(*legend)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "QCD", "f")
        leg1.AddEntry(histos[3], signal, "f")
        leg1.Draw()
        return (leg1)  # persistent


    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
    sel_trig1 = "(hltCaloGlobal.dijet_mindphi_j30>%.2f && hltCaloGlobal.dijet_mindphi_j40>%.2f && hltCaloGlobal.dijet_mindphi_cj30>%.2f && hltCaloGlobal.dijet_mindphi_cj40>%.2f && hltCaloGlobal.dijet_mindphi_2cj>%.2f && hltCaloGlobal.dijet_mindphi_3cj>%.2f)"
    sel_trig2 = "(hltPFGlobal.dijet_mindphi_j30>%.2f && hltPFGlobal.dijet_mindphi_j40>%.2f && hltPFGlobal.dijet_mindphi_cj30>%.2f && hltPFGlobal.dijet_mindphi_cj40>%.2f && hltPFGlobal.dijet_mindphi_2cj>%.2f && hltPFGlobal.dijet_mindphi_3cj>%.2f)"

    addsel = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>50 && hltPFMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>40)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1)"
    kColor = tones[2]

    # offline mindphi
    variable = ("patmindphi_j30", "patGlobal.dijet_mindphi_j30")
    binning = ("#scale[0.7]{RECO} min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("patmindphi_2cj", "patGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{RECO} min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    # offline alphaT
    variable = ("patalphat", "alphaT(patJets[0].pt, patJets[0].px, patJets[0].py, patJets[1].pt, patJets[1].px, patJets[1].py)")
    binning = ("#scale[0.7]{RECO} #alpha_{T}", 32, 0, 3.2)
    plotting.append((variable, binning))

    # all mindphi
    variable = ("hltCaloJetsL1Fast_mindphi_j30", "hltCaloGlobal.dijet_mindphi_j30")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltCaloJetsL1Fast_mindphi_j40", "hltCaloGlobal.dijet_mindphi_j40")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{j40}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltCaloJetsL1Fast_mindphi_cj30", "hltCaloGlobal.dijet_mindphi_cj30")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{cj30}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltCaloJetsL1Fast_mindphi_cj40", "hltCaloGlobal.dijet_mindphi_cj40")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{cj40}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltCaloJetsL1Fast_mindphi_2cj", "hltCaloGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltCaloJetsL1Fast_mindphi_3cj", "hltCaloGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} Calo min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    # (PF)
    variable = ("hltPFJetsL1FastL2L3_mindphi_j30", "hltPFGlobal.dijet_mindphi_j30")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{j30}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltPFJetsL1FastL2L3_mindphi_j40", "hltPFGlobal.dijet_mindphi_j40")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{j40}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltPFJetsL1FastL2L3_mindphi_cj30", "hltPFGlobal.dijet_mindphi_cj30")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{cj30}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltPFJetsL1FastL2L3_mindphi_cj40", "hltPFGlobal.dijet_mindphi_cj40")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{cj40}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltPFJetsL1FastL2L3_mindphi_2cj", "hltPFGlobal.dijet_mindphi_2cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{2cj}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))

    variable = ("hltPFJetsL1FastL2L3_mindphi_3cj", "hltPFGlobal.dijet_mindphi_3cj")
    binning = ("#scale[0.7]{HLT} PF min #Delta#phi_{3cj}(MET,jet)", 32, 0, 3.2)
    plotting.append((variable, binning))


    for p in plotting:
        (variable, binning) = p

        kVar = "mindphi_j30"
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise2])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise2])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD1])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD1])),
            ]
        #histos = book(params, binning)
        #project(params, histos)
        #draw_rate(params, histos)
        #legs = label_var(histos)
        #save(imgdir, "diagnosis_effvar_" + kVar + "_" + variable[0])

        kVar = "mindphi_2cj"
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise2])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise2])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD])),
            ]
        #histos = book(params, binning)
        #project(params, histos)
        #draw_rate(params, histos)
        #legs = label_var(histos)
        #save(imgdir, "diagnosis_effvar_" + kVar + "_" + variable[0])

        kVar = "alphaT"
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise2])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_noNoise2])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD2])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD2])),
            ]
        #histos = book(params, binning)
        #project(params, histos)
        #draw_rate(params, histos)
        #legs = label_var(histos)
        #save(imgdir, "diagnosis_effvar_" + kVar + "_" + variable[0])



    #___________________________________________________________________________
    # recoPFMETT0T1
    if plotting: del plotting[:]

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 40, 320)
    thresholds = [-99., 0.1, 0.2, 0.3, 0.4, 0.5]
    kHist = thresholds.index(0.3)

    kTrig = "hltCaloJetsL1Fast_mindphi_j30"
    trig = "HLT Calo min #Delta#phi_{j30}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_j40"
    trig = "HLT Calo min #Delta#phi_{j40}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_cj30"
    trig = "HLT Calo min #Delta#phi_{cj30}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_cj40"
    trig = "HLT Calo min #Delta#phi_{cj40}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_2cj"
    trig = "HLT Calo min #Delta#phi_{2cj}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_3cj"
    trig = "HLT Calo min #Delta#phi_{3cj}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_j30"
    trig = "HLT PF min #Delta#phi_{j30}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_j40"
    trig = "HLT PF min #Delta#phi_{j40}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_cj30"
    trig = "HLT PF min #Delta#phi_{cj30}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_cj40"
    trig = "HLT PF min #Delta#phi_{cj40}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_2cj"
    trig = "HLT PF min #Delta#phi_{2cj}>%.2f"
    plotting.append((kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_3cj"
    trig = "HLT PF min #Delta#phi_{3cj}>%.2f"
    plotting.append((kTrig, trig))


    for ip, p in enumerate(plotting):
        (kTrig, trig) = p
        #kTrig += "_2"

        if   ip<6   :  sel_trig = sel_trig1
        else        :  sel_trig = sel_trig2

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD1])),
            #(variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD])),
            ]
        for i, x in enumerate(thresholds):
            dphi1, dphi2, dphi3, dphi4, dphi5, dphi6 = -99., -99., -99., -99., -99., -99.
            if   ip%6==0:  dphi1 = x
            elif ip%6==1:  dphi2 = x
            elif ip%6==2:  dphi3 = x
            elif ip%6==3:  dphi4 = x
            elif ip%6==4:  dphi5 = x
            elif ip%6==5:  dphi6 = x

            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig % (dphi1, dphi2, dphi3, dphi4, dphi5, dphi6), sel_noNoise2, sel_noQCD1])),
                #(variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig % (dphi1, dphi2, dphi3, dphi4, dphi5, dphi6), sel_noNoise2, sel_noQCD])),
                ]
            trigs = [trig % -99.]+[trig % x for x in thresholds]

        #histos = book_ratio(params, binning)
        #project(params, histos, drawOverflow=False)
        #legs = draw_effnum(params, histos, kHist, trigs)
        #save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        #legs = draw_eff(params, histos, kHist, trigs)
        #save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # recoPFMETT0T1
    if plotting: del plotting[:]

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 40, 320)
    thresholds = [-99., 0.1, 0.2, 0.3, 0.4, 0.5]
    kHist = thresholds.index(0.3)
    dphi1, dphi2, dphi3, dphi4, dphi5, dphi6 = -99., -99., -99., -99., -99., -99.

    kTrig = "hltCaloJetsL1Fast_mindphi_2cj"
    trig = "HLT Calo min #Delta#phi_{2cj}>%.2f"
    addsel = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>50 && hltPFMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>40)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1)"
    plotting.append((addsel, kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_2cj"
    trig = "HLT PF min #Delta#phi_{2cj}>%.2f"
    addsel = addsel
    plotting.append((addsel, kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_2cj_CJ100CJ30"
    trig = "HLT Calo min #Delta#phi_{2cj}>%.2f"
    addsel = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>50 && hltPFMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>80)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>1)"
    plotting.append((addsel, kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_2cj_CJ100CJ30"
    trig = "HLT PF min #Delta#phi_{2cj}>%.2f"
    addsel = addsel
    plotting.append((addsel, kTrig, trig))

    kTrig = "hltCaloJetsL1Fast_mindphi_2cj_CJ60CJ30"
    trig = "HLT Calo min #Delta#phi_{2cj}>%.2f"
    addsel = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && hltPFMET.pt>100 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>40)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1)"
    plotting.append((addsel, kTrig, trig))

    kTrig = "hltPFJetsL1FastL2L3_mindphi_2cj_CJ60CJ30"
    trig = "HLT PF min #Delta#phi_{2cj}>%.2f"
    addsel = addsel
    plotting.append((addsel, kTrig, trig))


    for ip, p in enumerate(plotting):
        (addsel, kTrig, trig) = p
        #kTrig += "_2"

        if ip==0 or ip==2 or ip==4:  sel_trig = sel_trig1
        else                      :  sel_trig = sel_trig2

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD1])),
            #(variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise2, sel_noQCD])),
            ]
        for i, x in enumerate(thresholds):
            dphi5 = x
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig % (dphi1, dphi2, dphi3, dphi4, dphi5, dphi6), sel_noNoise2, sel_noQCD1])),
                #(variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig % (dphi1, dphi2, dphi3, dphi4, dphi5, dphi6), sel_noNoise2, sel_noQCD])),
                ]
            trigs = [trig % -99.]+[trig % x for x in thresholds]

        #histos = book_ratio(params, binning)
        #project(params, histos, drawOverflow=False)
        #legs = draw_effnum(params, histos, kHist, trigs)
        #save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        #legs = draw_eff(params, histos, kHist, trigs)
        #save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


if sections["multijet_eff"]:
    if plotting: del plotting[:]

    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
    sel_trig1 = "(hltCaloHTMHT.sumEt>%i && hltPFHTMHTNoPU.sumEt>%i)"
    addsel = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && hltPFMET.pt>100)"

    # HT
    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 150, 700)
    thresholds = [0, 200, 225, 250, 275, 300]
    kHist = thresholds.index(250)
    calo, pf = 0, 300
    plotting.append((variable, binning, thresholds, calo, pf, kHist))

    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 200, 750)
    thresholds = [0, 250, 275, 300, 325, 350]
    kHist = thresholds.index(300)
    calo, pf = 0, 350
    plotting.append((variable, binning, thresholds, calo, pf, kHist))

    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 250, 800)
    thresholds = [0, 300, 325, 350, 375, 400]
    kHist = thresholds.index(350)
    calo, pf = 0, 400
    plotting.append((variable, binning, thresholds, calo, pf, kHist))

    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 300, 850)
    thresholds = [0, 300, 350, 375, 400, 450]
    kHist = thresholds.index(350)
    calo, pf = 0, 450
    plotting.append((variable, binning, thresholds, calo, pf, kHist))

    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 450, 950)
    thresholds = [0, 350, 400, 450, 500, 550]
    kHist = thresholds.index(450)
    calo, pf = 0, 550
    plotting.append((variable, binning, thresholds, calo, pf, kHist))

    variable = ("patHT", "patHTMHT.sumEt")
    binning = ("#scale[0.7]{RECO} PFHT [GeV]", 500, 1000)
    thresholds = [0, 450, 500, 550, 600, 650]
    kHist = thresholds.index(600)
    calo, pf = 0, 650
    plotting.append((variable, binning, thresholds, calo, pf, kHist))


    for p in plotting:
        (variable, binning, thresholds, calo, pf, kHist) = p

        params = [
            (variable[0]+"_filt", toneC, toneC, variable[1], "*".join([sel, addsel, sel_noNoise2])),
            ]
        for i, x in enumerate(thresholds):
            params += [
                (variable[0]+("_trig%i" % i), tones[i], tones[i], variable[1], "*".join([sel, addsel, sel_trig1 % (x, pf), sel_noNoise2])),
                ]
        kTrig = "PFNoPUHT%i" % pf
        trigs = ["HLT PFNoPUHT<%i" % pf]+[("HLT PFNoPUHT>%i && CaloHT>%i" % (pf,x)) for x in thresholds]

        histos = book_ratio(params, binning)
        project(params, histos, drawOverflow=False)
        legs = draw_effnum(params, histos, kHist, trigs, legend=(0.48,0.695,0.94,0.94))
        save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
        legs = draw_eff(params, histos, kHist, trigs)
        save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


if sections["future_topology_hlt"]:
    if plotting: del plotting[:]

    pass


if sections["future_triggers"]:
    if plotting: del plotting[:]

    blinds = [
        TColor.GetColor("#332288"),
        TColor.GetColor("#44AA99"),
        TColor.GetColor("#117733"),
        TColor.GetColor("#999933"),
        TColor.GetColor("#882255"),
        TColor.GetColor("#AA4499"),
        ]

    def draw_future_A(variable, binning, selections, entries, sel_sig, addsel="(1)", benchmark=False):
        # Set colors
        if len(selections) == 2:
            colors = [blinds[0], blinds[4]]
        elif len(selections) == 3:
            colors = [blinds[0], blinds[2], blinds[4]]
        elif len(selections) == 4:
            colors = [blinds[0], blinds[1], blinds[2], blinds[4]]
        else:
            colors = blinds[:len(selections)]

        # Book/Make histograms
        params = []
        for i, s in enumerate(selections):
            params += [
                (variable[0]+"_trig%i" % i, colors[i], colors[i], variable[1], "*".join([sel, addsel, s])),
                ]
        for i, s in enumerate(selections):
            params += [
                (variable[0]+"_trig%i" % (i+len(selections)), colors[i], colors[i], variable[1], "*".join([sel, addsel, s, sel_sig])),
                ]
        histos = book(params, binning)
        project(params, histos)

        # Draw histograms
        ytitle = "Events"
        logy = False
        legend = (0.18,0.94-0.04*len(selections),0.96,0.94)

        ymax = max(histos[0].GetMaximum(), histos[1].GetMaximum())
        histos[0].SetMaximum(ymax * 1.25)
        histos[0].SetMinimum(0)
        histos[0].GetYaxis().SetTitle(ytitle)

        histos1 = []
        for h in histos[:len(histos)/2]:
            hclone = h.Clone(h.GetName()+"_clone")
            h.SetFillStyle(3004)
            hclone.SetFillStyle(0)
            histos1.append(hclone)

        for h in histos[len(histos)/2:]:
            hclone = h.Clone(h.GetName()+"_clone")
            hclone.SetLineColor(kGray+1)
            hclone.SetFillStyle(0)
            hclone.SetLineWidth(1)
            histos1.append(hclone)

        histos[0].Draw("hist")
        for h in histos:
            h.Draw("hist same")
        for h in histos1:
            h.Draw("hist same")

        assert(len(histos)/2 == len(entries))
        leg1 = TLegend(legend[0], legend[1], legend[0]+0.6*(legend[2]-legend[0]), legend[3])
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        for h, e in zip(histos[:len(histos)/2], entries):
            leg1.AddEntry(h, e+"  %.0f [%.0f%%]" % (h.Integral(), h.Integral() / histos[0].Integral() * 100), "f")
        leg1.Draw()

        leg2 = TLegend(legend[0]+0.6*(legend[2]-legend[0]), legend[1], legend[2], legend[3])
        leg2.SetFillStyle(0); leg2.SetLineColor(0); leg2.SetShadowColor(0); leg2.SetBorderSize(0)
        filt = "(no noise)" if not benchmark else "(benchmark)"
        for h, e in zip(histos[len(histos)/2:], entries):
            leg2.AddEntry(h, filt+"  %.0f [%.0f%%]" % (h.Integral(), h.Integral() / histos[len(histos)/2].Integral() * 100), "f")
        leg2.Draw()

        gPad.SetLogy(logy)
        CMS_label()
        return (histos, histos1, leg1, leg2)  # persistent


    def draw_future_B(selections, entries, numerator=False):
        colors = [blinds[0], blinds[4]]
        variable = ("nGoodPV", "event.nGoodPV")
        binning = ("#scale[0.7]{RECO} # good PV", 2, 40)

        # Book/make histograms
        assert(len(selections)==2)
        params = [(variable[0]+"_all", kBlack, kBlack, variable[1], "*".join([sel, addsel])),]
        for i, s in enumerate(selections):
            params += [
                (variable[0]+"_trig%i" % (i), colors[i], colors[i], variable[1], "*".join([sel, addsel, s])),
                ]
        histos = []
        xlow, xup = binning[1], binning[2]
        tmpbins = [2,4,6,8]+range(10,60+2,2)
        tmpbins2 = [xlow] + [x for x in tmpbins if (x > xlow and x < xup)] + [xup]
        for i, p in enumerate(params):
            h = TH1F("h_"+p[0], "; "+binning[0], len(tmpbins2)-1, numpy.array(tmpbins2, dtype=float))
            h.SetLineColor(p[1])
            h.SetMarkerColor(p[1])
            histos.append(h)
        project(params, histos, drawOverflow=False)

        # Divide
        if not numerator:
            for h in histos[1:]:
                h.Divide(h, histos[0], 1, 1, "b")

        ## Normalize
        #for h in histos:
        #    h.SetBinContent(1, h.GetBinContent(1)/3)
        #    h.SetBinError(1, h.GetBinError(1)/3)
        #    h.Scale(1.0/h.GetBinContent(1))

        # Draw histograms
        ytitle = "arbitrary unit" if not numerator else "Events"
        logy = False
        legend = (0.50,0.94-0.05*len(selections),0.96,0.94)

        ymax = max(histos[1].GetMaximum(), histos[2].GetMaximum())
        histos[1].SetMaximum(ymax * 1.25)
        histos[1].SetMinimum(0)
        histos[1].GetYaxis().SetTitle(ytitle)

        histos[1].Draw("p")
        for h in histos[1:]:
            h.Draw("p same")
        if numerator:
            scale = int(histos[0].GetSumOfWeights()/histos[1].GetSumOfWeights())
            for i in xrange(10):
                if scale / pow(10,i) <10:
                    scale = scale / pow(10,i) * pow(10,i)
                    break
            histos[0].Scale(1.0/scale)
        histos[0].Draw("hist same")

        leg1 = TLegend(*legend)
        leg1.SetFillStyle(0); leg1.SetLineColor(0); leg1.SetShadowColor(0); leg1.SetBorderSize(0)
        if numerator:
            leg1.SetY1NDC(leg1.GetY1NDC()-0.05)
            leg1.AddEntry(histos[0], "all (x1/%i)" % scale, "l")
        for h, e in zip(histos[1:], entries):
            leg1.AddEntry(h, e, "lp")
        leg1.Draw()

        gPad.SetLogy(logy)
        CMS_label()
        #latex.DrawLatex(0.58, 0.79, "(normalized to 1st bin in 2012D)")
        return (histos, leg1)  # persistent


    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} T0T1 PFMET [GeV]", 50, 0, 250)
    addsel = "(1)"

    #___________________________________________________________________________
    # Inclusive
    kColor = kBlue2
    kTrig = "PFMET150"
    selections = [
        "(hltCaloMET.pt>80 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>-99 && hltCaloMETCleanUsingJetID.pt>-99 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>-99 && hltPFMET.pt>150)",
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)",
        ]
    entries = [
        "2012D",
        "2015#alpha (CaloMET>90)",
        "2015#alpha (HBHE only)",
        "2015#alpha (HBHE && JetID)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_noNoise0, benchmark=False)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # MonoCentralJet
    kColor = kYellow3
    sel_bench = "((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && abs(deltaPhi(patJets[0].phi,patJets[1].phi))<2.5) || (Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1))"
    sel_bench = sel_noNoise + "*" + sel_bench

    kTrig = "PFMET100CJ100"
    selections = [
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)",
        ]
    entries = [
        "2012D",
        "2015#alpha (HBHE)",
        "2015#alpha (HBHE && CJ100)",
        "2015#alpha (HBHE && CJ100 && nch)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])

    kTrig = "PFMET120CJ80"
    selections = [
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>105||hltPFMET.pt>105) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>80 && hltCaloMETClean.pt>70 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>120||hltPFMET.pt>120) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
        "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>80 && hltCaloMETClean.pt>70 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>120||hltPFMET.pt>120) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)",
        ]
    entries = [
        "2012D",
        "2015#alpha (HBHE)",
        "2015#alpha (HBHE && PFMET120)",
        "2015#alpha (HBHE && PFMET120 && nch)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # DiCentralJet
    kColor = kOrange3
    kTrig = "PFMET100CJ60CJ30"
    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>60 && abs(patJets[0].eta)<2.5 && patJets[1].pt>30 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_maxpt>110 && patGlobal.dijet_maxpt_mjj<250 && patGlobal.dijet_mindphi_3cj>0.5)"
    sel_bench = sel_noNoise + "*" + sel_bench

    selections = [
        "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)",
        "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFGlobal.dijet_mindphi_2cj>-99)",
        "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFGlobal.dijet_mindphi_2cj>-99)",
        "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFGlobal.dijet_mindphi_2cj>0.5)",
        ]
    entries = [
        "2012D",
        "2015#alpha (relaxed)",
        "2015#alpha (relaxed && CJ60,30)",
        "2015#alpha (relaxed && CJ60,30 && mdphi)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])


    # DiCentralJetSUS
    kColor = kOrange3
    kTrig = "PFMET80CJ60x2"
    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>70 && abs(patJets[0].eta)<2.5 && patJets[1].pt>70 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.dijet_mindphi_2cj>0.5 && ((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>2 && abs(patJets[2].eta)<2.5 && patJets[2].jetID==1 && patGlobal.dijet_mindphi_3cj>0.3)||Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==2))"
    sel_bench = sel_noNoise + "*" + sel_bench

    selections = [
        "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )",
        "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) && hltPFGlobal.dijet_mindphi_2cj>-99.)",
        "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) && hltPFGlobal.dijet_mindphi_2cj>-99.)",
        "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) && hltPFGlobal.dijet_mindphi_2cj>0.3)",
        ]
    entries = [
        "2012D",
        "2015#alpha (CaloMET>65)",
        "2015#alpha (CaloMET>65 && CJ60x2)",
        "2015#alpha (CaloMET>65 && CJ60x2 && mdphi)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # HT

    kColor = kCocoa2
    sel_bench = "(Sum$(patJets.pt>50 && abs(patJets.eta)<2.5)>2 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[2].pt>50 && abs(patJets[2].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patJets[2].jetID==1 && patHTMHT.sumEt>450)"
    sel_bench = sel_noNoise + "*" + sel_bench

    kTrig = "PFMET100HT400"
    selections = [
        "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))",
        "(hltCaloHTMHT.sumEt>350 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))",
        "(hltCaloHTMHT.sumEt>350 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100) && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)",
        "(hltCaloHTMHT.sumEt>350 && hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && hltPFHTMHTNoPU.sumEt>400 && hltPFMET.pt>100 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)",
        ]
    entries = [
        "2012D",
        "2015#alpha (HT>400)",
        "2015#alpha (HT>400 && CJ40x2)",
        "2015#alpha (HT>400 && CJ40x2 && CaloMET>65)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])

    ## Change variable
    #variable = ("patMHT", "patHTMHT.pt")
    #binning = ("#scale[0.7]{RECO} PFMHT [GeV]", 50, 0, 250)
    #addsel = "(1)"
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])

    kTrig = "PFMET80HT400"
    selections = [
        "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))",
        "(hltCaloHTMHT.sumEt>350 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))",
        "(hltCaloHTMHT.sumEt>350 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100) && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)",
        "(hltCaloHTMHT.sumEt>350 && hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && hltPFHTMHTNoPU.sumEt>400 && hltPFMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)",
        ]
    entries = [
        "2012D",
        "2015#alpha (HT>400)",
        "2015#alpha (HT>400 && CJ40x2)",
        "2015#alpha (HT>400 && CJ40x2 && PFMET>80)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,3)(selections), entries=itemgetter(0,3)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])

    #___________________________________________________________________________
    # btag

    kColor = kBoson
    kTrig = "btag"

    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patGlobal.bjet_maxcsv>0.898)"
    sel_bench = sel_noNoise + "*" + sel_bench

    kTrig = "PFMET80CJ30x2CSV07"
    selections = [
        "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)",
        "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)",
        "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && hltPFMET.pt>80)",
        ]
    entries = [
        "2012D",
        "2015#alpha (HBHE)",
        "2015#alpha (HBHE && PFNoPU)",
        ]
    #legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    #save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,2)(selections), entries=itemgetter(0,2)(entries), numerator=True)
    #save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    #legsB = draw_future_B(selections=itemgetter(0,2)(selections), entries=itemgetter(0,2)(entries), numerator=False)
    #save(imgdir, "future_B_" + kTrig + "_" + variable[0])

    # mono
    sel_bench = "(Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>0 && patJets[0].pt>100 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1 && patGlobal.bjet_maxcsv>0.898)"
    sel_bench = sel_noNoise + "*" + sel_bench

    kTrig = "PFMET100CJ80CSV07"
    selections = [
        "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)",
        #"(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)",
        #"(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && hltPFMET.pt>80)",
        "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>80)",
        "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100)",
        ]
    entries = [
        "2012D",
        "2015#alpha (CJ80)",
        "2015#alpha (CJ80 && PFMET100)",
        ]
    legsA = draw_future_A(variable, binning, selections, entries, sel_bench, benchmark=True)
    save(imgdir, "future_A_" + kTrig + "_" + variable[0])
    legsB = draw_future_B(selections=itemgetter(0,2)(selections), entries=itemgetter(0,2)(entries), numerator=True)
    save(imgdir, "future_B2_" + kTrig + "_" + variable[0])
    legsB = draw_future_B(selections=itemgetter(0,2)(selections), entries=itemgetter(0,2)(entries), numerator=False)
    save(imgdir, "future_B_" + kTrig + "_" + variable[0])



