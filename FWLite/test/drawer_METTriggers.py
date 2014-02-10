#!/usr/bin/env python

from ROOT import TH1, TH1F, TH2F, TProfile, TFile, TChain, TCanvas, TLegend, TLatex, TLine, gROOT, gInterpreter, gStyle, gSystem, gPad
from rootcolors import *
from math import sqrt
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


# ______________________________________________________________________________
# Configurations
drawerInit = DrawerInit()
chain = TChain("tree", "tree")
infiles = [
    "../bin/compactified.L1ETM40.1.root"
    ]
for f in infiles:
    chain.Add(f)
tree = chain


sections = {}
sections["overview"]        = False
sections["overview_prof"]   = False
sections["topology"]        = False
sections["topology_hlt"]    = False
sections["puremet"]         = False
sections["puremet_eff"]     = False
sections["monojet"]         = False
sections["monojet_eff"]     = True
sections["higdijet"]        = False
sections["higdijet_eff"]    = False
sections["susdijet"]        = False
sections["susdijet_eff"]    = False
sections["multijet"]        = False
sections["multijet_eff"]    = False
sections["bjet"]            = False
sections["bjet_eff"]        = False
sections["vbf"]             = False
sections["vbf_eff"]         = False
sections["vbf2"]            = False
sections["vbf2_eff"]        = False
sections["inthepast"]       = False
sections["inthefuture"]     = False
plotting = []

#imgdir = "figures_20131130/"  # for Torino workshop
imgdir = "figures_20140206/"
if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)

wait = False
if not wait:  gROOT.SetBatch(1)


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
    tmpbins = [0,10,20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,280,320,360,400,480,600]
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

# Project
def fixOverflow(h):
    nbins = h.GetNbinsX()
    if h.GetBinContent(nbins+1) > 0:
        h.SetBinContent(nbins, h.GetBinContent(nbins) + h.GetBinContent(nbins+1))
        h.SetBinError(nbins, sqrt(h.GetBinError(nbins)**2 + h.GetBinError(nbins+1)**2))
        h.SetBinContent(nbins+1, 0)
        h.SetBinError(nbins+1, 0)
        h.SetEntries(h.GetEntries() - 2)  # SetBinEntries() automatically increases # entries by one

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

def draw_rate(params, histos, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", logy=False):
    ymax = histos[3].GetMaximum()
    histos[0].SetMaximum(ymax * 2.5)
    histos[0].SetMinimum(0)
    histos[0].GetYaxis().SetTitle(ytitle)

    histos[0].Draw("hist")
    histos[2].Draw("hist same")
    histos[3].Draw("hist same")

    gPad.SetLogy(logy)
    latex.DrawLatex(0.17, 0.97, text)
    return

def label_rate(histos, trig="Trigger", legend=(0.52,0.82,0.96,0.94)):
    leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
    leg1.SetFillStyle(0)
    leg1.SetLineColor(0)
    leg1.SetShadowColor(0)
    leg1.SetBorderSize(0)
    leg1.AddEntry(histos[0], "Noise", "f")
    leg1.AddEntry(histos[2], "Uncollected by %s" % trig, "f")
    leg1.AddEntry(histos[3], "%s" % trig, "f")
    leg1.Draw()
    latex.DrawLatex(0.56, 0.78, "All HLT filters except this one")
    return (leg1)

def draw_effnum(params, histos, k, trigs, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", ymin=0, ymax=-2.0, logy=False, legend=(0.46,0.625,0.96,0.94), normalize=False):
    if normalize:
        ytitle = "arbitrary unit"
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

    leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
    leg1.SetFillStyle(0)
    leg1.SetLineColor(0)
    leg1.SetShadowColor(0)
    leg1.SetBorderSize(0)
    for h, t in zip(histos, trigs):
        leg1.AddEntry(h, t, "f")
    leg1.Draw()
    latex.DrawLatex(0.17, 0.97, text)

    tmphistos = []
    if normalize:
        for h in histos:
            hclone = h.Clone(h.GetName() + "_clone")
            hclone.SetFillStyle(0)
            hclone.SetLineWidth(2)
            hclone.Draw("hist same")
            tmphistos.append(hclone)  # persistent
    return (leg1, tmphistos)  # persistent

def draw_eff(params, histos, k, trigs, ytitle="HLT Efficiency", text="Run2012D HLT_L1ETM40_v2", ymin=0, ymax=1.5, logy=False, legend=(0.48,0.695,0.96,0.94)):
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

    leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
    leg1.SetFillStyle(0)
    leg1.SetLineColor(0)
    leg1.SetShadowColor(0)
    leg1.SetBorderSize(0)
    for h, t in zip(histos[1:], trigs[1:]):
        leg1.AddEntry(h, t, "f")
    leg1.Draw()
    latex.DrawLatex(0.17, 0.97, text)
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

    def draw(params, histos, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", logy=False):
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
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def draw_norm(params, histos, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", logy=False):
        ymax = getMaximum(histos)
        histos[0].SetMaximum(ymax * 1.5)
        histos[0].GetYaxis().SetTitle(ytitle)

        histos[0].SetFillStyle(3004)
        histos[1].SetFillStyle(3003)
        histos[0].Draw("hist")
        histos[1].Draw("hist same")

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

    def label_norm(histos, legend=(0.26,0.84,0.90,0.94)):
        leg1 = TLegend(0.26,0.84,0.96,0.94)
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "1 #leq #PV #leq 10 (no noise) #mu,#sigma = %.1f,%.1f" % (histos[0].GetMean(), histos[0].GetRMS()), "f")
        leg1.AddEntry(histos[1], "25 #leq #PV #leq 35(no noise) #mu,#sigma = %.1f,%.1f" % (histos[1].GetMean(), histos[1].GetRMS()), "f")
        leg1.Draw()
        return (leg1)

    # nGoodPV
    variable = ("nGoodPV", "event.nGoodPV")
    binning = ("#scale[0.7]{RECO} # good PVs", 40, 0, 40)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # nJets
    variable = ("nJets", "Sum$(patJets.pt > 30 && abs(patJets.eta) < 2.5)")
    binning = ("#scale[0.7]{RECO} # jets_{p_{T} > 30 GeV, |#eta| < 2.5}", 5, 0, 5)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # Jet 1 pT
    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # Jet 2 pT
    variable = ("ptj2", "Alt$(patJets[1].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 2 p_{T} [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # Jet 1 eta
    variable = ("etaj1", "Alt$(patJets[0].eta, -99)")
    binning = ("#scale[0.7]{RECO} jet 1 #eta", 20, -5, 5)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # Jet 2 eta
    variable = ("etaj2", "Alt$(patJets[1].eta, -99)")
    binning = ("#scale[0.7]{RECO} jet 2 #eta", 20, -5, 5)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # HLT MET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # RECO MET
    variable = ("recoPFMET", "recoPFMET.pt")
    binning = ("#scale[0.7]{RECO} PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    variable = ("recoPFMETT1", "recoPFMETT1.pt")
    binning = ("#scale[0.7]{RECO} Type-1 PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # RECO MET variants
    variable = ("patMPT", "patMPT.pt")
    binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    variable = ("recoPFMETNoPU", "recoPFMETNoPU.pt")
    binning = ("#scale[0.7]{RECO} NoPU PFMET [GeV]", 30, 0, 150)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    save(imgdir, "overview_"+variable[0])

    # Two pileup scenarios
    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_lowPU = "(1<=event.nGoodPV && event.nGoodPV<=10)"
    sel_highPU = "(25<=event.nGoodPV && event.nGoodPV<=35)"


    # HLT MET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 0, 150)
    params = [
        (variable[0]+"_lowPU", kRed, kRed, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
        (variable[0]+"_highPU", kMaroon, kMaroon, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos, normalize=1)
    draw_norm(params, histos, ytitle="(normalized)")
    legs = label_norm(histos)
    save(imgdir, "overview_norm_"+variable[0])

    variable = ("hltTrackMET", "hltTrackMET.pt")
    binning = ("#scale[0.7]{HLT} TrackMET [GeV]", 30, 0, 150)
    params = [
        (variable[0]+"_lowPU", kBlue, kBlue, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
        (variable[0]+"_highPU", kNavy, kNavy, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos, normalize=1)
    draw_norm(params, histos, ytitle="(normalized)")
    legs = label_norm(histos)
    save(imgdir, "overview_norm_"+variable[0])

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 0, 150)
    params = [
        (variable[0]+"_lowPU", kRed, kRed, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
        (variable[0]+"_highPU", kMaroon, kMaroon, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos, normalize=1)
    draw_norm(params, histos, ytitle="(normalized)")
    legs = label_norm(histos)
    save(imgdir, "overview_norm_"+variable[0])

    variable = ("patMPT", "patMPT.pt")
    binning = ("#scale[0.7]{RECO} TrackMET [GeV]", 30, 0, 150)
    params = [
        (variable[0]+"_lowPU", kBlue, kBlue, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
        (variable[0]+"_highPU", kNavy, kNavy, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos, normalize=1)
    draw_norm(params, histos, ytitle="(normalized)")
    legs = label_norm(histos)
    save(imgdir, "overview_norm_"+variable[0])

    variable = ("recoPFMETMVA", "recoPFMETMVA.pt")
    binning = ("#scale[0.7]{RECO} MVA PFMET [GeV]", 30, 0, 150)
    params = [
        (variable[0]+"_lowPU", kMagenta, kMagenta, variable[1], "*".join([sel, sel_lowPU, sel_noNoise])),
        (variable[0]+"_highPU", kPurple2, kPurple2, variable[1], "*".join([sel, sel_highPU, sel_noNoise])),
    ]
    histos = book(params, binning)
    project(params, histos, normalize=1)
    draw_norm(params, histos, ytitle="(normalized)")
    legs = label_norm(histos)
    save(imgdir, "overview_norm_"+variable[0])


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

    def draw(params, histos, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", logy=False, ymin=1e-6, ymax=1e2):
        histos[0].SetMinimum(ymin)
        histos[0].SetMaximum(ymax)
        #histos[0].GetYaxis().SetTitle(ytitle)

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

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)

    # Mean
    params = [
        ("hltCaloMET"    , kBlack, kBlack, "hltCaloMET.pt:event.nGoodPV" , "*".join([sel, sel_noNoise])),
        ("hltPFMET"      , kRed  , kRed  , "hltPFMET.pt:event.nGoodPV"   , "*".join([sel, sel_noNoise])),
        ("hltTrackMET"   , kBlue , kBlue , "hltTrackMET.pt:event.nGoodPV", "*".join([sel, sel_noNoise])),
    ]
    binning = ("#scale[0.7]{RECO} # good PVs; #scale[0.7]{HLT} #mu_{MET} [GeV]", 20, 0, 40, 0, 200)
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
    binning = ("#scale[0.7]{RECO} # good PVs; #scale[0.7]{RECO} #mu_{MET} [GeV]", 20, 0, 40, 0, 200)
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
    binning = ("#scale[0.7]{RECO} # good PVs; #scale[0.7]{HLT} #sigma_{MET} [GeV]", 20, 0, 40, 0, 200)
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
    binning = ("#scale[0.7]{RECO} # good PVs; #scale[0.7]{RECO} #sigma_{MET} [GeV]", 20, 0, 40, 0, 200)
    histos = book_prof(params, binning, option="s")
    project_prof(params, histos)
    hhistos = rmserror(binning, histos)
    draw(params, hhistos, ymax=40)
    legs = label_RECO(histos)
    save(imgdir, "overview_prof_sigma_recoMETs")


# ______________________________________________________________________________
# Topology
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

    def draw(params, histos, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", logy=False, zoom=False):
        ymax = getMaximum(histos)
        if zoom:
            histos[0].SetMaximum(ymax * 1.5 / 100)
        else:
            histos[0].SetMaximum(ymax * 1.5)
        histos[0].GetYaxis().SetTitle(ytitle)

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
        leg1.AddEntry(histos[3], "DiCentralJet", "f")
        leg1.AddEntry(histos[4], "HT", "f")
        leg1.Draw()
        return (leg1)

    # HLT MET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 30, 50, 200)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save(imgdir, "topology_"+variable[0])

    draw(params, histos, zoom=True)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save(imgdir, "topology_"+variable[0]+"_zoom")

    # RECO MET
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 50, 200)
    params = prepare(variable)
    histos = book(params, binning)
    project(params, histos)
    draw(params, histos)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save(imgdir, "topology_"+variable[0])

    draw(params, histos, zoom=True)
    legs = label(histos)
    latex.DrawLatex(0.56, 0.66, "categorized by jet topology (offline)")
    save(imgdir, "topology_"+variable[0]+"_zoom")


if sections["topology_hlt"]:
    if plotting: del plotting[:]

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1  = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
    sel_trig2a = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
    sel_trig2b = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) )
    sel_trig2  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4'), triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')))
    sel_trig3  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
    sel_trig4a = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) )
    sel_trig4b = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) )
    sel_trig4  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5'), triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')))
    sel_trig5a = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) )
    sel_trig5b = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) )
    sel_trig5  = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9'), triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')))

    def draw(params, histos, ytitle="Events", text="Run2012D HLT_L1ETM40_v2", logy=False, zoom=False):
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
        gPad.RedrawAxis()
        latex.DrawLatex(0.17, 0.97, text)
        return

    def label(histos, trig="MET trigger", legend=(0.52,0.82,0.96,0.94), noNoise=False):
        rate0 = (histos[0].Integral() - histos[1].Integral()) / histos[0].Integral() * 100.
        rate1 = (histos[1].Integral() - histos[2].Integral()) / histos[0].Integral() * 100.
        rate2 = (histos[2].Integral()) / histos[0].Integral() * 100.
        leg1 = TLegend(legend[0], legend[1], legend[2], legend[3])
        leg1.SetFillStyle(0)
        leg1.SetLineColor(0)
        leg1.SetShadowColor(0)
        leg1.SetBorderSize(0)
        leg1.AddEntry(histos[0], "A && !B (%.0f%%)" % rate0, "f")
        leg1.AddEntry(histos[1], "A && B (%.0f%%)" % rate1, "f")
        leg1.AddEntry(histos[2], "!A && B (%.0f%%)" % rate2, "f")
        leg1.Draw()
        latex.DrawLatex(0.56, 0.78, "A=%s  B=existing" % trig)
        if noNoise:
            latex.DrawLatex(0.56, 0.75, "after offline noise cleaning")
        return (leg1)

    kColor = kYellow
    kTrig = "MonoCentralJet"
    sel_trigC  = "(%s &&  %s)" % (sel_trig0, sel_trig1)
    sel_trigD  = "(%s ||  %s)" % (sel_trig0, sel_trig1)
    sel_trigB  = "(%s && !%s)" % (sel_trig0, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trig0, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kOrange3
    kTrig = "DiCentralJetSUS"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig2a)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig2a)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kOrange3
    kTrig = "DiCentralJetHIG"
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

    kColor = kCyan
    kTrig = "btag"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig4a)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig4a)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kMagenta
    kTrig = "VBFAll"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig5a)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig5a)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))

    kColor = kMagenta
    kTrig = "VBFLead"
    sel_trigC  = "(%s &&  %s)" % (sel_trigA, sel_trig5b)
    sel_trigD  = "(%s ||  %s)" % (sel_trigA, sel_trig5b)
    sel_trigB  = "(%s && !%s)" % (sel_trigA, sel_trigC)  # bottom
    sel_trigAB = "(%s ||  %s)" % (sel_trigA, sel_trigC)  # middle
    sel_trigA  = "(%s)"        % (sel_trigD)  # top
    plotting.append((sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig))


    for p in plotting:
        # HLT PFMET
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        variable = ("hltPFMET", "hltPFMET.pt")
        binning = ("#scale[0.7]{HLT} PFMET [GeV]", 50, 0, 250)
        addsel = "(1)"
        params = [
            (variable[0]+"_nofilt_trigA" , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])),
            (variable[0]+"_nofilt_trigAB", kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB])),
            (variable[0]+"_nofilt_trigB" , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw(params, histos)
        legs = label(histos, kTrig)
        save(imgdir, "topology_hlt_nofilt_"+kTrig+"_"+variable[0])

        del params[:]
        params = [
            (variable[0]+"_filt_trigA"   , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA , sel_noNoise])),
            (variable[0]+"_filt_trigAB"  , kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB, sel_noNoise])),
            (variable[0]+"_filt_trigB"   , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB , sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw(params, histos)
        legs = label(histos, kTrig, noNoise=True)
        save(imgdir, "topology_hlt_"+kTrig+"_"+variable[0])

        # RECO PFMET
        (sel_trigA, sel_trigAB, sel_trigB, kColor, kTrig) = p
        variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
        binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 50, 0, 250)
        addsel = "(1)"
        params = [
            (variable[0]+"_nofilt_trigA" , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA])),
            (variable[0]+"_nofilt_trigAB", kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB])),
            (variable[0]+"_nofilt_trigB" , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw(params, histos)
        legs = label(histos, kTrig)
        save(imgdir, "topology_hlt_nofilt_"+kTrig+"_"+variable[0])

        del params[:]
        params = [
            (variable[0]+"_filt_trigA"   , kBlack, kColor   , variable[1], "*".join([sel, addsel, sel_trigA , sel_noNoise])),
            (variable[0]+"_filt_trigAB"  , kBlack, nBlue    , variable[1], "*".join([sel, addsel, sel_trigAB, sel_noNoise])),
            (variable[0]+"_filt_trigB"   , kBlack, kBlue    , variable[1], "*".join([sel, addsel, sel_trigB , sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw(params, histos)
        legs = label(histos, kTrig, noNoise=True)
        save(imgdir, "topology_hlt_"+kTrig+"_"+variable[0])


# ______________________________________________________________________________
if sections["puremet"]:
    if plotting: del plotting[:]

    kColor = kBlue
    kTrig = "PFMET150"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"
    sel_trig2 = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>70 && hltPFMET.pt>150)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 130, 300)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 34, 130, 300)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} PFMET [GeV]", 34, 130, 300)
    addsel = "(hltCaloMET.pt>80 && hltPFMET.pt>-99)"
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
    addsel = "(hltCaloMET.pt>-99 && hltPFMET.pt>150)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && hltPFMET.pt>150)"
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

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = sel_trig1
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 300)
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
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["monojet"]:
    if plotting: del plotting[:]

    kColor = kYellow
    kTrig = "MonoCentralJet"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
    sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 80, 250)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 34, 80, 250)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} PFMET [GeV]", 34, 80, 250)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>-99)"
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
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>-99 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>-99 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_pt2", "MaxIf$(hltPFJetsL1FastL2L3.pt, hltPFJetsL1FastL2L3.pt!=Max$(hltPFJetsL1FastL2L3.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt2", "MaxIf$(hltPFJetsL1FastL2L3NoPU.pt, hltPFJetsL1FastL2L3NoPU.pt!=Max$(hltPFJetsL1FastL2L3NoPU.pt))")
    binning = ("#scale[0.7]{HLT} Subleading PFNoPUJet p_{T} [GeV]", 30, 0, 150)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_eta1", "Sum$(hltPFJetsL1FastL2L3.eta * (hltPFJetsL1FastL2L3.pt==Max$(hltPFJetsL1FastL2L3.pt)))")
    binning = ("#scale[0.7]{HLT} Leading PFJet #eta", 20, -5, 5)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_pt2", "MaxIf$(hltCaloJetsL1Fast.pt, hltCaloJetsL1Fast.pt!=Max$(hltCaloJetsL1Fast.pt))")
    binning = ("#scale[0.7]{HLT} Subleading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
    plotting.append((variable, binning, addsel))

    variable = ("hltCaloJetsL1Fast_eta1", "Sum$(hltCaloJetsL1Fast.eta * (hltCaloJetsL1Fast.pt==Max$(hltCaloJetsL1Fast.pt)))")
    binning = ("#scale[0.7]{HLT} Leading CaloJet #eta", 20, -5, 5)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"
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

    # nef, chf, cef, nch, ntot
    variable = ("hltPFJetsL1FastL2L3_nef1", "hltPFJetsL1FastL2L3[0].nef")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(neutral EM)", 40, 0, 1)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>40 && hltCaloMETCleanUsingJetID.pt>40)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_chf1", "hltPFJetsL1FastL2L3[0].chf")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(charged HAD)", 40, 0, 1)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>40 && hltCaloMETCleanUsingJetID.pt>40)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_cef1", "hltPFJetsL1FastL2L3[0].cef")
    binning = ("#scale[0.7]{HLT} Leading PFJet f_{en}(charged EM)", 40, 0, 1)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>40 && hltCaloMETCleanUsingJetID.pt>40)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_nch1", "hltPFJetsL1FastL2L3[0].nch")
    binning = ("#scale[0.7]{HLT} Leading PFJet # charged constits", 40, 0, 40)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>40 && hltCaloMETCleanUsingJetID.pt>40)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3_ntot1", "hltPFJetsL1FastL2L3[0].ntot")
    binning = ("#scale[0.7]{HLT} Leading PFJet # constits", 40, 0, 80)
    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105 && hltCaloMETClean.pt>40 && hltCaloMETCleanUsingJetID.pt>40)"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["higdijet"]:
    if plotting: del plotting[:]

    kColor = kOrange3
    kTrig = "DiCentralJetHIG"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) )
    sel_trig1 = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETCleanUsingJetID", "hltPFMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT} PFMETCleanUsingJetID [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && hltCaloMETClean.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && hltCaloMETClean.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && hltCaloMETClean.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["susdijet"]:
    if plotting: del plotting[:]

    kColor = kOrange3
    kTrig = "DiCentralJetSUS"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
    sel_trig1 = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>-99||hltPFMETNoMu.pt>-99) )"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>-99||hltPFMETNoMu.pt>-99) )"
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>-99)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>-99)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["multijet"]:
    if plotting: del plotting[:]

    kColor = kCocoa2
    kTrig = "HT"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
    sel_trig1 = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 34, 80, 250)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>-99))"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 34, 80, 250)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>-99))"
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 40, 0, 400)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 40, 0, 400)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 40, 0, 400)
    addsel = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["bjet"]:
    if plotting: del plotting[:]

    kColor = kCyan
    kTrig = "btag"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) )
    sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    #sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) )
    #sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 28, 60, 200)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>-99)"
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>-99)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>-99)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["vbf"]:
    if plotting: del plotting[:]

    kColor = kMagenta
    kTrig = "VBFAll"

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i])" % len(metfilters)
    #sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) )
    #sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_leadmjj>600 && hltPFGlobal.vbf_leadmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) )
    sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"

    # hltPFMET
    variable = ("hltPFMET", "hltPFMET.pt")
    binning = ("#scale[0.7]{HLT} PFMET [GeV]", 26, 50, 180)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>-99)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFMETNoMu", "hltPFMETNoMu.pt")
    binning = ("#scale[0.7]{HLT} PFMETNoMu [GeV]", 26, 50, 180)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>-99)"
    plotting.append((variable, binning, addsel))

    # hltCaloMET
    variable = ("hltCaloMET", "hltCaloMET.pt")
    binning = ("#scale[0.7]{HLT} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETClean
    variable = ("hltCaloMETClean", "hltCaloMETClean.pt")
    binning = ("#scale[0.7]{HLT HBHENoiseCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    # hltCaloMETCleanUsingJetID
    variable = ("hltCaloMETCleanUsingJetID", "hltCaloMETCleanUsingJetID.pt")
    binning = ("#scale[0.7]{HLT JetIDCleaned} CaloMET [GeV]", 40, 0, 200)
    addsel = "(hltCaloMET.pt>-99 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    # hltPFJet
    variable = ("hltPFJetsL1FastL2L3_pt1", "Max$(hltPFJetsL1FastL2L3.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    variable = ("hltPFJetsL1FastL2L3NoPU_pt1", "Max$(hltPFJetsL1FastL2L3NoPU.pt)")
    binning = ("#scale[0.7]{HLT} Leading PFNoPUJet p_{T} [GeV]", 30, 0, 300)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    # hltCaloJet
    variable = ("hltCaloJetsL1Fast_pt1", "Max$(hltCaloJetsL1Fast.pt)")
    binning = ("#scale[0.7]{HLT} Leading CaloJet p_{T} [GeV]", 30, 0, 150)
    addsel = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>-99)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"
    plotting.append((variable, binning, addsel))

    for p in plotting:
        (variable, binning, addsel) = p
        params = [
            (variable[0]+"_nofilt"      , kBlack, kWhite , variable[1], "*".join([sel, addsel])),
            (variable[0]+"_nofilt_trig1", kBlack, kWhite , variable[1], "*".join([sel, addsel, sel_trig1])),
            (variable[0]+"_filt"        , kBlack, kGray  , variable[1], "*".join([sel, addsel, sel_noNoise])),
            (variable[0]+"_filt_trig1"  , kBlack, kColor , variable[1], "*".join([sel, addsel, sel_trig1, sel_noNoise])),
            ]
        histos = book(params, binning)
        project(params, histos)
        draw_rate(params, histos)
        legs = label_rate(histos, trig=kTrig)
        save(imgdir, "diagnosis_" + kTrig + "_" + variable[0])


if sections["puremet_eff"]:
    if plotting: del plotting[:]

    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise = "(metfilterFlags[%i] && (Sum$(patJets.pt>30)>0 && patJets[0].jetID==1 && ((Sum$(patJets.pt>30)>1 && patJets[1].jetID==1) || Sum$(patJets.pt>30)==1)) )" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
    sel_trig1 = "(hltCaloMET.pt>%i && hltPFMET.pt>%i)"
    sel_trig2 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltPFMET.pt>%i)"
    sel_trig3 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltCaloMETCleanUsingJetID.pt>%i && hltPFMET.pt>%i)"
    sel_trig4 = "(hltCaloMET.pt>%i && hltCaloMETClean.pt>%i && hltCaloMETCleanUsingJetID.pt>%i && hltPFMET.pt>%i && hltTrackMET.pt>%i)"

    colors = [
        TColor.GetColor("#9CE0EB"),  #? use TColor.GetColor("#CCF5FF"),
        TColor.GetColor("#66D6E5"),
        TColor.GetColor("#20ABBF"),
        TColor.GetColor("#0D8799"),
        TColor.GetColor("#055A66"),
        TColor.GetColor("#002E33"),
        ]
    kComplete = TColor.GetColor("#D67F20")
    addsel = "(1)"

    # recoPFMETT0T1
    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 80, 440)
    thresholds = [0, 80, 90, 100, 110, 120]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 200, 0
    kHist = thresholds.index(90)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 80, 400)
    thresholds = [0, 65, 70, 80, 90, 100]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 150, 0
    kHist = thresholds.index(90)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 70, 400)
    thresholds = [0, 55, 65, 70, 80, 90]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 120, 0
    kHist = thresholds.index(80)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 50, 360)
    thresholds = [0, 50, 55, 65, 70, 80]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 100, 0
    kHist = thresholds.index(70)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 40, 320)
    thresholds = [0, 50, 55, 60, 65, 70]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 80, 0
    kHist = thresholds.index(60)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 240)
    thresholds = [0, 45, 50, 55, 60, 65]
    calomet, calomet1, calomet2, pfmet, trackmet = 0, 0, 0, 65, 0
    kHist = thresholds.index(45)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    #for p in plotting:
    #    (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p
    #
    #    params = [
    #        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
    #        ]
    #    for i, x in enumerate(thresholds):
    #        params += [
    #            (variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig4 % (x, calomet1, calomet2, pfmet, trackmet), sel_noNoise])),
    #            ]
    #    kTrig = "PFMET%i" % pfmet
    #    trigs = ["HLT PFMET<%i" % pfmet]+[("HLT PFMET>%i && CaloMET>%i" % (pfmet,x)) for x in thresholds]
    #
    #    histos = book_ratio(params, binning)
    #    project(params, histos, drawOverflow=False)
    #    legs = draw_effnum(params, histos, kHist, trigs)
    #    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #    legs = draw_eff(params, histos, kHist, trigs)
    #    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # CaloMETClean or CaloMETCleanUsingJetID optimization
    if plotting: del plotting[:]

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 80, 440)
    thresholds = [0, 60, 70, 80, 90, 100]
    calomet, calomet1, calomet2, pfmet, trackmet = 100, 0, 0, 220, 0
    kHist = thresholds.index(70)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 80, 400)
    thresholds = [0, 50, 60, 70, 80, 90]
    calomet, calomet1, calomet2, pfmet, trackmet = 90, 0, 0, 150, 0
    kHist = thresholds.index(70)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 70, 400)
    thresholds = [0, 40, 50, 60, 70, 80]
    calomet, calomet1, calomet2, pfmet, trackmet = 80, 0, 0, 120, 0
    kHist = thresholds.index(50)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 50, 360)
    thresholds = [0, 40, 45, 50, 60, 70]
    calomet, calomet1, calomet2, pfmet, trackmet = 70, 0, 0, 100, 0
    kHist = thresholds.index(50)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 40, 320)
    thresholds = [0, 30, 40, 45, 50, 60]
    calomet, calomet1, calomet2, pfmet, trackmet = 60, 0, 0, 80, 0
    kHist = thresholds.index(50)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 240)
    thresholds = [0, 30, 35, 40, 45, 50]
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 0, 0, 65, 0
    kHist = thresholds.index(50)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    #for p in plotting:
    #    (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p
    #
    #    params = [
    #        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
    #        ]
    #    for i, x in enumerate(thresholds):
    #        params += [
    #            (variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, x, calomet2, pfmet, trackmet), sel_noNoise])),
    #            #(variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet, x, pfmet, trackmet), sel_noNoise])),
    #            ]
    #    kTrig = "CaloMET%i" % calomet
    #    #kTrig = "CaloMET%i_2" % calomet
    #    trigs = ["HLT CaloMET<%i" % pfmet]+[("HLT CaloMET>%i && CaloMETClean>%i" % (calomet,x)) for x in thresholds]
    #
    #    histos = book_ratio(params, binning)
    #    project(params, histos, drawOverflow=False)
    #    legs = draw_effnum(params, histos, kHist, trigs)
    #    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #    legs = draw_eff(params, histos, kHist, trigs)
    #    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])

    #___________________________________________________________________________
    # TrackMET optimization
    if plotting: del plotting[:]

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 80, 440)
    thresholds = [0, 10, 15, 20, 30, 40]
    calomet, calomet1, calomet2, pfmet, trackmet = 90, 0, 0, 200, 0
    kHist = thresholds.index(20)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 80, 400)
    thresholds = [0, 10, 15, 20, 25, 30]
    calomet, calomet1, calomet2, pfmet, trackmet = 90, 0, 0, 150, 0
    kHist = thresholds.index(15)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 70, 400)
    thresholds = [0, 10, 12, 15, 20, 25]
    calomet, calomet1, calomet2, pfmet, trackmet = 80, 0, 0, 120, 0
    kHist = thresholds.index(15)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 50, 360)
    thresholds = [0, 5, 8, 10, 12, 15]
    calomet, calomet1, calomet2, pfmet, trackmet = 70, 0, 0, 100, 0
    kHist = thresholds.index(10)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 40, 320)
    thresholds = [0, 5, 8, 10, 12, 15]
    calomet, calomet1, calomet2, pfmet, trackmet = 60, 0, 0, 80, 0
    kHist = thresholds.index(10)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    variable = ("recoPFMETT0T1", "recoPFMETT0T1.pt")
    binning = ("#scale[0.7]{RECO} Type-0+1 PFMET [GeV]", 30, 240)
    thresholds = [0, 5, 8, 10, 12, 15]
    calomet, calomet1, calomet2, pfmet, trackmet = 50, 0, 0, 65, 0
    kHist = thresholds.index(10)
    plotting.append((variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist))

    #for p in plotting:
    #    (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p
    #
    #    params = [
    #        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
    #        ]
    #    for i, x in enumerate(thresholds):
    #        params += [
    #            (variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet1, calomet2, pfmet, x), sel_noNoise])),
    #            ]
    #    kTrig = "PFMET%i_2" % pfmet
    #    trigs = ["HLT (PFMET<%i || CaloMET<%i)" % (pfmet, calomet)]+[("HLT PFMET>%i && CaloMET>%i && TrackMET>%i" % (pfmet,calomet,x)) for x in thresholds]
    #
    #    histos = book_ratio(params, binning)
    #    project(params, histos, drawOverflow=False)
    #    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.23,0.625,0.96,0.94))
    #    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #    legs = draw_eff(params, histos, kHist, trigs, legend=(0.24,0.695,0.96,0.94))
    #    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


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

    #for p in plotting:
    #    (variable, binning, thresholds, calomet, calomet1, calomet2, pfmet, trackmet, kHist) = p
    #
    #    params = []
    #    for i, x in enumerate(thresholds):
    #        params += [
    #            (variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig4 % (calomet, calomet1, calomet2, x, trackmet), sel_noNoise])),
    #            ]
    #    kTrig = "PFMET%i" % pfmet
    #    trigs = [("HLT CaloMET>%i && PFMET>%i" % (calomet,x)) for x in thresholds]
    #
    #    histos = book_ratio(params, binning)
    #    project(params, histos, drawOverflow=False)
    #    legs = draw_effnum(params, histos, kHist, trigs, normalize=True)
    #    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])


if sections["monojet_eff"]:
    if plotting: del plotting[:]

    sel = "(triggerFlags[%i])" % ireftrig
    #sel_noNoise = "(metfilterFlags[%i] && (Sum$(patJets.pt>0)>0 && patJets[0].jetID==1 && ((Sum$(patJets.pt>30)>1 && patJets[1].jetID==1) || Sum$(patJets.pt>30)==1)) )" % len(metfilters)
    sel_noNoise = "(metfilterFlags[%i] && (Sum$(patJets.pt>0)>0 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1) )" % len(metfilters)
    sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
    sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>%i)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>%i)>0 && hltPFMETNoMu.pt>105)"
    sel_trig2 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>%i)>0 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>%i)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>%i)>0)"
    sel_trig3 = "(hltPFJetsL1FastL2L3[0].nhf<%.2f && ((abs(hltPFJetsL1FastL2L3[0].eta)<2.4 && hltPFJetsL1FastL2L3[0].nch>%i) || abs(hltPFJetsL1FastL2L3[0].eta)>2.4))"  # use |eta|<2.4

    colors = [
        TColor.GetColor("#9CE0EB"),  #? use TColor.GetColor("#CCF5FF"),
        TColor.GetColor("#66D6E5"),
        TColor.GetColor("#20ABBF"),
        TColor.GetColor("#0D8799"),
        TColor.GetColor("#055A66"),
        TColor.GetColor("#002E33"),
        ]
    kComplete = TColor.GetColor("#D67F20")
    addsel = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && hltPFMET.pt>100)"

    # ptj1
    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 320)
    thresholds = [0, 90, 100, 105, 110, 120]
    calojet, pfjet, pfnopujet = 0, 120, 0
    kHist = thresholds.index(100)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 280)
    thresholds = [0, 70, 80, 85, 90, 100]
    calojet, pfjet, pfnopujet = 0, 100, 0
    kHist = thresholds.index(80)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 240)
    thresholds = [0, 50, 60, 65, 70, 80]
    calojet, pfjet, pfnopujet = 0, 80, 0
    kHist = thresholds.index(60)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 40, 50, 55, 60, 70]
    calojet, pfjet, pfnopujet = 0, 70, 0
    kHist = thresholds.index(50)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 30, 40, 45, 50, 60]
    calojet, pfjet, pfnopujet = 0, 60, 0
    kHist = thresholds.index(40)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 20, 30, 35, 40, 50]
    calojet, pfjet, pfnopujet = 0, 50, 0
    kHist = thresholds.index(30)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    #for p in plotting:
    #    (variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist) = p
    #
    #    params = [
    #        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
    #        ]
    #    for i, x in enumerate(thresholds):
    #        params += [
    #            (variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig2 % (x, pfjet, pfnopujet), sel_noNoise])),
    #            ]
    #    kTrig = "PFJet%i" % pfjet
    #    trigs = ["HLT PFJet1 p_{T}<%i" % pfjet]+[("HLT PFJet1 p_{T}>%i && CaloJet1 p_{T}>%i" % (pfjet,x)) for x in thresholds]
    #
    #    histos = book_ratio(params, binning)
    #    project(params, histos, drawOverflow=False)
    #    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.48,0.695,0.94,0.94))
    #    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #    legs = draw_eff(params, histos, kHist, trigs)
    #    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])

    #___________________________________________________________________________
    # PFnoPUJet
    if plotting: del plotting[:]

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 320)
    thresholds = [0, 90, 100, 105, 110, 120]
    calojet, pfjet, pfnopujet = 0, 0, 120
    kHist = thresholds.index(100)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 280)
    thresholds = [0, 70, 80, 85, 90, 100]
    calojet, pfjet, pfnopujet = 0, 0, 100
    kHist = thresholds.index(80)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 240)
    thresholds = [0, 50, 60, 65, 70, 80]
    calojet, pfjet, pfnopujet = 0, 0, 80
    kHist = thresholds.index(60)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 40, 50, 55, 60, 70]
    calojet, pfjet, pfnopujet = 0, 0, 70
    kHist = thresholds.index(50)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 30, 40, 45, 50, 60]
    calojet, pfjet, pfnopujet = 0, 0, 60
    kHist = thresholds.index(40)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 20, 200)
    thresholds = [0, 20, 30, 35, 40, 50]
    calojet, pfjet, pfnopujet = 0, 0, 50
    kHist = thresholds.index(30)
    plotting.append((variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist))

    #for p in plotting:
    #    (variable, binning, thresholds, calojet, pfjet, pfnopujet, kHist) = p
    #
    #    params = [
    #        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
    #        ]
    #    for i, x in enumerate(thresholds):
    #        params += [
    #            (variable[0]+("_trig%i" % i), colors[i], colors[i], variable[1], "*".join([sel, addsel, sel_trig2 % (x, pfjet, pfnopujet), sel_noNoise])),
    #            ]
    #    kTrig = "PFNoPUJet%i" % pfnopujet
    #    trigs = ["HLT PFNoPUJet1 p_{T}<%i" % pfnopujet]+[("HLT PFNoPUJet1 p_{T}>%i && CaloJet1 p_{T}>%i" % (pfnopujet,x)) for x in thresholds]
    #
    #    histos = book_ratio(params, binning)
    #    project(params, histos, drawOverflow=False)
    #    legs = draw_effnum(params, histos, kHist, trigs, legend=(0.48,0.695,0.94,0.94))
    #    save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #    legs = draw_eff(params, histos, kHist, trigs)
    #    save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


    #___________________________________________________________________________
    # nhf & nch
    if plotting: del plotting[:]

    addsel = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>60)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>0)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && hltPFMET.pt>100)"

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 320)
    thresholds = [99, 0.95, 0.90]
    nhf, nch = 99, -99
    kHist = thresholds.index(0.95)
    params = [
        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
        ]
    for i, x in enumerate(thresholds):
        params += [
            (variable[0]+("_trig%i" % i), colors[1+i*2], colors[1+i*2], variable[1], "*".join([sel, addsel, sel_trig3 % (x, nch), sel_noNoise])),
            ]
    kTrig = "nhf"
    trigs = ["HLT PFJet1 nch<%i" % nch]+[("HLT PFJet1 nch>%i && nhf<%.2f" % (nch,x)) for x in thresholds]

    #histos = book_ratio(params, binning)
    #project(params, histos, drawOverflow=False)
    #legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    #save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94), ymin=0.8, ymax=1.1)
    #save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])

    variable = ("ptj1", "Alt$(patJets[0].pt, 0)")
    binning = ("#scale[0.7]{RECO} jet 1 p_{T} [GeV]", 30, 320)
    thresholds = [-1, 0, 1]
    nhf, nch = 0.95, -99
    kHist = thresholds.index(0)
    params = [
        (variable[0]+"_filt"        , kComplete, kComplete, variable[1], "*".join([sel, addsel, sel_noNoise])),
        ]
    for i, x in enumerate(thresholds):
        params += [
            (variable[0]+("_trig%i" % i), colors[1+i*2], colors[1+i*2], variable[1], "*".join([sel, addsel, sel_trig3 % (nhf, x), sel_noNoise])),
            ]
    kTrig = "nch"
    trigs = ["HLT PFJet1 nhf>%.2f" % nhf]+[("HLT PFJet1 nhf>%.2f && nch>%i" % (nhf,x)) for x in thresholds]

    #histos = book_ratio(params, binning)
    #project(params, histos, drawOverflow=False)
    #legs = draw_effnum(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94))
    #save(imgdir, "diagnosis_effnum_" + kTrig + "_" + variable[0])
    #legs = draw_eff(params, histos, kHist, trigs, legend=(0.40,0.79,0.96,0.94), ymin=0.8, ymax=1.1)
    #save(imgdir, "diagnosis_eff_" + kTrig + "_" + variable[0])


