from drawer_init import *

drawerInit = DrawerInit()
wait = True
if not wait:  gROOT.SetBatch(1)

#tfile = TFile.Open("../bin/compactified.root")
tfile = TFile.Open("../bin/compactified.0.root")
#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
tree = tfile.Events


def book(varname, xtitle, nbins, xlow, xup):
    histos = []
    histos.append(TH1F("ha_"+varname, "; "+xtitle, nbins, xlow, xup))
    histos.append(TH1F("hb_"+varname, "; "+xtitle, nbins, xlow, xup))
    histos.append(TH1F("hc_"+varname, "; "+xtitle, nbins, xlow, xup))

    histos[0].SetLineWidth(2)
    histos[1].SetLineWidth(2)
    histos[2].SetLineWidth(2)
    histos[0].SetLineColor(kBlue)
    histos[1].SetLineColor(kGreen)
    histos[2].SetLineColor(kRed)
    histos[0].SetFillStyle(3003)
    histos[1].SetFillStyle(3003)
    histos[2].SetFillStyle(3003)
    histos[0].SetFillColor(kBlue2)
    histos[1].SetFillColor(kGreen2)
    histos[2].SetFillColor(kRed2)
    return histos

def project(var, histos, addsel="1", normalize=False):
    tree.Project(histos[0].GetName(), var, runsel1 + " * triggerFlags[%i] * (%s)" % (i_reftrig, addsel), "goff")
    tree.Project(histos[1].GetName(), var, runsel2 + " * triggerFlags[%i] * (%s)" % (i_reftrig, addsel), "goff")
    tree.Project(histos[2].GetName(), var, runsel3 + " * triggerFlags[%i] * (%s)" % (i_reftrig, addsel), "goff")

    if normalize:
        histos[0].Scale(1.0 / histos[0].GetSumOfWeights())
        histos[1].Scale(1.0 / histos[1].GetSumOfWeights())
        histos[2].Scale(1.0 / histos[2].GetSumOfWeights())
        histos[0].GetYaxis().SetTitle("(normalized)")
    else:
        histos[0].GetYaxis().SetTitle("count")

def draw(histos, imgname, fnum="%.1f"):
    print "(A) N=%i mu=%.3f sig=%.3f" % (histos[0].GetEntries(), histos[0].GetMean(), histos[0].GetRMS())
    print "(B) N=%i mu=%.3f sig=%.3f" % (histos[1].GetEntries(), histos[1].GetMean(), histos[1].GetRMS())
    print "(C) N=%i mu=%.3f sig=%.3f" % (histos[2].GetEntries(), histos[2].GetMean(), histos[2].GetRMS())

    histos[0].SetMaximum(max(histos[0].GetMaximum(), histos[1].GetMaximum(), histos[2].GetMaximum()) * 1.1)
    histos[0].Draw()
    histos[1].Draw("same")
    histos[2].Draw("same")

    latex.DrawLatex(0.74, 0.90, "#color[4]{#mu,#sigma} = " + fnum % histos[0].GetMean() + ", " + fnum % histos[0].GetRMS() )
    latex.DrawLatex(0.74, 0.86, "#color[3]{#mu,#sigma} = " + fnum % histos[1].GetMean() + ", " + fnum % histos[1].GetRMS() )
    latex.DrawLatex(0.74, 0.82, "#color[2]{#mu,#sigma} = " + fnum % histos[2].GetMean() + ", " + fnum % histos[2].GetRMS() )

    gPad.Print(imgdir + imgname + ".png")
    gPad.Print(imgdir + imgname + ".pdf")


# NPV
histos_NPV = book("NPV", "# good PVs", 40, 0, 40)
project("event.nGoodPV", histos_NPV, "1", True)
draw(histos_NPV, "NPV_nofilt_MET_R207454", "%.1f")
project("event.nGoodPV", histos_NPV, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_NPV, "NPV_MET_R207454", "%.1f")

# NJets
histos_NJets = book("NJets", "# jets_{p_{T} > 30 GeV, |#eta| < 2.5}", 10, 0, 10)
project("Sum$(patJets.pt > 30 && abs(patJets.eta) < 2.5)", histos_NJets, "1", True)
draw(histos_NJets, "NJets_nofilt_MET_R207454", "%.1f")
project("Sum$(patJets.pt > 30 && abs(patJets.eta) < 2.5)", histos_NJets, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_NJets, "NJets_MET_R207454", "%.1f")

# HLT MET
histos_hltCaloMET = book("hltCaloMET", "HLT Calo MET [GeV]", 100, 0, 200)
project("hltCaloMET.pt", histos_hltCaloMET, "1", True)
draw(histos_hltCaloMET, "hltCaloMET_nofilt_MET_R207454", "%.0f")
project("hltCaloMET.pt", histos_hltCaloMET, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_hltCaloMET, "hltCaloMET_MET_R207454", "%.0f")
project("hltCaloMET.pt", histos_hltCaloMET, "metfilterFlags[%i] && hltPFMET.pt>100" % len(metfilters), True)
draw(histos_hltCaloMET, "hltCaloMET_if_hltPFMETgt100_MET_R207454", "%.0f")

histos_hltPFMET = book("hltPFMET", "HLT PF MET [GeV]", 100, 0, 200)
project("hltPFMET.pt", histos_hltPFMET, "1", True)
draw(histos_hltPFMET, "hltPFMET_nofilt_MET_R207454", "%.1f")
project("hltPFMET.pt", histos_hltPFMET, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_hltPFMET, "hltPFMET_MET_R207454", "%.1f")

# RECO MET
histos_recoPFMET = book("recoPFMET", "RECO PF MET [GeV]", 100, 0, 200)
project("recoPFMET.pt", histos_recoPFMET, "1", True)
draw(histos_recoPFMET, "recoPFMET_nofilt_MET_R207454", "%.1f")
project("recoPFMET.pt", histos_recoPFMET, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_recoPFMET, "recoPFMET_MET_R207454", "%.1f")

histos_recoPFMETT1 = book("recoPFMETT1", "RECO Type-1 PF MET [GeV]", 100, 0, 200)
project("recoPFMETT1.pt", histos_recoPFMETT1, "1", True)
draw(histos_recoPFMETT1, "recoPFMETT1_nofilt_MET_R207454", "%.1f")
project("recoPFMETT1.pt", histos_recoPFMETT1, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_recoPFMETT1, "recoPFMETT1_MET_R207454", "%.1f")

histos_recoPFMETT0T1 = book("recoPFMETT0T1", "RECO Type-0+1 PF MET [GeV]", 100, 0, 200)
project("recoPFMETT0T1.pt", histos_recoPFMETT0T1, "1", True)
draw(histos_recoPFMETT0T1, "recoPFMETT0T1_nofilt_MET_R207454", "%.1f")
project("recoPFMETT0T1.pt", histos_recoPFMETT0T1, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_recoPFMETT0T1, "recoPFMETT0T1_MET_R207454", "%.1f")

histos_patMPT = book("patMPT", "RECO TRK MET [GeV]", 100, 0, 200)
project("patMPT.pt", histos_patMPT, "1", True)
draw(histos_patMPT, "patMPT_nofilt_MET_R207454", "%.1f")
project("patMPT.pt", histos_patMPT, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_patMPT, "patMPT_MET_R207454", "%.1f")

histos_recoPFMETMVA = book("recoPFMETMVA", "RECO MVA PF MET [GeV]", 100, 0, 200)
project("recoPFMETMVA.pt", histos_recoPFMETMVA, "1", True)
draw(histos_recoPFMETMVA, "recoPFMETMVA_nofilt_MET_R207454", "%.1f")
project("recoPFMETMVA.pt", histos_recoPFMETMVA, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_recoPFMETMVA, "recoPFMETMVA_MET_R207454", "%.1f")

histos_recoPFMETNoPU = book("recoPFMETNoPU", "RECO NoPU PF MET [GeV]", 100, 0, 200)
project("recoPFMETNoPU.pt", histos_recoPFMETNoPU, "1", True)
draw(histos_recoPFMETNoPU, "recoPFMETNoPU_nofilt_MET_R207454", "%.1f")
project("recoPFMETNoPU.pt", histos_recoPFMETNoPU, "metfilterFlags[%i]" % len(metfilters), True)
draw(histos_recoPFMETNoPU, "recoPFMETNoPU_MET_R207454", "%.1f")