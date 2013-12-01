from drawer_init import *

drawerInit = DrawerInit()
wait = True
if not wait:  gROOT.SetBatch(1)

#tfile = TFile.Open("../bin/compactified.root")
tfile = TFile.Open("../bin/compactified.0.root")
#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
tree = tfile.Events


mets = [
    ("hltCaloMET", "HLT Calo MET [GeV]", "hltCaloMET.pt", 1),
    ("hltPFMET", "HLT PF MET [GeV]", "hltPFMET.pt", 1),
    ("recoPFMET", "RECO PF MET [GeV]", "recoPFMET.pt", kBlue),
    ("recoPFMETT1", "RECO Type-1 PF MET [GeV]", "recoPFMETT1.pt", kGreen),
    ("recoPFMETT0T1", "RECO Type-0+1 PF MET [GeV]", "recoPFMETT0T1.pt", kRed),
    ("patMPT", "RECO TRK MET [GeV]", "patMPT.pt", kYellow2),
    ("recoPFMETMVA", "RECO MVA PF MET [GeV]", "recoPFMETMVA.pt", kMagenta2),
    ("recoPFMETNoPU", "RECO NoPU PF MET [GeV]", "recoPFMETNoPU.pt", kCyan2),
]

def book(mets, nbins, xlow, xup, ylow, yup):
    histos = []
    for met in mets:
        histos.append(TProfile("p_" + met[0], "; # good PVs; " + met[1], nbins, xlow, xup, ylow, yup) )
        histos[-1].SetLineWidth(2)
        histos[-1].SetLineColor(met[3])
        histos[-1].SetMarkerColor(met[3])
    return histos

def project(mets, histos, addsel="1"):
    for i, met in enumerate(mets):
        tree.Project(histos[i].GetName(), met[2] + ":event.nGoodPV", runsel + " * triggerFlags[%i] * (%s)" % (i_reftrig, addsel), "prof goff")

def draw(histos, imgname, fnum="%.1f"):
    h = histos[0].Clone("slave")
    h.GetYaxis().SetTitle("<MET> [GeV]")
    h.Draw()

    histos[1].Draw("same")
    histos[2].Draw("same")
    histos[3].Draw("same")
    histos[4].Draw("same")

    leg = TLegend(0.20, 0.70, 0.6, 0.94)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.SetShadowColor(0)
    leg.AddEntry(histos[0], histos[0].GetYaxis().GetTitle()[:-6])
    leg.AddEntry(histos[1], histos[1].GetYaxis().GetTitle()[:-6])
    leg.AddEntry(histos[2], histos[2].GetYaxis().GetTitle()[:-6])
    leg.AddEntry(histos[3], histos[3].GetYaxis().GetTitle()[:-6])
    leg.AddEntry(histos[4], histos[4].GetYaxis().GetTitle()[:-6])
    leg.Draw()

    gPad.Print(imgdir + imgname + ".png")
    gPad.Print(imgdir + imgname + ".pdf")


mets1 = mets[1:6]
histos1 = book(mets1, 20, 0, 40, 0, 200)
project(mets1, histos1, "1")
draw(histos1, "pileup_prof_nofilt_MET_R207454")
project(mets1, histos1, "metfilterFlags[%i]" % len(metfilters))
draw(histos1, "pileup_prof_MET_R207454")
del mets1[:]

mets2 = [mets[1]] + mets[4:]
histos2 = book(mets2, 20, 0, 40, 0, 200)
project(mets2, histos2, "1")
draw(histos2, "pileup2_prof_nofilt_MET_R207454")
project(mets2, histos2, "metfilterFlags[%i]" % len(metfilters))
draw(histos2, "pileup2_prof_MET_R207454")

