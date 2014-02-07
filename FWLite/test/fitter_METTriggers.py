#!/usr/bin/env python

from ROOT import TFile, TH1, TH1F, TChain, TTree, TTreeFormula, TCanvas, TLatex, TLine, gPad, gROOT, gSystem, gInterpreter, gStyle, kGray, RooFit, RooWorkspace, RooRealVar, RooDataSet, RooArgSet, RooArgList, RooFormulaVar, RooCategory, RooEfficiency, RooBinning
import numpy


# For init
class FitterInit:
    def __init__(self):
        # ROOT
        gROOT.LoadMacro("tdrstyle.C")
        #gROOT.LoadMacro("HelperFunctions.h")
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
latex.SetTextSize(0.044)

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


#_______________________________________________________________________________
# Configurations
fitterinit = FitterInit()
chain = TChain("tree", "tree")
infiles = [
    "../bin/compactified.L1ETM40.1.root"
    ]
for f in infiles:
    chain.Add(f)

imgdir = "figures_20140206/"
if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)

wait = True
if not wait:  gROOT.SetBatch(1)


#_______________________________________________________________________________
# Classes/Functions
class Properties:
    """Stores stuff"""
    def __init__(self, x, y):
        self.x = x
        self.y = y

def MakeTree(chain, selection, categories, variables, outtreename="tree", outfilename="dataset.root"):

    outfile = TFile.Open(outfilename, "RECREATE")
    outtree = TTree(outtreename, outtreename)

    ttfsel = TTreeFormula("selection", selection, chain)  #preselection
    outcats = []
    ttfcats = []
    for i, t in enumerate(categories):
        b = numpy.zeros(1, dtype=int)
        outcats.append(b)
        outtree.Branch("cat%i" % i, b, 'cat%i/I' % i)
        ttfcats.append(TTreeFormula("cat%i" % i, t, chain))
    outvars = []
    ttfvars = []
    for i, t in enumerate(variables):
        b = numpy.zeros(1, dtype=float)
        outvars.append(b)
        outtree.Branch("var%i" % i, b, 'var%i/D' % i)
        ttfvars.append(TTreeFormula("var%i" % i, t, chain))

    nevt = chain.GetEntries()
    tnumber = chain.GetTreeNumber()

    for ievt in xrange(0, nevt):
        localEntry = chain.LoadTree(ievt)  # used by TTreeFormula
        if localEntry < 0:  break

        if tnumber != chain.GetTreeNumber():
            tnumber = chain.GetTreeNumber()
            ttfsel.UpdateFormulaLeaves()
            for t in ttfcats:
                t.UpdateFormulaLeaves()
            for t in ttfvars:
                t.UpdateFormulaLeaves()

        ndata = ttfsel.GetNdata()
        keep = False
        for idata in xrange(ndata):
            keep |= (bool(ttfsel.EvalInstance(idata)) != 0)
            if keep:  break
        if not keep:  continue
        #keep = bool(ttfsel.EvalInstance())
        #if not keep:  continue

        for i, t in enumerate(ttfcats):
            t.GetNdata()
            outcats[i][0] = int(t.EvalInstance())
        for i, t in enumerate(ttfvars):
            t.GetNdata()
            outvars[i][0] = float(t.EvalInstance())

        #chain.GetEntry(ievt)
        outtree.Fill()

    outtree.Write()
    outfile.Close()

# FIXME: allow init values as input
def RunFactory(ws, eff, fitparams, var="var0", cat="cat0", xlow=0, xup=400):

    ws.factory("var0[%.2f,%.2f]" % (xlow,xup))
    ws.factory("var1[%.2f,%.2f]" % (xlow,xup))
    ws.factory("var2[%.2f,%.2f]" % (xlow,xup))
    ws.factory("var3[%.2f,%.2f]" % (xlow,xup))
    ws.factory("cat0[reject=0,accept=1]")

    for i, f in enumerate(fitparams):
        ws.factory("fit%i[%.2f,%.2f,%.2f]" % (i, f[0], f[1], f[2]))

    ws.factory("expr::x('(%s-fit0)/fit1', %s, fit0, fit1)" % (var,var))
    ws.factory("expr::effLogistic('fit2 / (1+exp(-x))', x, fit2)")
    ws.factory("expr::effTanh('fit2/2 * (1 + TMath::TanH(x/2))', x, fit2)")
    ws.factory("expr::effArctan('fit2/2 * (1 + (2/TMath::Pi())*TMath::ATan(x))', x, fit2)")
    #ws.factory("expr::effArctan('fit2/2 * (1 + (2/TMath::Pi())*TMath::ATan(exp(TMath::Pi()/2 * x)))', x, fit2)")
    ws.factory("expr::effAlgebra('fit2/2 * (1 + x/(TMath::Sqrt(1+TMath::Power(x,2))))', x, x, fit2)")
    # NOTE: effError won't converge properly if eff=0 expected but eff>0 observed. Need to truncate data range
    ws.factory("expr::effError('fit2/2 * (1 + TMath::Erf(x/TMath::Sqrt2()))', x, fit2)")
    ws.factory("Efficiency::model(%s, %s, 'accept')" % (eff, cat))
    ws.Print()

def save(imgdir, imgname):
    gPad.Print(imgdir+imgname+".pdf")
    gPad.Print(imgdir+imgname+".png")


#_______________________________________________________________________________
# Make dataset

selection = "(triggerFlags[%i] && metfilterFlags[%i])" % (ireftrig, len(metfilters))
categories = [
    ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) ),
    ]
variables = [
    "recoPFMETT0T1.pt",
    "hltPFMET.pt",
    "hltCaloMET.pt",
    "hltCaloMETClean.pt",
    ]
rootfilename = "dataset.root"
#MakeTree(chain, selection, categories, variables, outfilename=rootfilename)


#_______________________________________________________________________________
# Fit results


#eff = "effLogistic"
#eff = "effTanh"
#eff = "effArctan"
#eff = "effAlgebra"
eff = "effError"
fitparams = [
    # Used in effLogistic, effTanh, effArctan, effAlgebra, effError
    (170,  50, 250),  # mean
    ( 20,   1,  50),  # sigma
    (0.9, 0.3, 1.0),  # plateau
    ]
nfitparams = len(fitparams)

#xlow, xup = 50, 400
xlow, xup = 80, 400
var = "var2"
cat = "cat0"
ws = RooWorkspace("workspace", "workspace")
RunFactory(ws, eff, fitparams, var, cat, xlow=xlow, xup=xup)

dsfile = TFile.Open(rootfilename)
ds = RooDataSet("data", "data", RooArgSet(ws.var(var), ws.cat(cat)), RooFit.Import(dsfile.tree)); ds.Print()
res = ws.pdf("model").fitTo(ds, RooFit.ConditionalObservables(RooArgSet(ws.var(var))), RooFit.Minimizer("Minuit"), RooFit.Save()); res.Print()

wsfilename = "workspace.root"
ws.writeToFile(wsfilename)


#_______________________________________________________________________________
# Make plots

if xlow > (120-10):  print "ERROR: This binning won't work!"
bins = RooBinning(xlow, xup)
bins.addUniform((120-xlow)/10,xlow,120)
bins.addUniform(6,120,240)
bins.addUniform(2,240,320)
bins.addUniform(1,320,xup)

wsfile = TFile.Open(wsfilename)
ws = wsfile.workspace

kColor = 4
frame1 = ws.var(var).frame(RooFit.Bins(40), RooFit.Title("Data (all, accepted)"))
ds.plotOn(frame1)
ds.plotOn(frame1, RooFit.Cut("%s==%s::accept" % (cat,cat)), RooFit.MarkerColor(kColor), RooFit.LineColor(kColor))

#frame2 = ws.var(var).frame(RooFit.Bins(40), RooFit.Title("Fitted efficiency"))
frame2 = ws.var(var).frame(RooFit.Title("; #scale[0.7]{RECO} Type-0+1 PFMET [GeV]; HLT efficiency"))
ds.plotOn(frame2, RooFit.Binning(bins), RooFit.Efficiency(ws.cat(cat)))
ws.function(eff).plotOn(frame2, RooFit.LineColor(kColor))

#c1 = TCanvas("c1", "c1", 800, 400)
#c1.Divide(2);
#c1.cd(1); gPad.SetLeftMargin(0.15); frame1.GetYaxis().SetTitleOffset(1.6)
#frame1.Draw()
#c1.cd(2); gPad.SetLeftMargin(0.15); frame2.GetYaxis().SetTitleOffset(1.4)
#frame2.Draw()

c1 = TCanvas("c1", "c1")
gPad.SetLeftMargin(0.15); frame2.GetYaxis().SetTitleOffset(1.25); frame2.GetYaxis().SetRangeUser(0,1.1)
frame2.Draw()

chisq = frame2.chiSquare(nfitparams)
print "chi2/ndof = %.2f" % chisq
latex.DrawLatex(0.68, 0.22, "#chi^{2}/dof = %.2f" % chisq)
line.DrawLine(xlow, 1, xup, 1)
save(imgdir, "fiteff_PFMET150")
