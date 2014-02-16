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


#_______________________________________________________________________________
# Configurations
fitterinit = FitterInit()
chain = TChain("tree", "tree")
infiles = [
    "../bin/compactified.L1ETM40.3.root",
    ]
for f in infiles:
    chain.Add(f)


sections = {}
sections["makedataset"]     = False
sections["makefits"]        = False
sections["makenewfits"]     = True
plotting = []
writing = []

#rootfilename = "dataset.root"
rootfilename = "dataset_PFMET150.root"

imgdir = "figures_20140206/"
if not imgdir.endswith("/"):  imgdir += "/"
if gSystem.AccessPathName(imgdir):
    gSystem.mkdir(imgdir)


#_______________________________________________________________________________
# Classes/Functions
class Properties:
    """Stores stuff"""
    def __init__(self, x, y):
        self.x = x
        self.y = y

def CMS_label():
    old = (latex.GetTextFont(), latex.GetTextSize())
    latex.SetTextFont(42); latex.SetTextSize(0.026)
    latex.DrawLatex(0.665, 0.968, "Run2012D HLT_L1ETM40_v2")
    latex.SetTextFont(62); latex.SetTextSize(0.028)
    latex.DrawLatex(0.445, 0.968, "CMS Preliminary")
    latex.SetTextFont(old[0]); latex.SetTextSize(old[1])
    return

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

def RunFactory(ws, eff, fitparams, var="var0", cats=["cat0"], xlow=0, xup=400):
    ws.factory("%s[%.2f,%.2f]" % (var,xlow,xup))
    for cat in cats:
        ws.factory("%s[reject=0,accept=1]" % cat)

    for i, f in enumerate(fitparams):
        ws.factory("fit%i[%.2f,%.2f,%.2f]" % (i, f[0], f[1], f[2]))

    ws.factory("expr::x('(%s-fit0)/fit1', %s, fit0, fit1)" % (var,var))
    ws.factory("expr::effLogistic('fit2 / (1+exp(-x*2))', x, fit2)")
    ws.factory("expr::effTanh('fit2/2 * (1 + TMath::TanH(x))', x, fit2)")
    ws.factory("expr::effArctan('fit2/2 * (1 + (2/TMath::Pi())*TMath::ATan(x))', x, fit2)")
    #ws.factory("expr::effArctan('fit2/2 * (1 + (2/TMath::Pi())*TMath::ATan(exp(TMath::Pi()/2 * x)))', x, fit2)")
    ws.factory("expr::effAlgebra('fit2/2 * (1 + x/(TMath::Sqrt(1+TMath::Power(x,2))))', x, x, fit2)")
    # NOTE: effError won't converge properly if eff=0 expected but eff>0 observed. Need to truncate data range
    ws.factory("expr::effError('fit2/2 * (1 + TMath::Erf(x/TMath::Sqrt2()))', x, fit2)")
    for cat in cats:
        ws.factory("Efficiency::%s(%s, %s, 'accept')" % (cat.replace("cat","model"), eff, cat))
    #ws.Print()

def draw(ws, res, eff, var, cat, xtitle, trig, nfitparams, xlow=0, xup=400, ymax=1.1, new=False):
    # Set variable-width binning
    tmpbins = [0,10,20,30,40,50,60,70,80,90,100,110,120,140,160,180,200,220,240,280,320,360,400,480,600,750,1000]
    tmpbins2 = [xlow] + [x for x in tmpbins if (x > xlow and x < xup)] + [xup]
    bins = RooBinning(len(tmpbins2)-1, numpy.array(tmpbins2, dtype=float))

    kColor = 4  # kBlue
    kColorErr = 432  # kCyan
    if new:
        kColor = 618  # kMagenta+2
        kColorErr = 609  # kMagenta-7
    frame1 = ws.var(var).frame(RooFit.Bins(40), RooFit.Title("Data (all, accepted)"))
    ds.plotOn(frame1)
    ds.plotOn(frame1, RooFit.Cut("%s==%s::accept" % (cat,cat)), RooFit.MarkerColor(kColor), RooFit.LineColor(kColor))

    frame2 = ws.var(var).frame(RooFit.Title("; %s; HLT efficiency" % xtitle))
    ds.plotOn(frame2, RooFit.Binning(bins), RooFit.Efficiency(ws.cat(cat)))
    ws.function(eff).plotOn(frame2, RooFit.LineColor(kColor))
    nchisq = frame2.chiSquare(nfitparams)
    ndof = bins.numBins()-nfitparams
    chisq = nchisq * ndof
    ws.function(eff).plotOn(frame2, RooFit.VisualizeError(res), RooFit.FillColor(kColorErr), RooFit.MoveToBack())

    gPad.SetLeftMargin(0.15); frame2.GetYaxis().SetTitleOffset(1.25); frame2.GetYaxis().SetRangeUser(0,ymax)
    #frame1.Draw()
    frame2.Draw()
    CMS_label()

    latex.SetTextSize(0.044)
    latex.DrawLatex(0.60, 0.30, "#chi^{2}/ndof = %4.2f/%.0f" % (chisq,ndof))
    text = ["** "+trig+" **"]
    text.append("chi2/ndof      %6.2f  =  %4.2f/%.0f" % (nchisq, chisq, ndof))
    print "chi2/ndof      %6.2f  =  %4.2f/%.0f" % (nchisq, chisq, ndof)

    latex.SetTextSize(0.026)
    latex.DrawLatex(0.60, 0.26, trig)
    for i in xrange(3):
        v = ws.var("fit%i" % i)
        if v.getVal() > 1:
            text.append("p_%i            %6.2f  +/- %3.2f" % (i, v.getVal(), v.getError()))
            latex.DrawLatex(0.60, 0.26-(i+1)*0.03, "p_{%i}            %6.2f  +/- %3.2f" % (i, v.getVal(), v.getError()))
        else:
            text.append("p_%i             %6.3f +/- %4.3f" % (i, v.getVal(), v.getError()))
            latex.DrawLatex(0.60, 0.26-(i+1)*0.03, "p_{%i}             %6.3f +/- %4.3f" % (i, v.getVal(), v.getError()))
    line.DrawLine(xlow, 1, xup, 1)
    text.append("")
    return text

def save(imgdir, imgname):
    gPad.RedrawAxis()
    gPad.Print(imgdir+imgname+".pdf")
    gPad.Print(imgdir+imgname+".png")


#_______________________________________________________________________________
# Make dataset

if sections["makedataset"]:
    sel = "(triggerFlags[%i])" % ireftrig
    sel_noNoise  = "(metfilterFlags[%i] && event.json)" % len(metfilters)
    sel_noNoise0 = "(metfilterFlags[%i] && event.json && (Sum$(patJets.pt>20)>0 && patJets[0].jetID==1 && ((Sum$(patJets.pt>20)>1 && patJets[1].jetID==1) || Sum$(patJets.pt>20)==1)) )" % len(metfilters)
    sel_noNoise1 = "(metfilterFlags[%i] && event.json && Sum$(patJets.pt>20)>0 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1)" % len(metfilters)
    sel_noNoise2 = "(metfilterFlags[%i] && event.json && Sum$(patJets.pt>20 && abs(patJets.eta)<2.5)>1 && patJets[0].jetID==1 && ((abs(patJets[1].eta)<2.5 && patJets[1].jetID==1)||abs(patJets[1].eta)>=2.5) )" % len(metfilters)
    selection = sel + "*" + sel_noNoise0
    categories = [
        ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) ),
        ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) ),
        ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) ),
        ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) ),
        #("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4'), triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5'))),
        ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) ),
        ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) ),
        #("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) ),
        #("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5'), triggers.index('HLT_DiCentralPFJet30_PFMET80_v6'))),
        ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) ),
        ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) ),
        #("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9'), triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9'))),
        # dummy x2
        ("(!triggerFlags[%i])" % ireftrig),
        ("(!triggerFlags[%i])" % ireftrig),

        # proposal
        "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)",
        ]
    variables = [
        "recoPFMETT0T1.pt",
        "patMPT.pt",
        "hltPFMET.pt",
        "hltCaloMET.pt",
        "hltCaloMETClean.pt",
        "hltTrackMET.pt",
        ]
    MakeTree(chain, selection, categories, variables, outfilename=rootfilename)


#_______________________________________________________________________________
# Make fits
if sections["makefits"]:
    if plotting: del plotting[:]
    if writing: del writing[:]

    eff = "effError"
    #eff = "effLogistic"
    #eff = "effTanh"
    #eff = "effArctan"
    #eff = "effAlgebra"

    xlow, xup = 80, 400
    cat = "cat0"
    kTrig = "PFMET150"
    fitparams = [
        (170,  80, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.92, 0.3, 1.0),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 50, 360
    cat = "cat1"
    kTrig = "MonoCentralJet"
    fitparams = [
        (130,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.9, 0.2, 1.0),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat2"
    kTrig = "DiCentralJetSUS"
    fitparams = [
        (120,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.5, 0.2, 0.7),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat3"
    kTrig = "DiCentralJetHIG"
    fitparams = [
        (120,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.5, 0.2, 0.7),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat4"
    kTrig = "HT"
    fitparams = [
        (140,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.4, 0.2, 0.7),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat5"
    kTrig = "btag"
    fitparams = [
        (110,  40, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.2,0.05,0.25),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 30, 240
    cat = "cat6"
    kTrig = "VBFAll"
    fitparams = [
        ( 100,  40, 220),  # mean
        (  20,   5,  60),  # sigma
        (0.02,0.01,0.04),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))

    xlow, xup = 30, 240
    cat = "cat7"
    kTrig = "VBFLead"
    fitparams = [
        ( 100,  40, 220),  # mean
        (  20,   5,  60),  # sigma
        (0.02,0.01,0.04),  # plateau
        ]
    plotting.append((cat, kTrig, fitparams, xlow, xup))


    # Make fits and plots
    c1 = TCanvas("c1", "c1")
    for p in plotting:
        (cat, kTrig, fitparams, xlow, xup) = p
        nfitparams = len(fitparams)

        var = "var0"
        xtitle = "#scale[0.7]{RECO} T0T1 PFMET [GeV]"
        ws = RooWorkspace("workspace", "workspace")
        RunFactory(ws, eff, fitparams, var, [cat], xlow=xlow, xup=xup)
        #ws.var("fit0").setConstant(1); ws.var("fit2").setConstant(1); nfitparams -= 2
        #wsfilename = "workspace.root"
        #ws.writeToFile(wsfilename)

        dsfile = TFile.Open(rootfilename)
        ds = RooDataSet("data", "data", RooArgSet(ws.var(var), ws.cat(cat)), RooFit.Import(dsfile.tree)); ds.Print()

        model = cat.replace("cat","model")
        res = ws.pdf(model).fitTo(ds, RooFit.ConditionalObservables(RooArgSet(ws.var(var))), RooFit.Minimizer("Minuit"), RooFit.Save());

        ymax = fitparams[2][2]+0.1
        if ymax < 0.2:  ymax = fitparams[2][2]+0.01
        text = draw(ws, res, eff, var, cat, xtitle, kTrig, nfitparams, xlow=xlow, xup=xup, ymax=ymax)

        for t in text:  writing.append(t)
        save(imgdir, "fiteff_" + kTrig + "_" + eff)

    # Printout
    print "-" * 40
    for t in writing:  print t


#_______________________________________________________________________________
# Make new fits
if sections["makenewfits"]:
    if plotting: del plotting[:]
    if writing: del writing[:]

    eff = "effError"
    #eff = "effLogistic"
    #eff = "effTanh"
    #eff = "effArctan"
    #eff = "effAlgebra"

    xlow, xup = 80, 400
    cat = "cat0"
    newcat = "cat10"
    kTrig = "PFMET150"
    fitparams = [
        (170,  80, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.92, 0.3, 1.0),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 50, 360
    cat = "cat1"
    newcat = "cat11"
    kTrig = "MonoCentralJet"
    fitparams = [
        (130,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.9, 0.2, 1.0),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat2"
    newcat = "cat12"
    kTrig = "DiCentralJetSUS"
    fitparams = [
        (120,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.5, 0.2, 0.7),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat3"
    newcat = "cat13"
    kTrig = "DiCentralJetHIG"
    fitparams = [
        (120,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.5, 0.2, 0.7),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat4"
    newcat = "cat14"
    kTrig = "HT"
    fitparams = [
        (140,  50, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.4, 0.2, 0.7),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 40, 320
    cat = "cat5"
    newcat = "cat15"
    kTrig = "btag"
    fitparams = [
        (110,  40, 250),  # mean
        ( 20,   1,  50),  # sigma
        (0.2,0.05,0.25),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 30, 240
    cat = "cat6"
    newcat = "cat16"
    kTrig = "VBFAll"
    fitparams = [
        ( 100,  40, 220),  # mean
        (  20,   5,  60),  # sigma
        (0.02,0.01,0.04),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))

    xlow, xup = 30, 240
    cat = "cat7"
    newcat = "cat17"
    kTrig = "VBFLead"
    fitparams = [
        ( 100,  40, 220),  # mean
        (  20,   5,  60),  # sigma
        (0.02,0.01,0.04),  # plateau
        ]
    plotting.append((cat, newcat, kTrig, fitparams, xlow, xup))


    # Make fits and plots
    c1 = TCanvas("c1", "c1")
    for p in plotting[:1]:
        (cat, newcat, kTrig, fitparams, xlow, xup) = p
        nfitparams = len(fitparams)

        var = "var0"
        xtitle = "#scale[0.7]{RECO} T0T1 PFMET [GeV]"
        ws = RooWorkspace("workspace", "workspace")
        RunFactory(ws, eff, fitparams, var, [cat, newcat], xlow=xlow, xup=xup)
        #ws.var("fit0").setConstant(1); ws.var("fit2").setConstant(1); nfitparams -= 2
        #wsfilename = "workspace.root"
        #ws.writeToFile(wsfilename)

        dsfile = TFile.Open(rootfilename)
        ds = RooDataSet("data", "data", RooArgSet(ws.var(var), ws.cat(cat), ws.cat(newcat)), RooFit.Import(dsfile.tree)); ds.Print()

        model = cat.replace("cat","model")
        res = ws.pdf(model).fitTo(ds, RooFit.ConditionalObservables(RooArgSet(ws.var(var))), RooFit.Minimizer("Minuit"), RooFit.Save());

        ymax = fitparams[2][2]+0.1
        if ymax < 0.2:  ymax = fitparams[2][2]+0.01
        text = draw(ws, res, eff, var, cat, xtitle, kTrig, nfitparams, xlow=xlow, xup=xup, ymax=ymax)
        for t in text:  writing.append(t)

        text = draw(ws, res, eff, var, newcat, xtitle, kTrig, nfitparams, xlow=xlow, xup=xup, ymax=ymax, new=True)
        for t in text:  writing.append(t)

        save(imgdir, "fitneweff_" + kTrig + "_" + eff)

    # Printout
    print "-" * 40
    for t in writing:  print t

