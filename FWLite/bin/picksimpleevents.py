from ROOT import TFile, TTreeFormula, gInterpreter, gROOT
gInterpreter.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gROOT.LoadMacro("../src/SimpleCandidateLinkDef.h")
import warnings
warnings.filterwarnings( action='ignore', category=RuntimeWarning, message='creating converter.*' )

ibit = 0

def dimmi(kTrig, selection, rootfile="compactified.PFMETplusX.0.root"):
    global ibit
    textfile = "pickevents_%s.log" % kTrig
    pickevents = []
    with open(textfile) as f:
        for line in f.readlines():
            if ":" in line:
                run, lumi, event = line.strip().split(":")
                pickevents.append((long(run), long(lumi), long(event)))
    if len(pickevents) != 200:
        print "Expected 200 events!"

    tree = TFile.Open(rootfile).tree
    ttf = TTreeFormula("ttf", selection, tree)
    passed, failed = 0, 0
    for event in tree:
        run, lumi, event = long(tree.event.run), long(tree.event.lumi), long(tree.event.event)
        if (run, lumi, event) in pickevents:
            ttf.GetNdata()
            selected = bool(ttf.EvalInstance())
            if selected:
                passed += 1
            else:
                failed += 1

    if ibit == 0:
        print "TrigReport  Trig Bit#        Run     Passed     Failed      Error Name"
    print "TrigReport     1  %3i        %3i        %3i        %3i          0 %s  " % (ibit, passed+failed, passed, failed, kTrig)
    ibit += 1


if __name__ == "__main__":

    kTrig = "PFMET150"
    selection = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)"
    dimmi(kTrig, selection)

    kTrig = "PFMET150"
    selection = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>180)"
    dimmi(kTrig, selection)

    kTrig = "MonoCentralJet"
    selection = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)"
    dimmi(kTrig, selection)

    kTrig = "DiCentralJetHIG"
    selection = "(hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>0 && (hltPFMET.pt>100||hltPFMETNoMu.pt>100) && hltPFGlobal.dijet_mindphi_2cj>0.5)"
    dimmi(kTrig, selection)

    kTrig = "DiCentralJetSUS"
    selection = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>40)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>60)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) && hltPFGlobal.dijet_mindphi_2cj>0.3)"
    dimmi(kTrig, selection)

    kTrig = "HTMET"
    selection = "(hltCaloHTMHT.sumEt>350 && hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && hltPFHTMHTNoPU.sumEt>400 && (hltPFMET.pt>80||hltPFMET.pt>80) && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)"
    dimmi(kTrig, selection)

    kTrig = "btagMET"
    selection = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && hltPFMET.pt>80)"
    dimmi(kTrig, selection)

    kTrig = "btagMET"
    selection = "(hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>-9999. && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>30)>1 && hltPFMET.pt>80)"
    dimmi(kTrig, selection)

    kTrig = "PFMET150"
    selection = "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltTrackMET.pt>90)"
    dimmi(kTrig, selection)
