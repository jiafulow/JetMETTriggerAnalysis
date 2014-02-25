#!/usr/bin/env python

from ROOT import TFile, TH1F, gROOT, gSystem, gInterpreter

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

#gSystem.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gInterpreter.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gROOT.LoadMacro("../src/SimpleCandidateLinkDef.h")
gROOT.LoadMacro("HelperFunctions.h")

#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
#tfile.tree = tfile.Events
tfile = TFile.Open("../bin/compactified.L1ETM40.4.root")


# puremet
# 9-3115-5  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
#sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"

# monojet
# 71-4917-23  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
#sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"

# susdijet
# 267-3094-32  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
#sel_trig1 = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"

# higdijet
# 30-1351-6  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) )
#sel_trig1 = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"

# multijet
# 55-923-1  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
#sel_trig1 = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))"

# bjet
# 21-1185-13  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) )
#sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"

# bjet service
# 6991-898-8  #< (only rereco)-(common)-(only real data)
# HLT_L1ETM40 PS=1750, HLT_DiCentralPFJet30_PFMET80 PS=150
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) )
#sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"

# VBF lead
# 7-831-27  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) )
#sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_leadmjj>600 && hltPFGlobal.vbf_leadmjj_deta>3.5 && hltPFMETNoMu.pt>65)"

# VBF max
# 5-614-41  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ800VBF_AllJets_v9')) )
#sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)>1 && hltPFGlobal.vbf_maxmjj>800 && hltPFGlobal.vbf_maxmjj_deta>3.5 && hltPFMETNoMu.pt>65)"


h1 = TH1F("h1", "h1", 5, 0, 5)
h2 = TH1F("h2", "h2", 5, 0, 5)

#tfile.tree.Project(h1.GetName(), sel_trig0, ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1, "goff")
#tfile.tree.Project(h2.GetName(), sel_trig1, ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "goff")

#tfile.tree.Draw("hltPFGlobal.vbf_leadmjj", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0)
#tfile.tree.Draw("hltPFGlobal.vbf_leadmjj", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltCaloJetsL1Fast[1].pt >> h(100,0,200)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltCaloJetsL1Fast[1].pt", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("hltPFJetsL1FastL2L3[1].pt >> h(100,0,200)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltPFJetsL1FastL2L3[1].pt", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30) >> h(10,0,10)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("Sum$(abs(hltCaloJetsL1Fast.eta)<5.0 && hltCaloJetsL1Fast.pt>30)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40) >> h(10,0,10)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("Sum$(abs(hltPFJetsL1FastL2L3.eta)<5.0 && hltPFJetsL1FastL2L3.pt>40)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50) >> h(10,0,10)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("hltPFHTMHTNoPU.pt >> h(100,0,500)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltPFHTMHTNoPU.pt", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("hltCaloHTMHT.pt >> h(100,0,500)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltCaloHTMHT.pt", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")
#tfile.tree.Draw("hltPFGlobal.vbf_maxmjj >> h(100,0,1000)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltPFGlobal.vbf_maxmjj", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")

h1.Print("all")
h2.Print("all")


sel_noNoise  = "(metfilterFlags[%i] && event.json)" % 8
sel_noNoise0 = "(metfilterFlags[%i] && event.json && (Sum$(patJets.pt>20)>0 && patJets[0].jetID==1 && ((Sum$(patJets.pt>20)>1 && patJets[1].jetID==1) || Sum$(patJets.pt>20)==1)) )" % 8
#selections = [
#    "(hltCaloMET.pt>80 && hltPFMET.pt>150)",
#    "(hltCaloMET.pt>90 && hltCaloMETClean.pt>-99 && hltCaloMETCleanUsingJetID.pt>-99 && hltPFMET.pt>150)",
#    "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>-99 && hltPFMET.pt>150)",
#    "(hltCaloMET.pt>90 && hltCaloMETClean.pt>80 && hltCaloMETCleanUsingJetID.pt>80 && hltPFMET.pt>150)",
#    ]
#sel_bench = "((Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)>1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && abs(patJets[1].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && abs(deltaPhi(patJets[0].phi,patJets[1].phi))<2.5) || (Sum$(patJets.pt>30 && abs(patJets.eta)<2.5)==1 && patJets[0].pt>110 && abs(patJets[0].eta)<2.5 && patJets[0].jetID==1))"
#sel_bench = sel_noNoise + "*" + sel_bench
sel_bench = "(Sum$(patJets.pt>50 && abs(patJets.eta)<2.5)>2 && patJets[0].pt>50 && abs(patJets[0].eta)<2.5 && patJets[1].pt>50 && abs(patJets[1].eta)<2.5 && patJets[2].pt>50 && abs(patJets[2].eta)<2.5 && patJets[0].jetID==1 && patJets[1].jetID==1 && patJets[2].jetID==1 && patHTMHT.sumEt>450)"
sel_bench = sel_noNoise + "*" + sel_bench

selections = [
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)",
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>85)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>100)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)",
    ]
selections = [
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)",
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>70 && hltCaloMETClean.pt>60 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>100||hltPFMET.pt>100) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>80 && hltCaloMETClean.pt>70 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>120||hltPFMET.pt>120) && hltPFJetsL1FastL2L3[0].nhf<0.95)",
    "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>70)>0 && hltCaloMET.pt>80 && hltCaloMETClean.pt>70 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>80)>0 && (hltPFMETNoMu.pt>120||hltPFMET.pt>120) && hltPFJetsL1FastL2L3[0].nhf<0.95 && hltPFJetsL1FastL2L3[0].nch>0)",
    ]
selections = [
    "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))",
    "(hltCaloHTMHT.sumEt>350 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100))",
    "(hltCaloHTMHT.sumEt>350 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>400 && (hltCaloHTMHT.pt>150 || hltPFMET.pt>100) && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)",
    "(hltCaloHTMHT.sumEt>350 && hltCaloMET.pt>65 && hltCaloMETClean.pt>55 && hltPFHTMHTNoPU.sumEt>400 && hltPFMET.pt>100 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>40)>1)",
    ]


tfile.tree.Draw("recoPFMETT0T1.pt", sel_bench + "*" + selections[1] + "* !" + selections[3])
#tfile.tree.Scan("recoPFMETT0T1.pt:patJets[0].pt", sel_bench + "*" + selections[2] + "* !" + selections[3])
