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


#_______________________________________________________________________________

#gSystem.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gInterpreter.AddIncludePath("/uscms_data/d2/jiafu/Trigger/CMSSW_5_3_11/src/")
gROOT.LoadMacro("../src/SimpleCandidateLinkDef.h")

#tfile = TFile.Open("../bin/compactified.L1ETM40.0.root")
#tfile.tree = tfile.Events
tfile = TFile.Open("../bin/compactified_13.root")


# puremet
# 9-3109-5  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFMET150_v7'), triggers.index('HLT_PFMET180_v7')) )
#sel_trig1 = "(hltCaloMET.pt>80 && hltPFMET.pt>150)"

# monojet
# 111-4907-23  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4')) )
#sel_trig1 = "(Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>65)>0 && hltCaloMET.pt>65 && hltPFJetsL1FastL2L3[0].nhf<0.95 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>80)>0 && hltPFMETNoMu.pt>105)"

# susdijet
# XX-YY-ZZ  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFNoPUJet50_PFMETORPFMETNoMu80_v4')) )
#sel_trig1 = "(hltCaloMET.pt>80 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3NoPU.eta)<2.6 && hltPFJetsL1FastL2L3NoPU.pt>50)>1 && (hltPFMET.pt>80||hltPFMETNoMu.pt>80) )"

# higdijet
# XX-YY-ZZ  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralJetSumpT100_dPhi05_DiCentralPFJet60_25_PFMET100_HBHENoiseCleaned_v5')) )
#sel_trig1 = "(hltCaloMET.pt>50 && hltCaloMETClean.pt>25 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>15)>1 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>50)>0 && hltCaloGlobal.dijet_maxpt>100 && hltCaloGlobal.dijet_mindphi_j40>0.5 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>25)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>60)>0 && hltPFMET.pt>100)"

# multijet
# XX-YY-ZZ  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i] || triggerFlags[%i])" %(triggers.index('HLT_PFNoPUHT350_PFMET100_v4'), triggers.index('HLT_PFNoPUHT400_PFMET100_v4')) )
#sel_trig1 = "(hltCaloHTMHT.sumEt>300 && hltCaloHTMHT.pt>75 && hltPFHTMHTNoPU.sumEt>350 && (hltPFHTMHTNoPU.pt>150 || hltPFMET.pt>100))"

# bjet
# XX-YY-ZZ  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_BTagCSV07_v5')) )
#sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && hltCaloGlobal.bjet_maxcsv>0.7 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"

# bjet service
# XX-YY-ZZ  #< (only rereco)-(common)-(only real data)
#sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiCentralPFJet30_PFMET80_v6')) )
#sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>20)>1 && Sum$(abs(hltPFJetsL1FastL2L3.eta)<2.6 && hltPFJetsL1FastL2L3.pt>30)>1 && hltPFMET.pt>80)"

# VBF lead
# XX-YY-ZZ  #< (only rereco)-(common)-(only real data)
sel_trig0 = ("(triggerFlags[%i])" %(triggers.index('HLT_DiPFJet40_PFMETnoMu65_MJJ600VBF_LeadingJets_v9')) )
sel_trig1 = "(hltCaloMET.pt>65 && Sum$(abs(hltCaloJetsL1Fast.eta)<2.6 && hltCaloJetsL1Fast.pt>30)>1 && hltPFGlobal.vbf_leadmjj>600 && hltPFGlobal.vbf_leadmjj_deta>3.5 && hltPFMETNoMu.pt>65)"


h1 = TH1F("h1", "h1", 5, 0, 5)
h2 = TH1F("h2", "h2", 5, 0, 5)

tfile.tree.Project(h1.GetName(), sel_trig0, ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1, "goff")
tfile.tree.Project(h2.GetName(), sel_trig1, ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "goff")

#tfile.tree.Draw("hltPFGlobal.vbf_leadmjj", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0)
#tfile.tree.Draw("hltPFGlobal.vbf_leadmjj", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltPFJetsL1FastL2L3[1].pt >> h(100,0,500)", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig1)
#tfile.tree.Draw("hltPFJetsL1FastL2L3[1].pt", ("(triggerFlags[%i]) * " % ireftrig) + sel_trig0, "same")

h1.Print("all")
h2.Print("all")