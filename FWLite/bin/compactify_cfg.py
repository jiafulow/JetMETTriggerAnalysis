import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
import re


# Parse command line arguments
options = VarParsing('analysis')
options.register('njobs',
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Number of jobs (default: 1)."
)
options.register('jobid',
    0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Job ID (default: 0)."
)
options.register('isData',
    True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run on data or Monte Carlo simulation (default: True)."
)
options.register('verbose',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Enable verbose outputs (default: False)"
)
options.parseArguments()
if options.jobid < 0:  options.jobid = 0
if options.jobid >= options.njobs:  options.jobid = options.njobs - 1


# Get JSON file correctly parsed
import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
JSONfile = 'Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt'
lumilist = LumiList.LumiList(filename = JSONfile).getCMSSWString().split(',')

process = cms.Process("p")
process.input = cms.PSet(
    #fileNames   = cms.vstring("../../MyProducts.MET.50k.root"),
    #fileNames   = cms.vstring("../../../HLTrigger/HLTanalyzers/test/openHLT/MyProducts.MET.root"),
    fileNames   = cms.vstring("../../../HLTrigger/HLTanalyzers/test/openHLT/MyProducts.MET.root"),
    #fileNames   = cms.vstring("/eos/uscms/store/user/jiafu/METTriggers/skimHLTPFMET150_eos_20130909/TT_CT10_TuneZ2star_8TeV-powheg-tauola-PU25bx50_2/MyProducts.MC_2_1_ukP.root"),
    njobs       = cms.int32(options.njobs),
    jobid       = cms.int32(options.jobid),
    maxEvents   = cms.int32(options.maxEvents),
    runMin      = cms.int32(-1),   ## optional
    runMax      = cms.int32(-1),   ## optional
    skipEvents  = cms.int32(0),    ## optional
    reportEvery = cms.int32(1000), ## optional
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange(lumilist)),
)

if True:
    import os
    #dirname = "/pnfs/cms/WAX/resilient/jiafulow/METTriggers/skimHLTPFMET150_20130909/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-PU25bx50/"
    #dirname = "/eos/uscms/store/user/jiafu/METTriggers/skimHLTPFMET150_eos_20130909/TT_CT10_TuneZ2star_8TeV-powheg-tauola-PU25bx50_2/"
    #dirname = "/pnfs/cms/WAX/resilient/jiafulow/METTriggers/skimHLTPFMET150_20130909/MET2012D/"
    #dirname = "/pnfs/cms/WAX/resilient/jiafulow/METTriggers/skimHLTL1ETM40_20131007/MET-Run2012D-PromptReco-v1/"
    #dirname = "/eos/uscms/store/user/jiafu/METTriggers/skimHLTPFMETplusX_20131104/MET-Run2012D-PromptReco-v1/"
    dirname = "/eos/uscms/store/user/jiafu/METTriggers/openHLT2PAT_20131128/MET-Run2012D-skimHLTL1ETM40/"
    if not dirname.endswith("/"):  dirname += "/"
    basenamelist = sorted(os.listdir(dirname))
    if dirname.startswith("/pnfs/cms/"):
        dirname = "dcache:" + dirname
    if dirname.startswith("/eos/uscms"):  # some files have size = 0
        basenamelist = [basename for basename in basenamelist if os.path.getsize(dirname + basename) > 1000]
    process.input.fileNames = cms.vstring()
    #process.input.fileNames.extend([dirname + basename for basename in basenamelist[:100]])
    process.input.fileNames.extend([dirname + basename for basename in basenamelist])


process.output = cms.PSet(
    fileName = cms.string("compactified.root"),
)

process.analyzer = cms.PSet(
    hltJetID = cms.PSet(
        hltCaloJetID = cms.PSet(
            min_N90 = cms.int32( -2 ),
            min_N90hits = cms.int32( 2 ),
            min_EMF = cms.double( 1.0E-6 ),
            max_EMF = cms.double( 999.0 ),
            JetIDParams = cms.PSet(
              useRecHits = cms.bool( True ),
              hbheRecHitsColl = cms.InputTag( "hltHbhereco" ),
              hoRecHitsColl = cms.InputTag( "hltHoreco" ),
              hfRecHitsColl = cms.InputTag( "hltHfreco" ),
              ebRecHitsColl = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEB' ),
              eeRecHitsColl = cms.InputTag( 'hltEcalRecHitAll','EcalRecHitsEE' )
            ),
        ),
        hltPFJetID = cms.PSet(
            min_CEEF = cms.double( -99.0 ),
            max_CEEF = cms.double( 99.0 ),
            min_NEEF = cms.double( -99.0 ),
            max_NEEF = cms.double( 99.0 ),
            min_CHEF = cms.double( -99.0 ),
            max_CHEF = cms.double( 99.0 ),
            min_NHEF = cms.double( -99.0 ),
            max_NHEF = cms.double( 0.95 ),
            triggerType = cms.int32( 85 ),
            nJet = cms.uint32( 1 )
        ),
    ),
    jetID = cms.PSet(
        caloJetID = cms.PSet(
            version = cms.string('PURE09'),
            quality = cms.string('LOOSE'),
        ),
        pfJetID = cms.PSet(
            version = cms.string('FIRSTDATA'),
            quality = cms.string('LOOSE'),
        ),
    ),
    triggers = cms.vstring(
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
    ),
    metfilters = cms.vstring(
        "p_HBHENoiseFilter",
        "p_CSCTightHaloFilter",
        "p_hcalLaserEventFilter",
        "p_EcalDeadCellTriggerPrimitiveFilter",
        "p_trackingFailureFilter",
        "p_eeBadScFilter",
        "p_ecalLaserCorrFilter",
        "p_trkPOGFilters",
    ),
    optmetfilters = cms.vstring(
        "p_goodVerticesFilter",
        "p_noscraping",
        "p_hcallasereventfilter2012",
        "p_EcalDeadCellBoundaryEnergyFilter",
        "p_tobtecfakesFilters",
    ),
    hltCaloJetPtMin = cms.double(10),
    hltCaloJetEtaMax = cms.double(5),
    hltPFJetPtMin = cms.double(15),
    hltPFJetEtaMax = cms.double(5),
    patJetPtMin = cms.double(15),
    patJetEtaMax = cms.double(5),
    isData = cms.bool(options.isData),
    isGolden = cms.bool(False),
    verbose = cms.bool(options.verbose),
)

process.handler = cms.PSet(
    # L1
    l1METs = cms.string("hltL1extraParticles:MET"),
    l1MHTs = cms.string("hltL1extraParticles:MHT"),
    l1CenJets = cms.string("hltL1extraParticles:Central"),
    l1ForJets = cms.string("hltL1extraParticles:Forward"),
    l1TauJets = cms.string("hltL1extraParticles:Tau"),
    l1NoIsoEGs = cms.string("hltL1extraParticles:NonIsolated"),
    l1IsoEGs = cms.string("hltL1extraParticles:Isolated"),
    l1Mus = cms.string("hltL1extraParticles"),
    l1HFRings = cms.string("hltL1extraParticles"),
    # HLT
    hltCaloJets = cms.string("hltAntiKT5CaloJets"),
    #hltCaloJets = cms.string("hltAntiKT5CaloJetsPF"),
    hltCaloJetsIDPassed = cms.string(""),
    #hltCaloJetsIDPassed = cms.string("hltCaloJetIDPassed"),
    hltCaloJetsL1Fast = cms.string("hltCaloJetL1FastJetCorrected"),
    hltCaloMETs = cms.string("hltMet"),
    hltCaloMETCleans = cms.string("hltMetClean"),
    hltCaloMETCleansUsingJetID = cms.string("hltMetCleanUsingJetID"),
    hltPFMETs = cms.string("hltPFMETProducer2"),
    #hltPFMETs = cms.string("hltPFMETProducer"),
    hltPFMETsNoMu = cms.string("hltPFMETNoMuProducer2"),
    #hltPFMETsNoMu = cms.string("hltPFMETnoMu"),
    hltPFMETCleansUsingJetID = cms.string("hltPFMETCleanUsingJetID"),
    hltTrackMETs = cms.string("hltTrackMETProducer2"),
    #hltTrackMETs = cms.string("hltTrackMetProducer"),
    hltCaloHTMHTs = cms.string("hltHtMht"),
    hltPFHTMHTs = cms.string(""),
    hltPFHTMHTsNoPU = cms.string("hltPFHTMHTNoPUProducer2"),
    #hltPFHTMHTsNoPU = cms.string("hltPFHTNoPU"),
    hltMuons = cms.string(""),
    hltPFCandidates = cms.string("hltParticleFlow"),
    hltPFJets = cms.string("hltAntiKT5PFJets"),
    hltPFJetsNoPU = cms.string("hltAntiKT5PFJetsNoPU"),
    hltPFJetsL1FastL2L3 = cms.string("hltAK5PFJetL1FastL2L3Corrected"),
    hltPFJetsL1FastL2L3NoPU = cms.string("hltAK5PFJetL1FastL2L3CorrectedNoPU"),
    hltPFRecTracks = cms.string("hltLightPFTracks"),
    hltVertices = cms.string("hltOnlinePrimaryVertices"),
    hltCSVBJetTags = cms.string("hltL3CombinedSecondaryVertexBJetTags"),
    hltRho_kt6CaloJets = cms.string("hltKT6CaloJets:rho"),
    hltRho_kt6PFJets = cms.string("hltKT6PFJets:rho"),
    # RECO
    recoTracks = cms.string(""),
    recoCaloJets = cms.string("ak5CaloJets"),
    recoCaloMETs = cms.string("met"),
    recoPFCandidates = cms.string("particleFlow"),
    recoPFJets = cms.string("ak5PFJets"),
    recoPFMETs = cms.string("pfMet"),
    recoPFMETT1s = cms.string("pfMetT1"),
    recoPFMETT0T1s = cms.string("pfMetT0pcT1"),
    recoPFMETMVAs = cms.string("pfMEtMVA"),
    recoPFMETNoPUs = cms.string("noPileUpPFMEt"),
    #recoVertices = cms.string(""),
    recoVertices = cms.string("offlinePrimaryVertices"),
    recoGoodVertices = cms.string("goodOfflinePrimaryVertices"),
    recoRho_kt6CaloJets = cms.string("kt6CaloJets:rho"),
    recoRho_kt6PFJets = cms.string("kt6PFJets:rho"),
    # PAT
    patElectrons = cms.string("selectedPatElectronsPFlow"),
    patJets = cms.string("selectedPatJetsPFlow"),
    patMETs = cms.string("patMETsPFlow"),
    patMuons = cms.string("selectedPatMuonsPFlow"),
    patPhotons = cms.string(""),
    patTaus = cms.string(""),
    # Trigger results
    triggerResults = cms.string("TriggerResults::HLT"),
    triggerEvent = cms.string("hltTriggerSummaryAOD"),
    #triggerPrescaleTable = cms.string("hltPrescaleRecorder"),
    # GEN
    genJets = cms.string("ak5GenJets"),
    genMETs = cms.string("genMetTrue"),
    genParticles = cms.string("genParticles"),
    genEventInfo = cms.string("generator"),
    simPileupInfo = cms.string("addPileupInfo"),
    # User
    hltPFPileUpFlags = cms.string("hltPFPileUpFlag"),
    patPFPileUpFlags = cms.string("pfPileUpFlag"),
    hltCaloJetIDs = cms.string("hltAK5JetID"),
    )

import csv
lumiA = []
lumiB = []
lumiC = []
with open('lumibyls_HLT_MET80_v5.csv') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if not row[7][0].isdigit():  continue
        runlumi = int(row[0].split(":")[0]) * 100000 + int(row[1].split(":")[0])
        effective = float(row[7])
        if effective > 1200:
            lumiC.append(runlumi)
        elif effective > 800:
            lumiB.append(runlumi)
        else:
            lumiA.append(runlumi)

process.lumicalc = cms.PSet(
    lumiA = cms.vuint64(lumiA),
    lumiB = cms.vuint64(lumiB),
    lumiC = cms.vuint64(lumiC),
    )
