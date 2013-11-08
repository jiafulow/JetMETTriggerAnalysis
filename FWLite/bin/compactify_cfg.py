import FWCore.ParameterSet.Config as cms
import re

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
    maxEvents   = cms.int32(-1),                            ## optional
    runMin      = cms.int32(-1),                            ## optional
    runMax      = cms.int32(-1),                            ## optional
    skipEvents  = cms.int32(0),                             ## optional
    reportEvery = cms.int32(1000),                          ## optional
    lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange(lumilist)),
)

if False:
    import os
    #dirname = "/pnfs/cms/WAX/resilient/jiafulow/METTriggers/skimHLTPFMET150_20130909/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball-PU25bx50/"
    #dirname = "/eos/uscms/store/user/jiafu/METTriggers/skimHLTPFMET150_eos_20130909/TT_CT10_TuneZ2star_8TeV-powheg-tauola-PU25bx50_2/"
    #dirname = "/pnfs/cms/WAX/resilient/jiafulow/METTriggers/skimHLTPFMET150_20130909/MET2012D/"
    dirname = "/pnfs/cms/WAX/resilient/jiafulow/METTriggers/skimHLTL1ETM40_20131007/MET-Run2012D-PromptReco-v1/"
    if not dirname.endswith("/"):  dirname += "/"
    basenamelist = os.listdir(dirname)
    if dirname.startswith("/pnfs/cms/"):
        dirname = "dcache:" + dirname
    if dirname.startswith("/eos/uscms"):  # some files have size = 0
        basenamelist = [basename for basename in basenamelist if os.path.getsize(dirname + basename) > 1000]
    process.input.fileNames = cms.vstring()
    process.input.fileNames.extend([dirname + basename for basename in basenamelist[:1]])
    #process.input.fileNames.extend([dirname + basename for basename in basenamelist])


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
        "HLT_PFMET150_v7",
        #"HLT_MET200_v12",
        #"HLT_MET200_HBHENoiseCleaned_v5",
    ),
    metfilters = cms.vstring(
        "p_HBHENoiseFilter",
        "p_CSCTightHaloFilter",
        "p_trackingFailureFilter",
        "p_EcalDeadCellTriggerPrimitiveFilter",
        "p_hcalLaserEventFilter",
        "p_ecalLaserCorrFilter",
        "p_eeBadScFilter",
        "p_trkPOGFilters",
    ),
    pfjetPtMin = cms.double(30),
    pfjetEtaMax = cms.double(5),
    pfjetEtaMaxCtr = cms.double(2.5),
    calojetPtMin = cms.double(30),
    calojetEtaMax = cms.double(5),
    calojetEtaMaxCtr = cms.double(2.6),
    calometcleanPtMin = cms.double(60),
    calometjetidPtMin = cms.double(60),
    isData = cms.bool(True),
    #verbose = cms.bool(True),
    verbose = cms.bool(False),
)

process.handler = cms.PSet(
    # L1
    l1METs = cms.string("hltL1extraParticles", "MET"),
    l1MHTs = cms.string("hltL1extraParticles", "MHT"),
    l1CentralJets = cms.string("hltL1extraParticles", "Central"),
    l1ForwardJets = cms.string("hltL1extraParticles", "Forward"),
    l1TauJets = cms.string("hltL1extraParticles", "Tau"),
    # HLT
    hltCaloJets = cms.string("hltAntiKT5CaloJets"),
    #hltCaloJets = cms.string("hltAntiKT5CaloJetsPF"),
    hltCaloJetIDPasseds = cms.string(""),
    #hltCaloJetIDPasseds = cms.string("hltCaloJetIDPassed"),
    hltCaloJetL1Fasts = cms.string("hltCaloJetL1FastJetCorrected"),
    hltCaloMETs = cms.string("hltMet"),
    hltCaloMETCleans = cms.string("hltMetClean"),
    hltCaloMETJetIDCleans = cms.string("hltMetJetIDClean"),
    hltPFMETs = cms.string("hltPFMETProducer"),
    hltPFMETNoMus = cms.string("hltPFMETnoMu"),
    hltTrackMETs = cms.string("hltTrackMetProducer"),
    hltMuons = cms.string(""),
    hltPFCandidates = cms.string("hltParticleFlow"),
    hltPFJets = cms.string("hltAntiKT5PFJets"),
    hltPFJetsNoPU = cms.string("hltAntiKT5PFJetsNoPU"),
    hltPFJetsL1FastL2L3 = cms.string("hltAK5PFJetL1FastL2L3Corrected"),
    hltPFJetsL1FastL2L3NoPU = cms.string("hltAK5PFJetL1FastL2L3CorrectedNoPU"),
    hltPFRecTracks = cms.string("hltLightPFTracks"),
    hltVertices = cms.string("hltOnlinePrimaryVertices"),
    hltRho_kt6CaloJets = cms.string("hltKT6CaloJets:rho"),
    hltRho_kt6PFJets = cms.string("hltKT6PFJets:rho"),
    # RECO
    recoTracks = cms.string(""),
    recoCaloJets = cms.string("ak5CaloJets"),
    recoCaloMETs = cms.string("met"),
    recoPFCandidates = cms.string("particleFlow"),
    recoPFJets = cms.string("ak5PFJets"),
    recoPFMETs = cms.string("pfMet"),
    recoVertices = cms.string(""),
    #recoVertices = cms.string("offlinePrimaryVertices"),
    recoGoodVertices = cms.string("goodOfflinePrimaryVertices"),
    recoRho_kt6CaloJets = cms.string("kt6CaloJets:rho"),
    recoRho_kt6PFJets = cms.string("kt6PFJets:rho"),
    # PAT
    patElectrons = cms.string("selectedPatElectronsPFlow"),
    patJets = cms.string("selectedPatJetsPFlow"),
    patMETs = cms.string("patPFMetPFlow"),
    patMETTypeIs = cms.string("patMETsPFlow"),
    patMuons = cms.string("selectedPatMuonsPFlow"),
    patPhotons = cms.string(""),
    patTaus = cms.string(""),
    # Trigger results
    triggerResults = cms.string("TriggerResults"),
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

