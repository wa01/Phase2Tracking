# Imports
import FWCore.ParameterSet.Config as cms

# Create a new CMS process
process = cms.Process('hittree')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')


# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_15_0_0/RelValSingleMuPt10/GEN-SIM-RECO/141X_mcRun4_realistic_v3_STD_RegeneratedGS_Run4D110_noPU-v2/2580000/0b0d313e-56e1-4e64-aa58-8a2e61767bf5.root',
        #'root://eoscms.cern.ch//eos/cms/store/relval/CMSSW_15_0_0/RelValSingleMuPt10/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v2/2580000/012df10b-6e3c-4497-94be-c382898f2416.root',
        '/store/relval/CMSSW_15_0_0_pre3/RelValTTbar_14TeV/GEN-SIM-RECO/PU_141X_mcRun4_realistic_v3_STD_Run4D110_PU-v2/2580000/008977ec-54ed-42c7-a69f-7df26681d0c6.root',
    )
)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('rechits_tree.root')
)

process.load('RecoLocalTracker.SiPhase2Clusterizer.phase2TrackerClusterizer_cfi')
process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEESProducer_cfi')
#process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEGeometricESProducer_cfi')
process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2TrackerRecHits_cfi')
#process.siPhase2RecHits.Phase2StripCPE = cms.ESInputTag("phase2StripCPEESProducer", "Phase2StripCPE")
#process.siPhase2RecHits.Phase2StripCPE = cms.ESInputTag("phase2StripCPEGeometricESProducer", "Phase2StripCPEGeometric")


# Analyzer
process.analysis = cms.EDAnalyzer('RecHitTree',
    rechits = cms.InputTag("siPhase2RecHits"),
    clusters = cms.InputTag("siPhase2Clusters"),
    links = cms.InputTag("simSiPixelDigis", "Tracker"),
    simhitsbarrel = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
    simhitsendcap = cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof"),
    simtracks = cms.InputTag("g4SimHits"),
    genparticles = cms.InputTag("genParticles"),
    #ECasRings = cms.bool(True),
    SimTrackMinPt = cms.double(0.5),
    #MakeEtaPlots = cms.bool(False),
    #MinEta = cms.double(0.),
    #MaxEta = cms.double(10.)
)

# Processes to run
#process.rechits_step = cms.Path(process.siPhase2Clusters + process.siPhase2RecHits)
process.rechits_step = cms.Path(process.siPhase2RecHits)
process.analyze_step = cms.Path(process.analysis)

process.schedule = cms.Schedule(process.rechits_step, process.analyze_step)
