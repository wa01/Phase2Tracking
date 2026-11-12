# Imports
import FWCore.ParameterSet.Config as cms

# Create a new CMS process
from Configuration.Eras.Era_Phase2C22I13M9_cff import Phase2C22I13M9
process = cms.Process('tracktree',Phase2C22I13M9)

#process = cms.Process('tracktree')

# Import all the necessary files
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D121Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T35', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic_T33', '')

#process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
#from TrackingTools.TransientTrack.TransientTrackBuilder_cfi import TransientTrackBuilderESProducer
#process.TransientTrackBuilderESProducer = cms.ESProducer( "TransientTrackBuilderESProducer",
#  ComponentName = cms.string( "TransientTrackBuilder" ),
#  appendToDataLabel = cms.string( "" )
#)

# Number of events (-1 = all)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(10)
)
# Dataset: /RelValSingleMuPt10/CMSSW_15_1_0-150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_noPU-v1/GEN-SIM-RECO
# Input file
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring( \
'/store/relval/CMSSW_20_0_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_D121_RecycledGEN_noPU_16Aug26-v4/2590000/06b5209d-2ee5-4406-88a0-20610cdb7f07.root', \
'/store/relval/CMSSW_20_0_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_D121_RecycledGEN_noPU_16Aug26-v4/2590000/385d0125-1061-4587-b011-60b38fdd299b.root', \
'/store/relval/CMSSW_20_0_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_D121_RecycledGEN_noPU_16Aug26-v4/2590000/9fe4df4e-3095-4b53-b6de-7ca59cbdecf9.root', \
'/store/relval/CMSSW_20_0_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_D121_RecycledGEN_noPU_16Aug26-v4/2590000/aacffc32-4e49-4580-8d7e-74d4fea0c8b0.root', \
'/store/relval/CMSSW_20_0_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_D121_RecycledGEN_noPU_16Aug26-v4/2590000/ba0e16b2-36ca-43fe-b99b-ac92599ab49a.root'
#      'file:0166e3ce-e26d-483c-8226-cac472af699c.root'
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_noPU-v1/2590000/0166e3ce-e26d-483c-8226-cac472af699c.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/0562962b-66c8-488b-89d5-6df21d71f0b6.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/319ecc54-5d83-4e72-8d4a-3f73563e6d37.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/34e25999-54a2-4134-8703-241898703fed.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/37360339-d4dc-4270-85d1-8b89448a6504.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/49967669-571d-4a2b-aca2-6a148ebdb729.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/a88d679b-a909-47a4-8006-2455d1db2d5c.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/ae60ec5b-9009-4a17-8bcb-760c22430149.root', \
#      '/store/relval/CMSSW_15_1_0/RelValSingleMuPt10/GEN-SIM-RECO/150X_mcRun4_realistic_v1_STD_RegeneratedGS_Run4D110_n#oPU-v1/2590000/f931e711-1c11-47ce-9de3-49bd3966438d.root'
      )
)

# Output
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('file:track_tree_mu10.root')
)

#process.load('RecoLocalTracker.SiPhase2Clusterizer.phase2TrackerClusterizer_cfi')
#process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEESProducer_cfi')
##process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2StripCPEGeometricESProducer_cfi')
#process.load('RecoLocalTracker.Phase2TrackerRecHits.Phase2TrackerRecHits_cfi')
##process.siPhase2RecHits.Phase2StripCPE = cms.ESInputTag("phase2StripCPEESProducer", "Phase2StripCPE")
##process.siPhase2RecHits.Phase2StripCPE = cms.ESInputTag("phase2StripCPEGeometricESProducer", "Phase2StripCPEGeometric")


# Analyzer
process.analysis = cms.EDAnalyzer('TrackTree',
    rectracks = cms.InputTag("generalTracks"),
    simtracks = cms.InputTag("g4SimHits"),
    simvertices = cms.InputTag("g4SimHits"),
    simhitsbarrel = cms.InputTag("g4SimHits", "TrackerHitsPixelBarrelLowTof"),
    simhitsendcap = cms.InputTag("g4SimHits", "TrackerHitsPixelEndcapLowTof")
)

# Processes to run
#process.rechits_step = cms.Path(process.siPhase2Clusters + process.siPhase2RecHits)
#process.rechits_step = cms.Path(process.siPhase2RecHits)
process.analyze_step = cms.Path(process.analysis)

process.schedule = cms.Schedule(process.analyze_step)
