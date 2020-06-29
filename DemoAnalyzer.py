import FWCore.ParameterSet.Config as cms

process = cms.Process("DemoAnalyzer")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.INFO = cms.untracked.PSet(
                                                     limit = cms.untracked.int32(0)
                                                     )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False),
                                     SkipEvent = cms.untracked.vstring('ProductNotFound'))

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring())

process.source.fileNames =[#'file:./pat.root',
'/store/mc/RunIIFall17MiniAODv2/ZToMuMu_NNPDF31_13TeV-powheg_M_50_120/MINIAODSIM/MUOTrackFix_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/90000/50A3D4D2-C6FE-E811-AE32-24BE05C6D731.root'
                           ]

process.demo = cms.EDAnalyzer('DemoAnalyzer')


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histos1.root')
                                   )

process.p = cms.Path(process.demo)

