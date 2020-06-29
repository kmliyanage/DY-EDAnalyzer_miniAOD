import FWCore.ParameterSet.Config as cms

process = cms.Process('DemoAnalyzer')
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/b/bkailasa/public/DarkSusy-RunIIAutumn18DR_step3-MiniAOD.root'))
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histos1.root')
                                   )

#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#rom Configuration.AlCa.GlobalTag import GlobalTag
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
