import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:/afs/cern.ch/user/h/hongyi/public/powheg_nnlops.root'
    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/90000/FE30D7AE-081B-E911-B31F-0CC47AB35C3C.root')
)
process.demo = cms.EDAnalyzer('genAnalyzer',
genParticles = cms.untracked.InputTag('prunedGenParticles')
)


process.p = cms.Path(process.demo)
