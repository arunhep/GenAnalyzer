import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
    #'file:/afs/cern.ch/user/h/hongyi/public/powheg_nnlops.root'
    #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/WWTo2L2Nu_13TeV-powheg/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/90000/FE30D7AE-081B-E911-B31F-0CC47AB35C3C.root')
    #'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/WWJTo2L2Nu_NNLOPS_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/130000/02F65DD6-D6C8-144F-A9C9-ED87856487D3.root')
    'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/WWJTo2L2Nu_NNLOPS_TuneCUEP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/10000/0291D0AC-D7D6-E911-8745-246E96D74858.root')
)
process.demo = cms.EDAnalyzer('genAnalyzer',
genParticles = cms.untracked.InputTag('prunedGenParticles')
)


process.p = cms.Path(process.demo)
