import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
#        'file:/afs/cern.ch/work/a/arun/pythiaValidation/CMSSW_6_2_3/src/HIG-Fall13-00013_pythia6.root'
#        'file:/afs/cern.ch/work/n/nishu/public/HIG-Fall13-00013_pythia6.root'
	'file:/afs/cern.ch/work/n/nishu/public/HWW_MCNLO.root'
#         'file:/afs/cern.ch/work/n/nishu/public/HIG-Fall13-00013_Pythia8.root'
#         'file:/afs/cern.ch/work/a/arun/pythiaValidation/CMSSW_6_2_3/src/HIG-Fall13-00013_Pythia8.root'
    )
)

process.demo = cms.EDAnalyzer('genAnalyzer'
)

process.printTree = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(10),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("genParticles")
)

#process.p = cms.Path(process.demo + process.printTree)
process.p = cms.Path(process.demo)

