from CRABClient.UserUtilities import config
config = config()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'WW_NNLOPS'
config.General.transferLogs = False
#config.General.workarea = 'ww_analysis'
config.General.workArea = 'crab_WW_NNLOPS'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = './ConfFile_cfg.py'
#config.JobType.psetName    = '../gendumper_cfg.py'
config.JobType.outputFiles = ['WW_pair_production.root']
#config.JobType.pyCfgParams = ['outputFile=GEN_MC_ggHww.root']
config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDBS = 'global'
config.Data.inputDataset = '/WWJTo2L2Nu_NNLOPS_TuneCUEP8M1_13TeV-powheg-pythia8/RunIISummer16MiniAODv3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1/MINIAODSIM'
config.Data.splitting    = 'FileBased'  #'LumiBased'
config.Data.unitsPerJob  = 10  # Since files based, 10 files per job
config.Data.publication = False
#config.Data.inputDBS     = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader/'
#config.Data.outLFN       = '/store/user/hongyi/ww_data'_
config.Data.outLFNDirBase       = '/store/group/phys_higgs/cmshww/arun/WW/MCStudies/WW_NNLOPS'

config.section_('Site')
config.Site.storageSite = 'T2_CH_CERN'
