# GenAnalyzer
Generator Level code used for MC Validation of HWW


cmsrel CMSSW\_10\_6\_4
cd CMSSW\_10\_6\_4/src
git clone git@github.com:arunhep/GenAnalyzer.git
git checkout WW\_GenStudies
scramv1 b
cd GenAnalyzer/genAnalyzer/test
cmsRun ConfFile\_cfg.py
source /cvmfs/cms.cern.ch/common/crab-setup.sh(csh) : depends on the shell being used
voms-proxy-init --voms cms
crab submit -c crab\_cfg.py

More about crab :
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial#CRAB_commands
https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile

