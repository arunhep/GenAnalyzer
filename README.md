# GenAnalyzer
Generator Level code used SMP HWW

### Commands to checkout and run the package
```
cmsrel CMSSW_10_6_4`
cd CMSSW_10_6_4/src
git clone git@github.com:arunhep/GenAnalyzer.git
git checkout WW_GenStudies
scramv1 b
cd GenAnalyzer/genAnalyzer/test
cmsRun ConfFile_cfg.py
```
```
source /cvmfs/cms.cern.ch/common/crab-setup.sh(csh) : depends on the shell being used
voms-proxy-init --voms cms
crab submit -c crab_cfg.py
```
### More about crab :
(https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial#CRAB_commands)
(https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile)

