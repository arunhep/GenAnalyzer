// -*- C++ -*-
//
// Package:    genAnalyzer
// Class:      genAnalyzer
// 
/**\class genAnalyzer genAnalyzer.cc GenAnalyzer/genAnalyzer/plugins/genAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Arun Kumar
//         Created:  Tue, 20 May 2014 12:10:21 GMT
// $Id$
//
//

// system include files
//#include <memory>
//#include <vector>
//#include <string>

#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TBranch.h"
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
//#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
//#include "DataFormats/TrackReco/interface/Track.h"
//#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "TH1.h"
#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
//#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//---- for LHE information
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEXMLStringProduct.h"
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "utils.h"

using namespace std;
using namespace reco;
using namespace edm;
//
// class declaration
//

class genAnalyzer : public edm::EDAnalyzer {
   public:
      explicit genAnalyzer(const edm::ParameterSet&);
      ~genAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  // ----------member data ---------------------------

  void Clear();
  void AddBranch(double* x, std::string name);
  void AddBranch(int* x, std::string name);
  void AddBranch(std::vector<double>*, std::string name);
  void AddBranch(std::vector<std::string>*, std::string name);
  void AddBranch(std::vector<int>*, std::string name);
  void SetBranches();
  bool isFromHiggs(const reco::GenParticle& p)const;

//  std::vector<double>& evtWeights = NULL;

  edm::InputTag genPartLabel_;

  TTree* tree_;
  TFile* file;

//This is a simple number
  int ngenMu_,ngenEle_,ngenTau_,ngenW_,ngenH_,ngenJet_,ngenNu_,ngenLept_;
  double startWeight_,numWeight_,denomWeight_,finalWeight_;

 std::vector<double> Event_;
 std::vector<double> Run_;
 std::vector<double> baseWeight_;
 std::vector<double> genNuPt_;
 std::vector<double> genNuEta_;
 std::vector<double> genNuPhi_;
 std::vector<double> genNuQ_;
 std::vector<int> genNustatus_;
 std::vector<int> genNuMother_;
 std::vector<int> genNuGrandMother_;
 std::vector<int> genNuPdgId_;

/*
  std::vector<double> genMuPt_;
  std::vector<double> genMuEta_;
  std::vector<double> genMuPhi_;
  std::vector<double> genMuQ_;
  std::vector<int> genMustatus_;
  std::vector<int> genMuMother_;
  std::vector<int> genMuGrandMother_;


  std::vector<double> genElePt_;
  std::vector<double> genEleEta_;
  std::vector<double> genElePhi_;
  std::vector<double> genEleQ_;
  std::vector<int> genElestatus_;
  std::vector<int> genEleMother_;
  std::vector<int> genEleGrandMother_;

  std::vector<double> genTauPt_;
  std::vector<double> genTauEta_;
  std::vector<double> genTauPhi_;
  std::vector<double> genTauQ_;
  std::vector<int> genTaustatus_;
  std::vector<int> genTauMother_;
  std::vector<int> genTauGrandMother_;
*/

  std::vector<double> genWPt_;
  std::vector<double> genWEta_;
  std::vector<double> genWPhi_;
  std::vector<double> genWMass_;
  std::vector<double> genWQ_;
  std::vector<int> genWstatus_;
  std::vector<int> genWMother_;

  std::vector<double> genHPt_;
  std::vector<double> genHEta_;
  std::vector<double> genHPhi_;
  std::vector<double> genHMass_;
  std::vector<double> genHQ_;
  std::vector<int> genHstatus_;
  std::vector<int> genHdaughter_;

  std::vector<double> genLeptPt_;
  std::vector<double> genLeptEta_;
  std::vector<double> genLeptPhi_;
  std::vector<double> genLeptM_;
  std::vector<int> genLeptId_;
  std::vector<int> genLeptStatus_;
  std::vector<double> genLeptMother_;
  std::vector<int> genLeptGrandMother_;
  std::vector<int> genLeptQ_;

  std::vector<double> genJetPt_;
  std::vector<double> genJetEta_;
  std::vector<double> genJetPhi_;
  std::vector<double> genJetMass_;
  std::vector<double> genCaloMET_;
  std::vector<double> genCaloMETPhi_;
  std::vector<double> genTrueMET_;
  std::vector<double> genTrueMETPhi_;
  std::vector<double> genNPMET_;
  std::vector<double> genNPMETPhi_;
};


//---------------------------------------------------------------
// Add Branches to the Tree
//---------------------------------------------------------------

bool genAnalyzer::isFromHiggs(const reco::GenParticle &p) const {
     for (unsigned int i = 0, n = p.numberOfMothers(); i < n; ++i) {
         const reco::GenParticleRef & mom = p.motherRef(i);
         if (abs(mom->pdgId()) == 25) return true;
         if ((11 <= abs(mom->pdgId()) && abs(mom->pdgId()) <= 24) && isFromHiggs(*mom)) return true;
     }
     return false;
}

void
genAnalyzer::AddBranch(double* x, std::string name){
  tree_->Branch(name.c_str(),x,(name+"/D").c_str());

}

void
genAnalyzer::AddBranch(int* x, std::string name){
    tree_->Branch(name.c_str(),x,(name+"/I").c_str());
    }

void
genAnalyzer::AddBranch(std::vector<double>* vec, std::string name){
    tree_->Branch(name.c_str(),vec);
   }

void
genAnalyzer::AddBranch(std::vector<int>* vec, std::string name){
  tree_->Branch(name.c_str(),vec);
  }

void
genAnalyzer::AddBranch(std::vector<std::string>* vec, std::string name){
   tree_->Branch(name.c_str(),vec);
  }

void  
genAnalyzer::SetBranches(){

  AddBranch(&Event_,"Event");
  AddBranch(&Run_,"Run");

  AddBranch(&genNuPdgId_,"genNuPdgId");
  AddBranch(&ngenNu_,"ngenNu");
  AddBranch(&genNuPt_, "genNuPt");
  AddBranch(&genNuEta_,"genNuEta");
  AddBranch(&genNuPhi_,"genNuPhi");
  AddBranch(&genNuQ_,"genNuQ");
  AddBranch(&genNustatus_,"genNustatus");
  AddBranch(&genNuMother_,"genNuMother");
  AddBranch(&genNuGrandMother_,"genNuGrandMother");

   AddBranch(&startWeight_,"startWeight");
   AddBranch(&numWeight_,"numWeight");
   AddBranch(&denomWeight_,"denomWeight");
   AddBranch(&finalWeight_,"finalWeight");
   AddBranch(&baseWeight_,"baseWeight");


/*
  AddBranch(&ngenMu_,"ngenMu");
  AddBranch(&genMuPt_, "genMuPt");
  AddBranch(&genMuEta_,"genMuEta");
  AddBranch(&genMuPhi_,"genMuPhi");
  AddBranch(&genMuQ_,"genMuQ");
  AddBranch(&genMustatus_,"genMustatus");
  AddBranch(&genMuMother_,"genMuMother");
  AddBranch(&genMuGrandMother_,"genMuGrandMother");

  AddBranch(&ngenEle_,"ngenEle");
  AddBranch(&genElePt_, "genElePt");
  AddBranch(&genEleEta_,"genEleEta");
  AddBranch(&genElePhi_,"genElePhi");
  AddBranch(&genEleQ_,"genEleQ");
  AddBranch(&genElestatus_,"genElestatus");
  AddBranch(&genEleMother_,"genEleMother");
  AddBranch(&genEleGrandMother_,"genEleGrandMother");

  AddBranch(&ngenTau_,"ngenTau");
  AddBranch(&genTauPt_, "genTauPt");
  AddBranch(&genTauEta_,"genTauEta");
  AddBranch(&genTauPhi_,"genTauPhi");
  AddBranch(&genTauQ_,"genTauQ");
  AddBranch(&genTaustatus_,"genTaustatus");
  AddBranch(&genTauMother_,"genTauMother");
  AddBranch(&genTauGrandMother_,"genTauGrandMother");
*/

  AddBranch(&ngenW_,"ngenW");
  AddBranch(&genWPt_, "genWPt");
  AddBranch(&genWEta_,"genWEta");
  AddBranch(&genWPhi_,"genWPhi");
  AddBranch(&genWMass_,"genWMass");
  AddBranch(&genWQ_,"genWQ");
  AddBranch(&genWstatus_,"genWstatus");
  AddBranch(&genWMother_,"genWMother");

  AddBranch(&ngenH_,"ngenH");
  AddBranch(&genHPt_, "genHPt");
  AddBranch(&genHEta_,"genHEta");
  AddBranch(&genHPhi_,"genHPhi");
  AddBranch(&genHMass_,"genHMass");
  AddBranch(&genHQ_,"genHQ");
  AddBranch(&genHstatus_,"genHstatus");
  AddBranch(&genHdaughter_,"genHdaughter");

  AddBranch(&ngenLept_, "ngenLept");
  AddBranch(&genLeptPt_, "genLeptPt");
  AddBranch(&genLeptEta_,"genLeptEta");
  AddBranch(&genLeptPhi_,"genLeptPhi");
  AddBranch(&genLeptM_,"genLeptM");
  AddBranch(&genLeptQ_,"genLeptQ");
  AddBranch(&genLeptStatus_,"genLeptStatus");
  AddBranch(&genLeptId_,"genLeptId");
  AddBranch(&genLeptMother_,"genLeptMother");
  AddBranch(&genLeptGrandMother_,"genLeptGrandMother");


  AddBranch(&genJetPt_, "genJetPt");
  AddBranch(&genJetEta_, "genJetEta");
  AddBranch(&genJetPhi_, "genJetPhi");
  AddBranch(&genJetMass_, "genJetMass");
  AddBranch(&ngenJet_, "ngenJet");

  AddBranch(&genCaloMET_, "genCaloMET");
  AddBranch(&genCaloMETPhi_, "genCaloMETPhi");
  AddBranch(&genTrueMET_, "genTrueMET");
  AddBranch(&genTrueMETPhi_, "genTrueMETPhi");
  AddBranch(&genNPMET_,"genNPMET");
  AddBranch(&genNPMETPhi_,"genNPMETPhi");


}

void  
genAnalyzer::Clear(){

//ngenMu_ = -999;
//ngenEle_ = -999;
//ngenTau_ = -999;
ngenW_ = -999;
ngenH_ = -999;
ngenJet_ = -999;
ngenNu_ = -999;
ngenLept_ = -999;
startWeight_ = 1.;
finalWeight_ = 1.;
numWeight_ = 1.;
denomWeight_ = 1.;

baseWeight_.clear();

 Event_.clear();
 Run_.clear();

  genLeptPt_.clear();
  genLeptEta_.clear();
  genLeptPhi_.clear();
  genLeptQ_.clear();
  genLeptStatus_.clear();
  genLeptMother_.clear();
  genLeptGrandMother_.clear();
  genLeptId_.clear();
  genLeptM_.clear();

  genNuPt_.clear();
  genNuEta_.clear();
  genNuPhi_.clear();
  genNuQ_.clear();
  genNustatus_.clear();
  genNuMother_.clear();
  genNuGrandMother_.clear();
  genNuPdgId_.clear();
/*
  genMuPt_.clear();
  genMuEta_.clear();
  genMuPhi_.clear();
  genMuQ_.clear();
  genMustatus_.clear();
  genMuMother_.clear();
  genMuGrandMother_.clear();

  genElePt_.clear();
  genEleEta_.clear();
  genElePhi_.clear();
  genEleQ_.clear();
  genElestatus_.clear();
  genEleMother_.clear();
  genEleGrandMother_.clear();

  genTauPt_.clear();
  genTauEta_.clear();
  genTauPhi_.clear();
  genTauQ_.clear();
  genTaustatus_.clear();
  genTauMother_.clear();
  genTauGrandMother_.clear();
*/
  genWPt_.clear();
  genWEta_.clear();
  genWPhi_.clear();
  genWMass_.clear();
  genWQ_.clear();
  genWstatus_.clear();
  genWMother_.clear();

  genHPt_.clear();
  genHEta_.clear();
  genHPhi_.clear();
  genHMass_.clear();
  genHQ_.clear();
  genHstatus_.clear();
  genHdaughter_.clear();

genJetPt_.clear();
genJetEta_.clear();
genJetPhi_.clear();
genJetMass_.clear();

genCaloMET_.clear();
genCaloMETPhi_.clear();
genTrueMET_.clear();
genTrueMETPhi_.clear();
genNPMET_.clear();
genNPMETPhi_.clear();

}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
genAnalyzer::genAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed

  file = new TFile("HWW_125GeV.root","recreate");
//  file = new TFile("WW.root","recreate");
  tree_ = new TTree("tree","tree");

// genPartLabel_ = iConfig.getParameter<edm::InputTag>("genParticles");

 SetBranches();


}


genAnalyzer::~genAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  delete tree_;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
genAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
Clear();
   using namespace edm;

     Event_.push_back(iEvent.id().event());
     Run_.push_back(iEvent.id().run());

   edm::Handle<reco::GenParticleCollection> genParticleHandle;
  // if(not iEvent.getByLabel(genPartLabel_, genParticleHandle))
   if(not iEvent.getByLabel("genParticles", genParticleHandle))
     {
       std::cout<<
	 "GenAnalyzer: Generator Level Information not found\n"
		<<std::endl;
     }

//   std::cout<<"Running\n"<<std::endl;

//int nMu = 0;
//int nEle = 0;
//int nTau = 0;
int nZ = 0;
int nH = 0;
int nNu = 0;
int nLept = 0;
int njets = 0;

   const reco::GenParticleCollection* genColl= &(*genParticleHandle);
//   std::vector<reco::GenParticle> sorted(*genColl);
//   reco::GenParticleCollection genColl(*(genParticleHandle.product()));  //Obtain the electron collection.
//   for(reco::GenParticleCollection::const_iterator geni = genColl.begin(); geni != genColl.end(); ++geni){

     std::vector<const reco::GenParticle *> sortedPtrs;
     sortedPtrs.reserve(genColl->size());
     for (const reco::GenParticle &g : *genColl) { sortedPtrs.push_back(&g); }
     std::sort(sortedPtrs.begin(), sortedPtrs.end(),PtGreater());

for (auto const & genPtr : sortedPtrs) {
   auto const & geni = *genPtr;


//if(isFromHiggs(geni) == 1) {

//if(abs(geni.pdgId()) ==11 && geni.status() == 1){
//if(((abs(geni.pdgId())==11) || (abs(geni.pdgId())==13) || (abs(geni.pdgId())==15)) && ((geni.mother()->pdgId() == 24) || (geni.mother()->pdgId() == -24)) && (geni.mother()->mother()->pdgId() == 25)){
//if((abs(geni.pdgId())==11) || (abs(geni.pdgId())==13) || (abs(geni.pdgId())==15)){
if((abs(geni.pdgId())==11 || abs(geni.pdgId())==13) && (abs(geni.mother()->pdgId()) == 24 || (abs(geni.mother()->pdgId()) == 15 && abs(geni.mother()->mother()->pdgId()) == 24))){
//if((abs(geni.pdgId())==11 || abs(geni.pdgId())==13) && (abs(geni.mother()->pdgId()) == 24 || abs(geni.mother()->pdgId()) == 15)){
//genLeptPt_.push_back(geni.pt());
genLeptPt_.push_back(geni.pt());
genLeptEta_.push_back(geni.eta());
genLeptPhi_.push_back(geni.phi());
genLeptM_.push_back(geni.mass());
genLeptId_.push_back(geni.pdgId());
genLeptStatus_.push_back(geni.status());
genLeptMother_.push_back(geni.mother()->pdgId());
genLeptGrandMother_.push_back(geni.mother()->mother()->pdgId());
//cout << "geni.pdgId() = " << geni.pdgId() << " Mother = " <<  geni.mother()->pdgId() << " GrandMother = " << geni.mother()->mother()->pdgId() << endl;
//cout << "Lept Status = " << geni.status() << endl;
nLept++;
}
ngenLept_ = nLept;

//if(((abs(geni.pdgId())==12) || (abs(geni.pdgId())==14) || (abs(geni.pdgId())==16)) && ((geni.mother()->pdgId() == 24) || (geni.mother()->pdgId() == -24)) && (geni.mother()->mother()->pdgId() == 25)){
if((abs(geni.pdgId())==12 || abs(geni.pdgId())==14 || abs(geni.pdgId())==16) && (geni.status() == 1)){

       genNuPt_.push_back(geni.pt());
       genNuEta_.push_back(geni.eta());
       genNuPhi_.push_back(geni.phi());
       genNuQ_.push_back(geni.charge());
       genNustatus_.push_back(geni.status());
       genNuMother_.push_back(geni.mother()->pdgId());
       genNuGrandMother_.push_back(geni.mother()->mother()->pdgId());
       genNuPdgId_.push_back(geni.pdgId());
nNu++;
}
ngenNu_ = nNu;

//if((abs(geni.pdgId())==16)) {
//cout << "geni.status = " << geni.status() << endl;
//}

/*
if((abs(geni.pdgId())==13) && ((geni.mother()->pdgId() == 24) || (geni.mother()->pdgId() == -24)) && (geni.mother()->mother()->pdgId() == 25)){
       genMuPt_.push_back(geni.pt());
       genMuEta_.push_back(geni.eta());
       genMuPhi_.push_back(geni.phi());
       genMuQ_.push_back(geni.charge());
       genMustatus_.push_back(geni.status());
       genMuMother_.push_back(geni.mother()->pdgId());
       genMuGrandMother_.push_back(geni.mother()->mother()->pdgId());

 nMu++;
}
ngenMu_ = nMu;

//if(nMu == 2){
//cout << "nMu = " << nMu << endl;
//}

if((abs(geni.pdgId())==11) && ((geni.mother()->pdgId() == 24) || (geni.mother()->pdgId() == -24)) && (geni.mother()->mother()->pdgId() == 25)){
//if((abs(geni.pdgId())==11)){
       genElePt_.push_back(geni.pt());
       genEleEta_.push_back(geni.eta());
       genElePhi_.push_back(geni.phi());
       genEleQ_.push_back(geni.charge());
       genElestatus_.push_back(geni.status());
       genEleMother_.push_back(geni.mother()->pdgId());
       genEleGrandMother_.push_back(geni.mother()->mother()->pdgId());

 nEle++;
}
ngenEle_ = nEle;

if((abs(geni.pdgId())==15) && ((geni.mother()->pdgId() == 24) || (geni.mother()->pdgId() == -24)) && (geni.mother()->mother()->pdgId() == 25)){
//if((abs(geni.pdgId())==15)){
       genTauPt_.push_back(geni.pt());
       genTauEta_.push_back(geni.eta());
       genTauPhi_.push_back(geni.phi());
       genTauQ_.push_back(geni.charge());
       genTaustatus_.push_back(geni.status());
       genTauMother_.push_back(geni.mother()->pdgId());
       genTauGrandMother_.push_back(geni.mother()->mother()->pdgId());

 nTau++;
}
ngenTau_ = nTau;
*/

//if(abs(geni.pdgId())==24 && geni.mother()->pdgId() != 25) {
//cout << "geni.mother()->pdgId() = " << geni.mother()->pdgId() << endl;
//}

if((abs(geni.pdgId())==24) && (geni.mother()->pdgId() == 25)){
//if((abs(geni.pdgId())==23)){
       genWPt_.push_back(geni.pt());
       genWEta_.push_back(geni.eta());
       genWPhi_.push_back(geni.phi());
       genWMass_.push_back(geni.mass());
       genWQ_.push_back(geni.charge());
       genWstatus_.push_back(geni.status());
       genWMother_.push_back(geni.mother()->pdgId());

 nZ++;
//cout << geni.pt() << endl;
}
ngenW_ = nZ;
//cout << "nH (before) = " << nH << endl;
if(geni.pdgId()==25 && geni.numberOfDaughters() == 2){
//if(geni.pdgId()==25){
//cout << "geni.pdgId() = " << geni.pdgId() << endl;
       genHPt_.push_back(geni.pt());
       genHEta_.push_back(geni.eta());
       genHPhi_.push_back(geni.phi());
       genHMass_.push_back(geni.mass());
       genHQ_.push_back(geni.charge());
       genHstatus_.push_back(geni.status());
//if(geni.status() != 62) {
//cout << geni.status() << endl;
//}
     int n = geni.numberOfDaughters();
     for(int j = 0; j < n; ++ j) {
       const Candidate *d = geni.daughter(j);
       genHdaughter_.push_back(d->pdgId());
//cout << "daughter id = " << d->pdgId() << endl;

     }

 nH++;
//cout << nH << endl;
}
//cout << "nH (After) = " << nH << endl;
ngenH_ = nH;
//cout << "ngenH_ = " << ngenH_ << endl;

}


//cout << "nLept = " << nLept << endl;


//cout << "**************************************************" << endl;

///////////////////////////////////////// looking for genjets //////////////////////////////////////



   edm::Handle<reco::GenJetCollection> genJetsHandle;
   if( not iEvent.getByLabel("ak4GenJets",genJetsHandle)){
     edm::LogInfo("GenAnalyzer") << "genJets not found, "
       "skipping event";
     return;
   }


   const reco::GenJetCollection* genJetColl= &(*genJetsHandle);
//   std::vector<reco::GenParticle> sorted(*genColl);
//   reco::GenParticleCollection genColl(*(genParticleHandle.product()));  //Obtain the electron collection.
//   for(reco::GenParticleCollection::const_iterator geni = genColl.begin(); geni != genColl.end(); ++geni){

     std::vector<const reco::GenJet *> sortedjets;
     sortedjets.reserve(genJetColl->size());
     for (const reco::GenJet &g : *genJetColl) { sortedjets.push_back(&g); }
     std::sort(sortedjets.begin(),sortedjets.end(),PtGreater());

//ngenJet_ = sortedjets.size();
//cout << "Number of jets Initial = " << sortedjets.size() << endl;
for (auto const & genPtrjets : sortedjets) {
   auto const & gjets = *genPtrjets;
//cout << "Jets PT = " << gjets.pt() << "   Jets Eta = " << gjets.eta() << endl;
if(gjets.pt() > 30. && fabs(gjets.eta()) < 5.0) {
     genJetPt_.push_back(gjets.pt());
     genJetEta_.push_back(gjets.eta());
     genJetPhi_.push_back(gjets.phi());
     genJetMass_.push_back(gjets.mass());
njets++;
}
//cout << gjets.mass() << endl;
}
ngenJet_ = njets;
//cout << "number of jets After = " << njets << endl;
//cout <<"*****************************"<< endl;

//	edm::Handle<std::vector<reco::PFMET>> genCaloMETHandle;
        edm::Handle<reco::GenMETCollection> genCaloMETHandle;	
	if( not iEvent.getByLabel("genMetCalo", genCaloMETHandle)){
	edm::LogInfo("GenAnalyzer") << "genCaloMET not found, "
					"skipping event";
return;
}

const reco::GenMETCollection* genMETColl = &(*genCaloMETHandle);
reco::GenMETCollection::const_iterator i;
for(i=genMETColl->begin(); i!=genMETColl->end(); i++){

genCaloMET_.push_back(i->pt());
genCaloMETPhi_.push_back(i->phi());

}

        edm::Handle<reco::GenMETCollection> genTrueMETHandle;
        if( not iEvent.getByLabel("genMetTrue", genTrueMETHandle)){
        edm::LogInfo("GenAnalyzer") << "genTrueMET not found, "
                                        "skipping event";
return;
}

const reco::GenMETCollection* genTrueMETColl = &(*genTrueMETHandle);
reco::GenMETCollection::const_iterator j;
for(j=genTrueMETColl->begin(); j!=genTrueMETColl->end(); j++){

genTrueMET_.push_back(j->pt());
genTrueMETPhi_.push_back(j->phi());

}

        edm::Handle<reco::GenMETCollection> genNPMETHandle;
        if( not iEvent.getByLabel("genMetCaloAndNonPrompt", genNPMETHandle)){
        edm::LogInfo("GenAnalyzer") << "genNPMET not found, "
                                        "skipping event";
return;
}

const reco::GenMETCollection* genNPMETColl = &(*genNPMETHandle);
reco::GenMETCollection::const_iterator k;
for(k=genNPMETColl->begin(); k!=genNPMETColl->end(); k++){

genNPMET_.push_back(k->pt());
genNPMETPhi_.push_back(k->phi());

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++ FOR MCNLO WEIGHTS +++++++++++++++++++++++++++++++++++++++++++


edm::Handle<GenEventInfoProduct> genEvtInfo;
if( not iEvent.getByLabel("generator",genEvtInfo)){
        edm::LogInfo("GenAnalyzer") << "GenEventInfo not found, "
                                        "skipping event";
return;
}
//double qScale = genEvtInfo->qScale();  // in case of Pythia6, this will be pypars/pari(23)
//const std::vector<double>& binningValues = genEvtInfo->binningValues(); 

std::vector<double> evtWeights = genEvtInfo->weights();
double theWeight = genEvtInfo->weight();

startWeight_ = theWeight;

//cout << "Base Weight =  " << theWeight << endl;

 for (unsigned int iWeight = 0; iWeight < evtWeights.size(); iWeight++) {
//  std::cout << " evtWeights[" << iWeight << "] = " << evtWeights.at(iWeight) << std::endl;
  baseWeight_.push_back(evtWeights.at(iWeight));
 }


edm::Handle<LHEEventProduct> EvtHandle ;
if (not iEvent.getByLabel("externalLHEProducer",EvtHandle)){
edm::LogInfo("GenAnalyzer") << "GenEventInfo not found, "
                                        "skipping event";
return;
}

 unsigned int num_whichWeight = EvtHandle->weights().size();
//cout << "number of weights = " << num_whichWeight << endl;
 for (unsigned int iWeight = 0; iWeight < num_whichWeight; iWeight++) {
//  _weightsLHE.push_back( productLHEHandle->weights()[iWeight].wgt/productLHEHandle->originalXWGTUP() ); 
//  std::cout << " weightLHE[" << iWeight << "] = " << (EvtHandle->weights()[iWeight].wgt/EvtHandle->originalXWGTUP()) << std::endl;
//  cout << "productLHEHandle->originalXWGTUP() = " << EvtHandle->originalXWGTUP() << endl;
 }

 std::vector<std::string> comments_LHE;
 std::vector<std::string>::const_iterator it_end = (*(EvtHandle.product())).comments_end();
  std::vector<std::string>::const_iterator it = (*(EvtHandle.product())).comments_begin();
  for(; it != it_end; it++) {
   comments_LHE.push_back (*it);
  }

 for (unsigned int iComm = 0; iComm<comments_LHE.size(); iComm++) {
//cout << "comments_LHE = " << comments_LHE.at(iComm) << endl;
}

//edm::Handle<LHERunInfoProduct> run; 
//typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;
 
//iRun.getByLabel("externalLHEProducer",run);
//LHERunInfoProduct myLHERunInfoProduct = *(run.product());
 
//for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
//  std::cout << iter->tag() << std::endl;
//  std::vector<std::string> lines = iter->lines();
//  for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
//   std::cout << lines.at(iLine);
//  }
//}


int whichWeight = 1;

//cout << "theWeight before = " << theWeight << endl;

//cout << "Numerator =  " << EvtHandle->weights()[1].wgt << endl;
//cout << "Denominator = " << EvtHandle->originalXWGTUP() << endl;

numWeight_ = EvtHandle->weights()[whichWeight].wgt;
denomWeight_ = EvtHandle->originalXWGTUP();

theWeight *= EvtHandle->weights()[whichWeight].wgt/EvtHandle->originalXWGTUP(); 

finalWeight_ = theWeight;

//cout << " theWeight after = " << theWeight << endl;

//int nAllPartons = genEvtInfo->nMEPartons();
//int nPartonsEnteringMatching = genEvtInfo->nMEPartonsFiltered();

tree_->Fill();
//cout << "*****************************************************" <<endl;

}


// ------------ method called once each job just before starting event loop  ------------
void 
genAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
genAnalyzer::endJob() 
{
  file->cd();
  file->Write();
}

// ------------ method called when starting to processes a run  ------------
/*
void 
genAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
genAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
genAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
genAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
genAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(genAnalyzer);
