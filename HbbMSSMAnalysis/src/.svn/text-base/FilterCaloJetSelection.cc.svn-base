#include <memory>
#include <vector>
#include <map>
#include <set>

// user include files

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "../interface/FilterCaloJetSelection.h"


using namespace edm;
using namespace std;

FilterCaloJetSelection::FilterCaloJetSelection(const edm::ParameterSet& iConfig)
{
  
  applyfilter = iConfig.getUntrackedParameter<bool>("applyfilter",true);
  debugOn     = iConfig.getUntrackedParameter<bool>("debugOn",false);
  ptminJet1 =  iConfig.getUntrackedParameter<double>("ptminJet1",46);
  ptminJet2 =  iConfig.getUntrackedParameter<double>("ptminJet2",38);
  jets_ = iConfig.getUntrackedParameter<edm::InputTag>("JetSource",edm::InputTag("ak5CaloJets"));
}

FilterCaloJetSelection::~FilterCaloJetSelection()
{
}

bool FilterCaloJetSelection::filter( edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool accepted = false;
  // get GeneralTracks collection

  edm::Handle<reco::CaloJetCollection> jets;
  iEvent.getByLabel(jets_,jets);    

//   std::cout << "Number of jets = " << jets->size() << std::endl;
  
  // check pt ordering of jets
  double lastpt = 99999999;
  const int maxJet = 2;
  double jetPtMin[maxJet];
  jetPtMin[0] = ptminJet1;
  jetPtMin[1] = ptminJet2;
  int nJet = 0;

  for (reco::CaloJetCollection::const_iterator jit=jets->begin(); jit != jets->end(); ++jit) {
//     std::cout << "Jet " << jit-jets->begin()
// 	      << " pt " << jit->pt()
// 	      << " eta " << jit->eta()
// 	      << std::endl;
    if (jit->pt() > lastpt) {
      std::cout << "Bad pt ordering, " << jit->pt() << " " << lastpt << std::endl;
    }
    lastpt = jit->pt();

    if (! (fabs(jit->eta())<2.6) ) continue;
    if (nJet<2) {
      if (jit->pt() > jetPtMin[nJet]) {
// 	std::cout << "Leading jet " << nJet << " found!" << std::endl;
	++nJet;
      }
    }
  }
  if (nJet >= 2) accepted= true;
    
  if (debugOn) {
    int ievt = iEvent.id().event();
    int irun = iEvent.id().run();
    int ils = iEvent.luminosityBlock();
    int bx = iEvent.bunchCrossing();
    
    std::cout << "FilterCaloJetSelection_debug: Run " << irun << " Event " << ievt << " Lumi Block " << ils << " Bunch Crossing " << bx << " NJetss " << jets->size() << " Accepted " << accepted << std::endl;
  }
 
  if (applyfilter)
    return accepted;
  else
    return true;

}

//define this as a plug-in
DEFINE_FWK_MODULE(FilterCaloJetSelection);
