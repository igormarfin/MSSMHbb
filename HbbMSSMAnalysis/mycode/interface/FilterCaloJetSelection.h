#ifndef FilterCaloJetSelection_H
#define FilterCaloJetSelection_H

// system include files
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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//
// class declaration
//


class FilterCaloJetSelection : public edm::EDFilter {
public:
  explicit FilterCaloJetSelection( const edm::ParameterSet & );
  ~FilterCaloJetSelection();
  
private:
  virtual bool filter ( edm::Event &, const edm::EventSetup & );
  
  bool applyfilter;
  bool debugOn;
  double ptminJet1;
  double ptminJet2;
  edm::InputTag jets_;

};

#endif


