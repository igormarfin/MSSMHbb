#ifndef B_TAG_SF_H
#define B_TAG_SF_H

#include "BTag.h"
#include "JetFlavor.h"
#include "HbbSyst.h"


// Interface to b-tag efficiency data-MC scale factors as implemented
// in Utilities/interface/HbbSyst.h
class BTagSF {
public:
  // Return b-tag scale factor (Data/MC)
  double operator()(double jetPt, double jetEta, JetFlavor::Code code, BTag::WP wp) const;

  // Return uncertainty on b-tag scale factor
  double uncertaintyUp(double jetPt, double jetEta, JetFlavor::Code code, BTag::WP wp) const;
  double uncertaintyDown(double jetPt, double jetEta, JetFlavor::Code code, BTag::WP wp) const;


private:
  // Used internally to obtain scale factors and uncertainties
  mutable HbbSyst hbbSyst_;

  // To interface methods in HbbSyst
  int bTagToInt(BTag::WP wp) const;
};


// Return b-tag scale factor (Data/MC)
double BTagSF::operator()(double jetPt, double jetEta, JetFlavor::Code code, BTag::WP wp) const {
  double SF = 1.;

  if     ( code == JetFlavor::UDSG )
    SF = hbbSyst_.getSFlight(jetPt,jetEta,bTagToInt(wp));
  else if( code == JetFlavor::C || code == JetFlavor::CC || code == JetFlavor::B || code == JetFlavor::BB )
    SF = hbbSyst_.getSFbc(jetPt,jetEta,bTagToInt(wp));

  return SF;
}


// Return uncertainty on b-tag scale factor
double BTagSF::uncertaintyUp(double jetPt, double jetEta, JetFlavor::Code code, BTag::WP wp) const {
  double uncert = 0.;

  if     ( code == JetFlavor::UDSG )
    uncert = hbbSyst_.getSFlightUncertaintyUp(jetPt,jetEta,bTagToInt(wp));
  else if( code == JetFlavor::C || code == JetFlavor::CC )
    uncert = hbbSyst_.getSFcUncertaintyUp(jetPt,jetEta,bTagToInt(wp));
  else if( code == JetFlavor::B || code == JetFlavor::BB )
    uncert = hbbSyst_.getSFbUncertaintyUp(jetPt,jetEta,bTagToInt(wp));
  
  return uncert;
}

double BTagSF::uncertaintyDown(double jetPt, double jetEta, JetFlavor::Code code, BTag::WP wp) const {
  double uncert = 0.;

  if     ( code == JetFlavor::UDSG )
    uncert = hbbSyst_.getSFlightUncertaintyDown(jetPt,jetEta,bTagToInt(wp));
  else if( code == JetFlavor::C || code == JetFlavor::CC )
    uncert = hbbSyst_.getSFcUncertaintyDown(jetPt,jetEta,bTagToInt(wp));
  else if( code == JetFlavor::B || code == JetFlavor::BB )
    uncert = hbbSyst_.getSFbUncertaintyDown(jetPt,jetEta,bTagToInt(wp));
  
  return uncert;
}


// To interface methods in HbbSyst
int BTagSF::bTagToInt(BTag::WP wp) const {
  int tagger = -1;

  if     ( wp == BTag::TCHPT  ) tagger = 0;
  else if( wp == BTag::TCHP6  ) tagger = 1;
  else if( wp == BTag::CSVT   ) tagger = 2;
  else if( wp == BTag::SSVHPT ) tagger = 3;

  return tagger;
}
#endif
