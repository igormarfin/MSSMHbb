#ifndef TRIPLEBTAGANALYSIS_SF_H
#define TRIPLEBTAGANALYSIS_SF_H

#include <string>
#include <vector>

#include "TFile.h"
#include "TH1.h"

#include "Analysis/Utilities/interface/BTag.h"
#include "Analysis/Utilities/interface/JetFlavor.h"
#include "Analysis/Utilities/interface/FlavorClass.h"
#include "Analysis/Utilities/interface/TrigHistArray.h"
#include "Analysis/Utilities/interface/TrigHistArray2D.h"


void prepareHistograms(const std::vector<std::string> &genericTriggerList,
		       const std::vector<std::string> &triggerFilterList,
		       const std::vector<JetFlavor::Code> &jetFlavorCodes,
		       const BTag::WPs &bTagWPs,
		       unsigned int nSelJet,
		       const std::vector<FlavorClass::Dijet> &flavorClassesDijet,
		       const std::vector<FlavorClass::Jet3> &flavorClasses3rdJet,
		       const std::vector<FlavorClass::Trijet> &flavorClassesTrijet,
		       int ncorr, int ncateg, int ntpat, bool _doMC);


/// output file
TFile* hout;

/// other variable declarations


/// histogram definitions
TrigHistArray* nJetA;
TrigHistArray* nJetPostselA;
TrigHistArray* matchPatternA;
TrigHistArray* dPhiJet1Jet2A;
// for monitoring the trigger weights
TrigHistArray* triggerWeightA;

std::vector< TrigHistArray* > mDijetBtA; // [nbtag]
std::vector< std::vector< TrigHistArray* > > ptJetBtA; // [nbtag][nSelJet]
std::vector< std::vector< TrigHistArray* > > mDijetFcBtA; // [nbtag][nfcTrip]
std::vector< std::vector< std::vector< TrigHistArray* > > > ptJetFcBtA; // [nbtag][nfcTrip][nSelJet]

std::vector< TrigHistArray* > mDijetBtwA; // [nbtag]
std::vector< std::vector< TrigHistArray* > > ptJetBtwA; // [nbtag][nSelJet]
std::vector< std::vector< TrigHistArray* > > mDijetFcBtwA; // [nbtag][nfcTrip]
std::vector< std::vector< std::vector< TrigHistArray* > > > ptJetFcBtwA; // [nbtag][nfcTrip][nSelJet]

std::vector< TrigHistArray* > dPhiJet1Jet3BtwA; // [nbtag]
std::vector< TrigHistArray* > dPhiJet2Jet3BtwA; // [nbtag]

std::vector< TrigHistArray* > svMassA; // [nSelJet]
std::vector< std::vector< TrigHistArray* > > svMassBtA; // [nbtag][nSelJet]

std::vector< TrigHistArray* > evBtagBtA; // [nbtag]
std::vector< TrigHistArray2D* > massEvBtagBtA; // [nbtag]

std::vector< TrigHistArray2D* > massEvBtagBtTWA; // [nbtag]

std::vector< TH1* > hfc; // [nbtag]
std::vector< TH1* > hfcm; // [nbtag]
std::vector< TH1* > hfcww; // [nbtag]
std::vector< TH1* > hfcmww; // [nbtag]

TH1* hfctrip;

std::vector< std::vector< std::vector< std::vector< std::vector< TrigHistArray* > > > > > massTemplateA; // [nfc][nbtag][ncateg][ncorr][ntpat]

std::vector< std::vector< std::vector< std::vector< TrigHistArray* > > > > bTagTemplateA; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > massBTagTemplateA; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > errorMassBTagTemplateA; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > massBTagTemplateUncldA; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > massBTagTemplateCldA; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > massBTagTemplateNonbbRA; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > massBTagTemplateCldRA; // [nfc][nbtag][ncateg][ntpat]

std::vector< std::vector< TrigHistArray* > > tpatA; // [nbtag][ncateg]
std::vector< TrigHistArray* > tpatAllA; // [nbtag]
std::vector< TrigHistArray* > atpattripall; // [nbtag]

std::vector< std::vector< TrigHistArray* > > mDibBtcutA; // [nbtag][ncateg]
std::vector< std::vector< TrigHistArray* > > mDibBtweightA; // [nbtag][ncateg]
std::vector< std::vector< std::vector< TrigHistArray* > > > ptDibBtcutA; // [nbtag][ncateg][nSelJet]
std::vector< std::vector< std::vector< TrigHistArray* > > > ptDibBtweightA; // [nbtag][ncateg][nSelJet]

std::vector< std::vector< std::vector< TH1* > > > mDibBtcutFcH; // [nfcDijet][nbtag][ncateg]
std::vector< std::vector< std::vector< TH1* > > > mDibBtweightFcH; // [nfcDijet][nbtag][ncateg]
std::vector< std::vector< std::vector< std::vector< TH1* > > > > hptdibbt; // [nfc][nAB][nbtag][ncateg]
  
std::vector< std::vector< std::vector< std::vector< TrigHistArray* > > > > massPred; // [nfc][nbtag][ncateg][ntpat]
std::vector< std::vector< std::vector< std::vector< TrigHistArray2D* > > > > massBTagPred; // [nfc][nbtag][ncateg][ntpat]

TH1* hptL3Objects_btagtrigger;
TH1* hptL3Objects_matchedtoL2_btagtrigger;
TH1* h_distance_L3Objects_vs_L2Objects_btagtrigger;


#endif
