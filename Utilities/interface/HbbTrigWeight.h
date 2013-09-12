#ifndef HbbTrigWeight_h
#define HbbTrigWeight_h

#include "Analysis/Utilities/interface/HbbTrigger.h"


class HbbTrigWeight {
  // helper class to interface parameterized trigger efficiencies from HbbTrigger.h
  // to TripleBtagAnalysis macro
 public:
  HbbTrigger*** triggers;        // triggers[nGtl][nPeriod]
  std::vector<std::string>* gtl;  // genericTriggerList
  unsigned int filterPattern;
  unsigned int nPeriod;

  HbbTrigWeight() {
    std::cout << "HbbTrigWeight: Standard constructor called, this should not happen" << std::endl;
  }

  HbbTrigWeight(std::vector<std::string>* theGtl,std::vector<std::string>* theFilter) {
    gtl = theGtl;
    nPeriod = 2;
    // first set the filter pattern
    filterPattern = 0;
    if (theFilter != NULL) {
      for (std::vector<std::string>::iterator fit=theFilter->begin(); fit != theFilter->end(); ++fit) {
	std::vector<std::string>::iterator tSlotBtag = std::find(theGtl->begin(), theGtl->end(), *fit);
	if (tSlotBtag != theGtl->end()) {
	  filterPattern = filterPattern | (1<<(tSlotBtag - theGtl->begin()));
	} else {
	  std::cout << "HbbTrigWeight: Filter trigger " << *fit << " not found in any slot" << std::endl;
	}
      }
    } else {
      filterPattern = ~0;  // set all filter bits to one
    }
    
    // create array that will hold pointer to one HbbTrigger object for each generic trigger
    triggers = new HbbTrigger**[gtl->size()];
    for (unsigned int iGtl=0; iGtl<gtl->size(); ++iGtl) {
      triggers[iGtl] = new HbbTrigger*[nPeriod];
    }
    for (unsigned int ib=0; ib<gtl->size(); ++ib) {
      for (unsigned int iPeriod=0; iPeriod<nPeriod; ++iPeriod) {
	triggers[ib][iPeriod] = NULL;
	int robervalTNumber = -1;
	if (filterPattern & (1<<ib)) {
	  if ( ((*gtl)[ib] == "HLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D")
	       || ((*gtl)[ib] == "HLT_CentralJet46_CentralJet38_DiBTagIP3D") ) {
	    robervalTNumber = 0;
	  } else if ((*gtl)[ib] == "HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D") {
	    robervalTNumber = 1;
	  } else if ((*gtl)[ib] == "HLT_CentralJet60_CentralJet53_DiBTagIP3D") {
	    robervalTNumber = 2;
	  }
	  std::cout << "HbbTrigWeight: assign robervalTNumber= " << robervalTNumber
		    << " for trigger " << (*gtl)[ib] << std::endl;
	}
	if ( robervalTNumber >= 0 ) {
	  triggers[ib][iPeriod] = new HbbTrigger( robervalTNumber );
	  if (  robervalTNumber == 1 ) triggers[ib][iPeriod]->setPeriod( iPeriod );
	}
      }
    }
  }

  void getTrigWeight(bTagEff* theBTagReleffOnline,const char* theBtag,int* theFc,
		     float* thePt,float* theEta,float* weightArray,const bool L1L2Only=false) {
    // first compute the online btag relative efficiency
    float obeff[3];
    
    for (int iJ=0; iJ<3; ++iJ) {
      // correct for the fact that MinTag4 has only 94% of efficiency wrt. MinTag3
      obeff[iJ] = theBTagReleffOnline->eff(theFc[iJ],theBtag,thePt[iJ],theEta[iJ]) / 0.94;
    }
    float obtRel = obeff[0] * obeff[1] + obeff[0] * obeff[2] + obeff[1] * obeff[2]
      - 2 * obeff[0] * obeff[1] * obeff[2];
    if ( (obtRel < 0) || (obtRel > 1.) ) {
      std::cout << "HbbTrigWeight::GetTrigWeight: bad obtRel = " << obtRel << std::endl;
      obtRel = 0;
    }
    
    if (L1L2Only) obtRel = 1;
    
    for (unsigned int ib=0; ib<gtl->size(); ++ib) {
      weightArray[ib] = 0;
      if ( (triggers[ib][0] != NULL) && (triggers[ib][1] != NULL) ) {
	// mix according to period fractions
	float frac0 = 246.039 / (246.039 + 731.275);
	float frac1 = 1. - frac0;
	weightArray[ib] = obtRel * ( frac0 * triggers[ib][0]->getEfficiency(thePt[0],thePt[1],thePt[2])
				     + frac1 * triggers[ib][1]->getEfficiency(thePt[0],thePt[1],thePt[2]) );
      }
/*       if (ib == 0) { */
/* 	std::cout << "flav= " << theFc[0] << " " << theFc[1] << " " << theFc[2] */
/* 		  << " pt= " <<  thePt[0] << " " << thePt[1] */
/* 		  << " L1L2= " << triggers[ib]->getEfficiency(thePt[0],thePt[1],thePt[2]) */
/* 		  << " obtRel= " << obtRel << std::endl; */
/*       } */
    }
  }
};
#endif
