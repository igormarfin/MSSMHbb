#ifndef TriggerRunSelection_h
#define TriggerRunSelection_h

#include <vector>
#include <string>
#include <algorithm>
#include <exception>

class TriggerRunSelection {
 public:
  std::string scenario;
  std::vector<unsigned> trgMask;
  std::vector<int> firstRun;
  std::vector<int> lastRun;
  std::vector<std::string>* gtl;
  std::vector<std::string>* filter;
  std::vector<std::string> usedtrigger;

  TriggerRunSelection(){}
  
  TriggerRunSelection(std::string theScenario,std::vector<std::string>* theGtl, const bool doMC = false) { //default is to run on data => doMC = false
    scenario = theScenario;
    gtl = theGtl;

   /*  if (theScenario == "LowMass2011") { */
/*       addRange(makePat("HLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D"),165970,166967); */
/*       addRange(makePat("HLT_CentralJet46_CentralJet38_DiBTagIP3D"),167039,168437); */
/*       addRange(makePat("HLT_CentralJet46_CentralJet38_DiBTagIP3D"),170826,173198); */
/*       addRange(makePat("HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D"),173236,175770); */
/*       addRange(makePat("HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D"),175832,178380); */
/*       addRange(makePat("HLT_CentralJet46_CentralJet38_DiBTagIP3D"),178420,180252); */
/*     } else if (theScenario == "MediumMass2011") { */
/*       addRange(makePat("HLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D"),165970,166967); */
/*       addRange(makePat("HLT_CentralJet46_CentralJet38_DiBTagIP3D"),167039,168437); */
/*       addRange(makePat("HLT_CentralJet60_CentralJet53_DiBTagIP3D"),170826,175770); */
/*       addRange(makePat("HLT_CentralJet60_CentralJet53_DiBTagIP3D"),175832,180252); */
/*     } else { */
/*       std::cout << "TriggerRunSelection Ctor: bad scenario " << theScenario << std::endl; */
/*     } */
    if (theScenario == "MediumMass2012") {
      addRangeTriggerName("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v",doMC ? 1 : 186785, doMC ? 9999999 : 209465);
  
    } else if (theScenario == "HighMass2012") {
      addRangeTriggerName("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV_v",doMC ? 1 : 186785, doMC ? 9999999 : 209465);
    
    } else if (theScenario == "VeryHighMass2012") {
      addRangeTriggerName("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose_v",doMC ? 1 : 186785, doMC ? 9999999 : 209465);
    
    } else {
      std::cout << "error: TriggerRunSelection Ctor: bad scenario " << theScenario << std::endl;
      throw std::exception();
    }

    //check trgMask and usedtrigger
    if(usedtrigger.size() != trgMask.size())
      {
        std::cout << "error: TriggerRunSelection::ERROR: usedtrigger.size() != trgMask.size() " << usedtrigger.size() << " " <<  trgMask.size() << " exit." << std::endl;
          throw std::exception();
      }

    print();
  }
 
  //wrapper for addRange which in parallel also stores the trigger name
  void addRangeTriggerName(std::string name, int theFirstRun, int theLastRun) {
    addRange(makePat(name),theFirstRun, theLastRun);
    usedtrigger.push_back(name);
  }
  //get the name of the trigger used in run theRun
  std::string getTriggerName(int theRun)
    {
      std::string returnname("");
      for (unsigned int iRange=0; iRange<trgMask.size(); ++iRange) {
        if(theRun >= firstRun[iRange] && theRun <= lastRun[iRange])
          {
            returnname = usedtrigger.at(iRange);
            break;
          }
      }
      if(returnname == "")
        {
          std::cout << "TriggerRunSelection::ERROR: no trigger name found in run " << theRun << "? exit." << std::endl;
          throw std::exception();
        }
      return returnname;
    }
  void addRange(unsigned int thePat, int theFirstRun, int theLastRun) {
    if (thePat != 0) {
      trgMask.push_back( thePat );
      firstRun.push_back( theFirstRun );
      lastRun.push_back( theLastRun );
    } else {
      std::cout << "TriggerRunSelection::addRange: empty trigger mask rejected" << std::endl;
    }
  }

  unsigned int makePat(std::string theTrig) {
    // create bit pattern mask for this trigger
    std::vector<std::string>::iterator tSlot = std::find(gtl->begin(), gtl->end(), theTrig);
    if (tSlot != gtl->end()) {
      unsigned int tNumber = tSlot - gtl->begin();
      std::cout << "TriggerRunSelection::makePat: trigger " << theTrig << " found at slot "
		<< tNumber << std::endl;
      return (1<<tNumber);
    } else {
      std::cout << "error TriggerRunSelection::makePat: trigger " << theTrig << " not found" << std::endl;
      throw std::exception();
      return 0;
    }
  }
  
  unsigned int mask(int theRun) {
    for (unsigned int iRange=0; iRange<trgMask.size(); ++iRange) {
      if ( (theRun >= firstRun[iRange]) && (theRun <= lastRun[iRange]) ) {
	return trgMask[iRange];
      }
    }
    return 0;
  }

  void print() {
    std::cout << "=========== TriggerRunSelection for scenario: " 
	      << scenario << " ============" << std::endl;
    std::cout << std::setw(5) << "iRg"  
	      << std::setw(10) << "firstRun"  
	      << std::setw(8) << "lastRun" 
	      << std::setw(8) << "mask" << std::endl;
    for (unsigned int iRange=0; iRange<trgMask.size(); ++iRange) {
      std::cout << std::setw(5) << iRange 
		<< std::setw(10) << firstRun[iRange] 
		<< std::setw(8) << lastRun[iRange] 
		<< std::setw(8) << trgMask[iRange] << std::endl;
    }
  }
};

    

#endif // #ifdef
