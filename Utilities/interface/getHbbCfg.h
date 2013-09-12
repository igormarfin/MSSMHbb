#ifndef getHbbCfg_H
#define getHbbCfg_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

int getHbbCfg(std::string& scenario,bool& useTemplateError,bool & useNP) {
  ifstream cfgFile;
  std::string theScenario;
  std::string theNPMode;
  std::string theBTMode;
  cfgFile.open("hbb.cfg");
  cfgFile >> theScenario >>  theNPMode >> theBTMode;
  cfgFile.close();

  if ( (theScenario == "LowMass2011") || (theScenario == "MediumMass2011") ||  (theScenario == "MediumMass2012") ) {
    scenario = theScenario;
  } else {
    std::cout << " Bad scenario " << std::endl;
    return -1;
  }

  if ( theNPMode == "NPon" ) {
    useNP = true;
  } else if (theNPMode == "NPoff" ) {
    useNP = false;
  } else {
    std::cout << " Bad NPmode " << std::endl;
    return -1;
  }

  if (theBTMode == "BTon") {
    useTemplateError = true;
  } else if (theBTMode == "BToff") {
    useTemplateError = false;
  } else {
    std::cout << " Bad BTmode " << std::endl;
    return -1;
  }
  std::cout << "getHbbCfg: scenario " << scenario 
	    << " useNP " << useNP
	    << " useTemplateError " << useTemplateError
	    << std::endl;
  return 0;
}



#endif
