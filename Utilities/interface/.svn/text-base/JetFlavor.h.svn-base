#ifndef JET_FLAVOR_H
#define JET_FLAVOR_H

#include <vector>
#include <string>

#include "Analysis/Utilities/interface/HbbNtuple.h"


// Encapsulates the jet flavor code definition
// and returns the flavor code of a given jet
class JetFlavor {
public:

  // Flavor codes defined in the analysis
  // Assign fixed values for backward compatibility
  enum Code { UNDEFINED = -1,
	      UDSG      =  0,
	      C         =  1,
	      B         =  2,
	      CC        =  3,
	      BB        =  4  };


  // For convenience
  typedef std::vector<Code> Codes;
  typedef std::vector<Code>::const_iterator CodesIt;


  // Flavor codes used in the analysis
  static std::vector<JetFlavor::Code> getList();

  // std::string representation of JetFlavor::Code
  static std::string toString(JetFlavor::Code flav);

  // Return flavor code of jet
  static JetFlavor::Code code(const int iJet);

  // Return dijet flavor code. The dijet flavor code maps
  // the flavor code of both jets to a unique number
  //
  // TODO if still needed, this should also become a readable enum
  static int diJetCode(const int iJet0, const int iJet1);


  // for backward compatibility; if you have to use
  // any of the following methods your code is outdated
  static JetFlavor::Code fromInt(int flav);
  static JetFlavor::Code fromString(const std::string &flav);



private:

  // Returns true for udscb quarks and gluons,
  // false otherwise
  static bool validPdgID(int pdgID);
};


// Set of flavor codes used in the analysis
std::vector<JetFlavor::Code> JetFlavor::getList() {
  std::vector<JetFlavor::Code> list;
  list.push_back(UDSG);
  list.push_back(C);
  list.push_back(B);
  list.push_back(CC);
  list.push_back(BB);

  return list;
}


// std::string representation of JetFlavor::Code
std::string JetFlavor::toString(JetFlavor::Code flav) {
  std::string str = "undefined";

  if     ( flav == UDSG ) str = "udsg";
  else if( flav == C    ) str = "c";
  else if( flav == B    ) str = "b";
  else if( flav == CC   ) str = "cc";
  else if( flav == BB   ) str = "bb";

  return str;
}


// Return flavor code of jet
JetFlavor::Code JetFlavor::code(const int iJet) {
  // The flavor code of the jet
  JetFlavor::Code flav = UNDEFINED;

  // Check if parton was matched to jet and
  // if this parton has valid pdg id and status
  //    if( isJetMatchedPartonExist[ iJet ] && validPdgID( flavorJetMatchedParton[ iJet ] ) ) {
  if( ( isJetMatchedPartonExist[iJet] )
      &&
      ( statusJetMatchedParton[iJet] == 2 || statusJetMatchedParton[iJet] == 3 )
      &&
      ( std::abs(partonFlavorJet[iJet]) != 0 ) ) {
      
    // A valid parton has been matched. Now apply our flavor-code definition.
    // First, check if jet contains 1 or more b quarks
    if( static_cast<int>(BContentJet3[ iJet ]) > 0 ) {
      // Exactly 1 or >1 b quark(s)?
      if( static_cast<int>(BContentJet3[ iJet ]) == 1 ) flav = B;
      else                                              flav = BB;
    }
    // Otherwise, check if jet contains 1 or more c quarks
    else if( static_cast<int>(CContentJet3[ iJet ]) > 0 ) {
      // Exactly 1 or >1 c quark(s)?
      if( static_cast<int>(CContentJet3[ iJet ]) == 1 ) flav = C;
      else                                              flav = CC;
    }
    // Everything else is assumed to be light flavor
    else {
      flav = UDSG;
    }
  }
    
  return flav;
}



// Return dijet flavor code. The dijet flavor code maps
// the flavor code of both jets to a unique number

// TODO if still needed, this should also become a readable enum
int JetFlavor::diJetCode(const int iJet0, const int iJet1) {
  // determine dijet flavor from MC info
  // we set a flavor pair code that does not distinguish ordering
  int theFlav[2];
  // dijetCode uses this (outdated?) definition...
  for(unsigned int i = 0; i < 2; ++i) {
    unsigned int iJet = (i==0 ? iJet0 : iJet1);
    int theFlavCode = -1;
    int theFlavor = partonFlavorJet[iJet];
    if (hflContentJet[iJet] != 0) theFlavor = hflContentJet[iJet];
    switch (abs(theFlavor)) {
    case 0:
    case 1:
    case 2:
    case 3:
      theFlavCode = 0;
      break;
    case 4:
      theFlavCode = 1;
      break;
    case 5:
      theFlavCode = 2;
      break;
    case 21:
      theFlavCode = 0;
      break;
    default:
      std::cout << "bad flavor " << theFlavor << std::endl;
    }
    theFlav[i] = theFlavCode;
  }

  // now set the dijet flavor
  int theFcDijet = -1;
  int flavPatternDijet = 10*theFlav[0]+theFlav[1];
  switch (flavPatternDijet) {
  case 0:
    theFcDijet = 5;
    break;
  case 11:
    theFcDijet = 3;
    break;
  case 22:
    theFcDijet = 0;
    break;
  case 1:
  case 10:
    theFcDijet = 4;
    break;
  case 2:
  case 20:
    theFcDijet = 2;
    break;
  case 12:
  case 21:
    theFcDijet = 1;
    break;
  default:
    std::cout << "diJetFlavorCode: Bad flavor code " << theFlav[0] << " " << theFlav[1] << std::endl;
  }
    
  return theFcDijet;
}



// for backward compatibility; if you have to use this
// method your code is outdated
JetFlavor::Code JetFlavor::fromInt(int flav) {
  if     ( flav == -1 ) return UNDEFINED;
  else if( flav ==  0 ) return UDSG;
  else if( flav ==  1 ) return C;
  else if( flav ==  2 ) return B;
  else if( flav ==  3 ) return CC;
  else if( flav ==  4 ) return BB;
  else                  return UNDEFINED;
}

// for backward compatibility; if you have to use this
// method your code is outdated
JetFlavor::Code JetFlavor::fromString(const std::string &flav) {
  if     ( flav == "udsg" ) return UDSG;
  else if( flav ==  "c"   ) return C;
  else if( flav ==  "b"   ) return B;
  else if( flav ==  "cc"  ) return CC;
  else if( flav ==  "bb"  ) return BB;
  else                      return UNDEFINED;
}


// Returns true for udscb quarks and gluons,
// false otherwise
bool JetFlavor::validPdgID(int pdgID) {
  int id = std::abs(pdgID);
  return id==1 || id==2 || id==3 || id==4 || id==5 || id==21;
}


#endif
