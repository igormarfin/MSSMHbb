#ifndef B_TAG_H
#define B_TAG_H

#include <iostream>
#include <vector>
#include <string>
#include <exception>



// Defines offline b-tag discriminants and working points and
// stores related information such as the corresponding discriminants
// and cut values.
class BTag {
public:

  // Offline b-tag discriminants
  enum Discriminant { TCHP, CSV, SSVHP };
  // Offline b-tag working points
  enum WP { TCHPT, TCHP6, CSVT, SSVHPT };

  
  // For convenience
  typedef std::vector<WP> WPs;
  typedef std::vector<WP>::const_iterator WPsIt;


  // Returns vector of available offline b-tag working points
  static BTag::WPs getList();

  // Returns vector of offline b-tag working points as
  // specified in the string
  static BTag::WPs getList(const std::vector<std::string> &vecstr);


  // Returns cut value on discriminant for the
  // given working point.
  static double cut(BTag::WP wp);
  
  // Returns the discriminant of the given
  // working point
  static BTag::Discriminant discriminant(BTag::WP wp);

  // std::string representation of working point
  static std::string toString(BTag::WP wp);

  // std::string representation of discriminant
  static std::string toString(BTag::Discriminant dscr);

  // Return working point for given string representation.
  // This method is public only for backward compatibility;
  // if you have to use it your code is outdated!
  static BTag::WP fromString(const std::string &str);

};



// Returns vector of available offline b-tag working points
BTag::WPs BTag::getList() { 
  WPs wps; 
  wps.push_back(TCHPT); 
  wps.push_back(TCHP6); 
  wps.push_back(CSVT); 
  wps.push_back(SSVHPT); 
  return wps;
}

// Returns vector of offline b-tag working points as
// specified in the string
BTag::WPs BTag::getList(const std::vector<std::string> &vecstr) {
  WPs wps; 
  for(std::vector<std::string>::const_iterator it = vecstr.begin(); 
      it != vecstr.end(); it++) { 
    wps.push_back(BTag::fromString((*it))); 
  } 
  return wps; 
} 

    
// Returns cut value on discriminant for the
// given working point.
double BTag::cut(BTag::WP wp) {
  double val = -99999.;

  if     ( wp == TCHPT  ) val = 3.41;
  else if( wp == TCHP6  ) val = 6.;
  else if( wp == CSVT   ) val = 0.898;
  else if( wp == SSVHPT ) val = 2.;

  return val;
}

// Returns the discriminant of the given
// working point
BTag::Discriminant BTag::discriminant(BTag::WP wp) {
  Discriminant discr = TCHP;

  if     ( wp == TCHPT  ) discr = TCHP;
  else if( wp == TCHP6  ) discr = TCHP;
  else if( wp == CSVT   ) discr = CSV;
  else if( wp == SSVHPT ) discr = SSVHP;
  else {
    throw std::exception();
  }

  return discr;
}


// std::string representation of working point
std::string BTag::toString(BTag::WP wp) {
  std::string str = "undefined";

  if     ( wp == TCHPT  ) str = "TCHPT";
  else if( wp == TCHP6  ) str = "TCHP6";
  else if( wp == CSVT   ) str = "CSVT";
  else if( wp == SSVHPT ) str = "SSVHPT";

  return str;
}

// std::string representation of discriminant
std::string BTag::toString(BTag::Discriminant dscr) {
  std::string str = "undefined";

  if     ( dscr == TCHP  ) str = "TCHP";
  else if( dscr == CSV   ) str = "CSV";
  else if( dscr == SSVHP ) str = "SSVHP";
  else throw std::exception();

  return str;
}

// Return working point for given string representation.
// This method is public only for backward compatibility;
// if you have to use it (publicly) your code is outdated!
BTag::WP BTag::fromString(const std::string &str) {
  if     ( str == "TCHPT"  ) return TCHPT;
  else if( str == "TCHP6"  ) return TCHP6;
  else if( str == "CSVT"   ) return CSVT;
  else if( str == "SSVHPT" ) return SSVHPT;
  else {
    std::cout << "@BTag: error: BTagger '"<< str << "' not known!" << std::endl;
    throw std::exception();
    return SSVHPT;
  }
}

#endif
