#ifndef FLAVOR_CLASS_H
#define FLAVOR_CLASS_H

#include <string>


// Definitions of flavor classes for double and triple b-tag samples
class FlavorClass {
public:

  // Flavor classes for double b-tag samples
  enum Dijet { B1B1, B1C1, B1Q, C1C1, C1Q, QQ, C2B1, C2C1, C2Q, C2C2, B2B1, B2C1, B2Q, B2B2, B2C2 };
  enum Jet3 { Q, C, B };
  // Flavor classes for triple b-tag sample
  enum Trijet { BBB, BBC, BBQ, BCB, BQB, NON_BB };

  // Return list of FlavorClass values used in the analysis
  static std::vector<FlavorClass::Dijet> getListDijet();
  static std::vector<FlavorClass::Jet3> getListJet3();
  static std::vector<FlavorClass::Trijet> getListTrijet();

  // std::string representations of FlavorClass values
  static std::string toString(FlavorClass::Dijet fc);
  static std::string toString(FlavorClass::Jet3 fc);
  static std::string toString(FlavorClass::Trijet fc);
};



std::vector<FlavorClass::Dijet> FlavorClass::getListDijet() {
  std::vector<FlavorClass::Dijet> fcs;
  fcs.push_back(B1B1);
  fcs.push_back(B1C1);
  fcs.push_back(B1Q);
  fcs.push_back(C1C1);
  fcs.push_back(C1Q);
  fcs.push_back(QQ);
  fcs.push_back(C2B1);
  fcs.push_back(C2C1);
  fcs.push_back(C2Q);
  fcs.push_back(C2C2);
  fcs.push_back(B2B1);
  fcs.push_back(B2C1);
  fcs.push_back(B2Q);
  fcs.push_back(B2B2);
  fcs.push_back(B2C2);

  return fcs;
}


std::vector<FlavorClass::Jet3> FlavorClass::getListJet3() {
  std::vector<FlavorClass::Jet3> fcs;
  fcs.push_back(Q);
  fcs.push_back(C);
  fcs.push_back(B);

  return fcs;
}


std::vector<FlavorClass::Trijet> FlavorClass::getListTrijet() {
  std::vector<FlavorClass::Trijet> fcs;
  fcs.push_back(BBB);
  fcs.push_back(BBC);
  fcs.push_back(BBQ);
  fcs.push_back(BCB);
  fcs.push_back(BQB);
  fcs.push_back(NON_BB);
  
  return fcs;
}


std::string FlavorClass::toString(FlavorClass::Dijet fc) {
  std::string str = "undefined";
  if(      fc == B1B1 ) str = "b1b1";
  else if( fc == B1C1 ) str = "b1c1";
  else if( fc == B1Q  ) str = "b1q";
  else if( fc == C1C1 ) str = "c1c1";
  else if( fc == C1Q  ) str = "c1q";
  else if( fc == QQ   ) str = "qq";
  else if( fc == C2B1 ) str = "c2b1";
  else if( fc == C2C1 ) str = "c2c1";
  else if( fc == C2Q  ) str = "c2q";
  else if( fc == C2C2 ) str = "c2c2";
  else if( fc == B2B1 ) str = "b2b1";
  else if( fc == B2C1 ) str = "b2c1";
  else if( fc == B2Q  ) str = "b2q";
  else if( fc == B2B2 ) str = "b2b2";
  else if( fc == B2C2 ) str = "b2c2";

  return str;
}


std::string FlavorClass::toString(FlavorClass::Jet3 fc) {
  std::string str = "undefined";
  if(      fc == Q ) str = "q";
  else if( fc == C ) str = "c";
  else if( fc == B ) str = "b";

  return str;
}


std::string FlavorClass::toString(FlavorClass::Trijet fc) {
  std::string str = "undefined";
  if(      fc == BBB    ) str = "bbb";
  else if( fc == BBC    ) str = "bbc";
  else if( fc == BBQ    ) str = "bbq";
  else if( fc == BCB    ) str = "bcb";
  else if( fc == BQB    ) str = "bqb";
  else if( fc == NON_BB ) str = "non-bb";

  return str;
}    

#endif
