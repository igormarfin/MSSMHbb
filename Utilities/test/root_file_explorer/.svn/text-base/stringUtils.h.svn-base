#ifndef STRING_UTILS_H
#define STRING_UTILS_H

//Utility Class for lsRoot and hadd
#include <TString.h>
#include <vector>
#include <set>

bool splitString(const TString& full, TString& part1, TString& part2, const char split);//returns true if split was found                                                                                                                 
bool splitString(const TString& full, TString& part1, TString& part2, const TString &split);//returns true if split was found                                                                                                             
std::vector<TString> parseList(const TString& str, char c);
std::set<TString> parseSet(const TString& str, char c);
std::vector<TString> parseList(const TString& str, const TString& c);

void truncateString( TString& full, TString& residual, const char point);
void truncateString( TString& full, const char point);
void extractStringBeforePoint( TString& full, TString& residual, const char point);

TString returnString( int );
TString returnString( double, int );

//String formatting
void printPad(int level, int space=10, TString symbol= " ", TString *color = 0);
TString makePad(int level, int space=10);
void addPad(TString &s, int space);

void insertString( TString & target, TString insert, TString expr, int pad = 1, int level = 1 );

//Class to assist in wild-carding directories
class stringResiduals {
 public:
  stringResiduals(const TString& full, bool m_initialPath_regex = false);
  ~stringResiduals();
  void print();

  bool hasWildCard;
  bool isDirectory;
  //Path before and after wildcard
  TString firstPart;
  TString secondPart;

  TString wildcard;
};

#endif
