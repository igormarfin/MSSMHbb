#include "stringUtils.h"
#include <iostream>
#include <sstream>
#include <iomanip>

using namespace std;

TString colorNorm = "\033[0m";

stringResiduals::stringResiduals(const TString& full, bool m_initialPath_regex){
  TString part1, part2, res1, res2;

  TString splitter = "";
  if(!m_initialPath_regex){
    int star = full.Index("*");
    int question = full.Index("?");
    if(star < 0 ) star = 1000;
    if(question < 0) question = 1000;
    if( star < question ) { hasWildCard = splitString(full, part1, part2, '*'); splitter = "*"; }
    else { hasWildCard = splitString(full, part1, part2, '?'); splitter = '?'; }
  }
  else {
    int dot = (full.Index(".") >= 0) ? full.Index(".") : 1000;
    int star = (full.Index("*") >= 0) ? full.Index("*") : 1000;
    int braket = (full.Index("[") >= 0) ? full.Index("[") : 1000;
    int plus = (full.Index("+") >= 0) ? full.Index("+") : 1000;
    int pipe = (full.Index("|") >= 0) ? full.Index("|") : 1000;
    int slash_w = (full.Index("\\w") >= 0) ? full.Index("\\w") : 1000;
    if( dot < star && dot < braket && dot < plus && dot < pipe && dot < slash_w) { hasWildCard = splitString(full,part1,part2,'.'); splitter = "."; }
    else if( star < braket && star < plus && star < pipe && star < slash_w ) { hasWildCard = splitString(full,part1,part2,'*'); splitter = "*"; }
    else if( braket < plus && braket < pipe && braket < slash_w) { hasWildCard = splitString(full,part1,part2,'['); splitter = "["; }
    else if( plus < pipe && plus < slash_w) { hasWildCard = splitString(full,part1,part2,'+'); splitter = "+"; }
    else if( pipe < slash_w) { hasWildCard = splitString(full,part1,part2,'|'); splitter = "|"; }
    else if( slash_w >= 0) { hasWildCard = splitString(full,part1,part2,"\\w"); splitter = "\\w"; }
  }

  if(full.EndsWith("/")) isDirectory = true;
  else isDirectory = false;

  truncateString(part1, res1, '/');
  extractStringBeforePoint(part2, res2, '/');

  if( !full.Contains(splitter) && res2.IsNull() && !isDirectory) splitter = "";

  firstPart = part1;
  secondPart = res2;
  wildcard = res1+splitter+part2;

  if(!m_initialPath_regex){
    wildcard.ReplaceAll("*",".*");
    wildcard.ReplaceAll("?",".");
    if(wildcard==".") wildcard = "";
    //  while(wildcard.Contains("..")) wildcard.ReplaceAll("..",".");//protect against bad wild-cards
  }
  
}

stringResiduals::~stringResiduals(){}

void stringResiduals::print(){
  cout << "******************************" << endl;
  cout << "hasWildCard    " << hasWildCard << endl;
  cout << "firstPart:     " << firstPart << endl;
  cout << "secondPart:    " << secondPart << endl;
  cout << "wildcard:      " << wildcard << endl;
  cout << "isDirectory:   " << isDirectory << endl;
  cout << "******************************\n" << endl;
}

void truncateString( TString& full, TString& residual, const char point){
  int i = full.Sizeof();
  bool foundPoint = false;
  for(; i != 0; --i)  if (full(i,1)==point) {foundPoint = true; break;}
  if(foundPoint){
    residual = full(i+1,full.Sizeof()-i-1);
    full = full(0, i+1);
  }
  else {
    residual = full(i,full.Sizeof()-i);
    full = full(0, i);
  }    
}
void truncateString( TString& full, const char point){
  int i = full.Sizeof();
  bool foundPoint = false;
  for(; i != 0; --i)  if (full(i,1)==point) {foundPoint = true; break;}
  if(foundPoint){
    full = full(0, i+1);
  }
  else {
    full = full(0, i);
  }      
}
void extractStringBeforePoint( TString& full, TString& residual, const char point){
  int i = 0;
  for(; i != full.Sizeof(); ++i)  if (full(i,1)==point) break;
  residual = full(i+1, full.Sizeof()-i-1);
  full = full(0, i);
}

bool splitString(const TString& full, TString& part1, TString& part2, const char split){
  int i = full.Index(split);
  if(i >= 0) {
    part1 = full(0, full.Index(split));
    part2 = full(full.Index(split)+1,full.Sizeof() - full.Index(split));
    return true;
  }
  else {
    part1 = full;
    part2 = "";
    return false;
  }
}

bool splitString(const TString& full, TString& part1, TString& part2, const TString & split){
  int i = full.Index(split);
  if(i >= 0) {
    part1 = full(0, full.Index(split));
    part2 = full(full.Index(split)+split.Length(),full.Sizeof() - full.Index(split) - split.Length() );
    return true;
  }
  else {
    part1 = full;
    part2 = "";
    return false;
  }
}

vector<TString> parseList(const TString& str, char c){
  vector<TString> list;
  TString full = str;
  TString part1;
  TString part2;

  while( splitString(full, part1, part2, c) ){
    list.push_back(part1);
    full = part2;
  }
  list.push_back(full);
  return list;
}

set<TString> parseSet(const TString& str, char c){
  set<TString> list;
  TString full = str;
  TString part1;
  TString part2;

  while( splitString(full, part1, part2, c) ){
    list.insert(part1);
    full = part2;
  }
  list.insert(full);
  return list;
}

vector<TString> parseList(const TString& str, const TString& c){
  vector<TString> list;
  TString full = str;
  TString part1;
  TString part2;

  while( splitString(full, part1, part2, c) ){
    list.push_back(part1);
    full = part2;
  }
  list.push_back(full);
  return list;
}

TString returnString(int a){
  stringstream ss; ss << a;
  TString b; ss >> b;
  return b;
}

TString returnString(double a, int precision){
  stringstream ss; ss << setprecision(precision) << a;
  TString b; ss >> b;
  return b;
}

void printPad(int level, int space, TString symbol, TString* color){
  if(color){
    for(int i = 0; i!=level*space; ++i) cout << *color << symbol << colorNorm;   
  }
  else {
    for(int i = 0; i!=level*space; ++i) cout << symbol ;   
  }

}

TString makePad(int level, int space){
  TString pad;
  if (space < 0) return "";
  for(int i = 0; i!=level*space; ++i) pad+= " ";   
  return pad;
}

void addPad(TString &s, int space){ for(int i = 0; i!=space; ++i) s+=" ";}

void insertString( TString & line, TString name, TString expr, int pad, int level ){
  name+=makePad(level, std::max(20-name.Length(), pad));
  line.ReplaceAll("\\c", name );  
}


