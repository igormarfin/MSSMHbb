#ifndef LSROOT_UTILS_H
#define LSROOT_UTILS_H

class TKey;
class TH1;
class TH2;

#include <TObject.h>

#include <vector>

struct histDisplayOpt {
  histDisplayOpt(int loglv, int stat){
    loglevel = loglv;
    statOpt = stat;
  }
  int loglevel;
  int statOpt;
};

struct asciDrawOptions {
  asciDrawOptions(TH1* h, double ratio);
  asciDrawOptions();
  asciDrawOptions* clone();
  double minx;
  double miny;
  double minz;
  double maxx;
  double maxy;
  double maxz;
  bool logy;
  bool logz;  
  void print();
};

class JKey {
 public:
  JKey( TObject* o, TKey* k, TDirectory* d);
  ~JKey();
  TDirectory* dir;
  TObject* obj;
  TKey* key;
};

void getFloat(const TString& type, double &min, double &max);
TH1* extractHist(TH1* h, asciDrawOptions* opt);
TH2* extractHist2(TH2* h, asciDrawOptions* opt);

bool jkeyCompTime(const JKey* a, const JKey* b);
bool jkeyComp(const JKey* a, const JKey* b);
bool jkeyCompR(const JKey* a, const JKey* b);
bool jkeyCompSize(const JKey* a, const JKey* b);
bool jkeyCompSizeR(const JKey* a, const JKey* b);
bool jkeyCompN(const JKey* a, const JKey* b);
bool jkeyCompNR(const JKey* a, const JKey* b);
bool jkeyCompNBins(const JKey* a, const JKey* b);
bool jkeyCompNBinsR(const JKey* a, const JKey* b);

#endif
