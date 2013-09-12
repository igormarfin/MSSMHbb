#ifndef HbbTrigger_h
#define HbbTrigger_h

#include <vector>
#include <TF1.h>
#include <TMath.h>

class HbbTrigger {

public:
   HbbTrigger(int );
   // path = 0, HLT_Jet46_Jet38
   // path = 1, HLT_Jet46_Jet38_Jet20
   // path = 2, HLT_Jet60_Jet53
 
   ~HbbTrigger();

   float getEfficiency(float , float );
   float getEfficiency(float , float , float );
   void  setPeriod(int);   // Run2011A = 0; Run2011B = 1

private:
   TF1 * _f_eff[5]; // 0 = HLT20; 1 = HLT38; 2 = HLT46; 3 = HLT53; 4 = HLT60

   float _path;
   int _period;
 
};
#endif
