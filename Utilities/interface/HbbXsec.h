#ifndef HbbXsec_h
#define HbbXsec_h

#include <vector>
#include <string>
#include <TF1.h>
#include <TMath.h>

#include "Analysis/Utilities/interface/mssm_xs_tools.h"

class HbbXsec {

public:
   HbbXsec(std::string);
 
   ~HbbXsec();

   double getXsecBR(double , double );   // returns the cross section * BR
   double getXsecBRUnc_muup(double , double );  // returns the cross section* BR scale uncertainty up 
   double getXsecBRUnc_mudown(double , double );  // returns the cross section* BR scale uncertainty down
   double getXsecBRUnc_pdfalphas68up(double , double );  // returns the cross section* BR pdf+alphas uncertainty up
   double getXsecBRUnc_pdfalphas68down(double , double ); // returns the cross section* BR pdf+alphas uncertainty down
   double getXsecBRUnc_UEPS(double , double );  // returns the cross section* BR underlying even+parton shower uncertainty

private:

   mssm_xs_tools mssm_xs;

   double _calculateUncertainty(double, double, double, double);

};
#endif
