#include "Analysis/Utilities/interface/HbbXsec.h"

HbbXsec::HbbXsec(std::string file) {
//   std::string file = "/data/group/higgs/bbbar/HbbMSSMAnalysis/TheoryXsectionScans/out.mhmax_mu200_7_nnlo.tanBeta_gte1.root";
   char * cfile;
   cfile = new char [file.size()+1];
   strcpy (cfile, file.c_str());
   mssm_xs.SetInput(cfile);
}

HbbXsec::~HbbXsec() {
}

double HbbXsec::getXsecBR(double mA, double tanb)
{
   // Branching ratio   
   double br_Abb  = mssm_xs.Give_BR_A_bb(mA,tanb);
   double br_hbb  = mssm_xs.Give_BR_h_bb(mA,tanb);
   double br_Hbb  = mssm_xs.Give_BR_H_bb(mA,tanb);
   
   // Cross section (5-flavour scheme)
   double xSecbbA = mssm_xs.Give_Xsec_bbA5f(mA,tanb);
   double xSecbbh = mssm_xs.Give_Xsec_bbh5f(mA,tanb);
   double xSecbbH = mssm_xs.Give_Xsec_bbH5f(mA,tanb);
   
   // 
   
   double xs = 0;
   
   if ( mA < 130 ) 
   {
      xs = xSecbbA*br_Abb + xSecbbh*br_hbb;
   } else if ( mA == 130 )
   {
      xs = xSecbbA*br_Abb + xSecbbh*br_hbb + xSecbbH*br_Hbb;
   } else
   {
      xs = xSecbbA*br_Abb + xSecbbH*br_Hbb;
   }
   
   
   return xs;
}
double HbbXsec::getXsecBRUnc_muup(double mA, double tanb)
{
   // Branching ratio   
   double br_Abb  = mssm_xs.Give_BR_A_bb(mA,tanb);
   double br_hbb  = mssm_xs.Give_BR_h_bb(mA,tanb);
   double br_Hbb  = mssm_xs.Give_BR_H_bb(mA,tanb);
   
   // Total scale uncertainty
   double xSecUnc_bbA = mssm_xs.Give_XsecUnc_muup_bbA5f(mA,tanb)*br_Abb;
   double xSecUnc_bbh = mssm_xs.Give_XsecUnc_muup_bbh5f(mA,tanb)*br_hbb;
   double xSecUnc_bbH = mssm_xs.Give_XsecUnc_muup_bbH5f(mA,tanb)*br_Hbb;
   
   // 
   double xsUnc = _calculateUncertainty(mA,xSecUnc_bbA,xSecUnc_bbh,xSecUnc_bbH);

   return xsUnc;
}

double HbbXsec::getXsecBRUnc_mudown(double mA, double tanb)
{
   // Branching ratio   
   double br_Abb  = mssm_xs.Give_BR_A_bb(mA,tanb);
   double br_hbb  = mssm_xs.Give_BR_h_bb(mA,tanb);
   double br_Hbb  = mssm_xs.Give_BR_H_bb(mA,tanb);
   
   // Total scale uncertainty
   double xSecUnc_bbA = mssm_xs.Give_XsecUnc_mudown_bbA5f(mA,tanb)*br_Abb;
   double xSecUnc_bbh = mssm_xs.Give_XsecUnc_mudown_bbh5f(mA,tanb)*br_hbb;
   double xSecUnc_bbH = mssm_xs.Give_XsecUnc_mudown_bbH5f(mA,tanb)*br_Hbb;
   
   // 
   
   double xsUnc = _calculateUncertainty(mA,xSecUnc_bbA,xSecUnc_bbh,xSecUnc_bbH);
   
    return xsUnc;
}

double HbbXsec::getXsecBRUnc_pdfalphas68up(double mA, double tanb)
{
   // Branching ratio   
   double br_Abb  = mssm_xs.Give_BR_A_bb(mA,tanb);
   double br_hbb  = mssm_xs.Give_BR_h_bb(mA,tanb);
   double br_Hbb  = mssm_xs.Give_BR_H_bb(mA,tanb);
   
   // Total pdf+alphas uncertainty
   double xSecUnc_bbA = mssm_xs.Give_XsecUnc_pdfalphas68up_bbA5f(mA,tanb)*br_Abb;
   double xSecUnc_bbh = mssm_xs.Give_XsecUnc_pdfalphas68up_bbh5f(mA,tanb)*br_hbb;
   double xSecUnc_bbH = mssm_xs.Give_XsecUnc_pdfalphas68up_bbH5f(mA,tanb)*br_Hbb;
   
   // 
   
   double xsUnc = _calculateUncertainty(mA,xSecUnc_bbA,xSecUnc_bbh,xSecUnc_bbH);
   
   
   return xsUnc;
}

double HbbXsec::getXsecBRUnc_pdfalphas68down(double mA, double tanb)
{
   // Branching ratio   
   double br_Abb  = mssm_xs.Give_BR_A_bb(mA,tanb);
   double br_hbb  = mssm_xs.Give_BR_h_bb(mA,tanb);
   double br_Hbb  = mssm_xs.Give_BR_H_bb(mA,tanb);
   
   // Total pdf+alphas uncertainty
   double xSecUnc_bbA = mssm_xs.Give_XsecUnc_pdfalphas68down_bbA5f(mA,tanb)*br_Abb;
   double xSecUnc_bbh = mssm_xs.Give_XsecUnc_pdfalphas68down_bbh5f(mA,tanb)*br_hbb;
   double xSecUnc_bbH = mssm_xs.Give_XsecUnc_pdfalphas68down_bbH5f(mA,tanb)*br_Hbb;
   
   // 
   
   double xsUnc = _calculateUncertainty(mA,xSecUnc_bbA,xSecUnc_bbh,xSecUnc_bbH);
   
   return xsUnc;
}

double HbbXsec::getXsecBRUnc_UEPS(double mA, double tanb)
{
   double xsUnc = this -> getXsecBR(mA,tanb)*0.04;
   
   return xsUnc;
}

double HbbXsec::_calculateUncertainty(double mA, double uncA, double unch, double uncH )
{
   double xsUnc = 0;
   
   if ( mA < 130 ) 
   {
      xsUnc = sqrt(uncA*uncA + unch*unch);
   } else if ( mA == 130 )
   {
      xsUnc = sqrt(uncA*uncA + unch*unch + uncH*uncH);
   } else
   {
      xsUnc = sqrt(uncA*uncA + uncH*uncH);
   }
   
   return xsUnc;
}
