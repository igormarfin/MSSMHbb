#include "Analysis/Utilities/interface/HbbXsec.h"

//void TripleBtagAnalysis_SF() {
int main(int narg,char** varg) {

   HbbXsec hbb_xs("/data/group/higgs/bbbar/HbbMSSMAnalysis/TheoryXsectionScans/out.mhmax_mu200_7_nnlo.tanBeta_gte1.root");
   
   std::cout << "Total xsec*BR                   = " << hbb_xs.getXsecBR(120,20)                    << std::endl;
   std::cout << "Total xsec*BR Unc MuUP          = " << hbb_xs.getXsecBRUnc_muup(120,20)            << std::endl;
   std::cout << "Total xsec*BR Unc MuDW          = " << hbb_xs.getXsecBRUnc_mudown(120,20)          << std::endl;
   std::cout << "Total xsec*BR Unc PdfAlphas68UP = " << hbb_xs.getXsecBRUnc_pdfalphas68up(120,20)   << std::endl;
   std::cout << "Total xsec*BR Unc PdfAlphas68DW = " << hbb_xs.getXsecBRUnc_pdfalphas68down(120,20) << std::endl;
   std::cout << "Total xsec*BR Unc UE-PS         = " << hbb_xs.getXsecBRUnc_UEPS(120,20)            << "    -  This is simply the xsec*BR*0.04. " << std::endl;

}
