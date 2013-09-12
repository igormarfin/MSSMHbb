#ifndef RelativeOnlineSFandUncertainties_h
#define RelativeOnlineSFandUncertainties_h

#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <exception>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <TF1.h>
#include <TROOT.h>

/* fit result for scaling factor and uncertainties: control_btageff_tagandprobe_normal_btag_eff_tagprobe_data_vs_mc_merged.eps */
/* offset: 1.046605471563083e+00 */
/* slope: 9.138538654640538e-05 */
/* cov00 2.713675092615117e-03 */
/* cov11 1.409894277554858e-07 */
/* cov01 -1.730917795757378e-05 */
/* end of fit result */


/* fit result for scaling factor and uncertainties: control_btageff_tagandprobe_loose_btag_eff_tagprobe_data_vs_mc_merged.eps */
/* offset: 1.049022182824710e+00 */
/* slope: -1.779993975647544e-04 */
/* cov00 4.810413688794324e-04 */
/* cov11 1.574249187396761e-08 */
/* cov01 -2.418520028587507e-06 */
/* end of fit result */


class RelativeOnlineSFandUncertainties {
public:
  RelativeOnlineSFandUncertainties() : 
    def_error_function("sqrt([0]+x*x*[1]+2.0*x*[2])"),
    f_error_b_normal(TF1("f_error_b_normal", def_error_function.c_str(), 0.0, 3500.0)),
    f_error_b_loose(TF1("f_error_b_loose", def_error_function.c_str(), 0.0, 3500.0)),
    SF_b_normal(TF1("SF_b_normal", "pol1", 0.0, 3500.0)),
    SF_b_loose(TF1("SF_b_loose", "pol1", 0.0, 3500.0))
  {
    f_error_b_normal.SetParameter(0, 2.713675092615117e-03);
    f_error_b_normal.SetParameter(1, 1.409894277554858e-07);
    f_error_b_normal.SetParameter(2, -1.730917795757378e-05);
    
    f_error_b_loose.SetParameter(0, 4.810413688794324e-04);
    f_error_b_loose.SetParameter(1, 1.574249187396761e-08);
    f_error_b_loose.SetParameter(2, -2.418520028587507e-06);
    

    SF_b_normal.SetParameter(0,1.046605471563083e+00);
    SF_b_normal.SetParameter(1,9.138538654640538e-05);
    
    SF_b_loose.SetParameter(0,1.049022182824710e+00);
    SF_b_loose.SetParameter(1,-1.779993975647544e-04);
  }
  //flav is a string for the flavour category: e.g. 'b'
  //MC relative online btagging efficiency has to be multiplied with the result of GetSF
  double GetSF(double pt, std::string flav, bool normal)  const //if normal is equal to true, then the 'normal' btagging is used, otherwise the 'loose' btagging (i.e. for the very high mass trigger)
  {
    checkflavpt(pt,flav);  

    const double sf = normal ? SF_b_normal.Eval(pt) : SF_b_loose.Eval(pt);
    if(sf <= 0.0) {
      std::cout << "error: unreasonable relative online scale factor: SF = " << sf << std::endl;
      throw std::exception();
    }
    return 1.0/sf;
  }
  //flav is a string for the flavour category: e.g. 'b'
  double GetUncertainty(double pt, std::string flav, bool normal) const //if normal is equal to true, then the 'normal' btagging is used, otherwise the 'loose' btagging (i.e. for the very high mass trigger)
  {
    checkflavpt(pt,flav);  
    const double error = normal ? f_error_b_normal.Eval(pt) : f_error_b_loose.Eval(pt);
    if(error<=0.0) {
      std::cout << "error: unreasonable uncertainty: uncerainty = " << error << std::endl;
      throw std::exception();
    }
    

    return error;
  }
private:
  void checkflavpt(double pt, std::string &flav) const
  {
    if(flav != "b") {
      std::cout << "error: only b flavour is supported for the relative online SF (uncertainty).exit." << std::endl;
      throw std::exception();
    }
    if(pt <= 0.0 || pt >= 3500.0) {
      std::cout << "error: negative or large (pt>3500 GeV) p_T: pt = " << pt << std::endl;
      throw std::exception();
    }
  }


  const std::string def_error_function;
  TF1 f_error_b_normal, f_error_b_loose;
  TF1 SF_b_normal, SF_b_loose;

};

#endif
