#include "Analysis/Utilities/interface/HbbTrigger.h"

HbbTrigger::HbbTrigger(int path)
{
   if ( path < 0 || path > 2 )
   {
      std::cout << "Path not recognised. Setting path to 0 (HLT_46_38)." << std::endl;
      path = 0;
   }
   
   _f_eff[0] = new TF1("f_eff_0","[0]*0.5*(1+TMath::Erf((x-[2])/[1]))", 20.,10000.);
   _f_eff[1] = new TF1("f_eff_1","[0]*0.5*(1+TMath::Erf((x-[2])/[1]))", 20.,10000.);
   _f_eff[2] = new TF1("f_eff_2","[0]*0.5*(1+TMath::Erf((x-[2])/[1]))", 20.,10000.);
   _f_eff[3] = new TF1("f_eff_3","[0]*0.5*(1+TMath::Erf((x-[2])/[1]))", 20.,10000.);
   _f_eff[4] = new TF1("f_eff_4","[0]*0.5*(1+TMath::Erf((x-[2])/[1]))", 20.,10000.);
   
   _period = 1;
   this -> setPeriod(_period);
//   std::cout << " *** HbbTrigger::HbbTrigger(int)" << std::endl;
//   std::cout << "     Period set to default value = 1 (Run2011B)" << std::endl;
//   std::cout << "     To set for Run2011A use the method: setPeriod(0)" << std::endl;
   
   _path = path;
  
}

HbbTrigger::~HbbTrigger() {

}

float HbbTrigger::getEfficiency(float jet1Pt, float jet2Pt)
{
   float eff1 = 0;
   float eff2 = 0;
   float eff  = 0;
   if ( _path == 0 ) // HLT_46_38
   {
      eff1 = _f_eff[2] -> Eval(jet1Pt);
      eff2 = _f_eff[1] -> Eval(jet2Pt);
//       std::cout << "HbbTrigger: pt1= " << jet1Pt << " eff1=" << eff1
// 		<< " pt2=" << jet2Pt 
// 		<< " eff2=" << eff2 << std::endl;
   }
   else if ( _path == 2 ) // HLT_60_53
   {
      eff1 = _f_eff[4] -> Eval(jet1Pt);
      eff2 = _f_eff[3] -> Eval(jet2Pt);
   }
   
   eff = eff1*eff2;
   
   return eff;
}

float HbbTrigger::getEfficiency(float jet1Pt, float jet2Pt, float jet3Pt)
{
   float eff1 = 0;
   float eff2 = 0;
   float eff3 = 0;
   float eff  = 0;
   if ( _path == 1 ) // HLT_46_38_20
   {
      eff1 = _f_eff[2] -> Eval(jet1Pt);
      eff2 = _f_eff[1] -> Eval(jet2Pt);
      eff3 = _f_eff[0] -> Eval(jet3Pt);
      eff  = eff1*eff2*eff3;
//       std::cout << "HbbTrigger: pt1= " << jet1Pt << " eff1=" << eff1
// 		<< " pt2=" << jet2Pt <<  " eff2=" << eff2 
//      << " pt3=" << jet3Pt <<  " eff3=" << eff3 << std::endl;
//       std::cout << "PERIOD: " << _period << std::endl;
   }
   else
   {
      eff = this -> getEfficiency(jet1Pt,jet2Pt);
      
   }
   
   return eff;
}
void HbbTrigger::setPeriod(int newPeriod)
{
   if ( _period != newPeriod )
   {
      std::cout << "*** HbbTrigger::setPeriod(int) : " << "  Data period change from " << _period << " to " << newPeriod << std::endl;
   }
   
   _period = newPeriod;
   if ( _period < 0 || _period > 1 ) 
   {
      std::cout << "Unknown period " << _period << ". Setting to default period = 1 (Run2011B)" << std::endl;
      _period = 1;
   }
   float  par[3][5]; 
   float epar[3][5];

    par[0][0] = 0.981504; 
    par[1][0] = 19.5766; 
    par[2][0] = 11.9996; 
    if ( _period == 0 )
    {
       par[0][0] = 0.982952; 
       par[1][0] = 14.86; 
       par[2][0] = 26.9275; 
    }
    par[0][1] = 0.984353; 
    par[1][1] = 15.162; 
    par[2][1] = 42.5712; 
    par[0][2] = 0.98709; 
    par[1][2] = 14.9056; 
    par[2][2] = 44.6875; 
    par[0][3] = 0.983103; 
    par[1][3] = 14.7521; 
    par[2][3] = 52.9569; 
    par[0][4] = 0.985096; 
    par[1][4] = 15.746; 
    par[2][4] = 58.3774; 
    
    if ( _period == 0 )
    {
       epar[0][0] = 0.00234037; 
       epar[1][0] = 0.15269; 
       epar[2][0] = 0.060582; 
    }
   epar[0][0] = 0.000955436; 
   epar[1][0] = 0.217865; 
   epar[2][0] = 0.110966; 
   epar[0][1] = 0.00198478; 
   epar[1][1] = 0.227409; 
   epar[2][1] = 0.0778294; 
   epar[0][2] = 0.0019397; 
   epar[1][2] = 0.430618; 
   epar[2][2] = 0.207463; 
   epar[0][3] = 0.00281911; 
   epar[1][3] = 0.425887; 
   epar[2][3] = 0.152465; 
   epar[0][4] = 0.00322333; 
   epar[1][4] = 0.857611; 
   epar[2][4] = 0.426579; 
   
   for ( int i = 0; i < 5 ; ++i )
   {
      _f_eff[i] -> SetParameter(0,par[0][i]);
      _f_eff[i] -> SetParameter(1,par[1][i]);
      _f_eff[i] -> SetParameter(2,par[2][i]);
   }

}
