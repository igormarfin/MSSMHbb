#ifndef FindLine_C_H
#define FindLine_C_H
///run me
/// findLine.C(TString arg1,TString arg2,TString arg3,TString pathArg="./",Bool_t verbose=kFALSE)


#include "TString.h"
#include "TSystem.h"
#include <iostream>


static TString findLine(TString arg1,TString arg2,TString arg3,TString pathArg="./",Bool_t verbose=kFALSE){
TString cmd=TString("source ");
cmd+=pathArg+TString("utils.sh; ");
cmd+=TString("findLine " );
cmd+=TString(' ')+arg1+TString(' ');
cmd+=TString(' ')+arg2+TString(' ');
cmd+=TString(' ')+arg3+TString(' ');
TString res =gSystem->GetFromPipe(cmd.Data()) ;
if (verbose) std::cout<<res<<std::endl;
return res;}
#endif
