#ifndef Field_C_H
#define Field_C_H
///run me
/// field.C(TString arg1,TString arg2,TString pathArg="./",Bool_t verbose=kFALSE)


#include "TString.h"
#include "TSystem.h"
#include <iostream>


static TString field(TString arg1,TString arg2,TString pathArg="./",Bool_t verbose=kFALSE){
TString cmd=TString("source ");
cmd+=pathArg+TString("utils.sh; ");
cmd+=TString("field " );
cmd+=TString(' ')+arg1+TString(' ');
cmd+=TString(' ')+arg2+TString(' ');
TString res =gSystem->GetFromPipe(cmd.Data()) ;
if (verbose) std::cout<<res<<std::endl;
return res;}
#endif
