///run me
/// HowManyFieldsInFile.C(TString arg1)


#include "TString.h"
#include "TSystem.h"
#include <iostream>


TString HowManyFieldsInFile(TString arg1){
TString cmd=TString("source utils.sh; ");
cmd+=TString("HowManyFieldsInFile " );
cmd+=TString(' ')+arg1+TString(' ');
TString res =gSystem->GetFromPipe(cmd.Data()) ;
std::cout<<res<<std::endl;
return res;}
