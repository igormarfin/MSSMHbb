#ifndef SystReadInput_h
#define SystReadInput_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TString.h>


void SystReadInput(int & systype,double & variation)
{

int isFile3 = (int) gSystem->Exec(TString("[ -f systematics")+TString(" ]"));

if (isFile3==0) {
	ifstream ff3("systematics"); ///must exists
	TString str3;
	str3.ReadLine(ff3);


	while (str3.Length()>0 )
	{
		TObjArray * arr3 = str3.Tokenize(TString(':'));
		if (arr3->GetEntries()==2 && !str3.Contains("#")) {
			systype=TString(arr3->At(0)->GetName()).Atoi();
			variation=TString(arr3->At(1)->GetName()).Atof();
		}

	str3.ReadLine(ff3);

	} ///while

} ///if isFile3





return;

}




#endif

