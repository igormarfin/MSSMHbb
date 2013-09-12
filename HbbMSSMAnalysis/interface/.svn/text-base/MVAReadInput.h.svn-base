#ifndef MVAReadInput_h
#define MVAReadInput_h

#include "Analysis/Utilities/interface/MVAComputer.h"
#include "Analysis/Utilities/interface/MVATrainer.h"
#include "Analysis/Utilities/interface/MathRoot.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <string>
#include <vector>

#include <TSystem.h>
#include <TString.h>

#include <TROOT.h>
#include <map>
#include <TObjArray.h>
#include <TLorentzVector.h>


std::map<TString,MVAComputer *> _mvaComputersMap;
std::map<TString,MVATrainer *>   _mvaTrainersMap;
std::map<TString,TString >   _mvaCutsMap;
std::map<TString,TString >   _mvaMethodsMap;


std::map<TString,TString> _exclude_trainer_from_selection;


void ExcludeVariableFromMVA()
{

int isFile2 = (int) gSystem->Exec(TString("[ -f exclude_trainer")+TString(" ]"));

if (isFile2==0) {
	ifstream ff2("exclude_trainer"); ///must exists
	TString str2;
	str2.ReadLine(ff2);

	while (str2.Length()>0 )
	{
		if (!str2.Contains("#"))
			_exclude_trainer_from_selection[str2]=str2;
		str2.ReadLine(ff2);

	} ///while
} ///if isFile2


return; 
}


void ReadMVASettings(TString _mvasettings)
{

ifstream ff( _mvasettings.Data());

TString str;
str.ReadLine(ff);

while (str.Length()>0 )
{


///it should be of the form     NAME:METHODNAME:FILE.xml:
///


	TObjArray * arr1 = str.Tokenize(TString(':'));

	if (arr1->GetEntries()>2 && arr1->GetEntries()<5 && !str.Contains("#")) {

	TString _name = arr1->At(0)->GetName();
	TString _methodname = arr1->At(1)->GetName();
	TString _fileDir;


	TString __str=gSystem->GetFromPipe("echo $CMSSW_BASE");
	__str+="/src/";
	_fileDir = __str+TString(arr1->At(2)->GetName());

	TString _selection;
	if (arr1->GetEntries()>3) _selection=arr1->At(3)->GetName();


/*
cout<<_name<<endl;
cout<<_fileDir<<endl;
cout<<_methodname<<endl;
cout<<_selection<<endl;
*/

///check: is this trainer or computer real? Does weight folder exists?
	int isFile = (int) gSystem->Exec(TString("[ -f ")+_fileDir+TString(" ]"));

	bool _isOk = false;

///trainer doesn't have selection criteria
	if (isFile==0 && _selection.Length()==0 ) {
	 _mvaTrainersMap[_name] = new MVATrainer(_name,_name,std::string(_fileDir.Data()),std::string(""));

//cout<<"I'm trainer"<<endl;
//cout<<_name<<endl;

	} ///end of trainer ini

///mvacomputer requires to have _selection
	if (isFile==0 &&  _selection.Length()>0) {

	std::string _dir="";

	TString __dir = gSystem->GetFromPipe(TString("dirname ")+ _fileDir);
	_dir=std::string(__dir.Data());
	std::vector<std::string>  _mvaMethods;
	_mvaMethods.push_back(_methodname.Data());
	_mvaComputersMap[_name] = new MVAComputer(_name,_name,_dir,_mvaMethods);
	_mvaCutsMap[_name]= _selection;

	} ///end of computer ini

	if (isFile==0) _isOk = true;

///if settings are correct and weight folder or xml files exist, then save method in _mvaMethodsMap.
/// _mvaMethodsMap is a cenral storage of mva computers or trainers

	if (_isOk)  _mvaMethodsMap[_name]=_methodname;


} /// if (arr1->GetEntries()>2 ...


str.ReadLine(ff);


} ///end of while

return; 
} 

///  Do select events using MVAComputer or fill variables for MVATrainer

/// _isOk is true if event passed offline selection down.

/// _doexclude  to exclude some trainers from MVA selection!  if _doexclude is false,
///  then MVATrainer fill events selected by previous MVAComputer 

bool MVASelectAndFillEvent(bool _isOk,bool _doexclude,std::vector<math::XYZTLorentzVector> & _vecJets,std::vector<int> & _theFlav,std::vector<int> & _leading,double _weight)
{

bool _passMVA=true;
std::vector<double> * res=0;

TString _lastComp="";
for (std::map<TString,TString>::iterator it =  _mvaMethodsMap.begin(); it !=  _mvaMethodsMap.end(); it++)
{
bool _exclude=false;

//if mva computer
	if ( _mvaComputersMap.find(it->first)!= _mvaComputersMap.end() )
	{

	//      cout<<"I'm computer "<<_mvaComputersMap[it->first]->GetName()<<endl;
		 _mvaComputersMap[it->first]->ProcessEvent(_vecJets,_leading) ; /// calculate mva response
		  res= _mvaComputersMap[it->first]->GetOutPut(); /// get response
	 	double mvaVal =  res->at(0); /// in principle, one mva computer can carry several discriminators: res->at(1), res->at(2) etc
	//      cout<<"MVA response: "<<mvaVal<<endl;

        	if ( _mvaCutsMap.find(it->first)!= _mvaCutsMap.end() ) {
	        	TString _val = Form("%.5f",mvaVal);
		        TString _sel = _mvaCutsMap[it->first];
		        _sel.ReplaceAll(TString("CUT"),_val);
		        double xx=1.;
		        _passMVA=(bool) TFormula("fm",_sel.Data()).Eval(xx);  // if event pass mva, then...
		        _lastComp = _mvaComputersMap[it->first]->GetName();
		        }
	} ///if mva computer

	if (_doexclude)
	if (_exclude_trainer_from_selection.find(it->first)!=_exclude_trainer_from_selection.end()) _exclude=true;

	if ( _mvaTrainersMap.find(it->first)!= _mvaTrainersMap.end() && (_passMVA || _exclude) && _isOk)
	{
	//      cout<<"I'm trainer "<< _mvaTrainersMap[it->first]->GetName() <<endl;

	        if (res!=0) {
	//      cout<<"mva response from previous computer = "<<res->at(0)<<endl;
	        _mvaTrainersMap[it->first]->ProcessEvent(_vecJets,_theFlav,_leading,_weight,kFALSE,res,0);  /// fill int$
	        }
        	else  _mvaTrainersMap[it->first]->ProcessEvent(_vecJets,_theFlav,_leading,_weight);


	} ///if mva trainer

} /// end of iteration over all available computers and trainers



return _passMVA;
}

	///_MVA_BEGIN_ clean up memory and calculation of events went trought mva computers

void MVAClean()
{

for (std::map<TString,TString>::iterator it =  _mvaMethodsMap.begin(); it !=  _mvaMethodsMap.end(); it++)
{
	if ( _mvaComputersMap.find(it->first)!= _mvaComputersMap.end() ) delete  _mvaComputersMap[it->first];
	if ( _mvaTrainersMap.find(it->first)!= _mvaTrainersMap.end())  delete  _mvaTrainersMap[it->first];
}

return; 
}
#endif

