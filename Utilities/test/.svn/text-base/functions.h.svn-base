#ifndef FUNCTIONS_SAMPLE_H
#define FUNCTIONS_SAMPLE_H

#include "TPRegexp.h"
#include "TString.h"
#include "TObjString.h"
#include <vector>

#include <fstream>
#include <string>
#include <iostream>

using namespace std;

TString GetPartTString(TString str,TPRegexp reg,Ssiz_t part)
{

if (part ==0 ) return TString("error");
TString str2;
Ssiz_t ind=0;


while (part>0) {
ind=str.Index(reg,ind);
str2 = str(reg,ind);

ind++;
part--;
}
return str2;

}

Int_t GetTStringSize(TString str,TPRegexp reg)
{

Int_t part=0;
TString str2;
Ssiz_t ind=1;



while (ind>0) {
ind=str.Index(reg,ind);
str2 = str(reg,ind);
ind++;
part++;
}
return part-1;

}

///Read file containig processed by TSelector 
void Readinput( char *nameFile, char ** listName, int &Number) {
ifstream  infile;
infile.open(nameFile,std::ios::in);
if (!infile) {std::cout << "Couldn't read the file"; return;}
Number=0;
 
//while (!infile.eof()) {
while (infile) {
        char *buf = new char[256];
        infile.getline(buf,256);
        listName[Number]=buf;

        Number++;
}
infile.close();
return;
}


///Read file containig processed by TSelector 
void Readinput_v1( char *nameFile, std::vector<std::string> & listFiles, int &Number) {
ifstream  infile;
infile.open(nameFile,std::ios::in);
if (!infile) {std::cout << "Couldn't read the file"; return;}
Number=0;

std::string str;
int num;
double dbl;

while (1) {

	infile>>str;

cout<<str<<endl;	

if (!infile.eof()) {
        listFiles.push_back(str);
} else  break;


if (infile.eof()) break;
Number++; 
}
infile.close();
return;
}



///Read file containig processed by TSelector 
void Readinput_v2( char *nameFile, std::vector<std::string> & listFiles, std::vector<std::string> & process,std::vector<std::string> & type, std::vector<double> & xsects, std::vector<double> & effs ,std::vector<double> & wgts ,std::vector<int> & mergeList,int &Number) {
ifstream  infile;
infile.open(nameFile,std::ios::in);
if (!infile) {std::cout << "Couldn't read the file"; return;}
Number=0;

std::string str;
int num;
double dbl;

while (1) {
//while (infile) {
//        char *buf = new char[256];
//        infile.getline(buf,256);
//        listName[Number]=buf;


	infile>>str;


//cout<<str<<endl;

if (!infile.eof()) {
	listFiles.push_back(str);
} else  break;

	infile>>str;

//cout<<str<<endl;

if (!infile.eof()) {
	process.push_back(str);
} else  break;

	infile>>str;

if (!infile.eof()) {
	type.push_back(str);
} else  break;

	infile>>dbl;

if (!infile.eof()) {
	xsects.push_back(dbl);
} else  break;

	infile>>dbl;


if (!infile.eof()) {
	effs.push_back(dbl);
} else  break;

///It's optional
	infile>>dbl;



if (!infile.eof()) {
	wgts.push_back(dbl);
} else  break;


	infile>>num;
if (!infile.eof()) {
	mergeList.push_back(num);
} else  break;



if (infile.eof()) break;
Number++; 
}
infile.close();
return;
}


///Read file containig processed by TSelector
/// 
void Readinput_v3( char *nameFile, std::vector<std::string> & listFiles,  std::vector<std::string> & configs,std::vector<int> & mergeList,std::vector<double> & wgtList, std::vector<std::string> & nameList,std::vector<std::string> & typeList,int &Number) 
{
ifstream  infile;
infile.open(nameFile,std::ios::in);
if (!infile) {std::cout << "Couldn't read the file"; return;}
Number=0;

std::string str;
int num;
double dbl;

while (1) {
//while (infile) {
//        char *buf = new char[256];
//        infile.getline(buf,256);
//        listName[Number]=buf;


	infile>>str;

if (!infile.eof()) {
	listFiles.push_back(str);
} else  break;

	infile>>str;

if (!infile.eof()) {
	nameList.push_back(str);
} else  break;


	infile>>str;

if (!infile.eof()) {
	typeList.push_back(str);
} else  break;


	infile>>str;

if (!infile.eof()) {
	configs.push_back(str);
} else  break;


        infile>>dbl;

if (!infile.eof()) {
        wgtList.push_back(dbl);
} else  break;


	infile>>num;
if (!infile.eof()) {
	mergeList.push_back(num);
} else  break;



///Count the number of valid lines
if (infile.eof()) break;
Number++; 
}
infile.close();
return;
}

///Remove step0 from name of sample
TString RemoveStep0(TString samplename)
{

	TString tmpstr =  GetPartTString(samplename,TPRegexp("\\w+_\\d+\\.\\d+:\\d+\\.\\d+_"),1);
	TString tmpstr2 =  GetPartTString(tmpstr,TPRegexp("\\w+_"),1);
	tmpstr2.Chop();
	TString  tmpstr3 = samplename;
	tmpstr3.Remove(0,tmpstr.Length()-1);
	TString  tmpstr4=tmpstr2+tmpstr3; ///First part of each histo belonging to the concrete sampleName


///return either Hbb_10.0000:1000000.0000_30.00000:100000.00000
///	or Hbb 
return tmpstr4;
}

///Remove all 'vars' from the name of histo:  Hbb_10.0000:1000000.0000_30.00000:100000.00000__cp_eta1_step0 ---> Hbb_cp_eta1_step0
///type 0 --   return Hbb_cp_eta1_step0
///    1  -- return _cp_eta1_step0
/// 2 -- return Hbb
/// 3 -- return _10.0000:1000000.0000_30.00000:100000.00000_
TString RemoveAllVars(TString samplename, Int_t type=0)
{
TString tmpstr =  GetPartTString(samplename,TPRegexp("\\w+_\\d+\\.\\d+:\\d+\\.\\d+_"),1);
TString tmpstr2 =  GetPartTString(tmpstr,TPRegexp("\\w+_"),1);
tmpstr2.Chop(); ///Contains 'Hbb'



TString varCombination;
TString AllvarCombination;
UInt_t numVars = GetTStringSize(samplename,TPRegexp("_\\d+\\.\\d+:\\d+\\.\\d+_"));
for (UInt_t i=1;i<=numVars;i++) {
	varCombination = GetPartTString(samplename,TPRegexp("_\\d+\\.\\d+:\\d+\\.\\d+_"),i);
	if (i<numVars) varCombination.Chop();
	AllvarCombination+=varCombination;
}
TString tmpstr3=tmpstr2+AllvarCombination;



TString tmpstr4=samplename;

tmpstr4.Remove(0,tmpstr3.Length());


if (type==1) return tmpstr4;

if (type==2) return tmpstr2;

if (type==3) return AllvarCombination;

tmpstr4=tmpstr2+tmpstr4;



return tmpstr4;
}

///Get Name of sample
TString GetTypeName(const char * samplename, vector<Double_t> & mins, vector<Double_t> & maxes) {

TString stroutput(samplename);
stroutput+="_";

                for (unsigned int i=0;i<mins.size();i++) {      


                        TString stradd;
                        stradd.Form("%.3f:%.3f_",mins[i],maxes[i]);
                        stroutput+=stradd;


                }

//std::cout<<"TypeName  "<<stroutput<<std::endl;
return stroutput;
}

Int_t GetIDNumber(vector<int> & vec)
{
int mult=10;
Int_t res=0;
for (unsigned int i=0;i<vec.size();i++) {
	
	res+=mult*vec[i];
	mult*=mult;
}
return res;
}

#endif


/**

Useful application

root [9]  GetPartTString(GetPartTString(str,TPRegexp("\\w+"),1),TPRegexp("_\\w+_"),1)
(class TString)"_header_SGN_SemiMuonic_"
root [10]  GetPartTString(GetPartTString(TString("_header_SGN_SemiMuonic_"),TPRegexp("\\w+"),1),TPRegexp("_\\w+_"),1)
(class TString)"_header_SGN_SemiMuonic_"


**/
