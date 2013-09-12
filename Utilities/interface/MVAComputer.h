#ifndef MVAComputer_H
#define MVAComputer_H

//#define _test_


#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TClass.h"
#include "TList.h"
#include "TROOT.h"


#include "Analysis/Utilities/interface/MathRoot.h"
#include "Analysis/Utilities/interface/MVAInputVars.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iostream>

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

///from utils.h --> text processing support
#include "Analysis/Utilities/interface/field.C.h"
#include "Analysis/Utilities/interface/findLine.C.h"


class MVAComputer : public TNamed {

public:

MVAComputer ():TNamed("MVAComputer","MVAComputer"),
_isOk(false)
{

output_=0;
reader_=0;

///_test_, remove the block  after testing

#ifdef _test_


#endif

////_test_
};



///jets is  an array of TLorentzVectors

MVAComputer (TString name, TString title, std::string  mvaPath, std::vector<std::string> & mvaMethod ):TNamed(name,title),
_isOk(true)
{
mvaPath_ = mvaPath;
mvaMethod_ = mvaMethod;

if  (mvaPath_.empty() || mvaMethod_.size() ==0 ) { std::cout<<"Please, provide me proper settings, return"<<std::endl; _isOk=false ;} 



if (!_isOk) return;


output_ = new std::vector<double>();
reader_  = new TMVA::Reader( "!Color:Silent" );

///Extract all variables used by MVA!
/**
        idea of parsing *.weights.xml
        1) read section

        <Variables >    </Variables>

Example :

 <Variables NVar="6">
    <Variable VarIndex="0" Expression="detajet2jet3" Label="detajet2jet3" Title="detajet2jet3" Unit="F" Internal="detajet2jet3" Type="F" Min="2.13906e-05$
    <Variable VarIndex="1" Expression="detajet3jet4" Label="detajet3jet4" Title="detajet3jet4" Unit="F" Internal="detajet3jet4" Type="F" Min="3.42131e-05$
    <Variable VarIndex="2" Expression="dphiMETlepton" Label="dphiMETlepton" Title="dphiMETlepton" Unit="F" Internal="dphiMETlepton" Type="F" Min="6.91414$
    <Variable VarIndex="3" Expression="lepeta" Label="lepeta" Title="lepeta" Unit="F" Internal="lepeta" Type="F" Min="3.16916e-05" Max="2.49505"/>
    <Variable VarIndex="4" Expression="maxdijetmass" Label="maxdijetmass" Title="maxdijetmass" Unit="F" Internal="maxdijetmass" Type="F" Min="0.00531079"$
    <Variable VarIndex="5" Expression="mindijetmass" Label="mindijetmass" Title="mindijetmass" Unit="F" Internal="mindijetmass" Type="F" Min="0.00186576"$
  </Variables>


2) extract from section: NVars, VarIndex && Expression

        It can be done with a help of utilites from utils.h --> findLine && field

        So, do it : #include "findLine.C.h";
                     #include "field.C.h"
                this *.h can be made by calling function createCPPcode 'findLine' or  createCPPcode 'field' from utils.h



**/
TString prefix;
TString   methodName;
TString weightfile;

   std::vector<std::string>::const_iterator it2 = mvaMethod_.begin();
   prefix = "TMVAClassification";
   methodName = TString(*it2)+TString(" method");
   weightfile = TString(mvaPath_) + TString("/") + prefix + TString("_") +  TString(*it2) + TString(".weights.xml");


TString res=findLine("Variable","NVar",weightfile,TString(mvaPath_) + "/");
res.ReplaceAll("<"," ");  res.ReplaceAll("="," ");
res.ReplaceAll(">"," "); res.ReplaceAll("\""," ");
res=field(TString("\"")+res+TString("\""),"3",TString(mvaPath_) + "/");
Int_t nVars=res.Atoi();

res=findLine("Variable","VarIndex",weightfile,TString(mvaPath_) + "/");
res.ReplaceAll("<"," ");  res.ReplaceAll("="," ");  res.ReplaceAll(">"," ");
res.ReplaceAll("\""," ");
res=field(TString("\"")+res+TString("\""),"5",TString(mvaPath_) + "/");
res.ReplaceAll("\n"," ");

for (Int_t i=0;i<nVars;i++)
{
	TString res2=res;
	TString varNum=Form("%d",i+1);
	res2=field(TString("\"")+res2+TString("\""),varNum,TString(mvaPath_) + "/");

//std::cout<<"my var is "<<res2<<std::endl;
vars_[res2.Data()]=new float(-1e10);


//std::cout<<"I'm adding variables"<<std::endl;
// reader_->AddVariable(res2 , &vars_[i] );
 reader_->AddVariable(res2 , vars_[res2.Data()] );
//std::cout<<"I'm adding variables 2"<<std::endl;

}


for (std::vector<std::string>::const_iterator it = mvaMethod_.begin();it!=mvaMethod_.end();++it)
{
        if (it->empty()) continue;

    prefix = "TMVAClassification";
     methodName = TString(*it)+" method";
     weightfile = TString(mvaPath_) + "/" + prefix + "_" +  TString(*it) + ".weights.xml";

 std::cout<<"reading file"<<weightfile<<std::endl;
 std::cout<<"reading method"<<methodName<<std::endl;
 reader_->BookMVA( methodName, weightfile );
}

};

///Not supported now, cal is internaly created object
//void ProcessEvent(MVACalculator * & cal_)
//void ProcessEvent(TObjArray * jets_)
///New Code to test

void ProcessEvent(std::vector<math::XYZTLorentzVector> & jets_, std::vector<int> & leading_)
{
output_->clear();

if (!_isOk ) {std::cout <<"Something wrong with initiallization, return"<<std::endl; return;}
//if (!cal_ ) {std::cout <<"Something wrong with MVACalculator, return"<<std::endl; return;}
if (jets_.size() == 0 ) {
///std::cout <<"Something wrong with jets, return"<<std::endl; 
return;
}

 std::vector<int> flav_;
MVAInputVars KinVars("KinVars","KinVars",jets_,flav_,leading_);

///to store all available current kinematic variables
std::map<std::string,float> vals_;



        vals_["sumEt"]=KinVars.sumEt();
        vals_["Et1"]=KinVars.Et1();
//        vals_["lepeta"]=KinVars.lepeta();
//        vals_["MET"]=KinVars.MET();
        vals_["dphiMETlepton"]=KinVars.dphiMETlepton();
        vals_["detajet2jet3"]=KinVars.detajet2jet3();
        vals_["detajet3jet4"]=KinVars.detajet3jet4();
        vals_["mindijetmass"]=KinVars.mindijetmass();
        vals_["maxdijetmass"]=KinVars.maxdijetmass();
        vals_["sphericity"]=KinVars.sphericity();
        vals_["aplanarity"]=KinVars.aplanarity();
        vals_["Et2"]=KinVars.Et2();
        vals_["Et3"]=KinVars.Et3();
        vals_["Eta1"]=KinVars.Eta1();
        vals_["Eta2"]=KinVars.Eta2();
        vals_["Eta3"]=KinVars.Eta3();
        vals_["dptjet1jet2"]=KinVars.dptjet1jet2();
        vals_["dptjet1jet3"]=KinVars.dptjet1jet3();
        vals_["djet1jet2pt"]=KinVars.djet1jet2pt();
        vals_["detajet1jet2"]=KinVars.detajet1jet2();
        vals_["dthetajet1jet2_boost"]=KinVars.dthetajet1jet2_boost();
        vals_["dphijet1jet2_boost"]=KinVars.dphijet1jet2_boost();
        vals_["dphijet1jet2"]=KinVars.dphijet1jet2();
        vals_["dphijet2jet3"]=KinVars.dphijet2jet3();
        vals_["dphijet1jet3"]=KinVars.dphijet1jet3();
        vals_["Et2byEt1"]=KinVars.Et2byEt1();
        vals_["Et3byEt2"]=KinVars.Et3byEt2();
        vals_["Et3byEt1"]=KinVars.Et3byEt1();
        vals_["sphericity_boost"]=KinVars.sphericity_boost();


///new vars: boosted and flav content


 vals_["thetajet1_boost12"]=KinVars.thetajet1_boost12();
 vals_["thetajet3_boost12"]=KinVars.thetajet3_boost12();
// vals_["thetajet2_boost12"]=KinVars.thetajet2_boost12();
// vals_["thetajet1_boost123"]=KinVars.thetajet1_boost123();
// vals_["thetajet2_boost123"]=KinVars.thetajet2_boost123();
// vals_["dphijet1jet2_boost12"]=KinVars.dphijet1jet2_boost12();
 vals_["dphijet2jet3_boost12"]=KinVars.dphijet2jet3_boost12();
// vals_["dphijet1jet3_boost12"]=KinVars.dphijet1jet3_boost12();



 vals_["flav_cont"] = KinVars.flav_cont();



///new event shape variables

 vals_["D"] = KinVars.D();
 vals_["isotropy"] = KinVars.isotropy();


///new variable on jet area

    vals_["JA1byJA2"]=KinVars.JA1byJA2();
    vals_["JA2byJA3"]=KinVars.JA2byJA3();
    vals_["JA1byJA3"]=KinVars.JA1byJA3();

    vals_["ptd1byptd2"]=KinVars.ptd1byptd2();
    vals_["ptd2byptd3"]=KinVars.ptd2byptd3();
    vals_["ptd1byptd3"]=KinVars.ptd1byptd3();


    vals_["JA1"]=KinVars.JA1();
    vals_["JA2"]=KinVars.JA2();
    vals_["JA3"]=KinVars.JA3();

    vals_["ptd1"]=KinVars.ptd1();
    vals_["ptd2"]=KinVars.ptd2();
    vals_["ptd3"]=KinVars.ptd3();

/// boosted variables
    vals_["phijet3_boost12"]=KinVars.phijet3_boost12();
    vals_["phijet3_boost13"]=KinVars.phijet3_boost13();
    vals_["thetajet1_boost13"]=KinVars.thetajet1_boost13();
    vals_["thetajet2_boost13"]=KinVars.thetajet2_boost13();
    vals_["dphijet2jet3_boost13"]=KinVars.dphijet2jet3_boost13();



/// new event shape and jet broadening variables
    vals_["thrust"]=KinVars.thrust();
    vals_["major"]=KinVars.major(); 
    vals_["minor"]=KinVars.minor();
    vals_["oblateness"]=KinVars.oblateness();
    vals_["Bmax"]=KinVars.Bmax();
    vals_["Bmin"]=KinVars.Bmin();
    vals_["Bdiff"]=KinVars.Bdiff();



std::map<std::string,float *>::iterator it;
for ( it=vars_.begin() ; it != vars_.end(); it++ )
{
std::string _variableName = (*it).first;
float * _variableVal = (*it).second;

* _variableVal = vals_[_variableName];

//std::cout<<_variableName<<" = "<< * _variableVal<<std::endl;


}




for (std::vector<std::string>::const_iterator it = mvaMethod_.begin();it!=mvaMethod_.end();++it)
{
        if (it->empty()) continue;

TString    prefix = "TMVAClassification";
TString      methodName = TString(*it)+" method";
//TString      weightfile = TString(mvaPath_) + "/" + prefix + "_" +  TString(*it) + ".weights.xml";

output_->push_back(reader_->EvaluateMVA(methodName));


}

};

std::vector<double> * GetOutPut() { return output_;}

~MVAComputer() {

if (reader_) delete reader_; if (output_) delete output_; 

std::map<std::string,float *>::iterator it;
for ( it=vars_.begin() ; it != vars_.end(); it++ ) if ((*it).second) delete (*it).second;


 }




private:

std::string  mvaPath_;
std::vector<std::string>  mvaMethod_;

bool _isOk;
TMVA::Reader * reader_;
std::vector<double> *  output_; 

std::map<std::string,float *> vars_;		

public:
        ClassDef(MVAComputer,2)
};











#ifdef __MAKECINT__
ClassImp(MVAComputer)
#endif


#endif
