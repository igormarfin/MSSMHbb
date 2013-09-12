#ifndef MVATrainer_H
#define MVATrainer_H


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

#include "TFile.h"
#include "TTree.h"



#include "Analysis/Utilities/interface/MathRoot.h"
#include "Analysis/Utilities/interface/MVAInputVars.h"

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <iostream>


///from utils.h --> text processing support
#include "Analysis/Utilities/interface/field.C.h"
#include "Analysis/Utilities/interface/findLine.C.h"


class MVATrainer : public TNamed {

public:

MVATrainer ():TNamed("MVATrainer","MVATrainer")
{

tree_=0;
file_=0;

_isOk=false;

///_test_, remove the block  after testing

#ifdef _test_


#endif

////_test_
};



///jets is  an array of TLorentzVectors

MVATrainer (TString name, TString title, std::string  _treexml, std::string _file="ntuple.root" ) : TNamed(name,title),
tree_(0),_isOk(true),file_(0),nCounter(0),nCounterWeight(0)
{

if  (_treexml.empty() ) { std::cout<<"Please, provide me proper settings, return"<<std::endl; _isOk=false ;} 

treexml =TString( _treexml);

///External file must be opened
//file_= TFile::Open(_file.c_str(),"RECREATE");




if (!_isOk) return;


///Extract all variables used by MVA!
/**
        idea of parsing *.weights.xml
        1) read section

        <Variables >    </Variables>

Example :

<Variables NVar="25" TreeName="KinVars">
        <Variable VarIndex="0" Expression="detajet2jet3"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="1" Expression="detajet3jet4" Type="F" BINS="100"  MIN="0" MAX="3.5"/> 
        <Variable VarIndex="2" Expression="maxdijetmass" Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="3" Expression="mindijetmass"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="4" Expression="aplanarity"  Type="F" BINS="100"  MIN="0" MAX="0.5"/>
        <Variable VarIndex="5" Expression="sphericity"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="6" Expression="Et1"  Type="F" BINS="100"  MIN="0" MAX="1000"/>
        <Variable VarIndex="7" Expression="Et2"  Type="F" BINS="100"  MIN="0" MAX="1000"/>
        <Variable VarIndex="8" Expression="Et3"  Type="F" BINS="100"  MIN="0" MAX="1000"/>
        <Variable VarIndex="9" Expression="Eta1"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="10" Expression="Eta2"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="11" Expression="Eta3"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="12" Expression="dptjet1jet2"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="13" Expression="dptjet1jet3"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="14" Expression="djet1jet2pt"  Type="F" BINS="100"  MIN="0" MAX="1"/>


2) extract from section: NVars, VarIndex && Expression

        It can be done with a help of utilites from utils.h --> findLine && field

        So, do it : #include "findLine.C.h";
                     #include "field.C.h"
                this *.h can be made by calling function createCPPcode 'findLine' or  createCPPcode 'field' from utils.h



**/

TString treexmlPath_="";
TObjArray* arr = treexml.Tokenize("/");
for(int i=0;i<arr->GetEntries()-1;i++) treexmlPath_+=TString("/")+arr->At(i)->GetName();



TString res=findLine("Variable","NVar",treexml,TString(treexmlPath_) + "/");
res.ReplaceAll("<"," ");  res.ReplaceAll("="," ");
res.ReplaceAll(">"," "); res.ReplaceAll("\""," ");
res=field(TString("\"")+res+TString("\""),"3",TString(treexmlPath_) + "/");
Int_t nVars=res.Atoi();


res=findLine("Variable","NVar",treexml,TString(treexmlPath_) + "/");
res.ReplaceAll("<"," ");  res.ReplaceAll("="," ");
res.ReplaceAll(">"," "); res.ReplaceAll("\""," ");
res=field(TString("\"")+res+TString("\""),"5",TString(treexmlPath_) + "/");
TString nameTree = res;

if (tree_==0) {

//std::cout<<"Tree with name "<<nameTree<<std::endl;

tree_ = new TTree(nameTree,nameTree);
//tree_->SetDirectory(gROOT);
}


res=findLine("Variable","VarIndex",treexml,TString(treexmlPath_) + "/");
res.ReplaceAll("<"," ");  res.ReplaceAll("="," ");  res.ReplaceAll(">"," ");
res.ReplaceAll("\""," ");
res=field(TString("\"")+res+TString("\""),"5",TString(treexmlPath_) + "/");
res.ReplaceAll("\n"," ");

TString type=findLine("Variable","VarIndex",treexml,TString(treexmlPath_) + "/");
type.ReplaceAll("<"," ");  type.ReplaceAll("="," ");  type.ReplaceAll(">"," ");
type.ReplaceAll("\""," ");
type.ReplaceAll("/"," ");

type=field(TString("\"")+type+TString("\""),"7",TString(treexmlPath_) + "/");
type.ReplaceAll("\n"," ");

//cout<<"res = "<<res<<std::endl;
//cout<<"2 type = "<<type<<std::endl;


for (int i=0;i<nVars;i++)
{
TString res2=res;
TString type2=type;
TString varNum=Form("%d",i+1);
res2=field(TString("\"")+res2+TString("\""),varNum,TString(treexmlPath_) + "/");
type2=field(TString("\"")+type2+TString("\""),varNum,TString(treexmlPath_) + "/");

//std::cout<<"my var is "<<res2<<std::endl;
//std::cout<<"my type  is "<<type2<<std::endl;


//std::cout<<"I'm adding variable "<<res2+"/"+type2<<std::endl;

vars_[res2.Data()]=new float(-1e10);

 tree_->Branch(res2,vars_[res2.Data()],res2+"/"+type2);


//std::cout<<"I'm adding variables 2"<<std::endl;

_usedvars[i]=res2.Data();


}

vars_["weight"]=new float(-1e10);

tree_->Branch("weight",vars_["weight"],"weight/F");


};

///Not supported now, cal is internaly created object
//void ProcessEvent(MVACalculator * & cal_)
//void ProcessEvent(TObjArray * jets_)
///New Code to test

void ProcessEvent(std::vector<math::XYZTLorentzVector> & jets_, std::vector<int> & flav_, std::vector<int> & leadingjets_, double weight_=1e0, Bool_t _print = kFALSE, std::vector<double> *  mvaOutPut=(std::vector<double> *)0, int iNdexMva=0) 
{

if (!_isOk ) {std::cout <<"Something wrong with initiallization, return"<<std::endl; return;}
if (jets_.size() == 0 ) {
//std::cout <<"Something wrong with jets, return"<<std::endl; 
return;
}


nCounter++;
nCounterWeight+=weight_;

//if (_print) {
//
//cout<<"weight = "<<weight_<<endl;
//cout<<"weight events = "<<nCounterWeight<<endl;
//}

 MVAInputVars KinVars("KinVars","KinVars",jets_, flav_,leadingjets_);

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

	vals_["Pt1"]=KinVars.Pt1();
	vals_["Pt2"]=KinVars.Pt2();
	vals_["Pt3"]=KinVars.Pt3();

	vals_["Pt1_b"]=KinVars.Pt1_b();
	vals_["Pt2_b"]=KinVars.Pt2_b();
	vals_["Pt3_b"]=KinVars.Pt3_b();

//MC
	vals_["Pt1_nonb"]=KinVars.Pt1_nonb();
	vals_["Pt2_nonb"]=KinVars.Pt2_nonb();
	vals_["Pt3_nonb"]=KinVars.Pt3_nonb();


	vals_["Eta1_b"]=KinVars.Eta1_b();
	vals_["Eta2_b"]=KinVars.Eta2_b();
	vals_["Eta3_b"]=KinVars.Eta3_b();

	vals_["Eta1_nonb"]=KinVars.Eta1_nonb();
	vals_["Eta2_nonb"]=KinVars.Eta2_nonb();
	vals_["Eta3_nonb"]=KinVars.Eta3_nonb();
	vals_["M12"]=KinVars.M12();


	vals_["Pt1_c"]=KinVars.Pt1_c();
	vals_["Pt2_c"]=KinVars.Pt2_c();
	vals_["Pt1_q"]=KinVars.Pt1_q();
	vals_["Pt2_q"]=KinVars.Pt2_q();

	vals_["Eta1_c"]=KinVars.Eta1_c();
	vals_["Eta2_c"]=KinVars.Eta2_c();
	vals_["Eta1_q"]=KinVars.Eta1_q();
	vals_["Eta2_q"]=KinVars.Eta2_q();

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


///new shape vars

 vals_["D"] = KinVars.D();
 vals_["isotropy"] = KinVars.isotropy();


//PU variables

	vals_["Loose_id"]=KinVars.Loose_id();
	vals_["Medium_id"]=KinVars.Medium_id();
	vals_["Tight_id"]=KinVars.Tight_id();
	vals_["nPV"]=KinVars.nPV();


/// Jet Area vars

    vals_["JA1byJA2"]=KinVars.JA1byJA2();
    vals_["JA2byJA3"]=KinVars.JA2byJA3();
    vals_["JA1byJA3"]=KinVars.JA1byJA3();

    vals_["ptd1byptd2"]=KinVars.ptd1byptd2();     
    vals_["ptd2byptd3"]=KinVars.ptd2byptd3();     
    vals_["ptd1byptd3"]=KinVars.ptd1byptd3();     

    vals_["JA1"]=KinVars.JA1();
    vals_["JA2"]=KinVars.JA2();
    vals_["JA3"]=KinVars.JA1();

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



	vals_["weight"]=weight_;

if( mvaOutPut)  {
	vals_["mva"]=mvaOutPut->at(iNdexMva);
//	std::cout<<"MVA output "<<vals_["mva"]<<std::endl;
}
else vals_["mva"]=-1e10;



//if (nCounter%10000==0) 
//if ( _print)  cout<<"================"<<endl;


std::map<std::string,float *>::iterator it;
std::map<std::string,float >::iterator it2;

float * _variableVal;
for ( it=vars_.begin() ; it != vars_.end(); it++ )
{
std::string _variableName = (*it).first;
float * _variableVal = (*it).second;
it2 =  vals_.find(_variableName);
if ( it2!= vals_.end() )  * _variableVal = vals_[_variableName];

//if (nCounter%10000==0)  
//if ( _print) std::cout<<GetName()<<" : "<<_variableName<<" = "<< * _variableVal<<std::endl;

}
//if (nCounter%10000==0) 
//if ( _print)  cout<<"================"<<endl;

if (tree_!=0) 
{
tree_->Fill();

}

};


~MVATrainer() {

//if (reader_) delete reader_; 
/*

if (file_) {

     const char* filename = file_->GetName();
        TString str(filename);
	str+=TString(":/");

file_->cd(str);
}


if (tree_) tree_->Write();
if (tree_) delete tree_;
if (file_) file_->Close();
if (file_) delete file_;

*/


std::map<std::string,float *>::iterator it;
for ( it=vars_.begin() ; it != vars_.end(); it++ ) if ((*it).second) delete (*it).second;


std::cout<<GetName()<<" : "<<"processed # events --> "<<nCounter<<std::endl;
std::cout<<GetName()<<" : "<<"processed # weighted events --> "<<nCounterWeight<<std::endl;

 };



TTree * tree_;

private:
bool _isOk;
TFile * file_;
std::map<int,std::string> _usedvars;
TString  treexml;
Int_t nCounter;
Double_t nCounterWeight;
std::map<std::string,float *> vars_;		

public:
        ClassDef(MVATrainer,2)
};











#ifdef __MAKECINT__
ClassImp(MVATrainer)
#endif


#endif
