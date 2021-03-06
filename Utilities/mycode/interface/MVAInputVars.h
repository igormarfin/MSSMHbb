#ifndef MVAInputVars_h
#define MVAInputVars_h



#include "EventShapeVariables.h"
#include "EventShape.h"

#include "TObject.h"
#include "TNamed.h"

#include "TString.h"
#include "TList.h"
#include <vector>
#include "TMath.h"
#include "TLorentzVector.h" // to perform boost
#include <iostream>
#include "TObjArray.h"

#include "Analysis/Utilities/interface/MathRoot.h"

///to get access to PU variables
///extern  int nPV;

//#include "Analysis/Utilities/interface/HbbNtuple.h"

class MVAInputVars :public TNamed 
{

  // common calculator class for likelihood
  // variables in semi leptonic ttbar decays

public:

  MVAInputVars():TNamed("MVAInputVars",(const char *)("")) { }
  MVAInputVars(TString name, TString header,const std::vector<math::XYZTLorentzVector>&, std::vector<int> & , std::vector<int> &);

  ~MVAInputVars();
  
  double sumEt() const { return var_sumEt; }
  double Et1() const { if (var_sumEt>0  &&  var_Et1>0) return var_Et1/var_sumEt; return -1e10; }
  double lepeta() const { return fabs(var_lepeta); }
  double MET() const { return var_MET; }
    
double dphiMETlepton() const { return var_dphiMETlepton; }
  
double detajet2jet3() const { return var_detajet2jet3; }
double detajet3jet4() const { return var_detajet3jet4; }

double mindijetmass() const { return var_mindijetmass/massalljets; }
double maxdijetmass() const { return var_maxdijetmass/massalljets; }

double mindRjetlepton() const { return var_mindRjetlepton; }

   
double DeltaPhi(math::XYZTLorentzVector, math::XYZTLorentzVector);
double DeltaR(math::XYZTLorentzVector, math::XYZTLorentzVector);

double sphericity() const { return var_sphericity;}
double aplanarity() const { return var_aplanarity;}

double Et2 () const { if (var_sumEt>0) return var_Et2/var_sumEt; return -1e10;}
double Et3 () const { if (var_sumEt>0) return var_Et3/var_sumEt; return -1e10;}
double Eta1() const {return var_Eta1;}
double Eta2() const {return var_Eta2;}
double Eta3() const {return var_Eta3;}
double dptjet1jet2 () const  {if (var_sumEt>0) return var_dptjet1jet2/var_sumEt; return -1e10;} 
double dptjet1jet3 () const  {if (var_sumEt>0) return var_dptjet1jet3/var_sumEt; return -1e10;} 
double djet1jet2pt () const  {if (var_sumEt>0) return var_djet1jet2pt/var_sumEt; return -1e10;}


double detajet1jet2() const { return var_detajet1jet2; }
double dthetajet1jet2_boost() const { return var_dthetajet1jet2_boost;}
double thetajet1_boost12() const { return var_thetajet1_boost12;}
double thetajet3_boost12() const { return var_thetajet3_boost12;}
//double thetajet2_boost12() const { return var_thetajet2_boost12;}
//double thetajet1_boost123() const { return var_thetajet1_boost123;}
//double thetajet2_boost123() const { return var_thetajet2_boost123;}
//double dphijet1jet2_boost12() const { return var_dphijet1jet2_boost12;}
double dphijet2jet3_boost12() const { return var_dphijet2jet3_boost12;}
//double dphijet1jet3_boost12() const { return var_dphijet1jet3_boost12;}

double dphijet1jet2_boost() const { return var_dphijet1jet2_boost;}

double Et2byEt1 () const {if (var_Et1>0 && var_Et2>0 ) return var_Et2/var_Et1; return -1e10;}
double Et3byEt1 () const {if (var_Et1>0 && var_Et3>0) return var_Et3/var_Et1; return -1e10;}
double Et3byEt2 () const {if (var_Et2>0 && var_Et3>0)  return var_Et3/var_Et2; return -1e10;}

double dphijet1jet2 () const { return var_dphijet1jet2;}
double dphijet1jet3 () const { return var_dphijet1jet3;}
double dphijet2jet3 () const { return var_dphijet2jet3;}

double sphericity_boost() const {return var_sphericity_boost;}

double Pt1() const { return var_Et1;}
double Pt2() const { return var_Et2;}
double Pt3() const { return var_Et3;}


double Pt1_b() const { return var_Pt1_b;}
double Pt2_b() const { return var_Pt2_b;}
double Pt3_b() const { return var_Pt3_b;}


double Pt1_nonb() const { return var_Pt1_nonb;}
double Pt2_nonb() const { return var_Pt2_nonb;}
double Pt3_nonb() const { return var_Pt3_nonb;}


double Eta1_b() const { return var_Eta1_b;}
double Eta2_b() const { return var_Eta2_b;}
double Eta3_b() const { return var_Eta3_b;}


double Eta1_nonb() const { return var_Eta1_nonb;}
double Eta2_nonb() const { return var_Eta2_nonb;}
double Eta3_nonb() const { return var_Eta3_nonb;}

double M12() const { return var_M12;}

double Pt1_c() const { return var_Pt1_c;}
double Pt2_c() const { return var_Pt2_c;}
double Pt1_q() const { return var_Pt1_q;}
double Pt2_q() const { return var_Pt2_q;}

double Eta1_c() const { return var_Eta1_c;}
double Eta2_c() const { return var_Eta2_c;}
double Eta1_q() const { return var_Eta1_q;}
double Eta2_q() const { return var_Eta2_q;}

double flav_cont () const {return var_flav_cont;}


///new event shapes variables

///and measures the 4-jet structure of the event (D vanishes for a planar event)
double D()  const { return var_D;}

/// the return value is 1 for spherical events and 0 for events linear in r-phi. This function 
  /// needs the number of steps to determine how fine the granularity of the algorithm in phi 
  /// should be

double isotropy() const { return var_isotropy;}

///PU variables
double Loose_id() const {return var_Loose_id;}
double Medium_id() const {return var_Medium_id;}
double Tight_id() const {return var_Tight_id;}
double nPV() const { return var_nPV;}

double JA1byJA2 () const {return var_JA1byJA2;}
double JA2byJA3 () const {return var_JA2byJA3;}
double JA1byJA3 () const {return var_JA1byJA3;}

double ptd1byptd2 () const {return var_ptd1byptd2;}
double ptd2byptd3 () const {return var_ptd2byptd3;}
double ptd1byptd3 () const {return var_ptd1byptd3;}


double JA1 () const {return var_JA1;}
double JA2 () const {return var_JA2;}
double JA3 () const {return var_JA3;}

double ptd1 () const {return var_ptd1;}
double ptd2 () const {return var_ptd2;}
double ptd3 () const {return var_ptd3;}


double phijet3_boost12() const { return var_phijet3_boost12;}
double phijet3_boost13() const {return var_phijet3_boost13;}
double thetajet1_boost13() const {return var_thetajet1_boost13;}
double thetajet2_boost13() const {return var_thetajet2_boost13 ;}
double dphijet2jet3_boost13() const {return var_dphijet2jet3_boost13;}

double thrust()   const {return var_thrust;}
double major() const {return var_major;}
double minor() const {return var_minor;}
double oblateness() const {return var_oblateness;}
double Bmax() const {return var_Bmax;}
double Bmin() const {return var_Bmin;}
double Bdiff() const {return var_Bdiff;}


private:

double var_isotropy;
double var_D;


double var_thetajet1_boost12;
double var_thetajet3_boost12;
//double var_thetajet2_boost12;
//double var_thetajet1_boost123;
//double var_thetajet2_boost123;
//double var_dphijet1jet2_boost12;
double var_dphijet2jet3_boost12;
//double var_dphijet1jet3_boost12;



double var_flav_cont; ///flavor content of first 3 jets

double var_Pt1_c;
double var_Pt2_c;
double var_Pt1_q;
double var_Pt2_q;


double var_Eta1_c;
double var_Eta2_c;
double var_Eta1_q;
double var_Eta2_q;

double var_Pt1_b;
double var_Pt2_b;
double var_Pt3_b;

double var_Pt1_nonb;
double var_Pt2_nonb;
double var_Pt3_nonb;

double var_Eta1_b;
double var_Eta2_b;
double var_Eta3_b;

double var_Eta1_nonb;
double var_Eta2_nonb;
double var_Eta3_nonb;

double var_M12;



double var_dphijet1jet2;
double var_dphijet2jet3;
double var_dphijet1jet3;

double var_sphericity_boost;

double var_Et2;
double var_Et3;
double var_Eta1;
double var_Eta2;
double var_Eta3;
double var_dptjet1jet2;
double var_dptjet1jet3;
double var_djet1jet2pt;
double var_detajet1jet2;
double var_dthetajet1jet2_boost;
double var_dphijet1jet2_boost;

double var_sphericity;
double var_aplanarity;

  
  double var_sumEt;
  double var_Et1;
  double var_lepeta;
  double var_MET;
  
  double var_dphiMETlepton; 
  
  double var_detajet2jet3;
  double var_detajet3jet4;

  double var_mindijetmass;
  double var_maxdijetmass;

  double var_mindRjetlepton;
  

  double massalljets;

///PU variables
  double var_Loose_id;
  double var_Medium_id;
  double var_Tight_id;
  double var_nPV;


//Jet area related variables
double var_JA1byJA2;
double var_JA2byJA3;
double var_JA1byJA3;

double var_ptd1byptd2;
double var_ptd2byptd3;
double var_ptd1byptd3;


double var_JA1;
double var_JA2;
double var_JA3;

double var_ptd1;
double var_ptd2;
double var_ptd3;

///new boosted variables

double var_phijet3_boost12;
double var_phijet3_boost13;
double var_thetajet1_boost13;
double var_thetajet2_boost13;
double var_dphijet2jet3_boost13;

/// new event shape and jet broadening variables

double var_thrust;
double var_major;
double var_minor;
double var_oblateness;
double var_Bmax;
double var_Bmin;
double var_Bdiff;

public:
 	ClassDef(MVAInputVars,2);
};


 #ifdef __MAKECINT__
ClassImp(MVAInputVars)
#endif


#if defined(MVAInputVars_CC)


MVAInputVars::MVAInputVars(TString name, TString header,const std::vector<math::XYZTLorentzVector>& topJets, std::vector<int> & flav,std::vector<int>  & leadingJets):TNamed(name,header)
{ 



  unsigned int nJetsMax = topJets.size();


var_isotropy=-1e10;
var_D=-1e10;
 
//var_MET = -1e10;
var_sumEt = 0;

var_Pt1_b=-1e10;
var_Pt2_b=-1e10;
var_Pt3_b=-1e10;

var_Pt1_nonb=-1e10;
var_Pt2_nonb=-1e10;
var_Pt3_nonb=-1e10;

var_Eta1_b=-1e10;
var_Eta2_b=-1e10;
var_Eta3_b=-1e10;

var_Eta1_nonb=-1e10;
var_Eta2_nonb=-1e10;
var_Eta3_nonb=-1e10;

var_M12=-1e10;

 var_Pt1_c=-1e10;
 var_Pt2_c=-1e10;
 var_Pt1_q=-1e10;
 var_Pt2_q=-1e10;


 var_Eta1_c=-1e10;
 var_Eta2_c=-1e10;
 var_Eta1_q=-1e10;
 var_Eta2_q=-1e10;

 var_thetajet1_boost12=-1e10;
 var_thetajet3_boost12=-1e10;

// var_thetajet2_boost12=-1e10;
// var_thetajet1_boost123=-1e10;
// var_thetajet2_boost123=-1e10;

// var_dphijet1jet2_boost12=-1e10;
 var_dphijet2jet3_boost12=-1e10;
// var_dphijet1jet3_boost12=-1e10;



var_flav_cont=-1e10;

//PU variables
var_Loose_id=-1e10;
var_Medium_id=-1e10;
var_Tight_id=-1e10;
var_nPV=-1e10;


if (flav.size()>2)
var_flav_cont = abs(flav[0]-2)*3*3 +  abs(flav[1]-2)*3 + abs(flav[2]-2);
///  fabs(flav[i]-2) to change notation
/// b was having flav =2 , now it has 0, c was 1, now 1, q was 0 now it's 2


///PU variables

/*

  bool jetIDLoose[100];
  bool jetIDMedium[100];
  bool jetIDTight[100];
*/
  extern int nPV;

  var_nPV=nPV;



//Jet area related variables
var_JA1byJA2=-1e10;
var_JA2byJA3=-1e10;
var_JA1byJA3=-1e10;

var_ptd1byptd2=-1e10;
var_ptd2byptd3=-1e10;
var_ptd1byptd3=-1e10;

var_JA1=-1e10;
var_JA2=-1e10;
var_JA3=-1e10;

var_ptd1=-1e10;
var_ptd2=-1e10;
var_ptd3=-1e10;



extern    float areaJet[100];
extern    float ptdJet[100];//see CMS AN-2012/074

if (leadingJets.size()>0) 
{
var_JA1=areaJet[leadingJets[0]];
var_ptd1=ptdJet[leadingJets[0]];
}

if (leadingJets.size()>1) {
if (abs(areaJet[leadingJets[1]])>0e0) 	 var_JA1byJA2 = areaJet[leadingJets[0]]/areaJet[leadingJets[1]];
if (abs(ptdJet[leadingJets[1]])>0e0)    var_ptd1byptd2=ptdJet[leadingJets[0]]/ptdJet[leadingJets[1]];

var_JA2=areaJet[leadingJets[1]];
var_ptd2=ptdJet[leadingJets[1]];


//std::cout<<"JetArae[0]= "<< areaJet[leadingJets[0]]<<std::endl;
}

if (leadingJets.size()>2) {
if (abs(areaJet[leadingJets[2]])>0e0) {
	 var_JA1byJA3 = areaJet[leadingJets[0]]/areaJet[leadingJets[2]];
	 var_JA2byJA3 = areaJet[leadingJets[1]]/areaJet[leadingJets[2]];
}
if (abs(ptdJet[leadingJets[2]])>0e0)   {
 var_ptd1byptd3=ptdJet[leadingJets[0]]/ptdJet[leadingJets[2]];
 var_ptd2byptd3=ptdJet[leadingJets[1]]/ptdJet[leadingJets[2]];
}

var_JA3=areaJet[leadingJets[2]];
var_ptd3=ptdJet[leadingJets[2]];

}


 



if (flav.size()>0 && topJets.size()>0)  {
if (flav[0] == 2 ) {
var_Pt1_b = topJets[0].Pt(); 
var_Eta1_b = topJets[0].Eta(); 
}
else {  
var_Pt1_nonb=topJets[0].Pt();
var_Eta1_nonb=topJets[0].Eta();
}

if (flav[0] == 0 ) {
var_Pt1_q = topJets[0].Pt();
var_Eta1_q=topJets[0].Eta();
}


if (flav[0] == 1 ) {
var_Pt1_c = topJets[0].Pt();
var_Eta1_c=topJets[0].Eta();
}

}

if (topJets.size()>1)
{

 math::XYZTLorentzVector Jet0 = topJets[0];
 math::XYZTLorentzVector Jet1 = topJets[1];
var_M12 =   (Jet0+Jet1).M();


}

if (flav.size()>1 && topJets.size()>1)  {

if (flav[1] == 2 ) {
var_Pt2_b = topJets[1].Pt(); 
var_Eta2_b = topJets[1].Eta(); 
}
else {  
var_Pt2_nonb=topJets[1].Pt();
var_Eta2_nonb=topJets[1].Eta();
}

if (flav[1] == 0 ) {
var_Pt2_q = topJets[1].Pt();
var_Eta2_q=topJets[1].Eta();
}


if (flav[1] == 1 ) {
var_Pt2_c = topJets[1].Pt();
var_Eta2_c=topJets[1].Eta();
}


}


if (flav.size()>2 && topJets.size()>2)  {
if (flav[2] == 2 ) {
var_Pt3_b = topJets[2].Pt(); 
var_Eta3_b = topJets[2].Eta(); 
}
else {  
var_Pt3_nonb=topJets[2].Pt();
var_Eta3_nonb=topJets[2].Eta();
}
}


//std::cout<<"MVAInputVars!!! "<< var_MET<<std::endl;  

//std::cout<<nJetsMax<<std::endl;

  math::XYZTLorentzVector Jetsum(0.,0.,0.,0.);
  
  for(unsigned int i=0; i<nJetsMax; i++) {
    math::XYZTLorentzVector aJet = topJets[i];
    Jetsum += aJet;
    var_sumEt += topJets[i].Pt();
  }
  massalljets = Jetsum.M();

//std::cout<<"MVAInputVars!!! "<<massalljets<<std::endl;  

//std::cout<<"var_sumEt = "<<var_sumEt<<std::endl;
 
//  var_lepeta = -1e10;

//std::cout<<"MVAInputVars!!! "<<var_lepeta<<std::endl;  


//  math::XYZTLorentzVector Met = MET->p4();
  std::vector<double> Etjet;
  std::vector<double> Jetjet;
  double dijetmass;
  var_mindijetmass = 99999.;
  var_maxdijetmass = -1e10;

  for(unsigned int i=0; i<nJetsMax; i++) {
    math::XYZTLorentzVector aJet = topJets[i];
    Etjet.push_back( aJet.Pt());
    for(unsigned int j=i+1; j<nJetsMax; j++) {
      math::XYZTLorentzVector asecJet = topJets[j];
      dijetmass = (aJet+asecJet).M();
      if(dijetmass<var_mindijetmass) var_mindijetmass = dijetmass;
      if(dijetmass>var_maxdijetmass) var_maxdijetmass = dijetmass;
    }
  }

//std::cout<<"var_mindijetmass = "<<var_mindijetmass<<std::endl;


if (Etjet.size()>0)  {
	var_Et1 = Etjet[0];
//   std::cout<<"var_Et1 is "<<var_Et1<<std::endl;
}
   else var_Et1=-1e10;
  var_dphiMETlepton =-1e10;


std::vector<math::XYZTLorentzVector> _vec1;
std::vector<math::XYZTLorentzVector> _vec2;

for(unsigned int i=0; i<nJetsMax; i++ ) _vec1.push_back(topJets[i]);
for(unsigned int i=0; i<nJetsMax; i++ ) _vec2.push_back(topJets[i]);

if (nJetsMax>2)  var_detajet2jet3 = fabs(_vec1[1].Eta()-_vec2[2].Eta());
else var_detajet2jet3=-1e10;
if (nJetsMax>3)  var_detajet3jet4 = fabs(_vec1[2].Eta()-_vec2[3].Eta());
else var_detajet3jet4=-1e10;


//std::cout<<"var_detajet2jet3 = "<<var_detajet2jet3<<std::endl;

 
  var_mindRjetlepton = -1e+10;


///Event shape variables

std::vector<math::XYZVector> _vecXYZ;
for ( unsigned int i=0;i<nJetsMax;i++) _vecXYZ.push_back(topJets[i].Vect());
EventShapeVariables _etsh(_vecXYZ);

///New code! to test
///sphericity, sphericity_boost
/**
	TObjArray *  _jets = new TObjArray();
	_jets->SetOwner(kTRUE); // To 'autodelete' all daughters

	for ( unsigned int i=0;i<nJetsMax;i++) _jets->Add(new TLorentzVector(topJets[i].X(),topJets[i].Y(), topJets[i].Z(), topJets[i].T() ));
        Sphericity * spher = new Sphericity();
        spher->analyze(_jets,2);
	delete _jets;
**/		

if (topJets.size()>0) {
var_sphericity = _etsh.sphericity();
///var_sphericity = spher->sphericity();
var_aplanarity = _etsh.aplanarity();
var_isotropy= _etsh.isotropy();
var_D=_etsh.D();

} 
else {

var_sphericity = -1e10;
var_aplanarity = -1e10;
var_D=-1e10;
}

if (Etjet.size()>1)  var_Et2 = Etjet[1];
   else var_Et2=-1e10;
if (Etjet.size()>2)  var_Et3 = Etjet[2];
   else var_Et3=-1e10;

//   std::cout<<"var_Et2 is "<<var_Et2<<std::endl;
//   std::cout<<"var_Et3 is "<<var_Et3<<std::endl;


if (topJets.size()>0)  var_Eta1 = topJets[0].Eta();
   else var_Eta1=-1e10;
if (topJets.size()>1)  var_Eta2 = topJets[1].Eta();
   else var_Eta2=-1e10;
if (topJets.size()>2)  var_Eta3 = topJets[2].Eta();
   else var_Eta3=-1e10;

//std::cout<<" 1111!!!"<<std::endl;

if (topJets.size()>1)  var_dptjet1jet2 = topJets[0].Pt() - topJets[1].Pt();
else var_dptjet1jet2=-1e10;
if (topJets.size()>2)  var_dptjet1jet3 = topJets[0].Pt() - topJets[2].Pt();
else var_dptjet1jet3=-1e10;
if (topJets.size()>1)  var_djet1jet2pt = (topJets[0]- topJets[1]).Pt();
else var_djet1jet2pt=-1e10;
if (nJetsMax>1)  var_detajet1jet2 = fabs(_vec1[0].Eta()-_vec2[1].Eta());
else var_detajet1jet2=-1e10;



//std::cout<<"var_sphericity = "<<var_sphericity<<std::endl;






var_phijet3_boost12 = -1e10;
var_phijet3_boost13 = -1e10;
var_thetajet1_boost13 = -1e10;
var_thetajet2_boost13 = -1e10;
var_dphijet2jet3_boost13=-1e10;

std::vector<TLorentzVector>_vecs12;
std::vector<TLorentzVector>_vecs13;

if (topJets.size()>1) 
{
///How many objects to be used to construct CM frame
Int_t boostNum=topJets.size();

//for(unsigned int i=0; i<2; i++ ) _vecs12.push_back(TLorentzVector(topJets[i].X(),topJets[i].Y(),topJets[i].Z(),topJets[i].T()));
//TLorentzVector _Higgs = _vecs12[0]+_vecs12[1];
for( int i=0; i< boostNum; i++ ) _vecs12.push_back(TLorentzVector(topJets[i].X(),topJets[i].Y(),topJets[i].Z(),topJets[i].T()));
TLorentzVector _Higgs = _vecs12[0];

for(int i=1; i< boostNum; i++ ) _Higgs+=_vecs12[i];

//TVector3 _CM = -_Higgs.BoostVector();
_vecs12[0].Boost( -_Higgs.BoostVector());
_vecs12[1].Boost( -_Higgs.BoostVector());

var_dthetajet1jet2_boost = fabs(ROOT::Math::VectorUtil::Angle<TLorentzVector,TLorentzVector>(_vecs12[0],_vecs12[1])); 
var_dphijet1jet2_boost = fabs(ROOT::Math::VectorUtil::DeltaPhi<TLorentzVector,TLorentzVector>(_vecs12[0],_vecs12[1])); 

//std::cout<<" _vecs12[0].Boost.X = "<< _vecs12[0].X()<<std::endl; 

_vecs12.clear();
for(unsigned int i=0; i<2; i++ ) _vecs12.push_back(TLorentzVector(topJets[i].X(),topJets[i].Y(),topJets[i].Z(),topJets[i].T()));
_Higgs = _vecs12[0];
for(unsigned int i=1; i<2; i++ ) _Higgs+=_vecs12[i];

_vecs12[0].Boost( -_Higgs.BoostVector());
_vecs12[1].Boost( -_Higgs.BoostVector());


var_thetajet1_boost12 = _vecs12[0].Theta();
//var_thetajet2_boost12 = _vecs12[1].Theta();

if (topJets.size()>2) {
_vecs12.push_back(TLorentzVector(topJets[2].X(),topJets[2].Y(),topJets[2].Z(),topJets[2].T()));
_vecs12[2].Boost( -_Higgs.BoostVector());
var_thetajet3_boost12 = _vecs12[2].Theta();

//var_dphijet1jet3_boost12 =   fabs(ROOT::Math::VectorUtil::DeltaPhi<TLorentzVector,TLorentzVector>(_vecs12[0],_vecs12[2]));
var_dphijet2jet3_boost12 =   fabs(ROOT::Math::VectorUtil::DeltaPhi<TLorentzVector,TLorentzVector>(_vecs12[1],_vecs12[2]));

/// new variables, added 20.07.2013

var_phijet3_boost12 =  _vecs12[2].Phi();
for( int i=0; i< 3; i++ ) _vecs13.push_back(TLorentzVector(topJets[i].X(),topJets[i].Y(),topJets[i].Z(),topJets[i].T()));
TLorentzVector _Higgs13 = _vecs13[0];
_Higgs13+=_vecs13[2];
_vecs13[0].Boost( -_Higgs13.BoostVector());
_vecs13[1].Boost( -_Higgs13.BoostVector());
_vecs13[2].Boost( -_Higgs13.BoostVector());

var_thetajet1_boost13 = _vecs13[0].Theta();
var_thetajet2_boost13 = _vecs13[1].Theta();
var_dphijet2jet3_boost13 =   fabs(ROOT::Math::VectorUtil::DeltaPhi<TLorentzVector,TLorentzVector>(_vecs13[1],_vecs13[2]));
var_phijet3_boost13 =  _vecs13[2].Phi();


}
else 
{
 var_thetajet3_boost12=-1e10;
 var_dphijet2jet3_boost12=-1e10;
}
//var_dphijet1jet2_boost12 =   fabs(ROOT::Math::VectorUtil::DeltaPhi<TLorentzVector,TLorentzVector>(_vecs12[0],_vecs12[1]));





} else {
var_dthetajet1jet2_boost=-1e10;
var_dphijet1jet2_boost = -1e10;
}

/*
std::vector<TLorentzVector>_vecs123;
if (topJets.size()>2) {
for(unsigned int i=0; i< 3; i++ ) _vecs123.push_back(TLorentzVector(topJets[i].X(),topJets[i].Y(),topJets[i].Z(),topJets[i].T()));

TLorentzVector _Higgs  = _vecs123[0];
for(unsigned int i=1; i<3; i++ ) _Higgs+=_vecs123[i];

_vecs123[0].Boost( -_Higgs.BoostVector());
_vecs123[1].Boost( -_Higgs.BoostVector());
_vecs123[2].Boost( -_Higgs.BoostVector());

var_thetajet1_boost123 = _vecs123[0].Theta();
var_thetajet2_boost123 =  _vecs123[1].Theta();


}

*/


if (topJets.size()>1)
var_dphijet1jet2= fabs(ROOT::Math::VectorUtil::DeltaPhi<math::XYZTLorentzVector,math::XYZTLorentzVector>(topJets[0],topJets[1]));
else var_dphijet1jet2=-1e10;

if (topJets.size()>2){
var_dphijet1jet3= fabs(ROOT::Math::VectorUtil::DeltaPhi<math::XYZTLorentzVector,math::XYZTLorentzVector>(topJets[0],topJets[2]));
var_dphijet2jet3= fabs(ROOT::Math::VectorUtil::DeltaPhi<math::XYZTLorentzVector,math::XYZTLorentzVector>(topJets[1],topJets[2]));

}else {
var_dphijet1jet3=-1e10;
var_dphijet2jet3=-1e10;
}

//std::cout<<"22323 !!!"<<std::endl;

_vecs12.clear();
_vecXYZ.clear();
for(unsigned int i=0; i<nJetsMax;i++) _vecs12.push_back(TLorentzVector(topJets[i].X(),topJets[i].Y(),topJets[i].Z(),topJets[i].T()));
TLorentzVector _CM (0.,0.,0.,0.); 
for(unsigned int i=0; i<nJetsMax;i++) _CM+=_vecs12[i];

//cout<<"Before boost"<<endl;

for(unsigned int i=0; i<nJetsMax;i++) _vecs12[i].Boost(-_CM.BoostVector());


for(unsigned int i=0; i<nJetsMax;i++) _vecXYZ.push_back(math::XYZVector(_vecs12[i].X(),_vecs12[i].Y(),_vecs12[i].Z()));
EventShapeVariables _etsh_boost(_vecXYZ);

//cout<<"After boost"<<endl;

///New code! to test
///sphericity, sphericity_boost
/**
        TObjArray *  _jets2 = new TObjArray();
        _jets2->SetOwner(kTRUE); // To 'autodelete' all daughters
        for ( unsigned int i=0;i<nJetsMax;i++) _jets2->Add(new TLorentzVector(_vecs12[i].X(),_vecs12[i].Y(), _vecs12[i].Z(), _vecs12[i].T() ));
  
        spher->analyze(_jets2,2);
	delete _jets2;
**/


if (topJets.size()>0)
var_sphericity_boost = _etsh_boost.sphericity();
///var_sphericity_boost = spher->sphericity();
else
var_sphericity_boost = -1e10;

//std::cout<<"END "<<std::endl;

///if (spher) delete spher;



/// new event shape variables, added 20.07.13
/// about variables http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node214.html
/// http://home.fnal.gov/~mrenna/lutp0613man2/node235.html
var_thrust=-1e10;
var_major=-1e10;
var_minor=-1e10;
var_oblateness=-1e10;

/// about jet broadening  http://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=5&cad=rja&ved=0CEgQFjAE&url=http%3A%2F%2Filcagenda.linearcollider.org%2FmaterialDisplay.py%3FcontribId%3D53%26sessionId%3D10%26materialId%3Dslides%26confId%3D5468&ei=soPqUbXBHIbHtAbC54GABg&usg=AFQjCNE9MCePpZKQ3tcNrnOEU9_lPljvvQ&sig2=Xhd5tQ2H1nfzo6BNBg3aqg
/**
* Jet broadening related event shapes
*/

/**
*  The wide jet broadening
*/
var_Bmax=-1e10;
/**
 *  The narrow jet broadening
*/
var_Bmin=-1e10;
/**
*  The difference of the jet broadenings
*/
var_Bdiff=-1e10;


TObjArray* veclist = new TObjArray();
veclist->SetOwner(kTRUE); // To 'autodelete' all daughters
Float_t Pvec3[3];

for ( unsigned int i=0;i<nJetsMax;i++) {
Pvec3[0]=topJets[i].X();
Pvec3[1]=topJets[i].Y();
Pvec3[2]=topJets[i].Z();
TVector3* vec3 = new TVector3(Pvec3);
veclist->Add(vec3);
}
EventShape eshape;
//  Do Thrust finding
eshape.setPartList(veclist); // Input particle 3-vector list

var_thrust=eshape.thrust().X();
var_major=eshape.thrust().Y();
var_minor=eshape.thrust().Z();

var_oblateness=eshape.oblateness();

/// jet broadening taken from 
/// http://herwig.hepforge.org/svn/tags/release-2-4-0/Analysis/EventShapes.h
///http://herwig.hepforge.org/svn/tags/release-2-4-0/Analysis/EventShapes.cc

TVector3  thrustaxis = eshape.thrustAxis();
//TLorentzVector pos, neg;
double epos=0e0;
double eneg=0e0;
double pden=0e0;
for ( unsigned int i=0;i<nJetsMax;i++) {
Pvec3[0]=topJets[i].X();
Pvec3[1]=topJets[i].Y();
Pvec3[2]=topJets[i].Z();
TVector3  _vec3(Pvec3);
if(_vec3 * thrustaxis > 0e0)
	{
//      pos += topJets[i];
      epos += _vec3.Cross(thrustaxis).Mag();

	}
      else
	{
//      neg += topJets[i];
      eneg += _vec3.Cross(thrustaxis).Mag();

   }
      pden += _vec3.Mag();
} /// for nJetsMax


if (epos<eneg) 
{
double tmp_var;
tmp_var=epos;
epos=eneg;
eneg=tmp_var;
} 


if ( pden!=0e0) {
  var_Bmax  = 0.5*epos/pden;
  var_Bmin  = 0.5*eneg/pden;
}


if ( var_Bmax>-1e10 && var_Bmin>-1e10) var_Bdiff=var_Bmax-var_Bmin;


delete veclist;


}




double MVAInputVars::DeltaPhi(math::XYZTLorentzVector v1, math::XYZTLorentzVector v2)
{
  double dPhi = fabs(v1.Phi() - v2.Phi());
  if (dPhi > TMath::Pi()) dPhi =  2*TMath::Pi() - dPhi;
  return dPhi;
}

double MVAInputVars::DeltaR(math::XYZTLorentzVector v1, math::XYZTLorentzVector v2)
{
  double dPhi = DeltaPhi(v1,v2);
  double dR = TMath::Sqrt((v1.Eta()-v2.Eta())*(v1.Eta()-v2.Eta())+dPhi*dPhi);
  return dR;
}

MVAInputVars::~MVAInputVars() 
{
}


#endif

#endif

