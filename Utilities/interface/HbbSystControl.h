#ifndef HbbSystControl_h
#define HbbSystControl_h



#include "TObject.h"
#include "TNamed.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <vector>

#include "HbbSyst.h"
#include "HbbNtuple.h"



/**

	Performs changing values of some ntuple's variable in accordance to taken systematics


**/


/**
	type :
		-1 --nothing to do, to get cenrtal value
		0 -- JSE uncert,
		1 -- SF central,
        2 -- SF uncert b/c,
		3 -- SF uncert udsg,
		4 -- JER uncert 

/// added to 2012 fat jets
        5 -- SF uncert cc/bb
	    

	shift :	
		>0 -- UP change (+1sigma*shift)
		<0 -- DOWN change (-1sigma*shift)
		==0 returns central value (for example, SF)

	tagger :

	      // 0 : TCHPT
	      // 1 : TCHP6
	      // 2 : CSVT
	      // 3 : SSVHPT


**/


using namespace std;

class HbbSystControl:public HbbSyst
{
	

public:
	HbbSystControl():_wasJESapplied(false),_wasJERapplied(false) { 
	 nbtag=4;
	sbtag = new TString [nbtag];
	sbtag[0]="TCHPT";
	sbtag[1]="TCHP6";
	sbtag[2]="CSVT";
	sbtag[3]="SSVHPT";

	}

	~HbbSystControl() { delete [] sbtag;}

	void operator()(int type, double shift, int itagger, int jetFlavor, double jetPtX, double jetEtaX, double & _SF);
       std::vector<float> genJetMatched(float jetPt, float jetEta, float jetPhi, float dRMin);


void SysJES( double shift);
void SysJER( double shift);

double GetSF(TString tagger, int jetFlavor,  double jetPtX, double jetEtaX);
double SysSFbc(TString tagger, int jetFlavor, double shift , double jetPtX, double jetEtaX);
double SysSFudsg(TString tagger, int jetFlavor, double shift, double jetPtX, double jetEtaX);
double SysSFbbcc(TString tagger, int jetFlavor, double shift , double jetPtX, double jetEtaX);
void Reset() { /* _wasJESapplied=false; _wasJERapplied=false; */ }

private:

bool _wasJESapplied; // to know if JES was applied
bool _wasJERapplied; // to know if JER was applied

int nbtag;
TString * sbtag;
};





#if defined(HbbSystControl_CC)




void  HbbSystControl::operator() (int type, double shift, int itagger, int jetFlavor,  double jetPtX, double jetEtaX,double & _SF )
	{

_SF=1.0;


		if (type<0 ) return;


		if (type == 0) {
/*
			if (_wasJESapplied) return;
			else _wasJESapplied=true;
*/



//			cout<<"JES shift is started"<<endl;
//			cout<<"JES shift = "<<shift<<endl;

			    for (int i=0;i<numberOfJets;i++)
			   {
			           double jesUncertUp = 0;

			           if (shift>0)  jesUncertUp = getJESuncertaintyUp(ptJet[i],etaJet[i]);
  		  	           else jesUncertUp = getJESuncertaintyDown(ptJet[i],etaJet[i]);

//				cout<<"jesUncertScale ["<<i<<"] = "<<jesUncertUp<<endl;

				 if (shift!=0) {	
				  ptJet[i]*=1+shift*jesUncertUp;
				  pxJet[i]*=1+shift*jesUncertUp;
				  pyJet[i]*=1+shift*jesUncertUp;
				  pzJet[i]*=1+shift*jesUncertUp;		  
				  energyJet[i]*=1+shift*jesUncertUp;
				}
		   	   }
		} ///type 0

	
		 if (type == 1 && jetFlavor>=0 && itagger>=0) {

                 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);

                }

		


		if (type == 2 && jetFlavor>=1 && jetFlavor<3 && itagger>=0) {
		 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);

		  if (shift!=0)  {
       	        _SF += shift*getSFuncertainty(jetPtX,jetEtaX,jetFlavor,itagger,shift>0); 
		}

		}

		if (type == 3 && jetFlavor==0 && itagger>=0) {
			
		 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);

		  if (shift!=0)  {
       	        _SF += shift*getSFuncertainty(jetPtX,jetEtaX,jetFlavor,itagger,shift>0); 
		}
		}

		if (type == 5 && jetFlavor>2 && itagger>=0) {
			
		 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);

		  if (shift!=0)  {
       	        _SF += shift*getSFuncertainty(jetPtX,jetEtaX,jetFlavor,itagger,shift>0); 
		}
		}

      
      
		if (type == 4) 
      {
/*
			if (_wasJERapplied) return;
			else _wasJERapplied=true;
*/

//cout<<"JER  is started"<<endl;
//cout<<"JER shift = "<<shift<<endl;
         for (int i=0;i<numberOfJets;i++)
         {
            float jerUncert = 0.;
            if (shift>0)  jerUncert = getJERuncertaintyUp(ptJet[i],etaJet[i]);
            else jerUncert = getJERuncertaintyDown(ptJet[i],etaJet[i]);
            
            float jerSF = getSFJER(ptJet[i],etaJet[i]);
//std::cout << " *** HbbSystControl:  JER scale factor = " << jerSF << " for a jet["<<i<<"] with pt="<<ptJet[i]<<" , with eta = " << etaJet[i] << std::endl;
            
            std::vector<float> genJet = genJetMatched( ptJet[i], etaJet[i] , phiJet[i], 0.5 );
            
            if ( genJet.size() == 0 ) 
            {
//               std::cout << " No match!!! Event  = " << event << std::endl;
               continue; // no match, no change
            }
            
            float genJetPt = genJet[0];
            float genJetPx = genJet[1];
            float genJetPy = genJet[2];
            float genJetPz = genJet[3];
            float genJetEnergy = genJet[4];
            
            float modPtJet;
            float modPxJet;
            float modPyJet;
            float modPzJet;
            float modEnergyJet;

	
//cout<<"ptJet-genJetPt ["<<i<<"] = "<<ptJet[i]-genJetPt<<endl;
//cout<<"fbas(pxJet)-fabs(genJetPt) [" <<i<<"]= "<<fabs(pxJet[i])-fabs(genJetPx)<<endl;
	                  
            if (shift!=0) {
//std::cout << " *** HbbSystControl: jet PT = " << ptJet[i] << ", gen PT = " << genJetPt  << ", jet ETA = " << etaJet[i] << ", JER uncertainty = " << jerUncert << std::endl;
//std::cout << "SF variation" << std::endl;
               modPtJet     = genJetPt      + (jerSF+shift*jerUncert)*(ptJet[i]-genJetPt);
               modPxJet     = fabs(genJetPx)+ (jerSF+shift*jerUncert)*(fabs(pxJet[i])-fabs(genJetPx));
               modPyJet     = fabs(genJetPy)+ (jerSF+shift*jerUncert)*(fabs(pyJet[i])-fabs(genJetPy));
               modPzJet     = fabs(genJetPz)+ (jerSF+shift*jerUncert)*(fabs(pzJet[i])-fabs(genJetPz));
               modEnergyJet = genJetEnergy  + (jerSF+shift*jerUncert)*(energyJet[i]-genJetEnergy);
            } else {
//std::cout << "SF only" << std::endl;
               modPtJet     = genJetPt      + jerSF*(ptJet[i]-genJetPt);
               modPxJet     = fabs(genJetPx)+ jerSF*(fabs(pxJet[i])-fabs(genJetPx));
               modPyJet     = fabs(genJetPy)+ jerSF*(fabs(pyJet[i])-fabs(genJetPy));
               modPzJet     = fabs(genJetPz)+ jerSF*(fabs(pzJet[i])-fabs(genJetPz));
               modEnergyJet = genJetEnergy+jerSF*(energyJet[i]-genJetEnergy);
//std::cout << "mod jet pt = " << modPtJet << ",   " << ptJet[i] << std::endl;
            }
            float zero = 0.;
            float pxSign = pxJet[i]/fabs(pxJet[i]);
            float pySign = pyJet[i]/fabs(pyJet[i]);
            float pzSign = pzJet[i]/fabs(pzJet[i]);
            ptJet[i] = std::max(zero, modPtJet);
            pxJet[i] = pxSign*modPxJet;
            pyJet[i] = pySign*modPyJet;
            pzJet[i] = pzSign*modPzJet;
            energyJet[i] = std::max(zero, modEnergyJet);
//std::cout << " *** HbbSystControl: mod PT = " << ptJet[i] << std::endl;
        
            
         }
      }


return;
}


 

std::vector<float> HbbSystControl::genJetMatched( float jetPt, float jetEta, float jetPhi, float dRMin )
{
   std::vector<float> genJet;
   float pi = acos(-1);
   for ( int j = 0; j < numberOfGenJets; ++j )
   {
      float dEta = jetEta - etaGenJet[j];
      float dPhi = jetPhi - phiGenJet[j];
      if ( dPhi > pi )  dPhi -= 2.*pi;
      if ( dPhi <= -pi ) dPhi += 2.*pi;
      float dR = sqrt(dEta*dEta + dPhi*dPhi);
      if ( dR < dRMin )
      {
         genJet.push_back(ptGenJet[j]);
         genJet.push_back(pxGenJet[j]);
         genJet.push_back(pyGenJet[j]);
         genJet.push_back(pzGenJet[j]);
         genJet.push_back(energyGenJet[j]);
         break;
      }
   }
   return genJet;
   
}

void HbbSystControl::SysJES( double shift){
/*
                        if (_wasJESapplied) return;
                        else _wasJESapplied=true;
*/

                            for (int i=0;i<numberOfJets;i++)
                           {
                                   double jesUncertUp = 0;

                                   if (shift>0)  jesUncertUp = getJESuncertaintyUp(ptJet[i],etaJet[i]);
                                   else jesUncertUp = getJESuncertaintyDown(ptJet[i],etaJet[i]);

                                 if (shift!=0) {
                                  ptJet[i]*=1+shift*jesUncertUp;
                                  pxJet[i]*=1+shift*jesUncertUp;
                                  pyJet[i]*=1+shift*jesUncertUp;
                                  pzJet[i]*=1+shift*jesUncertUp;
                                  energyJet[i]*=1+shift*jesUncertUp;
                                }
                           }

return;
}



void HbbSystControl::SysJER( double shift)
{
/*
		if (_wasJERapplied) return;
		else _wasJERapplied=true;
*/




//std::cout<<"JER  is started"<<endl;
//std::cout<<"JER shift = "<<shift<<endl;


         for (int i=0;i<numberOfJets;i++)
         {
            float jerUncert = 0.;
            if (shift>0)  jerUncert = getJERuncertaintyUp(ptJet[i],etaJet[i]);
            else jerUncert = getJERuncertaintyDown(ptJet[i],etaJet[i]);
            
            float jerSF = getSFJER(ptJet[i],etaJet[i]);
//std::cout << " *** HbbSystControl:  JER scale factor = " << jerSF << " for a jet["<<i<<"] with pt="<<ptJet[i]<<" , with eta = " << etaJet[i] << std::endl;
            
            std::vector<float> genJet = genJetMatched( ptJet[i], etaJet[i] , phiJet[i], 0.5 );
            
            if ( genJet.size() == 0 ) 
            {
//               std::cout << " No match!!! Event  = " << event << std::endl;
               continue; // no match, no change
            }
            
            float genJetPt = genJet[0];
            float genJetPx = genJet[1];
            float genJetPy = genJet[2];
            float genJetPz = genJet[3];
            float genJetEnergy = genJet[4];
            
            float modPtJet;
            float modPxJet;
            float modPyJet;
            float modPzJet;
            float modEnergyJet;


//std::cout<<"ptJet-genJetPt ["<<i<<"] = "<<ptJet[i]-genJetPt<<endl;
//std::cout<<"fbas(pxJet)-fabs(genJetPt) [" <<i<<"]= "<<fabs(pxJet[i])-fabs(genJetPx)<<endl;
	
	                  
            if (shift!=0) {
//std::cout << "SF variation" << std::endl;
//std::cout << " *** HbbSystControl: jet PT = " << ptJet[i] << ", gen PT = " << genJetPt  << ", jet ETA = " << etaJet[i] << ", JER uncertainty = " << jerUncert << std::endl;

               modPtJet     = genJetPt      + (jerSF+shift*jerUncert)*(ptJet[i]-genJetPt);
               modPxJet     = fabs(genJetPx)+ (jerSF+shift*jerUncert)*(fabs(pxJet[i])-fabs(genJetPx));
               modPyJet     = fabs(genJetPy)+ (jerSF+shift*jerUncert)*(fabs(pyJet[i])-fabs(genJetPy));
               modPzJet     = fabs(genJetPz)+ (jerSF+shift*jerUncert)*(fabs(pzJet[i])-fabs(genJetPz));
               modEnergyJet = genJetEnergy  + (jerSF+shift*jerUncert)*(energyJet[i]-genJetEnergy);
            } else {
               modPtJet     = genJetPt      + jerSF*(ptJet[i]-genJetPt);
               modPxJet     = fabs(genJetPx)+ jerSF*(fabs(pxJet[i])-fabs(genJetPx));
               modPyJet     = fabs(genJetPy)+ jerSF*(fabs(pyJet[i])-fabs(genJetPy));
               modPzJet     = fabs(genJetPz)+ jerSF*(fabs(pzJet[i])-fabs(genJetPz));
               modEnergyJet = genJetEnergy+jerSF*(energyJet[i]-genJetEnergy);
            }
            float zero = 0.;
            float pxSign = pxJet[i]/fabs(pxJet[i]);
            float pySign = pyJet[i]/fabs(pyJet[i]);
            float pzSign = pzJet[i]/fabs(pzJet[i]);
            ptJet[i] = std::max(zero, modPtJet);
            pxJet[i] = pxSign*modPxJet;
            pyJet[i] = pySign*modPyJet;
            pzJet[i] = pzSign*modPzJet;
            energyJet[i] = std::max(zero, modEnergyJet);
//std::cout << " *** HbbSystControl: mod PT = " << ptJet[i] << std::endl;
        
            
         } ///for 



return;

}


double HbbSystControl::GetSF(TString tagger, int jetFlavor,  double jetPtX, double jetEtaX)
{


int itagger=(int) std::distance(sbtag,std::find( sbtag, sbtag+nbtag,tagger));


double _SF=1e0;
                 if ( jetFlavor>=0 && itagger>=0) 
                 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);

                
return _SF;


}


double HbbSystControl::SysSFbc(TString tagger, int jetFlavor, double shift , double jetPtX, double jetEtaX){

int itagger=(int) std::distance(sbtag,std::find( sbtag, sbtag+nbtag,tagger));


double _SF=1e0;
                if (jetFlavor>=1 && jetFlavor<3 && itagger>=0) {
                 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);
                 if (shift!=0) _SF += shift*getSFuncertainty(jetPtX,jetEtaX,jetFlavor,itagger,shift>0);

                }


return _SF;


}

double HbbSystControl::SysSFudsg(TString tagger, int jetFlavor, double shift, double jetPtX, double jetEtaX){

int itagger=(int) std::distance(sbtag,std::find( sbtag, sbtag+nbtag,tagger));

double _SF=1e0;
                if (jetFlavor==0 && itagger>=0) {

                 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);
                if (shift!=0) _SF += shift*getSFuncertainty(jetPtX,jetEtaX,jetFlavor,itagger,shift>0);
                
                }


return _SF;
}



double HbbSystControl::SysSFbbcc(TString tagger, int jetFlavor, double shift , double jetPtX, double jetEtaX){

int itagger=(int) std::distance(sbtag,std::find( sbtag, sbtag+nbtag,tagger));


double _SF=1e0;
                if (jetFlavor>2 && itagger>=0) {
                 _SF = getSFtag(jetPtX,jetEtaX,jetFlavor,itagger);
                 if (shift!=0) _SF += shift*getSFuncertainty(jetPtX,jetEtaX,jetFlavor,itagger,shift>0);

                }


return _SF;


}



#endif

#endif
