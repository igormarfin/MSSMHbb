#ifndef PileUpProducer_H
#define PileUpProducer_H

//#define _debug
//#define _use_seed

#include "TMath.h"
#include "TRandom.h"
#include "TObjString.h"
#include "TString.h"
#include "TClass.h"
#include "TList.h"

#include <iostream>
#include <map>
#include <vector>
#include <strstream>



#include "TObject.h"
#include "TNamed.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TH1F.h"

///Gamma and Landau distribution
#include "TF1.h" 
#include "TMath.h"
#include "TRandom.h"
#include "TROOT.h"

using namespace std;

class PileUpProducer : public TNamed 
{
public:


PileUpProducer():TNamed("PileUpProducer",(const char *)("")) {

dataHist_=0;
scenario_=-1;
}

PileUpProducer(TString name, TString header):TNamed(name,header) {
dataHist_=0;
scenario_=-1;
}
	

int shootNumVert();


std::vector<double> PileUpReweight(TH1* mcNPU, TH1 *dataNPU);


std::vector<double> PileUpReweight(TH1 *dataNPU);

PileUpProducer(TString name, TString title,TString dataHistPath,TString dataHistName,int scenario=0, double average=0,double sigma=0,int seed=0);

PileUpProducer(TString name, TString title,TString dataHistPath,TString MCHistPath,TString dataHistName,int scenario=0, double average=0,double sigma=0,int seed=0);

PileUpProducer(TString name, TString title,TH1F * dataHist, int scenario=0, double average=0,double sigma=0,int seed=0);




std::vector<double>  &  GetWeights() { return pileupwgt_;}  

double GetWeight(int nPU);

~PileUpProducer();


    



private:
TString dataHistPath_;
TString mcHistPath_;
TString dataHistName_;
TH1F * dataHist_;
TH1F * mcHist_;
std::vector<double> pileupwgt_;
int scenario_;
double average_;
double sigma_;
int seed_;


public:


 	ClassDef(PileUpProducer,2);
};



 #ifdef __MAKECINT__
ClassImp(PileUpProducer)
#endif


#if defined(PileUpProducer_CC)


PileUpProducer::PileUpProducer(TString name, TString title,TH1F * dataHist, int scenario, double average,double sigma,int seed):TNamed(name,title)
{


scenario_=scenario;
average_=average;
sigma_=sigma;

seed_=seed;

#ifdef _use_seed
gRandom->SetSeed(seed_);
#endif


dataHist_=dataHist;
if (dataHist_) dataHist_->SetDirectory(gROOT);

#ifdef _debug

if (dataHist_ ) std::cout<<"Histo Name_ = "<<dataHist_->GetName()<<std::endl;
std::cout<<" scenario_ = "<<scenario_<<std::endl;

#endif

pileupwgt_ = PileUpReweight(dataHist_);


}




PileUpProducer::PileUpProducer(TString name, TString title,TString dataHistPath,TString dataHistName,int scenario, double average,double sigma,int seed):TNamed(name,title)
{

dataHistPath_=dataHistPath;
dataHistName_=dataHistName;
scenario_=scenario;
average_=average;
sigma_=sigma;

seed_=seed;

#ifdef _use_seed
gRandom->SetSeed(seed_);
#endif

TFile file0(dataHistPath_.Data());	

dataHist_=(TH1F *)file0.Get(dataHistName_.Data());
if (dataHist_) dataHist_->SetDirectory(gROOT);

#ifdef _debug
if (dataHist_ ) std::cout<<"Histo Name_ = "<<dataHist_->GetName()<<std::endl;
std::cout<<" scenario_ = "<<scenario_<<std::endl;

#endif

if (dataHist_) pileupwgt_ = PileUpReweight(dataHist_);


}

PileUpProducer::PileUpProducer(TString name, TString title,TString dataHistPath,TString mcHistPath,TString dataHistName,int scenario, double average,double sigma,int seed):TNamed(name,title)
{

dataHistPath_=dataHistPath;
mcHistPath_=mcHistPath;
dataHistName_=dataHistName;
scenario_=scenario;
average_=average;
sigma_=sigma;

seed_=seed;

#ifdef _use_seed
gRandom->SetSeed(seed_);
#endif

TFile file0(dataHistPath_.Data());	

dataHist_=(TH1F *)file0.Get(dataHistName_.Data());
if (dataHist_) dataHist_->SetDirectory(gROOT);

TFile file1(mcHistPath_.Data());	

mcHist_=(TH1F *)file1.Get(dataHistName_.Data());
if (mcHist_) mcHist_->SetDirectory(gROOT);




#ifdef _debug
if (dataHist_ ) std::cout<<"Histo Name_ = "<<dataHist_->GetName()<<std::endl;
std::cout<<" scenario_ = "<<scenario_<<std::endl;

#endif

//if (dataHist_ && mcHist_) pileupwgt_ = PileUpReweight(mcHist_,dataHist_);


}


PileUpProducer::~PileUpProducer() { }



int PileUpProducer::shootNumVert()
{

int nzero_crossing = -1;


/// Taken from
///http://cmslxr.fnal.gov/lxr/source/Mixing/Base/src/BMixingModule.cc#028
if (scenario_==0) 
{

const double npu_probs[25] = {
0.0698146584,0.0698146584,0.0698146584,0.0698146584,
0.0698146584,0.0698146584,0.0698146584,0.0698146584,
0.0698146584,0.0698146584,0.0698146584,0.0630151648,
0.0526654164,0.0402754482,0.0292988928,0.0194384503,
0.0122016783,0.007207042,0.004003637,0.0020278322,
0.0010739954,0.0004595759,0.0002229748,0.0001028162,
4.58337152809607E-05
};

int xmin = (int) 0;
int xmax = (int) 25;  // need upper edge to be one beyond last value (npu_probs.size()-1)
int numBins = (int) 25;

TString nameHist="hprob";


TH1F *hprob = new TH1F(nameHist,"Histo from the user's probability function",numBins,xmin,xmax); 
for (int j=0; j < numBins ; j++)              hprob->Fill(j+0.5,npu_probs[j]); // assuming integer values for the bins, fill bin centers, not edges 
              


///Taken from
///http://cmslxr.fnal.gov/lxr/source/Mixing/Base/src/PileUp.cc

double d = hprob->GetRandom();
nzero_crossing =  int(d);

delete hprob;


}

///Taken from
///http://cmslxr.fnal.gov/lxr/source/Mixing/Base/src/PileUp.cc
///Poisson
if (scenario_==1)
{

nzero_crossing = (int) gRandom->Poisson(average_);

}

if (scenario_==2)
{

nzero_crossing = (int) gRandom->Uniform(0,25);

}


if (scenario_==3)
{

 
        TF1 * logNormFun=new TF1("logNormFun","TMath::LogNormal(x,[0],[1],[2])",0,25);
        logNormFun->SetParameters(sigma_,-average_,1.0);

        nzero_crossing = (int) logNormFun->GetRandom();  

        delete logNormFun;


}

if (scenario_==4)
{


        TF1 * gammaFun=new TF1("gammaFun","TMath::GammaDist(x,[0],[1],[2])",0,25);
        gammaFun->SetParameters(sigma_,0,1.0);

        nzero_crossing = (int) gammaFun->GetRandom();

        delete gammaFun;
}

if (scenario_==5)
{


        TF1 * landauFun=new TF1("landauFun","TMath::Landau(x,[0],[1],[2])",0,25);
        landauFun->SetParameters(average_,sigma_,kFALSE);

        nzero_crossing = (int) landauFun->GetRandom();

        delete landauFun;
}





return nzero_crossing;
}

/**

	The code was taken from Alexei's  Utilities/src/PUWeight.cc


**/


double PileUpProducer::GetWeight(int nPU)
{

double  result =1;
if (!(dataHist_ && mcHist_) ) {/*cout<<"Null pointers ; return"<<endl;*/ return result;}


  /*
      The applied rebinning is not really necessary as the
      reweighting function uses FindBin(). The same applies to the
      range (0,25) for the MC histogram. However, the histograms for
      
      both data and MC have the same binning to avoid edge/binning
      effects. A coarse binning should also limit the effect of the
      missing random seed during the pile-up creation process.
    */

Int_t _datanBins =dataHist_->GetNbinsX();
Int_t _mcnBins= mcHist_->GetNbinsX();

if (_datanBins ==0 || _mcnBins==0) { cout<<" Empty TH1s? ; return"<<endl; return result; }

Int_t _minNbins =-1;
if (_datanBins>_mcnBins) _minNbins=_mcnBins;
if (_mcnBins>_datanBins) _minNbins=_datanBins;

if (_minNbins>0)
if (!(_datanBins%_minNbins == 0 && _mcnBins%_minNbins==0)) {cout<<"Can't make rebining ; return"<<endl; return result;}


if (_datanBins>_mcnBins)
    dataHist_->Rebin(_datanBins/_minNbins);

if (_mcnBins>_datanBins)
     mcHist_->Rebin(_mcnBins/_minNbins);

    mcHist_->Scale(1./mcHist_->Integral());
    dataHist_->Scale(1./dataHist_->Integral());

//    for (int idx1 = 1; idx1 < mcHist_->GetNbinsX(); idx1++)
//      std::cout << idx1 << "\t" << mcHist_->GetBinContent(idx1) << "\t" << mcHist_->GetBinContent(idx1) << "\n";

  float mean = float(nPU);
  float weightMC = mcHist_->GetBinContent(mcHist_->FindBin(mean));
  float weightData = dataHist_->GetBinContent(dataHist_->FindBin(mean));

  if (!weightMC==0.)
    result = (double) weightData/ (double)weightMC;
  else result=1;	


return result;
}




/**
		mcNPU && dataNPU are TH1 histograms supposed to be non-normalized to 1

		They must be defined with equal parameters (xmin,xmax,numOfBins) and several conditions:
			1) xmin == 0
			2) xmax == numOfBins


**/




std::vector<double> PileUpProducer::PileUpReweight(TH1* mcNPU, TH1 *dataNPU)
{
std::vector<double> result;

if (!(mcNPU && dataNPU) ) {cout<<"Null pointers ; return"<<endl; return result;}

const int _datanPU =dataNPU->GetNbinsX();
const int _mcnPU =mcNPU->GetNbinsX();

int _dataxmin = (int) dataNPU->GetBinLowEdge(1); /// must be 0
int _dataxmax = (int) (dataNPU->GetBinLowEdge(_datanPU) + dataNPU->GetBinWidth(_datanPU));  /// must be maximal value of PU verticies

int _mcxmin = (int) mcNPU->GetBinLowEdge(1); /// must be 0
int _mcxmax = (int) (mcNPU->GetBinLowEdge(_mcnPU) + mcNPU->GetBinWidth(_mcnPU));  /// must be maximal value of PU verticies

if ( _datanPU != _dataxmax ) {cout<<"Problem in dataPU: NumOfBins!=MaxNPU; return"<<endl; return result;}
if ( _mcnPU != _mcxmax ) {cout<<"Problem in mcPU: NumOfBins!=MaxNPU; return"<<endl; return result;}

if ( 0 != _dataxmin ) {cout<<"Problem in dataPU: MinNPU != 0; return"<<endl; return result;}
if ( 0 != _mcxmin ) {cout<<"Problem in mcPU: MinNPU != 0; return"<<endl; return result;}

if ( _datanPU != _mcnPU ) { cout <<"Problem in consistency: NumOfBinsDataPU!= NumOfBinsMCPU; return"<<endl;return result;}

///Normalize MCPU

mcNPU->Scale(1e0/mcNPU->Integral());

 double _s = 0.0;

 for(int _npu=0; _npu<_dataxmax; ++_npu){
double _npu_data = dataNPU->GetBinContent(dataNPU->GetXaxis()->FindBin(_npu));
double _npu_mc = mcNPU->GetBinContent(mcNPU->GetXaxis()->FindBin(_npu));

 result.push_back(_npu_data/_npu_mc);
        _s += _npu_data;

}



 for(int _npu=0; _npu<_datanPU; ++_npu) result[_npu] /= _s;

return result; 
}


std::vector<double> PileUpProducer::PileUpReweight(TH1 *dataNPU)
{

///http://cmslxr.fnal.gov/lxr/source/SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_50ns_PoissonOOT.py
const double npu_probs[25] = {
0.0698146584,0.0698146584,0.0698146584,0.0698146584,
0.0698146584,0.0698146584,0.0698146584,0.0698146584,
0.0698146584,0.0698146584,0.0698146584,0.0630151648,
0.0526654164,0.0402754482,0.0292988928,0.0194384503,
0.0122016783,0.007207042,0.004003637,0.0020278322,
0.0010739954,0.0004595759,0.0002229748,0.0001028162,
4.58337152809607E-05      
};

    std::vector<double> result(25);
    double s = 0.0;
    for(int npu=0; npu<25; ++npu){

        double npu_estimated ;
        if (dataNPU)    npu_estimated = dataNPU->GetBinContent(dataNPU->GetXaxis()->FindBin(npu));                              
        else npu_estimated=npu_probs[npu];

        result[npu] = npu_estimated / npu_probs[npu];
        s += npu_estimated;
    }
    // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
    for(int npu=0; npu<25; ++npu) result[npu] /= s;
    
    return result;

}





#endif

#endif
