#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <THStack.h>
#include <TMinuit.h>
#include <TMath.h>

/// systematics read-out
#include <map>
#include <TString.h>
#include <TObjArray.h>
#include <TSystem.h>

TCanvas* canvas;
class templateFitter;
templateFitter* gTemplateFitter;  // pointer needs to be global, to be called in fcn

 const int ncateg = 3;
  const int ncorr=2;
  const int nfc = 3;
  const int nbtag = 4;

///_SYST_BEGIN
typedef struct{ TH2F * _array[nfc][nbtag][ncateg];} mjjEbtTemplateType;
typedef struct{ TH2F * _array[nbtag];} mjjEbtTemplateTypeSgn;
///_SYST_END


void setTDRStyle() {

   TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

// For the canvas:
   tdrStyle->SetCanvasBorderMode(0);
   tdrStyle->SetCanvasColor(kWhite);
   tdrStyle->SetCanvasDefH(600); //Height of canvas
   tdrStyle->SetCanvasDefW(600); //Width of canvas
   tdrStyle->SetCanvasDefX(0);   //POsition on screen
   tdrStyle->SetCanvasDefY(0);

// For the Pad:
   tdrStyle->SetPadBorderMode(0);
   // tdrStyle->SetPadBorderSize(Width_t size = 1);
   tdrStyle->SetPadColor(kWhite);
   tdrStyle->SetPadGridX(false);
   tdrStyle->SetPadGridY(false);
   tdrStyle->SetGridColor(0);
   tdrStyle->SetGridStyle(3);
   tdrStyle->SetGridWidth(1);

// For the frame:
   tdrStyle->SetFrameBorderMode(0);
   tdrStyle->SetFrameBorderSize(1);
   tdrStyle->SetFrameFillColor(0);
   tdrStyle->SetFrameFillStyle(0);
   tdrStyle->SetFrameLineColor(1);
   tdrStyle->SetFrameLineStyle(1);
   tdrStyle->SetFrameLineWidth(1);

// For the histo:
   // tdrStyle->SetHistFillColor(1);
   // tdrStyle->SetHistFillStyle(0);
   tdrStyle->SetHistLineColor(1);
   tdrStyle->SetHistLineStyle(0);
   tdrStyle->SetHistLineWidth(1);
   // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
   // tdrStyle->SetNumberContours(Int_t number = 20);

//   tdrStyle->SetEndErrorSize(2);
//   tdrStyle->SetErrorMarker(20);
//   tdrStyle->SetErrorX(0.);

   tdrStyle->SetMarkerStyle(20);

//For the fit/function:
   tdrStyle->SetOptFit(1);
   tdrStyle->SetFitFormat("5.4g");
   tdrStyle->SetFuncColor(2);
   tdrStyle->SetFuncStyle(1);
   tdrStyle->SetFuncWidth(1);

//For the date:
   tdrStyle->SetOptDate(0);
   // tdrStyle->SetDateX(Float_t x = 0.01);
   // tdrStyle->SetDateY(Float_t y = 0.01);

// // For the statistics box:
//   tdrStyle->SetOptFile(0);
//   tdrStyle->SetOptStat(0); // To display the mean and RMS:    
   tdrStyle->SetOptStat("mrei");
   tdrStyle->SetStatColor(kWhite);
   tdrStyle->SetStatFont(42);
    tdrStyle->SetStatFontSize(0.025);
   tdrStyle->SetStatTextColor(1);
   tdrStyle->SetStatFormat("6.4g");
   tdrStyle->SetStatBorderSize(1);
   tdrStyle->SetStatH(0.1);
   tdrStyle->SetStatW(0.15);
   // tdrStyle->SetStatStyle(Style_t style = 1001);
   // tdrStyle->SetStatX(Float_t x = 0);
   // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
   tdrStyle->SetPadTopMargin(0.05);
   tdrStyle->SetPadBottomMargin(0.13);
   tdrStyle->SetPadLeftMargin(0.16);
   tdrStyle->SetPadRightMargin(0.02);

// For the Global title:

   tdrStyle->SetOptTitle(0);
   tdrStyle->SetTitleFont(42);
   tdrStyle->SetTitleColor(1);
   tdrStyle->SetTitleTextColor(1);
   tdrStyle->SetTitleFillColor(10);
   tdrStyle->SetTitleFontSize(0.05);
   // tdrStyle->SetTitleH(0); // Set the height of the title box
   // tdrStyle->SetTitleW(0); // Set the width of the title box
   // tdrStyle->SetTitleX(0); // Set the position of the title box
   // tdrStyle->SetTitleY(0.985); // Set the position of the title box
   // tdrStyle->SetTitleStyle(Style_t style = 1001);
   // tdrStyle->SetTitleBorderSize(2);

// For the axis titles:

   tdrStyle->SetTitleColor(1, "XYZ");
   tdrStyle->SetTitleFont(42, "XYZ");
   //tdrStyle->SetTitleSize(0.06, "XYZ");
   tdrStyle->SetTitleSize(0.052, "XYZ");
   // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to  set the size?
   // tdrStyle->SetTitleYSize(Float_t size = 0.02);
   //tdrStyle->SetTitleXOffset(0.9);
   tdrStyle->SetTitleXOffset(0.875);
   //tdrStyle->SetTitleYOffset(1.25);
   tdrStyle->SetTitleYOffset(1.05);
   // tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the  Offset

// For the axis labels:

   tdrStyle->SetLabelColor(1, "XYZ");
   tdrStyle->SetLabelFont(42, "XYZ");
   tdrStyle->SetLabelOffset(0.007, "XYZ");
   tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

   tdrStyle->SetAxisColor(1, "XYZ");
   tdrStyle->SetStripDecimals(kTRUE);
   tdrStyle->SetTickLength(0.03, "XYZ");
   tdrStyle->SetNdivisions(510, "XYZ");
   tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
   tdrStyle->SetPadTickY(1);

// Change for log plots:
   tdrStyle->SetOptLogx(0);
   tdrStyle->SetOptLogy(0);
   tdrStyle->SetOptLogz(0);

// Postscript options:
   tdrStyle->SetPaperSize(20.,20.);
   // tdrStyle->SetLineScalePS(Float_t scale = 3);
   // tdrStyle->SetLineStyleString(Int_t i, const char* text);
   // tdrStyle->SetHeaderPS(const char* header);
   // tdrStyle->SetTitlePS(const char* pstitle);

   // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
   // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
   // tdrStyle->SetPaintTextFormat(const char* format = "g");
   // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
   // tdrStyle->SetTimeOffset(Double_t toffset);
   // tdrStyle->SetHistMinimumZero(kTRUE);

   tdrStyle->cd();

}

class templateId {
public:
  int flav;
  int categ;
  templateId(int theFlav,int theCateg) : flav(theFlav), categ(theCateg) {}
  string name() {
    string sFlav[3] = {"Q", "C", "B"};
    string sLegend("bbb");
    sLegend.replace(categ,1,sFlav[flav]);
    return sLegend;
  }
};


class templateRef {
public:
  string name;
  TH2F* hist;
  double inival;
  int color;
  double sumW; // total content before renormalization
  int nNuisance;  // number of nuisance parameters
  vector<TH2F*> tempPlus; // vector of templates for +nsig of nuisance parameter
  vector<TH2F*> tempMinus; // vector of templates for -nsig of nuisance parameter
  vector<float> nSigma; // vector of nsig for nuisance parameters

  templateRef() {
    std::cout << "Warning: templateRef default constructor called" << std::endl;
  }
  templateRef(const char* theName,TH2F* theHist,const double theInival = 1,const int theColor=-1) : name(theName), hist(theHist), 
												    inival( theInival ), 
												    color( theColor), nNuisance(0) {
    std::cout << "templateRef constructor: theHist=" << theHist << " theColor=" << theColor << std::endl;
    if (theColor == -1) {
      color = theHist->GetLineColor();
    } else {
      theHist->SetLineColor(color);
    }
    sumW = theHist->GetSumOfWeights();
  }
  void normalize() {
    if (hist->GetSumOfWeights() >0) {
      float scaleFac = 1. / hist->GetSumOfWeights();
      hist->Scale( scaleFac );
      for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
	if (tempPlus[iNuisance]->GetSumOfWeights() >0) {
	  tempPlus[iNuisance]->Scale( 1. / tempPlus[iNuisance]->GetSumOfWeights() );
	} else {
	  std::cout << "templateRef::normalize: Error: Plus template non-normalizable" << std::endl;
	  return;
	}
	if (tempMinus[iNuisance]->GetSumOfWeights() >0) {
	  tempMinus[iNuisance]->Scale( 1. / tempMinus[iNuisance]->GetSumOfWeights() );
	} else {
	  std::cout << "templateRef::normalize: Error: Minus template non-normalizable" << std::endl;
	  return;
	}
      }
    } else {
      std::cout << "templateRef::normalize(): bad sum of weights " << hist->GetSumOfWeights() << std::endl;
    }
  }
  void addNuisance(TH2F* theTempPlus, TH2F* theTempMinus, float theNSigma) {
    tempPlus.push_back( theTempPlus );
    tempMinus.push_back( theTempMinus );
    nSigma.push_back( theNSigma );
    ++nNuisance;
  }
  // copy constructor
  templateRef(const templateRef& otr) : name(otr.name), hist(otr.hist), inival(otr.inival),
					color (otr.color), sumW(otr.sumW), 
					nNuisance(otr.nNuisance), 
					tempPlus(otr.tempPlus),
					tempMinus(otr.tempMinus),
					nSigma(otr.nSigma)
  
{
    std::cout << "Here is the templateRef copy constructor" << std::endl;
//     tempPlus = otr.tempPlus;
//     tempMinus = otr.tempMinus;
//     nSigma = otr.nSigma;
  }
  
  static templateRef* merge(const char* theName, const char* theHistName, 
			    templateRef* tr1, templateRef* tr2, const double theInival = 1,const int theColor=-1) {
    // create a new template as merge of two existing ones
    TH2F* theHist = new TH2F( *(tr1->hist) );
    theHist->SetName( theHistName );
    theHist->SetTitle( theHistName );
    theHist->Add( tr2->hist );
    templateRef* mergedTemplate = new templateRef(theName,theHist,theInival,theColor);
    
    // add the nuisance parameter histos
    if ( tr1->nNuisance != tr2->nNuisance ) {
      std::cout << "templateRef::merge: nNuisance mismatch: " << tr1->nNuisance
		<< " " << tr2->nNuisance << std::endl;
      return mergedTemplate;
    }
    for (int iNuisance=0; iNuisance<tr1->nNuisance; ++iNuisance) {
      if ( tr1->nSigma[iNuisance] != tr2->nSigma[iNuisance] ) {
	std::cout << "templateRef::merge: nSigma mismatch: " << tr1->nSigma[iNuisance]
		  << " " << tr2->nSigma[iNuisance] << std::endl;
	return mergedTemplate;
      }
      // sum plus
      TH2F* theTempPlus = new TH2F( *(tr1->tempPlus[iNuisance]) );
      theTempPlus->SetName( Form("%s_Nuis%d_%s",theHistName,iNuisance,"Plus") );
      theTempPlus->SetTitle( Form("%s_Nuis%d_%s",theHistName,iNuisance,"Plus") );
      theTempPlus->Add( tr2->tempPlus[iNuisance] );
      // sum minus
      TH2F* theTempMinus = new TH2F( *(tr1->tempMinus[iNuisance]) );
      theTempMinus->SetName( Form("%s_Nuis%d_%s",theHistName,iNuisance,"Minus") );
      theTempMinus->SetTitle( Form("%s_Nuis%d_%s",theHistName,iNuisance,"Minus") );
      theTempMinus->Add( tr2->tempMinus[iNuisance] );
      // add nuisance
      mergedTemplate->addNuisance( theTempPlus, theTempMinus, tr1->nSigma[iNuisance] );
    }
    return mergedTemplate;
  }
};

double uplimit(double x0,double sigma) {
  double CL = 0.95;
  // CL = 0.835;

  //std::cout << " Confidence level: " << CL << std::endl;
  //double dxOverSig = TMath::ErfInverse( 2*CL - 1 );
  //std::cout << " dxOverSig = " << dxOverSig << std::endl;

  //std::cout << "Upper limit  (non-constrained) : " << (x0 + sigma * dxOverSig) << std::endl;

  double dxOverSig = TMath::ErfInverse( CL + (1 - CL) * TMath::Erf( -x0/(sigma*sqrt(2)) ) );
  
  //std::cout << "Upper limit (zero constrained) : " << (x0 + sigma * dxOverSig) << std::endl;

  return (x0 + sqrt(2) * sigma * dxOverSig);
}


void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);  // declaration is needed in templateFitter

class templateFitter {
  
public:
  TH2F* data;
  std::vector<TH2F*> mc;
  std::vector<templateRef>* templates;
  std::vector<double> minConstraint;
  std::vector<double> maxConstraint;
  int fitLow;
  int fitHigh;
  double fitChi2;
  int fitNdf;
  int fitStat;
  int nNuisance;

  templateFitter(TH2F* theData,TObjArray* theMc,std::vector<templateRef>* theTemplates,int theNNuisance=0) {
    data = theData;
    TIter next(theMc);
    TH2F* theTemp;
    while ( (theTemp = (TH2F*) next()) ) {
      mc.push_back( theTemp );
      minConstraint.push_back(0);
      maxConstraint.push_back(0);
    }
    templates = theTemplates;
    fitLow = 1;  // TH convention, count from 1...nbins
    fitHigh = theData->GetNbinsX();
    nNuisance = theNNuisance;
    // no constraints for nuisance parameters
    for (unsigned int i=mc.size(); i<mc.size()+nNuisance; ++i) {
      minConstraint.push_back(0);
      maxConstraint.push_back(0);
    }
  }
  void SetRangeX(int theLow,int theHigh) {
    fitLow = theLow;
    fitHigh = theHigh;
  }
  void Constrain(unsigned int jpar,double lowVal,double highVal) {
    unsigned int ipar = jpar - 1;
    if (ipar<minConstraint.size() && ipar<maxConstraint.size()) {
      minConstraint[ipar] = lowVal;
      maxConstraint[ipar] = highVal;
    } else {
      std::cout << "Constrain: bad ipar=" << ipar << std::endl;
    }
  }
  TH2F* GetPlot() {
    return GetData();
  }
  TH2F* GetData() {
    return data;
  }
  TH2F* GetMc(unsigned int imc) {
    if (imc < mc.size()) {
      return mc[imc];
    } else {
      std::cout << "GetMc: Bad imc = " << imc << std::endl;
      return 0;
    }
  }    
  void GetResult(int iPar, double& theValue, double& theError){
    // parameter numbering from zero
    double value;
    double error;
    int getpar = gMinuit->GetParameter(iPar,value,error);
    theValue = value;
    theError = error;
  };
  double GetChisquare() {
    return fitChi2;
  };
  int GetNDF() {
    return fitNdf;
  };
  double GetProb() { 
    return TMath::Prob(fitChi2,fitNdf);
  };
  int Fit(){
    int nFitPar = mc.size() + nNuisance;
    TMinuit* gMinuit = new TMinuit(nFitPar);
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    arglist[0] = 1;
    Int_t ierflg = 0;
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
    Double_t* vstart = new Double_t[nFitPar];
    Double_t* vstep = new Double_t[nFitPar];
    // for normalization of the initial fractions
    double iniTot = 0;
    // initialize the template fraction parameters
    for (unsigned int i=0; i<mc.size(); ++i) {
      iniTot += (*templates)[i].inival;
    }
    for (unsigned int i=0; i<mc.size(); ++i) {
      vstart[i] = (*templates)[i].inival / iniTot;
      vstep[i] = vstart[i] / 100;
      //gMinuit->mnparm(i,"par",vstart[i],vstep[i],0,2*data->GetSumOfWeights() / mc[i]->GetSumOfWeights(),ierflg);
      gMinuit->mnparm(i,"par",vstart[i],vstep[i],minConstraint[i],maxConstraint[i],ierflg);
    }
    // initialize the nuisance parameters
    for (unsigned int i=mc.size(); i<mc.size()+nNuisance; ++i) {
      vstart[i] = 0;
      vstep[i] = 0.1;
      minConstraint[i] = 0;
      maxConstraint[i] = 0;
      gMinuit->mnparm(i,"par",vstart[i],vstep[i],minConstraint[i],maxConstraint[i],ierflg);
    }
      

//     // tell minuit to use gradient
//     arglist[0] = 1;
//     arglist[1] = 1;
//     gMinuit->mnexcm("SET GRA",arglist,2,ierflg);
#ifdef INIHESSE
    arglist[0] = 0;
    arglist[1] = 0;
    gMinuit->mnexcm("HESSE",arglist,0,ierflg);
#endif
    arglist[0] = 50000;
    arglist[1] = .1;
    gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
    // get statistics from minuit
    Double_t fmin; 
    Double_t fedm; 
    Double_t errdef; 
    Int_t npari; 
    Int_t nparx; 
    Int_t istat;
    gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
    fitChi2 = fmin;
    fitNdf = (gTemplateFitter->fitHigh - gTemplateFitter->fitLow + 1) * 6 - npari;
    fitStat = istat;
    std::cout << " fmin=" << fmin << " fedm=" << fedm << " errdef=" << errdef
	      << " npari=" << npari << " nparx=" << nparx << " istat=" << istat << std::endl;
    std::cout << " ndf=" << fitNdf << std::endl;

    return 0;
  };
};

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {    
  // this is the FCN defining the objective function
  double T = gTemplateFitter->GetData()->GetSumOfWeights();
  double chisq = 0;
  for (int binx=gTemplateFitter->fitLow; binx<=gTemplateFitter->fitHigh; ++binx) {
    for (int biny=1; biny<=6; ++biny) {
      double errData = gTemplateFitter->GetData()->GetBinError(binx,biny);
//       if (errData == 0) {
// 	std::cout << " errData=0 for binx, biny=" << binx << " " << biny << std::endl;
//       }
      if (errData == 0) errData = 1;
      double delta = gTemplateFitter->GetData()->GetBinContent(binx,biny);
      double errMcSq = 0;
//      for (int imc=0; imc<npar; ++imc) {
      for (int imc=0; imc<npar-gTemplateFitter->nNuisance; ++imc) {
	double templateVal = gTemplateFitter->GetMc(imc)->GetBinContent(binx,biny);
	// template shifts according to nuisance parameters
	for (int iNuisance=0; iNuisance<gTemplateFitter->nNuisance; ++iNuisance) {
	  double theNuisancePar = par[ gTemplateFitter->mc.size() + iNuisance ];
	  if (theNuisancePar>0) {
	    templateVal += ( (*gTemplateFitter->templates)[imc].tempPlus[iNuisance]->GetBinContent(binx,biny) - templateVal )
	      * theNuisancePar / (*gTemplateFitter->templates)[imc].nSigma[iNuisance];
	  } else {
	    templateVal -= ( (*gTemplateFitter->templates)[imc].tempMinus[iNuisance]->GetBinContent(binx,biny) - templateVal )
	      * theNuisancePar / (*gTemplateFitter->templates)[imc].nSigma[iNuisance];
	  }
	}

	delta -= par[imc] * T * templateVal;
	double errMCi = par[imc] * T * gTemplateFitter->GetMc(imc)->GetBinError(binx,biny);
	errMcSq += errMCi * errMCi;
      }
//       std::cout << "binx " << binx << " biny " << biny 
// 		<< " data " <<  gTemplateFitter->GetData()->GetBinContent(binx,biny)
// 		<< " edata " << errData 
// 		<< " mc0 " << par[0] * T *  gTemplateFitter->GetMc(0)->GetBinContent(binx,biny)
// 		<< " emc " << par[0] * T * gTemplateFitter->GetMc(0)->GetBinError(binx,biny)
// 		<< " delta " << delta << " edelta " << sqrt( errData * errData + errMcSq)
// 		<< " dchisq " << delta*delta / (errData * errData + errMcSq) << std::endl;
      chisq += delta*delta / (errData * errData + errMcSq);
    }
  }
  // add contributions from nuisance parameters to chi2
  for (int iNuisance=0; iNuisance<gTemplateFitter->nNuisance; ++iNuisance) {
    double theNuisancePar = par[ gTemplateFitter->mc.size() + iNuisance ];
    chisq += theNuisancePar * theNuisancePar;
  }

  f = chisq;

  if (iflag == 2) {
    // compute derivative
    for (int k=0; k<npar; ++k) {
      double deriK = 0;
      for (int binx=gTemplateFitter->fitLow; binx<=gTemplateFitter->fitHigh; ++binx) {
	for (int biny=1; biny<=6; ++biny) {
	  double delta = gTemplateFitter->GetData()->GetBinContent(binx,biny);
	  double errData = gTemplateFitter->GetData()->GetBinError(binx,biny);
	  if (errData == 0) errData = 1;
	  double errSq = errData*errData;
	  for (int imc=0; imc<npar; ++imc) {
	    delta -= par[imc] * T * gTemplateFitter->GetMc(imc)->GetBinContent(binx,biny);
	    double errMCi = par[imc] * T * gTemplateFitter->GetMc(imc)->GetBinError(binx,biny);
	    errSq += errMCi * errMCi;
	  }
	  deriK += (- 2 * T * gTemplateFitter->GetMc(k)->GetBinContent(binx,biny) * delta) / errSq;
	  deriK += (- 2 * par[k] * T * T * gTemplateFitter->GetMc(k)->GetBinError(binx,biny) * gTemplateFitter->GetMc(k)->GetBinError(binx,biny)
		    * delta * delta ) / (errSq * errSq);
	}
      }
      gin[k] = deriK;
    }
  }
}


TH2F* getTrigsAndMerge(TFile* fa,char* hname,const int nTCombData,string tCombData[]) {
  TH2F* theData[nTCombData];
  if (fa == NULL) {
    std::cout << "getTrigsAndMerge: no hfile opened" << std::endl;
    return NULL;
  }
  std::cout << "getTrigsAndMerge: hname=" << hname << std::endl;
  for (int iTCombData=0; iTCombData<nTCombData; ++iTCombData) {
    std::cout << "getTrigsAndMerge: iTCombData= " << iTCombData << " tCombData=" << tCombData[iTCombData] << std::endl;

    theData[iTCombData] =  (TH2F*) fa->Get(Form("%s%s",hname,tCombData[iTCombData].c_str()));
    if (theData[iTCombData] == NULL) {
      std::cout << "getTrigsAndMerge: histogram " << Form("%s%s",hname,tCombData[iTCombData].c_str()) << " not found" << std::endl;
      return 0;
    }
  }
  string theName( hname );
  // strip off any path in histogram name
  long unsigned int nsl = theName.find_last_of('/');
  if (nsl != theName.npos) {
    theName = theName.substr(nsl+1,theName.size()-nsl-1);
  }
  string theTitle( theData[0]->GetTitle() );
  TH2F* mergedData = new TH2F( *theData[0] );
  for (int iTCombData=0; iTCombData<nTCombData; ++iTCombData) {
    theName += tCombData[iTCombData];
    if (iTCombData != 0) theTitle += tCombData[iTCombData];
  }
  mergedData->SetName( theName.c_str() );
  for (int iTCombData=1; iTCombData<nTCombData; ++iTCombData) {
    mergedData->Add( theData[iTCombData] );
  }
  return mergedData;
}




void customFitNuisance() {
  // fit with corrected templates using the minuit fitter
  setTDRStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  canvas = new TCanvas ("cg1","PadX",10,10,800,600);
  gStyle->SetPadColor(0);
  canvas->SetFillColor(0);

  const int maxfa = 2;
  TFile* fa[maxfa];

  // this is for triple tag mass spectrum and background templates
  fa[0] = new TFile("../TripleBtagAnalysis.root");
  fa[1] = new TFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/CMSSW424p1/v0.2/mc/Fall11/SUSYBBHToBB_M-120_7TeV-pythia6-tauola/TripleBtagAnalysisV3/TripleBtagAnalysis.root");

  for (int ifa=0; ifa<maxfa; ++ifa) {
    if (fa[ifa] == NULL) {
      std::cout << "Problem opening file " << ifa << std::endl;
      return;
    }
  }


/// _SYST_BEGIN  read 'systematics' template from files

std::map<TString,TFile*>_tempPlusFiles;
std::map<TString,TFile*>_tempMinusFiles;
std::map<TString,TFile*>_tempCentralFiles;

std::map<TString,double>_tempSigmas;


int isFileOk = (int) gSystem->Exec(TString("[ -f syst_templates")+TString(" ]"));

if (isFileOk==0) {  
ifstream ff("syst_templates"); ///must exists
TString str;
str.ReadLine(ff);

while (str.Length()>0 )
{


TObjArray * arr = str.Tokenize(TString(':'));

if (arr->GetEntries()==4 && !str.Contains("#")) 
{

TString _id = arr->At(0)->GetName();
TString _type = arr->At(1)->GetName();
TString _sigma = arr->At(2)->GetName();
TString _fileName = arr->At(3)->GetName();

if (_type.Contains("up")) _tempPlusFiles[_id]=TFile::Open(_fileName);
if (_type.Contains("down")) _tempMinusFiles[_id]=TFile::Open(_fileName);
if (_type.Contains("central")) _tempCentralFiles[_id]=TFile::Open(_fileName);
if ( _tempSigmas.find(_id)==_tempSigmas.end() ) _tempSigmas[_id]=_sigma.Atof();


}

str.ReadLine(ff);

} ///while
} ///if isFileOk



/*
 // Example of iteration trough systematics

for (std::map<TString,double>::iterator it =  _tempSigmas.begin(); it !=   _tempSigmas.end(); it++)
{
	std::cout<<"The file for up shift of "<<it->first<<" is "<<_tempPlusFiles[it->first]->GetName()<<std::endl;	
	std::cout<<"The file for down shift of "<<it->first<<" is "<<_tempMinusFiles[it->first]->GetName()<<std::endl;	
}

*/


/// _SYST_END

  // now read the templates (taken from prottemplate2.C)
#define MC
#ifdef MC
  std::string tSel("AllTrigWeighted");
#else
  //std::string tSel("Trig2");
  std::string tSel("Trig0");
#endif

  // this is for combination of triggers
  const int nTCombData = 3;
  std::string tCombData[nTCombData] = {"Trig0", "Trig1", "Trig2"};
//   const int nTCombData = 1;
//   std::string tCombData[nTCombData] = {"Trig0"};

//  const int ncateg = 3;
//  const int ncorr=2;
//  const int nfc = 3;  
//  const int nbtag = 4;

  // mass range bins
  int lowBinX = 3;
  int highBinX = 22;
//   int lowBinX = 6;
//   int highBinX = 9;

  TH2F* mjjEbtTemplate[nfc][nbtag][ncateg];
  string tFlav[3] = {"Uds", "C", "B"};
  string sFlav[3] = {"Q", "C", "B"};

  string sfc[nfc] = { "q", "c", "b" };

  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
//   const std::string sbtag[nbtag] = { "TCHPT"};

  // prepare the drawing
  const int maxcol = 100;
  int tc[maxcol];
  tc[0] = kBlack;
  tc[1] = kBlue;
  tc[2] = kRed;
  tc[3] = kGreen;
  tc[4] = kMagenta;
  tc[5] = kCyan;
  tc[6] = kOrange;
  tc[7] = 38;  // pale blue
  tc[8] = 46;  // pale red
  tc[9] = 30;  // pale green
  for (int ic=10; ic<maxcol; ++ic) {
    tc[ic] = ic+10;
  }

  TH2F* mjjEbtdata[nbtag];
  TH1D* mjjEbtdataProX[nbtag];
  TH1D* mjjEbtdataProY[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtdata[ibtag] = getTrigsAndMerge(fa[0],Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()),nTCombData,tCombData);

    //mjjEbtdata[ibtag] = (TH2F*) fa[0]->Get(Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),tSel.c_str()));
    if (mjjEbtdata[ibtag] == NULL) {
      std::cout << "Histogram not found: " << Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),tSel.c_str()) << std::endl;
      return;
    }
    mjjEbtdataProX[ibtag] = mjjEbtdata[ibtag]->ProjectionX(Form("%sProX",mjjEbtdata[ibtag]->GetName()),0,-1,"e");
    mjjEbtdataProX[ibtag]->SetMarkerStyle(20);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetTitle("m(Jet1 Jet2) [GeV/c^{2}]");
    mjjEbtdataProX[ibtag]->GetYaxis()->SetTitle("N");
    mjjEbtdataProY[ibtag] = mjjEbtdata[ibtag]->ProjectionY(Form("%sProY",mjjEbtdata[ibtag]->GetName()),lowBinX,highBinX,"e");
    mjjEbtdataProY[ibtag]->SetMarkerStyle(20);
    mjjEbtdataProY[ibtag]->GetXaxis()->SetTitle("EvtBTag");
    mjjEbtdataProY[ibtag]->GetYaxis()->SetTitle("N");
  }

  //int itpat = 0;
  int itpat = 3;

  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	mjjEbtTemplate[ifc][ibtag][icateg] = getTrigsAndMerge(fa[0],Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d",
									 sfc[ifc].c_str(),
									 sbtag[ibtag].c_str(),icateg,itpat),nTCombData,tCombData);

// 	mjjEbtTemplate[ifc][ibtag][icateg] = (TH2F*) fa[0]->Get(Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d%s",
// 									    sfc[ifc].c_str(),
// 									    sbtag[ibtag].c_str(),icateg,itpat,tSel.c_str()));
	  
	if (mjjEbtTemplate[ifc][ibtag][icateg] == NULL) {
	  std::cout << "Histogram not found: " << Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d%s",
									    sfc[ifc].c_str(),
						       sbtag[ibtag].c_str(),icateg,itpat,tSel.c_str())
		    << std::endl;
	  return;
	}
#ifdef SCALE
	mjjEbtTemplate[ifc][ibtag][icateg][itpat]->Scale( mjjEbtTemplate[ifc][ibtag][icateg][itpat]->GetEntries() / mjjEbtTemplate[ifc][ibtag][icateg][itpat]->GetSumOfWeights() );
#endif
	mjjEbtTemplate[ifc][ibtag][icateg]->SetLineColor( tc[icateg+1+3*ifc] );
	mjjEbtTemplate[ifc][ibtag][icateg]->SetLineWidth(3);  
	
      }
    }
  }



///_SYST_BEGIN_
///Read systematics for Background

std::vector<mjjEbtTemplateType> _plusTemplatesBkg;
std::vector<mjjEbtTemplateType> _minusTemplatesBkg;
std::vector<double> _sigmaTemplatesBkg;


for (std::map<TString,double>::iterator it =  _tempSigmas.begin(); it !=   _tempSigmas.end(); it++)
{
        std::cout<<"The file for up shift of "<<it->first<<" is "<<_tempPlusFiles[it->first]->GetName()<<std::endl;
        std::cout<<"The file for down shift of "<<it->first<<" is "<<_tempMinusFiles[it->first]->GetName()<<std::endl;


if (it->first.Contains("_Sgn")) continue; ///skip signal shifted templates
        mjjEbtTemplateType _plusTemplate;
        mjjEbtTemplateType _minusTemplate;

  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {

        _plusTemplate._array[ifc][ibtag][icateg] = getTrigsAndMerge(_tempPlusFiles[it->first],Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d",
                                                                         sfc[ifc].c_str(),
                                                                         sbtag[ibtag].c_str(),icateg,itpat),nTCombData,tCombData);

///Change a little bit name, is it necesseary?
 _plusTemplate._array[ifc][ibtag][icateg]->SetName( TString(_plusTemplate._array[ifc][ibtag][icateg]->GetName())+"_"+it->first+"_plus");


     _minusTemplate._array[ifc][ibtag][icateg] = getTrigsAndMerge(_tempMinusFiles[it->first],Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d",
                                                                         sfc[ifc].c_str(),
                                                                         sbtag[ibtag].c_str(),icateg,itpat),nTCombData,tCombData);

///Change a little bit name, is it necesseary?
 _minusTemplate._array[ifc][ibtag][icateg]->SetName( TString(_minusTemplate._array[ifc][ibtag][icateg]->GetName())+"_"+it->first+"_minus");


        if (_minusTemplate._array[ifc][ibtag][icateg] == NULL || _plusTemplate._array[ifc][ibtag][icateg] == NULL) {
          std::cout << "Plus/Minus Histogram not found: " << Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d%s",
                                                                            sfc[ifc].c_str(),
                                                       sbtag[ibtag].c_str(),icateg,itpat,tSel.c_str())
                    << std::endl;
          return;


        } ///if

} ///icateg
} ///ibtag
} ///ifc

_sigmaTemplatesBkg.push_back(it->second);
_plusTemplatesBkg.push_back(_plusTemplate);
_minusTemplatesBkg.push_back(_minusTemplate);

} ///for  _tempSigmas


///_SYST_END_




  // now read the MC signal template
  TH2F* mjjEbtsignal[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtsignal[ibtag] = (TH2F*) fa[1]->Get(Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0"));
    if ( mjjEbtsignal[ibtag] == NULL ) {
      std::cout << "Signal histogram not found: " << Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0") << std::endl;
    }
    else {
      // change name to avoid duplication
      mjjEbtsignal[ibtag]->SetName( Form("%s_Signal",mjjEbtsignal[ibtag]->GetName()) );
    }
  }  


///_SYST_BEGIN
///Read systematics for Signal


std::vector<mjjEbtTemplateTypeSgn> _plusTemplatesSgn;
std::vector<mjjEbtTemplateTypeSgn> _minusTemplatesSgn;
std::vector<double >_sigmaTemplatesSgn;

for (std::map<TString,double>::iterator it =  _tempSigmas.begin(); it !=   _tempSigmas.end(); it++)
{
        std::cout<<"The file for up shift of "<<it->first<<" is "<<_tempPlusFiles[it->first]->GetName()<<std::endl;
        std::cout<<"The file for down shift of "<<it->first<<" is "<<_tempMinusFiles[it->first]->GetName()<<std::endl;


if (it->first.Contains("_Bkg")) continue; ///skip background shifted templates

  mjjEbtTemplateTypeSgn _plusTemplateSgn;
  mjjEbtTemplateTypeSgn _minusTemplateSgn;


 for (int ibtag=0; ibtag<nbtag; ++ibtag) {

        _plusTemplateSgn._array[ibtag]=(TH2F*)_tempPlusFiles[it->first]->Get(Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0"));
        _minusTemplateSgn._array[ibtag]=(TH2F*)_tempMinusFiles[it->first]->Get(Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0"));

	if ( _plusTemplateSgn._array[ibtag] == NULL ||  _minusTemplateSgn._array[ibtag]==NULL ) {
      		std::cout << "Syst Signal histogram not found: " << Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0") << std::endl;
	    }
	    else {
      // change name to avoid duplication
      _plusTemplateSgn._array[ibtag]->SetName( Form("%s_Signal_%s_plus", _plusTemplateSgn._array[ibtag]->GetName(),it->first.Data()) );
      _minusTemplateSgn._array[ibtag]->SetName( Form("%s_Signal_%s_minus", _minusTemplateSgn._array[ibtag]->GetName(),it->first.Data()) );
    	} /// if 

}   ///ibtag


_sigmaTemplatesSgn.push_back(it->second);
_plusTemplatesSgn.push_back(_plusTemplateSgn);
_minusTemplatesSgn.push_back(_minusTemplateSgn);

} /// _tempSigmas

///_SYST_END_



  // prepare legend
  TLatex legend;
  legend.SetNDC();
  float xLegend = 0.55;
  float yLegend0 = 0.85;
  float yIncrement = 0.055;
  legend.SetTextAlign(12);
  legend.SetTextSize(0.035);

  // build the set of MC histograms
  int icorr = 1;

  // build the nine basic templateRef's
  templateRef* mjjEbtTemplateRef[nfc][nbtag][ncateg];
  for (int ibtag = 0; ibtag < nbtag; ++ibtag) {
    mjjEbtTemplateRef[0][ibtag][0] = new templateRef("Qbb",mjjEbtTemplate[0][ibtag][0], 0.1);
    mjjEbtTemplateRef[1][ibtag][0] = new templateRef("Cbb",mjjEbtTemplate[1][ibtag][0], 0.2);
    mjjEbtTemplateRef[2][ibtag][0] = new templateRef("Bbb",mjjEbtTemplate[2][ibtag][0], 1);
    mjjEbtTemplateRef[0][ibtag][1] = new templateRef("bQb",mjjEbtTemplate[0][ibtag][1], 0.1);
    mjjEbtTemplateRef[1][ibtag][1] = new templateRef("bCb",mjjEbtTemplate[1][ibtag][1], 0.2);
    mjjEbtTemplateRef[2][ibtag][1] = new templateRef("bBb",mjjEbtTemplate[2][ibtag][1], 1);
    mjjEbtTemplateRef[0][ibtag][2] = new templateRef("bbQ",mjjEbtTemplate[0][ibtag][2], 0.1);
    mjjEbtTemplateRef[1][ibtag][2] = new templateRef("bbC",mjjEbtTemplate[1][ibtag][2], 0.2);
    mjjEbtTemplateRef[2][ibtag][2] = new templateRef("bbB",mjjEbtTemplate[2][ibtag][2], 0.3);
  }




  //
  // Add nuisance parameter information to templates HERE
  //

///_SYST_BEGIN
 
 for (int sys=0;sys<(int) _plusTemplatesBkg.size();sys++){
	  for (int ifc=0; ifc<nfc; ++ifc) {
		for (int ibtag=0; ibtag<nbtag; ++ibtag) {
		      for (int icateg=0; icateg<ncateg; ++icateg) {
			mjjEbtTemplateRef[ifc][ibtag][icateg]->addNuisance(_plusTemplatesBkg[sys]._array[ifc][ibtag][icateg],
							_minusTemplatesBkg[sys]._array[ifc][ibtag][icateg],
			 				 _sigmaTemplatesBkg[sys]);
			} ///icateg
		}///ibtag
	}///ifc
 } ///sys


///_SYST_END


  for (int ibtag = 0; ibtag < nbtag; ++ibtag) {

  // produce templates merging the first two categories
  templateRef* mjjEbtTemplateMerge01[nfc][nbtag];
  for (int ifc=0; ifc<nfc; ++ifc) {
    // templates to be merged
    mjjEbtTemplateMerge01[ifc][ibtag] = 
      templateRef::merge( Form("(%sb)b",sFlav[ifc].c_str()),
			  Form("MassBTagTemplate_%s_%s_Cat%sTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),
			       "01",itpat),
			  mjjEbtTemplateRef[ifc][ibtag][0],mjjEbtTemplateRef[ifc][ibtag][1],0.3);
  }

  



  // set up the templates configuration for this fit
  vector<templateRef> templates;
  templates.push_back( *mjjEbtTemplateRef[2][ibtag][0] ); // Bbb
  templates.push_back( *mjjEbtTemplateMerge01[1][ibtag] );  // (Cb)b
  templates.push_back( *mjjEbtTemplateMerge01[0][ibtag] );  // (Qb)b
  templates.push_back( *mjjEbtTemplateRef[2][ibtag][2] ); // bbB
  templates.push_back( *mjjEbtTemplateRef[0][ibtag][2] ); // bbQ
  

  // add the signal template
#define SIGNAL
#ifdef SIGNAL
  if (mjjEbtsignal[ibtag] != NULL) {
    templates.push_back(  templateRef("H",mjjEbtsignal[ibtag],0.001,2) );

///_SYST_BEGIN
 for (int sys=0;sys<(int) _plusTemplatesSgn.size();sys++)
 templates.back().addNuisance(_plusTemplatesSgn[sys]._array[ibtag],_minusTemplatesSgn[sys]._array[ibtag],_sigmaTemplatesSgn[sys]);
///_SYST_END


  } else {
    std::cout << "cannot add signal template for " << sbtag[ibtag] << std::endl;
  }


#endif




  bool constrainSignal = false;   // set to true if signal fraction must be non-negative
                                  // only works if signal template is labeled H

  // normalize the templates
  for (vector<templateRef>::iterator trit=templates.begin(); trit != templates.end(); ++trit) {
    (*trit).normalize();
  }

  TObjArray *mc = new TObjArray();
  for (vector<templateRef>::iterator tit=templates.begin(); tit != templates.end(); ++tit) {
    mc->Add( (*tit).hist );
  }
  std::cout << "Number of templates = " << mc->GetEntries() << std::endl;

  gTemplateFitter = new templateFitter(mjjEbtdata[ibtag], mc, &templates,  _plusTemplatesBkg.size() );

  for (int itemp=0; itemp<mc->GetEntries(); ++itemp) {
    if ((*gTemplateFitter->templates)[itemp].name != "H" || constrainSignal) {
      gTemplateFitter->Constrain(itemp+1,0.0,2*gTemplateFitter->data->GetSumOfWeights() / gTemplateFitter->mc[itemp]->GetSumOfWeights());
    }
  }
  gTemplateFitter->SetRangeX(lowBinX,highBinX);
  Int_t status = gTemplateFitter->Fit(); 
  cout << sbtag[ibtag] << " fit status: " << status << endl;
  TH2F* result = (TH2F*) gTemplateFitter->GetPlot();
  mjjEbtdata[ibtag]->SetMarkerStyle(20);
  mjjEbtdata[ibtag]->GetXaxis()->SetTitle("m(Jet1 Jet2) [GeV/c^{2}]");
  mjjEbtdata[ibtag]->GetYaxis()->SetTitle("EvtBTag");
  mjjEbtdata[ibtag]->GetZaxis()->SetTitle("N / 10 GeV/c^{2}");
  mjjEbtdata[ibtag]->ProjectionX()->Draw("Ep");
  //result->Draw("same");
  double value;
  double error;
  double fPenalty = -1;

  // loop over templates
  int icol = 1;
  int ii=0;
  TIter next(mc);
  TH2F* theTemp;
  THStack hsX("hsX","Stack of fitted templates (X projection)");
  THStack hsY("hsY","Stack of fitted templates (Y projection)");
  // find out how many BG components are inactive
  int nInactive = 0;
  string model = "";
  int iH = -1;  // pointer to Higgs template

  // vector to host projections
  std::vector<TH1D*> vProX;
  std::vector<TH1D*> vProY;

  while ( (theTemp = (TH2F*) next()) ) {
    gTemplateFitter->GetResult(ii, value, error);
    std::cout << "Template " << std::setw(8) << templates[ii].name << " fraction=" 
	      << std::setw(12) << std::setprecision(5) << value << " +- " 
	      << std::setw(12) << std::setprecision(5) << error << std::endl;
    //theTemp->Scale( value * mjjEbtdata[ibtag]->GetSumOfWeights() / theTemp->GetSumOfWeights() );
    theTemp->Scale( value * mjjEbtdata[ibtag]->GetSumOfWeights());
    theTemp->GetXaxis()->SetRange( gTemplateFitter->fitLow,gTemplateFitter->fitHigh );
    //theTemp->SetLineColor(tc[icol]);
    TH1D* theTempProX = theTemp->ProjectionX(Form("%sProX",theTemp->GetName()),0,-1,"e");
    theTempProX->SetLineColor( theTemp->GetLineColor() );
    theTempProX->SetFillColor( theTemp->GetLineColor() );
    theTempProX->SetMarkerStyle(20);
    vProX.push_back( theTempProX );
    std::cout << "Created projection of " << theTemp->GetName() << " as " << theTempProX->GetName() << std::endl;
    hsX.Add(theTempProX);
    //theTempProX->Draw("hist,same");

    TH1D* theTempProY = theTemp->ProjectionY(Form("%sProY",theTemp->GetName()),lowBinX,highBinX,"e");
    theTempProY->SetLineColor( theTemp->GetLineColor() );
    theTempProY->SetFillColor( theTemp->GetLineColor() );
    theTempProY->SetMarkerStyle(20);
    vProY.push_back( theTempProY );
    hsY.Add(theTempProY);

    string sLegend(Form("%6s %8.3f +- %8.3f",templates[ii].name.c_str(),value,error));
    legend.SetTextColor(templates[ii].color);
    legend.DrawLatex(xLegend+0.1,yLegend0-0.2-ii*yIncrement,sLegend.c_str());

    if (templates[ii].name == "H") {
      iH = ii;
    } else {
      if (value <1e-3) ++nInactive;
      if (model == "") {
	model += templates[ii].name;
      } else {
	model += ("+"+templates[ii].name);
      }
    }
    ++icol;
    ++ii;
  }

  // draw the projections
  int iT = 0;
  for (std::vector<TH1D*>::iterator vpit=vProX.begin(); vpit<vProX.end(); ++vpit) {
    (*vpit)->Draw();
    canvas->Print(Form("templateProX_%s_%d.png",sbtag[ibtag].c_str(),iT));
    ++iT;
  }

  // print Higgs statistics (if H template was fitted)
  if (iH>=0) {
    gTemplateFitter->GetResult(iH, value, error);
    // determine the fitting penalty
    fPenalty = error * mjjEbtdata[ibtag]->GetSumOfWeights() / sqrt( mjjEbtdata[ibtag]->GetSumOfWeights() );
    std::cout << "deltaS / sqrt(B) = " << error * mjjEbtdata[ibtag]->GetSumOfWeights()
	      << "  /  " << sqrt( mjjEbtdata[ibtag]->GetSumOfWeights() ) << "  =  "
	      << fPenalty << std::endl;
    // expected signal if MC Xsection is correct
    double HExpected = templates[iH].sumW * 0.501 * 111./38.;
    double HSeen = value * mjjEbtdata[ibtag]->GetSumOfWeights();
    double HSeenError = error * mjjEbtdata[ibtag]->GetSumOfWeights();
    std::cout << "Higgs scaling factor: " << HSeen / HExpected << " +- " << HSeenError / HExpected << std::endl;
    double HScal95 = uplimit( HSeen / HExpected, HSeenError / HExpected );
    std::cout << "Higgs scaling factor UL: " << HScal95 << std::endl;
    std::cout << "tan-beta             UL: " << 20 * sqrt( HScal95 ) << std::endl; 
    double HMCGen = 1051000;
    double HMCAccepted = mjjEbtsignal[ibtag]->GetEntries();
    double HEfficiency = HMCAccepted / HMCGen;
    double HCorr = HSeen / HEfficiency;
    double HCorrError = HSeenError / HEfficiency;
    double intLumi = 2.8; // in fb-1
    double HXSect = 0.001 * HCorr / intLumi;
    double HXSectError = 0.001 * HCorrError / intLumi;
    std::cout << "Efficiency = " << HEfficiency << std::endl
	      << "Higgs observed = " << HSeen << " +- " << HSeenError << std::endl
	      << "Higgs corrected = " << HCorr << " +- " << HCorrError << std::endl
	      << "Higgs XSection = " << HXSect << " +- " << HXSectError << " pb-1 " << std::endl
	      << "Higgs XSection < " << uplimit( HXSect, HXSectError) << " pb-1 " << std::endl
	      << "Assumed IntLumi = " << intLumi << "fb-1" << std::endl;

    // print the summary line
    std::cout << Form("| %6s | %28s+H | %10.2f | %d | %8.2f | %8.2f |",
		      sbtag[ibtag].c_str(),model.c_str(),gTemplateFitter->GetChisquare(),nInactive,fPenalty,HXSectError) << std::endl;

  }

  // print also the chi2
  legend.SetTextColor(1);
  legend.DrawLatex(xLegend,yLegend0,Form("#chi^{2}/N_{DF} = %6.2f / %d",gTemplateFitter->GetChisquare(),
						 gTemplateFitter->GetNDF()));
  legend.DrawLatex(xLegend,yLegend0-yIncrement,Form("Prob  = %6.4f ",gTemplateFitter->GetProb()));
  legend.DrawLatex(xLegend,yLegend0-2*yIncrement,sbtag[ibtag].c_str());

  //canvas->Print(Form("templateFit_%s.png",sbtag[ibtag].c_str()));

  std::cout << "ChiSquare   = " << gTemplateFitter->GetChisquare() << std::endl;
  std::cout << "NDF         = " << gTemplateFitter->GetNDF() << std::endl;
  std::cout << "Probability = " << gTemplateFitter->GetProb() << std::endl;

  // draw X stacked
  //mjjEbtdata[ibtag]->ProjectionX()->Draw("Ep");
  mjjEbtdataProX[ibtag]->Draw("Ep");
  std::cout << mjjEbtdataProX[ibtag]->GetName() << " Number of bins " << mjjEbtdataProX[ibtag]->GetXaxis()->GetNbins() << std::endl;
  //result->Draw("same");
  next.Reset();
  while ( (theTemp = (TH2F*) next()) ) {
    //theTemp->SetFillColor(kGreen);
    theTemp->SetFillColor( theTemp->GetLineColor() );
    theTemp->SetFillStyle(1001);
    std::cout << theTemp->GetName() << " Number of bins " << theTemp->GetXaxis()->GetNbins() << std::endl;
  }
  std::cout << "Draw stack hsX" << std::endl;
  hsX.Draw("HIST,SAME");
  //mjjEbtdata[ibtag]->ProjectionX()->Draw("EpSAME");
  std::cout << "Draw " << mjjEbtdataProX[ibtag]->GetName() << " Number of bins " << mjjEbtdataProX[ibtag]->GetXaxis()->GetNbins() << std::endl;
  mjjEbtdataProX[ibtag]->Draw("EpSAME");
  //result->Draw("same");
  for (int iTemplate=0; iTemplate<templates.size(); ++iTemplate) {
    gTemplateFitter->GetResult(iTemplate, value, error);
    string sLegend(Form("%6s %8.3f +- %8.3f",templates[iTemplate].name.c_str(),value,error));
    legend.SetTextColor(templates[iTemplate].color);
    legend.DrawLatex(xLegend+0.1,yLegend0-0.2-iTemplate*yIncrement,sLegend.c_str());
  }    
  // print also the chi2
  legend.SetTextColor(1);
  legend.DrawLatex(xLegend,yLegend0,Form("#chi^{2}/N_{DF} = %6.2f / %d",gTemplateFitter->GetChisquare(),
						gTemplateFitter->GetNDF()));
  legend.DrawLatex(xLegend,yLegend0-yIncrement,Form("Prob  = %6.4f ",gTemplateFitter->GetProb()));
  legend.DrawLatex(xLegend,yLegend0-2*yIncrement,sbtag[ibtag].c_str());

  canvas->Print(Form("templateFitStacked_%s_X.png",sbtag[ibtag].c_str()));
  std::cout << "Done drawing X" << std::endl;

  // draw Y stacked
  mjjEbtdataProY[ibtag]->Draw("Ep");
  next.Reset();
  std::cout << "Draw stack hsY" << std::endl;
  hsY.Draw("HIST,SAME");
  
  std::cout << "Draw " << mjjEbtdataProY[ibtag]->GetName() << std::endl;
  mjjEbtdataProY[ibtag]->Draw("EpSAME");
  for (int iTemplate=0; iTemplate<templates.size(); ++iTemplate) {
    gTemplateFitter->GetResult(iTemplate, value, error);
    string sLegend(Form("%6s %8.3f +- %8.3f",templates[iTemplate].name.c_str(),value,error));
    legend.SetTextColor(templates[iTemplate].color);
    legend.DrawLatex(xLegend+0.1,yLegend0-0.2-iTemplate*yIncrement,sLegend.c_str());
  }    
  // print also the chi2
  legend.SetTextColor(1);
  legend.DrawLatex(xLegend,yLegend0,Form("#chi^{2}/N_{DF} = %6.2f / %d",gTemplateFitter->GetChisquare(),
						gTemplateFitter->GetNDF()));
  legend.DrawLatex(xLegend,yLegend0-yIncrement,Form("Prob  = %6.4f ",gTemplateFitter->GetProb()));
  legend.DrawLatex(xLegend,yLegend0-2*yIncrement,sbtag[ibtag].c_str());

  canvas->Print(Form("templateFitStacked_%s_Y.png",sbtag[ibtag].c_str()));

  std::cout << "Done drawing Y" << std::endl;
  std::cout << "After printing YStack" << std::endl;

  // in case we fitted with signal, scale up
  next.Reset();
  while ( (theTemp = (TH2F*) next()) ) {
    for (int iTemplate=0; iTemplate<templates.size(); ++iTemplate) {
      gTemplateFitter->GetResult(iTemplate, value, error);
      if (templates[iTemplate].name == "H") {
	std::cout << "Higgs fit values : " << Form("%6s %8.3f +- %8.3f",templates[iTemplate].name.c_str(),value,error) << std::endl;
	double T = gTemplateFitter->GetData()->GetSumOfWeights();
	double nHiggsFitted = T * value;
	double enHiggsFitted = T * error;
	std::cout << "Higgs fitted number of events : " << Form("%6s %8.3f +- %8.3f",templates[iTemplate].name.c_str(),nHiggsFitted,enHiggsFitted) << std::endl;
	std::cout << "Higgs MC seen scaled to 1 fb-1 " << std::endl; 
      }
    }
  }

  // now draw X slices
  for (int biny=1; biny<=6; ++biny) {
    THStack hsSliX("hsSliX","Stack of fitted templates (Slices in X)");
    next.Reset();
    TH1D* mjjEbtdataSliX = mjjEbtdata[ibtag]->ProjectionX(Form("%sSliX%d",mjjEbtdata[ibtag]->GetName(),biny),biny,biny,"e");
    
    while ( (theTemp = (TH2F*) next()) ) {
      TH1D* mjjEbtTemplateSliX = theTemp->ProjectionX(Form("%sSliX%d",theTemp->GetName(),biny),biny,biny,"e");
      mjjEbtTemplateSliX->SetFillColor( theTemp->GetLineColor() );
      mjjEbtTemplateSliX->SetLineColor( theTemp->GetLineColor() );
      hsSliX.Add( mjjEbtTemplateSliX );
    }
    mjjEbtdataSliX->Draw("Ep");
    hsSliX.Draw("HIST,SAME");
    mjjEbtdataSliX->Draw("Ep,SAME");
    canvas->Print(Form("templateFitSliceX_%s_%d.png",sbtag[ibtag].c_str(),biny));
  }
}
}
