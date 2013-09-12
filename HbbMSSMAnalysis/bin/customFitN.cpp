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

TH2F* mergeSignal(TFile* fa,char* hname,const char* newname) {
  TH2F* theData[2];
  if (fa == NULL) {
    std::cout << "mergeSignal: no hfile opened" << std::endl;
    return NULL;
  }

//   // weights
//   double tWeight[2];
//   tWeight[0] = (524.90+265.75)*.94*0.94 + 193.58+251.12+453.33;
//   tWeight[1] = 246.53+732.73;
//   double sumWeight = tWeight[0]+tWeight[1];
//   for (int ii=0; ii<2; ++ii) {
//     tWeight[ii] = tWeight[ii] / sumWeight;
//   }
//   std::cout << "Trig weights: " << tWeight[0] << " " << tWeight[1] << std::endl;
  



  // only Trig0 is needed
  for (int ii=0; ii<1; ++ii) {
    std::cout << "mergeSignal: pick histogram " << Form("%sTrig%d",hname,ii) << std::endl;
    theData[ii] = (TH2F*) fa->Get(Form("%sTrig%d",hname,ii));
    if (theData[ii] == NULL) {
      std::cout << "mergeSignal: histogram " << Form("%sTrig%d",hname,ii) << " not found" << std::endl;
      return 0;
    }
    //theData[ii]->Draw();
  }
  std::string theName( hname );
  // strip off any path in histogram name
  long unsigned int nsl = theName.find_last_of('/');
  if (nsl != theName.npos) {
    theName = theName.substr(nsl+1,theName.size()-nsl-1);
  }
  std::string theTitle( theData[0]->GetTitle() );
  TH2F* mergedData = new TH2F( *theData[0] );

  //mergedData->Add(theData[0],theData[1],tWeight[0],tWeight[1]);

  if (newname == 0) {
    mergedData->SetName( Form("%s%s",theName.c_str(),"TrigMerged") );
  } else {
    mergedData->SetName( newname );
  }
  mergedData->SetTitle( mergedData->GetName() );
  std::cout << "mergeSignal: hist has name " << mergedData->GetName() << std::endl;
  //mergedData->Draw();
  //canvas->Print(Form("%s.png",newname));
  return mergedData;
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

    std::cout << "templateRef merge " << tr1->name << " " 
	      << tr1->hist->GetName() << " and " << tr2->name
	      << " " << tr2->hist->GetName() 
	      << " to " << theName << std::endl;
    
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

	delta -= (par[imc] * T * templateVal);
	double errMCi = par[imc] * T * gTemplateFitter->GetMc(imc)->GetBinError(binx,biny);
	errMcSq += (errMCi * errMCi);
      }
//       std::cout << "binx " << binx << " biny " << biny 
// 		<< " data " <<  gTemplateFitter->GetData()->GetBinContent(binx,biny)
// 		<< " edata " << errData 
// 		<< " mc0 " << par[0] * T *  gTemplateFitter->GetMc(0)->GetBinContent(binx,biny)
// 		<< " emc " << par[0] * T * gTemplateFitter->GetMc(0)->GetBinError(binx,biny)
// 		<< " delta " << delta << " edelta " << sqrt( errData * errData + errMcSq)
// 		<< " dchisq " << delta*delta / (errData * errData + errMcSq) << std::endl;
      chisq += (delta*delta / (errData * errData + errMcSq));
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
  std::cout << "getTrigsAndMerge: new name: " << mergedData->GetName() << std::endl;
  return mergedData;
}


void dump2D( TH2F* h ) {
  std::cout << " ===== dump2D for " << h->GetName() << std::endl;
  for (int ix=1; ix<h->GetXaxis()->GetNbins(); ++ix) {
    for (int iy=1; iy<h->GetYaxis()->GetNbins(); ++iy) {
       std::cout << ix << " " << iy << " "
		 << h->GetBinContent(ix,iy) << "  " << h->GetBinError(ix,iy) << std::endl;
     }
  }
}



void customFitN_Mass(int iSignal) {
  // fit with corrected templates using the minuit fitter
  setTDRStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  canvas = new TCanvas ("cg1","PadX",10,10,800,600);
  gStyle->SetPadColor(0);
  canvas->SetFillColor(0);

  const int ncateg = 3;
  const int ncorr=2;
  const int nfc = 3;
  const int nbtag = 1;
  // const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  const std::string sbtag[nbtag] = { "CSVT"};

  string sfc[nfc] = { "q", "c", "b" };

  std::string L1L2Mode("Weight");
  std::string signalMode("PU_WEIGHTED-NEW");
  //   std::string L1L2Mode("Cut");
  //   std::string signalMode("CUT_BASED");
  std::string scenario("LowMass2011");

  const int nSignal=7;
  int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };

  bool fitSignal = true;
  bool useTemplateError = false;

  bool activeNuisance[3] = { false, false, false };

  double efficiency[nSignal][nbtag];

  double intLumi = 0;
  string IgorScen("");
  string spacer("");
  string SashaPath("");

  if (scenario == "LowMass2011") {
    //intLumi = 2.66794; // in fb-1
    intLumi = 2.692643;  // with new method
    IgorScen.assign("low");
    spacer.assign("");
    SashaPath.assign("Data-Run2011AB");
  } else if (scenario == "MediumMass2011") {
    //intLumi = 3.99983; // in fb-1
    intLumi = 4.040802;
    IgorScen.assign("medium");
    spacer.assign("/MEDIUM");
    SashaPath.assign("Data-Run2011AB-Medium");
  } else {
    std::cout << "Bad scenario" << std::endl;
    return;
  }

  // summary variables across the masses
  double fPenalty[nbtag][nSignal];
  double HFraction[nbtag][nSignal];
  double HFractionError[nbtag][nSignal];
  double HSeen[nbtag][nSignal];
  double HSeenError[nbtag][nSignal];
  double nominalXSect[nbtag][nSignal];
  double HEfficiency[nbtag][nSignal];
  double HCorr[nbtag][nSignal];
  double HCorrError[nbtag][nSignal];
  double HXSect[nbtag][nSignal];
  double HXSectError[nbtag][nSignal];
  double HXSectUpLimit[nbtag][nSignal];
  double tanBetaUpLimit[nbtag][nSignal];

  string signalHistPattern("massEvBtag/mjjEvBTag_%s");
  if (L1L2Mode == "Weight") {
    signalHistPattern.assign("massEvBtagTW/mjjEvBTagTW_%s");
  }

  // arrays to hold the signal templates
  const int nSyst = 3;
  const int nUpDown = 2;
  std::string signalFile( Form("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-3/%s/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola%s/SF/job_1/TripleBtagAnalysisM-%d_%s.root",signalMode.c_str(),signalMass[iSignal],spacer.c_str(),signalMass[iSignal],IgorScen.c_str() ) );
  std::string signalSystFiles[nSyst][nUpDown];
  std::string systName[nSyst] = { "JES", "SFbc", "SFudsg" };
  std::string upDownName[nUpDown] = { "Up", "Down" };
  TH2F* hSignalSyst[nSyst][nUpDown][nbtag];

  string tFlav[3] = {"Uds", "C", "B"};
  string sFlav[3] = {"Q", "C", "B"};

  // this is for the combination of triggers
  const int nTCombData = 4;
  std::string tCombData[nTCombData] = {"Trig0", "Trig1", "Trig2", "Trig3"};



  // output file
  TFile* hout = new TFile(Form("customFitN-M-%d.root",signalMass[iSignal]),"recreate");
  hout->cd();
  TH2::AddDirectory(true);

  TFile* fSig = new TFile( signalFile.c_str() );
  if ( fSig == NULL ) {
    std::cout << "Could not open signal central file " << signalFile.c_str() << std::endl;
    return;
  } else {
    std::cout << "Open signal file " << signalFile.c_str() << std::endl;
  }

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

  TH2F* hSignalCentral[nbtag];

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hSignalCentral[ibtag] = mergeSignal(fSig,Form(signalHistPattern.c_str(),sbtag[ibtag].c_str()),
				       "bbH");
    if (hSignalCentral[ibtag] == NULL) {
      std::cout << "Reading of SignalCentral failed" << std::endl;
      return;
    }
    // rename
    hSignalCentral[ibtag]->SetName( Form("%s_M-%d",
					 hSignalCentral[ibtag]->GetName(),signalMass[iSignal]) );
    std::cout << "The signal hist has name " << hSignalCentral[ibtag]->GetName() << std::endl;
					     
    // read the efficiency
    TH1F* histEffMerged = (TH1F*) fSig->Get(Form("TrigEff/EffMerged%s",sbtag[ibtag].c_str()));
    if ( histEffMerged == NULL) {
      std::cout << "Efficiency histo not found" << std::endl;
      return;
    }
    double newEff = histEffMerged->GetBinContent(1);
    std::cout << "Mass= " << signalMass[iSignal] 
	      << " btag= " << sbtag[ibtag]
	      << " Efficiency = " << newEff << std::endl;
    efficiency[iSignal][ibtag] = newEff;

    // read the nominal cross section
    TH1F* histXSect = (TH1F*) fSig->Get("xsection/xsect");
    if ( histXSect == NULL) {
      std::cout << "xsection/xsect" << " not found" << std::endl;
      return;
    }
    nominalXSect[ibtag][iSignal] = histXSect->GetBinContent(1);


//     double normShould = 1000 * intLumi * efficiency[iSignal][ibtag];
//     double normIs = hSignalCentral[ibtag]->GetSumOfWeights();
    std::cout << hSignalCentral[ibtag]->GetName() << " TotalContents=" << hSignalCentral[ibtag]->GetSumOfWeights()
	      << std::endl;
//     fScal[iSignal][ibtag] = normShould / normIs;
//     std::cout << "normShould = " << normShould << " normIs " << normIs
// 	      << " rescale by " << fScal[iSignal][ibtag] << std::endl;
//     hSignalCentral[ibtag]->Scale( fScal[iSignal][ibtag] );
    hout->cd();
    hSignalCentral[ibtag]->Write();

    // create empty file just as marker
    ofstream markerFile;
    markerFile.open(Form("pack-%s-%s.txt",sbtag[ibtag].c_str(),scenario.c_str()),ios::app);
    markerFile << "Template for mass " << signalMass[iSignal] << std::endl;
    markerFile.close();
  }


  for (int iSyst=0; iSyst<nSyst; ++iSyst) {
    for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
      signalSystFiles[iSyst][iUpDown] = Form( "/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-3/%s/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola%s/%s_Sys%s/job_1/TripleBtagAnalysisM-%d_%s.root",
					      signalMode.c_str(),signalMass[iSignal],spacer.c_str(),systName[iSyst].c_str(),
					      upDownName[iUpDown].c_str(),signalMass[iSignal],IgorScen.c_str());
      std::cout << "Signal systematics file " << signalSystFiles[iSyst][iUpDown] << std::endl;
      TFile* fSigSys = new TFile( signalSystFiles[iSyst][iUpDown].c_str() );
      if ( fSigSys == NULL ) {
	std::cout << "Could not open signal syst file " << signalSystFiles[iSyst][iUpDown].c_str() << std::endl;
	return;
      }
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	hSignalSyst[iSyst][iUpDown][ibtag] = mergeSignal(fSigSys,Form(signalHistPattern.c_str(),sbtag[ibtag].c_str()),
							 Form("bbH_%s_%s",systName[iSyst].c_str(),upDownName[iUpDown].c_str()));
	if (hSignalSyst[iSyst][iUpDown][ibtag] == NULL) {
	  std::cout << "Reading of SignalSyst failed" << std::endl;
	  return;
	}
	// rename
	hSignalSyst[iSyst][iUpDown][ibtag]->SetName( Form("%s_%s_%s",
						     hSignalSyst[iSyst][iUpDown][ibtag]->GetName(),
							  systName[iSyst].c_str(),upDownName[iUpDown].c_str()) );
						     
	std::cout << "The merged hist has name " << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << std::endl;
	std::cout << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << " TotalContents=" << hSignalSyst[iSyst][iUpDown][ibtag]->GetSumOfWeights()
	      << std::endl;

// 	hSignalSyst[iSyst][iUpDown][ibtag]->Scale( fScal[iSignal][ibtag] );
	hout->cd();
	hSignalSyst[iSyst][iUpDown][ibtag]->Write();
      }
      fSigSys->Close();
    }
  }

  // here we should read the background templates
  std::string backgroundFile( Form("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/v1/%s/TripleBtagAnalysis_CR3_SF7/TripleBtagAnalysis_SF/TripleBtagAnalysis.root",SashaPath.c_str()) );
  std::cout << "Background central file : " << backgroundFile << std::endl;
  TFile* fBac = new TFile( backgroundFile.c_str() );
  if ( fBac == NULL ) {
    std::cout << "Could not open background central file " << signalFile.c_str() << std::endl;
    return;
  }
  TH2F* hBackgroundCentral[nbtag][ncateg][nfc];
  TH2F* hBackgroundCentralError[nbtag][ncateg][nfc];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	string templateCore("massBTagTemplatesCld/MassBTagTemplateCld");
	string hbSystName( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
				sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,3) );
	hBackgroundCentral[ibtag][icateg][ifc] = getTrigsAndMerge(fBac,Form("%s_%s_%s_Cat%dTpat%d",
									   templateCore.c_str(),
									   sfc[ifc].c_str(),
									   sbtag[ibtag].c_str(),icateg,3),nTCombData,tCombData);
	std::cout << "background central merged for " << templateCore << std::endl;
// 	hBackgroundCentral[ibtag][icateg][ifc] = (TH2F*) fBac->Get( hbSystName.c_str() );
	if ( hBackgroundCentral[ibtag][icateg][ifc] == NULL ) {
	  std::cout << "Hist not found: " << hbSystName << std::endl;
	  return;
	}
	templateCore.assign("errorMassBTagTemplates/ErrorMassBTagTemplate");
	string hbSystNameError( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
				sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,3) );
	hBackgroundCentralError[ibtag][icateg][ifc] = getTrigsAndMerge(fBac,Form("%s_%s_%s_Cat%dTpat%d",
									   templateCore.c_str(),
									   sfc[ifc].c_str(),
									   sbtag[ibtag].c_str(),icateg,3),nTCombData,tCombData);
// 	hBackgroundCentralError[ibtag][icateg][ifc] = (TH2F*) fBac->Get( hbSystNameError.c_str() );
	if ( hBackgroundCentralError[ibtag][icateg][ifc] == NULL ) {
	  std::cout << "Hist not found: " << hbSystNameError << std::endl;
	  return;
	}
	if (useTemplateError) {
	  // add the template error
	  std::cout << " ==== ErrorAdd ==== " << hBackgroundCentral[ibtag][icateg][ifc]->GetName() << std::endl;
	  for (int ibinx=1; ibinx<= (hBackgroundCentral[ibtag][icateg][ifc]->GetXaxis()->GetNbins()); ++ibinx) {
	    for (int ibiny=1; ibiny<= (hBackgroundCentral[ibtag][icateg][ifc]->GetYaxis()->GetNbins()); ++ibiny) {
	      float oldError = hBackgroundCentral[ibtag][icateg][ifc]->GetBinError(ibinx,ibiny);
	      float addError = hBackgroundCentralError[ibtag][icateg][ifc]->GetBinContent(ibinx,ibiny);
	      float newError = sqrt( oldError * oldError + addError * addError );
	      hBackgroundCentral[ibtag][icateg][ifc]->SetBinError(ibinx,ibiny,newError);
	    }
	  }
	}
	hout->cd();
	hBackgroundCentral[ibtag][icateg][ifc]->Write();
      }
    }
  }

  std::string backgroundSystFiles[nSyst][nUpDown];
  std::string systNameSasha[nSyst] = { "JES", "SFbc", "SFq" };
  std::string upDownNameSasha[nUpDown] = { "plus2", "minus2" };
  TH2F* hBackgroundSyst[nSyst][nUpDown][nbtag][ncateg][nfc];
  TH2F* hBackgroundSystError[nSyst][nUpDown][nbtag][ncateg][nfc];

  // in case we do not have an up/down, the default is always central value
  for (int iSyst=0; iSyst<nSyst; ++iSyst) {
    for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	for (int icateg=0; icateg<ncateg; ++icateg) {
	  for (int ifc=0; ifc<nfc; ++ifc) {
	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] = hBackgroundCentral[ibtag][icateg][ifc];
	  }
	}
      }
    }
  }

  // read explicitly for SFbc and SFq
  for (int iSyst=1; iSyst<nSyst; ++iSyst) {
    for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
      backgroundSystFiles[iSyst][iUpDown] = Form( "/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/v1/%s/TripleBtagAnalysis_CR3_SF7/TripleBtagAnalysis_%s%s/TripleBtagAnalysis.root",SashaPath.c_str(),systNameSasha[iSyst].c_str(),upDownNameSasha[iUpDown].c_str() );
      TFile* fBacSys = new TFile( backgroundSystFiles[iSyst][iUpDown].c_str() );
      std::cout << "Background systematics file " << backgroundSystFiles[iSyst][iUpDown] << std::endl;
      if ( fBacSys == NULL ) {
	std::cout << "Could not open background syst file " << backgroundSystFiles[iSyst][iUpDown] << std::endl;
	return;
      }
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	for (int icateg=0; icateg<ncateg; ++icateg) {
	  for (int ifc=0; ifc<nfc; ++ifc) {
	    string templateCore("massBTagTemplatesCld/MassBTagTemplateCld");
	    string hbSystName( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
				   sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,3) );
	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] = getTrigsAndMerge(fBac,Form("%s_%s_%s_Cat%dTpat%d",
									   templateCore.c_str(),
									   sfc[ifc].c_str(),
									   sbtag[ibtag].c_str(),icateg,3),nTCombData,tCombData);
// 	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] = 
// 	      (TH2F*) fBacSys->Get( hbSystName.c_str() );
	    if ( hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] == NULL ) {
	      std::cout << "Hist not found: " << hbSystName << std::endl;
	      return;
	    }
	    // rename
	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->SetName( Form("%s_%s_%s",
									hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName(),
									       systName[iSyst].c_str(),upDownName[iUpDown].c_str()) );  
	    std::cout << "The template syst name has name " << hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName()
		      << std::endl;

	    templateCore.assign("errorMassBTagTemplates/ErrorMassBTagTemplate");
	    string hbSystNameError( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
					 sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,3) );
	    hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc] = getTrigsAndMerge(fBac,Form("%s_%s_%s_Cat%dTpat%d",
									   templateCore.c_str(),
									   sfc[ifc].c_str(),
									   sbtag[ibtag].c_str(),icateg,3),nTCombData,tCombData);
// 	    hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc] = (TH2F*) fBac->Get( hbSystNameError.c_str() );
	    if ( hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc] == NULL ) {
	      std::cout << "Hist not found: " << hbSystNameError << std::endl;
	      return;
	    }
	    if (useTemplateError) {
	      // add the template error
	      std::cout << " ==== ErrorAdd ==== " << hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName() << std::endl;
	      for (int ibinx=1; ibinx<= (hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetXaxis()->GetNbins()); ++ibinx) {
		for (int ibiny=1; ibiny<= (hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetYaxis()->GetNbins()); ++ibiny) {
		  float oldError = hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetBinError(ibinx,ibiny);
		  float addError = hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc]->GetBinContent(ibinx,ibiny);
		  float newError = sqrt( oldError * oldError + addError * addError );
		  hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->SetBinError(ibinx,ibiny,newError);
		  if ( (icateg == 0) && (ifc == 0) ) {
		    std::cout << ibinx << " " << ibiny << " old error " << oldError
			      << " add error " << addError 
			      << " new error " << newError << std::endl;		    
		  }
		}
	      }
	    }

	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->SetName( 
									 Form("%s_syst%d_updown%d",
									      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName(),iSyst,iUpDown) );   

	    hout->cd();
	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->Write();
	  }
	}
      }
      //fBacSys->Close();
    }
  }

    // mass range bins
  int lowBinX = 3;
  int highBinX = 22;
//   int lowBinX = 6;
//   int highBinX = 9;

  TH2F* mjjEbtTemplate[nfc][nbtag][ncateg];
  // copy pointers to avoid messing up of code
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	mjjEbtTemplate[ifc][ibtag][icateg] = hBackgroundCentral[ibtag][icateg][ifc];
	std::cout << "ifc= " << ifc << " ibtag= " << ibtag
		  << " icateg= " << " copied from "
		  << hBackgroundCentral[ibtag][icateg][ifc]->GetName() << std::endl;
      }
    }
  }

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
  // read the hist-to-be-fitted
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtdata[ibtag] = getTrigsAndMerge(fBac,Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()),nTCombData,tCombData);

    if (mjjEbtdata[ibtag] == NULL) {
      std::cout << "Histogram not found: " << Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()) << std::endl;
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
    //dump2D( mjjEbtdata[ibtag] );
  }


  //int itpat = 0;
  int itpat = 3;

  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	mjjEbtTemplate[ifc][ibtag][icateg]->SetLineColor( tc[icateg+1+3*ifc] );
	mjjEbtTemplate[ifc][ibtag][icateg]->SetLineWidth(3);  
	
      }
    }
  }


  // set pointer for MC signal template
  TH2F* mjjEbtsignal[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtsignal[ibtag] = hSignalCentral[ibtag];
  }

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




#ifdef ADDNUISANCE
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
#endif


 // need a dummy tpat
 int itpatD = 4;

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
    // produce templates merging bbC and bbQ
    templateRef* mjjEbtTemplateMerge_bbCQ[nbtag];
//     mjjEbtTemplateMerge_bbCQ[ibtag] = templateRef::merge( "bbX",
// 							  Form("MassBTagTemplate_%s_%s_Cat%sTpat%d","CQ",sbtag[ibtag].c_str(),"2",itpatD),
// 							  mjjEbtTemplateRef[0][ibtag][2],mjjEbtTemplateRef[0][ibtag][2],0.1);
    mjjEbtTemplateMerge_bbCQ[ibtag] = templateRef::merge( "bbX",
							  Form("MassBTagTemplate_%s_%s_Cat%sTpat%d","CQ",sbtag[ibtag].c_str(),"2",itpatD),
							  mjjEbtTemplateRef[0][ibtag][2],mjjEbtTemplateRef[1][ibtag][2],0.1);
    //    std::cout << "Merge " << mjjEbtTemplateRef[0][ibtag][2]->hist->GetName()
    //      << " and " << mjjEbtTemplateRef[0][ibtag][2]->hist->GetName()



    // set up the templates configuration for this fit
    vector<templateRef> templates;
#undef ONETEMPLATE
#ifdef ONETEMPLATE
    templates.push_back( *mjjEbtTemplateRef[2][ibtag][2] ); // bbB
    //dump2D( mjjEbtTemplateRef[2][ibtag][2]->hist );
#else
//     templates.push_back( *mjjEbtTemplateRef[2][ibtag][0] ); // Bbb
    templates.push_back( *mjjEbtTemplateMerge01[2][ibtag] );  // (Bb)b
    templates.push_back( *mjjEbtTemplateMerge01[1][ibtag] );  // (Cb)b
    templates.push_back( *mjjEbtTemplateMerge01[0][ibtag] );  // (Qb)b
    templates.push_back( *mjjEbtTemplateRef[2][ibtag][2] ); // bbB
    //    templates.push_back( *mjjEbtTemplateRef[0][ibtag][2] ); // bbQ
    templates.push_back( *mjjEbtTemplateMerge_bbCQ[ibtag] ); // bbX
#endif

    // add the signal template
    if (fitSignal) {
      if (mjjEbtsignal[ibtag] != NULL) {
	templates.push_back(  templateRef("H",mjjEbtsignal[ibtag],0.001,2) );

#ifdef ADDNUISANCE
	///_SYST_BEGIN
	for (int sys=0;sys<(int) _plusTemplatesSgn.size();sys++)
	  templates.back().addNuisance(_plusTemplatesSgn[sys]._array[ibtag],_minusTemplatesSgn[sys]._array[ibtag],_sigmaTemplatesSgn[sys]);
	///_SYST_END
#endif
	

      } else {
	std::cout << "cannot add signal template for " << sbtag[ibtag] << std::endl;
      }

    }



    // list the templates
    std::cout << "Templates being fitted: " << std::endl;
    for (vector<templateRef>::iterator trit=templates.begin(); trit != templates.end(); ++trit) {
      std::cout << "   " << (*trit).name 
		<< "  inival " << (*trit).inival << " nNuisance " << (*trit).nNuisance
		<< std::endl;
    }

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

    // count number of active nuisance parameters
    int nActiveNuisance = 0;
    for (int iNuisance=0; iNuisance<3; ++iNuisance) {
      if ( activeNuisance[iNuisance] ) ++nActiveNuisance;
    }

    std::cout << "Number of active nuisance parameters: " << nActiveNuisance << std::endl
	      << "  " << activeNuisance[0]
	      << "  " << activeNuisance[1]
	      << "  " << activeNuisance[2]
	      << std::endl;
    
    gTemplateFitter = new templateFitter(mjjEbtdata[ibtag], mc, &templates,  nActiveNuisance );
    //gTemplateFitter = new templateFitter(mjjEbtdata[ibtag], mc, &templates,  0 );  // first test without any nuisance parameters
    
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
      gTemplateFitter->GetResult(iH, HFraction[ibtag][iSignal], HFractionError[ibtag][iSignal]);
      //gTemplateFitter->GetResult(iH, value, error);
      // determine the fitting penalty
      fPenalty =  HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights() / sqrt( mjjEbtdata[ibtag]->GetSumOfWeights() );
      std::cout << "deltaS / sqrt(B) = " << error * mjjEbtdata[ibtag]->GetSumOfWeights()
		<< "  /  " << sqrt( mjjEbtdata[ibtag]->GetSumOfWeights() ) << "  =  "
		<< fPenalty << std::endl;
      std::cout << " HFraction=" << HFraction[ibtag][iSignal] 
		<< " +- " << HFractionError[ibtag][iSignal] 
		<< " SumW= " << mjjEbtdata[ibtag]->GetSumOfWeights()
		<< std::endl;
      HSeen[ibtag][iSignal] = HFraction[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights();
      HSeenError[ibtag][iSignal] = HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights();

      HEfficiency[ibtag][iSignal] = efficiency[iSignal][ibtag];

//       // nominal cross section
//       TH1F* histXSect = (TH1F*) fSig->Get("xsection/xsect");
//       if ( histXSect == NULL) {
// 	std::cout << "xsection/xsect" << " not found" << std::endl;
// 	return;
//       }
//       nominalXSect[ibtag][iSignal] = histXSect->GetBinContent(1);

      HCorr[ibtag][iSignal] = HSeen[ibtag][iSignal] / HEfficiency[ibtag][iSignal];
      HCorrError[ibtag][iSignal] = HSeenError[ibtag][iSignal] / HEfficiency[ibtag][iSignal];
      HXSect[ibtag][iSignal] = 0.001 * HCorr[ibtag][iSignal] / intLumi;
      HXSectError[ibtag][iSignal] = 0.001 * HCorrError[ibtag][iSignal] / intLumi;
      HXSectUpLimit[ibtag][iSignal] = uplimit( HXSect[ibtag][iSignal], HXSectError[ibtag][iSignal]);
      tanBetaUpLimit[ibtag][iSignal] = sqrt( HXSectUpLimit[ibtag][iSignal] / nominalXSect[ibtag][iSignal] ) * 20;
      std::cout << "Efficiency = " << HEfficiency[ibtag][iSignal] << std::endl
		<< "Higgs observed = " << HSeen[ibtag][iSignal] << " +- " << HSeenError[ibtag][iSignal] << std::endl
		<< "Higgs corrected = " << HCorr[ibtag][iSignal] << " +- " << HCorrError[ibtag][iSignal] << std::endl
		<< "Higgs XSection = " << HXSect[ibtag][iSignal] << " +- " << HXSectError[ibtag][iSignal] << " pb " << std::endl
		<< "Higgs XSection < " << HXSectUpLimit[ibtag][iSignal] << " pb " << std::endl
		<< "tanBeta        < " << tanBetaUpLimit[ibtag][iSignal] << std::endl
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

    canvas->Print(Form("templateFitStacked_%s_M-%d_X.png",sbtag[ibtag].c_str(),signalMass[iSignal]));
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

#ifdef DONOTHING
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
#endif
  }  // loop over ibtag
}


void customFitN() {
  int iMass = 1;
  customFitN_Mass( iMass );
  return;
}
