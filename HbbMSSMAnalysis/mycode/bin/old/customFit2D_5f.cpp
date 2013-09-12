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
#include <cmath>

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
  templateRef() {
    std::cout << "Warning: templateRef default constructor called" << std::endl;
  }
  templateRef(const char* theName,TH2F* theHist,const double theInival = 1,const int theColor=-1) : name(theName), hist(theHist), 
												    inival( theInival ), color( theColor) {
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
      hist->Scale( 1./hist->GetSumOfWeights() );
    } else {
      std::cout << "templateRef::normalize(): bad sum of weights " << hist->GetSumOfWeights() << std::endl;
    }
  }
};

double uplimit(double x0,double sigma) {
  double CL = 0.95;

  double dxOverSig = TMath::ErfInverse( CL + (1 - CL) * TMath::Erf( -x0/(sigma*sqrt(2)) ) );
  
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

  templateFitter(TH2F* theData,TObjArray* theMc,std::vector<templateRef>* theTemplates) {
    data = theData;
    TIter next(theMc);
    TH2F* theTemp;
    while ( (theTemp = (TH2F*) next()) ) {
      mc.push_back( theTemp );
      minConstraint.push_back(0);
      maxConstraint.push_back(0);
    }
//     for (std::vector<templateRef>::iterator tit=theTemplates.begin(); tit != theTemplates.end(); ++tit) {
//       templates.push_back( *tit );
//     }
    templates = theTemplates;
    fitLow = 1;  // TH convention, count from 1...nbins
    fitHigh = theData->GetNbinsX();
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
    TMinuit* gMinuit = new TMinuit(mc.size());
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    arglist[0] = 1;
    Int_t ierflg = 0;
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
    Double_t* vstart = new Double_t[mc.size()];
    Double_t* vstep = new Double_t[mc.size()];
    // for normalization of the initial fractions
    double iniTot = 0;
    for (unsigned int i=0; i<mc.size(); ++i) {
      iniTot += (*templates)[i].inival;
    }
    for (unsigned int i=0; i<mc.size(); ++i) {
      vstart[i] = (*templates)[i].inival / iniTot;
      vstep[i] = vstart[i] / 100;
      //gMinuit->mnparm(i,"par",vstart[i],vstep[i],0,2*data->GetSumOfWeights() / mc[i]->GetSumOfWeights(),ierflg);
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
      for (int imc=0; imc<npar; ++imc) {
	delta -= par[imc] * T *  gTemplateFitter->GetMc(imc)->GetBinContent(binx,biny);
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


void findYMinMax(int ipos,int nPoints,double* theY,double* theEY,double& yMin,double& yMax) {
  double S = 0;
  for (int iPoint=0; iPoint<nPoints; ++iPoint) {
    double YHere = theY[iPoint] + theEY[iPoint];
    if (ipos == 0) YHere = abs(theY[iPoint]) + abs(theEY[iPoint]);
    YHere *= 1.2; // margin
    if (YHere>S) S = YHere;
  }
  if (ipos == 1) {
    yMin = 0;
    yMax = S;
  } else {
    yMin = -S;
    yMax = S;
  }
  return;
}

void customFit2D_5f() {
  // fit with corrected templates using the minuit fitter
  setTDRStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  canvas = new TCanvas ("cg1","PadX",10,10,800,600);
  gStyle->SetPadColor(0);
  canvas->SetFillColor(0);

  // this is for combination of triggers
  const int nTCombData = 3;
  std::string tCombData[nTCombData] = {"Trig0", "Trig1", "Trig2"};
  //   const int nTCombData = 1;
  //   std::string tCombData[nTCombData] = {"Trig0"};
  
  string bbPurity("DataDriven"); //  "MC", "DataDrivenR", "None"
  string Var("Var2");
  bool onlineBtagCorr = true;
  bool fitSignal = false;
  bool useTemplateError = true;

  const int ncateg = 3;
  const int ncorr=2;
  const int nfc = 3;  
  //const int nbtag = 4;
  const int nbtag = 1;

  //const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  const std::string sbtag[nbtag] = { "CSVT"};

  // mass range bins
  int lowBinX = 3;
  int highBinX = 22;
  //   int lowBinX = 6;
  //   int highBinX = 9;
  
  TH2F* mjjEbtTemplate[nfc][nbtag][ncateg];
  TH2F* mjjEbtTemplateError[nfc][nbtag][ncateg];
  string tFlav[3] = {"Uds", "C", "B"};
  string sFlav[3] = {"Q", "C", "B"};
  
  string sfc[nfc] = { "q", "c", "b" };
  
  // text summary files
  ofstream* sumFile[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    sumFile[ibtag] = new ofstream(Form("sumFile_%s.txt",sbtag[ibtag].c_str()));
  }

  // location of signal templates
  const int nSignal=7;
  //const int nSignal=1;
  int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };
  //int signalMass[nSignal] = { 120 };

  // summary variables across the masses
  double intLumi = 2.66794; // in fb-1
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
  

  // loop over the signal masses
  for (int iSignal=0; (fitSignal && (iSignal<nSignal)) || (iSignal<1); ++iSignal) {
    
    const int maxfa = 2;
    TFile* fa[maxfa];

    // this is for triple tag mass spectrum and background templates
    fa[0] = new TFile("../TripleBtagAnalysis.root");
//       fa[1] = new TFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/CMSSW424p1/v0.2/mc/Fall11/SUSYBBHToBB_M-120_7TeV-pythia6-tauola/TripleBtagAnalysisV3/TripleBtagAnalysis.root");
    fa[1] = new TFile(Form("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-2/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola/SF/job_1/TripleBtagAnalysisM-%d_low.root",signalMass[iSignal],signalMass[iSignal]));
    for (int ifa=0; ifa<maxfa; ++ifa) {
      if (fa[ifa] == NULL) {
	std::cout << "Problem opening file " << ifa << std::endl;
	return;
      }
    }

    // mass label (for plots)
    string cMass;
    if (fitSignal) {
      cMass.assign( Form("mH_%d",signalMass[iSignal]) );
    } else {
      cMass.assign( "BGonly" );
    }

    // now read the templates (taken from prottemplate2.C)
#define MC
#ifdef MC
    std::string tSel("AllTrigWeighted");
#else
    //std::string tSel("Trig2");
    std::string tSel("Trig0");
#endif

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
    int itpat = 3;  // this is the default (=use obt corrected templates)
    
    string templateCore;
    if (bbPurity == "DataDrivenR") {
      templateCore.assign("massBTagTemplatesNonbbR/MassBTagTemplateCldR");
    } else if (bbPurity == "DataDriven") {
      templateCore.assign("massBTagTemplatesCld/MassBTagTemplateCld");
    } else if (bbPurity == "MC") {
      templateCore.assign("massBTagTemplates/MassBTagTemplate");
    } else if (bbPurity == "None") {
      templateCore.assign("massBTagTemplatesCld/MassBTagTemplateUncld");
    }
    string templateCoreError("errorMassBTagTemplates/ErrorMassBTagTemplate");

    std::cout << "Var = " << Var << "   bbPurity = " << bbPurity << std::endl;

    for (int ifc=0; ifc<nfc; ++ifc) {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	for (int icateg=0; icateg<ncateg; ++icateg) {
	  if (onlineBtagCorr) {
	    itpat = 3;
	  } else {
	    itpat = icateg;
	  }
	  mjjEbtTemplate[ifc][ibtag][icateg] = getTrigsAndMerge(fa[0],Form("%s_%s_%s_Cat%dTpat%d",
									   templateCore.c_str(),
									   sfc[ifc].c_str(),
									   sbtag[ibtag].c_str(),icateg,itpat),nTCombData,tCombData);
	  
	  // 	mjjEbtTemplate[ifc][ibtag][icateg] = (TH2F*) fa[0]->Get(Form("massBTagTemplates/MassBTagTemplate_%s_%s_Cat%dTpat%d%s",
	  // 									    sfc[ifc].c_str(),
	  // 									    sbtag[ibtag].c_str(),icateg,itpat,tSel.c_str()));
	  
	  if (mjjEbtTemplate[ifc][ibtag][icateg] == NULL) {
	    std::cout << "Histogram not found: " << Form("%s_%s_%s_Cat%dTpat%d%s",
							 templateCore.c_str(),
							 sfc[ifc].c_str(),
							 sbtag[ibtag].c_str(),icateg,itpat,tSel.c_str())
		      << std::endl;
	    return;
	  }

	  if (useTemplateError) {
	    mjjEbtTemplateError[ifc][ibtag][icateg] = getTrigsAndMerge(fa[0],Form("%s_%s_%s_Cat%dTpat%d",
										  templateCoreError.c_str(),
										  sfc[ifc].c_str(),
										  sbtag[ibtag].c_str(),icateg,itpat),nTCombData,tCombData);
	    if (mjjEbtTemplateError[ifc][ibtag][icateg] == NULL) {
	      std::cout << "Histogram not found: " << Form("%s_%s_%s_Cat%dTpat%d%s",
							   templateCoreError.c_str(),
							   sfc[ifc].c_str(),
							   sbtag[ibtag].c_str(),icateg,itpat,tSel.c_str())
			<< std::endl;
	      return;
	    }
	    // add the template error
	    std::cout << " ==== ErrorAdd ==== " << mjjEbtTemplate[ifc][ibtag][icateg]->GetName() << std::endl;
	    for (int ibinx=1; ibinx<= (mjjEbtTemplate[ifc][ibtag][icateg]->GetXaxis()->GetNbins()); ++ibinx) {
	      for (int ibiny=1; ibiny<= (mjjEbtTemplate[ifc][ibtag][icateg]->GetYaxis()->GetNbins()); ++ibiny) {
		std::cout << ibinx << " " << ibiny << " old error " << mjjEbtTemplate[ifc][ibtag][icateg]->GetBinError(ibinx,ibiny)
			  << " add error " <<  mjjEbtTemplateError[ifc][ibtag][icateg]->GetBinContent(ibinx,ibiny)
			  << std::endl;
		float oldError = mjjEbtTemplate[ifc][ibtag][icateg]->GetBinError(ibinx,ibiny);
		float addError = mjjEbtTemplateError[ifc][ibtag][icateg]->GetBinContent(ibinx,ibiny);
		float newError = sqrt( oldError * oldError + addError * addError );
		mjjEbtTemplate[ifc][ibtag][icateg]->SetBinError(ibinx,ibiny,newError);
	      }
	    }
	  }

	  mjjEbtTemplate[ifc][ibtag][icateg]->SetLineColor( tc[icateg+1+3*ifc] );
	  mjjEbtTemplate[ifc][ibtag][icateg]->SetLineWidth(3);  
	  
	}
      }
    }

    // now read the MC signal template
    TH2F* mjjEbtsignal[nbtag];
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {

      if (sbtag[ibtag] != "CSVT") {
	std::cout << "Currently only CSVT is supported by Igor. Sorry." << std::endl;
	return;
      }
      mjjEbtsignal[ibtag] = (TH2F*) fa[1]->Get(Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0"));
      if ( mjjEbtsignal[ibtag] == NULL ) {
	std::cout << "Signal histogram not found: " << Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),"Trig0") << std::endl;
      }
      else {
	// change name to avoid duplication
	mjjEbtsignal[ibtag]->SetName( Form("%s_Signal",mjjEbtsignal[ibtag]->GetName()) );
      }
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

    // need a dummy tpat
    int itpatD = 4;

    for (int ibtag = 0; ibtag < nbtag; ++ibtag) {

      // produce templates merging the first two categories
      TH2F* mjjEbtTemplateMerge01[nfc][nbtag];
      for (int ifc=0; ifc<nfc; ++ifc) {
	std::cout << "Merge for" 
		  << " ifc=" << ifc << " ibtag=" << ibtag << " icorr=" << icorr << " itpatD=" << itpatD
		  << std::endl;
	mjjEbtTemplateMerge01[ifc][ibtag] = new TH2F( *mjjEbtTemplate[ifc][ibtag][0] );
	mjjEbtTemplateMerge01[ifc][ibtag]->SetName(  Form("MassBTagTemplate_%s_%s_Cat%sTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),"01",itpatD) );
	mjjEbtTemplateMerge01[ifc][ibtag]->SetTitle( Form("MassBTagTemplate_%s_%s_Cat%sTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),"01",itpatD) );
	mjjEbtTemplateMerge01[ifc][ibtag]->Add( mjjEbtTemplate[ifc][ibtag][1] );
      }

      // produce templates merging bbC and bbQ
      TH2F* mjjEbtTemplateMerge_bbCQ[nbtag];
      mjjEbtTemplateMerge_bbCQ[ibtag] = new TH2F( *mjjEbtTemplate[0][ibtag][2] );
      mjjEbtTemplateMerge_bbCQ[ibtag]->SetName(  Form("MassBTagTemplate_%s_%s_Cat%sTpat%d","CQ",sbtag[ibtag].c_str(),"2",itpatD) );
      mjjEbtTemplateMerge_bbCQ[ibtag]->SetTitle( Form("MassBTagTemplate_%s_%s_Cat%sTpat%d","CQ",sbtag[ibtag].c_str(),"2",itpatD) );
      mjjEbtTemplateMerge_bbCQ[ibtag]->Add( mjjEbtTemplate[1][ibtag][2] );
      
      vector<templateRef> templates;

      if (Var == "") {
	//  templates.push_back( templateRef("Qbb",mjjEbtTemplate[0][ibtag][0],1) );
	templates.push_back( templateRef("Bbb",mjjEbtTemplate[2][ibtag][0],1) );
	//  templates.push_back( templateRef("bBb",mjjEbtTemplate[2][ibtag][1],1) );
	//  templates.push_back( templateRef("(Bb)b",mjjEbtTemplateMerge01[2][ibtag],0.7) );
	templates.push_back( templateRef("(Cb)b",mjjEbtTemplateMerge01[1][ibtag],0.1) );
	templates.push_back( templateRef("(Qb)b",mjjEbtTemplateMerge01[0][ibtag],0.1) );

	//   templates.push_back( templateRef("bBb",mjjEbtTemplate[2][ibtag][1],1) );
	templates.push_back( templateRef("bbB",mjjEbtTemplate[2][ibtag][2],0.3) );
      
	templates.push_back( templateRef("bbQ",mjjEbtTemplate[0][ibtag][2],0.1) );
	//  templates.push_back( templateRef("bbX",mjjEbtTemplateMerge_bbCQ[ibtag],0.1) );
      } else if (Var == "Var2") {
	templates.push_back( templateRef("(Bb)b",mjjEbtTemplateMerge01[2][ibtag],0.7) );
	templates.push_back( templateRef("(Cb)b",mjjEbtTemplateMerge01[1][ibtag],0.1) );
	templates.push_back( templateRef("(Qb)b",mjjEbtTemplateMerge01[0][ibtag],0.1) );
	templates.push_back( templateRef("bbB",mjjEbtTemplate[2][ibtag][2],0.3) );
	templates.push_back( templateRef("bbX",mjjEbtTemplateMerge_bbCQ[ibtag],0.1) );
      } else {
	std::cout << "No such Var scenario " << Var << std::endl;
	return;
      }
      
      bool constrainSignal = false;   // set to true if signal fraction must be non-negative
      // only works if signal template is labeled H

      


      //   templates.push_back( templateRef("bbC",mjjEbtTemplate[1][ibtag][2],1) );
      //   templates.push_back( templateRef("bbQ",mjjEbtTemplate[0][ibtag][2],1) ); 
  

      // add the signal template
      if ( fitSignal ) {
	if (mjjEbtsignal[ibtag] != NULL) {
	  templates.push_back(  templateRef("H",mjjEbtsignal[ibtag],0.001,2) );
	} else {
	  std::cout << "cannot add signal template for " << sbtag[ibtag] << std::endl;
	}
      }
      
      // normalize the templates
      for (vector<templateRef>::iterator trit=templates.begin(); trit != templates.end(); ++trit) {
	(*trit).normalize();
      }
      
      TObjArray *mc = new TObjArray();
      for (vector<templateRef>::iterator tit=templates.begin(); tit != templates.end(); ++tit) {
	mc->Add( (*tit).hist );
      }
      std::cout << "Number of templates = " << mc->GetEntries() << std::endl;
      
      gTemplateFitter = new templateFitter(mjjEbtdata[ibtag], mc, &templates);
      
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
      fPenalty[ibtag][iSignal] = -1;
      
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
	canvas->Print(Form("%s_TemplateProX_%s_%d.png",cMass.c_str(),sbtag[ibtag].c_str(),iT));
	++iT;
      }

      // print Higgs statistics (if H template was fitted)
      if (iH>=0) {
	gTemplateFitter->GetResult(iH, HFraction[ibtag][iSignal], HFractionError[ibtag][iSignal]);
	// determine the fitting penalty
	fPenalty[ibtag][iSignal] = HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights() / sqrt( mjjEbtdata[ibtag]->GetSumOfWeights() );
	std::cout << "deltaS / sqrt(B) = " << HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights()
		  << "  /  " << sqrt( mjjEbtdata[ibtag]->GetSumOfWeights() ) << "  =  "
		  << fPenalty[ibtag][iSignal] << std::endl;
	// expected signal if MC Xsection is correct
	//double HExpected = templates[iH].sumW * 0.501 * 111./38.;
	HSeen[ibtag][iSignal] = HFraction[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights();
	HSeenError[ibtag][iSignal] = HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSumOfWeights();
	//std::cout << "Higgs scaling factor: " << HSeen / HExpected << " +- " << HSeenError / HExpected << std::endl;
	//double HScal95 = uplimit( HSeen / HExpected, HSeenError / HExpected );
	//std::cout << "Higgs scaling factor UL: " << HScal95 << std::endl;
	//std::cout << "tan-beta             UL: " << 20 * sqrt( HScal95 ) << std::endl; 
	//double HMCGen = 1051000;
	//double HMCAccepted = mjjEbtsignal[ibtag]->GetEntries();
	// use Igor's efficiency
	TH1F* histEff[2];
	for (int ite=0; ite<2; ++ite) {
	  histEff[ite] = (TH1F*) fa[1]->Get(Form("TrigEff/Eff%sTrig%d",sbtag[ibtag].c_str(),ite));
	  if ( histEff[ite] == NULL) {
	    std::cout << Form("TrigEff/Eff%sTrig%d",sbtag[ibtag].c_str(),ite) << " not found" << std::endl;
	    return;
	  }
	}
	//double HEfficiency = HMCAccepted / HMCGen;
	HEfficiency[ibtag][iSignal] = (1688.68/2667.94) * histEff[0]->GetBinContent(1)
	  + (979.26/2667.94) * histEff[1]->GetBinContent(1);

	// nominal cross section
	TH1F* histXSect = (TH1F*) fa[1]->Get("xsection/xsect");
	if ( histXSect == NULL) {
	  std::cout << "xsection/xsect" << " not found" << std::endl;
	  return;
	}
	nominalXSect[ibtag][iSignal] = histXSect->GetBinContent(1);

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
			  sbtag[ibtag].c_str(),model.c_str(),gTemplateFitter->GetChisquare(),nInactive,fPenalty[ibtag][iSignal],HXSectError[ibtag][iSignal]) << std::endl;
	
	// write statistics to a file
	if (iSignal == 0) {
	  *sumFile[ibtag] << std::setw(5) << "mH"
			  << " " << std::setw(12) << "HEff"
			  << " " << std::setw(12) << "nomXSect"
			  << " " <<  std::setw(12) << "fraction" << " " << std::setw(10) << "error" 
			  << " " << std::setw(10) << "HSeen" << " " << std::setw(10) << "HSeenErr" 
			  << " " << std::setw(10) << "HCorr" << " " << std::setw(10) << "HCorrErr"
			  << " " << std::setw(10) << "HXSect" << " " << std::setw(10) << "HXSectErr"
			  << " " << std::setw(10) << "HXSectUL"
			  << " " << std::setw(10) << "tanBetaUL"
			  << " " << std::setw(10) << "fPenalty"
			  << " " << std::setw(10) << "intLumi"
			  << std::endl;
	}
	*sumFile[ibtag] << std::setw(5) << signalMass[iSignal]
			<< " " << std::setw(12) << HEfficiency[ibtag][iSignal]
			<< " " << std::setw(12) << nominalXSect[ibtag][iSignal]
			<< " " <<  std::setw(12) << HFraction[ibtag][iSignal] << " " << std::setw(10) << HFractionError[ibtag][iSignal]
			<< " " << std::setw(10) << HSeen[ibtag][iSignal] << " " << std::setw(10) << HSeenError[ibtag][iSignal] 
			<< " " << std::setw(10) << HCorr[ibtag][iSignal] << " " << std::setw(10) << HCorrError[ibtag][iSignal]
			<< " " << std::setw(10) << HXSect[ibtag][iSignal] << " " << std::setw(10) << HXSectError[ibtag][iSignal]
			<< " " << std::setw(10) << HXSectUpLimit[ibtag][iSignal]
			<< " " << std::setw(10) << tanBetaUpLimit[ibtag][iSignal]
			<< " " << std::setw(10) << fPenalty[ibtag][iSignal]
			<< " " << std::setw(10) << intLumi
			<< std::endl;
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
      
      canvas->Print(Form("%s_templateFitStacked_%s_X.png",cMass.c_str(),sbtag[ibtag].c_str()));
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
      
      canvas->Print(Form("%s_templateFitStacked_%s_Y.png",cMass.c_str(),sbtag[ibtag].c_str()));

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
	canvas->Print(Form("%s_templateFitSliceX_%s_%d.png",cMass.c_str(),sbtag[ibtag].c_str(),biny));
      }
    }
    for (int ifa=0; ifa<maxfa; ++ifa) {
      fa[ifa]->Close();
      delete fa[ifa];
    }
  }

  if (fitSignal) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      // statistics histograms across masses
      const int nStati = 100;
      TH1D* statH[nStati];
      TGraphErrors* statEH[nStati];
      double Pos[nStati];
      string sName[nStati];
      string sTitle[nStati];
      double* theStat[nStati];
      double* theError[nStati];
      int jStati = 0;
      sName[jStati] = "Efficiency";
      sTitle[jStati] = "Higgs efficiency vs mass";
      theStat[jStati] = HEfficiency[ibtag];
      theError[jStati] = NULL;
      Pos[jStati] = 1;

      ++jStati;
      sName[jStati] = "nominalXSection";
      sTitle[jStati] = "Higgs theory XSect (tanBeta=20) vs mass";
      theStat[jStati] = nominalXSect[ibtag];
      theError[jStati] = NULL;
      Pos[jStati] = 1;
      
      ++jStati;
      sName[jStati] = "HFraction";
      sTitle[jStati] = "Higgs fraction in fit";
      theStat[jStati] = HFraction[ibtag];
      theError[jStati] = HFractionError[ibtag];
      Pos[jStati] = 0;

      ++jStati;
      sName[jStati] = "HSeen";
      sTitle[jStati] = "Higgs count seen in fit";
      theStat[jStati] = HSeen[ibtag];
      theError[jStati] = HSeenError[ibtag];
      Pos[jStati] = 0;

      ++jStati;
      sName[jStati] = "HCorr";
      sTitle[jStati] = "Higgs count corrected from fit";
      theStat[jStati] = HCorr[ibtag];
      theError[jStati] = HCorrError[ibtag];
      Pos[jStati] = 0;
      
      ++jStati;
      sName[jStati] = "HXSect";
      sTitle[jStati] = "Higgs cross section [pb] from fit";
      theStat[jStati] = HXSect[ibtag];
      theError[jStati] = HXSectError[ibtag];
      Pos[jStati] = 0;

      ++jStati;
      sName[jStati] = "HXSectUpLimit";
      sTitle[jStati] = "Higgs cross section upper limit [pb] from fit";
      theStat[jStati] = HXSectUpLimit[ibtag];
      theError[jStati] = NULL;
      Pos[jStati] = 1;

      ++jStati;
      sName[jStati] = "tanBetaUpLimit";
      sTitle[jStati] = "tan(beta) upper limit from fit";
      theStat[jStati] = tanBetaUpLimit[ibtag];
      theError[jStati] = NULL;
      Pos[jStati] = 1;

      ++jStati;
      sName[jStati] = "fPenalty";
      sTitle[jStati] = "Fitting penalty";
      theStat[jStati] = fPenalty[ibtag];
      theError[jStati] = NULL;
      Pos[jStati] = 1;

      int mStati = jStati + 1;  // number of positions filled
      
      // fill the TGraphError
      for (int iStati=0; iStati<mStati; ++iStati) {
	double theX[nSignal];
	double theY[nSignal];
	double theEX[nSignal];
	double theEY[nSignal];
    
	for (int iSignal=0; iSignal<nSignal; ++iSignal) {
	  theX[iSignal] = signalMass[iSignal];
	  theEX[iSignal] = 0;
	  theY[iSignal] = theStat[iStati][iSignal];
	  theEY[iSignal] = 0;
	  if (theError[iStati] != NULL) {
	    theEY[iSignal] = theError[iStati][iSignal];
	  }
	}
	statEH[iStati] = new TGraphErrors(nSignal,theX,theY,theEX,theEY);

	double theYMin = 0;
	double theYMax = 0;
	findYMinMax(Pos[iStati],nSignal,theY,theEY,theYMin,theYMax);
	TH2F hist("hdum","hdum",40,0,400,40,theYMin,theYMax);
	hist.GetXaxis()->SetTitle("M_{H} [GeV/c^{2}]");
	hist.GetYaxis()->SetTitle(sName[iStati].c_str());    
	hist.SetStats(kFALSE);
	hist.Draw();
	//       statEH[iStati]->GetXaxis()->SetTitle("M_{H} [GeV/c^{2}]");
	//       statEH[iStati]->GetYaxis()->SetTitle(sName[iStati].c_str());
	statEH[iStati]->Draw("EPL,SAME");
	TLatex statLegend;
	statLegend.SetNDC();
	statLegend.SetTextAlign(12);
	statLegend.SetTextSize(0.035);
	statLegend.DrawLatex(0.7,0.9,sbtag[ibtag].c_str());
	canvas->Print(Form("%s_%s.png",sName[iStati].c_str(),sbtag[ibtag].c_str()));
      }

      // pull histogram
      TH1D* pullFraction = new TH1D(Form("pullFraction_%s",sbtag[ibtag].c_str()),"pull of fraction",40,-4,4);
      pullFraction->GetXaxis()->SetTitle("HFraction / HFractionError");
      pullFraction->GetYaxis()->SetTitle("N");
      for (int iSignal=0; iSignal<nSignal; ++iSignal) {
	if ( HFractionError[ibtag][iSignal]  != 0) {
	  pullFraction->Fill( HFraction[ibtag][iSignal] /  HFractionError[ibtag][iSignal] );
	}
      }
      pullFraction->Draw("");
      canvas->Print(Form("pullFraction_%s.png",sbtag[ibtag].c_str()));
    }

    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      sumFile[ibtag]->close();
    }
  }
}
