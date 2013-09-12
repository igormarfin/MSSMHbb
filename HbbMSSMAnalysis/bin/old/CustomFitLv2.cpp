//
//  This is the looper version of customFitY.cpp, for evaluation of 
//  the expected limit
//
#define SUMMARYCARD 1

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
#include <TMatrixDSym.h>

/// systematics read-out
#include <map>
#include <TString.h>
#include <TObjArray.h>
#include <TSystem.h>
#include <TRandom3.h>

#include "Analysis/Utilities/interface/getHbbCfg.h"
#include "Analysis/Utilities/interface/HbbXsec.h"

TCanvas* canvas;
class templateFitter;
templateFitter* gTemplateFitter;  // pointer needs to be global, to be called in fcn

double tanBeta = 20;

using std::vector;
using std::cout;
using std::endl;

HbbXsec hbb_xs("/scratch/hh/current/cms/user/walsh/HbbMSSMAnalysis/TheoryXsectionScans/out.mhmax_mu200_7_nnlo.tanBeta_gte1.root");

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
    mergedData->SetName( Form("%s%s",theName.c_str(),newname) );
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
  std::string name() {
    std::string sFlav[3] = {"Q", "C", "B"};
    std::string sLegend("bbb");
    sLegend.replace(categ,1,sFlav[flav]);
    return sLegend;
  }
};


class templateRef {
public:
  std::string name;
  TH2F* hist;
  double inival;
  int color;
  double sumW; // total content before renormalization
  int nNuisance;  // number of nuisance parameters
  vector<TH2F*> tempPlus; // vector of templates for +nsig of nuisance parameter
  vector<TH2F*> tempMinus; // vector of templates for -nsig of nuisance parameter
  vector<float> nSigma; // vector of nsig for nuisance parameters
  vector<float> NormSlope;

  float getNormSlope(const unsigned int iNuisance) {
    if (iNuisance > NormSlope.size()) {
      std::cout << "templateRef::getNormSlope: bad iNuisance=" << iNuisance << std::endl;
      return -999999;
    }
    return NormSlope[iNuisance];
  }

  templateRef() {
    std::cout << "Warning: templateRef default constructor called" << std::endl;
  }
  templateRef(const char* theName,TH2F* theHist,const double theInival = 1,const int theColor=-1) : name(theName), 
												    inival( theInival ), 
												    color( theColor), nNuisance(0) {
//     std::cout << "templateRef constructor: theHist=" << theHist << " theColor=" << theColor << std::endl;

    // to be safe, make a copy of the histogram
    hist = new TH2F( *theHist );
    hist->SetName( Form("%s_tR", theHist->GetName()) );
    std::cout << "templateRef: Created a copy " << hist->GetName() << std::endl; 

    if (theColor == -1) {
      color = theHist->GetLineColor();
    } else {
      theHist->SetLineColor(color);
    }
    sumW = theHist->GetSum();
  }
  ~templateRef() {
    //std::cout << "This is the TemplateRef destructor for " << name << std::endl;
    //std::cout << "Delete " << hist->GetName() << std::endl;
    delete hist;
    for (unsigned int iNuisance=0; iNuisance<tempPlus.size(); ++iNuisance) {
      //std::cout << "Delete" << tempPlus[iNuisance]->GetName() << std::endl;
      delete tempPlus[iNuisance];
      //std::cout << "Delete" << tempMinus[iNuisance]->GetName() << std::endl;
      delete tempMinus[iNuisance];
    }
  }
  void normalize(const int centralMode = 1) {
    if (hist->GetSum() >0) {
      std::cout << hist->GetSum() << std::endl;
      float scaleFacCentral = 1. / hist->GetSum();
      hist->Scale( scaleFacCentral );
      for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
	float scaleFac = scaleFacCentral;
	if (centralMode == 0) {
	  if (tempPlus[iNuisance]->GetSum() >0) {
	    scaleFac = 1. / tempPlus[iNuisance]->GetSum();
	  } else {
	    std::cout << "templateRef::normalize: Error: Plus template non-normalizable" << std::endl;
	    return;
	  }
	}
	tempPlus[iNuisance]->Scale( scaleFac );

	scaleFac = scaleFacCentral;
	if (centralMode == 0) {
	  if (tempMinus[iNuisance]->GetSum() >0) {
	    scaleFac = 1. / tempMinus[iNuisance]->GetSum();
	  } else {
	    std::cout << "templateRef::normalize: Error: Minus template non-normalizable" << std::endl;
	    return;
	  }
	}
	tempMinus[iNuisance]->Scale( scaleFac );
      }
    } else {
      std::cout << "templateRef::normalize(): bad sum of weights " << hist->GetSum() << std::endl;
    }
  }
  void addNuisance(TH2F* theTempPlus, TH2F* theTempMinus, float theNSigma) {
    // create a copy, just to be safe
    TH2F* theTempPlusCopy = new TH2F( *theTempPlus );
    theTempPlusCopy->SetName( Form("%s_%s_COPY",theTempPlus->GetName(),name.c_str()) );
    TH2F* theTempMinusCopy = new TH2F( *theTempMinus );
    theTempMinusCopy->SetName( Form("%s_%s_COPY",theTempMinus->GetName(),name.c_str()) );
    tempPlus.push_back( theTempPlusCopy );
    tempMinus.push_back( theTempMinusCopy );
    nSigma.push_back( theNSigma );
    // create the norm slope
    if ( (theNSigma !=0) && (hist->GetSum() != 0) ) {
      float theNormSlope = (theTempPlus->GetSum() - theTempMinus->GetSum()) / 
	( 2 *  theNSigma * hist->GetSum() );
      std::cout << "addNuisance: " << name << " " << tempPlus.size() << " theNormSlope= " << theNormSlope << std::endl;
      NormSlope.push_back( theNormSlope );
    }

//     std::cout << "addNuisance: " << name << std::endl
// 	      << "   upper " << theTempPlus->GetName() << "  " << theTempPlusCopy->GetSum() << std::endl
// 	      << "   lower " << theTempMinus->GetName() << "  " << theTempMinusCopy->GetSum() << std::endl;
    ++nNuisance;
  }
  // copy constructor
  templateRef(const templateRef& otr) : name(otr.name), inival(otr.inival),
					color (otr.color), sumW(otr.sumW), 
					nNuisance(otr.nNuisance), 
					nSigma(otr.nSigma) {
    std::cout << "Here is the templateRef copy constructor for " << otr.name << std::endl;

    std::cout << "Print the original: " << std::endl;
    otr.Print();

    //name = name + "_cp";
    hist = new TH2F( *(otr.hist) );
    hist->SetName( Form("%s_tRc", otr.hist->GetName()) );
    std::cout << "   hist was " << (otr.hist)->GetName() << std::endl;
    std::cout << "   hist copied to " << hist->GetName() << std::endl;
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      std::cout << "Make a copy of hist " << otr.tempPlus[iNuisance]->GetName() << std::endl;
      std::cout << "Size of otr.tempPlus = " << otr.tempPlus.size() << std::endl;
      std::cout << "Size of otr.tempMinus = " <<  otr.tempMinus.size() << std::endl;
      TH2F* thePlus = new TH2F( *(otr.tempPlus[iNuisance]) );
      thePlus->SetName( Form("%s_tRc", otr.tempPlus[iNuisance]->GetName()) );
      std::cout << "  " << otr.tempPlus[iNuisance]->GetName() << "  --> " 
		<< thePlus->GetName() << std::endl;
      tempPlus.push_back( thePlus );
      TH2F* theMinus = new TH2F( *(otr.tempMinus[iNuisance]) );
      theMinus->SetName( Form("%s_tRc", otr.tempMinus[iNuisance]->GetName()) );
      std::cout << "  " << otr.tempMinus[iNuisance]->GetName() << "  --> " 
		<< theMinus->GetName() << std::endl;
      tempMinus.push_back( theMinus );
    }
    std::cout << "Print the result:" << std::endl;
    this->Print();
  }
  
  static templateRef* merge(const char* theName, const char* theHistName, 
			    templateRef* tr1, templateRef* tr2, const double theInival = 1,const int theColor=-1) {
    // create a new template as merge of two existing ones
    TH2F* theHist = new TH2F( *(tr1->hist) );
    theHist->SetName( theHistName );
    theHist->SetTitle( theHistName );
    theHist->Add( tr2->hist );
    std::cout << "templateRef::merge for " << theName << std::endl;
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
  int randomize(TRandom3& rdm ) {

    // check if all nuisance parameters are in sync
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      if ( ( int(tempPlus.size()) != nNuisance) || (int(tempMinus.size()) != nNuisance) ) {
	std::cout << "TemplateRef::randomize: Mismatch in number of nuisance parameters " << std::endl;
	return -1;
      }
    }

    for (int binx=1; binx<=hist->GetXaxis()->GetNbins(); ++binx) {
      for (int biny=1; biny<=hist->GetYaxis()->GetNbins(); ++biny) {
	// apply to central histogram
	float unified = rdm.Gaus(0.,hist->GetBinError(binx,biny));
	float oldContent = hist->GetBinContent(binx,biny);
	float newContent = oldContent + unified;
// 	std::cout << "rdmTempl: binx= " << binx << " biny= " << biny
// 		  << " oldContent= " << oldContent
// 		  << " newContent= " << newContent << std::endl;
	if (newContent>0) hist->SetBinContent(binx,biny,newContent);

	// apply the same random variation to all NP
	for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
	  oldContent = tempPlus[iNuisance]->GetBinContent(binx,biny);
	  newContent = oldContent + unified;
	  if (newContent>0) tempPlus[iNuisance]->SetBinContent(binx,biny,newContent);
	  //
	  oldContent = tempMinus[iNuisance]->GetBinContent(binx,biny);
	  newContent = oldContent + unified;
	  if (newContent>0) tempMinus[iNuisance]->SetBinContent(binx,biny,newContent);
	}
      }
    }
    return 0;
  }


  void Print() const {
    std::cout << "---- Template " << name << "  base hist " << hist->GetName()
	      << " norm " << hist->GetSum() << std::endl;
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      std::cout << "      " << "Nuisance parameter " << iNuisance << std::endl;
      float varPlus = -1;
      float varMinus = -1;
      if (hist->GetSum() != 0) {
	varPlus = tempPlus[iNuisance]->GetSum() / hist->GetSum();
	varMinus = tempMinus[iNuisance]->GetSum() / hist->GetSum();
      }
      std::cout << "          " << "Upper: " << std::setw(40) << tempPlus[iNuisance]->GetName() 
		<< "  " << tempPlus[iNuisance]->GetSum() 
		<< " ratio " << varPlus << std::endl;
      std::cout << "          " << "Lower: " << std::setw(40) << tempMinus[iNuisance]->GetName() 
		<< "  " << tempMinus[iNuisance]->GetSum() 
		<< " ratio " << varMinus << std::endl;
    }
    std::cout << std::endl;
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
  std::vector<templateRef>* templates;
  std::vector<double> minConstraint;
  std::vector<double> maxConstraint;
  int fitLow;
  int fitHigh;
  double fitChi2;
  int fitNdf;
  int fitStat;
  int nNuisance;
  TMinuit* theMinuit;
  unsigned int nTemplates() { 
    return templates->size(); 
  }

  templateFitter(TH2F* theData,std::vector<templateRef>* theTemplates,int theNNuisance=0) {

    data = theData;
    templates = theTemplates;
    fitLow = -1;
    fitHigh = -1;
    // initialize the constraints for the template parameters
    for (vector<templateRef>::iterator trit=templates->begin(); trit != templates->end(); ++trit) {
      minConstraint.push_back(0);
      maxConstraint.push_back(0);
    }
    
    fitLow = 1;  // TH convention, count from 1...nbins
    fitHigh = theData->GetNbinsX();
    nNuisance = theNNuisance;
    // no constraints for nuisance parameters
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      minConstraint.push_back(0);
      maxConstraint.push_back(0);
    }
    theMinuit = NULL;
  }
  ~templateFitter() {
    //std::cout << "TemplateFitter destructor called " << std::endl;
    if (theMinuit != NULL) delete theMinuit;
    //std::cout << "TemplateFitter destructor done " << std::endl;
  }
  void SetRangeX(int theLow,int theHigh) {
    fitLow = theLow;
    fitHigh = theHigh;
  }
  void Constrain(unsigned int jpar,double lowVal,double highVal) {
    unsigned int ipar = jpar - 1;  // jpar runs from 1...npar
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
  void GetResult(int iPar, double& theValue, double& theError){
    // parameter numbering from zero
    double value;
    double error;
    int getpar = theMinuit->GetParameter(iPar,value,error);
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
  TH2F* getTemplateHist(int iTemplate) {
    if ( (iTemplate<0) || (iTemplate>=int(nTemplates())) ) {
      std::cout << "getTemplateHist: Bad template number " << iTemplate << std::endl;
      return NULL;
    }
    return (*templates)[iTemplate].hist;
  }

  int Fit(){
    // first check if all templates have the required number of nuisance parameters
    for (vector<templateRef>::iterator trit=templates->begin(); trit != templates->end(); ++trit) {
      if ( (*trit).nNuisance != this->nNuisance ) {
	std::cout << "TemplateFitter::Fit: nuisance parameter mismatch !" << std::endl;
      }
    }
    int nFitPar = nTemplates() + nNuisance;
    theMinuit = new TMinuit(nFitPar);
    theMinuit->SetFCN(fcn);
    Double_t arglist[10];
    arglist[0] = 1;
    Int_t ierflg = 0;
    // set error definition to 1 (correct for chi2 fits)
    theMinuit->mnexcm("SET ERR",arglist,1,ierflg);
    Double_t* vstart = new Double_t[nFitPar];
    Double_t* vstep = new Double_t[nFitPar];
    // for normalization of the initial fractions
    double iniTot = 0;
    // initialize the template fraction parameters
    for (vector<templateRef>::iterator trit=templates->begin(); trit != templates->end(); ++trit) {
      iniTot += (*trit).inival;
    }
    for (unsigned int i=0; i<nTemplates(); ++i) {
      vstart[i] = (*templates)[i].inival / iniTot;
      vstep[i] = vstart[i] / 100;
      theMinuit->mnparm(i,Form("par%d",i),vstart[i],vstep[i],minConstraint[i],maxConstraint[i],ierflg);
    }
    // initialize the nuisance parameters
    for (unsigned int i=nTemplates(); i<nTemplates()+nNuisance; ++i) {
      vstart[i] = 0;
      vstep[i] = 0.1;
      minConstraint[i] = 0;
      maxConstraint[i] = 0;
      theMinuit->mnparm(i,Form("nupar%d",i),vstart[i],vstep[i],minConstraint[i],maxConstraint[i],ierflg);
    }
      

//     // tell minuit to use gradient
//     arglist[0] = 1;
//     arglist[1] = 1;
//     theMinuit->mnexcm("SET GRA",arglist,2,ierflg);
#ifdef INIHESSE
    arglist[0] = 0;
    arglist[1] = 0;
    theMinuit->mnexcm("HESSE",arglist,0,ierflg);
#endif
    arglist[0] = 50000;
    arglist[1] = .1;
    theMinuit->mnexcm("MIGRAD",arglist,2,ierflg);
    // get statistics from minuit
    Double_t fmin; 
    Double_t fedm; 
    Double_t errdef; 
    Int_t npari; 
    Int_t nparx; 
    Int_t istat;
    theMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
    fitChi2 = fmin;
    fitNdf = (gTemplateFitter->fitHigh - gTemplateFitter->fitLow + 1) * 6 - npari;
    fitStat = istat;
    std::cout << " fmin=" << fmin << " fedm=" << fedm << " errdef=" << errdef
	      << " npari=" << npari << " nparx=" << nparx << " istat=" << istat << std::endl;
    std::cout << " ndf=" << fitNdf << std::endl;

//     arglist[0] = 0;
//     arglist[1] = 0;
//     theMinuit->mnexcm("HESSE",arglist,0,ierflg);


    return fitStat;
  };

  void getBGTotalError() {
    Int_t npars = theMinuit->GetNumPars();
    TMatrixDSym _corrMatrix(npars); 
    cout<<"External correlation Matrix from TMinuit"<<"\n";
    theMinuit->mnemat(_corrMatrix.GetMatrixArray(),npars );
    for (int i=0;i<npars;i++){
      for (int j=0;j<npars;j++) {
	std::cout<<"a"<<i<<j<<"="<<_corrMatrix(i,j)<<"   ";
      }
      std::cout<< std::endl;
    }

    // fractions
    double frac[100];
    double theError;
    for (int i=0; i<npars;i++){
      GetResult(i,frac[i],theError);
    }

    double sigBSq = 0;
    double B = 0;
    for (int iTemplate=0; iTemplate<int(nTemplates()); ++iTemplate) {
      B += frac[iTemplate];
    }
    for (int iTemplate=0; iTemplate<int(nTemplates()); ++iTemplate) {
      for (int jTemplate=0; jTemplate<int(nTemplates()); ++jTemplate) {
	sigBSq += ( frac[iTemplate] * frac[jTemplate] * _corrMatrix(iTemplate,jTemplate) );
      }
    }
    double sigB = -1;
    if (sigBSq >= 0.) sigB = sqrt( sigBSq );

    std::cout << " getBGTotalError:  B = " << B << "  +- " << sigB << std::endl;

    return;
  }

  TH2F* getSummedTemplateHist() {
    TH2F* theSum = NULL;

    Int_t npars = theMinuit->GetNumPars();
    TMatrixDSym _corrMatrix(npars); 
    cout<<"External correlation Matrix from TMinuit"<<"\n";
    theMinuit->mnemat(_corrMatrix.GetMatrixArray(),npars );
    for (int i=0;i<npars;i++){
      for (int j=0;j<npars;j++) {
	std::cout<<"a"<<i<<j<<"="<<_corrMatrix(i,j)<<"   ";
      }
      std::cout<< std::endl;
    }

    // fractions
    double frac[100];
    double theError;
    for (int i=0; i<npars;i++){
      GetResult(i,frac[i],theError);
    }

    for (int iTemplate=0; iTemplate<int(nTemplates()); ++iTemplate) {
      TH2F* theTemp = this->getTemplateHist(iTemplate);
      
      TH2F theScaledTemplateHist( *theTemp);
      theScaledTemplateHist.SetName( Form("%sScld",theTemp->GetName()) );
      theScaledTemplateHist.Scale( frac[iTemplate] );

      if ( theSum == NULL ) {
	theSum = new TH2F( theScaledTemplateHist );
	theSum->SetName( Form("%sSum",theTemp->GetName()) );
      } else {
	theSum->Add( &theScaledTemplateHist );
      }
    }
    return theSum;
  }

  float getTemplateContWithNuisance(int npar,double* par,
				    int iTemplate,int binx,int biny,
				    float theError) {
    // get template content of bin ix,iy
    if ((iTemplate<0) || (iTemplate>= int(nTemplates()))) {
      std::cout << " getTemplateContWithNuisance: bad template " << std::endl;
      return -999999;
    }
    if (iTemplate >=npar) {
      std::cout << " getTemplateContWithNuisance: iTemplate out of range"
		<< std::endl;
    }
    TH2F* theTemplateHist = this->getTemplateHist(iTemplate);
    if ( (binx<1) || (binx>theTemplateHist->GetXaxis()->GetNbins()) ||
	 (biny<1) || (biny>theTemplateHist->GetYaxis()->GetNbins()) ) {
      std::cout << " getTemplateContWithNuisance: bad bin numbers " 
		<<	binx << " " << biny << std::endl;
      return -999998;
    }
    double templateVal = theTemplateHist->GetBinContent(binx,biny);
    // template shifts according to nuisance parameters
    for (int iNuisance=0; iNuisance<this->nNuisance; ++iNuisance) {
      double theNuisancePar = par[ this->nTemplates() + iNuisance ];
      if (theNuisancePar>0) {
	templateVal += ( (*this->templates)[iTemplate].tempPlus[iNuisance]->GetBinContent(binx,biny) - templateVal )
	  * theNuisancePar / (*this->templates)[iTemplate].nSigma[iNuisance];
      } else {
	templateVal -= ( (*this->templates)[iTemplate].tempMinus[iNuisance]->GetBinContent(binx,biny) - templateVal )
	  * theNuisancePar / (*this->templates)[iTemplate].nSigma[iNuisance];
      }
    }
    theError = theTemplateHist->GetBinError(binx,biny);
    return templateVal;
  };

  TH2F* getTemplateHistWithNuisance(const int iTemplate) {
    if ((iTemplate<0) || (iTemplate>= int(nTemplates()))) {
      std::cout << " getTemplateHistWithNuisance: bad template " << std::endl;
      return NULL;
    }
    // copy the central template
    TH2F* theTemplateHist = new TH2F( *this->getTemplateHist(iTemplate) );
    // rename
    theTemplateHist->SetName( Form("%s_WNuis",this->getTemplateHist(iTemplate)->GetName()) );
    std::cout << "getTemplateHistWithNuisance: created template hist "
	      << theTemplateHist->GetName() << std::endl;
    // get the parameter vector
    int nPar = this->nTemplates() + nNuisance;
    double* parVec = new double[nPar];
    for (int iPar=0; iPar<nPar; ++iPar) {
      double value = 0;
      double error = 0;
      this->GetResult( iPar, value, error );
//       std::cout << "getTemplateHistWithNuisance " << iPar << " "
// 		<< value << " " << error << std::endl;
      parVec[iPar] = value;
    }
    // set contents and bins
    for (int binx=1; binx<=theTemplateHist->GetXaxis()->GetNbins(); ++binx) {
      for (int biny=1; biny<=theTemplateHist->GetYaxis()->GetNbins(); ++biny) {
	float theError = 0;
	float theContent = getTemplateContWithNuisance(nPar,parVec,iTemplate,binx,biny,
						       theError);
      }
    }
    // cleanup dynamic vector
    delete [] parVec;
	
    return theTemplateHist;
  }

  void printNPPulls(const char* pullOutFile,const char* leadText) {
    ofstream pF;
    pF.open(pullOutFile  );
    std::cout << "Pulls of nuisance parameters" << std::endl;
    double valArray[100];
    for (int iNuisance=0; iNuisance<this->nNuisance; ++iNuisance) {
      double value;
      double error;
      double vOverE = -99999;
      if (error != 0) vOverE = value / error;
      int getpar = theMinuit->GetParameter(nTemplates()+iNuisance,value,error);
      valArray[iNuisance] = value;
      //pF << Form("%8s   %4d  %12.6f  +-  %12.6f  ratio %10.4f",leadText,iNuisance,value,error,vOverE) << std::endl;
    }
    pF << Form("%8s   ",leadText);
    for (int iNuisance=0; iNuisance<this->nNuisance; ++iNuisance) {
      pF << Form(" %12.6f ",valArray[iNuisance]);
    }
    pF << std::endl;
    pF.close();
  }
};

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {    
  // this is the FCN defining the objective function
  double T = gTemplateFitter->GetData()->GetSum();
  double chisq = 0;
  // loop over the histogram bins
  for (int binx=gTemplateFitter->fitLow; binx<=gTemplateFitter->fitHigh; ++binx) {
    for (int biny=1; biny<=6; ++biny) {
      double errData = gTemplateFitter->GetData()->GetBinError(binx,biny);
      if (errData == 0) errData = 1;
      double delta = gTemplateFitter->GetData()->GetBinContent(binx,biny);
      double errMcSq = 0;
      for (int iTemplate=0; iTemplate<int(gTemplateFitter->nTemplates()); ++iTemplate) {
	TH2F* theTemplateHist = gTemplateFitter->getTemplateHist(iTemplate);
	double templateVal = theTemplateHist->GetBinContent(binx,biny);
	// template shifts according to nuisance parameters
	for (int iNuisance=0; iNuisance<gTemplateFitter->nNuisance; ++iNuisance) {
	  double theNuisancePar = par[ gTemplateFitter->nTemplates() + iNuisance ];
	  if (theNuisancePar>0) {
	    templateVal += ( (*gTemplateFitter->templates)[iTemplate].tempPlus[iNuisance]->GetBinContent(binx,biny) - templateVal )
	      * theNuisancePar / (*gTemplateFitter->templates)[iTemplate].nSigma[iNuisance];
	  } else {
	    templateVal -= ( (*gTemplateFitter->templates)[iTemplate].tempMinus[iNuisance]->GetBinContent(binx,biny) - templateVal )
	      * theNuisancePar / (*gTemplateFitter->templates)[iTemplate].nSigma[iNuisance];
	  }
	}

	delta -= (par[iTemplate] * T * templateVal);
	double errMCi = par[iTemplate] * T * theTemplateHist->GetBinError(binx,biny);
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
    double theNuisancePar = par[ gTemplateFitter->nTemplates() + iNuisance ];
    chisq += (theNuisancePar * theNuisancePar);
  }

  f = chisq;

  if (iflag == 2) {
    // compute derivative
    std::cout << "Derivative currently not operational" << std::endl;
    for (int k=0; k<npar; ++k) {
      double deriK = 0;
      for (int binx=gTemplateFitter->fitLow; binx<=gTemplateFitter->fitHigh; ++binx) {
	for (int biny=1; biny<=6; ++biny) {
	  double delta = gTemplateFitter->GetData()->GetBinContent(binx,biny);
	  double errData = gTemplateFitter->GetData()->GetBinError(binx,biny);
	  if (errData == 0) errData = 1;
	  double errSq = errData*errData;
	  
	  for (int iTemplate=0; iTemplate<int(gTemplateFitter->nTemplates()); ++iTemplate) {
	    TH2F* theTemplateHist = gTemplateFitter->getTemplateHist(iTemplate);
	    delta -= par[iTemplate] * T * theTemplateHist->GetBinContent(binx,biny);
	    double errMCi = par[iTemplate] * T * theTemplateHist->GetBinError(binx,biny);
	    errSq += errMCi * errMCi;
	  }
	  TH2F* theTemplateHistK = gTemplateFitter->getTemplateHist(k);
	  deriK += (- 2 * T * theTemplateHistK->GetBinContent(binx,biny) * delta) / errSq;
	  deriK += (- 2 * par[k] * T * T * theTemplateHistK->GetBinError(binx,biny) 
		    * theTemplateHistK->GetBinError(binx,biny)
		    * delta * delta ) / (errSq * errSq);
	}
      }
      gin[k] = deriK;
    }
  }
}


TH2F* getTrigsAndMerge(TFile* fa,char* hname,const int nTCombData,std::string tCombData[]) {
  vector<TH2F *> theData(nTCombData);
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
  std::string theName( hname );
  // strip off any path in histogram name
  long unsigned int nsl = theName.find_last_of('/');
  if (nsl != theName.npos) {
    theName = theName.substr(nsl+1,theName.size()-nsl-1);
  }
  std::string theTitle( theData[0]->GetTitle() );
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

float makeSquare(const float x) { return (x*x); }

double getSystBbpurTpat( const std::string theScenario, const int theSignalMass ) {
  const int mSignal = 11;
  int masses[mSignal] = {90, 100, 120, 130, 140, 160, 180, 200, 250, 300, 350};
  double systLowMass[mSignal] = { 0.001154123,
				  0.004609176,	
				  0.002284732,	
				  0.001883083,	
				  0.001425553,	
				  0.00056921,	
				  0.00033541,	
				  0.000618466,
				  0.000492037,
                                  0.0000,
				  0.000212132};

  double systMediumMass[mSignal] = { 0.006302071,
				     0.002675892,
				     0.003208878,
				     0.002974693,
				     0.002928652,
				     0.001868689,
				     0.000524786,
				     0.00068542,
				     0.000764853,
                                     0.0000,
				     0.000212132 };

  // find the mass bin
  int iMass = -1;
  for (int i=0; i<mSignal; ++i) {
    if ( masses[i] == theSignalMass ) {
      iMass = i;
      break;
    }
  }
  if (iMass <0) {
    std::cout << "getSystBbpurTpat: bad signal mass " << theSignalMass << std::endl;
    return -1.;
  }

  if (theScenario == "LowMass2011") {
    return systLowMass[iMass];
  } else if (theScenario == "MediumMass2011") {
    return systMediumMass[iMass];
  } else {
    std::cout << "getSystBbpurTpat: bad scenario " << theScenario << std::endl;
    return -1.;
  }
}


int getSignalStats(const int theSignalMass) {
  const int mSignal = 11;
  int masses[mSignal] = {90, 100, 120, 130, 140, 160, 180, 200, 250, 300, 350};
  int signalStats[mSignal] = { 1088786, 1089696, 1093650,1061374, 550000,
			       550000, 550000, 549478, 550000, 550000, 550000 };
  // find the mass bin
  int iMass = -1;
  for (int i=0; i<mSignal; ++i) {
    if ( masses[i] == theSignalMass ) {
      iMass = i;
      break;
    }
  }
  if (iMass <0) {
    std::cout << "getSignalStats: bad signal mass " << theSignalMass << std::endl;
    return -1.;
  }
  return signalStats[iMass];
}


int getScale(const int theSignalMass,double& sFrac,double& sSeen,double& sXSect, double& sTb) {
  const int mSignal = 11;
  int masses[mSignal] = {90, 100, 120, 130, 140, 160, 180, 200, 250, 300, 350};
  // here we need only typical numbers to set the limits of the histograms
  double Fracs[mSignal] = { 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02 };
  double Seens[mSignal] = {8000., 4000., 1200., 1000., 1000., 800., 700., 700., 400., 400., 300.};
  double XSects[mSignal] = {2000., 2000., 200., 100., 80., 50., 60., 40., 20., 10., 10.};
  double Tbs[mSignal] = { 25., 25., 25., 25., 25., 25., 40., 40., 40., 40., 40. };

  // find the mass bin
  int iMass = -1;
  for (int i=0; i<mSignal; ++i) {
    if ( masses[i] == theSignalMass ) {
      iMass = i;
      break;
    }
  }
  if (iMass <0) {
    std::cout << "getSignalStats: bad signal mass " << theSignalMass << std::endl;
    return -1.;
  }
  sFrac = Fracs[iMass];
  sSeen = Seens[iMass];
  sXSect = XSects[iMass];
  sTb = Tbs[iMass];
  return 0;
}




void customFitS_Mass(const int iSignal,const bool fitSignal) {
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

  std::string sfc[nfc] = { "q", "c", "b" };

  //std::string scenario("LowMass2011");
  //std::string scenario("MediumMass2011");
  std::string scenario;
  bool useTemplateError;
  bool useNP;
  if (  getHbbCfg(scenario,useTemplateError,useNP) != 0 ) return;

//   const int nSignal=7;
//   int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };

#include "Analysis/Utilities/interface/HbbMass.h"

  if (iSignal >= nSignal) {
    std::cout << "Bad iMass=" << iSignal << std::endl;
    return;
  }


//   const int nSyst = 3;
//   std::string systName[nSyst] = { "JES", "SFbc", "SFudsg" };
  const int nSyst = 4;
  std::string systName[nSyst] = { "JES", "SFbc", "SFudsg", "JER" };

  bool activeNuisance[nSyst] = { true, true, true, true };
  //bool activeNuisance[nSyst] = { false, false, false };

  double efficiency[nSignal][nbtag];
  bool signalNuisanceCentral = true;

  bool includeSys = false;

  bool drawTemplateErrors = true;

  bool addBbpurTpat = true;
 
  double intLumi = 0;
  if (scenario == "LowMass2011") {
    //intLumi = 2.66794; // in fb-1
    intLumi = 2.692643;  // with new method
  } else if (scenario == "MediumMass2011") {
    //intLumi = 3.99983; // in fb-1
    intLumi = 4.040802;
  } else {
    std::cout << "Bad scenario" << std::endl;
    return;
  }

  // mass label (for plots)
  std::string cMass;
  if (fitSignal) {
    cMass.assign( Form("mH_%04d",signalMass[iSignal]) );
  } else {
    cMass.assign( "BGonly" );
  }

  // summary variables across the masses
  double fPenalty[nbtag][nSignal];
  double HFraction[nbtag][nSignal];
  double HFractionError[nbtag][nSignal];
  double HSeen[nbtag][nSignal];
  double HSeenError[nbtag][nSignal];
  double nominalXSect[nbtag][nSignal];
  double xSectUncUp[nbtag][nSignal];
  double xSectUncDown[nbtag][nSignal];
  double tanBUncUp[nbtag][nSignal];
  double tanBUncDown[nbtag][nSignal];
  double HEfficiency[nbtag][nSignal];
  double HCorr[nbtag][nSignal];
  double HCorrError[nbtag][nSignal];
  double HXSect[nbtag][nSignal];
  double HXSectError[nbtag][nSignal];
  double HXSectUpLimit[nbtag][nSignal];
  double tanBetaUpLimit[nbtag][nSignal];

  // and the syst limits
  double HXSectErrorSys[nbtag][nSignal];
  double HXSectUpLimitSys[nbtag][nSignal];
  double tanBetaUpLimitSys[nbtag][nSignal];
    


  // open the packed template file
  std::string ptFileName( Form("packedTemplates-M-%d.root",signalMass[iSignal]) );
  std::cout << "Open packed template file " << ptFileName << std::endl;
  TFile* ptFile = new TFile( ptFileName.c_str() );
  if (ptFile == NULL) {
    std::cout << "Error opening " << ptFileName << std::endl;
    return;
  }

  // arrays to hold the hist-to-be-fitted, the signal & background templates
  const int nUpDown = 2;
  std::string upDownName[nUpDown] = { "Up", "Down" };
  TH2F* hSignalCentral[nbtag];
  TH2F* hSignalSyst[nSyst][nUpDown][nbtag];

  TH2F* hBackgroundCentral[nbtag][ncateg][nfc];
  TH2F* hBackgroundSyst[nSyst][nUpDown][nbtag][ncateg][nfc];

  TH2F* mjjEbtdata[nbtag];
  TH2F* mjjEbtdataRdm[nbtag]; // array for randomized hist-to-be-fitted

  std::string tFlav[3] = {"Uds", "C", "B"};
  std::string sFlav[3] = {"Q", "C", "B"};

  // output file
  TFile* hout = new TFile(Form("customFitS-%s.root",cMass.c_str()),"recreate");
  hout->cd();
  TH2::AddDirectory(true);

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    // get the signal central template
    hSignalCentral[ibtag] = (TH2F*) ptFile->Get( Form("bbH_%s",sbtag[ibtag].c_str()) );
    if (hSignalCentral[ibtag] == NULL)  {
      std::cout << "Reading of SignalCentral failed" << std::endl;
      return;
    }
    std::cout << "Read " << hSignalCentral[ibtag]->GetName() << std::endl;
    std::cout << hSignalCentral[ibtag]->GetName() << " TotalContents=" << hSignalCentral[ibtag]->GetSum()
	      << std::endl;

    // read the efficiency
    TH1F* histEffMerged = (TH1F*) ptFile->Get(Form("EffMerged%s",sbtag[ibtag].c_str()));
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
    TH1F* histXSect = (TH1F*) ptFile->Get("xsect");
    if ( histXSect == NULL) {
      std::cout << "xsect" << " not found" << std::endl;
      return;
    }
//    nominalXSect[ibtag][iSignal] = histXSect->GetBinContent(1);
    nominalXSect[ibtag][iSignal] = hbb_xs.getXsecBR(signalMass[iSignal],tanBeta);
//     if (signalMass[iSignal] == 130) {
//       nominalXSect[ibtag][iSignal] = 93.44;
//       std::cout << "fix triple-merge problem for mH=" << signalMass[iSignal] << std::endl
// 		<< " Replace Xsection of " << histXSect->GetBinContent(1)
// 		<< " by " << nominalXSect[ibtag][iSignal] << std::endl;
//     }
    
    double muup    = hbb_xs.getXsecBRUnc_muup(signalMass[iSignal],tanBeta);
    double pdfup   = hbb_xs.getXsecBRUnc_pdfalphas68up(signalMass[iSignal],tanBeta);
    double mudown  = hbb_xs.getXsecBRUnc_mudown(signalMass[iSignal],tanBeta);
    double pdfdown = hbb_xs.getXsecBRUnc_pdfalphas68down(signalMass[iSignal],tanBeta);
    double ueps    = hbb_xs.getXsecBRUnc_UEPS(signalMass[iSignal],tanBeta);
    
    xSectUncUp[ibtag][iSignal]   = sqrt(muup*muup + pdfup*pdfup + ueps*ueps );
    xSectUncDown[ibtag][iSignal] = sqrt(mudown*mudown + pdfdown*pdfdown + ueps*ueps );
    
    // tanBetaError = 1/2 tanBeta * xSecError/xSec
    tanBUncUp[ibtag][iSignal]   = 0.5 * tanBeta * xSectUncUp[ibtag][iSignal]/nominalXSect[ibtag][iSignal];
    tanBUncDown[ibtag][iSignal] = 0.5 * tanBeta * xSectUncDown[ibtag][iSignal]/nominalXSect[ibtag][iSignal];
    
    
    std::cout << "Tan Beta "  << tanBeta << std::endl;
   std::cout << "Total xsec*BR                   = " << hbb_xs.getXsecBR(signalMass[iSignal],tanBeta)                    << std::endl;
   std::cout << "Total xsec*BR Unc MuUP          = " << hbb_xs.getXsecBRUnc_muup(signalMass[iSignal],tanBeta)            << std::endl;
   std::cout << "Total xsec*BR Unc MuDW          = " << hbb_xs.getXsecBRUnc_mudown(signalMass[iSignal],tanBeta)          << std::endl;
   std::cout << "Total xsec*BR Unc PdfAlphas68UP = " << hbb_xs.getXsecBRUnc_pdfalphas68up(signalMass[iSignal],tanBeta)   << std::endl;
   std::cout << "Total xsec*BR Unc PdfAlphas68DW = " << hbb_xs.getXsecBRUnc_pdfalphas68down(signalMass[iSignal],tanBeta) << std::endl;
   std::cout << "Total xsec*BR Unc UE-PS         = " << hbb_xs.getXsecBRUnc_UEPS(signalMass[iSignal],tanBeta)            << "    -  This is simply the xsec*BR*0.04. " << std::endl;
   


//     // create empty file just as marker
//     ofstream markerFile;
//     markerFile.open(Form("pack-%s-%s.txt",sbtag[ibtag].c_str(),scenario.c_str()),ios::app);
//     markerFile << "Template for mass " << signalMass[iSignal] << std::endl;
//     markerFile.close();

    // read the signal syst templates
    for (int iSyst=0; iSyst<nSyst; ++iSyst) {
      for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	hSignalSyst[iSyst][iUpDown][ibtag] 
	  = (TH2F*) ptFile->Get( Form("bbH_%s_%s_%s",systName[iSyst].c_str(),
				      upDownName[iUpDown].c_str(),sbtag[ibtag].c_str()) );
	if (hSignalSyst[iSyst][iUpDown][ibtag] == NULL) {
	  std::cout << "Reading of SignalSyst failed" << std::endl;
	  return;
	}
	std::cout << "Read " << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << std::endl;
      }
    }

    // read the background central templates
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	templateId tName(ifc,icateg);
	std::string bgCentralTemplateName( Form("%s_%s",tName.name().c_str(),sbtag[ibtag].c_str()) );
	hBackgroundCentral[ibtag][icateg][ifc] = 
	  (TH2F*) ptFile->Get( bgCentralTemplateName.c_str() );
	if ( hBackgroundCentral[ibtag][icateg][ifc] == NULL ) {
	  std::cout << "Hist not found: " << bgCentralTemplateName << std::endl;
	  return;
	}
	std::cout << "Read " << hBackgroundCentral[ibtag][icateg][ifc]->GetName() << std::endl;
      }
    }

    // read the background syst templates
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	for (int iSyst=0; iSyst<nSyst; ++iSyst) {
	  for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	    templateId tName(ifc,icateg);
	    std::string bgSystTemplateName( Form("%s_%s_%s_%s",tName.name().c_str(),systName[iSyst].c_str(),
					    upDownName[iUpDown].c_str(),sbtag[ibtag].c_str()) );
	    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]
	      = (TH2F*) ptFile->Get( bgSystTemplateName.c_str() );
	    if ( hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] == NULL ) {
	      std::cout << "Hist not found: " << bgSystTemplateName << std::endl;
	      return;
	    }
	    std::cout << "Read " << hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName()
		      << std::endl;
	  }
	}
      }
    }

    // read the hist-to-be-fitted
    std::string dataName( Form("Data_%s",sbtag[ibtag].c_str()) );
    mjjEbtdata[ibtag] = (TH2F*) ptFile->Get( dataName.c_str() );
    if (mjjEbtdata[ibtag] == NULL) {
      std::cout << "Histogram not found: " << dataName << std::endl;
      return;
    }
    std::cout << "Read " << mjjEbtdata[ibtag]->GetName() << std::endl;


    

  }

  // mass range bins
  int lowBinX = 3;
  int highBinX = 22;

  TH2F* mjjEbtTemplate[nfc][nbtag][ncateg];
  // create synonyms for background templates
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	mjjEbtTemplate[ifc][ibtag][icateg] = hBackgroundCentral[ibtag][icateg][ifc];
	std::cout << "ifc= " << ifc << " ibtag= " << ibtag
		  << " icateg= " << icateg << " mapped from "
		  << hBackgroundCentral[ibtag][icateg][ifc]->GetName() << std::endl;
      }
    }
  }
  // create synonyms for signal templates
  TH2F* mjjEbtsignal[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtsignal[ibtag] = hSignalCentral[ibtag];
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

  // create projections for later use
  TH1D* mjjEbtdataProX[nbtag];
  TH1D* mjjEbtdataProY[nbtag];
  // read the hist-to-be-fitted
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
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

  // set drawing colors
  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	mjjEbtTemplate[ifc][ibtag][icateg]->SetLineColor( tc[icateg+1+3*ifc] );
	mjjEbtTemplate[ifc][ibtag][icateg]->SetLineWidth(3);  
	
      }
    }
  }

  // build the set of template histograms

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

  // build the signal templates
  templateRef* mjjEbtSignalTemplateRef[nbtag];
  if (fitSignal) {
    for (int ibtag = 0; ibtag < nbtag; ++ibtag) {
      mjjEbtSignalTemplateRef[ibtag] = new templateRef("H",mjjEbtsignal[ibtag],0.001,2);
    }

    // add the nuisance parameters to the signal template
    for (int iNuisance=0; iNuisance<nSyst; ++iNuisance) {
      if ( activeNuisance[iNuisance] ) {
	for (int ibtag = 0; ibtag < nbtag; ++ibtag) {
	  std::cout << "Add signal nusances " << iNuisance << " " 
		    << " ibtag " << ibtag << systName[0] << std::endl;
	  std::cout << "        pointer " << hSignalSyst[iNuisance][0][ibtag] << std::endl;
	  std::cout << "        pointer " << hSignalSyst[iNuisance][1][ibtag] << std::endl;
	  std::cout << "Add signal nuisance " << hSignalSyst[iNuisance][0][ibtag]->GetName() << std::endl;
	  std::cout << "Add signal nuisance " << hSignalSyst[iNuisance][1][ibtag]->GetName() << std::endl;
	  if ( (hSignalSyst[iNuisance][0][ibtag] != NULL) && (hSignalSyst[iNuisance][1][ibtag] != NULL) ) {	  
	    mjjEbtSignalTemplateRef[ibtag]->addNuisance( hSignalSyst[iNuisance][0][ibtag],hSignalSyst[iNuisance][1][ibtag], 2.);
	  } else {
	    std::cout << "Signal syst template missing " << std::endl;
	  }
	}
      }
    }
    

    std::cout  << std::endl
	       << "===========================================" << std::endl
	       << "=                                         =" << std::endl
	       << "=   Signal     templates                  =" << std::endl
	       << "=                                         =" << std::endl
	       << "===========================================" << std::endl << std::endl;
    for (int ibtag = 0; ibtag < nbtag; ++ibtag) {
      mjjEbtSignalTemplateRef[ibtag]->Print();
    }
  }


  // add the nuisance parameters to the background templates
  for (int iNuisance=0; iNuisance<nSyst; ++iNuisance) {
    if ( activeNuisance[iNuisance] ) {
      for (int ibtag = 0; ibtag < nbtag; ++ibtag) {
	for (int icateg=0; icateg<ncateg; ++icateg) {
	  for (int ifc=0; ifc<nfc; ++ifc) {
	    if ( (hBackgroundSyst[iNuisance][0][ibtag][icateg][ifc] != NULL) && (hBackgroundSyst[iNuisance][1][ibtag][icateg][ifc] != NULL) ) {
	      mjjEbtTemplateRef[ifc][ibtag][icateg]->addNuisance(hBackgroundSyst[iNuisance][0][ibtag][icateg][ifc],hBackgroundSyst[iNuisance][1][ibtag][icateg][ifc],2.);
	      std::cout << "Add BG nuisance " << iNuisance << "  " <<  hBackgroundSyst[iNuisance][0][ibtag][icateg][ifc]->GetName()
			<< "  " << hBackgroundSyst[iNuisance][1][ibtag][icateg][ifc]->GetName() << std::endl;
	    } else {
	      std::cout << "hBackgroundSyst missing " << std::endl;
	      return;
	    }
	  }
	}
      }
    }
  }

  std::cout  << std::endl
	     << "===========================================" << std::endl
	     << "=                                         =" << std::endl
	     << "=   Background templates                  =" << std::endl
	     << "=                                         =" << std::endl
	     << "===========================================" << std::endl << std::endl;
  for (int ibtag = 0; ibtag < nbtag; ++ibtag) {
    for (int ifc=0; ifc<nfc; ++ifc) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
// 	templateId tName(ifc,icateg);
// 	std::cout << "Print BG template ifc=" << ifc << " icateg=" << icateg
// 		  << "  " << tName.name() << std::endl;
	mjjEbtTemplateRef[ifc][ibtag][icateg]->Print();
      }
    }
  }

  for (int ibtag = 0; ibtag < nbtag; ++ibtag) {

    // produce templates merging the first two categories
    templateRef* mjjEbtTemplateMerge01[nfc][nbtag];
    for (int ifc=0; ifc<nfc; ++ifc) {
      // templates to be merged
      std::string mergedName( Form("(%sb)b",sFlav[ifc].c_str()) );
      mjjEbtTemplateMerge01[ifc][ibtag] = 
	templateRef::merge(mergedName.c_str(),
			   Form("%s_%s",mergedName.c_str(),sbtag[ibtag].c_str()),
			   mjjEbtTemplateRef[ifc][ibtag][0],mjjEbtTemplateRef[ifc][ibtag][1],0.3);
      mjjEbtTemplateMerge01[ifc][ibtag]->Print();
    }
    // produce templates merging bbC and bbQ
    templateRef* mjjEbtTemplateMerge_bbCQ[nbtag];
    mjjEbtTemplateMerge_bbCQ[ibtag] = 
      templateRef::merge( "bbX",
			  Form("bbX_%s",sbtag[ibtag].c_str()),
			  mjjEbtTemplateRef[0][ibtag][2],mjjEbtTemplateRef[1][ibtag][2],0.1);
    
    mjjEbtTemplateMerge_bbCQ[ibtag]->Print();




    //
    // here the looper starts 
    //

    // get the typical scales
    double sFrac;
    double sSeen;
    double sCorr;
    double sXSect;
    double sTb;
    if ( getScale(signalMass[iSignal],sFrac,sSeen,sXSect,sTb) <0) {
      std::cout << "getScale fails for mass " << signalMass[iSignal] << std::endl;
      return;
    }
    sCorr = sSeen / efficiency[iSignal][ibtag];

    // make the histograms
    int nBins = 400;
    TH1D* FitStatus_H = new TH1D("FitStatus_H","Fit Status",10,-0.5,9.5);
    TH1D* HFraction_H = new TH1D("HFraction_H","Higgs Fraction",nBins,-5*sFrac,+5*sFrac);
    TH1D* HSeen_H = new TH1D("HSeen_H","Higgs Seen",nBins,-5*sSeen,+5*sSeen);
    TH1D* HCorr_H = new TH1D("HCorr_H","Higgs Corr",nBins,-5*sCorr,+5*sCorr);
    TH1D* HXSect_H = new TH1D("HXSect_H","Higgs XSect",nBins,-5*sXSect,+5*sXSect);

    TH1D* HFractionPull_H = new TH1D("HFractionPull_H","Higgs Fraction Pull",nBins,-5,5);
    //TH1D* HSeenPull_H = new TH1D("HSeenPull_H","Higgs Seen",nBins,-3,3);
    //TH1D* HCorrPull_H = new TH1D("HCorrPull_H","Higgs Corr",nBins,-3,3);
    TH1D* HXSectPull_H = new TH1D("HXSectPull_H","Higgs XSect Pull",nBins,-5,5);

    //TH1D* HSeenUL95_H = new TH1D("HSeenUL95_H","HSeenUL95",nBins,0,5*sSeen);
    TH1D* HXSectUL95_H = new TH1D("HXSectUL95_H","Higgs XSectUL95",nBins,0,5*sXSect);
    TH1D* HtanBetaUL95_H = new TH1D("HtanBetaUL95_H","Higgs tanBetaUL95",nBins,0,5*sTb);

    TRandom3 rdm;  // set up the random generator
    rdm.SetSeed(0);

    int nLoop = 5000;

    for (int iLoop=0; iLoop<nLoop; ++iLoop) {

      std::cout << "Start loop for iLoop=" << iLoop << std::endl;

      // set up the templates configuration for this fit
      vector<templateRef> templates;
      templates.reserve(10);
      //std::cout << "Push " << mjjEbtTemplateMerge01[2][ibtag]->name << std::endl;
      
      //mjjEbtTemplateMerge01[2][ibtag]->Print();

//       std::cout << "Just exercise the copy constructor" << std::endl;
//       templateRef testRef( * mjjEbtTemplateMerge01[2][ibtag] );
//       testRef.name = "Exercise";
//       std::cout << "Print the copy" << std::endl;
//       testRef.Print();

      mjjEbtTemplateMerge01[2][ibtag]->Print();
      templates.push_back( *mjjEbtTemplateMerge01[2][ibtag] );  // (Bb)b
      //templates.back().Print();
      //std::cout << "Push " << mjjEbtTemplateMerge01[1][ibtag]->name << std::endl;
      templates.push_back( *mjjEbtTemplateMerge01[1][ibtag] );  // (Cb)b
      //std::cout << "Push " << mjjEbtTemplateMerge01[0][ibtag]->name << std::endl;
      templates.push_back( *mjjEbtTemplateMerge01[0][ibtag] );  // (Qb)b
      //std::cout << "Push " << mjjEbtTemplateRef[2][ibtag][2]->name << std::endl;
      templates.push_back( *mjjEbtTemplateRef[2][ibtag][2] ); // bbB
      //std::cout << "Push " << mjjEbtTemplateMerge_bbCQ[ibtag]->name << std::endl;
      templates.push_back( *mjjEbtTemplateMerge_bbCQ[ibtag] ); // bbX

      //std::cout << "Looper: after building templates" << std::endl;

//       std::cout << "Master templates before randomize :" << std::endl;
//       std::cout << mjjEbtTemplateMerge01[1][ibtag]->hist->GetSum() << " "
// 	<< mjjEbtTemplateMerge01[1][ibtag]->hist->GetMean() << " "
// 		<< mjjEbtTemplateMerge01[1][ibtag]->hist->GetRMS() << std::endl;

      // randomize the background templates in the templates vector
      for (unsigned int iTemplate=0; iTemplate<templates.size(); ++iTemplate) {
	int rdm_stat = templates[iTemplate].randomize( rdm );
	if (rdm_stat < 0) {
	  std::cout << "Error from Randomize" << std::endl;
	  return;
	}
      }

//       std::cout << "Master templates after randomize :" << std::endl;
//       std::cout << mjjEbtTemplateMerge01[1][ibtag]->hist->GetSum() << " "
// 	<< mjjEbtTemplateMerge01[1][ibtag]->hist->GetMean() << " "
// 		<< mjjEbtTemplateMerge01[1][ibtag]->hist->GetRMS() << std::endl;

      
      // add the signal template
      if (fitSignal) {
	if (mjjEbtsignal[ibtag] != NULL) {
	  templates.push_back( *mjjEbtSignalTemplateRef[ibtag] );
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

      // print the templates before normalization
      std::cout << "Before normalization " << std::endl;
      for (vector<templateRef>::iterator trit=templates.begin(); trit != templates.end(); ++trit) {
	(*trit).Print();
      }

      // normalize the templates
      std::cout << "Looper: normalize " << std::endl;
      for (vector<templateRef>::iterator trit=templates.begin(); trit != templates.end(); ++trit) {
	if ( ! signalNuisanceCentral ) {
	  (*trit).normalize(0);
	} else {
	  if ((*trit).name == "H") {
	    (*trit).normalize(1);
	    std::cout << "Normalize central for signal" << std::endl;
	  } else {
	    (*trit).normalize(0);
	  }
	}
      }

      // print the templates after normalization
      std::cout << "Looper: after normalization " << std::endl;
      for (vector<templateRef>::iterator trit=templates.begin(); trit != templates.end(); ++trit) {
	(*trit).Print();
      }

      // count number of active nuisance parameters
      int nActiveNuisance = 0;
      for (int iNuisance=0; iNuisance<nSyst; ++iNuisance) {
	if ( activeNuisance[iNuisance] ) ++nActiveNuisance;
      }

      std::cout << "Number of active nuisance parameters: " << nActiveNuisance << std::endl
		<< "  " << activeNuisance[0]
		<< "  " << activeNuisance[1]
		<< "  " << activeNuisance[2]
		<< std::endl;

      // randomize the hist-to-be-fitted
      mjjEbtdataRdm[ibtag] = new TH2F( *mjjEbtdata[ibtag] );
      mjjEbtdataRdm[ibtag]->SetName( Form("%s_Rdm_%d",mjjEbtdata[ibtag]->GetName(),iLoop) );
      // Poissonian smearing
      for (int binx=1; binx<=mjjEbtdataRdm[ibtag]->GetXaxis()->GetNbins(); ++binx) {
	for (int biny=1; biny<=mjjEbtdataRdm[ibtag]->GetYaxis()->GetNbins(); ++biny) {
	  float newContent = rdm.Poisson( mjjEbtdata[ibtag]->GetBinContent(binx,biny) );
	  mjjEbtdataRdm[ibtag]->SetBinContent(binx,biny,newContent);
	}
      }

      std::cout << "Looper: build gTemplateFitter" << std::endl;
      gTemplateFitter = new templateFitter(mjjEbtdataRdm[ibtag], &templates,  nActiveNuisance );
    
      for (int iTemplate=0; iTemplate< int(gTemplateFitter->nTemplates()); ++iTemplate) {
	if ((*gTemplateFitter->templates)[iTemplate].name != "H" || constrainSignal) {
	  gTemplateFitter->Constrain(iTemplate+1,0.0,2*gTemplateFitter->data->GetSum() / 
				     gTemplateFitter->getTemplateHist(iTemplate)->GetSum());
	}
      }
      gTemplateFitter->SetRangeX(lowBinX,highBinX);
      std::cout << std::endl
		<< "   Fit for signal:  " << signalMass[iSignal] << std::endl << std::endl; 

      Int_t status = gTemplateFitter->Fit(); 
      cout << sbtag[ibtag] << " fit status: " << status << endl;
      //gTemplateFitter->printNPPulls(Form("NuisanceParPulls-%s.txt",cMass.c_str()),cMass.c_str());
      //gTemplateFitter->getBGTotalError();

      mjjEbtdata[ibtag]->SetMarkerStyle(20);
      mjjEbtdata[ibtag]->GetXaxis()->SetTitle("m(Jet1 Jet2) [GeV/c^{2}]");
      mjjEbtdata[ibtag]->GetYaxis()->SetTitle("EvtBTag");
      mjjEbtdata[ibtag]->GetZaxis()->SetTitle("N / 10 GeV/c^{2}");
      mjjEbtdata[ibtag]->ProjectionX()->Draw("Ep");
      double value = 0;
      double error = 0;
      double fPenalty = -1;

      int iH = -1;  // pointer to Higgs template
      for (int ii=0; ii<int(gTemplateFitter->nTemplates()); ++ii) {
	if (templates[ii].name == "H") {
	  iH = ii;
	}
      }

      // print Higgs statistics (if H template was fitted)
      if (iH>=0) {
	gTemplateFitter->GetResult(iH, HFraction[ibtag][iSignal], HFractionError[ibtag][iSignal]);
	// determine the fitting penalty
	fPenalty =  HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSum() / sqrt( mjjEbtdata[ibtag]->GetSum() );

	// add systematic errors
	if (addBbpurTpat) {
	  double theBbpurTpat = getSystBbpurTpat( scenario.c_str(), signalMass[iSignal] );
	  if (theBbpurTpat < 0) {
	    std::cout << "Error from getSystBbpurTpat" << std::endl;
	    return;
	  }
	  HFractionError[ibtag][iSignal] = 
	    sqrt( HFractionError[ibtag][iSignal] * HFractionError[ibtag][iSignal] + theBbpurTpat * theBbpurTpat );
	}

	HSeen[ibtag][iSignal] = HFraction[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSum();
	HSeenError[ibtag][iSignal] = HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSum();

	HEfficiency[ibtag][iSignal] = efficiency[iSignal][ibtag];

	// apply the systematics
	double deltaSigmaSysSq = 0;
      
	deltaSigmaSysSq += makeSquare(0.022); // luminosity
	deltaSigmaSysSq += makeSquare(0.10); // trigger
	deltaSigmaSysSq += makeSquare(0.3225); // online btag
	deltaSigmaSysSq += makeSquare(0.06); // btag efficiency topology dependence
	// statistics of signal sample
	int theSignalStatistics = getSignalStats(signalMass[iSignal]);
	if (theSignalStatistics < 0) {
	  std::cout << "Error from getSignalStats" << std::endl;
	  return;
	}
	float sysSignalStatistics = sqrt( (1 - efficiency[iSignal][ibtag]) 
					  / (efficiency[iSignal][ibtag] * theSignalStatistics) );
	//std::cout << " Signal stats rel error: " << sysSignalStatistics << std::endl;
	deltaSigmaSysSq += makeSquare(sysSignalStatistics);

	for (int iNuisance=0; iNuisance<mjjEbtSignalTemplateRef[ibtag]->nNuisance; ++iNuisance) {
	  float nuisNorm = mjjEbtSignalTemplateRef[ibtag]->getNormSlope(iNuisance);
	  std::cout<< " iTemplate=" << iNuisance << " nuisNorm=" << nuisNorm << std::endl;
	  deltaSigmaSysSq += makeSquare(nuisNorm); // SFbc effect on signal
	}
	//std::cout << "Relative deltaSigmaSys=" << sqrt( deltaSigmaSysSq ) << std::endl;

	HCorr[ibtag][iSignal] = HSeen[ibtag][iSignal] / HEfficiency[ibtag][iSignal];
	HCorrError[ibtag][iSignal] = HSeenError[ibtag][iSignal] / HEfficiency[ibtag][iSignal];
	HXSect[ibtag][iSignal] = 0.001 * HCorr[ibtag][iSignal] / intLumi;
	HXSectError[ibtag][iSignal] = 0.001 * HCorrError[ibtag][iSignal] / intLumi;
	
	
	HXSectErrorSys[ibtag][iSignal] = sqrt( HXSectError[ibtag][iSignal]*HXSectError[ibtag][iSignal]
					     + HXSect[ibtag][iSignal]*HXSect[ibtag][iSignal]*deltaSigmaSysSq );
	//HXSectError[ibtag][iSignal] * sqrt(1 + deltaSigmaSysSq); 

	HXSectUpLimit[ibtag][iSignal] = uplimit( HXSect[ibtag][iSignal], HXSectError[ibtag][iSignal]);

	HXSectUpLimitSys[ibtag][iSignal] = uplimit( HXSect[ibtag][iSignal], HXSectErrorSys[ibtag][iSignal] );

	tanBetaUpLimit[ibtag][iSignal] = sqrt( HXSectUpLimit[ibtag][iSignal] / nominalXSect[ibtag][iSignal] ) * 20;
   
   // tan beta aqui (find error in tan beta)
   
	std::cout << "Efficiency = " << HEfficiency[ibtag][iSignal] << std::endl
		  << "Higgs observed = " << HSeen[ibtag][iSignal] << " +- " << HSeenError[ibtag][iSignal] << std::endl
		  << "Higgs corrected = " << HCorr[ibtag][iSignal] << " +- " << HCorrError[ibtag][iSignal] << std::endl
		  << "Higgs XSection = " << HXSect[ibtag][iSignal] << " +- " << HXSectError[ibtag][iSignal] << " pb " << std::endl
		  << "Higgs XSection < " << HXSectUpLimit[ibtag][iSignal] << " pb " << std::endl
		  << "tanBeta        < " << tanBetaUpLimit[ibtag][iSignal] << std::endl
		  << "Assumed IntLumi = " << intLumi << "fb-1" << std::endl;

	std::cout << "Higgs XSection Sys                         = " << HXSectErrorSys[ibtag][iSignal] << " pb " << std::endl;
	std::cout << "Higgs XSection with Sys < " << HXSectUpLimitSys[ibtag][iSignal] << " pb " << std::endl;

	if (includeSys) {
	  std::cout << "Including syst errors " << std::endl;
	  HXSectError[ibtag][iSignal] = HXSectErrorSys[ibtag][iSignal];
	  HXSectUpLimit[ibtag][iSignal] = HXSectUpLimitSys[ibtag][iSignal];
	}

	// fill the histograms
	FitStatus_H->Fill(double(status));
	if (status == 3) {
	  HFraction_H->Fill( HFraction[ibtag][iSignal] );
	  HSeen_H->Fill( HSeen[ibtag][iSignal] );
	  HCorr_H->Fill( HCorr[ibtag][iSignal] );
	  HXSect_H->Fill( HXSect[ibtag][iSignal] );
	  HXSectUL95_H->Fill( HXSectUpLimit[ibtag][iSignal] );
	  HtanBetaUL95_H->Fill( tanBetaUpLimit[ibtag][iSignal] );
	
	  // fill the pulls 
	  if ( HFractionError[ibtag][iSignal] != 0 ) 
	    HFractionPull_H->Fill( HFraction[ibtag][iSignal] / HFractionError[ibtag][iSignal] );
	  if ( HXSectError[ibtag][iSignal] != 0 )
	    HXSectPull_H->Fill( HXSect[ibtag][iSignal] / HXSectError[ibtag][iSignal] );
	}

#ifdef SUMMARYCARD
	// write the summary card to a file
	ofstream headCard;
	headCard.open("sumCard-Header.txt");
	headCard << Form("| %4s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %5s |",
			 "Mass","Eff","chi2","frac","Efrac","XSect","EXSect","XSectUL95","tbUL95","theoXS","fstat")
		 << std::endl;
	ofstream sumCard;
	sumCard.open(Form("sumCard-M%04d-%s-%s.txt",signalMass[iSignal],
			  sbtag[ibtag].c_str(),scenario.c_str()));
	sumCard << Form("| %4d | %10.7f | %10.2f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %10.4f | %5d |",
			signalMass[iSignal],efficiency[iSignal][ibtag],
			gTemplateFitter->GetChisquare(),
			HFraction[ibtag][iSignal],HFractionError[ibtag][iSignal],
			HXSect[ibtag][iSignal],HXSectError[ibtag][iSignal],
			HXSectUpLimit[ibtag][iSignal],tanBetaUpLimit[ibtag][iSignal],
			nominalXSect[ibtag][iSignal],status)
		<< std::endl;
	sumCard.close();
	headCard.close();
#endif
      }
    
      // delete the TemplateFitter object
      std::cout << "End of looper: delete gTemplateFitter" << std::endl;
      delete gTemplateFitter;
      std::cout << "End of looper: delete mjjEbtdataRdm" << std::endl;
      delete mjjEbtdataRdm[ibtag];
      std::cout << "Scope of looper ends here" << std::endl;
    }  // end of the looper


    // write looper histograms to special root file
    std::cout << "Before opening looper output file " << std::endl;
    //TFile* fLooper = new TFile( Form("Looper_%s_%s.root",cMass,sbtag[ibtag].c_str()),"recreate");
    std::cout << "iSignal=" << iSignal 
	      << " signalMass=" << signalMass[iSignal]
	      << "  sbtag=" << sbtag[ibtag].c_str() << std::endl;
    std::string LooperOutputName( Form("Looper_M%04d_%s.root",signalMass[iSignal],sbtag[ibtag].c_str()) );
    std::cout << "Designed file name is " << LooperOutputName << std::endl;
    TFile* fLooper = new TFile( LooperOutputName.c_str(),"recreate");
    std::cout << "After opening looper output file " << std::endl;
    fLooper->cd();
    std::cout << "Before writing FitStatus " << std::endl;

    FitStatus_H->Write();
    HFraction_H->Write();
    HSeen_H->Write();
    HCorr_H->Write();
    HXSect_H->Write();
    HFractionPull_H->Write();
    HXSectPull_H->Write();
    HXSectUL95_H->Write();
    HtanBetaUL95_H->Write();
    std::cout << "Before closing " << std::endl;
    fLooper->Close();


  }  // loop over ibtag

  // close the packed template file
  ptFile->Close();
  //hout->Write();
  hout->Close();
}


int main(int narg,char** varg) {
//   int iMass = 2;
//   customFitS_Mass( iMass , true );
//   customFitS_Mass( 0 , false );

#include "Analysis/Utilities/interface/HbbMass.h"

   customFitS_Mass( 3 , true );

//   // bg-only fit
//   customFitS_Mass( 0 , false );
  // fits with signal
//  for (int iMass=0; iMass<nSignal; ++iMass) {
//    customFitS_Mass( iMass , true );
//  }
  return 0;
}
