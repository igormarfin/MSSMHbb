//
//  This is the customFitter for the "split mode", reading only packed template root file
//

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
  string name() {
    string sFlav[3] = {"Q", "C", "B"};
    string sLegend("bbb");
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

    // to be safe, make a copy of the histogram
    hist = new TH2F( *theHist );
    hist->SetName( Form("%s_tR", theHist->GetName()) );
    std::cout << "templateRef: Created a copy " << hist->GetName() << std::endl; 

    if (theColor == -1) {
      color = theHist->GetLineColor();
    } else {
      theHist->SetLineColor(color);
    }
//    sumW = theHist->GetSum();
    sumW = theHist->Integral();
  }
  ~templateRef() {
    delete hist;
    for (unsigned int iNuisance=0; iNuisance<tempPlus.size(); ++iNuisance) {
      delete tempPlus[iNuisance];
      delete tempMinus[iNuisance];
    }
  }
  void normalize(const int centralMode = 1) {
    if (hist->GetSum() >0) {
//      std::cout << hist->GetSum() << std::endl;
      std::cout << hist->Integral() << std::endl;
//      float scaleFacCentral = 1. / hist->GetSum();
      float scaleFacCentral = 1. / hist->Integral();
      hist->Scale( scaleFacCentral );
      for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
	float scaleFac = scaleFacCentral;
	if (centralMode == 0) {
	  if (tempPlus[iNuisance]->GetSum() >0) {
//	    scaleFac = 1. / tempPlus[iNuisance]->GetSum();
	    scaleFac = 1. / tempPlus[iNuisance]->Integral();
	  } else {
	    std::cout << "templateRef::normalize: Error: Plus template non-normalizable" << std::endl;
	    return;
	  }
	}
	tempPlus[iNuisance]->Scale( scaleFac );

	scaleFac = scaleFacCentral;
	if (centralMode == 0) {
	  if (tempMinus[iNuisance]->GetSum() >0) {
//	    scaleFac = 1. / tempMinus[iNuisance]->GetSum();
	    scaleFac = 1. / tempMinus[iNuisance]->Integral();

	  } else {
	    std::cout << "templateRef::normalize: Error: Minus template non-normalizable" << std::endl;
	    return;
	  }
	}
	tempMinus[iNuisance]->Scale( scaleFac );
      }
    } else {
//      std::cout << "templateRef::normalize(): bad sum of weights " << hist->GetSum() << std::endl;
      std::cout << "templateRef::normalize(): bad sum of weights " << hist->Integral() << std::endl;
    }
  }
  void scale(const float fScale) {
    hist->Scale( fScale );
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      tempPlus[iNuisance]->Scale( fScale );
      tempMinus[iNuisance]->Scale( fScale );
    }
  }
  void setNames(const int nSyst,std::string* systName,std::string* upDownName) {
    std::cout << "Here is setNames, nSyst=" << nSyst << std::endl;
    hist->SetName(Form("%s",this->name.c_str()));
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      std::cout << "setNames for iNuisance=" << iNuisance << std::endl;
      if (iNuisance < nSyst) {
	tempPlus[iNuisance]->SetName(Form("%s_%s_%s",this->name.c_str()
					  ,systName[iNuisance].c_str(),upDownName[0].c_str()));
	tempMinus[iNuisance]->SetName(Form("%s_%s_%s",this->name.c_str()
					   ,systName[iNuisance].c_str(),upDownName[1].c_str()));
      } else {
	std::cout << "templateRef::setNames: Bad  iNuisance=" << iNuisance << std::endl;
      }
    }
  }
  void write() {
    hist->Write();
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      tempPlus[iNuisance]->Write();
      tempMinus[iNuisance]->Write();
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
//      float theNormSlope = (theTempPlus->GetSum() - theTempMinus->GetSum()) / 
//	( 2 *  theNSigma * hist->GetSum() );
      float theNormSlope = (theTempPlus->Integral() - theTempMinus->Integral()) / 
	( 2 *  theNSigma * hist->Integral() );
      std::cout << "addNuisance: " << name << " " << tempPlus.size() << " theNormSlope= " << theNormSlope << std::endl;
      NormSlope.push_back( theNormSlope );
    }

    ++nNuisance;
  }
  // copy constructor
  templateRef(const templateRef& otr) : name(otr.name), inival(otr.inival),
					color (otr.color), sumW(otr.sumW), 
					nNuisance(otr.nNuisance), 
					nSigma(otr.nSigma) {
    //std::cout << "Here is the templateRef copy constructor for " << otr.name << std::endl;
    //std::cout << "Print the original: " << std::endl;
    //otr.Print();

    hist = new TH2F( *(otr.hist) );
    hist->SetName( Form("%s_tRc", otr.hist->GetName()) );
    for (unsigned int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      TH2F* thePlus = new TH2F( *(otr.tempPlus[iNuisance]) );
      thePlus->SetName( Form("%s_tRc", otr.tempPlus[iNuisance]->GetName()) );
      tempPlus.push_back( thePlus );
      TH2F* theMinus = new TH2F( *(otr.tempMinus[iNuisance]) );
      theMinus->SetName( Form("%s_tRc", otr.tempMinus[iNuisance]->GetName()) );
      tempMinus.push_back( theMinus );
    }
    //std::cout << "Print the result:" << std::endl;
    //this->Print();
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
      return mergedTemplate;
    }
    for (int iNuisance=0; iNuisance<tr1->nNuisance; ++iNuisance) {
      if ( tr1->nSigma[iNuisance] != tr2->nSigma[iNuisance] ) {
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
      if ( (tempPlus.size() != nNuisance) || (tempMinus.size() != nNuisance) ) {
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

  templateRef* posterior(const int mNuisance,  float* valNP, float* errNP ) {
    // create posterior templates
    std::string newName = this->name + "_Post";
    // first create copies
    templateRef* thePost = new templateRef( *this );
    thePost->name = newName;
    // then apply the nuisance parameter shifts
    for (int binx=1; binx<=hist->GetXaxis()->GetNbins(); ++binx) {
      for (int biny=1; biny<=hist->GetYaxis()->GetNbins(); ++biny) {
	// first shift the central histogram
	float oldContent = hist->GetBinContent(binx,biny);
	float theShift = 0;
	for (int iM=0; iM<mNuisance; ++iM) {  // loop over posterior NP
	  if ( valNP[iM]> 0) {
	    theShift += ( 0.5 * valNP[iM] * ( tempPlus[iM]->GetBinContent(binx,biny)
					      - hist->GetBinContent(binx,biny) ) );
	  } else {
	    theShift -= ( 0.5 * valNP[iM] * ( tempMinus[iM]->GetBinContent(binx,biny)
					      - hist->GetBinContent(binx,biny) ) );
	  }
	}
	float newContent = oldContent + theShift;
	thePost->hist->SetBinContent(binx,biny,newContent);
	// then create the new +/- 2 sigma contours
	for (int iM=0; iM<mNuisance; ++iM) {
	  if (iM < thePost->nNuisance) {
	    float theDeltaPlus = errNP[iM] * ( tempPlus[iM]->GetBinContent(binx,biny) 
					  - hist->GetBinContent(binx,biny) );
	    float thePostPlus = newContent + theDeltaPlus;
	    thePost->tempPlus[iM]->SetBinContent(binx,biny,thePostPlus); 
	    float theDeltaMinus = errNP[iM] * ( tempMinus[iM]->GetBinContent(binx,biny) 
					  - hist->GetBinContent(binx,biny) );
	    float thePostMinus = newContent + theDeltaMinus;
	    thePost->tempMinus[iM]->SetBinContent(binx,biny,thePostMinus);
// 	    std::cout << " np=" << iM << " binx=" << binx << " biny=" << biny
// 		      << " central=" << hist->GetBinContent(binx,biny)
// 		      << " theDeltaPlus=" << theDeltaPlus
// 		      << " theDeltaMinus=" << theDeltaMinus << std::endl;
	  } else {
	    std::cout << " Bad nuisance parameter " << iM << std::endl;
	    return NULL;
	  }
	}
      }
    }
    return thePost;
  }

  void Print() const {
    std::cout << "---- Template " << name << "  base hist " << hist->GetName()
//	      << " norm " << hist->GetSum() << std::endl;
	      << " norm " << hist->Integral() << std::endl;
    for (int iNuisance=0; iNuisance<nNuisance; ++iNuisance) {
      std::cout << "      " << "Nuisance parameter " << iNuisance << std::endl;
      float varPlus = -1;
      float varMinus = -1;
      if (hist->GetSum() != 0) {
//	varPlus = tempPlus[iNuisance]->GetSum() / hist->GetSum();
//	varMinus = tempMinus[iNuisance]->GetSum() / hist->GetSum();
	varPlus = tempPlus[iNuisance]->Integral() / hist->Integral();
	varMinus = tempMinus[iNuisance]->Integral() / hist->Integral();
      }
      std::cout << "          " << "Upper: " << std::setw(40) << tempPlus[iNuisance]->GetName() 
//		<< "  " << tempPlus[iNuisance]->GetSum() 
		<< "  " << tempPlus[iNuisance]->Integral() 
		<< " ratio " << varPlus << std::endl;
      std::cout << "          " << "Lower: " << std::setw(40) << tempMinus[iNuisance]->GetName() 
//		<< "  " << tempMinus[iNuisance]->GetSum() 
		<< "  " << tempMinus[iNuisance]->Integral() 
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
  TH2F* getTemplateHist(int iTemplate) {
    if ( (iTemplate<0) || (iTemplate>=nTemplates()) ) {
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
    TMinuit* gMinuit = new TMinuit(nFitPar);
    gMinuit->SetFCN(fcn);
    Double_t arglist[10];
    arglist[0] = 1;
    Int_t ierflg = 0;
    // set error definition to 1 (correct for chi2 fits)
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg);
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
      gMinuit->mnparm(i,Form("par%d",i),vstart[i],vstep[i],minConstraint[i],maxConstraint[i],ierflg);
    }
    // initialize the nuisance parameters
    for (unsigned int i=nTemplates(); i<nTemplates()+nNuisance; ++i) {
      vstart[i] = 0;
      vstep[i] = 0.1;
      minConstraint[i] = 0;
      maxConstraint[i] = 0;
      gMinuit->mnparm(i,Form("nupar%d",i),vstart[i],vstep[i],minConstraint[i],maxConstraint[i],ierflg);
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

//     arglist[0] = 0;
//     arglist[1] = 0;
//     gMinuit->mnexcm("HESSE",arglist,0,ierflg);


    return fitStat;
  };

  void getBGTotalError() {
    Int_t npars = gMinuit->GetNumPars();
    TMatrixDSym _corrMatrix(npars); 
    cout<<"External correlation Matrix from TMinuit"<<"\n";
    gMinuit->mnemat(_corrMatrix.GetMatrixArray(),npars );
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
    for (int iTemplate=0; iTemplate<nTemplates(); ++iTemplate) {
      B += frac[iTemplate];
    }
    for (int iTemplate=0; iTemplate<nTemplates(); ++iTemplate) {
      for (int jTemplate=0; jTemplate<nTemplates(); ++jTemplate) {
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

    Int_t npars = gMinuit->GetNumPars();
    TMatrixDSym _corrMatrix(npars); 
    cout<<"External correlation Matrix from TMinuit"<<"\n";
    gMinuit->mnemat(_corrMatrix.GetMatrixArray(),npars );
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

    for (int iTemplate=0; iTemplate<nTemplates(); ++iTemplate) {
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

  TH2F* getSummedTemplateHistNuisance(const int iNuisance, const int iUpDown) {
    TH2F* theSum = NULL;
    Int_t npars = gMinuit->GetNumPars();

    std::cout << "getSummedTemplateHistNuisance: iNuisance=" << iNuisance
	      << " iUpDown=" << iUpDown << std::endl;

    // fractions
    double frac[100];
    double theError;
    for (int i=0; i<npars;i++){
      GetResult(i,frac[i],theError);
    }

    for (int iTemplate=0; iTemplate<nTemplates(); ++iTemplate) {
      templateRef* theTemplate = &((*this->templates)[iTemplate]);
      TH2F* theHist;
      if (iUpDown ==0) {
	theHist = theTemplate->tempPlus[iNuisance];
      } else {
	theHist = theTemplate->tempMinus[iNuisance];
      }

      std::cout << "getSummedTemplateHistNuisance: iTemplate=" << iTemplate 
		<< " iUpDown=" << iUpDown << " theHist=" << theHist << std::endl;
//       std::cout << "Name=" << theHist->GetName() << std::endl;

      TH2F theScaledTemplateHist( *theHist);
      theScaledTemplateHist.SetName( Form("%sScld_NP%d_iup%d",theHist->GetName(),iNuisance,iUpDown) );
      theScaledTemplateHist.Scale( frac[iTemplate] );

      if ( theSum == NULL ) {
	theSum = new TH2F( theScaledTemplateHist );
	theSum->SetName( Form("%sSum_%d_%d",theHist->GetName(),iNuisance,iUpDown) );
      } else {
	theSum->Add( &theScaledTemplateHist );
      }
    }
//    std::cout << "getSummedTemplateHistNuisance: sum = " << theSum->GetSum() << std::endl;
    std::cout << "getSummedTemplateHistNuisance: sum = " << theSum->Integral() << std::endl;
    return theSum;
  }

  float getTemplateContWithNuisance(int npar,double* par,
				    unsigned int iTemplate,int binx,int biny,
				    float theError) {
    // get template content of bin ix,iy
    if ((iTemplate<0) || (iTemplate>= nTemplates())) {
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
    if ((iTemplate<0) || (iTemplate>= nTemplates())) {
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

    // here we support two formats, depending on the printErrors flag
    bool printErrors = true;

    pF.open(pullOutFile  );
    std::cout << "Pulls of nuisance parameters" << std::endl;
    double valArray[100];
    for (int iNuisance=0; iNuisance<this->nNuisance; ++iNuisance) {
      double value;
      double error;
      double vOverE = -99999;
      if (error != 0) vOverE = value / error;
      int getpar = gMinuit->GetParameter(nTemplates()+iNuisance,value,error);
      valArray[iNuisance] = value;
      //pF << Form("%8s   %4d  %12.6f  +-  %12.6f  ratio %10.4f",leadText,iNuisance,value,error,vOverE) << std::endl;
      if (printErrors) {
	pF << Form("%4d  %12.6f   %12.6f ",iNuisance,value,error) << std::endl;
      }
    }
    if (! printErrors) {
      pF << Form("%8s   ",leadText);
      for (int iNuisance=0; iNuisance<this->nNuisance; ++iNuisance) {
	pF << Form(" %12.6f ",valArray[iNuisance]);
      }
      pF << std::endl;
    }
    pF.close();
  }
  //
  templateRef* getPosterior(const int iTemplate) {
    if ((iTemplate<0) || (iTemplate>= nTemplates())) {
      std::cout << " getPosterior: bad template " << std::endl;
      return NULL;
    }
    // get the NP values
    float valArray[100];
    float errArray[100];
    for (int iNuisance=0; iNuisance<this->nNuisance; ++iNuisance) {
      double value;
      double error;
      int getpar = gMinuit->GetParameter(nTemplates()+iNuisance,value,error);
      valArray[iNuisance] = value;
      errArray[iNuisance] = error;
      std::cout << "templateFitter::getPosterior: np # " << iNuisance 
		<< " val=" << value << " err=" << error << std::endl;
    }
    templateRef* thePosterior = (*this->templates)[iTemplate].posterior(nNuisance,valArray,errArray);
    return thePosterior;
  }
};

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {    
  // this is the FCN defining the objective function
//  double T = gTemplateFitter->GetData()->GetSum();
  double T = gTemplateFitter->GetData()->Integral();
  double chisq = 0;
  // loop over the histogram bins
  for (int binx=gTemplateFitter->fitLow; binx<=gTemplateFitter->fitHigh; ++binx) {
    for (int biny=1; biny<=6; ++biny) {
      double errData = gTemplateFitter->GetData()->GetBinError(binx,biny);
      if (errData == 0) errData = 1;
      double delta = gTemplateFitter->GetData()->GetBinContent(binx,biny);
      double errMcSq = 0;
      for (int iTemplate=0; iTemplate<gTemplateFitter->nTemplates(); ++iTemplate) {
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
	  
	  for (int iTemplate=0; iTemplate<gTemplateFitter->nTemplates(); ++iTemplate) {
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

float makeSquare(const float x) { return (x*x); }

double getSystBbpurTpat( const string theScenario, const int theSignalMass ) {
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

  string sfc[nfc] = { "q", "c", "b" };

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
  int nSyst = 4;
  std::string systName[] = { "JES", "SFbc", "SFudsg", "JER" };

  bool activeNuisance[] = { true, true, true, true };

  double efficiency[nSignal][nbtag];
  bool signalNuisanceCentral = false;

  bool includeSys = true;

  bool drawTemplateErrors = true;

  bool addBbpurTpat = true;
 
  double intLumi = 0;
  if (scenario == "LowMass2011") {
    //intLumi = 2.66794; // in fb-1
    intLumi = 2.692643;  // with new method
  } else if (scenario == "MediumMass2011") {
    //intLumi = 3.99983; // in fb-1
    intLumi = 4.040802;
  } 

   else if (scenario == "MediumMass2012")
  {
   intLumi = 2.66;
  }

else {
    std::cout << "Bad scenario" << std::endl;
    return;
  }

//#if defined(MEDIUM2012)
//intLumi=2.66;
//nSyst=0; no syst
//#endif

  // mass label (for plots)
  string cMass;
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
  string ptFileName( Form("packedTemplates-M-%d.root",signalMass[iSignal]) );
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

  string tFlav[3] = {"Uds", "C", "B"};
  string sFlav[3] = {"Q", "C", "B"};

  // output file
  TFile* hout = new TFile(Form("customFitS-%s.root",cMass.c_str()),"recreate");
  hout->cd();
  TH2::AddDirectory(true);

  TRandom3 rdm;
  rdm.SetSeed(0);

 

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {


//#if !defined(MEDIUM2012)

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
    nominalXSect[ibtag][iSignal] = histXSect->GetBinContent(1);

    if (signalMass[iSignal] == 130) {
      nominalXSect[ibtag][iSignal] = 93.44;
      std::cout << "fix triple-merge problem for mH=" << signalMass[iSignal] << std::endl
		<< " Replace Xsection of " << histXSect->GetBinContent(1)
		<< " by " << nominalXSect[ibtag][iSignal] << std::endl;
    }

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

//#endif

    // read the background central templates
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	templateId tName(ifc,icateg);
	string bgCentralTemplateName( Form("%s_%s",tName.name().c_str(),sbtag[ibtag].c_str()) );
	hBackgroundCentral[ibtag][icateg][ifc] = 
	  (TH2F*) ptFile->Get( bgCentralTemplateName.c_str() );
	if ( hBackgroundCentral[ibtag][icateg][ifc] == NULL ) {
	  std::cout << "Hist not found: " << bgCentralTemplateName << std::endl;
	  return;
	}
	std::cout << "Read " << hBackgroundCentral[ibtag][icateg][ifc]->GetName() << std::endl;
      }
    }

//#if !defined(MEDIUM2012)

    // read the background syst templates
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	for (int iSyst=0; iSyst<nSyst; ++iSyst) {
	  for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	    templateId tName(ifc,icateg);
	    string bgSystTemplateName( Form("%s_%s_%s_%s",tName.name().c_str(),systName[iSyst].c_str(),
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

//#endif

    // read the hist-to-be-fitted
    string dataName( Form("Data_%s",sbtag[ibtag].c_str()) );
    mjjEbtdata[ibtag] = (TH2F*) ptFile->Get( dataName.c_str() );
    if (mjjEbtdata[ibtag] == NULL) {
      std::cout << "Histogram not found: " << dataName << std::endl;
      return;
    }
    std::cout << "Read " << mjjEbtdata[ibtag]->GetName() << std::endl;
	std::cout <<"Real data are "<< mjjEbtdata[ibtag]->Integral()<<std::endl;


    

  }

  // mass range bins
//  int lowBinX = 3;
//  int highBinX = 22;
  int lowBinX = 1;
  int highBinX = 25;

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
    mjjEbtdataProX[ibtag] = mjjEbtdata[ibtag]->ProjectionX(Form("%sProX",mjjEbtdata[ibtag]->GetName()),lowBinX,highBinX,"e");
    mjjEbtdataProX[ibtag]->SetMarkerStyle(20);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetTitle("m(Jet1 Jet2) [GeV/c^{2}]");
    mjjEbtdataProX[ibtag]->GetYaxis()->SetTitle("N");
//    mjjEbtdataProY[ibtag] = mjjEbtdata[ibtag]->ProjectionY(Form("%sProY",mjjEbtdata[ibtag]->GetName()),lowBinX,highBinX,"e");
    mjjEbtdataProY[ibtag] = mjjEbtdata[ibtag]->ProjectionY(Form("%sProY",mjjEbtdata[ibtag]->GetName()),0,6,"e");
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

  // prepare legend
  TLatex legend;
  legend.SetNDC();
  float xLegend = 0.55;
  float yLegend0 = 0.85;
  float yIncrement = 0.055;
  legend.SetTextAlign(12);
  legend.SetTextSize(0.035);

  // build the set of MC histograms

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
//#if !defined(MEDIUM2012)

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
  } ///fit signal


//#endif


//#if !defined(MEDIUM2012)

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

//#endif


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
      string mergedName( Form("(%sb)b",sFlav[ifc].c_str()) );
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
// 	templates.push_back(  templateRef("H",mjjEbtsignal[ibtag],0.001,2) );
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
    std::cout << "Normalize " << std::endl;
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
    std::cout << "After normalization " << std::endl;
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
    
    gTemplateFitter = new templateFitter(mjjEbtdata[ibtag], &templates,  nActiveNuisance );
    //gTemplateFitter = new templateFitter(mjjEbtdata[ibtag], &templates,  0 );  // first test without any nuisance parameters
    
    for (int iTemplate=0; iTemplate< gTemplateFitter->nTemplates(); ++iTemplate) {
      if ((*gTemplateFitter->templates)[iTemplate].name != "H" || constrainSignal) {
//	gTemplateFitter->Constrain(iTemplate+1,0.0,2*gTemplateFitter->data->GetSum() / 
//				   gTemplateFitter->getTemplateHist(iTemplate)->GetSum());
	gTemplateFitter->Constrain(iTemplate+1,0.0,2*gTemplateFitter->data->Integral() / 
				   gTemplateFitter->getTemplateHist(iTemplate)->Integral());
      }
    }
    gTemplateFitter->SetRangeX(lowBinX,highBinX);
    std::cout << std::endl
	      << "   Fit for signal:  " << signalMass[iSignal] << std::endl << std::endl; 

    Int_t status = gTemplateFitter->Fit(); 
    cout << sbtag[ibtag] << " fit status: " << status << endl;
    gTemplateFitter->printNPPulls(Form("NuisanceParPulls-%s.txt",cMass.c_str()),cMass.c_str());
    gTemplateFitter->getBGTotalError();

    TH2F* result = (TH2F*) gTemplateFitter->GetPlot();
    mjjEbtdata[ibtag]->SetMarkerStyle(20);
    mjjEbtdata[ibtag]->GetXaxis()->SetTitle("m(Jet1 Jet2) [GeV/c^{2}]");
    mjjEbtdata[ibtag]->GetYaxis()->SetTitle("EvtBTag");
    mjjEbtdata[ibtag]->GetZaxis()->SetTitle("N / 10 GeV/c^{2}");
    mjjEbtdata[ibtag]->ProjectionX()->Draw("Ep");
    //result->Draw("same");
    double value = 0;
    double error = 0;
    double fPenalty = -1;


    // build TObjArray of templates for plotting
    TObjArray *mc = new TObjArray();
    for (int iTemplate=0; iTemplate< gTemplateFitter->nTemplates(); ++iTemplate) {
      //mc->Add( (*gTemplateFitter->templates)[iTemplate].hist );
      TH2F* theTHist = gTemplateFitter->getTemplateHistWithNuisance(iTemplate);
      mc->Add( theTHist );
    }
    std::cout << "Number of templates = " << mc->GetEntries() << std::endl;




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
      //theTemp->Scale( value * mjjEbtdata[ibtag]->GetSum() / theTemp->GetSum() );
//      theTemp->Scale( value * mjjEbtdata[ibtag]->GetSum());
    std::cout<<"scale in stack "<<value * mjjEbtdata[ibtag]->Integral()<<std::endl;
    std::cout<<"integral before scale "<< theTemp->Integral()<<std::endl;

      theTemp->Scale( value * mjjEbtdata[ibtag]->Integral());
    std::cout<<"integral after scale "<< theTemp->Integral()<<std::endl;

      theTemp->GetXaxis()->SetRange( gTemplateFitter->fitLow,gTemplateFitter->fitHigh );
      //theTemp->SetLineColor(tc[icol]);
//      TH1D* theTempProX = theTemp->ProjectionX(Form("%sProX",theTemp->GetName()),0,-1,"e");
     TH1D* theTempProX = theTemp->ProjectionX(Form("%sProX",theTemp->GetName()),lowBinX,highBinX,"e");
       std::cout<<"integral of projection "<< theTempProX->Integral()<<std::endl;
      theTempProX->SetLineColor( theTemp->GetLineColor() );
      theTempProX->SetFillColor( theTemp->GetLineColor() );
      theTempProX->SetMarkerStyle(20);
      vProX.push_back( theTempProX );
      std::cout << "Created projection of " << theTemp->GetName() << " as " << theTempProX->GetName() << std::endl;
      hsX.Add(theTempProX);
      //theTempProX->Draw("hist,same");


    double stckint=0e0;
    TObjArray* arr= hsX.GetStack(); 
	if (arr){
    std::cout<<"Stack size "<<arr->GetSize()<<std::endl;
	for (int kkk=0;kkk<arr->GetSize();kkk++) { stckint+=((TH1D*)arr->At(kkk))->Integral();
    std::cout<<"Comp integral "<<((TH1D*)arr->At(kkk))->Integral()<<std::endl;}
    std::cout<<"stack integral "<< stckint<<std::endl;
    }

//      TH1D* theTempProY = theTemp->ProjectionY(Form("%sProY",theTemp->GetName()),lowBinX,highBinX,"e");
      TH1D* theTempProY = theTemp->ProjectionY(Form("%sProY",theTemp->GetName()),1,6,"e");
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

    TFile* fTot = new TFile(Form("TotalT_%s_%s.root",sbtag[ibtag].c_str(),cMass.c_str()),"recreate");
    fTot->cd();
    // draw the summed template hist
    TH2F* theSummedTemplateHist = gTemplateFitter->getSummedTemplateHist();
    theSummedTemplateHist->SetName(Form("TotalBG_%s_%s",cMass.c_str(),sbtag[ibtag].c_str()));
    theSummedTemplateHist->SetTitle( theSummedTemplateHist->GetName());
    // 
//    theSummedTemplateHist->Scale( mjjEbtdata[ibtag]->GetSum());
    theSummedTemplateHist->Scale( mjjEbtdata[ibtag]->Integral());
    theSummedTemplateHist->Draw();
    theSummedTemplateHist->Write();

    // make total template hist with error bars
    TH2F* theTotal = new TH2F( *theSummedTemplateHist );
    theTotal->SetName( Form("%s_SimErrors",theSummedTemplateHist->GetName()) );
    theTotal->SetTitle( theTotal->GetName() );
    // simulate Poissonian errors
    for (int binx=1; binx<=theTotal->GetXaxis()->GetNbins(); ++binx) {
      for (int biny=1; biny<=theTotal->GetYaxis()->GetNbins(); ++biny) {
	float oldError = theTotal->GetBinError(binx,biny);
	float newError = 0;
	if (theTotal->GetBinContent(binx,biny)>=0) {
	  newError = sqrt( theTotal->GetBinContent(binx,biny) );
	}
	theTotal->SetBinError(binx,biny,newError);
	std::cout << " old error: " << oldError << "  new error " << newError << std::endl;
      }
    }
    theTotal->Write();
    
    // make randomized total template hist
    TH2F* theRdm = new TH2F( *theTotal );
    theRdm->SetName( Form("%s_rdm",theTotal->GetName()) );
    theRdm->SetTitle( theRdm->GetName() );
    for (int binx=1; binx<=theRdm->GetXaxis()->GetNbins(); ++binx) {
      for (int biny=1; biny<=theRdm->GetYaxis()->GetNbins(); ++biny) {
	float oldContent = theRdm->GetBinContent(binx,biny);
	float newContent = rdm.Poisson( theRdm->GetBinContent(binx,biny) );
	theRdm->SetBinContent(binx,biny,newContent);
	std::cout << "old content " << oldContent
		  << " new content " << newContent << std::endl;
      }
    }
    theRdm->Write();

    // get also the nuisance parameter variants
    std::cout << "Loop until nNuisance=" << nSyst << std::endl;
    for (int iNuisance=0; iNuisance<nSyst; ++iNuisance) {
      for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	std::cout << "call getSummedTemplateHistNuisance iNuisance=" << iNuisance << std::endl;
	TH2F* theSummedTemplateHistNP = gTemplateFitter->getSummedTemplateHistNuisance(iNuisance,iUpDown);
	theSummedTemplateHistNP->SetName(Form("TotalBG_%s_%s_%s_%s",cMass.c_str(),sbtag[ibtag].c_str(),
					      systName[iNuisance].c_str(),upDownName[iUpDown].c_str()));
//	theSummedTemplateHistNP->Scale( mjjEbtdata[ibtag]->GetSum());
	theSummedTemplateHistNP->Scale( mjjEbtdata[ibtag]->Integral());
	theSummedTemplateHistNP->Write();
	
      }
    }

    if (! fitSignal) {
      // background posterior templates
      TFile* fPost = new TFile(Form("PosteriorBG_%s_%s.root",sbtag[ibtag].c_str(),cMass.c_str()),"recreate");
      fPost->cd();    // make the posterior templates
      ofstream fracOut;
      fracOut.open( Form("FractionsBG_%s_%s.txt",sbtag[ibtag].c_str(),cMass.c_str()) );
      for (int iTemplate=0; iTemplate< gTemplateFitter->nTemplates(); ++iTemplate) {
	templateRef* thePostTemplate = gTemplateFitter->getPosterior( iTemplate );
	gTemplateFitter->GetResult(iTemplate, value, error);
	std::cout << "Post template name " << thePostTemplate->name << std::endl;
//	thePostTemplate->scale( value * mjjEbtdata[ibtag]->GetSum());

    std::cout<<"scale in PosteriorBG_ "<<value * mjjEbtdata[ibtag]->Integral()<<std::endl;
    std::cout<<"integral before scale "<< (thePostTemplate->hist)->Integral()<<std::endl;

	thePostTemplate->scale( value * mjjEbtdata[ibtag]->Integral());
    std::cout<<"integral after scale "<< (thePostTemplate->hist)->Integral()<<std::endl;

	thePostTemplate->setNames(nSyst,systName,upDownName);
	thePostTemplate->write();
	fracOut << Form(" %10s  %12.6f  +-  %12.6f  ",thePostTemplate->name.c_str(),value,error) << std::endl;
      }
      fracOut.close();
    }

    if (fitSignal) {
      // signal posterior templates
      ifstream npIn;
      float valNp[100];
      float errNp[100];
      npIn.open("NuisanceParPulls-BGonly.txt");
      for (int iNuisance=0; iNuisance<4; ++iNuisance) {
	int theNp;
	npIn >> theNp >> valNp[iNuisance] >> errNp[iNuisance];
      }
      if (npIn.eof()) {
	std::cout << "Premature eof on npIn file" << std::endl;
      } else {
	valNp[0] = 0; // JES
	errNp[0] = 0;
	valNp[3] = 0; // JER
	errNp[3] = 0;

	templateRef* postSignalTemplate = mjjEbtSignalTemplateRef[ibtag]->posterior(4,valNp,errNp);
	postSignalTemplate->setNames(nSyst,systName,upDownName);
	TFile* fPostSignal = new TFile(Form("PosteriorSignal_%s_%s.root",
					    sbtag[ibtag].c_str(),cMass.c_str()),"recreate");
	fPostSignal->cd();
	postSignalTemplate->write();
      }
    }


    //fTot->Close();
    hout->cd();

    TH1D* theSummedTemplateHistProX = 
  //    theSummedTemplateHist->ProjectionX(Form("%sProX",theSummedTemplateHist->GetName()),0,-1,"e");
      theSummedTemplateHist->ProjectionX(Form("%sProX",theSummedTemplateHist->GetName()),lowBinX,highBinX,"e");
    TH1D* theSummedTemplateHistProY =
//      theSummedTemplateHist->ProjectionY(Form("%sProY",theSummedTemplateHist->GetName()),lowBinX,highBinX,"e");
      theSummedTemplateHist->ProjectionY(Form("%sProY",theSummedTemplateHist->GetName()),1,6,"e");
    theSummedTemplateHistProX->SetFillColor(7);
    theSummedTemplateHistProX->Draw("PE2");
    canvas->Print( Form("%s_SummedTemplateHistX_%s.png",cMass.c_str(),sbtag[ibtag].c_str()) );
    
    theSummedTemplateHistProY->SetFillColor(7);
    theSummedTemplateHistProY->Draw("PE2");
    canvas->Print( Form("%s_SummedTemplateHistY_%s.png",cMass.c_str(),sbtag[ibtag].c_str()) );

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
      //gTemplateFitter->GetResult(iH, value, error);
      // determine the fitting penalty
      fPenalty =  HFractionError[ibtag][iSignal] * mjjEbtdata[ibtag]->GetSum() / sqrt( mjjEbtdata[ibtag]->GetSum() );
      std::cout << "deltaS / sqrt(B) = " << error * mjjEbtdata[ibtag]->GetSum()
		<< "  /  " << sqrt( mjjEbtdata[ibtag]->GetSum() ) << "  =  "
		<< fPenalty << std::endl;
      std::cout << " HFraction=" << HFraction[ibtag][iSignal] 
		<< " +- " << HFractionError[ibtag][iSignal] 
		<< " SumW= " << mjjEbtdata[ibtag]->GetSum()
		<< std::endl;


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
      std::cout << " Signal stats rel error: " << sysSignalStatistics << std::endl;
      deltaSigmaSysSq += makeSquare(sysSignalStatistics);

      for (int iNuisance=0; iNuisance<mjjEbtSignalTemplateRef[ibtag]->nNuisance; ++iNuisance) {
	float nuisNorm = mjjEbtSignalTemplateRef[ibtag]->getNormSlope(iNuisance);
	std::cout<< " iTemplate=" << iNuisance << " nuisNorm=" << nuisNorm << std::endl;
	deltaSigmaSysSq += makeSquare(nuisNorm); // SFbc effect on signal
      }
      std::cout << "Relative deltaSigmaSys=" << sqrt( deltaSigmaSysSq ) << std::endl;

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
      

      // print the summary line
      std::cout << Form("| %6s | %28s+H | %10.2f | %d | %8.2f | %8.2f |",
			sbtag[ibtag].c_str(),model.c_str(),
			gTemplateFitter->GetChisquare(),nInactive,fPenalty,
			HXSectError[ibtag][iSignal]) << std::endl;

      // write the summary card to a file
      ofstream headCard;
      headCard.open("sumCard-Header.txt",ios::out);
      headCard << Form("| %4s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %10s | %5s |",
		       "Mass","Eff","chi2","frac","Efrac","XSect","EXSect","XSectUL95","tbUL95","theoXS","fstat")
	       << std::endl;
      ofstream sumCard;
      sumCard.open(Form("sumCard-M%04d-%s-%s.txt",signalMass[iSignal],
			sbtag[ibtag].c_str(),scenario.c_str()),ios::out);
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
//    std::cout<<"Final stack integral"<< hsX.Integral()<<std::endl;
    double stckint=0e0;
    TObjArray* arr= hsX.GetStack(); 
	if (arr){
	for (int kkk=0;kkk<arr->GetSize();kkk++) stckint+=((TH1D*)arr->At(kkk))->Integral();
    std::cout<<"stack final integral "<< stckint<<std::endl;
     }

    hsX.Draw("HIST,SAME");
    if (drawTemplateErrors) {
      theSummedTemplateHistProX->Draw("PE2,SAME");
//       theSummedTemplateHistProX->SetFillColor(0);
//       theSummedTemplateHistProX->Draw("HIST,SAME");
    }
    //mjjEbtdata[ibtag]->ProjectionX()->Draw("EpSAME");
    std::cout << "Draw " << mjjEbtdataProX[ibtag]->GetName() << " Number of bins " << mjjEbtdataProX[ibtag]->GetXaxis()->GetNbins() << std::endl;
    mjjEbtdataProX[ibtag]->Draw("EpSAME");
    //result->Draw("same");
    for (int iTemplate=0; iTemplate<templates.size(); ++iTemplate) {
      gTemplateFitter->GetResult(iTemplate, value, error);
      string sLegend(Form("%6s %8.3f +- %8.3f",templates[iTemplate].name.c_str(),value,error));
      legend.SetTextColor(templates[iTemplate].color);
      legend.DrawLatex(xLegend+0.1,yLegend0-0.2-iTemplate*yIncrement,sLegend.c_str());
    }        // print also the chi2
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
    if (drawTemplateErrors) {
      theSummedTemplateHistProY->Draw("PE2,SAME");
    }

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

#ifdef DONOTHING
    // in case we fitted with signal, scale up
    next.Reset();
    while ( (theTemp = (TH2F*) next()) ) {
      for (int iTemplate=0; iTemplate<templates.size(); ++iTemplate) {
	gTemplateFitter->GetResult(iTemplate, value, error);
	if (templates[iTemplate].name == "H") {
	  std::cout << "Higgs fit values : " << Form("%6s %8.3f +- %8.3f",templates[iTemplate].name.c_str(),value,error) << std::endl;
	  double T = gTemplateFitter->GetData()->GetSum();
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

  // close the packed template file
  ptFile->Close();
  hout->Write();
  hout->Close();
}


void customFitPostA() {
//   int iMass = 2;
//   customFitS_Mass( iMass , true );
//   customFitS_Mass( 0 , false );

#include "Analysis/Utilities/interface/HbbMass.h"


  // bg-only fit
  customFitS_Mass( 0 , false );
  // fits with signal
#if !defined(MEDIUM2012)
  for (int iMass=0; iMass<nSignal; ++iMass) {
    customFitS_Mass( iMass , true );
  }
#endif
  return;
}
