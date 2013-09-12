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




void findMedian() {
  // 
  setTDRStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  canvas = new TCanvas ("cg1","PadX",10,10,800,600);
  gStyle->SetPadColor(0);
  canvas->SetFillColor(0);

  const int nbtag = 1;
  // const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  const std::string sbtag[nbtag] = { "CSVT"};

  //string sfc[nfc] = { "q", "c", "b" };

  //std::string scenario("LowMass2011");
  std::string scenario("MediumMass2011");

//   const int nSignal=7;
//   //const int nSignal=1;
//   int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };
//   //int signalMass[nSignal] = { 120 };

#include "Analysis/Utilities/interface/HbbMass.h"

  //int iSignal = 6;

  std::cout << "Printout for scenario " << scenario << std::endl;
  for (int iSignal=0; iSignal<nSignal; ++iSignal) {

    for (int ibtag=0; ibtag<nbtag; ++ibtag) {

      // open the Root file
      std::string sName( Form( "LooperMerge_M%04d_%s.root",signalMass[iSignal],sbtag[ibtag].c_str()) );
      TFile* fLooper = new TFile( sName.c_str() );
      if (fLooper == NULL) {
	std::cout << "File not found " << sName << std::endl;
	continue;
      }
      TH1D* h = (TH1D*) fLooper->Get("HXSectUL95_H");
      if (h == NULL) {
	std::cout << "Histogram not found" << std::endl;
	break;
      }

      // find the observed limit
      ifstream sumCard;
//       sumCard.open(Form("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/fit-Expected/LowMass2011/NPon/BTon/sumCard-M%04d-%s-%s.txt",signalMass[iSignal],
// 			sbtag[ibtag].c_str(),scenario.c_str()));
      sumCard.open(Form("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/fit-V6-moreSyst/%s/NPon/BTon/sumCard-M%04d-%s-%s.txt",scenario.c_str(),signalMass[iSignal],
			sbtag[ibtag].c_str(),scenario.c_str()));
      double theMH = 0;
      double theEff = 0;
      double theChi2 = 0;
      double theFrac = 0;
      double theEfrac = 0;
      double theXSect = 0;
      double theEXSect = 0;
      double theXSectUL95 = 0;
      double theTbUL95 = 0;
      double theTheoXS = 0;
      double theFstat = 0;
      std::string sLine;
      getline(sumCard,sLine);
      if (! sumCard.eof()) {
	//std::cout << sLine << std::endl;
	replace(sLine.begin(),sLine.end(),'|',' ');
	//std::cout << sLine << std::endl;
	std::istringstream iss(sLine,std::istringstream::in);
	iss >> theMH >> theEff >> theChi2 >> theFrac >> theEfrac
	    >> theXSect >> theEXSect >> theXSectUL95 >> theTbUL95
	    >> theTheoXS >> theFstat;
      } else {
	std::cout << "Could not read observed limit " << std::endl;
      }
      sumCard.close();
      

      double xmin = h->GetXaxis()->GetXmin();
      double xmax = h->GetXaxis()->GetXmax();
      int nx = h->GetXaxis()->GetNbins();
      double dx = (xmax-xmin) / nx;
    
      double theMedian = -999999;
      
      double yTotal = 0;
      // first, determine the total content
      for (int binx=0; binx<=(nx+1); ++binx) {
	yTotal += h->GetBinContent( binx );
      }
      //std::cout << "yTotal=" << yTotal << " GetSum=" << h->GetSum() << std::endl;

      // find the median
      int binMed = -1;
      double yMed = 0;
      for (int binx=0; binx<=(nx+1); ++binx) {
	yMed += h->GetBinContent( binx );
	if (yMed > (0.5*yTotal)) {
	  binMed = binx;
	  break;
	}
    }

      if (binMed != nx) {
	theMedian = h->GetBinLowEdge( binMed+1 );
      } else {
	std::cout << "Funny median" << std::endl;
	theMedian = h->GetBinLowEdge( binMed );
      }

      //std::cout << "theMedian = " << theMedian << std::endl;
      
      double quant[2] = { 0.682, 0.954 };
    
      double crit[5];
      crit[0] = yTotal * 0.5;
      crit[1] =  yTotal * (0.5 - 0.5*quant[0]);
      crit[2] =  yTotal * (0.5 + 0.5*quant[0]);
      crit[3] =  yTotal * (0.5 - 0.5*quant[1]);
      crit[4] =  yTotal * (0.5 + 0.5*quant[1]);
      
      double val[5] = {0, 0, 0, 0, 0};
      for (int ii=0; ii<5; ++ii) {
	val[ii] = 0;
	double ySum = 0;
	int binSum = -1;
	for (int binx=0; binx<=(nx+1); ++binx) {
	  ySum += h->GetBinContent( binx );
	  if (ySum > crit[ii]) {
	    binSum = binx;
	    //std::cout << "ii=" << ii << "binSum=" << binSum << std::endl;
	    break;
	  }
	}
	
	val[ii] = xmin + binSum * dx;
	//std::cout << " ii=" << ii << " binSum=" << binSum << " val=" << val[ii] << std::endl;
	
      }

      std::cout << signalMass[iSignal] << "  " << theXSectUL95 << "  " ;
      for (int ii=0; ii<5; ++ii) {
	std::cout << val[ii] << "  " ;
      }
      std::cout << std::endl;

      fLooper->Close();
    }  // ibtag
  } // iSignal
}
      
