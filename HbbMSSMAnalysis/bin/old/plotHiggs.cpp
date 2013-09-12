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




void plotHiggs() {
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

  std::string scenario("LowMass2011");

  const int nSignal=7;
  //const int nSignal=1;
  int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };
  //int signalMass[nSignal] = { 120 };

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {

    // arrays to hold the numbers
    const int maxPos=500;
    double mH[maxPos];
    double eff[maxPos];
    double chi2[maxPos];
    double frac[maxPos];
    double efrac[maxPos];
    double XSect[maxPos];
    double EXSect[maxPos];
    double XSectUL95[maxPos];
    double tbUL95[maxPos];
    double theoXS[maxPos];
    double fstat[maxPos];

    // open the sumCard files
    int nPos = 0;
    for (int iSignal=0; iSignal<nSignal; ++iSignal) {
      ifstream sumCard;
      sumCard.open(Form("sumCard-M%04d-%s-%s.txt",signalMass[iSignal],
			sbtag[ibtag].c_str(),scenario.c_str()));
      std::string sLine;
      while ( (! sumCard.eof()) && (nPos < maxPos) ) {
	getline(sumCard,sLine);
	if (! sumCard.eof()) {
	  std::cout << sLine << std::endl;

	  replace(sLine.begin(),sLine.end(),'|',' ');
	  std::cout << sLine << std::endl;

	  istringstream iss(sLine,std::istringstream::in);
	  double theMH;
	  double theEff;
	  double theChi2;
	  double theFrac;
	  double theEfrac;
	  double theXSect;
	  double theEXSect;
	  double theXSectUL95;
	  double theTbUL95;
	  double theTheoXS;
	  double theFstat;
	  iss >> theMH >> theEff >> theChi2 >> theFrac >> theEfrac
	      >> theXSect >> theEXSect >> theXSectUL95 >> theTbUL95
	      >> theTheoXS >> theFstat;
	  mH[nPos] = theMH;
	  eff[nPos] = theEff;
	  chi2[nPos] = theChi2;
	  frac[nPos] = theFrac;
	  efrac[nPos] = theEfrac;
	  XSect[nPos] = theXSect;
	  EXSect[nPos] = theEXSect;
	  XSectUL95[nPos] = theXSectUL95;
	  tbUL95[nPos] = theTbUL95;
	  theoXS[nPos] = theTheoXS;
	  fstat[nPos] = theFstat;
	  ++nPos;
	}
      }
    }

    // prepare the plotting
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
    sTitle[jStati] = "#epsilon_{H}";
    theStat[jStati] = eff;
    theError[jStati] = NULL;
    Pos[jStati] = 1;

    ++jStati;
    sName[jStati] = "Chi2";
    sTitle[jStati] = "#chi^{2}";
    theStat[jStati] = chi2;
    theError[jStati] = NULL;
    Pos[jStati] = 1;
    
    ++jStati;
    sName[jStati] = "Fraction";
    sTitle[jStati] = "H fraction";
    theStat[jStati] = frac;
    theError[jStati] = efrac;
    Pos[jStati] = 0;

    ++jStati;
    sName[jStati] = "XSect";
    sTitle[jStati] = "#sigma_{H} [pb]";
    theStat[jStati] = XSect;
    theError[jStati] = EXSect;
    Pos[jStati] = 0;

    ++jStati;
    sName[jStati] = "XSectUL95";
    sTitle[jStati] = "UL(95CL) #sigma_{H} [pb]";
    theStat[jStati] = XSectUL95;
    theError[jStati] = NULL;
    Pos[jStati] = 1;
    
    ++jStati;
    sName[jStati] = "TbUL95";
    sTitle[jStati] = "UL(95CL) tan #beta";
    theStat[jStati] = tbUL95;
    theError[jStati] = NULL;
    Pos[jStati] = 1;
    
    ++jStati;
    sName[jStati] = "TheoXS";
    sTitle[jStati] = "MSSM Theory #sigma_{H} [pb] (tan #beta=20)";
    theStat[jStati] = theoXS;
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
	//theX[iSignal] = signalMass[iSignal];
	theX[iSignal] = mH[iSignal];
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
      TH2F hist("hdum","hdum",40,50,400,40,theYMin,theYMax);
      hist.GetXaxis()->SetTitle("M_{H} [GeV/c^{2}]");
      hist.GetYaxis()->SetTitle(sTitle[iStati].c_str());    
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
      canvas->Print(Form("plh_%s_%s.png",sName[iStati].c_str(),sbtag[ibtag].c_str()));
    }



    // here we should plot
    for (int ii=0; ii<nPos; ++ii) {
      std::cout << " mass= " << mH[ii] << std::endl;
    } 
  }
}
      
