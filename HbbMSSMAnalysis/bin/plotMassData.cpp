//
//  This is the mass spectrum plotter, based on TripleBtagAnalysis.root
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
#include <TPave.h>
#include <TLegend.h>

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
   tdrStyle->SetOptStat(0); // To display the mean and RMS:    
   //   tdrStyle->SetOptStat("mrei");
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
   //tdrStyle->SetPadRightMargin(0.02);
   tdrStyle->SetPadRightMargin(0.03);

// For the Global title:

   tdrStyle->SetOptTitle(0);
   //tdrStyle->SetOptTitle(1);
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


void cmsPrel(double intLumi){
   TLatex latex;
   latex.SetNDC();
   //latex.SetTextSize(0.04);
   latex.SetTextSize(0.03);

   latex.SetTextAlign(31); // align right
   //latex.DrawLatex(0.98,0.965,"#sqrt{s} = 7 TeV");
   if (intLumi > 0.) {
     latex.SetTextAlign(31); // align right
     //latex.DrawLatex(0.98,0.88,Form("CMS Preliminary #int #font[12]{L}dt = %.1f  fb^{-1}",intLumi));
     latex.DrawLatex(0.98,0.965,Form("CMS WorkInProgress 2012, L = %.1f  fb^{-1}, #sqrt{s} = 8 TeV",intLumi));
     //latex.DrawLatex(0.98,0.965,Form("CMS 2011, L_{int} = %.1f  fb^{-1}, #sqrt{s} = 7 TeV",intLumi));
     //latex.DrawLatex(0.98,0.965,Form("CMS Preliminary 2011, L = %.1f  fb^{-1}, #sqrt{s} = 7 TeV",intLumi));
   }
   latex.SetTextAlign(11); // align left
   //latex.DrawLatex(0.02,0.965,"#font[22]{CMS Preliminary, L=2.7 fb^{-1}");
   return;
}



float makeSquare(const float x) { return (x*x); }

void plotMassData(char* theScenario) {
  // plot data mass spectra
  setTDRStyle();
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  canvas = new TCanvas ("cg1","PadX",10,10,800,800);
  double yScale = 800./600.;
  gStyle->SetPadColor(0);
  canvas->SetFillColor(0);

  const int ncateg = 3;
  const int ncorr=2;
  const int nfc = 3;
  const int nbtag = 2;
  // const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  const std::string sbtag[nbtag] = { "CSVT", "TCHPT"};
  string strig("");

  string sfc[nfc] = { "q", "c", "b" };

  //std::string scenario("LowMass2011");
  //std::string scenario("MediumMass2011");
  std::string scenario;
  bool useTemplateError;
  bool useNP;
  //scenario.assign("MediumMass2012");
  scenario.assign(theScenario);

  double intLumi = 0;
  if (scenario == "MediumMass2012") {
    intLumi = 3.8757; // in fb-1
    strig.assign("Trig0");
  } else if (scenario == "HighMass2012") {
    intLumi = 18.6485; // in fb-1
    strig.assign("Trig1");
  } else if (scenario == "VeryHighMass2012") {
    intLumi = 18.6485; // in fb-1
    strig.assign("Trig2");
  } else if (scenario == "UCLAd755") {
    intLumi = 19.3; // in fb-1
    strig.assign("Trig7");
  } else if (scenario == "UCLAd753") {
    intLumi = 19.3; // in fb-1
    strig.assign("Trig7");
  } else if (scenario == "SUSYMC") {
    intLumi = 0; // in fb-1
    strig.assign("AllTrig");
  } else {
    std::cout << "Bad scenario" << std::endl;
    return;
  }


  // open the Root file
  string ptFileName( "TripleBtagAnalysis.root" );
  std::cout << "Open packed template file " << ptFileName << std::endl;
  TFile* ptFile = new TFile( ptFileName.c_str() );
  if (ptFile == NULL) {
    std::cout << "Error opening " << ptFileName << std::endl;
    return;
  }

  // arrays to hold the hist-to-be-fitted, the signal & background templates
  TH2F* mjjEbtdata[nbtag];
  TH1F* mjjdata[nbtag];

  string tFlav[3] = {"Uds", "C", "B"};
  string sFlav[3] = {"Q", "C", "B"};
  string sLogCode[2] = {"Lin","Log"};

  // output file
  TFile* hout = new TFile(Form("plotMassData.root"),"recreate");
  hout->cd();
  TH2::AddDirectory(true);

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {

    // read the data mass spectra (2D)
    string dataName( Form("mjjBtcut/mjjbt%s%s",sbtag[ibtag].c_str(),strig.c_str()) );
    mjjdata[ibtag] = (TH1F*) ptFile->Get( dataName.c_str() );
    if (mjjdata[ibtag] == NULL) {
      std::cout << "Histogram not found: " << dataName << std::endl;
      return;
    }
    std::cout << "Read " << mjjdata[ibtag]->GetName() << std::endl;

    dataName.assign( Form("massEvBtag/mjjEvBTag_%s%s",sbtag[ibtag].c_str(),strig.c_str()) );
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

  // prepare the drawing
  const int maxcol = 100;
  int tc[maxcol];
  tc[0] = kBlack;
  tc[1] = kBlue;
  tc[2] = kRed;
  tc[3] = kSpring+9; //kGreen;
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
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtdataProX[ibtag] = mjjEbtdata[ibtag]->ProjectionX(Form("%sProX",mjjEbtdata[ibtag]->GetName()),0,-1,"e");
    mjjEbtdataProX[ibtag]->SetMarkerStyle(20);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetTitle("M_{12} [GeV]");
    mjjEbtdataProX[ibtag]->GetXaxis()->SetTitleFont(42);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetLabelFont(42);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetLabelSize(0.05/yScale);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetTitleOffset(1.1);
    mjjEbtdataProX[ibtag]->GetXaxis()->SetTitleSize(0.05/yScale);
    mjjEbtdataProX[ibtag]->GetYaxis()->SetTitle("Events / 20 GeV");
    mjjEbtdataProX[ibtag]->GetYaxis()->SetLabelSize(0.05/yScale);
    mjjEbtdataProX[ibtag]->GetYaxis()->SetTitleSize(0.05/yScale);
    mjjEbtdataProX[ibtag]->GetYaxis()->SetTitleOffset(2.1);
    mjjEbtdataProY[ibtag] = mjjEbtdata[ibtag]->ProjectionY(Form("%sProY",mjjEbtdata[ibtag]->GetName()),lowBinX,highBinX,"e");
    mjjEbtdataProY[ibtag]->SetMarkerStyle(20);
    mjjEbtdataProY[ibtag]->GetXaxis()->SetTitle("X_{123}");
    mjjEbtdataProY[ibtag]->SetTitleFont(42);
    mjjEbtdataProY[ibtag]->SetLabelFont(42);
    mjjEbtdataProY[ibtag]->GetXaxis()->SetLabelSize(0.05/yScale);
    mjjEbtdataProY[ibtag]->GetXaxis()->SetTitleOffset(1.1);
    mjjEbtdataProY[ibtag]->GetXaxis()->SetTitleSize(0.05/yScale);
    mjjEbtdataProY[ibtag]->GetYaxis()->SetTitle("Events");
    mjjEbtdataProY[ibtag]->GetYaxis()->SetLabelSize(0.05/yScale);
    mjjEbtdataProY[ibtag]->GetYaxis()->SetTitleSize(0.05/yScale);
    mjjEbtdataProY[ibtag]->GetYaxis()->SetTitleOffset(2.1);
    //dump2D( mjjEbtdata[ibtag] );
  }

  // prepare legend
  TLatex legend;
  legend.SetNDC();
  float xLegend = 0.55;
  float yLegend0 = 0.85 + 0.05 * (yScale - 1);
  float yIncrement = 0.055 / yScale;
  legend.SetTextAlign(12);
  legend.SetTextSize(0.035);

  // do the drawing
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjdata[ibtag]->SetMarkerStyle(20);
    mjjdata[ibtag]->GetXaxis()->SetTitle("M_{12} [GeV]");
    mjjdata[ibtag]->GetXaxis()->SetTitleFont(42);
    mjjdata[ibtag]->GetXaxis()->SetLabelFont(42);
    mjjdata[ibtag]->GetXaxis()->SetLabelSize(0.05/yScale);
    mjjdata[ibtag]->GetXaxis()->SetTitleOffset(1.1);
    mjjdata[ibtag]->GetXaxis()->SetTitleSize(0.05/yScale);
    mjjdata[ibtag]->GetYaxis()->SetTitle("Events / 20 GeV");
    mjjdata[ibtag]->GetYaxis()->SetLabelSize(0.05/yScale);
    mjjdata[ibtag]->GetYaxis()->SetTitleSize(0.05/yScale);
    mjjdata[ibtag]->GetYaxis()->SetTitleOffset(2.1);

    // loop over lin/log
    for (int iLogCode=0; iLogCode<2; ++iLogCode) {
      canvas->SetLogy( iLogCode );
      mjjdata[ibtag]->Draw("Ep");
      legend.DrawLatex(xLegend,yLegend0,(scenario+"  "+sbtag[ibtag]).c_str());
      cmsPrel(intLumi);
      canvas->Print( Form("mjjData_%s_%s.png",sbtag[ibtag].c_str(),sLogCode[iLogCode].c_str()) );

      mjjEbtdataProX[ibtag]->Draw("Ep");
      legend.DrawLatex(xLegend,yLegend0,(scenario+"  "+sbtag[ibtag]).c_str());
      cmsPrel(intLumi);
      canvas->Print( Form("mjjDataPro_%s_%s.png",sbtag[ibtag].c_str(),sLogCode[iLogCode].c_str()) );
      if (scenario == "SUSYMC") {
	mjjEbtdataProX[ibtag]->Fit("gaus","","",650.,1050.);
	canvas->Print( Form("mjjDataPro_%s_%s_Fit.png",sbtag[ibtag].c_str(),sLogCode[iLogCode].c_str()) );
      }
    }
    // linear view range
    canvas->SetLogy( 0 );
    mjjEbtdataProX[ibtag]->GetXaxis()->SetRange(33,75);
    legend.DrawLatex(xLegend,yLegend0,(scenario+"  "+sbtag[ibtag]).c_str());
    cmsPrel(intLumi);
    canvas->Print( Form("mjjDataPro_%s_Range.png",sbtag[ibtag].c_str()) );

    // print the number of triple btag events
    int nEntries = mjjEbtdata[ibtag]->GetEntries();
    std::cout << Form("Number of events (%s): %d",sbtag[ibtag].c_str(),
		      nEntries) << std::endl;

    // histogram output
    //mjjEbtdataProX[ibtag]->Write();
  }


  // close the packed template file
  ptFile->Close();
  hout->Write();
  hout->Close();
}

