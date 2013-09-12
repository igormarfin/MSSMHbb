#include "lsDraw.h"
#include "lsUtils.h"
#include "stringUtils.h"
#include "lsConfig.h"

#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TTree.h>

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

void drawAsci(TH1* h, asciDrawOptions* drawOpts, bool originalCall){
  TH1* hOriginal = h;
  if(h->InheritsFrom("TH2")) { drawAsci2( (TH2*) h); return; }
  if(!h->InheritsFrom("TH1")) return;
  if(dynamic_cast<TH2*> (h) || h->InheritsFrom("TH2") ) { cerr << "ERROR    Will not Draw TH2 " << endl; return; }
  if(dynamic_cast<TTree*> (h) || h->InheritsFrom("TTree")) return; //not sure how this happens// perhaps unnecessary with first line

  bool glob_drawOpts = false;
  if(!drawOpts) drawOpts = new asciDrawOptions(h, conf.getColor("LogyThreshold"));
  else glob_drawOpts = true;
  //Create hClone from drawOpts
  h = extractHist(h, drawOpts);
  TH1* hClone = (TH1*) h->Clone();

  TString FillColor;  conf.getColor(FillColor,"FillColor");
  TString StatBox;  conf.getColor(StatBox,"StatBox");
  TString XAxis;  conf.getColor(XAxis,"XAxis");
  TString YAxis;  conf.getColor(YAxis,"YAxis");
  TString FillSymbol;  conf.getColor(FillSymbol,"FillSymbol");
  if (FillSymbol.IsNull()) FillSymbol = "XO";
  
  int binWidth = 1;
  if(h->GetNbinsX() + 16 < conf.m_nCols/2 - 1) binWidth = (conf.m_nCols)/(h->GetNbinsX()+16);

  TH1* hbinned = 0;
  map<int, std::vector<int> > asciHist; // map of value by bins
  if( h->GetNbinsX() + 16 > conf.m_nCols) {
    cout << "WARNING    Must Rebin hist cols " << conf.m_nCols << " Bins " << h->GetNbinsX() << endl;    
    if( h->GetNbinsX() + 16 < 2*conf.m_nCols) { hbinned = h->Rebin(2, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./2); } 
    else if( h->GetNbinsX() +16 < 3*conf.m_nCols) { hbinned = h->Rebin(3, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./3); }
    else if( h->GetNbinsX() +16 < 4*conf.m_nCols) { hbinned = h->Rebin(4, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./4); }
    else if( h->GetNbinsX() +16 < 5*conf.m_nCols) { hbinned = h->Rebin(5, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./5); }
    else if( h->GetNbinsX() +16 < 6*conf.m_nCols) { hbinned = h->Rebin(6, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./6); }
    else if( h->GetNbinsX() +16 < 7*conf.m_nCols) { hbinned = h->Rebin(7, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./7); }
    else if( h->GetNbinsX() +16 < 8*conf.m_nCols) { hbinned = h->Rebin(8, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./8); }
    else if( h->GetNbinsX() +16 < 9*conf.m_nCols) { hbinned = h->Rebin(9, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./9); }
    else if( h->GetNbinsX() +16 < 10*conf.m_nCols) { hbinned = h->Rebin(10, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(1./10); }
    else { hbinned = h->Rebin( (int) (h->GetNbinsX()+16)/conf.m_nCols, ((TString)h->GetName())+"_reBinned"); if(h->GetMaximum() < 1) hbinned ->Scale(conf.m_nCols/(h->GetNbinsX()+16)); }

    h = hbinned;
    if(h->GetNbinsX() + 16 > conf.m_nCols) {
      cerr << "ERROR    Could not rebing: " << h->GetName() << " cols " << conf.m_nCols << " Bins " << h->GetNbinsX() << endl;    

      delete hClone;
      delete hbinned;
      return;
    }
    delete hClone;
    hClone = (TH1*) h->Clone();
  }
  int nbins = h->GetNbinsX()*binWidth;
  nbins -= ((TString) h->GetTitle()).Length();
  nbins += 12;
  if(nbins >= 2 ) { printPad(1,nbins/2,"_"); cout << hOriginal->GetTitle(); printPad(1,nbins/2,"_"); cout << endl;}
  else cout << hOriginal->GetTitle() << endl;
  cout << StatBox << "                                                  Name    " << hOriginal->GetName() << conf.colorNorm << endl;
  cout << StatBox << "                                                  Entries " << hOriginal->GetEntries() << conf.colorNorm << endl;
  cout << StatBox << "                                                  Mean    " << hOriginal->GetMean() << conf.colorNorm << endl;
  cout << StatBox << "                                                  RMS     " << hOriginal->GetRMS() << conf.colorNorm << endl;

  // Implement StatBox here
//   h->Draw();//to get TPaveStats
//   TPaveStats *p = (TPaveStats*) h->GetListOfFunctions()->FindObject("stats");
//   TIter pItr(p->GetListOfLines()) ;
//   while (pItr()){
//     TString p_str = (TString) *pItr;
//   }


  double scale = (conf.m_nRows-14)/h->GetMaximum();
  h->Scale(scale) ;

  //Do Logy if max userDefined ratio non-Zero Min
  double ratio = conf.getColor("LogyThreshold");
  if( ratio < 1) ratio = 100;
  if( drawOpts->logy ) {
    cout << "Setting Logy " << hClone->GetMaximum() / hClone->GetMinimum(10) << " Logy CutOff " << ratio << endl;
    for( int ibin = 1; ibin != hClone->GetNbinsX()+1; ++ ibin){
      hClone->SetBinContent(ibin, std::log(hClone->GetBinContent(ibin))/std::log(10));
    }
    scale = std::abs((conf.m_nRows-14)/hClone->GetMaximum());
    hClone->Scale( scale );
    h = hClone;
  }
  
  for(int ibin = 1; ibin != h->GetNbinsX()+1; ++ibin){
    int value = (int) h->GetBinContent(ibin);
    map<int, std::vector<int> >::iterator mitr = asciHist.find(value);
    if( mitr != asciHist.end() ) {
      (*mitr).second.push_back(ibin);
    }
    else {
      vector<int> bins;
      bins.push_back(ibin);
      asciHist.insert(make_pair(value,bins));
    }
  }

  int counter = 0;
  std::set<int> filledBins;
  for(map<int, std::vector<int> >::reverse_iterator mitr = asciHist.rbegin(); mitr != asciHist.rend(); ++mitr,counter++){
    int content = (*mitr).first ;
    double value = content*1./scale;
    if( drawOpts->logy ) value = std::pow(10, value);

    if(counter%4==0) cout << YAxis << setprecision(3) << left << setw(9) << value << "|-- " ; 
    else cout << YAxis << "         |   " << conf.colorNorm;

    for(vector<int>::const_iterator itr = (*mitr).second.begin(); itr != (*mitr).second.end(); ++itr){
      filledBins.insert(*itr);
    }
    for(int i = 1; i != h->GetNbinsX(); ++i){
      if(filledBins.find(i)!=filledBins.end()) {
	if(FillSymbol.Length()>1) {
	  if(i%2==0) printPad(1, binWidth, FillSymbol(0,1), &FillColor);
	  else printPad(1, binWidth, FillSymbol(1,2), &FillColor);
	}
	else printPad(1, binWidth, FillSymbol(0,1), &FillColor);
      }
      else printPad(1, binWidth, " ");
    }
    cout << endl;
  }

  nbins = h->GetNbinsX()*binWidth;
  nbins -= ((TString) h->GetXaxis()->GetTitle()).Length();
  nbins += 12;
  if(nbins >= 2) {printPad(1,nbins/2,"_",&XAxis); cout << XAxis << h->GetXaxis()->GetTitle(); printPad(1,nbins/2,"_",&XAxis); cout << conf.colorNorm << endl;}
  else cout << XAxis << h->GetXaxis()->GetTitle() << conf.colorNorm << endl;

  cout << "             ";
  for( int ibin = 1; ibin != h->GetNbinsX()+1; ++ibin) {    
    if(h->GetNbinsX() > 10) {
      if(ibin%10==1) cout << XAxis << setprecision(3) << left << setw(10*binWidth) << h->GetBinLowEdge(ibin)+h->GetBinWidth(ibin)/2. ;
    }
    else cout << XAxis << setprecision(3) << left << setw(binWidth) << h->GetBinLowEdge(ibin)+h->GetBinWidth(ibin)/2. ;
  }
  cout << conf.colorNorm << endl << endl;
  if(conf.m_interactive) {
    bool doLoop = true;
    bool interactive = (!originalCall);
    cin.clear();
    while ( doLoop ) { 
      char myDec[100];
      if(!interactive) {
	cout << "Press 'return' to continue or 'i' for interactive mode " << endl;
	cin.getline(myDec,100);
	if( ((TString) myDec) == "i" ) interactive = true;
	else if( ((TString) myDec).Length()) {cout << "In valid Option Try again" << endl; continue; }
	else { cout << flush; break; }
      }
      if(interactive) {
	cout << "Interactive Mode" << endl;
	cout << "'x' for X-Axis" << endl;
	cout << "'y' for Y-Ayis" << endl;
	cout << "'l' for Toggle log scale " << drawOpts->logy << endl;
	cout << "'p' Print Current Hist Options" << endl;
	cout << "'d' Draw the hist again" << endl;
	cout << "'r' Reset Histogram To Original Size" << endl;
	cout << "'n' Go to next hist" << endl;
	cout << "'L' Lock Current Settings For Rest Of Program" << endl;
	cout << "'q' Quit Interactive Mode" << endl;
	cin >> myDec;
	TString opt = myDec;
	if( opt.Contains("p") ) drawOpts->print();
	if( opt == "n") { cin.getline(myDec,100); break; }
	if( opt == "q") { cin.getline(myDec,100); conf.m_interactive = false; break; } 
	if( opt == "L") { 
	  cin.getline(myDec,100); 
	  if(conf.m_drawOpt) delete conf.m_drawOpt;
	  conf.m_drawOpt = drawOpts->clone();
	}
	if( opt.Contains("x") ) getFloat("X", drawOpts->minx, drawOpts->maxx);
	if( opt.Contains("y") ) getFloat("Y", drawOpts->miny, drawOpts->maxy);
	if( opt.Contains("r") ) { delete drawOpts; drawOpts = new asciDrawOptions(hOriginal, conf.getColor("LogyThreshold")); }
	if( opt.Contains("l") ) {
	  if( drawOpts->logy ) drawOpts->logy = 0;
	  else drawOpts->logy = 1;
	}
	if( opt.Contains("d") ) { 
	  cin.getline(myDec, 100); //flush buffer
	  drawAsci( hOriginal, drawOpts, 0 ); 
	  interactive = false; doLoop = false; //ensure prompt exit w/ 'q' option
	  //delete allocated hists that'd other wise be lost
	  if(hbinned) delete hbinned;
	  delete hClone;
	  hbinned = 0;
	  hClone = 0;
	}
      }
    }
  }  
  if( originalCall ){
    //Do not continually delete these pointers
    if(hClone) delete hClone;
    if(hbinned) delete hbinned;
    if(!glob_drawOpts){
      delete drawOpts;
      drawOpts = 0;
    }
  }
}

void drawAsci2(TH2* h, asciDrawOptions* drawOpts, bool originalCall){
  TH2* hOriginal = h;

  bool glob_drawOpts = false;
  if(!drawOpts) drawOpts = new asciDrawOptions(h, conf.getColor("LogyThreshold"));
  else glob_drawOpts = true;
  //Create hClone from drawOpts
  h = extractHist2(h, drawOpts);
  TH2* hClone = (TH2*) h->Clone();

  TString FillColor;  conf.getColor(FillColor,"FillColor");
  TString StatBox;  conf.getColor(StatBox,"StatBox");
  TString XAxis;  conf.getColor(  XAxis,"XAxis");
  TString YAxis;  conf.getColor(  YAxis,"YAxis");

  if( h->GetNbinsX() + 16 > conf.m_nCols) {
    cerr << "ERROR     Cannot rebin 2D hist Exiting " << endl;
    return;
  }

  int binWidth = 1;
  if(h->GetNbinsX() + 16 < conf.m_nCols/2 - 1) binWidth = (conf.m_nCols)/(h->GetNbinsX()+25);
  TString FillSymbol;  conf.getColor(  FillSymbol,"FillSymbol");
  if (FillSymbol.IsNull() && conf.isTerminal) FillSymbol = "*";
  else if (!conf.isTerminal || !conf.m_doColor) FillSymbol = "";

  int nbins = h->GetNbinsX()*binWidth;
  nbins -= ((TString) h->GetTitle()).Length();
  nbins += 12;
  if(nbins >= 2 ) { printPad(1,nbins/2,"_"); cout << hOriginal->GetTitle(); printPad(1,nbins/2,"_"); cout << endl;}
  else cout << hOriginal->GetTitle() << endl;
  cout << StatBox << "                                                  Name    " << hOriginal->GetName() << conf.colorNorm << endl;
  cout << StatBox << "                                                  Entries " << hOriginal->GetEntries() << conf.colorNorm << endl;
  cout << StatBox << "                                                  Mean X  " << hOriginal->GetMean() << conf.colorNorm << endl;
  cout << StatBox << "                                                  RMS X   " << hOriginal->GetRMS() << conf.colorNorm << endl;
  cout << StatBox << "                                                  Mean Y  " << hOriginal->GetMean(2) << conf.colorNorm << endl;
  cout << StatBox << "                                                  RMS Y   " << hOriginal->GetRMS(2) << conf.colorNorm << endl;

  int nColors = conf.m_colors.size();
  
  map<int, int> invertedHist; //map of bin # by color index in colors array

  double maxVal = h->GetMaximum();
  double minVal = h->GetMinimum(0.000001);
  double binStep = (maxVal - minVal)/(nColors -1); // last color is null color
  bool logz = false;
  if( drawOpts->logz ) {
    cout << "Setting Logz" << endl;
    logz = true;
    for( int ibiny = 1; ibiny != h->GetNbinsY() +1; ++ibiny) {
      for ( int ibinx = 1; ibinx != h->GetNbinsX() +1; ++ibinx){
	int bin = h->GetBin( ibinx, ibiny);
	double contents = h->GetBinContent( ibinx, ibiny );
	if ( contents > 0.000001 ) h->SetBinContent(bin, std::log(contents)/std::log(10.));
	else h->SetBinContent(bin, 0);
      }
    }
    maxVal = h->GetMaximum();
    minVal = h->GetMinimum(0.00001);
    binStep = (maxVal - minVal)/(nColors -1); // last color is null color
  }

  for( int ibiny = 1; ibiny != h->GetNbinsY() +1; ++ibiny) {
    for ( int ibinx = 1; ibinx != h->GetNbinsX() +1; ++ibinx){
      int bin = h->GetBin( ibinx, ibiny);
      double contents = h->GetBinContent( ibinx, ibiny );
      //      cout << contents << endl;
      //determine color bin here
      int colorBin = 0;
      if( contents < minVal ) colorBin = nColors -1;
      else {
	colorBin = (int) ((contents-minVal)/binStep);  
	colorBin = nColors - colorBin - 1 - 1; //subtract null //subtract index starts at 0
      }
      //      cout << colorBin << " " << binStep << " " << contents << endl;
      if( colorBin >= nColors) colorBin = colorBin-1;
      if( colorBin < 0) colorBin = 0;
      invertedHist.insert(make_pair( bin, colorBin ));
    }
  }

  //Now draw!!!!
  int scaleBin = 0;
  int counter = 0;
  int scaleSize = h->GetNbinsY()/nColors;

  for( int ibiny = h->GetNbinsY(); ibiny != 0; --ibiny, ++counter) {
    if(ibiny%4==0) cout << YAxis << setprecision(3) << left << setw(9) << h->GetYaxis()->GetBinLowEdge(ibiny) << "|-- " << conf.colorNorm; 
    else cout << YAxis << "         |   " << conf.colorNorm;
    for ( int ibinx = 1; ibinx != h->GetNbinsX() +1; ++ibinx){
      int bin = h->GetBin( ibinx, ibiny);
      std::map< int, int>::const_iterator mitr = invertedHist.find(bin);
      if( mitr== invertedHist.end()) continue;
      if(ibiny%2==0 && FillSymbol.Length() == 2){
	if(ibinx%2==0) printPad(1, binWidth, FillSymbol(0,1), &conf.m_colors[(*mitr).second]);
	else printPad(1, binWidth, FillSymbol(1,1), &conf.m_colors[(*mitr).second]);
      }
      else if( FillSymbol.Length()==2){
	if(ibinx%2==1) printPad(1, binWidth, FillSymbol(0,1), &conf.m_colors[(*mitr).second]);
	else printPad(1, binWidth, FillSymbol(1,1), &conf.m_colors[(*mitr).second]);
      }
      else printPad(1, binWidth, FillSymbol(0,1), &conf.m_colors[(*mitr).second]);
    }
    //Draw Scale
    cout << "    ";
    if(counter == scaleSize) { scaleBin++; counter = 0;}
    if(counter == 0 && scaleBin < nColors) { 
      printPad(1, 3, FillSymbol(0,1), &conf.m_colors[scaleBin]); 
      if( scaleBin  == nColors -1) cout << " " << setprecision(2) << left << setw(10) << "0";	
      else if(logz) cout << " " << setprecision(2) << left << setw(10) << std::pow(10, maxVal - scaleBin*binStep);
      else cout << " " << setprecision(2) << left << setw(10) << maxVal - scaleBin*binStep;
    }
    else if ( scaleBin < nColors )  printPad(1, 3, FillSymbol(0,1), &conf.m_colors[scaleBin]); 
    cout << endl;
  }  

  cout << conf.colorNorm;// << endl;
  //
 
  nbins = h->GetNbinsX()*binWidth;
  nbins -= ((TString) h->GetXaxis()->GetTitle()).Length();
  nbins += 12;
  if(nbins >= 2) {printPad(1,nbins/2,"_",&XAxis); cout << XAxis << h->GetXaxis()->GetTitle(); printPad(1,nbins/2,"_",&XAxis); cout << endl;}
  else cout << XAxis << h->GetXaxis()->GetTitle() << conf.colorNorm << endl;

  cout << "            ";
  for( int ibin = 1; ibin != h->GetNbinsX()+1; ++ibin) {    
    if(h->GetNbinsX() > 10) {
      if(ibin%10==1) cout << XAxis << setprecision(3) << left << setw(10*binWidth) << h->GetBinLowEdge(ibin)+h->GetBinWidth(ibin)/2. ;
    }
    else cout << XAxis << setprecision(3) << left << setw(binWidth) << h->GetBinLowEdge(ibin)+h->GetBinWidth(ibin)/2. ;
  }
  cout << conf.colorNorm << endl << endl;

  if(conf.m_interactive) {
    bool doLoop = true;
    bool interactive = (!originalCall);
    cin.clear();
    while ( doLoop ) { 
      char myDec[100];
      if(!interactive) {
	cout << "Press 'return' to continue or 'i' for interactive mode " << endl;
	cin.getline(myDec,100);
	if( ((TString) myDec) == "i" ) interactive = true;
	else if( ((TString) myDec).Length()) {cout << "In valid Option Try again" << endl; continue; }
	else { cout << flush; break; }
      }
      if(interactive) {
	cout << "Interactive Mode" << endl;
	cout << "'x' for X-Axis" << endl;
	cout << "'y' for Y-Ayis" << endl;
	cout << "'z' for Z-Ayis" << endl;
	cout << "'l' for Toggle log scale " << drawOpts->logz << endl;
	cout << "'p' Print Current Hist Options" << endl;
	cout << "'d' Draw the hist again" << endl;
	cout << "'r' Reset Histogram To Original Size" << endl;
	cout << "'n' Go to next hist" << endl;
	cout << "'L' Lock Current Draw Options" << endl;
	cout << "'q' Quit Interactive Mode" << endl;
	cin >> myDec;
	TString opt = myDec;
	if( opt == "p") drawOpts->print();
	if( opt == "n") { cin.getline(myDec,100); break; }
	if( opt.Contains("x") ) getFloat("X", drawOpts->minx, drawOpts->maxx);
	if( opt.Contains("y") ) getFloat("Y", drawOpts->miny, drawOpts->maxy);
	if( opt.Contains("z") ) getFloat("Z", drawOpts->miny, drawOpts->maxz);
	if( opt.Contains("r") ) { delete drawOpts; drawOpts = new asciDrawOptions(hOriginal, conf.getColor("LogyThreshold")); }
	if( opt.Contains("l") ) {
	  if( drawOpts->logz ) drawOpts->logz = 0;
	  else drawOpts->logz = 1;
	}
	if( opt == "L") { 
	  cin.getline(myDec,100); 
	  if(conf.m_drawOpt2) delete conf.m_drawOpt2;
	  conf.m_drawOpt2 = drawOpts->clone();
	}
	if( opt == "q") {
	  cin.getline(myDec,100);
	  conf.m_interactive = false;
	  break;
	}
	if( opt.Contains("d") ) { 
	  cin.getline(myDec, 100); //flush buffer
	  drawAsci2( hOriginal, drawOpts, 0 ); 
	  interactive = false; doLoop = false; //ensure prompt exit w/ 'q' option
	  //delete allocated hists that'd other wise be lost
	  delete hClone;
	  hClone = 0;
	}
      }
    }
  }
  if( originalCall ){
    //Do not continually delete these pointers
    if(hClone) delete hClone;
    if(!glob_drawOpts) {
      delete drawOpts;
      drawOpts = 0;
    }
  }
}

void drawGraph(TGraph* g){

  double* x = g->GetX();
  double* y = g->GetY();

  if ( g->GetN() < 2) return;
  map< int , set<int> > invertedGraph;
  double maxX = 0;
  double minX = 0;
  double maxY = 0;
  double minY = 0;
  for(int i = 0; i != g->GetN(); ++i) {
    if( x[i] > maxX ) maxX = x[i];
    if( x[i] < minX ) minX = x[i];
    if( y[i] > maxY ) maxY = y[i];
    if( y[i] < minY ) minY = y[i];
  }
  double scale = conf.m_nRows/maxY;
  for(int i = 0; i != g->GetN(); ++i) y[i]/=scale;

  int binWidth = (int) (conf.m_nCols/(maxX + 25));
  cout << binWidth << " maxX " << maxX << " maxY " << maxY << " minX " << minX << " minY " << minY << endl;

  for(int i = 0; i!=g->GetN(); ++i){
    map<int, set<int> >::iterator itr = invertedGraph.find((int) y[i]);
    if( itr != invertedGraph.end() ) (*itr).second.insert((int) x[i]);
    else {
      set<int> vx; 
      vx.insert((int) x[i]);
      invertedGraph.insert(make_pair( (int) y[i], vx ));
    }
  }
  binWidth = 1;
  for(map<int, set<int> >::reverse_iterator itr = invertedGraph.rbegin(); itr != invertedGraph.rend(); ++itr){
    set<int> vx = (*itr).second;
    int currentCol = 0;
    for( set<int>::const_iterator sitr = vx.begin(); sitr != vx.end(); ++sitr ){
      int x_ = *sitr;
      printPad( 1, x_*binWidth - currentCol, " ");
      printPad( 1, (binWidth-1)/2, " ");
      printPad( 1, 1, "*");
      if ( (binWidth-1)%2==1 ) printPad( 1, binWidth/2 + 1, " ");
      else printPad( 1, binWidth/2, " ");
      currentCol = binWidth*x_;
    }
    cout << endl;
  }
}
