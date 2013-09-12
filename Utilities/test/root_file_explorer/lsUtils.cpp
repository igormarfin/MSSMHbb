#include "lsUtils.h"
#include "stringUtils.h"
#include <iostream>

#include <TH2.h>
#include <TKey.h>
#include <TH1.h>
#include <TTree.h>

using namespace std;

asciDrawOptions::asciDrawOptions(TH1* h, double ratio){
  minx = h->GetXaxis()->GetXmin();
  maxx = h->GetXaxis()->GetXmax();

  TH2* h2 = dynamic_cast<TH2*> (h);
  if(h2){
    miny = h->GetYaxis()->GetXmin();
    maxy = h->GetYaxis()->GetXmax();
    minz = h->GetMinimum(0.00001);
    maxz = h->GetMaximum();
  }
  else {
    miny = h->GetMinimum(0.00001);
    maxy = h->GetMaximum();
  }
  if(h2 && h2->GetMaximum() > 100 * h2->GetMinimum(0.000001)) { logz = true; logy = false; }
  else logz = false;
    
  if(!h2){
    if( ratio < 0.001 && ratio > -0.001 ) ratio = 100;
    if( h->GetMaximum() > ratio*h->GetMinimum(0.0001)) { cout << "Setting Logy " << ratio << " " << h->GetMaximum() << " " << h->GetMinimum(0.0001) << endl; logy = true; }
    else logy = false;
  }           
}

asciDrawOptions::asciDrawOptions(){}

asciDrawOptions* asciDrawOptions::clone(){
  asciDrawOptions* n = new asciDrawOptions();
  n->minx = this->minx;
  n->miny = this->miny;
  n->minz = this->minz;
  n->maxx = this->maxx;
  n->maxy = this->maxy;
  n->maxz = this->maxz;
  n->logy = this->logy;
  n->logz = this->logz; 
  return n;
}

void asciDrawOptions::print(){
    cout << "Min x " <<  minx << endl;
    cout << "Min y " <<  miny << endl;
    cout << "Min z " <<  minz << endl;
    cout << "Max x " <<  maxx << endl;
    cout << "Max y " <<  maxy << endl;
    cout << "Max z " <<  maxz << endl;
    cout << "log y " << logy << endl;
    cout << "log z " << logz << endl;
  }


JKey::JKey( TObject* o, TKey* k, TDirectory* d){obj = o; key = k; dir = d;}
JKey::~JKey() {  if(obj) delete obj; }


bool jkeyCompTime(const JKey* a, const JKey* b){
  TDatime aT = a->key->GetDatime();
  TDatime bT = b->key->GetDatime();
  return aT < bT;  
}

bool jkeyComp(const JKey* a, const JKey* b)
{
  TString aName = a->key->GetName();
  TString bName = b->key->GetName();
  return aName < bName;
}
bool jkeyCompR(const JKey* a, const JKey* b)
{
  TString aName = a->key->GetName();
  TString bName = b->key->GetName();
  return aName > bName;
}

bool jkeyCompSize(const JKey* a, const JKey* b){
  TTree* t_a = dynamic_cast<TTree*> (a->obj);
  TTree* t_b = dynamic_cast<TTree*> (b->obj);
  double m = a->key->GetNbytes();
  if(t_a) m = t_a->GetTotBytes();
  double n = b->key->GetNbytes();
  if(t_b) n = t_b->GetTotBytes();
  return m < n;
}
bool jkeyCompSizeR(const JKey* a, const JKey* b){
  TTree* t_a = dynamic_cast<TTree*> (a->obj);
  TTree* t_b = dynamic_cast<TTree*> (b->obj);
  double m = a->key->GetNbytes();
  if(t_a) m = t_a->GetTotBytes();
  double n = b->key->GetNbytes();
  if(t_b) n = t_b->GetTotBytes();
  return m > n;
}

bool jkeyCompN(const JKey* a, const JKey* b){
  TTree* t_a = dynamic_cast<TTree*> (a->obj);
  TTree* t_b = dynamic_cast<TTree*> (b->obj);
  TH1* h_a = dynamic_cast<TH1*> (a->obj);
  TH1* h_b = dynamic_cast<TH1*> (b->obj);

  double m = 0;
  if(t_a) m = t_a->GetEntries();
  if(h_a) m = h_a->GetEntries();

  double n = 0;
  if(t_b) n = t_b->GetEntries();
  if(h_b) n = h_b->GetEntries();
  return m < n;
}
bool jkeyCompNR(const JKey* a, const JKey* b){
  TTree* t_a = dynamic_cast<TTree*> (a->obj);
  TTree* t_b = dynamic_cast<TTree*> (b->obj);
  TH1* h_a = dynamic_cast<TH1*> (a->obj);
  TH1* h_b = dynamic_cast<TH1*> (b->obj);

  double m = 0;
  if(t_a) m = t_a->GetEntries();
  if(h_a) m = h_a->GetEntries();

  double n = 0;
  if(t_b) n = t_b->GetEntries();
  if(h_b) n = h_b->GetEntries();
  return m > n;
}

bool jkeyCompNBins(const JKey* a, const JKey* b){
  TH1* h_a = dynamic_cast<TH1*> (a->obj);
  TH1* h_b = dynamic_cast<TH1*> (b->obj);

  double m = 0;
  if(h_a) m = h_a->GetNbinsX()*h_a->GetNbinsY()*h_a->GetNbinsZ();

  double n = 0;
  if(h_b) n = h_b->GetNbinsX()*h_b->GetNbinsY()*h_b->GetNbinsZ();
  return m < n;
}
bool jkeyCompNBinsR(const JKey* a, const JKey* b){
  TH1* h_a = dynamic_cast<TH1*> (a->obj);
  TH1* h_b = dynamic_cast<TH1*> (b->obj);

  double m = 0;
  if(h_a) m = h_a->GetNbinsX()*h_a->GetNbinsY()*h_a->GetNbinsZ();

  double n = 0;
  if(h_b) n = h_b->GetNbinsX()*h_b->GetNbinsY()*h_b->GetNbinsZ();
  return m > n;
}

void getFloat(const TString& type, double &min, double &max){
  cout << type << "Min: " ;
  char myDec[100];
  cin.getline(myDec, 100);
  cin.getline(myDec, 100);
  TString opt = myDec;
  TString n1, n2;
  splitString(opt,n1,n2,' ');
  if(n1.Length() && !n1.IsFloat()) { cout << "Invalid Entry.  Hit 'return' for no change." << endl; getFloat(type, min, max); }
  if(!n2.Length()) {
    cout << type << "Max: " ;
    cin.getline(myDec, 100);
    n2 = myDec;
    if(n2.Length() && !n2.IsFloat()) { cout << "Invalid Entry. Hit 'return' for no change." << endl; getFloat(type, min, max);}
  }
  if( n1.Length() ) min = n1.Atof();
  if( n2.Length() ) max = n2.Atof();
}

TH1* extractHist(TH1* h, asciDrawOptions* opt){
  //ensure the x axis is in range
  if(opt->maxx > h->GetXaxis()->GetXmax()) opt->maxx = h->GetXaxis()->GetXmax();
  if(opt->minx < h->GetXaxis()->GetXmin()) opt->minx = h->GetXaxis()->GetXmin();
  //

  int xBin_low = h->FindBin(opt->minx+0.000001);
  if( xBin_low <= 0 ) xBin_low++;
  int xBin_high = h->FindBin(opt->maxx-0.000001);
  if( xBin_high >= h->GetNbinsX() + 1 ) xBin_high--;
  int nBins_new = xBin_high - xBin_low + 1;
  if( nBins_new <= 0 ) return (TH1*) h->Clone();
  TH1* hNew = 0;
  if(h->InheritsFrom("TH1D")) hNew = new TH1D(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nBins_new, h->GetBinLowEdge( xBin_low ), h->GetBinLowEdge(xBin_high) + h->GetBinWidth(xBin_high));
  else if(h->InheritsFrom("TH1I")) hNew = new TH1I(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nBins_new, h->GetBinLowEdge( xBin_low ), h->GetBinLowEdge(xBin_high) + h->GetBinWidth(xBin_high));
  else if(h->InheritsFrom("TH1S")) hNew = new TH1S(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nBins_new, h->GetBinLowEdge( xBin_low ), h->GetBinLowEdge(xBin_high) + h->GetBinWidth(xBin_high));
  else  hNew = new TH1F(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nBins_new, h->GetBinLowEdge( xBin_low ), h->GetBinLowEdge(xBin_high) + h->GetBinWidth(xBin_high));
  for( int ibin = 1 ; ibin <= nBins_new; ++ibin ) {
    double contents = h->GetBinContent(ibin+xBin_low-1);
    contents = std::min(contents, opt->maxy);
    contents = std::max(contents, opt->miny);
    hNew->SetBinContent(ibin, contents );
  }
  return hNew;
}

TH2* extractHist2(TH2* h, asciDrawOptions* opt){

  //ensure the x/y axis is in range
  if(opt->maxx > h->GetXaxis()->GetXmax()) opt->maxx = h->GetXaxis()->GetXmax();
  if(opt->minx < h->GetXaxis()->GetXmin()) opt->minx = h->GetXaxis()->GetXmin();
  if(opt->maxy > h->GetYaxis()->GetXmax()) opt->maxy = h->GetYaxis()->GetXmax();
  if(opt->miny < h->GetYaxis()->GetXmin()) opt->miny = h->GetYaxis()->GetXmin();
  //

  int xBin_low = h->GetXaxis()->FindBin(opt->minx+0.000001);
  if( xBin_low <= 0 ) xBin_low++;
  int xBin_high = h->GetXaxis()->FindBin(opt->maxx-0.00001);
  if( xBin_high >= h->GetNbinsX() + 1 ) xBin_high--;
  int nxBins_new = xBin_high - xBin_low + 1;
  if( nxBins_new <= 0 ) return (TH2*) h->Clone(); 

  int yBin_low = h->GetYaxis()->FindBin(opt->miny+0.00001);
  if( yBin_low <= 0 ) yBin_low++;
  int yBin_high = h->GetYaxis()->FindBin(opt->maxy-0.00001);
  if( yBin_high >= h->GetNbinsY() + 1 ) yBin_high--;
  int nyBins_new = yBin_high - yBin_low + 1;
  if( nyBins_new <= 0 ) return (TH2*) h->Clone(); 

  TH2* hNew = 0;
  if(h->InheritsFrom("TH2D")) hNew = new TH2D(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nxBins_new, 
					      h->GetXaxis()->GetBinLowEdge( xBin_low ), h->GetXaxis()->GetBinLowEdge(xBin_high) + h->GetXaxis()->GetBinWidth(xBin_high),
					      nyBins_new, h->GetYaxis()->GetBinLowEdge( yBin_low ), h->GetYaxis()->GetBinLowEdge(yBin_high) + h->GetYaxis()->GetBinWidth(yBin_high)
					      );
  else if(h->InheritsFrom("TH2I")) hNew = new TH2I(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nxBins_new, 
						   h->GetXaxis()->GetBinLowEdge( xBin_low ), h->GetXaxis()->GetBinLowEdge(xBin_high) + h->GetXaxis()->GetBinWidth(xBin_high),
						   nyBins_new, h->GetYaxis()->GetBinLowEdge( yBin_low ), h->GetYaxis()->GetBinLowEdge(yBin_high) + h->GetYaxis()->GetBinWidth(yBin_high)
					      );
  else if(h->InheritsFrom("TH2S")) hNew = new TH2S(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nxBins_new, 
						   h->GetXaxis()->GetBinLowEdge( xBin_low ), h->GetXaxis()->GetBinLowEdge(xBin_high) + h->GetXaxis()->GetBinWidth(xBin_high),
						   nyBins_new, h->GetYaxis()->GetBinLowEdge( yBin_low ), h->GetYaxis()->GetBinLowEdge(yBin_high) + h->GetYaxis()->GetBinWidth(yBin_high)
					      );
  else hNew = new TH2F(((TString) h->GetName())+"_temp", ((TString) h->GetTitle())+"_temp", nxBins_new, 
		       h->GetXaxis()->GetBinLowEdge( xBin_low ), h->GetXaxis()->GetBinLowEdge(xBin_high) + h->GetXaxis()->GetBinWidth(xBin_high),
		       nyBins_new, h->GetYaxis()->GetBinLowEdge( yBin_low ), h->GetYaxis()->GetBinLowEdge(yBin_high) + h->GetYaxis()->GetBinWidth(yBin_high)
					      );
  for( int ibinx = 1; ibinx <= hNew->GetNbinsX(); ++ibinx ){    
    for( int ibiny = 1; ibiny <= hNew->GetNbinsY(); ++ibiny ){
      double contents = h->GetBinContent(ibinx + xBin_low-1, ibiny + yBin_low-1);
      contents = std::min(opt->maxz, contents);
/*       contents = std::max(opt->minz, contents); */
      hNew->SetBinContent(ibinx, ibiny, contents);
    }
  }

  return hNew;

}
