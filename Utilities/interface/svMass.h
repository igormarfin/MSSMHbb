#include <TFile.h>
#include <TH1F.h>
#include <TString.h>

class svMass {

 public:

  // ***********************
  // CONSTRUCTOR
  // reads histograms, incapsulating
  // probabilities of udsg-, c- and b-jets
  // to fall into certain bin in SV Mass distribution
  // as a function of jet Pt
  // ***********************
  svMass(TString fileName) {
   
    TString massSv("MassSV");
    TString flavorStringVector[4] = {"Qjet",
				     "Cjet",
				     "Bjet",
				     "Gjet"};

    TString binNumberString[3] = {"Bin1","Bin2","Bin3"};

    // default values (re-written later)
    _nPtBins = 20;
    _nPt2Dbins = 20;
    _nEta2Dbins = 20;
    _nSVBins = 3;
    
    TFile * file = new TFile(fileName);

    TH1F * svBinsHist = (TH1F*)file->Get("massBins");

    _nSVBins = svBinsHist->GetNbinsX();

    for (int bin=0; bin<=_nSVBins; ++bin)
      _svBins[bin] = svBinsHist->GetBinLowEdge(bin+1);

    for (int iFlv=0; iFlv<3; ++iFlv) {
      for (int iBin=0; iBin<_nSVBins; ++iBin) {
	massSvPtBinsTchptH[iFlv][iBin] = (TH1F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_Tchpt_H");
	massSvPtBinsTchp6H[iFlv][iBin] = (TH1F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_Tchp6_H");
	massSvPtBinsCSVTH[iFlv][iBin] = (TH1F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_CSVT_H");
	massSvPtBinsSSVHPTH[iFlv][iBin] = (TH1F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_SSVHPT_H");
	massSvPtEtaBinsTchptH[iFlv][iBin] = (TH2F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_Eta_Tchpt_H");
	massSvPtEtaBinsTchp6H[iFlv][iBin] = (TH2F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_Eta_Tchp6_H");
	massSvPtEtaBinsCSVTH[iFlv][iBin] = (TH2F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_Eta_CSVT_H");
	massSvPtEtaBinsSSVHPTH[iFlv][iBin] = (TH2F*)file->Get(massSv+flavorStringVector[iFlv]+binNumberString[iBin]+"_Pt_Eta_SSVHPT_H");
      }
    }

    _nPtBins = massSvPtBinsTchptH[0][0]->GetNbinsX();
    for (int bin=0; bin<=_nPtBins; ++bin)
      _ptBins[bin] = massSvPtBinsTchptH[0][0]->GetBinLowEdge(bin+1);

    _nPt2Dbins = massSvPtEtaBinsTchptH[0][0]->GetNbinsX();
    for (int bin=0; bin<=_nPt2Dbins; ++bin)
      _pt2Dbins[bin] = massSvPtBinsTchptH[0][0]->GetXaxis()->GetBinLowEdge(bin+1);
    
    _nEta2Dbins = massSvPtEtaBinsTchptH[0][0]->GetNbinsY();
    for (int bin=0; bin<=_nEta2Dbins; ++bin)
      _eta2Dbins[bin] = massSvPtBinsTchptH[0][0]->GetYaxis()->GetBinLowEdge(bin+1);
    
    

  }
  int binNumber(float x, int nbins, float * bins) {

    int binN = 0;
  
    if (x<=bins[0])
      return 0;
    if (x>=bins[nbins])
      return nbins - 1;

    for (int iB=0; iB<nbins; ++iB) {
      binN = iB;
      if (x>=bins[iB]&&x<bins[iB+1]) 
      break;
    }
    
    return binN;

  }

  int getNSVbins(float * bins) {
    for (int iB=0; iB<=_nSVBins; ++iB)
      bins[iB] = _svBins[iB];
    return _nSVBins;
  }

  int getPtBins(float * bins) {
    for (int iB=0; iB<=_nPtBins; ++iB)
      bins[iB] = _ptBins[iB];
    return _nPtBins;
  }

  void getSVbins(int iFlav, int bTagCut, float pt, float * probs){
    // iFlav = 0 - udsg
    // iFlav = 1 - c
    // iFlav > 1 - b 
    // bTagCut = 0 - TCHPT
    // bTagCut = 1 - TCHP6
    // bTagCut = 2 - CSVT
    // bTagCut > 2 - SSVHPT
    // pt - jet Pt 
    // output : probs[3] - returned probabilities 
    //          for jet to fall into one of SV Mass bins: 0, 1 or 2

    if (pt>_ptBins[_nPtBins]) pt = _ptBins[_nPtBins] - 1;

    for (int iBin=0; iBin<3; ++iBin) {
      if (bTagCut==0) // TCHPT
	probs[iBin] = massSvPtBinsTchptH[iFlav][iBin]->Interpolate(pt);
      else if (bTagCut==1)// TCHP6
	probs[iBin] = massSvPtBinsTchp6H[iFlav][iBin]->Interpolate(pt);
      else if (bTagCut==2) // CSVT
	probs[iBin] = massSvPtBinsCSVTH[iFlav][iBin]->Interpolate(pt);
      else 
	probs[iBin] = massSvPtBinsSSVHPTH[iFlav][iBin]->Interpolate(pt);
    }

  }

  void getSVbins(int iFlav, int bTagCut, float pt, float eta, float * probs){
    // iFlav = 0 - udsg
    // iFlav = 1 - c
    // iFlav > 1 - b 
    // bTagCut = 0 - TCHPT
    // bTagCut = 1 - TCHP6
    // bTagCut = 2 - CSVT
    // bTagCut > 2 - SSVHPT
    // pt - jet Pt 
    // eta - jet Eta
    // output : probs[3] - returned probabilities 
    //          for jet to fall into one of SV Mass bins: 0, 1 or 2

    //    std::cout << "Pt = " << pt << "  Eta = " << eta <<std::endl;

    if (pt>_pt2Dbins[_nPt2Dbins]) pt = _pt2Dbins[_nPt2Dbins] - 1;
    if (eta>_eta2Dbins[_nEta2Dbins]) eta = _eta2Dbins[_nEta2Dbins] - 0.01;
    

    for (int iBin=0; iBin<3; ++iBin) {

      if (bTagCut==0) // TCHPT
	probs[iBin] = massSvPtEtaBinsTchptH[iFlav][iBin]->Interpolate(pt,eta);
      else if (bTagCut==1) // TCHP6
	probs[iBin] = massSvPtEtaBinsTchp6H[iFlav][iBin]->Interpolate(pt,eta);
      else if (bTagCut==2) // CVST
	probs[iBin] = massSvPtEtaBinsCSVTH[iFlav][iBin]->Interpolate(pt,eta);
      else 
	probs[iBin] = massSvPtEtaBinsSSVHPTH[iFlav][iBin]->Interpolate(pt,eta);

    }

  }

  int eventXBTag(float * svMassBinJets) {
    int nSVMassBins[3];
    for (int iB=0; iB<3; ++iB) {
      nSVMassBins[iB] = binNumber(svMassBinJets[iB],_nSVBins,_svBins);
    }
    int evBTag = eventBTag(nSVMassBins);
    return evBTag;
  }

  int getBinFromSvMass(float svMassValue) {
    int bin = binNumber(svMassValue,_nSVBins,_svBins);
    return bin;
  }

  int eventBTag(int * svMassBinJets) {

    int sum12 = svMassBinJets[0] + svMassBinJets[1];  
    int i12 = 0;
    int i3  = 0;

    if (sum12<2) 
      i12 = 0;
    else if (sum12<3) 
      i12 = 1;
    else 
      i12 = 2;

    if (svMassBinJets[2]<2)
      i3 = 0;
    else
      i3 = 3;

    return i12 + i3; 

 }

 private:

  TH1F * massSvPtBinsTchptH[3][3];
  TH1F * massSvPtBinsTchp6H[3][3];
  TH1F * massSvPtBinsCSVTH[3][3];
  TH1F * massSvPtBinsSSVHPTH[3][3];

  TH2F * massSvPtEtaBinsTchptH[3][3];
  TH2F * massSvPtEtaBinsTchp6H[3][3];
  TH2F * massSvPtEtaBinsCSVTH[3][3];
  TH2F * massSvPtEtaBinsSSVHPTH[3][3];

  int _nSVBins;
  float _svBins[10];

  int _nPtBins;
  float _ptBins[20];

  int _nPt2Dbins;
  float _pt2Dbins[20];

  int _nEta2Dbins;
  float _eta2Dbins[20];

};
