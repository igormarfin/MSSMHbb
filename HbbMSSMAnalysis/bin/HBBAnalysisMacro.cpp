#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h" 
#include "TH1.h" 
#include "TH2.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TRFIOFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"

#include "Analysis/Utilities/interface/HBBTo4B.h"
#include "Analysis/Utilities/interface/Trigger.h"
#include "Analysis/Utilities/interface/svMass.h"
#include "Analysis/Utilities/interface/bTagEff.h"


int binNumber(float x, int nbins, float * bins) {

  int binN = 0;
  
  for (int iB=0; iB<nbins; ++iB) {
    binN = iB;
    if (x>=bins[iB]&&x<bins[iB+1]) 
      break;
  }
    
  return binN;

}

int main(int argc, char * argv[])
{
  // usage 
  // > HBBAnalysisMacro [filelist]
  // example filelist is Analysis/HbbMSSMAnalysis/test/Run2011A_PR4

  // Enable Sumw2 for all histograms so we can access the correct errors
  // from the histograms as well.
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  

  // requires as a parameter the name of the filelist

  std::ifstream inputfile(argv[1]);

  std::string filelist(argv[1]); 
  TString FileName(filelist);

  bool isData = false;

  if (FileName.Contains("Run2011")) 
    isData = true;

  TString ChainName("hbbanalysis/HBBTo4B");
  TString SVMassFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_svMass.root");
  TString SVMassOnlineFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_BTagHLTmatched_svMass.root");

  svMass * sv = new svMass(SVMassFileName);
  svMass * svOnline = new svMass(SVMassOnlineFileName);

  TString OutputFileName = TString(filelist) + TString(".root");

  float weight = 1;
  int numberOfFiles = 0;
  inputfile >> numberOfFiles;

  if (!isData) {

    int nGenEvents = 0;
    float internalXsec = 0;
    float externalXsecLO = 0;
    float externalXsecNLO = 0;
    float filterEfficiency = 0;
    for (int iFile=0; iFile<numberOfFiles; ++iFile) {
      TString filename;
      inputfile >> filename;
      std::cout << "File " << filename ;
      TFile * file = new TFile(filename);
      if (file->IsZombie()) {
	std::cout << " does not exist... quitting... " << std::endl;
      }
      std::cout << std::endl;
      TH1F * events = (TH1F*)file->Get("InputEvents/EventCount");
      nGenEvents += events->GetEntries(); 
      TTree * genInfo = (TTree*)file->Get("hbbanalysis/GenInfo");
      float internalX;
      float externalXLO;
      float externalXNLO;
      float filterEff; 
      genInfo->SetBranchAddress("internalXsec",&internalX);
      genInfo->SetBranchAddress("externalXsecLO",&externalXLO);
      genInfo->SetBranchAddress("externalXsecNLO",&externalXNLO);
      genInfo->SetBranchAddress("filterEfficiency",&filterEff);
      genInfo->GetEntry(0);
      internalXsec += internalX;
      externalXsecLO += externalXLO;
      externalXsecNLO += externalXNLO;
      filterEfficiency += filterEff; 
    }

    internalXsec /= float(numberOfFiles);
    externalXsecLO /= float(numberOfFiles);
    externalXsecNLO /= float(numberOfFiles);
    filterEfficiency /= float(numberOfFiles);

    std::cout << std::endl;
    std::cout << "   Total number of generated events = " << nGenEvents << std::endl;
    std::cout << "   InternalXsec      = " << internalXsec << std::endl;
    std::cout << "   ExternalXsecLO    = " << externalXsecLO << std::endl;
    std::cout << "   ExternalXsecNLO   = " << externalXsecNLO << std::endl;
    std::cout << "   Filter efficiency = " << filterEfficiency << std::endl;
    std::cout << std::endl;
    
    weight = internalXsec*fabs(filterEfficiency)/float(nGenEvents);

    // float weight = externalXsecLO*filterEfficiency/float(nGenEvents);
    // float weight = externalXsecNLO*filterEfficiency/float(nGenEvents);
    //  weight = 1;

  }

  if (isData)
    std::cout << "Running on data ..." << std::endl;
  else
    std::cout << "Running on Monte Carlo ..." << std::endl;

  std::cout << std::endl;
  std::cout << "Event weight for 1 pb-1 = " << weight << std::endl;
  std::cout << std::endl;


  // TRIGGER LIST
  std::vector<TString> triggers;
  triggers.clear();
  triggers.push_back("CentralJet46_BTagIP3D_CentralJet38_BTagIP3D");
  triggers.push_back("CentralJet46_CentralJet38_DiBTagIP3D");

  // BTAG EFF CLASS
  std::string flavlabel[3] = {"udsg","c","b"};
  // define ovl labels
  std::string sovl[4] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };


  // corresponding TCHP cuts
  float btcut[4]    = { 3.41, 6., 0.898, 2.0};
  float jetPtMin[3] = { 46,    38,   15};
  float jetPtMax[3] = {1000, 1000, 1000};
  float jetEtaMax = 2.6;
  int minJet = 3;

  std::string offlineBtag("offline");
  std::string onlineBtag("online");

  std::string btagOfflineEffFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt-N14/btagMatrixOffline-csv/plotmistag-b.root");
  std::string btagOnlineEffFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt_bEnriched-T14b/btagMatrixOnline-csv/smooth/btagMatrixOnline-csv-smpf.root");


  bTagEff * bTagEffOffline = new bTagEff(btagOfflineEffFile,offlineBtag,flavlabel,sovl);
  bTagEff * bTagRelEffOnline  = new bTagEff(btagOnlineEffFile,onlineBtag,flavlabel,sovl);

  TFile * OutputFile = new TFile(FileName+".root","recreate");
  OutputFile->cd("");

  float massBins = 5;
  float binsMass[6] = {0,100,200,300,400,500};

  TH1F * bTmassTCHPTH[5];
  TH1F * bTmassTCHP6H[5];
  TH1F * bTmassCSVTH[5];
  TH1F * bTmassSSVHPTH[5];

  TString massBinsTStr[5] = {"M100","M200","M300","M400","M500"};

  TH1F * dijetMassTCHPTH = new TH1F("dijetMassTCHPTH","",24,0.,600.);
  TH2F * dijetMassEvBTagTCHPTH = new TH2F("dijetMassEvBTagTCHPTH","",20,0.,600.,6,-0.5,5.5);
  TH1F * evBTagTCHPTH = new TH1F("evBTagTCHPTH","",6,-0.5,5.5);

  TH1F * dijetMassTCHP6H = new TH1F("dijetMassTCHP6H","",24,0.,600.);
  TH2F * dijetMassEvBTagTCHP6H = new TH2F("dijetMassEvBTagTCHP6H","",20,0.,600.,6,-0.5,5.5);
  TH1F * evBTagTCHP6H = new TH1F("evBTagTCHP6H","",6,-0.5,5.5);

  TH1F * dijetMassCSVTH = new TH1F("dijetMassCSVTH","",24,0.,600.);
  TH2F * dijetMassEvBTagCSVTH = new TH2F("dijetMassEvBTagCSVTH","",20,0.,600.,6,-0.5,5.5);
  TH1F * evBTagCSVTH = new TH1F("evBTagCSVTH","",6,-0.5,5.5);

  TH1F * dijetMassSSVHPTH = new TH1F("dijetMassSSVHPTH","",24,0.,600.);
  TH2F * dijetMassEvBTagSSVHPTH = new TH2F("dijetMassEvBTagSSVHPTH","",20,0.,600.,6,-0.5,5.5);
  TH1F * evBTagSSVHPTH = new TH1F("evBTagSSVHPTH","",6,-0.5,5.5);

  for (int imass=0;imass<massBins;++imass) {
    bTmassTCHPTH[imass] = new TH1F("bTag"+massBinsTStr[imass]+"TH","",6,-0.5,5.5);
    bTmassTCHP6H[imass] = new TH1F("bTag"+massBinsTStr[imass]+"6H","",6,-0.5,5.5);
    bTmassCSVTH[imass] = new TH1F("bTag"+massBinsTStr[imass]+"CSVTH","",6,-0.5,5.5);
    bTmassSSVHPTH[imass] = new TH1F("bTag"+massBinsTStr[imass]+"SSVHPTH","",6,-0.5,5.5);
  }

  TH1F * svMassFirstJetH = new TH1F("svMassFirstJetH","",20,0.,10.);
  TH1F * svMassFirstJetTCHPTH = new TH1F("svMassFirstJetTCHPTH","",20,0.,10.);
  TH1F * svMassFirstJetTCHP6H = new TH1F("svMassFirstJetTCHP6H","",20,0.,10.);
  TH1F * svMassFirstJetCSVTH = new TH1F("svMassFirstJetCSVTH","",20,0.,10.);
  TH1F * svMassFirstJetSSVHPTH = new TH1F("svMassFirstJetSSVHPTH","",20,0.,10.);

  TH1F * svMassSecondJetH = new TH1F("svMassSecondJetH","",20,0.,10.);
  TH1F * svMassSecondJetTCHPTH = new TH1F("svMassSecondJetTCHPTH","",20,0.,10.);
  TH1F * svMassSecondJetTCHP6H = new TH1F("svMassSecondJetTCHP6H","",20,0.,10.);
  TH1F * svMassSecondJetCSVTH = new TH1F("svMassSecondJetCSVTH","",20,0.,10.);
  TH1F * svMassSecondJetSSVHPTH = new TH1F("svMassSecondJetSSVHPTH","",20,0.,10.);

  TH1F * svMassThirdJetH = new TH1F("svMassThirdJetH","",20,0.,10.);
  TH1F * svMassThirdJetTCHPTH = new TH1F("svMassThirdJetTCHPTH","",20,0.,10.);
  TH1F * svMassThirdJetTCHP6H = new TH1F("svMassThirdJetTCHP6H","",20,0.,10.);
  TH1F * svMassThirdJetCSVTH = new TH1F("svMassThirdJetCSVTH","",20,0.,10.);
  TH1F * svMassThirdJetSSVHPTH = new TH1F("svMassThirdJetSSVHPTH","",20,0.,10.);

  TH1F * ptFirstJetH = new TH1F("ptFirstJetH","",30,0.,300.);
  TH1F * ptSecondJetH = new TH1F("ptSecondJetH","",30,0.,300.);
  TH1F * ptThirdJetH = new TH1F("ptThirdJetH","",30,0.,300.);

  // Templates ====>

  TString templateNames[3][3];

  // names [categ][flavor];
  // categ  = 0 (Xbb) 
  //        = 1 (bXb)
  //        = 2 (bbX)
  // flavor = 0 (udsg)
  //        = 1 (c)
  //        = 2 (b)  
  templateNames[0][0] = TString("Qbb");
  templateNames[0][1] = TString("Cbb");
  templateNames[0][2] = TString("Bbb");
  templateNames[1][0] = TString("bQb");
  templateNames[1][1] = TString("bCb");
  templateNames[1][2] = TString("bBb");
  templateNames[2][0] = TString("bbQ");
  templateNames[2][1] = TString("bbC");
  templateNames[2][2] = TString("bbB");

  // templateH[tagger][categ][flavor]
  //
  // tagger    = 0 : TCHPT
  //           = 1 : TCHP6
  //           = 2 : CSVT
  //           = 3 : SSVHPT
  //
  // categ     = 0 : Xbb
  //           = 1 : bXb
  //           = 2 : bbX
  // 
  // flavor    = 0 : udsg
  //           = 1 : c
  //           = 2 : b

  TH1F * bTagTemplateH[4][3][3];
  TH2F * massBTagTemplateH[4][3][3];

  TH1F * bTagOnlineCorrectedTemplateH[4][3][3];
  TH2F * massBTagOnlineCorrectedTemplateH[4][3][3];

  for (int isovl=0; isovl<4; ++isovl) {
    for (int iFlav=0; iFlav<3; ++iFlav) {
      for (int iCat=0; iCat<3; ++iCat) {
        bTagTemplateH[isovl][iCat][iFlav] = new TH1F("bTagTemplate"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",6,-0.5,5.5);
        massBTagTemplateH[isovl][iCat][iFlav] = new TH2F("MassBTagTemplate_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",20,0.,600.,6,-0.5,5.5);
	bTagOnlineCorrectedTemplateH[isovl][iCat][iFlav] = new TH1F("bTagTemplate_OnlCorr_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",6,-0.5,5.5); 
	massBTagOnlineCorrectedTemplateH[isovl][iCat][iFlav] = new TH2F("MassBTagTemplate_OnlCorr_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",20,0.,600.,6,-0.5,5.5);
      }
    }
  }

  //   tPatTemplateH[tagger][categ]
  TH1F * tPatternH[4][3];
  TString CategNames[3] = {"Xbb","bXb","bbX"};
  TString tPatNames[8] = {"xxx","Txx","xTx","xTT","xTT","TxT","TTx","TTT"};

  for (int iTagger=0; iTagger<4; ++iTagger) {
    for (int iCat=0; iCat<3; ++iCat) {
      tPatternH[iTagger][iCat] = new TH1F("tBattern_"+TString(sovl[iTagger])+"_"+CategNames[iCat]+"_H",CategNames[iCat]+" "+TString(sovl[iTagger]),8,-0.5,7.5);
      for (int iB=0; iB<8; ++iB) 
	tPatternH[iTagger][iCat]->GetXaxis()->SetBinLabel(iB+1,tPatNames[iB]);
    }
  }

  std::ifstream inputFile(argv[1]);
  inputFile >> numberOfFiles;

  int nEvents = 0;
  for (int iFile=0; iFile<numberOfFiles; ++iFile) {
    
    TString filename;
    inputFile >> filename;
    TFile * file = new TFile(filename);
    TTree * tree = (TTree*)file->Get(ChainName);
    int numberOfEntries = tree->GetEntries();
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;
    std::cout << filename << " : " << numberOfEntries << std::endl;

    HBBTo4B * hbbTo4b = new HBBTo4B();
    hbbTo4b->Init(tree);

    Trigger * trigger = new Trigger(file,triggers);
    trigger->PrintOut();

    for (int iE=0; iE<numberOfEntries; ++iE) {

      nEvents++;
      if (nEvents%100000==0)
    	std::cout << "Processed " << nEvents << " events " << std::endl;

      hbbTo4b->GetEntry(iE);

      int TriggerAccept = hbbTo4b->trgAccept;
      bool accept = trigger->accept(TriggerAccept);
      
      if (!accept) continue;

      int nJets = 0;
      int numberOfJets = hbbTo4b->NumberOfJets;

      numberOfJets = TMath::Min(numberOfJets,int(28));

      std::vector<int> leadingJets;
      leadingJets.clear();

      for (int iJ=0; iJ<numberOfJets; ++iJ) {
	
	if (nJets>=minJet) break;

	bool jetIDLoose = hbbTo4b->JetIDLoose[iJ];
	if (!jetIDLoose) continue;

	float absJetEta = TMath::Abs(hbbTo4b->JetEta[iJ]);
	if (absJetEta>jetEtaMax) continue;

	float jetPt = hbbTo4b->JetPt[iJ];   

	if (jetPt<jetPtMin[nJets]) continue;
	if (jetPt>jetPtMax[nJets]) continue;

	leadingJets.push_back(iJ);
	nJets++;


      }

      if (nJets<minJet) continue;

      int j1 = leadingJets[0];
      int j2 = leadingJets[1];
      int j3 = leadingJets[2];

      // Triggers

      bool bTagHlt1 = hbbTo4b->IsJetWithHltBtag[j1];
      bool bTagHlt2 = hbbTo4b->IsJetWithHltBtag[j2];
      bool bTagHlt3 = hbbTo4b->IsJetWithHltBtag[j3];

      bool tPatTTx = bTagHlt1 && bTagHlt2;
      bool tPatTxT = bTagHlt1 && bTagHlt3;
      bool tPatxTT = bTagHlt2 && bTagHlt3;
      bool tPatTTT = bTagHlt1 && bTagHlt2 && bTagHlt3;

      int thisTpat = -1;

      if (tPatTTT) {
	thisTpat = 3;
      }
      else {
	if (tPatxTT)
	  thisTpat = 0;
	if (tPatTxT)
	  thisTpat = 1;
	if (tPatTTx)
	  thisTpat = 2;
	
      }

      int tpat = 0;

      if (bTagHlt1 && !bTagHlt2 && !bTagHlt3)
	tpat = 1;
      else if (!bTagHlt1 && bTagHlt2 && !bTagHlt3)
	tpat = 2;
      else if (!bTagHlt1 && !bTagHlt2 && bTagHlt3)
	tpat = 3;
      else if (!bTagHlt1 && bTagHlt2 && bTagHlt3)
	tpat = 4;
      else if (bTagHlt1 && !bTagHlt2 && bTagHlt3)
	tpat = 5;
      else if (bTagHlt1 && bTagHlt2 && !bTagHlt3)
	tpat = 6;
      else if (bTagHlt1 && bTagHlt2 && bTagHlt3)
	tpat = 7;

      //      if (thisTpat<0) continue; //  

      // **********************************************
      // ***** Kinematics of three leading jets ======>
      // **********************************************

      TLorentzVector firstJet(hbbTo4b->JetPx[j1],
			      hbbTo4b->JetPy[j1],
			      hbbTo4b->JetPz[j1],
			      hbbTo4b->JetEnergy[j1]);

      TLorentzVector secondJet(hbbTo4b->JetPx[j2],
			       hbbTo4b->JetPy[j2],
			       hbbTo4b->JetPz[j2],
			       hbbTo4b->JetEnergy[j2]);
      

      TLorentzVector diJet = firstJet + secondJet;

      // invariant mass of teh two leading jets 
      float diJetMass = float(diJet.M());

      int jetIndex[3];
      jetIndex[0] = j1;
      jetIndex[1] = j2;
      jetIndex[2] = j3;


      float jetPtArray[3];
      float jetEtaArray[3];
      for (int index=0; index<3; ++index) {
	jetPtArray[index]  = hbbTo4b->JetPt[jetIndex[index]];
	jetEtaArray[index] = TMath::Abs(hbbTo4b->JetEta[jetIndex[index]]);
      }

      // TCHP tagger
      float btag1 = hbbTo4b->TCHPBJetTag[j1];
      float btag2 = hbbTo4b->TCHPBJetTag[j2];
      float btag3 = hbbTo4b->TCHPBJetTag[j3];

      // CSV tagger
      float dCsv1 = hbbTo4b->CombSVBJetTag[j1];
      float dCsv2 = hbbTo4b->CombSVBJetTag[j2];
      float dCsv3 = hbbTo4b->CombSVBJetTag[j3];

      // SSVHP tagger
      float dSsvhp1 = hbbTo4b->SVHPBJetTag[j1];
      float dSsvhp2 = hbbTo4b->SVHPBJetTag[j2];
      float dSsvhp3 = hbbTo4b->SVHPBJetTag[j3];

      float svMass1 = hbbTo4b->JetSvMass[j1];
      float svMass2 = hbbTo4b->JetSvMass[j2];
      float svMass3 = hbbTo4b->JetSvMass[j3];

      float pt1 = hbbTo4b->JetPt[j1];
      float pt2 = hbbTo4b->JetPt[j2];
      float pt3 = hbbTo4b->JetPt[j3];

      //      float eta1 = TMath::Abs(hbbTo4b->JetEta[j1]);
      //      float eta2 = TMath::Abs(hbbTo4b->JetEta[j2]);
      //      float eta3 = TMath::Abs(hbbTo4b->JetEta[j3]);

      if (svMass1<=0) svMass1 = 0.01;
      if (svMass2<=0) svMass2 = 0.01;
      if (svMass3<=0) svMass3 = 0.01;

      if (svMass1>=6) svMass1 = 5.99;
      if (svMass2>=6) svMass2 = 5.99;
      if (svMass3>=6) svMass3 = 5.99;

      // TCHPT tagger

      bool b1TCHPT = btag1>btcut[0];
      bool b2TCHPT = btag2>btcut[0];
      bool b3TCHPT = btag3>btcut[0];

      bool evBtagTCHPT = b1TCHPT && b2TCHPT && b3TCHPT;

      // TCHP6 tagger

      bool b1TCHP6 = btag1>btcut[1];
      bool b2TCHP6 = btag2>btcut[1];
      bool b3TCHP6 = btag3>btcut[1];

      bool evBtagTCHP6 = b1TCHP6 && b2TCHP6 && b3TCHP6;

      // CSVT tagger

      bool b1CSVT  = dCsv1>btcut[2];
      bool b2CSVT  = dCsv2>btcut[2];
      bool b3CSVT  = dCsv3>btcut[2];
      
      bool evBtagCSVT = b1CSVT && b2CSVT && b3CSVT;

      // SSVHPT tagger

      bool b1SSVHPT  = dSsvhp1>btcut[3];
      bool b2SSVHPT  = dSsvhp2>btcut[3];
      bool b3SSVHPT  = dSsvhp3>btcut[3];

      bool evBtagSSVHPT = b1SSVHPT && b2SSVHPT && b3SSVHPT;

      // bTagger[tagger][categ]; 
      // tagger = 0 : TCHPT
      //        = 1 : TCHP6
      //        = 2 : CSVT
      //        = 3 : SSVHPT
      // categ  = 0 : xBB;
      //        = 1 : BxB;
      //        = 2 : BBx;

      bool bTagger[4][3];
      bool jetsTag[4][3];

      jetsTag[0][0] = b1TCHPT;
      jetsTag[0][1] = b2TCHPT;
      jetsTag[0][2] = b3TCHPT;

      jetsTag[1][0] = b1TCHP6;
      jetsTag[1][1] = b2TCHP6;
      jetsTag[1][2] = b3TCHP6;

      jetsTag[2][0] = b1CSVT;
      jetsTag[2][1] = b2CSVT;
      jetsTag[2][2] = b3CSVT;

      jetsTag[2][0] = b1SSVHPT;
      jetsTag[2][1] = b2SSVHPT;
      jetsTag[2][2] = b3SSVHPT;

      bTagger[0][0] = b2TCHPT && b3TCHPT;
      bTagger[0][1] = b1TCHPT && b3TCHPT;
      bTagger[0][2] = b1TCHPT && b2TCHPT;

      bTagger[1][0] = b2TCHP6 && b3TCHP6;
      bTagger[1][1] = b1TCHP6 && b3TCHP6;
      bTagger[1][2] = b1TCHP6 && b2TCHP6;

      bTagger[2][0] = b2CSVT && b3CSVT;
      bTagger[2][1] = b1CSVT && b3CSVT;
      bTagger[2][2] = b1CSVT && b2CSVT;

      bTagger[3][0] = b2SSVHPT && b3SSVHPT;
      bTagger[3][1] = b1SSVHPT && b3SSVHPT;
      bTagger[3][2] = b1SSVHPT && b2SSVHPT;

      for (int iTagger=0; iTagger<4; ++iTagger) {
	for (int iCat=0; iCat<3; ++iCat) {
	  if (bTagger[iTagger][iCat]) {
	    tPatternH[iTagger][iCat]->Fill(float(tpat),weight);
	  }
	}
      }

      if (thisTpat<0) // select events where at least two out of 
	              // three leading jets are matched to BTagHlt object
	continue;

      // SV Mass array

      float svMassJets[3] ;
      svMassJets[0] = svMass1;
      svMassJets[1] = svMass2;
      svMassJets[2] = svMass3;

      int nEvBTag = sv->eventXBTag(svMassJets); 

      //      std::cout << "Pt 1 = " << pt1 << std::endl;
      //      std::cout << "Pt 2 = " << pt2 << std::endl;
      //      std::cout << "Pt 3 = " << pt3 << std::endl;

      //      std::cout << std::endl;
      //      std::cout << "d(CSV)   1 = " << dCsv1 << std::endl;
      //      std::cout << "d(CSV)   2 = " << dCsv2 << std::endl;
      //      std::cout << "d(CSV)   3 = " << dCsv3 << std::endl;
      //      std::cout << "CSVT       = " << btcut[2] << std::endl;
      //      std::cout << std::endl;
      //      std::cout << "d(SSVHP) 1 = " << dSsvhp1 << std::endl;
      //      std::cout << "d(SSVHP) 2 = " << dSsvhp2 << std::endl;
      //      std::cout << "d(SSVHP) 3 = " << dSsvhp3 << std::endl;
      //      std::cout << "SSVHPT     = " << btcut[3] << std::endl;
      //      std::cout << std::endl;

      // ********************************************
      // **** Filling some control distributions ****
      // ********************************************

      ptFirstJetH->Fill(pt1,weight);
      ptSecondJetH->Fill(pt2,weight);
      ptThirdJetH->Fill(pt3,weight);

      svMassFirstJetH->Fill(svMass1,weight);
      if (b1TCHPT)
	svMassFirstJetTCHPTH->Fill(svMass1,weight);
      if (b1TCHP6)
	svMassFirstJetTCHP6H->Fill(svMass1,weight);
      if (b1CSVT)
	svMassFirstJetCSVTH->Fill(svMass1,weight);
      if (b1SSVHPT)
        svMassFirstJetSSVHPTH->Fill(svMass1,weight);

      svMassSecondJetH->Fill(svMass2,weight);
      if (b2TCHPT)
	svMassSecondJetTCHPTH->Fill(svMass2,weight);
      if (b2TCHP6)
	svMassSecondJetTCHP6H->Fill(svMass2,weight);
      if (b2CSVT)
	svMassSecondJetCSVTH->Fill(svMass2,weight);
      if (b2SSVHPT)
	svMassSecondJetSSVHPTH->Fill(svMass2,weight);

      svMassThirdJetH->Fill(svMass3,weight);
      if (b3TCHPT)
	svMassThirdJetTCHPTH->Fill(svMass3,weight);
      if (b3TCHP6)
	svMassThirdJetTCHP6H->Fill(svMass3,weight);
      if (b3CSVT)
	svMassThirdJetCSVTH->Fill(svMass3,weight);
      if (b3SSVHPT)
	svMassThirdJetSSVHPTH->Fill(svMass3,weight);


      if (evBtagTCHPT) {
	dijetMassTCHPTH->Fill(diJetMass,weight);
	evBTagTCHPTH->Fill(float(nEvBTag),weight);
	dijetMassEvBTagTCHPTH->Fill(diJetMass,float(nEvBTag),weight);
	int massBin = sv->binNumber(diJetMass,massBins,binsMass);
	bTmassTCHPTH[massBin]->Fill(float(nEvBTag),weight);
      }

      if (evBtagTCHP6) {
	dijetMassTCHP6H->Fill(diJetMass,weight);
	evBTagTCHP6H->Fill(float(nEvBTag),weight);
	dijetMassEvBTagTCHP6H->Fill(diJetMass,float(nEvBTag),weight);
	int massBin = sv->binNumber(diJetMass,massBins,binsMass);
	bTmassTCHP6H[massBin]->Fill(float(nEvBTag),weight);
      }

      if (evBtagCSVT) {
	//	std::cout << "Wow! Triple-tag with CSVT tagger..." << std::endl;
	dijetMassCSVTH->Fill(diJetMass,weight);
	evBTagCSVTH->Fill(float(nEvBTag),weight);
	dijetMassEvBTagCSVTH->Fill(diJetMass,float(nEvBTag),weight);
	int massBin = sv->binNumber(diJetMass,massBins,binsMass);
        bTmassTCHP6H[massBin]->Fill(float(nEvBTag),weight);
      }

      if (evBtagSSVHPT) {
	//	std::cout << "Wow! Triple-tag with SSVHPT tagger..." << std::endl;
        dijetMassSSVHPTH->Fill(diJetMass,weight);
        evBTagSSVHPTH->Fill(float(nEvBTag),weight);
        dijetMassEvBTagSSVHPTH->Fill(diJetMass,float(nEvBTag),weight);
        int massBin = sv->binNumber(diJetMass,massBins,binsMass);
        bTmassSSVHPTH[massBin]->Fill(float(nEvBTag),weight);
      }


      // Creating templates ------>

      int nB1 = sv->getBinFromSvMass(svMass1);
      int nB2 = sv->getBinFromSvMass(svMass2);
      int nB3 = sv->getBinFromSvMass(svMass3);

      //      std::cout << "isvMass 1 = " << nB1 << std::endl;
      //      std::cout << "isvMass 2 = " << nB2 << std::endl;
      //      std::cout << "isvMass 3 = " << nB3 << std::endl;

      float probs[3];
      float probsOnline[3];
      float eff[3];
      float relEffOnl[3];

      int nBTags[3];
      nBTags[0] = nB1;
      nBTags[1] = nB2;
      nBTags[2] = nB3;

      // Templates

      for (int iTagger=0; iTagger<4; ++iTagger) {
	std::vector<float*> probsTpat;
	probsTpat.clear();
	for (int itpat=0; itpat<3; ++itpat) {
	  float probsP[3];
	  svOnline->getSVbins(2,iTagger,jetPtArray[itpat],jetEtaArray[itpat],probsP);
	  probsTpat.push_back(probsP);
	}
	for (int iCat=0; iCat<3; ++iCat) {
	  if (thisTpat==iCat || thisTpat == 3) { // online btag pattern must match offline
	    if (bTagger[iTagger][iCat]) {
	      for (int iFlv=0; iFlv<3; ++iFlv) {
		eff[iFlv] = bTagEffOffline->eff(iFlv,sovl[iTagger],jetPtArray[iCat],jetEtaArray[iCat]);
		relEffOnl[iFlv] = bTagRelEffOnline->eff(iFlv,sovl[iTagger],jetPtArray[iCat],jetEtaArray[iCat]);
		sv->getSVbins(iFlv,iTagger,jetPtArray[iCat],jetEtaArray[iCat],probs);
		svOnline->getSVbins(iFlv,iTagger,jetPtArray[iCat],jetEtaArray[iCat],probsOnline);
		int nTags[3];
		for (int iC=0; iC<3; iC++)
		  nTags[iC] = nBTags[iC];
		for (int iTag=0; iTag<3; ++iTag) {
		  nTags[iCat] = iTag;
		  int evtBTag = sv->eventBTag(nTags);
		  float templateWeight = probs[iTag]*weight*eff[iFlv];
		  bTagTemplateH[iTagger][iCat][iFlv]->Fill(evtBTag,templateWeight);
		  massBTagTemplateH[iTagger][iCat][iFlv]->Fill(diJetMass,evtBTag,templateWeight);
		  float addWeight = 0;
		  for (int itpat=0; itpat<3; ++itpat) {
		    float wtpat = 1;
		    if ((itpat != iCat) && (thisTpat != 3)) { // weight correction if offline 
		                                              // & online btag pattern are different (except TTT)
		      float effOnlineTpat = bTagRelEffOnline->eff(2,sovl[iTagger],jetPtArray[itpat],jetEtaArray[itpat]);
		      if (relEffOnl[iFlv]>0&&effOnlineTpat>0)
			wtpat = relEffOnl[iFlv] / effOnlineTpat;
		      
		      if (probsOnline[iTag]>0&&probsTpat.at(itpat)[iTag]>0)
			wtpat *= probsOnline[iTag]/probsTpat.at(itpat)[iTag];
		      
		    }
		    addWeight += wtpat;
		  }
		  bTagOnlineCorrectedTemplateH[iTagger][iCat][iFlv]->Fill(evtBTag,templateWeight * addWeight);
		  massBTagOnlineCorrectedTemplateH[iTagger][iCat][iFlv]->Fill(diJetMass,evtBTag,templateWeight * addWeight);
		}
	      }
	    }
	  }
	}
      }

    }

    delete file;
    delete hbbTo4b;
    delete trigger;
  }

  OutputFile->cd("");
  OutputFile->Write();
  OutputFile->Close();

}
