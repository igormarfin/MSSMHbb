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

#include "Analysis/Utilities/interface/HBBTo4B_MC.h"
#include "Analysis/Utilities/interface/Trigger.h"
#include "Analysis/Utilities/interface/svMass.h"
#include "Analysis/Utilities/interface/bTagEff.h"
#include "Analysis/Utilities/interface/BBPurity.h"
#include "Analysis/Utilities/interface/HbbSyst.h"

float dPhiFrom2P(float Px1, float Py1,
		 float Px2, float Py2) {
  

  float prod = Px1*Px2 + Py1*Py2;
  float mod1 = TMath::Sqrt(Px1*Px1+Py1*Py1);
  float mod2 = TMath::Sqrt(Px2*Px2+Py2*Py2);
  
  float cosDPhi = prod/(mod1*mod2);

  return TMath::ACos(cosDPhi);
  
}

float dPhiFromPhi(float Phi1, 
		  float Phi2) {

  
  float Px1 = TMath::Cos(Phi1);
  float Py1 = TMath::Sin(Phi1);

  float Px2 = TMath::Cos(Phi2);
  float Py2 = TMath::Sin(Phi2);

  float dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);

  return dPhi;

}

float deltaR(float Eta1, float Phi1,
	     float Eta2, float Phi2) {

  float Px1 = TMath::Cos(Phi1);
  float Py1 = TMath::Sin(Phi1);

  float Px2 = TMath::Cos(Phi2);
  float Py2 = TMath::Sin(Phi2);

  float dPhi = dPhiFrom2P(Px1,Py1,Px2,Py2);
  float dEta = Eta1 - Eta2;

  float dR = TMath::Sqrt(dPhi*dPhi+dEta*dEta);

  return dR;

}


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

  // ****************************
  // some control variables ---->
  // ****************************
  bool printOut = false;
  bool applyTrigger = true;
  bool applyHltBTagMatch = true;
  bool applyWeightOnTriple = false;
  bool applyWeightOnDouble = false;
  bool applyBTagHltWeightOnDouble = false;
  bool applyBTagHltWeightOnTriple = false;

  // bH->bbb analysis
  // ****************
  // steering cards
  // after changing these steering cards
  // program should be recompiled
  // please do not set applyJESShift
  // and applyBTag to true simultenously
  // one should set applyJESShift = true in order to derive
  // templates corresponding to systematic shifts 
  // in JES
  // one should set applyBTagShift = true in order
  // to derive templates corresponding to systematic shifts
  // in b-tag efficiency / mistag rate

  bool applyJESShift = false; // apply systematic shift in JES (uncertainty)
  float jesShift = 2; // JES shift value is sigmas (+2sigma)
  // float jesUncertainty = -2; // JES shift value in sigmas (-2sigma)

  bool applyBTagSF = true; // apply Data/MC scale factors for b-tagging
  bool applyBTagShift = true; // apply systematic shift in BTag SF (uncertainty)  

  bool btagSFlighShift = 2; // light flavor mistag shift in sigmas (+2sigma)
  // float btagSFlightUncertainty = -2; // light flavor mistag shift in sigmas (-2sigma)
  bool btagSFbcShift = 2; // b- and c-jet tagging efficiency shift in sigmas (+2sigma)
  // float btagSFbcShift = 2; // b- and c-jet tagging efficiency shift in sigmas (+2sigma)

  // integrated luminosity
  float intLumi = 1000; // 1 fb-1

  bool jesUp = false;
  bool btagLFUp = false;
  bool btagBCUp = false;

  if (jesShift>0)
    jesUp = true;
  if (btagSFlighShift>0)
    btagLFUp = true;
  if (btagSFbcShift>0)
    btagBCUp = true;


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
  TString BBPurityFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/LFQCD_bbPurity.root");

  svMass * sv = new svMass(SVMassFileName);
  svMass * svOnline = new svMass(SVMassOnlineFileName);
  BBPurity * bbpur = new BBPurity(BBPurityFileName);
  HbbSyst * hbbSyst = new HbbSyst();

  TString OutputFileName = TString(filelist) + TString(".root");

  float eventWeight = 1;
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
    
    eventWeight = internalXsec*fabs(filterEfficiency)/float(nGenEvents);

    // float weight = externalXsecLO*filterEfficiency/float(nGenEvents);
    // float weight = externalXsecNLO*filterEfficiency/float(nGenEvents);
    //  weight = 1;

    std::cout << std::endl;
    std::cout << "Event weight for 1 fb-1 = " << eventWeight << std::endl;
    std::cout << "Integrated luminosity   = " << 1e-3*float(nGenEvents)/(internalXsec*fabs(filterEfficiency)) << " fb-1" << std::endl;
    std::cout << std::endl;

    eventWeight *= intLumi; 

  }

  if (isData)
    std::cout << "Running on data ..." << std::endl;
  else
    std::cout << "Running on Monte Carlo ..." << std::endl;

  // exit(1);

  // TRIGGER LIST
  std::vector<TString> triggers;
  triggers.clear();

  // Jet Pt triggers ----->

  // auxiliarly triggers
  //  triggers.push_back("HLT_L1DoubleJet36Central");

  triggers.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D"); // only light flavor
  //  triggers.push_back("HLT_bbPhi_CentralJet60_CentralJet53_v1"); // only light flavor

  // Single BTag triggers --->
  //  triggers.push_back("CentralJet46_BTagIP3D_CentralJet38_BTagIP3D");

  // bbPhi triggers ---->
  //  triggers.push_back("HLT_bbPhi_CentralJet46_CentralJet38_DiBTagIP3D_L25MinTag3_v1"); // light flavor QCD
  //  triggers.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D_v2"); // b-enriched QCD               
  //  triggers.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D"); // Alpgen

  // BTAG EFF CLASS
  std::string flavlabel[3] = {"udsg","c","b"};
  // define ovl labels
  std::string sovl[4] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };

  TString flavors[3] = {"Qquark","Cquark","Bquark"};

  // corresponding TCHP cuts
  float btcut[4]    = { 3.41, 6., 0.898, 2.0};
  float jetPtMin[3] = { 46,    38,   20};
  float jetPtMax[3] = {1000, 1000, 1000};
  float jetEtaMax = 2.2;
  int minJet = 3;

  std::string offlineBtag("offline");
  std::string onlineBtag("online");

  std::string btagOfflineEffFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt-N14/btagMatrixOffline-csv/plotmistag-b.root");
  std::string btagOnlineEffFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt_bEnriched-T14b/btagMatrixOnline-csv/smooth/btagMatrixOnline-csv-smpf.root");


  bTagEff * bTagEffOffline = new bTagEff(btagOfflineEffFile,offlineBtag,flavlabel,sovl);
  bTagEff * bTagRelEffOnline  = new bTagEff(btagOnlineEffFile,onlineBtag,flavlabel,sovl);

  TFile * OutputFile = new TFile(FileName+".root","recreate");
  OutputFile->cd("");


  TH1F * diJetMassBeforeTriggerH = new TH1F("diJetMassBeforeTriggerH","",120,0.,600.);
  TH1F * diJetMassAfterTriggerH = new TH1F("diJetMassAfterTriggerH","",120,0.,600.);


  TString tripleBsamples[6] = {"bbb",
			       "bbc",
			       "bbq",
			       "bcb",
			       "bqb",
			       "nonbb"};

  TH1F * dijetMassSamplesH[4][6];
  TH1F * evBTagSamplesH[4][6];
  TH2F * dijetMassEvBTagSamplesH[4][6];

  TH1F * deltaR12SamplesH[4][6];
  TH1F * deltaR13SamplesH[4][6];
  TH1F * deltaR23SamplesH[4][6];

  TH1F * deltaEta12SamplesH[4][6];
  TH1F * deltaEta13SamplesH[4][6];
  TH1F * deltaEta23SamplesH[4][6];

  TH1F * deltaPhi12SamplesH[4][6];
  TH1F * deltaPhi13SamplesH[4][6];
  TH1F * deltaPhi23SamplesH[4][6];

  TH1F * ptFirstSamplesH[4][6];
  TH1F * ptSecondSamplesH[4][6];
  TH1F * ptThirdSamplesH[4][6];

  TH1F * etaFirstSamplesH[4][6];
  TH1F * etaSecondSamplesH[4][6];
  TH1F * etaThirdSamplesH[4][6];

  TH1F * nJets20InTripleSamplesH[4][6];
  TH1F * nJets30InTripleSamplesH[4][6];

  TH1F * dijetMassH[4];
  TH1F * evBTagH[4];
  TH2F * dijetMassEvBTagH[4];

  TH1F * deltaR12H[4];
  TH1F * deltaR13H[4];
  TH1F * deltaR23H[4];

  TH1F * deltaEta12H[4];
  TH1F * deltaEta13H[4];
  TH1F * deltaEta23H[4];

  TH1F * deltaPhi12H[4];
  TH1F * deltaPhi13H[4];
  TH1F * deltaPhi23H[4];

  TH1F * ptFirstH[4];
  TH1F * ptSecondH[4];
  TH1F * ptThirdH[4];

  TH1F * etaFirstH[4];
  TH1F * etaSecondH[4];
  TH1F * etaThirdH[4];

  TH1F * nJets20InTripleH[4];
  TH1F * nJets30InTripleH[4];
  
  for (int iTagger=0; iTagger<4; ++iTagger) {

    TString TaggerName(sovl[iTagger]);

    dijetMassH[iTagger] = new TH1F("dijetMass"+TaggerName+"H","",30,0.,600.);
    evBTagH[iTagger] = new TH1F("evBTag"+TaggerName+"H","",6,-0.5,5.5);
    dijetMassEvBTagH[iTagger] = new TH2F("dijetMassEvBTag"+TaggerName+"H","",20,0.,600.,6,-0.5,5.5);

    deltaR12H[iTagger] = new TH1F("deltaR12"+TaggerName+"H","",100,0.,5.); 
    deltaR13H[iTagger] = new TH1F("deltaR13"+TaggerName+"H","",100,0.,5.); 
    deltaR23H[iTagger] = new TH1F("deltaR23"+TaggerName+"H","",100,0.,5.); 

    deltaEta12H[iTagger] = new TH1F("deltaEta12"+TaggerName+"H","",100,0.,5.);
    deltaEta13H[iTagger] = new TH1F("deltaEta13"+TaggerName+"H","",100,0.,5.);
    deltaEta23H[iTagger] = new TH1F("deltaEta23"+TaggerName+"H","",100,0.,5.);

    deltaPhi12H[iTagger] = new TH1F("deltaPhi12"+TaggerName+"H","",100,0.,TMath::Pi());
    deltaPhi13H[iTagger] = new TH1F("deltaPhi13"+TaggerName+"H","",100,0.,TMath::Pi());
    deltaPhi23H[iTagger] = new TH1F("deltaPhi23"+TaggerName+"H","",100,0.,TMath::Pi());

    ptFirstH[iTagger] = new TH1F("ptFirst"+TaggerName+"H","",100,0.,500.);
    ptSecondH[iTagger] = new TH1F("ptSecond"+TaggerName+"H","",100,0.,500.);
    ptThirdH[iTagger] = new TH1F("ptThird"+TaggerName+"H","",100,0.,500.);

    etaFirstH[iTagger] = new TH1F("etaFirst"+TaggerName+"H","",100,-2.5,2.5);
    etaSecondH[iTagger] = new TH1F("etaSecond"+TaggerName+"H","",100,-2.5,2.5);
    etaThirdH[iTagger] = new TH1F("etaThird"+TaggerName+"H","",100,-2.5,2.5);

    nJets20InTripleH[iTagger] = new TH1F("nJets20InTriple"+TaggerName+"H","",10,-0.5,9.5);
    nJets30InTripleH[iTagger] = new TH1F("nJets30InTriple"+TaggerName+"H","",10,-0.5,9.5);

    for (int isample=0; isample<6; ++isample) {

      dijetMassSamplesH[iTagger][isample] = new TH1F("dijetMass"+TaggerName+tripleBsamples[isample]+"H","",30,0.,600.);
      evBTagSamplesH[iTagger][isample] = new TH1F("evBTag"+TaggerName+tripleBsamples[isample]+"H","",6,-0.5,5.5);
      dijetMassEvBTagSamplesH[iTagger][isample] = new TH2F("dijetMassEvBTag"+TaggerName+tripleBsamples[isample]+"H","",20,0.,600.,6,-0.5,5.5);

      deltaR12SamplesH[iTagger][isample] = new TH1F("deltaR12"+TaggerName+tripleBsamples[isample]+"H","",100,0.,5.); 
      deltaR13SamplesH[iTagger][isample] = new TH1F("deltaR13"+TaggerName+tripleBsamples[isample]+"H","",100,0.,5.); 
      deltaR23SamplesH[iTagger][isample] = new TH1F("deltaR23"+TaggerName+tripleBsamples[isample]+"H","",100,0.,5.); 

      deltaEta12SamplesH[iTagger][isample] = new TH1F("deltaEta12"+TaggerName+tripleBsamples[isample]+"H","",100,0.,5.);
      deltaEta13SamplesH[iTagger][isample] = new TH1F("deltaEta13"+TaggerName+tripleBsamples[isample]+"H","",100,0.,5.);
      deltaEta23SamplesH[iTagger][isample] = new TH1F("deltaEta23"+TaggerName+tripleBsamples[isample]+"H","",100,0.,5.);

      deltaPhi12SamplesH[iTagger][isample] = new TH1F("deltaPhi12"+TaggerName+tripleBsamples[isample]+"H","",100,0.,TMath::Pi());
      deltaPhi13SamplesH[iTagger][isample] = new TH1F("deltaPhi13"+TaggerName+tripleBsamples[isample]+"H","",100,0.,TMath::Pi());
      deltaPhi23SamplesH[iTagger][isample] = new TH1F("deltaPhi23"+TaggerName+tripleBsamples[isample]+"H","",100,0.,TMath::Pi());

      ptFirstSamplesH[iTagger][isample] = new TH1F("ptFirst"+TaggerName+tripleBsamples[isample]+"H","",100,0.,500.);
      ptSecondSamplesH[iTagger][isample] = new TH1F("ptSecond"+TaggerName+tripleBsamples[isample]+"H","",100,0.,500.);
      ptThirdSamplesH[iTagger][isample] = new TH1F("ptThird"+TaggerName+tripleBsamples[isample]+"H","",100,0.,500.);

      etaFirstSamplesH[iTagger][isample] = new TH1F("etaFirst"+TaggerName+tripleBsamples[isample]+"H","",100,-2.5,2.5);
      etaSecondSamplesH[iTagger][isample] = new TH1F("etaSecond"+TaggerName+tripleBsamples[isample]+"H","",100,-2.5,2.5);
      etaThirdSamplesH[iTagger][isample] = new TH1F("etaThird"+TaggerName+tripleBsamples[isample]+"H","",100,-2.5,2.5);

      nJets20InTripleSamplesH[iTagger][isample] = new TH1F("nJets20InTriple"+TaggerName+tripleBsamples[isample]+"H","",10,-0.5,9.5);
      nJets30InTripleSamplesH[iTagger][isample] = new TH1F("nJets30InTriple"+TaggerName+tripleBsamples[isample]+"H","",10,-0.5,9.5);
    }
  }


  // ***************************************************
  // to-do : indexing of jet rank and btag working point 
  // ***************************************************

  TH1F * svMassFirstJetH = new TH1F("svMassFirstJetH","",11,-1.,10.);
  TH1F * svMassFirstJetTCHPTH = new TH1F("svMassFirstJetTCHPTH","",11,-1.,10.);
  TH1F * svMassFirstJetTCHP6H = new TH1F("svMassFirstJetTCHP6H","",11,-1.,10.);
  TH1F * svMassFirstJetCSVTH = new TH1F("svMassFirstJetCSVTH","",11,-1.,10.);
  TH1F * svMassFirstJetSSVHPTH = new TH1F("svMassFirstJetSSVHPTH","",11,-1.,10.);

  TH1F * svMassSecondJetH = new TH1F("svMassSecondJetH","",11,-1.,10.);
  TH1F * svMassSecondJetTCHPTH = new TH1F("svMassSecondJetTCHPTH","",11,-1.,10.);
  TH1F * svMassSecondJetTCHP6H = new TH1F("svMassSecondJetTCHP6H","",11,-1.,10.);
  TH1F * svMassSecondJetCSVTH = new TH1F("svMassSecondJetCSVTH","",11,-1.,10.);
  TH1F * svMassSecondJetSSVHPTH = new TH1F("svMassSecondJetSSVHPTH","",11,-1.,10.);

  TH1F * svMassThirdJetH = new TH1F("svMassThirdJetH","",11,-1.,10.);
  TH1F * svMassThirdJetTCHPTH = new TH1F("svMassThirdJetTCHPTH","",11,-1.,10.);
  TH1F * svMassThirdJetTCHP6H = new TH1F("svMassThirdJetTCHP6H","",11,-1.,10.);
  TH1F * svMassThirdJetCSVTH = new TH1F("svMassThirdJetCSVTH","",11,-1.,10.);
  TH1F * svMassThirdJetSSVHPTH = new TH1F("svMassThirdJetSSVHPTH","",11,-1.,10.);

  TH1F * ptFirstJetH = new TH1F("ptFirstJetH","",30,0.,300.);
  TH1F * ptSecondJetH = new TH1F("ptSecondJetH","",30,0.,300.);
  TH1F * ptThirdJetH = new TH1F("ptThirdJetH","",30,0.,300.);

  TH1F * svMassJetFlavorH[3][3][4];

  TString jetStr[3] = {"FirstJet","SecondJet","ThirdJet"};

  for (int iJet=0; iJet<3; ++iJet) {
    for (int iFlv=0; iFlv<3; ++iFlv) {
      for (int iTagger=0; iTagger<4; ++iTagger) {
	TString histoName = "svMass"+jetStr[iJet] + flavors[iFlv] + TString(sovl[iTagger])+"H";
	svMassJetFlavorH[iJet][iFlv][iTagger] = new TH1F(histoName,"",11,-1.,10.);
      }
    }
  }


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

  TH1F * bTagTemplateBBimpurityH[4][3][3];
  TH2F * massBTagTemplateBBimpurityH[4][3][3];

  TH1F * bTagOnlineCorrectedTemplateH[4][3][3];
  TH2F * massBTagOnlineCorrectedTemplateH[4][3][3];

  TH1F * bTagOnlineCorrectedTemplateBBimpurityH[4][3][3];
  TH2F * massBTagOnlineCorrectedTemplateBBimpurityH[4][3][3];

  for (int isovl=0; isovl<4; ++isovl) {
    for (int iFlav=0; iFlav<3; ++iFlav) {
      for (int iCat=0; iCat<3; ++iCat) {
        bTagTemplateH[isovl][iCat][iFlav] = new TH1F("bTagTemplate_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",6,-0.5,5.5);
        massBTagTemplateH[isovl][iCat][iFlav] = new TH2F("MassBTagTemplate_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",20,0.,600.,6,-0.5,5.5);

        bTagTemplateBBimpurityH[isovl][iCat][iFlav] = new TH1F("bTagTemplateBBimpurity_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",6,-0.5,5.5);
        massBTagTemplateBBimpurityH[isovl][iCat][iFlav] = new TH2F("MassBTagTemplateBBimpurity_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",20,0.,600.,6,-0.5,5.5);

	bTagOnlineCorrectedTemplateH[isovl][iCat][iFlav] = new TH1F("bTagTemplate_OnlCorr_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",6,-0.5,5.5); 
	massBTagOnlineCorrectedTemplateH[isovl][iCat][iFlav] = new TH2F("MassBTagTemplate_OnlCorr_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",20,0.,600.,6,-0.5,5.5);

	bTagOnlineCorrectedTemplateBBimpurityH[isovl][iCat][iFlav] = new TH1F("bTagTemplateBBimpurity_OnlCorr_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",6,-0.5,5.5); 
	massBTagOnlineCorrectedTemplateBBimpurityH[isovl][iCat][iFlav] = new TH2F("MassBTagTemplateBBimpurity_OnlCorr_"+TString(sovl[isovl])+"_"+templateNames[iCat][iFlav]+"H","",20,0.,600.,6,-0.5,5.5);


      }
    }
  }

  //   tPatTemplateH[tagger][categ]
  TH1F * tPatternH[4][3];
  TString CategNames[3] = {"Xbb","bXb","bbX"};
  TString tPatNames[8] = {"xxx","Txx","xTx","xxT","xTT","TxT","TTx","TTT"};
  TString evBTag12names[6] = {"Sum0","Sum1","Sum2",
			      "xxx0","xxx1","xxx2"};
  TString evBTag3names[6]  = {"00","10","20",
			      "01","11","21"};
  
  TH1F * dijetMassDoubleBBH[4][3][6];
  TH1F * dijetMassDoubleBBXH[4][3][6];

  TH1F * dijetMassDoubleTotBBH[4][3];
  TH1F * dijetMassDoubleTotBBXH[4][3];

  for (int iTagger=0; iTagger<4; ++iTagger) {
    for (int iCat=0; iCat<3; ++iCat) {
      tPatternH[iTagger][iCat] = new TH1F("tBattern_"+TString(sovl[iTagger])+"_"+CategNames[iCat]+"_H",CategNames[iCat]+" "+TString(sovl[iTagger]),8,-0.5,7.5);
      for (int iB=0; iB<8; ++iB) 
	tPatternH[iTagger][iCat]->GetXaxis()->SetBinLabel(iB+1,tPatNames[iB]);
      for (int iConfig=0; iConfig<6; ++iConfig) {
	TString histName = "diJetMass_" + TString(sovl[iTagger])+"_" + CategNames[iCat] + "_" + evBTag12names[iConfig] + "_BB";
	if (iCat<2)
	  histName = "diJetMass_" + TString(sovl[iTagger])+"_" + CategNames[iCat] + "_" + evBTag3names[iConfig] + "_BB";

	dijetMassDoubleBBH[iTagger][iCat][iConfig] = new TH1F(histName+"H","",60,0.,600.);
	dijetMassDoubleBBXH[iTagger][iCat][iConfig] = new TH1F(histName+"XH","",60,0.,600.);
      }
      TString histName = "diJetMass_" + TString(sovl[iTagger])+"_" + CategNames[iCat] + "_Tot_BB";
      dijetMassDoubleTotBBH[iTagger][iCat] = new TH1F(histName+"H","",60,0.,600.);
      dijetMassDoubleTotBBXH[iTagger][iCat] = new TH1F(histName+"XH","",60,0.,600.);
    }
  }

  // Opening filelist
  // reading out number of files
  std::ifstream inputFile(argv[1]);
  inputFile >> numberOfFiles;

  // *************************************
  // *** L O O P   O V E R   F I L E S ***
  // *************************************
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

    HBBTo4B_MC * hbbTo4b = new HBBTo4B_MC(tree);
    hbbTo4b->Init(tree);

    Trigger * trigger = new Trigger(file,triggers,true);
    //    trigger->PrintOut();

    // ***************************************
    // *** L O O P   O V E R   E V E N T S ***
    // ***************************************
    for (int iE=0; iE<numberOfEntries; ++iE) {

      nEvents++;
      if (nEvents%100000==0)
    	std::cout << "Processed " << nEvents << " events " << std::endl;

      hbbTo4b->GetEntry(iE);

      int TriggerAccept = hbbTo4b->trgAccept;
      bool accept = trigger->acceptAnd(TriggerAccept);

      if (!applyTrigger)
	accept = true;

      if (!accept) continue;

      int nJets = 0;
      int numberOfJets = hbbTo4b->NumberOfJets;

      numberOfJets = TMath::Min(numberOfJets,int(28));

      std::vector<int> leadingJets;
      leadingJets.clear();

      int nJets20 = 0;
      int nJets30 = 0;

      
      // bH->bbb analysis
      // *************************************************
      // L O O P    O V E R    J E T S
      // *************************************************

      for (int iJ=0; iJ<numberOfJets; ++iJ) {

	// bH->bbb analysis
	// Here we apply JES shifts
	// to jet's Pt 
	if (applyJESShift) {
	  float jetEtaX = TMath::Abs(hbbTo4b->JetEta[iJ]);
	  float jetPtX  = hbbTo4b->JetPt[iJ];
	    
	  if (jesShift>0) {
	    float jesUncertUp = hbbSyst->getJESuncertaintyUp(jetPtX,jetEtaX);
	    hbbTo4b->JetPt[iJ] = (1+jesShift*jesUncertUp)*hbbTo4b->JetPt[iJ];
	  }
	  else {
	    float jesUncertDown = hbbSyst->getJESuncertaintyUp(jetPtX,jetEtaX);
	    hbbTo4b->JetPt[iJ] = (1+jesShift*jesUncertDown)*hbbTo4b->JetPt[iJ];
	  }
	}

	bool jetIDLoose = hbbTo4b->JetIDLoose[iJ];
	if (!jetIDLoose) continue;

	float absJetEta = TMath::Abs(hbbTo4b->JetEta[iJ]);
	if (absJetEta>jetEtaMax) continue;

	float jetPt = hbbTo4b->JetPt[iJ];  

	if (jetPt>20)
	  nJets20++;

	if (jetPt>30)
	  nJets30++;

      }

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

      bool refTrigger = trigger->acceptAnd(TriggerAccept);


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

      diJetMassBeforeTriggerH->Fill(diJetMass,eventWeight);

      if (refTrigger)
	diJetMassAfterTriggerH->Fill(diJetMass,eventWeight);

      int jetIndex[3];
      jetIndex[0] = j1;
      jetIndex[1] = j2;
      jetIndex[2] = j3;


      float jetPtArray[3];
      float jetEtaArray[3];
      float jetEtaSignedArray[3];
      float jetPhiArray[3];
      int jetFlavor[3];
      float jetOfflineBTagEff[3][4];
      float jetRelOnlineBTagEff[3][4];
      for (int index=0; index<3; ++index) {
	jetPtArray[index]  = hbbTo4b->JetPt[jetIndex[index]];
	jetEtaArray[index] = TMath::Abs(hbbTo4b->JetEta[jetIndex[index]]);
	jetEtaSignedArray[index] = hbbTo4b->JetEta[jetIndex[index]];
	jetPhiArray[index] = hbbTo4b->JetPhi[jetIndex[index]];
	//	float jesUncertUp = hbbSyst->getJESuncertaintyUp(jetPtArray[index],jetEtaArray[index]);
	//	float jesUncertDown = hbbSyst->getJESuncertaintyUp(jetPtArray[index],jetEtaArray[index]);
	//	std::cout << index << " : pt = " << jetPtArray[index] << " ; eta = " << jetEtaArray[index]
	//		  << "  : up = " << jesUncertUp << "  : down = " << jesUncertDown << std::endl;
	int partonFlavorJet = hbbTo4b->PartonFlavorJet[jetIndex[index]];
	int jetFlavorAbs = TMath::Abs(partonFlavorJet);
	int HflContentInJet = TMath::Abs(hbbTo4b->HflContentJet[jetIndex[index]]);
	jetFlavor[index] = 0;
	if (jetFlavorAbs==4)
	  jetFlavor[index] = 1;
	else if (jetFlavorAbs==5)
	  jetFlavor[index] = 2;
	else {
	  if (HflContentInJet==4)
	    jetFlavor[index] = 1;
	  if (HflContentInJet==5)
	    jetFlavor[index] = 2;
	}
	for (int itagger=0; itagger<4;++itagger) {
	  jetOfflineBTagEff[index][itagger] = bTagEffOffline->eff(jetFlavor[index],sovl[itagger],jetPtArray[index],jetEtaArray[index]);
	  jetRelOnlineBTagEff[index][itagger] = bTagRelEffOnline->eff(jetFlavor[index],sovl[itagger],jetPtArray[index],jetEtaArray[index]);
	}
      }

      // *******************************************************
      // to-do : indexing of taggers, flavors ; iCat (if needed)
      // pointer to taggers --------------->
      // *******************************************************


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
      //
      if (printOut) {
	std::cout << "Jet 1 : dTCHP = " << btag1
		  << "  dCSV = " << dCsv1
		  << "  dSSVHP = " << dSsvhp1 
		  << "  hltB = " << bTagHlt1 
		  << "  Flavor = " << jetFlavor[0] << std::endl;
	
	std::cout << "Jet 2 : dTCHP = " << btag2
		  << "  dCSV = " << dCsv2
		  << "  dSSVHP = " << dSsvhp2 
		  << "  hltB = " << bTagHlt2 
		  << "  Flavor = " << jetFlavor[1] << std::endl;
      
	std::cout << "Jet 3 : dTCHP = " << btag3
		  << "  dCSV = " << dCsv3
		  << "  dSSVHP = " << dSsvhp3 
		  << "  hltB = " << bTagHlt3 
		  << "  Flavor = " << jetFlavor[2] << std::endl;
	std::cout << std::endl;
      }
      //      float eta1 = TMath::Abs(hbbTo4b->JetEta[j1]);
      //      float eta2 = TMath::Abs(hbbTo4b->JetEta[j2]);
      //      float eta3 = TMath::Abs(hbbTo4b->JetEta[j3]);

      float svMassFull1 = svMass1;
      if (svMassFull1<0) svMassFull1 = -0.5;

      float svMassFull2 = svMass2;
      if (svMassFull2<0) svMassFull2 = -0.5;

      float svMassFull3 = svMass3;
      if (svMassFull3<0) svMassFull3 = -0.5;
      

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
	    tPatternH[iTagger][iCat]->Fill(float(tpat),eventWeight);
	  }
	}
      }

      // *************************************************
      // select events where at least two out of 
      // three leading jets are matched to BTag Hlt object
      // *************************************************

      if (thisTpat<0 && applyHltBTagMatch) 
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

      float weight12 = jetOfflineBTagEff[0][1]*jetRelOnlineBTagEff[0][1]*jetOfflineBTagEff[1][1]*jetRelOnlineBTagEff[1][1];
      float weight13 = jetOfflineBTagEff[0][1]*jetRelOnlineBTagEff[0][1]*jetOfflineBTagEff[2][1]*jetRelOnlineBTagEff[2][1];
      float weight23 = jetOfflineBTagEff[1][1]*jetRelOnlineBTagEff[1][1]*jetOfflineBTagEff[2][1]*jetRelOnlineBTagEff[2][1];

      float svMassFullJets[3];
      svMassFullJets[0] = svMassFull1;
      svMassFullJets[1] = svMassFull2;
      svMassFullJets[2] = svMassFull3;

      bool jetTags[3][4];

      jetTags[0][0] = b1TCHPT;
      jetTags[0][1] = b1TCHP6;
      jetTags[0][2] = b1CSVT;
      jetTags[0][3] = b1SSVHPT;

      jetTags[1][0] = b2TCHPT;
      jetTags[1][1] = b2TCHP6;
      jetTags[1][2] = b2CSVT;
      jetTags[1][3] = b2SSVHPT;
     
      jetTags[2][0] = b3TCHPT;
      jetTags[2][1] = b3TCHP6;
      jetTags[2][2] = b3CSVT;
      jetTags[2][3] = b3SSVHPT;

      float jetWeights[3];
      jetWeights[0] = weight23;
      jetWeights[1] = weight13;
      jetWeights[2] = weight12;

      float weight = eventWeight;

      for (int iJet=0; iJet<3; ++iJet) {
	for (int itagger=0; itagger<4; ++itagger) {
	  if (jetTags[iJet][itagger])
	    svMassJetFlavorH[iJet][jetFlavor[iJet]][itagger]->Fill(svMassFullJets[iJet],weight*jetWeights[iJet]);
	}
      }

      //      std::cout << "weight = " << weight << "  weight23 = " << weight23 << std::endl;
      
      ptFirstJetH->Fill(pt1,weight*weight23);
      svMassFirstJetH->Fill(svMassFull1,weight*weight23);
      if (b1TCHPT)
	svMassFirstJetTCHPTH->Fill(svMassFull1,weight*weight23);
      if (b1TCHP6)
	svMassFirstJetTCHP6H->Fill(svMassFull1,weight*weight23);
      if (b1CSVT)
	svMassFirstJetCSVTH->Fill(svMassFull1,weight*weight23);
      if (b1SSVHPT)
	svMassFirstJetSSVHPTH->Fill(svMassFull1,weight*weight23);
      
      ptSecondJetH->Fill(pt2,weight*weight13);
      svMassSecondJetH->Fill(svMassFull2,weight*weight13);
      if (b2TCHPT)
	svMassSecondJetTCHPTH->Fill(svMassFull2,weight*weight13);
      if (b2TCHP6)
	svMassSecondJetTCHP6H->Fill(svMassFull2,weight*weight13);
      if (b2CSVT)
	svMassSecondJetCSVTH->Fill(svMassFull2,weight*weight13);
      if (b2SSVHPT)
	svMassSecondJetSSVHPTH->Fill(svMassFull2,weight*weight13);

      ptThirdJetH->Fill(pt3,weight*weight12);
      svMassThirdJetH->Fill(svMassFull3,weight*weight12);
      if (b3TCHPT)
	svMassThirdJetTCHPTH->Fill(svMassFull3,weight*weight12);
      if (b3TCHP6)
	svMassThirdJetTCHP6H->Fill(svMassFull3,weight*weight12);
      if (b3CSVT)
	svMassThirdJetCSVTH->Fill(svMassFull3,weight*weight12);
      if (b3SSVHPT)
	svMassThirdJetSSVHPTH->Fill(svMassFull3,weight*weight12);

      // *****************************
      // indexing of b-taggers ----->
      // *****************************

      float trigWeight[4] = {1,1,1,1};

      if (applyBTagHltWeightOnTriple) {
	for (int iTagger=0; iTagger<4; ++iTagger) {
	  trigWeight[iTagger] = 
	    jetRelOnlineBTagEff[1][iTagger]*jetRelOnlineBTagEff[2][iTagger] +
	    jetRelOnlineBTagEff[0][iTagger]*jetRelOnlineBTagEff[1][iTagger]*(1-jetRelOnlineBTagEff[2][iTagger]) +
	    jetRelOnlineBTagEff[0][iTagger]*(1-jetRelOnlineBTagEff[1][iTagger])*jetRelOnlineBTagEff[2][iTagger];
	}      
      }

      if (applyWeightOnTriple) {
	evBtagTCHPT  = true;
	evBtagTCHP6  = true;
	evBtagCSVT   = true;
	evBtagSSVHPT = true;
      }


      bool tripleTag[4];
      tripleTag[0] = evBtagTCHPT;
      tripleTag[1] = evBtagTCHP6;
      tripleTag[2] = evBtagCSVT;
      tripleTag[3] = evBtagSSVHPT;

      float offlTripleTagEff[4] = {1,1,1,1}; 
      float relOnlTripleTagEff[4] = {1,1,1,1};

      if (applyWeightOnTriple) {
	for (int itagger=0; itagger<4; ++itagger) {
	  for (int index=0; index<3; ++index) {
	    offlTripleTagEff[itagger] *= jetOfflineBTagEff[index][itagger]; 
	  }
	}
      }

      for (int itagger=0; itagger<4; ++itagger) {
	for (int index=0; index<3; ++index) {
	  relOnlTripleTagEff[itagger] *= jetRelOnlineBTagEff[index][itagger]; 
	}
      }

      float offlDoubleTagEff[4][3];
      float relOnlDoubleTagEff[4][3];

      for (int itagger=0; itagger<4; ++itagger) {
	for (int icat=0; icat<3; ++icat) {
	  if (applyWeightOnDouble) {
	    bTagger[itagger][icat] = true;
	    offlDoubleTagEff[itagger][icat] = offlTripleTagEff[itagger]/jetOfflineBTagEff[icat][itagger];
	  }
	  else {
	    offlDoubleTagEff[itagger][icat] = 1;
	  }    
	  if (applyBTagHltWeightOnDouble) 
	    relOnlDoubleTagEff[itagger][icat] = relOnlTripleTagEff[itagger]/jetRelOnlineBTagEff[icat][itagger];
	  else 
	    relOnlDoubleTagEff[itagger][icat] = 1;
	}
      }

      int isample = 5;

      if ((jetFlavor[0]==2)&&(jetFlavor[1]==2)&&(jetFlavor[2]==2))
	isample = 0;
      if ((jetFlavor[0]==2)&&(jetFlavor[1]==2)&&(jetFlavor[2]==1))
	isample = 1;
      if ((jetFlavor[0]==2)&&(jetFlavor[1]==2)&&(jetFlavor[2]==0))
	isample = 2;
      if ((jetFlavor[0]==2)&&(jetFlavor[1]==1)&&(jetFlavor[2]==2))
	isample = 3;
      if ((jetFlavor[0]==1)&&(jetFlavor[1]==2)&&(jetFlavor[2]==2))
        isample = 3;
      if ((jetFlavor[0]==2)&&(jetFlavor[1]==0)&&(jetFlavor[2]==2))
        isample = 4;
      if ((jetFlavor[0]==0)&&(jetFlavor[1]==2)&&(jetFlavor[2]==2))
        isample = 4;

      float dR12 = deltaR(jetEtaSignedArray[0],jetPhiArray[0],
			  jetEtaSignedArray[1],jetPhiArray[1]);
      float dR13 = deltaR(jetEtaSignedArray[0],jetPhiArray[0],
                          jetEtaSignedArray[2],jetPhiArray[2]);
      float dR23 = deltaR(jetEtaSignedArray[1],jetPhiArray[1],
                          jetEtaSignedArray[2],jetPhiArray[2]);

      float dEta12 = TMath::Abs(jetEtaSignedArray[0]-jetEtaSignedArray[1]);
      float dEta13 = TMath::Abs(jetEtaSignedArray[0]-jetEtaSignedArray[2]);
      float dEta23 = TMath::Abs(jetEtaSignedArray[1]-jetEtaSignedArray[2]);

      float dPhi12 = TMath::Abs(dPhiFromPhi(jetPhiArray[0],jetPhiArray[1]));
      float dPhi13 = TMath::Abs(dPhiFromPhi(jetPhiArray[0],jetPhiArray[2]));
      float dPhi23 = TMath::Abs(dPhiFromPhi(jetPhiArray[1],jetPhiArray[2]));

      float eta1 = jetEtaSignedArray[0]; 
      float eta2 = jetEtaSignedArray[1]; 
      float eta3 = jetEtaSignedArray[2]; 


      // bH->bbb analysis
      // ****************************
      // Triple tag selection ======>
      // ****************************
      // run over 4 taggers 
      // 0 : TCHPT
      // 1 : TCHP6
      // 2 : CSVT
      // 3 : SSVHPT
      //

      for (int itagger=0; itagger<4; ++itagger) {

	if (tripleTag[itagger]) {

	  // bH->bbb analysis
	  // *********************
	  // inclusive samples
	  // *********************
	  // When running on MC 
	  // additional weight, corresponding to product of BTag scale factors
	  // for each jet, should be applied 
	  float SF[3]; // Scale Factors for the flavors of three leading jets
	  float btagSFweight = 1;
	  if (applyBTagSF) {
	    for (int iJ=0; iJ<3; ++iJ) {

	      SF[iJ] = hbbSyst->getSFtag(jetPtArray[iJ],jetEtaArray[iJ],jetFlavor[iJ],itagger);

	      std::cout << "flavor : " << jetFlavor[iJ] << "   ;  Scale factor = " << SF[iJ] << std::endl;

	      if (applyBTagShift) { // systematic shift in the btag scale factor 
		if (jetFlavor[iJ]==0)
		  SF[iJ] += btagSFlighShift*hbbSyst->getSFuncertainty(jetPtArray[iJ],jetEtaArray[iJ],jetFlavor[iJ],itagger,btagLFUp); 
		else 
		  SF[iJ] += btagSFbcShift*hbbSyst->getSFuncertainty(jetPtArray[iJ],jetEtaArray[iJ],jetFlavor[iJ],itagger,btagLFUp);
	      }

	      std::cout << " Syst shift : " << hbbSyst->getSFuncertainty(jetPtArray[iJ],jetEtaArray[iJ],jetFlavor[iJ],itagger,btagLFUp)
			<< std::endl;
		
	      // additional weight for triple-btag in Monte Carlo
	      btagSFweight *= SF[iJ];
	      
	    }
	  }

	  dijetMassH[itagger]->Fill(diJetMass, weight * btagSFweight * offlTripleTagEff[itagger] * trigWeight[itagger] );
	  evBTagH[itagger]->Fill(float(nEvBTag), weight * btagSFweight * offlTripleTagEff[itagger] * trigWeight[itagger] );
	  dijetMassEvBTagH[itagger]->Fill(diJetMass,float(nEvBTag), weight * btagSFweight * offlTripleTagEff[itagger] * trigWeight[itagger] );

	  deltaR12H[itagger]->Fill(dR12,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaR13H[itagger]->Fill(dR13,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaR23H[itagger]->Fill(dR23,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  deltaEta12H[itagger]->Fill(dEta12,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaEta13H[itagger]->Fill(dEta13,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaEta23H[itagger]->Fill(dEta23,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  deltaPhi12H[itagger]->Fill(dPhi12,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaPhi13H[itagger]->Fill(dPhi13,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaPhi23H[itagger]->Fill(dPhi23,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  ptFirstH[itagger]->Fill(pt1,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  ptSecondH[itagger]->Fill(pt2,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  ptThirdH[itagger]->Fill(pt3,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  etaFirstH[itagger]->Fill(eta1,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  etaSecondH[itagger]->Fill(eta2,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  etaThirdH[itagger]->Fill(eta3,weight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  nJets20InTripleH[itagger]->Fill(float(nJets20),weight*offlTripleTagEff[itagger]* trigWeight[itagger]);
	  nJets30InTripleH[itagger]->Fill(float(nJets30),weight*offlTripleTagEff[itagger]* trigWeight[itagger]);

	  // ****************************************
	  // break-down into different flavor content
	  // ****************************************

	  dijetMassSamplesH[itagger][isample]->Fill(diJetMass, weight * btagSFweight * offlTripleTagEff[itagger] * trigWeight[itagger] );
	  evBTagSamplesH[itagger][isample]->Fill(float(nEvBTag), weight * btagSFweight * offlTripleTagEff[itagger] * trigWeight[itagger] );
	  dijetMassEvBTagSamplesH[itagger][isample]->Fill(diJetMass,float(nEvBTag), weight * btagSFweight * offlTripleTagEff[itagger] * trigWeight[itagger] );

	  deltaR12SamplesH[itagger][isample]->Fill(dR12,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaR13SamplesH[itagger][isample]->Fill(dR13,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaR23SamplesH[itagger][isample]->Fill(dR23,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  deltaEta12SamplesH[itagger][isample]->Fill(dEta12,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaEta13SamplesH[itagger][isample]->Fill(dEta13,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaEta23SamplesH[itagger][isample]->Fill(dEta23,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  deltaPhi12SamplesH[itagger][isample]->Fill(dPhi12,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaPhi13SamplesH[itagger][isample]->Fill(dPhi13,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  deltaPhi23SamplesH[itagger][isample]->Fill(dPhi23,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  ptFirstSamplesH[itagger][isample]->Fill(pt1,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  ptSecondSamplesH[itagger][isample]->Fill(pt2,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  ptThirdSamplesH[itagger][isample]->Fill(pt3,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  etaFirstSamplesH[itagger][isample]->Fill(eta1,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  etaSecondSamplesH[itagger][isample]->Fill(eta2,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );
	  etaThirdSamplesH[itagger][isample]->Fill(eta3,weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger] );

	  nJets20InTripleSamplesH[itagger][isample]->Fill(float(nJets20),weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger]);
	  nJets30InTripleSamplesH[itagger][isample]->Fill(float(nJets30),weight*btagSFweight*offlTripleTagEff[itagger]* trigWeight[itagger]);

	}

      }

      // **************************
      // Creating templates 
      // from double tag samples
      // **************************

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

      float svMass[3];
      svMass[0] = svMass1;
      svMass[1] = svMass2;
      svMass[2] = svMass3;

      // Templates

      bool doubleBflv[3] = {false, false, false};
      doubleBflv[0] = (jetFlavor[1]==2) && (jetFlavor[2]==2);
      doubleBflv[1] = (jetFlavor[0]==2) && (jetFlavor[2]==2);
      doubleBflv[2] = (jetFlavor[0]==2) && (jetFlavor[1]==2);
      

      for (int iTagger=0; iTagger<4; ++iTagger) {
	std::vector<float*> probsOfflineTpat;
	std::vector<float*> probsOnlineTpat;
	probsOfflineTpat.clear();
	probsOnlineTpat.clear();
	for (int itpat=0; itpat<3; ++itpat) {
	  float probsOnl[3];
	  float probsOffl[3];
	  svOnline->getSVbins(2,iTagger,jetPtArray[itpat],jetEtaArray[itpat],probsOnl);
	  probsOnlineTpat.push_back(probsOnl);
	  sv->getSVbins(2,iTagger,jetPtArray[itpat],jetEtaArray[itpat],probsOffl);
	  probsOfflineTpat.push_back(probsOffl);
	}
	for (int iCat=0; iCat<3; ++iCat) {
	  bool isDoubleHltMatch = (thisTpat==iCat) || (thisTpat == 3) ;
	  if (applyBTagHltWeightOnDouble) 
	    isDoubleHltMatch = true;
	  if (isDoubleHltMatch) { // online btag pattern must match offline
	    if (bTagger[iTagger][iCat]) {
	      int nJ1 = 1;
	      int nJ2 = 2;
	      if (iCat==0) {
		nJ1 = 1;
		nJ2 = 2;
	      }
	      if (iCat==1) {
		nJ1 = 0;
		nJ2 = 2;
	      }
	      if (iCat==2) {
		nJ1 = 0;
		nJ2 = 1;
	      }

	      float purity = bbpur->getBBPurity(iTagger,iCat,diJetMass,svMass[nJ1],svMass[nJ2]);

	      //	      if ( svMass[nJ1]>0.5 && svMass[nJ2]>0.5 )
	      //		std::cout << sovl[iTagger] << " " << CategNames[iCat] << " : M12 = " 
	      //			  << diJetMass << " svM1 = " << svMass[nJ1] << "  svM2 = " << svMass[nJ2] 
	      //			  << "  bb(purity) = " << purity << std::endl;

	      int nConfig = 0;
	      if (iCat==2) { // bbX
		int T12 = nB1 + nB2;
		if (T12<2) 
		  nConfig = 0;
		else if (T12<3)
		  nConfig = 1;
		else 
		  nConfig = 2;
	      }
	      else { // bXb or Xbb
		int T12 = nB1;
		if (iCat==1)
		  T12 = nB2;
		int T3 = 0;
		if (nB3==2)
		  T3 = 3;
		nConfig = T3 + T12; 
	      }
	      float weightTag = offlDoubleTagEff[iTagger][iCat] * relOnlDoubleTagEff[iTagger][iCat];
	      if (doubleBflv[iCat]) {
		dijetMassDoubleTotBBH[iTagger][iCat]->Fill(diJetMass,weight*weightTag);
		dijetMassDoubleBBH[iTagger][iCat][nConfig]->Fill(diJetMass,weight*weightTag);
	      }
	      else {
		dijetMassDoubleTotBBXH[iTagger][iCat]->Fill(diJetMass,weight*weightTag);
		dijetMassDoubleBBXH[iTagger][iCat][nConfig]->Fill(diJetMass,weight*weightTag);
	      }
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
		  float templateWeight = weight*eff[iFlv]*offlDoubleTagEff[iTagger][iCat]*relOnlDoubleTagEff[iTagger][iCat];
		  bTagTemplateH[iTagger][iCat][iFlv]->Fill(evtBTag,templateWeight * probs[iTag]);
		  massBTagTemplateH[iTagger][iCat][iFlv]->Fill(diJetMass,evtBTag,templateWeight * probs[iTag]);
		  if (!doubleBflv[iCat]) {
		    bTagTemplateBBimpurityH[iTagger][iCat][iFlv]->Fill(evtBTag,templateWeight * probs[iTag]);
		    massBTagTemplateBBimpurityH[iTagger][iCat][iFlv]->Fill(diJetMass,evtBTag,templateWeight * probs[iTag]);
		  }
		  float addWeight = 0;
		  for (int itpat=0; itpat<3; ++itpat) {
		    float wtpat = 1;
		    if ((itpat != iCat) && (thisTpat != 3)) { // weight correction if offline 
		                                              // & online btag pattern are different (except TTT)
		      float effOnlineTpat = bTagRelEffOnline->eff(2,sovl[iTagger],jetPtArray[itpat],jetEtaArray[itpat]);
		      if ( effOnlineTpat>0 )
			wtpat = relEffOnl[iFlv] / effOnlineTpat;
		      
		      if ( probsOnlineTpat.at(itpat)[iTag]>0 )
			wtpat *= probsOnline[iTag]*probsOfflineTpat.at(itpat)[iTag]/probsOnlineTpat.at(itpat)[iTag];
		      
		    }
		    else 
		      wtpat *= probs[iTag];
		    addWeight += wtpat;
		  }
		  bTagOnlineCorrectedTemplateH[iTagger][iCat][iFlv]->Fill(evtBTag,templateWeight * addWeight);
		  massBTagOnlineCorrectedTemplateH[iTagger][iCat][iFlv]->Fill(diJetMass,evtBTag,templateWeight * addWeight);
		  if (!doubleBflv[iCat]) {
		    bTagOnlineCorrectedTemplateBBimpurityH[iTagger][iCat][iFlv]->Fill(evtBTag,templateWeight * addWeight);
		    massBTagOnlineCorrectedTemplateBBimpurityH[iTagger][iCat][iFlv]->Fill(diJetMass,evtBTag,templateWeight * addWeight);
		  }
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
