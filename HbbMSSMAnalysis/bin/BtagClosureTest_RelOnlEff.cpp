#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TH1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <TFileCollection.h>
#include <TGraph2D.h>
#include <TLatex.h>
#include <TMarker.h>
#include <math.h>
#include <TMath.h>
#include <TFileInfo.h>

//#include "Analysis/Utilities/interface/HBBTo4B.h"
//#include "Analysis/Utilities/src/HBBTo4B.cc"
#include "Analysis/Utilities/interface/svMass.h"
#include "Analysis/Utilities/interface/bTagEff.h"
#include "Analysis/Utilities/interface/TrigHistArray.h"
#include "Analysis/Utilities/interface/TrigHistArray2D.h"
#include "Analysis/Utilities/interface/getHbbMetadata.h"
#include "Analysis/Utilities/interface/checkHLTBTagMatch.h"
#include "Analysis/Utilities/interface/jetFlavorCode.h"
#include "Analysis/HbbMSSMAnalysis/test/NegativeBTag.C"

#include "Analysis/Utilities/interface/HbbNtuple.h"
#include "Analysis/Utilities/interface/HbbSelect.h"
#include "Analysis/Utilities/interface/TriggerRunSelection.h"

#include "Analysis/Utilities/interface/HbbSystControl.h"
#include "Analysis/Utilities/interface/HbbTrigger.h"
#include "Analysis/Utilities/interface/HbbTrigWeight.h"
//#include "Analysis/Utilities/interface/testJetIsolation.h"
#include "Analysis/Utilities/interface/deltaR.h"



bool RedoHLTBtagMatch(int iJet)
{
  bool matched = false;
  for(int iL3jet = 0; iL3jet < l3NumberOfBJets; iL3jet++) {
        
            
    const float delta_R_L3_RecoJet = CalculateDeltaR(l3PhiBJet[iL3jet], phiJet[iJet], l3EtaBJet[iL3jet], etaJet[iJet]);

 

    double mindeltaR_L3_L2 = 1000.0;


    for (int iJetL2=0; iJetL2 < l2NumberOfJets; ++iJetL2) {
            
      const double deltaR = CalculateDeltaR(l2PhiJet[iJetL2] ,l3PhiBJet[iL3jet], l2EtaJet[iJetL2] , l3EtaBJet[iL3jet]);
          
      if(deltaR < mindeltaR_L3_L2) {
        mindeltaR_L3_L2 = deltaR;
      }
    }
            
    if(delta_R_L3_RecoJet < 0.5 && mindeltaR_L3_L2 < 0.5) {
              
      matched = true;
    }
  }
          
  return matched;
}



//void TripleBtagAnalysis_SF() {
int main(int narg,char** varg) {


  const bool doTrigSelect = false;//true;  // perform trigger selection according to scenario
  const std::string Scenario("MediumMass2012");  //not used// Higgs mass scenario // available: LowMass2012, MediumMass2012, HighMass2012
  const std::string onlinescenario("normal");
  // open an ntuple file

  const bool _doMC = true;

  const  bool doLumiWeighting = true;  // weight MC events according to intgrated lumi

  // chain mode
  TChain f("hbbanalysis/HBBTo4B");
  TFileCollection theFileColl("myNtFileCollection","","theMergeList.txt");
  if (theFileColl.GetNFiles() <=0) {
    std::cout << "Problem opening list of input files " << std::endl;
    return 0;
  }
  TIter next( (TCollection*) theFileColl.GetList() );
  TFileInfo* tfi = (TFileInfo*) next();
  TString firstFile( tfi->GetFirstUrl()->GetFile() );
  if (firstFile.Contains("Run2011") || firstFile.Contains("Run2012")) {
    std::cout << "Assume this is RealData. works only for MC. exit" << std::endl;
    return 1;
  } else {
    std::cout << "Assume this is MonteCarlo" << std::endl;
   
  }  

  // create the trigger filter list. Only histograms for these triggers will be created
  std::vector<std::string> genericTriggerList;
  std::vector<std::string> HLTBtagTriggerObjectList;
  std::vector<std::string> L25BtagTriggerObjectList;
  std::vector<std::string> triggerFilterList;
  std::vector<std::string> HLTBtagTriggerObjectFilterList;
  std::vector<std::string> L25BtagTriggerObjectFilterList;
  
  // extract generic trigger list and number of input events from ntuple files
  float lumiScaleFac = 1;
  getHbbMetadata(theFileColl,genericTriggerList,HLTBtagTriggerObjectList,1000.,lumiScaleFac,_doMC);
  getHbbL25Metadata(theFileColl,L25BtagTriggerObjectList);

  std::string theTrigNameBtag("");
  if(onlinescenario == "normal") {
    theTrigNameBtag.assign("HLT_DiJet40Eta2p6_BTagIP3DFastPV_v");
  } else if(onlinescenario == "loose") {
    theTrigNameBtag.assign("HLT_DiJet80Eta2p6_BTagIP3DFastPVLoose_v");
  }

  const std::string theTrigNameRef("HLT_L1DoubleJet36Central_v");
  
  triggerFilterList.push_back(theTrigNameBtag);

  if(onlinescenario == "normal") {
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BL1FastJetFastPV");
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhi1BL1FastJetFastPV");
    
  } else if(onlinescenario == "loose") {
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BLooseL1FastJetFastPV");
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhi1B1stTrackL1FastJetFastPV");
  }
  
  f.AddFileInfoList((TCollection*) theFileColl.GetList());
  TTree* hbbtree = &f;

  // create trigger run selection object
  TriggerRunSelection trigSelector( Scenario, &genericTriggerList );
  checkHLTBtagMatch HLTBtagMatchObj(
                                    triggerFilterList, 
                                    HLTBtagTriggerObjectFilterList, 
                                    HLTBtagTriggerObjectList,
                                    doTrigSelect,
                                    &trigSelector
                                    );
  HLTBtagMatchObj.SetL25(L25BtagTriggerObjectFilterList,L25BtagTriggerObjectList);
  
  const std::vector<std::string>::const_iterator tSlotBtag = std::find(genericTriggerList.begin(), genericTriggerList.end(), theTrigNameBtag);
  if (tSlotBtag != genericTriggerList.end()) {
    std::cout << "Btag trigger found at slot " << tSlotBtag - genericTriggerList.begin() << std::endl;
  } else {
    std::cout << "Btag trigger not found in any slot" << std::endl;
    return 1;
  }

  const unsigned int tNumberBtag = tSlotBtag - genericTriggerList.begin();
 

  const std::vector<std::string>::const_iterator tSlotRef = std::find(genericTriggerList.begin(), genericTriggerList.end(), theTrigNameRef);
  if (tSlotRef != genericTriggerList.end()) {
    std::cout << "Ref trigger found at slot " << tSlotRef - genericTriggerList.begin() << std::endl;
  } else {
    std::cout << "Ref trigger not found in any slot" << std::endl;
    return 1;
  }
  const unsigned int tNumberRef = tSlotRef - genericTriggerList.begin();
 
  //switch off TTree branches that are currently not needed?
  //disable everything
  hbbtree->SetBranchStatus("*",0);
  hbbtree->SetBranchStatus("MatchedPartonFlavor",1);
  //enable certainly needed branches
  hbbtree->SetBranchStatus("nJetConstituents",1);
  hbbtree->SetBranchStatus("nJetChargedConstituents",1);
  hbbtree->SetBranchStatus("NumberOfJets",1);
  hbbtree->SetBranchStatus("JetEta",1);
  hbbtree->SetBranchStatus("JetPhi",1);
  hbbtree->SetBranchStatus("JetPt",1);

  hbbtree->SetBranchStatus("IsJetMatchedPartonExist",1);
  hbbtree->SetBranchStatus("HflContentJet",1);
  hbbtree->SetBranchStatus("PartonFlavorJet",1);

  hbbtree->SetBranchStatus("NumberOfPU",1);
  hbbtree->SetBranchStatus("NumberOfPUInTime",1);

  hbbtree->SetBranchStatus("CombSVBJetTag",1);
  hbbtree->SetBranchStatus("TCHEBJetTag",0);
  hbbtree->SetBranchStatus("TCHPBJetTag",0);
  hbbtree->SetBranchStatus("nTCHEBJetTag",0);
  hbbtree->SetBranchStatus("nTCHPBJetTag",0);
  hbbtree->SetBranchStatus("SVHPBJetTag",0);

 
  hbbtree->SetBranchStatus("L3*",1);
 
  hbbtree->SetBranchStatus("*L2*",0);
  //   hbbtree->SetBranchStatus("L2NumberOfJets",1);
  //   hbbtree->SetBranchStatus("L2JetPhi",1);
  hbbtree->SetBranchStatus("*L2*",1);
  hbbtree->SetBranchStatus("IsJetWithHltBtagBitPattern",1);
  hbbtree->SetBranchStatus("JetIDLoose",1);

  hbbtree->SetBranchStatus("PUJetMVA",0);
  hbbtree->SetBranchStatus("PUJetIDLoose",1);
  hbbtree->SetBranchStatus("PUJetIDMedium",0);
  hbbtree->SetBranchStatus("PUJetIDTight",0);

  hbbtree->SetBranchStatus("BContentJet5", 0);
  hbbtree->SetBranchStatus("CContentJet5", 0);
  hbbtree->SetBranchStatus("BContentJet3", 1);
  hbbtree->SetBranchStatus("CContentJet3", 1);

  //end of switch off TTree branches that are currently not needed?
  
  const double jetPtMax = 3500.;
  double maxEta = 1.7;

  const int nbtag = 4;  // number of offline btag working points

  const int nflav=5;  // number of flavor groups to be distinguished
  const std::string flavlabel[nflav] = {"udsg","c","b", "cc", "bb"};  // labels for flavor groups


  // define names of btag variants
  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };

  // corresponding discriminant cuts
  double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags




  // make pointers to appropriate btag discriminants
  float* theBJetTag[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    if (btdiscr[ibtag] == "TCHP") {
      theBJetTag[ibtag] = tcHPBJetTag;
    }  else if (btdiscr[ibtag] == "CSV") {
      theBJetTag[ibtag] = combSVBJetTag;
    }  else if (btdiscr[ibtag] == "SSVHP") {
      theBJetTag[ibtag] = svHPBJetTag;
    } else {
      std::cout << "invalid btag discriminant " << btdiscr[ibtag] << std::endl;
    }
  }


  // create the btag efficiency objects
  // offline btag object

  bTagEff* bTagEffOffline;

  bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_NOPU.root","offline",flavlabel,sbtag,nflav,nbtag);

  //online btag efficiencies 
  bTagEff* bTagEffOnline;
  
  if(onlinescenario == "normal")
    {
    //   bTagEffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V18_normal.root","online",flavlabel,sbtag,nflav,nbtag); 


    //   if(firstFile.Contains("52X")) {
//         bTagEffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V53_52X_normal.root","online",flavlabel,sbtag,nflav,nbtag); 
//       } else {
        bTagEffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_normal_NOPU.root","online",flavlabel,sbtag,nflav,nbtag); 
        // }
    }
  else { 
    std::cout << "not yet implemented." << std::endl;
    return 0;
  }
  
  // create the Root output file
  TFile* hout = new TFile("BtagClosureTest_RelOnlEff.root","recreate");

  TH1::AddDirectory(true);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();  // not really needed (we inherit from TH1)

  Int_t nentries = (Int_t) hbbtree->GetEntries();
  std::cout << "Number of events in ntuple: " << nentries << std::endl;

#include "Analysis/Utilities/interface/HbbNtuple.cc"

  // book the general histograms (before btag selection)
 
  std::cout << "book histograms" << std::endl;

  //binning in et
  const int n_bins_et = 33;
  const float binning_et[n_bins_et+1] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300, 400., 500.0, 800.0};
  




  const int nSelJet = 3;
  //online closure test first
  //Single offline tag
  hout->mkdir("RelOnlBtag_PtJetBtcut_SingleTag");
  hout->mkdir("RelOnlBtag_PtJetBtweight_SingleTag");

  

  hout->mkdir("RelOnlBtag_PtJetBtcut_FC_SingleTag");
  hout->mkdir("RelOnlBtag_PtJetBtweight_FC_SingleTag");


  hout->mkdir("RelOnlBtag_nChargedConstituentsBtcut_SingleTag");
  hout->mkdir("RelOnlBtag_nChargedConstituentsBtweight_SingleTag");


  hout->mkdir("RelOnlBtag_nConstituentsBtcut_SingleTag");
  hout->mkdir("RelOnlBtag_nConstituentsBtweight_SingleTag");



  hout->mkdir("RelOnlBtag_MjjBtcut_SingleTag");
  hout->mkdir("RelOnlBtag_MjjBtweight_SingleTag");
  hout->mkdir("RelOnlBtag_EtaJetBtcut_SingleTag");
  hout->mkdir("RelOnlBtag_EtaJetBtweight_SingleTag");

  hout->mkdir("RelOnlBtag_PhiJetBtcut_SingleTag");
  hout->mkdir("RelOnlBtag_PhiJetBtweight_SingleTag");


  hout->mkdir("RelOnlBtag_nL3Btcut_SingleTag");
  hout->mkdir("RelOnlBtag_nL3Btweight_SingleTag");
  
  TH1D* RelOnlBtag_PtJetBtcut_FC_SingleTag[nbtag][nSelJet][nflav];
  TH1D* RelOnlBtag_PtJetBtweight_FC_SingleTag[nbtag][nSelJet][nflav];


  TH1D* RelOnlBtag_PtJetBtcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_PtJetBtweight_SingleTag[nbtag][nSelJet];

  TH1D* RelOnlBtag_PhiJetBtcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_PhiJetBtweight_SingleTag[nbtag][nSelJet];


  TH1D* RelOnlBtag_MjjBtcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_MjjBtweight_SingleTag[nbtag][nSelJet];

  TH1D* RelOnlBtag_nChargedConstituentsBtcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_nChargedConstituentsBtweight_SingleTag[nbtag][nSelJet];
  
  TH1D* RelOnlBtag_nConstituentsBtcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_nConstituentsBtweight_SingleTag[nbtag][nSelJet];


  TH1D* RelOnlBtag_EtaJetBtcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_EtaJetBtweight_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_nL3Btcut_SingleTag[nbtag][nSelJet];
  TH1D* RelOnlBtag_nL3Btweight_SingleTag[nbtag][nSelJet];


  for (int ibtag=0; ibtag<nbtag; ++ibtag) {


    for (int iJ=0; iJ<nSelJet; ++iJ) {



      hout->cd("RelOnlBtag_nChargedConstituentsBtcut_SingleTag");
      RelOnlBtag_nChargedConstituentsBtcut_SingleTag[ibtag][iJ] = new TH1D(
                                                                           Form("RelOnlBtag_nChargedConstituentsBtcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                                           Form("nChargedConstituents onl. btag on jet %d %s",iJ,sbtag[ibtag].c_str()), 60,0.0, 60.0);
      
      
      hout->cd("RelOnlBtag_nChargedConstituentsBtweight_SingleTag");
      RelOnlBtag_nChargedConstituentsBtweight_SingleTag[ibtag][iJ] = new TH1D(
                                                                              Form("RelOnlBtag_nChargedConstituentsBtweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                                              Form("nChargedConstituents onl. btag on jet %d weighted %s",iJ,sbtag[ibtag].c_str()), 60,0.0, 60.0);



      hout->cd("RelOnlBtag_nConstituentsBtcut_SingleTag");
      RelOnlBtag_nConstituentsBtcut_SingleTag[ibtag][iJ] = new TH1D(
                                                                           Form("RelOnlBtag_nConstituentsBtcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                                           Form("nConstituents onl. btag on jet %d %s",iJ,sbtag[ibtag].c_str()), 80,0.0, 80.0);
      
      
      hout->cd("RelOnlBtag_nConstituentsBtweight_SingleTag");
      RelOnlBtag_nConstituentsBtweight_SingleTag[ibtag][iJ] = new TH1D(
                                                                              Form("RelOnlBtag_nConstituentsBtweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                                              Form("nConstituents onl. btag on jet %d weighted %s",iJ,sbtag[ibtag].c_str()), 80,0.0, 80.0);





    
      hout->cd("RelOnlBtag_MjjBtcut_SingleTag");
      RelOnlBtag_MjjBtcut_SingleTag[ibtag][iJ] = new TH1D(
                                                          Form("RelOnlBtag_MjjBtcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                          Form("Mjj onl. btag on jet %d %s",iJ,sbtag[ibtag].c_str()), 50,0.0, 500.0);
         
          
      hout->cd("RelOnlBtag_MjjBtweight_SingleTag");
      RelOnlBtag_MjjBtweight_SingleTag[ibtag][iJ] = new TH1D(
                                                             Form("RelOnlBtag_MjjBtweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                             Form("Mjj onl. btag on jet %d weighted %s",iJ,sbtag[ibtag].c_str()), 50,0.0, 500.0);
         
      for(int iflav = 0; iflav < nflav; iflav++) {
      hout->cd("RelOnlBtag_PtJetBtcut_FC_SingleTag");
      RelOnlBtag_PtJetBtcut_FC_SingleTag[ibtag][iJ][iflav] = new TH1D(
                                                            Form("RelOnlBtag_PtJetBtcut_FC_SingleTag_j%d_%s_%s",iJ,sbtag[ibtag].c_str(),flavlabel[iflav].c_str()),
                                                            Form("pt leading jet %d %s %s",iJ,sbtag[ibtag].c_str(),flavlabel[iflav].c_str()),n_bins_et, binning_et);
         
          
      hout->cd("RelOnlBtag_PtJetBtweight_FC_SingleTag");
      RelOnlBtag_PtJetBtweight_FC_SingleTag[ibtag][iJ][iflav] = new TH1D(
                                                                         Form("RelOnlBtag_PtJetBtweight_FC_SingleTag_j%d_%s_%s",iJ,sbtag[ibtag].c_str(),flavlabel[iflav].c_str()),
                                                               Form("pt leading jet %d weighted %s %s",iJ,sbtag[ibtag].c_str(),flavlabel[iflav].c_str()),n_bins_et, binning_et);
      
      }



      hout->cd("RelOnlBtag_PtJetBtcut_SingleTag");
      RelOnlBtag_PtJetBtcut_SingleTag[ibtag][iJ] = new TH1D(
                                                            Form("RelOnlBtag_PtJetBtcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                            Form("pt leading jet %d %s",iJ,sbtag[ibtag].c_str()),n_bins_et, binning_et);
         
          
      hout->cd("RelOnlBtag_PtJetBtweight_SingleTag");
      RelOnlBtag_PtJetBtweight_SingleTag[ibtag][iJ] = new TH1D(
                                                               Form("RelOnlBtag_PtJetBtweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                               Form("pt leading jet %d weighted %s",iJ,sbtag[ibtag].c_str()),n_bins_et, binning_et);
      




      
      hout->cd("RelOnlBtag_EtaJetBtcut_SingleTag");
      RelOnlBtag_EtaJetBtcut_SingleTag[ibtag][iJ] = new TH1D(
                                                             Form("RelOnlBtag_EtaJetBtcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                             Form("#eta leading jet %d %s",iJ,sbtag[ibtag].c_str()),30,-2.6,2.6);
         
          
      hout->cd("RelOnlBtag_EtaJetBtweight_SingleTag");
      RelOnlBtag_EtaJetBtweight_SingleTag[ibtag][iJ] = new TH1D(
                                                                Form("RelOnlBtag_EtaJetBtweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                                Form("#eta leading jet %d weighted %s",iJ,sbtag[ibtag].c_str()),30,-2.6,2.6);
         
          
    hout->cd("RelOnlBtag_PhiJetBtcut_SingleTag");
      RelOnlBtag_PhiJetBtcut_SingleTag[ibtag][iJ] = new TH1D(
                                                             Form("RelOnlBtag_PhiJetBtcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                             Form("#phi leading jet %d %s",iJ,sbtag[ibtag].c_str()),32,-1.0*TMath::Pi(), TMath::Pi());
         
          
      hout->cd("RelOnlBtag_PhiJetBtweight_SingleTag");
      RelOnlBtag_PhiJetBtweight_SingleTag[ibtag][iJ] = new TH1D(
                                                                Form("RelOnlBtag_PhiJetBtweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                                Form("#phi leading jet %d weighted %s",iJ,sbtag[ibtag].c_str()),32,-1.0*TMath::Pi(), TMath::Pi());
         
          





      hout->cd("RelOnlBtag_nL3Btcut_SingleTag");
      RelOnlBtag_nL3Btcut_SingleTag[ibtag][iJ] = new TH1D(
                                                             Form("RelOnlBtag_nL3Btcut_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                             Form(">=nL3, jet %d %s",iJ,sbtag[ibtag].c_str()),10, 0.0, 10.0);
         
          
      hout->cd("RelOnlBtag_nL3Btweight_SingleTag");
      RelOnlBtag_nL3Btweight_SingleTag[ibtag][iJ] = new TH1D(
                                                             Form("RelOnlBtag_nL3Btweight_SingleTag_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                             Form(">=nL3, jet %d weighted %s",iJ,sbtag[ibtag].c_str()),10, 0.0, 10.0);
      
          





    }
  }

  hout->cd();









  //
  hout->mkdir("RelOnlBtag_PtJetBtcut");
  hout->mkdir("RelOnlBtag_PtJetBtweight");
  hout->mkdir("RelOnlBtag_MjjBtcut");
  hout->mkdir("RelOnlBtag_MjjBtweight");
  hout->mkdir("RelOnlBtag_EtaJetBtcut");
  hout->mkdir("RelOnlBtag_EtaJetBtweight");


  
  TH1D* RelOnlBtag_PtJetBtcut[nbtag][nSelJet];
  TH1D* RelOnlBtag_PtJetBtweight[nbtag][nSelJet];
  TH1D* RelOnlBtag_MjjBtcut[nbtag][nSelJet];
  TH1D* RelOnlBtag_MjjBtweight[nbtag][nSelJet];
  TH1D* RelOnlBtag_EtaJetBtcut[nbtag][nSelJet];
  TH1D* RelOnlBtag_EtaJetBtweight[nbtag][nSelJet];
 

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int iJ=0; iJ<nSelJet; ++iJ) {
      
      hout->cd("RelOnlBtag_MjjBtcut");
      RelOnlBtag_MjjBtcut[ibtag][iJ] = new TH1D(
                                                Form("RelOnlBtag_MjjBtcut_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                Form("Mjj onl. btag on jet %d %s",iJ,sbtag[ibtag].c_str()), 50,0.0, 500.0);
      
      
      hout->cd("RelOnlBtag_MjjBtweight");
      RelOnlBtag_MjjBtweight[ibtag][iJ] = new TH1D(
                                                   Form("RelOnlBtag_MjjBtweight_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                   Form("Mjj onl. btag on jet %d weighted %s",iJ,sbtag[ibtag].c_str()), 50,0.0, 500.0);
      
      
      hout->cd("RelOnlBtag_PtJetBtcut");
      RelOnlBtag_PtJetBtcut[ibtag][iJ] = new TH1D(
                                                  Form("RelOnlBtag_PtJetBtcut_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                  Form("pt leading jet %d %s",iJ,sbtag[ibtag].c_str()),n_bins_et, binning_et);
      
      
      hout->cd("RelOnlBtag_PtJetBtweight");
      RelOnlBtag_PtJetBtweight[ibtag][iJ] = new TH1D(
                                                     Form("RelOnlBtag_PtJetBtweight_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                     Form("pt leading jet %d weighted %s",iJ,sbtag[ibtag].c_str()),n_bins_et, binning_et);
                
      hout->cd("RelOnlBtag_EtaJetBtcut");
      RelOnlBtag_EtaJetBtcut[ibtag][iJ] = new TH1D(
                                                   Form("RelOnlBtag_EtaJetBtcut_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                   Form("#eta leading jet %d %s",iJ,sbtag[ibtag].c_str()),30,-2.6,2.6);
         
          
      hout->cd("RelOnlBtag_EtaJetBtweight");
      RelOnlBtag_EtaJetBtweight[ibtag][iJ] = new TH1D(
                                                      Form("RelOnlBtag_EtaJetBtweight_j%d_%s",iJ,sbtag[ibtag].c_str()),
                                                      Form("#eta leading jet %d weighted %s",iJ,sbtag[ibtag].c_str()),30,-2.6,2.6);
        




 
    }
  }

  hout->cd();

  //end of online btagging closure test

  hout->mkdir("debug");
  hout->cd("debug");
  TH1D *h_check_HLTbtagMatch = new TH1D("h_check_HLTbtagMatch","h_check_HLTbtagMatch",2, 0.0, 2.0);
  TH1D *h_check_n_HLTbtagMatch = new TH1D("h_check_n_HLTbtagMatch","h_check_n_HLTbtagMatch",10, 0.0, 10.0);
  TH1D *h_check_n_HLTbtagMatch_selevent = new TH1D("h_check_n_HLTbtagMatch_selevent","h_check_n_HLTbtagMatch_selevent",10, 0.0, 10.0);
  
  TH1D *h_L2_consistency = new TH1D("h_L2_consistency","h_L2_consistency",2, 0.0, 2.0);
  

  const int required_n_jets = 3;
  double jetPtMinRelOnl[required_n_jets] = { 60.0, 53.0, 20.0};
  if(onlinescenario == "loose") {
    jetPtMinRelOnl[0] = 160.0;
    jetPtMinRelOnl[1] = 120.0;
    jetPtMinRelOnl[2] = 20.0;
    maxEta = 2.4;
  }
  
  // loop over the tree
  for (Int_t iE=0; iE<nentries; iE++) {

    hbbtree->GetEntry(iE);
   

    if(iE % 200000 == 0)
      //  if(iE % 200 == 0)
      {
        const double percent = 100.0*(double)iE/((double)nentries);
        std::cout << "process event " << iE << " - " << percent << " % done." << std::endl;
      }
  
    //test trigger object bit pattern
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      HLTBtagMatchObj.testinternalconsistency(isJetWithHltBtagBitPattern[iJet],isJetWithL25JetBitPattern[iJet],etaJet[iJet], phiJet[iJet],l25NumberOfJets,l25EtaJet,l25PhiJet,l3NumberOfBJets,l3EtaBJet,l3PhiBJet);
    }
    

    // compute the weight for this event
    float weight = 1;
    if ( _doMC && doLumiWeighting )  weight = lumiScaleFac;

    // trigger selection
    unsigned trgSelect = trgAccept;






    //check HLT btag match

    for(int iL3jet = 0; iL3jet < l3NumberOfBJets; iL3jet++) {
      int matched = 0;
      int iJtested = -1;
      for (int iJ=0; iJ<numberOfJets; ++iJ) {
      
       
     
        float delta_phi = l3PhiBJet[iL3jet] - phiJet[iJ];
              
        if(delta_phi > TMath::Pi()) delta_phi -= 2.0*TMath::Pi();
        if (delta_phi< -1.0*TMath::Pi()) delta_phi += 2.0*TMath::Pi();
              
        const float delta_eta = l3EtaBJet[iL3jet] - etaJet[iJ];
              
        const float delta_R = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
        if(delta_R < 0.5) {
          matched++;
          iJtested = iJ;
        }
       
      }
      h_check_n_HLTbtagMatch->Fill(matched,weight);
      bool matchedBOOL = matched >= 1 ? true: false;
      if(matched>= 1)
        {
          if(HLTBtagMatchObj.checkL3(isJetWithHltBtagBitPattern[iJtested],run) != matchedBOOL) {
            h_check_HLTbtagMatch->Fill(0.0, weight);
          }
          else {
            h_check_HLTbtagMatch->Fill(1.0, weight);
          }
        }

    }
    //end of check HLT btag match
























    //closure test for relative online btag efficiencies
    std::vector<int> leadingJetsRelOnl;


    // find set of leading jets
    int nJetRelOnl = 0;
    // loop over the jets
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      if ( ptJet[iJet] < 20. ) continue;
//       if(phiJet[iJet] > -0.5 && phiJet[iJet] < 0.5) 
//         continue;
      if(!puJetIDLoose[iJet]) continue;
     
      if(!jetIDLoose[iJet] ) continue;


      if (nJetRelOnl >= required_n_jets) break;
      if (! (fabs(etaJet[iJet])<maxEta) ) continue;
          
      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
          
      if ( (ptJet[iJet] > jetPtMinRelOnl[nJetRelOnl] ) && (ptJet[iJet] < jetPtMax) ) {
        leadingJetsRelOnl.push_back(iJet);
        ++nJetRelOnl;
      }
    }
    
    if (nJetRelOnl >= required_n_jets)
      {
        float deltaR12 = -1.0f;
        float deltaR23 = -1.0f;
        float deltaR13 = -1.0f;
        float deltaEta = 999.0f; // between the two leading jets

        if(nJetRelOnl >=2)
          {
            float dphij =  phiJet[leadingJetsRelOnl[1]] - phiJet[leadingJetsRelOnl[0]];
            if (dphij>3.1415926) dphij -= 2.0*3.1415926;
            if (dphij<-3.1415926) dphij += 2.0*3.1415926;
            
            float detaj = etaJet[leadingJetsRelOnl[1]] - etaJet[leadingJetsRelOnl[0]];
            deltaR12 = sqrt( dphij*dphij + detaj*detaj );
            deltaEta = TMath::Abs(detaj);
          }

        if(nJetRelOnl >= 3)
          {
            float dphij =  phiJet[leadingJetsRelOnl[2]] - phiJet[leadingJetsRelOnl[1]];                                                                                        
            if (dphij>3.1415926) dphij -= 2.0*3.1415926;                                                                                                           
            if (dphij<-3.1415926) dphij += 2.0*3.1415926;                                                                                                
            float detaj = etaJet[leadingJetsRelOnl[2]] - etaJet[leadingJetsRelOnl[1]];                                                                                         
            deltaR23 = sqrt( dphij*dphij + detaj*detaj );  
          }
        if(nJetRelOnl >= 3)
          {
            float dphij =  phiJet[leadingJetsRelOnl[0]] - phiJet[leadingJetsRelOnl[2]];                                                                                        
            if (dphij>3.1415926) dphij -= 2.0*3.1415926;                                                                                                           
            if (dphij<-3.1415926) dphij += 2.0*3.1415926;                                                                                                
            float detaj = etaJet[leadingJetsRelOnl[0]] - etaJet[leadingJetsRelOnl[2]];                                                                                         
            deltaR13 = sqrt( dphij*dphij + detaj*detaj );  
          }

        double deltaEtaCut12 = 1.7;
        if(onlinescenario == "loose")
          deltaEtaCut12 = 999999.0;

        if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) {

          //get status of reference and btag trigger bits
          bool referenceTriggerFired = true;
          bool btagTriggerFired = true;
          
          if (! (trgAccept & (1<<tNumberRef)))
            referenceTriggerFired = false;

          if(! (trgAccept & (1<<tNumberBtag)))
            btagTriggerFired = false;
 
          if(referenceTriggerFired)
            {

              // hltDoubleBJet40Eta2p6L1FastJet = HLT1CaloJet {
              //     bool saveTags = true
              //     double MinPt = 40.0
              //     int32 MinN = 2
              //     double MaxEta = 2.6
              //     double MinMass = -1.0
              //     InputTag inputTag = hltCaloJetL1FastJetCorrected
              //     double MinE = -1.0
              //     int32 triggerType = 86
              // }



              //now L2 selection
              int nSelectedL2Jets = 0;
                
              for (int iJet=0; iJet < l2NumberOfJets; ++iJet) {
                 
                //the selection
                if(TMath::Abs(l2EtaJet[iJet]) < 2.6 && l2PtJet[iJet] > (onlinescenario == "normal" ? 40.0 : 80.0 ))
                  nSelectedL2Jets++;
              }
              
              //we require at least 2 selected jets on L2

              int theFlav[required_n_jets];
              bool flavmatchingOK = true;
              for (int iJ=0; iJ<required_n_jets; ++iJ) {
                theFlav[iJ] = jetFlavorCodeNewCat(leadingJetsRelOnl[iJ]);
                if(theFlav[iJ] == -1)
                  flavmatchingOK = false; //only select events where the flavour of all three jets is defined
                
               //  if(!(isJetMatchedPartonExist[iJ]))
//                    flavmatchingOK = false;
              }
              if(btagTriggerFired && nSelectedL2Jets < 2) {
                h_L2_consistency->Fill(0.0,weight);
              } else {
                 h_L2_consistency->Fill(1.0,weight);
              }
              //apply L2 selection and check whether all jets have defined flavours
              if(flavmatchingOK && nSelectedL2Jets >= 2)
                {
                  //calculate mjj using the two leading jets
                  const float energyTot = energyJet[leadingJetsRelOnl[0]] + energyJet[leadingJetsRelOnl[1]];
                  const float pxTot = pxJet[leadingJetsRelOnl[0]] + pxJet[leadingJetsRelOnl[1]];
                  const float pyTot = pyJet[leadingJetsRelOnl[0]] + pyJet[leadingJetsRelOnl[1]];
                  const float pzTot = pzJet[leadingJetsRelOnl[0]] + pzJet[leadingJetsRelOnl[1]];
                
                  const float dijet_mass_sq = energyTot * energyTot - pxTot * pxTot - pyTot * pyTot - pzTot * pzTot;
                
                  const float dijet_mass = TMath::Sqrt( dijet_mass_sq );



                  //check HLT btag match

                  for(int iL3jet = 0; iL3jet < l3NumberOfBJets; iL3jet++) {
                    int matched = 0;
  
                    for (int iJ=0; iJ<numberOfJets; ++iJ) {
      
       
     
                      float delta_phi = l3PhiBJet[iL3jet] - phiJet[iJ];
              
                      if(delta_phi > TMath::Pi()) delta_phi -= 2.0*TMath::Pi();
                      if (delta_phi< -1.0*TMath::Pi()) delta_phi += 2.0*TMath::Pi();
              
                      const float delta_eta = l3EtaBJet[iL3jet] - etaJet[iJ];
              
                      const float delta_R = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
                      if(delta_R < 0.5) {
                        matched++;
                      }
       
                    }
                    h_check_n_HLTbtagMatch_selevent->Fill(matched,weight);
                           
                  }
                  //end of check HLT btag match












                  for (int ibtag=0; ibtag<nbtag; ++ibtag)
                    {
                      //check btag trigger
                      if(btagTriggerFired)
                        {
                          //check offline btagging for two leading jets
                          if( (theBJetTag[ibtag][leadingJetsRelOnl[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJetsRelOnl[1]]>btcut[ibtag])
                              && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJetsRelOnl[0]],isJetWithL25JetBitPattern[leadingJetsRelOnl[0]],run)
                              && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJetsRelOnl[1]],isJetWithL25JetBitPattern[leadingJetsRelOnl[1]],run)
                              && isJetMatchedPartonExist[leadingJetsRelOnl[0]] && isJetMatchedPartonExist[leadingJetsRelOnl[1]]
                              )
                            {
                              for (int iJ=0; iJ<required_n_jets; ++iJ) 
                                {
                                  if(
                                     1 == 1
                                     //HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJetsRelOnl[iJ]],isJetWithL25JetBitPattern[leadingJetsRelOnl[iJ]],run) 
                                     //&& isJetMatchedPartonExist[leadingJetsRelOnl[iJ]]
                                     //RedoHLTBtagMatch(leadingJetsRelOnl[iJ])
                                     )
                                    {
                                      //if third jet require also btagging for third jet
                                      if(iJ < 2 || (iJ == 2 && (theBJetTag[ibtag][leadingJetsRelOnl[2]] > btcut[ibtag])) ) 
                                        {
                                          RelOnlBtag_PtJetBtcut[ibtag][iJ]->Fill(ptJet[leadingJetsRelOnl[iJ]],weight);
                                          RelOnlBtag_EtaJetBtcut[ibtag][iJ]->Fill(etaJet[leadingJetsRelOnl[iJ]],weight);
                                          RelOnlBtag_MjjBtcut[ibtag][iJ]->Fill(dijet_mass,weight); 

                                        }
                                    }
                                } 
                            }



                       





                        }
                  
                      //now the weighting
                      if (_doMC) 
                        {
                          // determine jet flavor from MC info
                          double weight_offline[required_n_jets];
                          for (int iJ=0; iJ<required_n_jets; ++iJ) {
                            weight_offline[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJetsRelOnl[iJ]],etaJet[leadingJetsRelOnl[iJ]]);
                          }
                          double weight_online[required_n_jets];
                          for (int iJ=0; iJ<required_n_jets; ++iJ) {
                            weight_online[iJ] = bTagEffOnline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJetsRelOnl[iJ]],etaJet[leadingJetsRelOnl[iJ]]);
                          }

                       
                          //   for (int iJ=0; iJ<required_n_jets; ++iJ) {

                          //                           double totalweight = weight *
                          //                             (
                          //                              weight_offline[0] * weight_offline[1]
                          //                              * weight_online[iJ]
                          //                              );
                          
                          //                           if(iJ ==2)
                          //                             {
                          //                               totalweight = weight *
                          //                                 (
                          //                                  weight_offline[0] * weight_offline[1] * weight_offline[2]
                          //                                  * weight_online[iJ]
                          //                                  );
                          //                             }
                          
                          //                           RelOnlBtag_MjjBtweight[ibtag][iJ]->Fill(dijet_mass,totalweight); 
                          //                           RelOnlBtag_PtJetBtweight[ibtag][iJ]->Fill(ptJet[leadingJetsRelOnl[iJ]],totalweight);
                          //                           RelOnlBtag_EtaJetBtweight[ibtag][iJ]->Fill(etaJet[leadingJetsRelOnl[iJ]],totalweight);
                          //                         } 

                        
                          for (int iJ=0; iJ < 2; ++iJ) {
                            if( (theBJetTag[ibtag][leadingJetsRelOnl[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJetsRelOnl[1]]>btcut[ibtag])) {
                              
                              if(isJetMatchedPartonExist[leadingJetsRelOnl[0]] && isJetMatchedPartonExist[leadingJetsRelOnl[1]]) {


                              const double totalweight = weight *
                                (
                                 // ( iJ == 0 ? weight_offline[1] : weight_offline[0] )
                                 //1.0 * weight_online[leadingJetsRelOnl[iJ]]
                                 //( iJ == 2 ?  weight_offline[0] * weight_offline[1] * weight_offline[2] * weight_online[iJ] : weight_offline[0] * weight_offline[1] * weight_online[iJ] )

                                 
                                 weight_online[0] *  weight_online[1] 
                                 );
                              RelOnlBtag_MjjBtweight[ibtag][iJ]->Fill(dijet_mass,totalweight); 
                              RelOnlBtag_PtJetBtweight[ibtag][iJ]->Fill(ptJet[leadingJetsRelOnl[iJ]],totalweight);
                              RelOnlBtag_EtaJetBtweight[ibtag][iJ]->Fill(etaJet[leadingJetsRelOnl[iJ]],totalweight);
                            } 
                          }
                          }
                          
                        }
                    }















                  for (int ibtag=0; ibtag<nbtag; ++ibtag)
                    {
                      //check btag trigger
                      if(btagTriggerFired)
                        {
                          for (int iJ=0; iJ<required_n_jets; ++iJ) 
                            {
                              //check offline btagging
                              if( (theBJetTag[ibtag][leadingJetsRelOnl[iJ]]>btcut[ibtag]))
                                {
                                  if(
                                     HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJetsRelOnl[iJ]],isJetWithL25JetBitPattern[leadingJetsRelOnl[iJ]],run) 
                                     && isJetMatchedPartonExist[leadingJetsRelOnl[iJ]]
                                     //RedoHLTBtagMatch(leadingJetsRelOnl[iJ])
                                     )
                                    {
                                      RelOnlBtag_PtJetBtcut_FC_SingleTag[ibtag][iJ][theFlav[iJ]]->Fill(ptJet[leadingJetsRelOnl[iJ]],weight);
                                      RelOnlBtag_PtJetBtcut_SingleTag[ibtag][iJ]->Fill(ptJet[leadingJetsRelOnl[iJ]],weight);
                                      RelOnlBtag_EtaJetBtcut_SingleTag[ibtag][iJ]->Fill(etaJet[leadingJetsRelOnl[iJ]],weight);
                                      RelOnlBtag_PhiJetBtcut_SingleTag[ibtag][iJ]->Fill(phiJet[leadingJetsRelOnl[iJ]],weight);


                                      RelOnlBtag_MjjBtcut_SingleTag[ibtag][iJ]->Fill(dijet_mass,weight);
                                      if(ptJet[leadingJetsRelOnl[iJ]] > 53.0 && ptJet[leadingJetsRelOnl[iJ]] < 120.0) {
                                        RelOnlBtag_nChargedConstituentsBtcut_SingleTag[ibtag][iJ]->Fill(numberOfChargedConstituentsInJet[leadingJetsRelOnl[iJ]],weight);
                                        RelOnlBtag_nConstituentsBtcut_SingleTag[ibtag][iJ]->Fill(numberOfConstituentsInJet[leadingJetsRelOnl[iJ]],weight);
                                        
                                      }

                                      for(int ibin = 1; ibin <= RelOnlBtag_nL3Btcut_SingleTag[ibtag][iJ]->GetXaxis()->GetLast(); ibin++) {
                                        if(l3NumberOfBJets >= (ibin-1))
                                        RelOnlBtag_nL3Btcut_SingleTag[ibtag][iJ]->Fill(ibin-1,weight); 
                                      }
                                    }
                                } 
                            }
                        }
                  
                      //now the weighting
                      if (_doMC) 
                        {
                          // determine jet flavor from MC info
                          double weight_offline[required_n_jets];
                          for (int iJ=0; iJ<required_n_jets; ++iJ) {
                            weight_offline[iJ] = bTagEffOffline->eff(theFlav[iJ]
                                                                     ,sbtag[ibtag].c_str()
                                                                     ,ptJet[leadingJetsRelOnl[iJ]]
                                                                     ,etaJet[leadingJetsRelOnl[iJ]]);
                          }
                          double weight_online[required_n_jets];
                          for (int iJ=0; iJ<required_n_jets; ++iJ) {
                            weight_online[iJ] = bTagEffOnline->eff(theFlav[iJ]
                                                                   ,sbtag[ibtag].c_str()
                                                                   ,ptJet[leadingJetsRelOnl[iJ]]
                                                                   ,etaJet[leadingJetsRelOnl[iJ]]);
                          }

                        
                          for (int iJ=0; iJ < required_n_jets; ++iJ) {
                            if( (theBJetTag[ibtag][leadingJetsRelOnl[iJ]]>btcut[ibtag])   && isJetMatchedPartonExist[leadingJetsRelOnl[iJ]] ) {
                              const double totalweight = weight *
                                (
                                 1.0 * weight_online[iJ]// * weight_offline[iJ] 
                                 );
                              
                              if(ptJet[leadingJetsRelOnl[iJ]] > 53.0 && ptJet[leadingJetsRelOnl[iJ]] < 120.0) {
                                RelOnlBtag_nConstituentsBtweight_SingleTag[ibtag][iJ]->Fill(numberOfConstituentsInJet[leadingJetsRelOnl[iJ]],totalweight);
                                RelOnlBtag_nChargedConstituentsBtweight_SingleTag[ibtag][iJ]->Fill(numberOfChargedConstituentsInJet[leadingJetsRelOnl[iJ]],totalweight);
                                
                              }
                              

                              RelOnlBtag_MjjBtweight_SingleTag[ibtag][iJ]->Fill(dijet_mass,totalweight); 
                              RelOnlBtag_PtJetBtweight_SingleTag[ibtag][iJ]->Fill(ptJet[leadingJetsRelOnl[iJ]],totalweight);
                              RelOnlBtag_PtJetBtweight_FC_SingleTag[ibtag][iJ][theFlav[iJ]]->Fill(ptJet[leadingJetsRelOnl[iJ]],totalweight);
                              RelOnlBtag_EtaJetBtweight_SingleTag[ibtag][iJ]->Fill(etaJet[leadingJetsRelOnl[iJ]],totalweight);
                              RelOnlBtag_PhiJetBtweight_SingleTag[ibtag][iJ]->Fill(phiJet[leadingJetsRelOnl[iJ]],totalweight);
                              
                              for(int ibin = 1; ibin <= RelOnlBtag_nL3Btcut_SingleTag[ibtag][iJ]->GetXaxis()->GetLast(); ibin++) {
                                if(l3NumberOfBJets >= (ibin-1))
                                  RelOnlBtag_nL3Btweight_SingleTag[ibtag][iJ]->Fill((ibin-1),weight); 
                              }
                            } 
                          }
                        }
                    }
































                }
            }
        
        }
      }


  
    //end of closure test for relative online btag efficiencies


  } // loop over the tree
  // termination
  delete bTagEffOnline;
  delete bTagEffOffline;
  

  hout->Write();
  hout->Close();
}
