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
//#include "Analysis/Utilities/interface/FatSlimMixture.h"
#include "Analysis/Utilities/interface/TrigHistArray.h"
#include "Analysis/Utilities/interface/TrigHistArray2D.h"
#include "Analysis/Utilities/interface/getHbbMetadata.h"
#include "Analysis/Utilities/interface/checkHLTBTagMatch.h"
#include "Analysis/Utilities/interface/jetFlavorCode.h"
#include "Analysis/HbbMSSMAnalysis/test/NegativeBTag.C"
#include "Analysis/Utilities/interface/deltaR.h"

#include "Analysis/Utilities/interface/HbbNtuple.h"
#include "Analysis/Utilities/interface/HbbSelect.h"
#include "Analysis/Utilities/interface/TriggerRunSelection.h"

#include "Analysis/Utilities/interface/HbbSystControl.h"
#include "Analysis/Utilities/interface/HbbTrigger.h"
#include "Analysis/Utilities/interface/HbbTrigWeight.h"
//#include "Analysis/Utilities/interface/testJetIsolation.h"






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
  std::vector<std::string> triggerFilterList;
  std::vector<std::string> HLTBtagTriggerObjectFilterList;
  



  // extract generic trigger list and number of input events from ntuple files
  float lumiScaleFac = 1;
  getHbbMetadata(theFileColl,genericTriggerList,HLTBtagTriggerObjectList,1000.,lumiScaleFac,_doMC);
   
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
  } else if(onlinescenario == "loose") {
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BLooseL1FastJetFastPV");
  }
  
 
  
  f.AddFileInfoList((TCollection*) theFileColl.GetList());
  TTree* hbbtree = &f;

  // create trigger run selection object
  TriggerRunSelection trigSelector( Scenario, &genericTriggerList );
 //  checkHLTBtagMatch HLTBtagMatchObj(
//                                     triggerFilterList, 
//                                     HLTBtagTriggerObjectFilterList, 
//                                     HLTBtagTriggerObjectList,
//                                     doTrigSelect,
//                                     &trigSelector
//                                     );
  
 
  


 
 
  //switch off TTree branches that are currently not needed?
  //disable everything
  hbbtree->SetBranchStatus("*",0);
  hbbtree->SetBranchStatus("MatchedPartonFlavor",1);
  //enable certain needed branches
  hbbtree->SetBranchStatus("nJetConstituents",1);
  hbbtree->SetBranchStatus("NumberOfJets",1);
  hbbtree->SetBranchStatus("JetEta",1);
  hbbtree->SetBranchStatus("JetPhi",1);
  hbbtree->SetBranchStatus("JetPt",1);
  hbbtree->SetBranchStatus("JetPx",1);
  hbbtree->SetBranchStatus("JetPy",1);
  hbbtree->SetBranchStatus("JetPz",1);

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

 
  hbbtree->SetBranchStatus("L3*",0);
 
  hbbtree->SetBranchStatus("*L2*",0);
  //   hbbtree->SetBranchStatus("L2NumberOfJets",1);
  //   hbbtree->SetBranchStatus("L2JetPhi",1);
  hbbtree->SetBranchStatus("*L2*",0);
  hbbtree->SetBranchStatus("IsJetWithHltBtagBitPattern",0);


  hbbtree->SetBranchStatus("PUJetMVA",0);
  hbbtree->SetBranchStatus("PUJetIDLoose",1);
  hbbtree->SetBranchStatus("PUJetIDMedium",0);
  hbbtree->SetBranchStatus("PUJetIDTight",0);
  hbbtree->SetBranchStatus("JetIDLoose",1);
  hbbtree->SetBranchStatus("BContentJet5", 1);
  hbbtree->SetBranchStatus("CContentJet5", 1);
  hbbtree->SetBranchStatus("BContentJet3", 1);
  hbbtree->SetBranchStatus("CContentJet3", 1);

  //end of switch off TTree branches that are currently not needed?


  double jetPtMax = 3500.;
  double maxEta = 1.7;//2.6;
  // here we define dimensions and labels

//   const int nbtag = 4;  // number of offline btag working points

  
 
 

//   const int nflav=5;  // number of flavor groups to be distinguished
//   const std::string flavlabel[nflav] = {"udsg","c","b", "cc", "bb"};  // labels for flavor groups

//   // define names of btag variants
//   const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
//   string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };

//   // corresponding discriminant cuts
//   double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags





  const int nbtag = 1;  // number of offline btag working points

  
 
 

  const int nflav=5;  // number of flavor groups to be distinguished
  const std::string flavlabel[nflav] = {"udsg","c","b", "cc", "bb"};  // labels for flavor groups

  // define names of btag variants
  const std::string sbtag[nbtag] = {  "CSVT"};
  string btdiscr[nbtag] = {  "CSV"};

  // corresponding discriminant cuts
  double btcut[nbtag] = { 0.898 };  // working points for btags





  //

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

 
    //bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v13_partonmatchingmode3_bincentercorrected.root","offline",flavlabel,sbtag,nflav,nbtag);
    //bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v13_partonmatchingmode3_cone3.root","offline",flavlabel,sbtag,nflav,nbtag);
    // bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v14_partonmatchingmode3_cone5.root","offline",flavlabel,sbtag,nflav,nbtag);
    //bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v14_partonmatchingmode3_cone3.root","offline",flavlabel,sbtag,nflav,nbtag);
    // bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v14_partonmatchingmode3_cone3_bincentercorrected.root","offline",flavlabel,sbtag,nflav,nbtag);

  //   bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v14_partonmatchingmode3_cone3_jetPUidloose.root","offline",flavlabel,sbtag,nflav,nbtag);
   //  bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v15_partonmatchingmode3_cone3_jetPUidloose.root","offline",flavlabel,sbtag,nflav,nbtag);
//  bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v15_partonmatchingmode3_cone3_jetPUidloose.root","offline",flavlabel,sbtag,nflav,nbtag);
  
  // if(firstFile.Contains("53X")) {
//     bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v104_53X.root","offline",flavlabel,sbtag,nflav,nbtag);
//   }
//   else {
//     bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v104_52X.root","offline",flavlabel,sbtag,nflav,nbtag);
//   }

  //bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v104_53X_NOPU.root","offline",flavlabel,sbtag,nflav,nbtag);
  // bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v111_53X_NOPU.root","offline",flavlabel,sbtag,nflav,nbtag);
  bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_NOPU.root","offline",flavlabel,sbtag,nflav,nbtag);
  //FatSlimMixture flavmix("mediummass", "/afs/naf.desy.de/user/j/jbehr/public/cms/btag/flavourmixture.txt", bTagEffOffline);
  
  
  




  // create the Root output file
  TFile* hout = new TFile("BtagClosureTest.root","recreate");

  TH1::AddDirectory(true);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();  // not really needed (we inherit from TH1)

  Int_t nentries = (Int_t) hbbtree->GetEntries();
  std::cout << "Number of events in ntuple: " << nentries << std::endl;

#include "Analysis/Utilities/interface/HbbNtuple.cc"

  // book the general histograms (before btag selection)
 
  std::cout << "book histograms" << std::endl;


  //binning in et
  const int n_bins_et = 24;
  float binning_et[n_bins_et+1] = {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,270, 370., 500.0, 800.0};
  



  const int nSelJet = 3;
  
  hout->cd();


  TH1D *h_csv_discr = new TH1D("h_csv_discr","h_csv_discr", 1000, -20.0, 1.1);


  const int nhistoscuts = 3; //number of pt cut scenarios used for the 'general' histograms
  const double HistosJetCuts[3*nhistoscuts] = {
 //    20.0, -1.0, -1.0,
//     20.0, 20.0, -1.0,
//     20.0, 20.0, 20.0,
//     60.0, 20.0, 20.0,
    60.0, 53.0, 20.0,
    80, 70, 20,
    160, 120, 20
  //   80.0, 20.0, -1.0,
//     80.0, 70.0, -1.0,
//     80.0, 70.0, 20.0
  };

  const double etaCut[nhistoscuts] = {
    1.7,
    1.7,
    2.4
  };
  const double deltaEtaCut[nhistoscuts] = {
    1.7,
    1.7,
    99999999.0
  };

  const std::string sc[nhistoscuts] = {
    "MediumMass2012",
    "HighMass2012",
    "VeryHighMass2012"
  };
 
  hout->mkdir("SingleTag_OffBtagPtJetBtcut");
  hout->mkdir("SingleTag_OffBtagPtJetBtweight");

  hout->mkdir("SingleTag_OffBtagPtJetBtcutFc");
  hout->mkdir("SingleTag_OffBtagPtJetBtweightFc");

  const int nstat = 3;
  TrigHistArray* SingleTag_OffBtagptJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagptJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* SingleTag_OffBtagptJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagptJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* SingleTag_OffBtagptJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagptJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];
  
 

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("SingleTag_OffBtagPtJetBtcut");
          SingleTag_OffBtagptJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                            Form("OffBtag_ptj%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                            Form("pt leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          SingleTag_OffBtagptJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                      Form("OffBtag_ptj%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                      Form("pt leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          SingleTag_OffBtagptJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                               Form("OffBtag_ptj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                               Form("pt leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          
          hout->cd("SingleTag_OffBtagPtJetBtweight");
          for(int istat = 0; istat< nstat; istat++)
            {
              SingleTag_OffBtagptJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                        Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                        Form("pt leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);
              SingleTag_OffBtagptJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                  Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                  Form("pt leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);
          
              SingleTag_OffBtagptJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                           Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                           Form("pt leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);
          

            }
        }
      }
    }


  //
  
  

  //now mjj

  hout->mkdir("DoubleTag_OffBtagMjjBtcut");
  hout->mkdir("DoubleTag_OffBtagMjjBtweight");

 
  TrigHistArray* DoubleTag_OffBtagMjjBt_deltaRcut_deltaEta[nbtag][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagMjjBtw_deltaRcut_deltaEta[nbtag][nhistoscuts];

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
       

        hout->cd("DoubleTag_OffBtagMjjBtcut");
                                                                                    
        DoubleTag_OffBtagMjjBt_deltaRcut_deltaEta[ibtag][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_Mjj_bt%s_JetCuts%d_deltaRcut_deltaEta",sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("Mjj 2*%s,JetCuts%d",sbtag[ibtag].c_str(), iJetCuts),60,0.0, 500.0);
          
        

        hout->cd("DoubleTag_OffBtagMjjBtweight");
        DoubleTag_OffBtagMjjBtw_deltaRcut_deltaEta[ibtag][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                        Form("OffBtag_Mjj_ww%s_JetCuts%d_deltaRcut_deltaEta",sbtag[ibtag].c_str(), iJetCuts),
                                                                                        Form("Mjj 2*%s weighted,JetCuts%d",sbtag[ibtag].c_str(), iJetCuts),60,0.0, 500.0);
          
             
      }
      
    }


  //nJetConstituens
  
  hout->mkdir("DoubleTag_OffBtag_nJetConstituents_BtcutFc");
  hout->mkdir("DoubleTag_OffBtag_nJetConstituents_BtweightFc");

  TrigHistArray* DoubleTag_OffBtag_nJetConstituents_Bt_Fc_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][5];
  TrigHistArray* DoubleTag_OffBtag_nJetConstituents_Btw_Fc_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat][5];

  for(int iflav = 0; iflav < 5; iflav++)
    {
      for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
        {
          for (int ibtag=0; ibtag<nbtag; ++ibtag) {
            for (int iJ=0; iJ<nSelJet; ++iJ) {

              hout->cd("DoubleTag_OffBtag_nJetConstituents_BtcutFc");
          
              DoubleTag_OffBtag_nJetConstituents_Bt_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][iflav] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                                          Form("OffBtag_nJetConstituents_j%dbt%s_JetCuts%d_flav%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,iflav),
                                                                                                                          Form("nJetConstituents leading jet %d 2*%s,JetCuts%d,%d",iJ,sbtag[ibtag].c_str(), iJetCuts,iflav),60, 0.0, 60.0);
          
          
              for(int istat = 0; istat< nstat; istat++)
                {

                  hout->cd("DoubleTag_OffBtag_nJetConstituents_BtweightFc");
            
          
                  DoubleTag_OffBtag_nJetConstituents_Btw_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat][iflav] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                                                      Form("OffBtag_nJetConstituents_j%dww%s_JetCuts%d_stat%d_flav%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat,iflav),
                                                                                                                                      Form("nJetConstituents leading jet %d weighted as 2*%s,JetCuts%d_stat%d,%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat,iflav),60, 0.0, 60.0);
                }
            }
          }
        }

    }


  //

  hout->mkdir("DoubleTag_OffBtagPtJetBtcut");
  hout->mkdir("DoubleTag_OffBtagPtJetBtweight");

  hout->mkdir("DoubleTag_OffBtagPtJetBtcutFc");
  hout->mkdir("DoubleTag_OffBtagPtJetBtweightFc");

  TrigHistArray* DoubleTag_OffBtagptJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagptJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* DoubleTag_OffBtagptJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagptJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];

  TrigHistArray* DoubleTag_OffBtagptJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagptJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("DoubleTag_OffBtagPtJetBtcut");
          DoubleTag_OffBtagptJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                            Form("OffBtag_ptj%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                            Form("pt leading jet %d 2*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          DoubleTag_OffBtagptJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                      Form("OffBtag_ptj%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                      Form("pt leading jet %d 2*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          DoubleTag_OffBtagptJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                               Form("OffBtag_ptj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                               Form("pt leading jet %d 2*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          
          
          for(int istat = 0; istat< nstat; istat++)
            {

              hout->cd("DoubleTag_OffBtagPtJetBtweight");
              DoubleTag_OffBtagptJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                        Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                        Form("pt leading jet %d weighted as 2*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);
              DoubleTag_OffBtagptJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                  Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                  Form("pt leading jet %d weighted as 2*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);
          
              DoubleTag_OffBtagptJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                           Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                           Form("pt leading jet %d weighted as 2*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);
            }
        }
      }
    }

  TrigHistArray* DoubleTag_OffBtagptJetBt_Fc_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][5];
  TrigHistArray* DoubleTag_OffBtagptJetBtw_Fc_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat][5];

  for(int iflav = 0; iflav < 5; iflav++)
    {
      for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
        {
          for (int ibtag=0; ibtag<nbtag; ++ibtag) {
            for (int iJ=0; iJ<nSelJet; ++iJ) {

              hout->cd("DoubleTag_OffBtagPtJetBtcutFc");
          
              DoubleTag_OffBtagptJetBt_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][iflav] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                             Form("OffBtag_ptj%dbt%s_JetCuts%d_flav%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,iflav),
                                                                                                             Form("pt leading jet %d 2*%s,JetCuts%d,%d",iJ,sbtag[ibtag].c_str(), iJetCuts,iflav),n_bins_et,binning_et);
          
          
              for(int istat = 0; istat< nstat; istat++)
                {

                  hout->cd("DoubleTag_OffBtagPtJetBtweightFc");
            
          
                  DoubleTag_OffBtagptJetBtw_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat][iflav] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                                         Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_flav%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat,iflav),
                                                                                                                         Form("pt leading jet %d weighted as 2*%s,JetCuts%d_stat%d,%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat,iflav),n_bins_et,binning_et);
                }
            }
          }
        }

    }

  //---
  hout->mkdir("TripleTag_OffBtagPtJetBtcut");
  hout->mkdir("TripleTag_OffBtagPtJetBtweight");

  hout->mkdir("TripleTag_OffBtagPtJetBtcutFc");
  hout->mkdir("TripleTag_OffBtagPtJetBtweightFc");

  TrigHistArray* TripleTag_OffBtagptJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagptJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagptJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagptJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagptJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagptJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagptJetBtw_Fc_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nflav];

  TrigHistArray* TripleTag_OffBtagptJetBt_Fc_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nflav];



  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("TripleTag_OffBtagPtJetBtcut");
          TripleTag_OffBtagptJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                            Form("OffBtag_ptj%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                            Form("pt leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          TripleTag_OffBtagptJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                      Form("OffBtag_ptj%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                      Form("pt leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
          TripleTag_OffBtagptJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                               Form("OffBtag_ptj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                               Form("pt leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),n_bins_et,binning_et);
           hout->cd("TripleTag_OffBtagPtJetBtcutFc");
          

           for(int iflav = 0; iflav < nflav; iflav++)
            {
              TripleTag_OffBtagptJetBt_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][iflav] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                              Form("OffBtag_ptj%dbt%s_JetCuts%d_flav%s_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,flavlabel[iflav].c_str()),
                                                                                                             Form("pt leading jet %d %s,JetCuts%d,%s",iJ,sbtag[ibtag].c_str(), iJetCuts,flavlabel[iflav].c_str()),n_bins_et,binning_et);
            }


          hout->cd("TripleTag_OffBtagPtJetBtweightFc");
          

          for(int iflav = 0; iflav < nflav; iflav++)
            {
              
              TripleTag_OffBtagptJetBtw_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][iflav] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                              Form("OffBtag_ptj%dww%s_JetCuts%d_flav%s_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,flavlabel[iflav].c_str()),
                                                                                                              Form("pt leading jet %d weighted as 3*%s,JetCuts%d,%s",iJ,sbtag[ibtag].c_str(), iJetCuts,flavlabel[iflav].c_str()),n_bins_et,binning_et);
              
            }


          for(int istat = 0; istat< nstat; istat++)
            {
              hout->cd("TripleTag_OffBtagPtJetBtweight");
              TripleTag_OffBtagptJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                        Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                        Form("pt leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);

              TripleTag_OffBtagptJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                  Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                  Form("pt leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);

              TripleTag_OffBtagptJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                           Form("OffBtag_ptj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                           Form("pt leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),n_bins_et,binning_et);



            }}
      }
    }

 







  //----------- now eta and delta R



  TrigHistArray* SingleTag_OffBtagetaJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagetaJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* SingleTag_OffBtagetaJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagetaJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* SingleTag_OffBtagetaJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];

  TrigHistArray* SingleTag_OffBtagdeltaRJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagdeltaRJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];
 

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("SingleTag_OffBtagPtJetBtcut");
          SingleTag_OffBtagetaJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                             Form("OffBtag_etaj%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                             Form("eta leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          SingleTag_OffBtagetaJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_etaj%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("eta leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          SingleTag_OffBtagetaJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                Form("OffBtag_etaj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                Form("eta leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          
          SingleTag_OffBtagdeltaRJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_deltaRj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                   Form("#Delta R leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, 0.0, 3.0);
          
          for(int istat = 0; istat< nstat; istat++)
            {

              hout->cd("SingleTag_OffBtagPtJetBtweight");
              SingleTag_OffBtagetaJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                         Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                         Form("eta leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);
              SingleTag_OffBtagetaJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                   Form("eta leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);
              SingleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                            Form("eta leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);



              SingleTag_OffBtagdeltaRJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_deltaRj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                               Form("#deltaR leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, 0.0, 3.0);
              



            }
        }}
    }


  //



  TrigHistArray* DoubleTag_OffBtagetaJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagetaJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* DoubleTag_OffBtagetaJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagetaJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* DoubleTag_OffBtagetaJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("DoubleTag_OffBtagPtJetBtcut");
          DoubleTag_OffBtagetaJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                             Form("OffBtag_etaj%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                             Form("eta leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          DoubleTag_OffBtagetaJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_etaj%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("eta leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          
          DoubleTag_OffBtagetaJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                Form("OffBtag_etaj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                Form("eta leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          for(int istat = 0; istat< nstat; istat++)
            {
              hout->cd("DoubleTag_OffBtagPtJetBtweight");
              DoubleTag_OffBtagetaJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                         Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                         Form("eta leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);
              DoubleTag_OffBtagetaJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                   Form("eta leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);
              DoubleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                            Form("eta leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);
          
            }
        }
      }
    }





  TrigHistArray* TripleTag_OffBtagetaJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagetaJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagetaJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagetaJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagetaJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("TripleTag_OffBtagPtJetBtcut");
          TripleTag_OffBtagetaJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                             Form("OffBtag_etaj%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                             Form("eta leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          TripleTag_OffBtagetaJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_etaj%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("eta leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          TripleTag_OffBtagetaJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                Form("OffBtag_etaj%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                Form("eta leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -2.6,2.6);
          for(int istat = 0; istat< nstat; istat++)
            {

              hout->cd("TripleTag_OffBtagPtJetBtweight");
              TripleTag_OffBtagetaJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                         Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                         Form("eta leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);

              TripleTag_OffBtagetaJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                   Form("eta leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);
              TripleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_etaj%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                            Form("eta leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -2.6,2.6);



            }
        }
      }



    }
 

  //now phi



  TrigHistArray* SingleTag_OffBtagphiJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagphiJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* SingleTag_OffBtagphiJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagphiJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* SingleTag_OffBtagphiJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* SingleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];



  TrigHistArray* DoubleTag_OffBtagphiJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagphiJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* DoubleTag_OffBtagphiJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagphiJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* DoubleTag_OffBtagphiJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* DoubleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];



 

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("SingleTag_OffBtagPtJetBtcut");
          SingleTag_OffBtagphiJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                             Form("OffBtag_phij%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                             Form("#phi leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          SingleTag_OffBtagphiJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_phij%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("#phi leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          SingleTag_OffBtagphiJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                Form("OffBtag_phij%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                Form("#phi leading jet %d 1*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          for(int istat = 0; istat < nstat; istat++)
            {
              hout->cd("SingleTag_OffBtagPtJetBtweight");
              SingleTag_OffBtagphiJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                         Form("OffBtag_phij%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                         Form("#phi leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());
              SingleTag_OffBtagphiJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_phij%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                   Form("#phi leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());
              SingleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_phij%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                            Form("#phi leading jet %d weighted as 1*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());

            }   }
      }
    }





  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("DoubleTag_OffBtagPtJetBtcut");
          DoubleTag_OffBtagphiJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                             Form("OffBtag_phij%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                             Form("#phi leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          DoubleTag_OffBtagphiJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_phij%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("#phi leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          
          DoubleTag_OffBtagphiJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                Form("OffBtag_phij%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                Form("#phi leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          for(int istat = 0; istat< nstat; istat++)
            {
              hout->cd("DoubleTag_OffBtagPtJetBtweight");
              DoubleTag_OffBtagphiJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                         Form("OffBtag_phij%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                         Form("#phi leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());
              DoubleTag_OffBtagphiJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_phij%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                   Form("#phi leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());
              DoubleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_phij%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                            Form("#phi leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());
          

            }
        }
      }

    }



  TrigHistArray* TripleTag_OffBtagphiJetBt[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagphiJetBtw[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagphiJetBt_deltaRcut[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagphiJetBtw_deltaRcut[nbtag][nSelJet][nhistoscuts][nstat];
  TrigHistArray* TripleTag_OffBtagphiJetBt_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts];
  TrigHistArray* TripleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[nbtag][nSelJet][nhistoscuts][nstat];

  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        for (int iJ=0; iJ<nSelJet; ++iJ) {

          hout->cd("TripleTag_OffBtagPtJetBtcut");
          TripleTag_OffBtagphiJetBt[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                             Form("OffBtag_phij%dbt%s_JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                             Form("#phi leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          TripleTag_OffBtagphiJetBt_deltaRcut[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                       Form("OffBtag_phij%dbt%s_JetCuts%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                       Form("#phi leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          TripleTag_OffBtagphiJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                Form("OffBtag_phij%dbt%s_JetCuts%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts),
                                                                                                Form("#phi leading jet %d 3*%s,JetCuts%d",iJ,sbtag[ibtag].c_str(), iJetCuts),32, -1.0*TMath::Pi(), TMath::Pi());
          
          for(int istat = 0; istat< nstat; istat++)
            {
              hout->cd("TripleTag_OffBtagPtJetBtweight");
              TripleTag_OffBtagphiJetBtw[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                         Form("OffBtag_phij%dww%s_JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                         Form("#phi leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());

              TripleTag_OffBtagphiJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                   Form("OffBtag_phij%dww%s_JetCuts%d_stat%d_deltaRcut",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                   Form("#phi leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());
              TripleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                                                            Form("OffBtag_phij%dww%s_JetCuts%d_stat%d_deltaRcut_deltaEta",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),
                                                                                                            Form("#phi leading jet %d weighted as 3*%s,JetCuts%d_stat%d",iJ,sbtag[ibtag].c_str(), iJetCuts,istat),32, -1.0*TMath::Pi(), TMath::Pi());

            }

        }
      }
    }
  
  hout->cd();
  hout->mkdir("fat_jet_efficiencies");
  TH1D *h_fat_cc_jet_efficiencies[nhistoscuts];
  TH1D *h_fat_bb_jet_efficiencies[nhistoscuts];
  
  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      hout->cd("fat_jet_efficiencies");
      std::string s(Form("h_fat_cc_jet_efficiencies_%d",iJetCuts));
      h_fat_cc_jet_efficiencies[iJetCuts] = new TH1D(s.c_str(), s.c_str(), 300, 0.0, 4.0);
    }
   
  for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
    {
      hout->cd("fat_jet_efficiencies");
      std::string s(Form("h_fat_bb_jet_efficiencies_%d",iJetCuts));
      h_fat_bb_jet_efficiencies[iJetCuts] = new TH1D(s.c_str(), s.c_str(), 300, 0.0, 4.0);
    }
   
   


  // loop over the tree
  for (Int_t iE=0; iE<nentries; iE++) {

    hbbtree->GetEntry(iE);
   

    if(iE % 200000 == 0)
      //  if(iE % 200 == 0)
      {
        const double percent = 100.0*(double)iE/((double)nentries);
        std::cout << "process event " << iE << " - " << percent << " % done."<< std::endl;
      }
    
  

    // compute the weight for this event
    float weight = 1;
    if ( _doMC && doLumiWeighting )  weight = lumiScaleFac;

    // trigger selection
    unsigned trgSelect = trgAccept;

   

    for(int iJetCuts = 0; iJetCuts < nhistoscuts; iJetCuts++)
      {
        std::vector<int> leadingJets;

        // find set of leading jets
        int nJet = 0;
        // loop over the jets
        for (int iJet=0; iJet<numberOfJets; ++iJet) {
          if (nJet >= nSelJet) break;
          //if( combSVBJetTag[iJet] < 0.0 ) continue; //REMOVE! ITS JUST FOR DEBUGGING!
          if (! (fabs(etaJet[iJet])< etaCut[iJetCuts]) ) continue;
          if ( ptJet[iJet] < 20. ) continue;
          if(!puJetIDLoose[iJet])
            continue;
          if(!jetIDLoose[iJet] ) continue;
          

          if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
          
          if ( (ptJet[iJet] > HistosJetCuts[(unsigned int)(3*iJetCuts+nJet)] ) && (ptJet[iJet] < jetPtMax) ) {
            leadingJets.push_back(iJet);
            ++nJet;
          }
        }
        int nJetsCut = 0;

        for (int inJetsCut=0; inJetsCut<nSelJet; ++inJetsCut) {
          if(HistosJetCuts[(unsigned int)(3*iJetCuts+inJetsCut)] > 0.0)
            nJetsCut++;
        }

       


        if (nJet < nJetsCut) continue;
        const int nSelJet = nJetsCut;
        
        float deltaR12 = -1.0f;
        float deltaR23 = -1.0f;
        float deltaR13 = -1.0f;
        float deltaEta = 999.0f; // between the two leading jets
        if(nJet >=2)
          {
            float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
            if (dphij>3.1415926) dphij -= 2.0*3.1415926;
            if (dphij<-3.1415926) dphij += 2.0*3.1415926;
            
            float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[0]];
            deltaR12 = sqrt( dphij*dphij + detaj*detaj );
            deltaEta = TMath::Abs(detaj);
          }

        if(nJet >= 3)
          {
            float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[1]];                                                                                        
            if (dphij>3.1415926) dphij -= 2.0*3.1415926;                                                                                                           
            if (dphij<-3.1415926) dphij += 2.0*3.1415926;                                                                                                
            float detaj = etaJet[leadingJets[2]] - etaJet[leadingJets[1]];                                                                                         
            deltaR23 = sqrt( dphij*dphij + detaj*detaj );  
          }
        if(nJet >= 3)
          {
            float dphij =  phiJet[leadingJets[0]] - phiJet[leadingJets[2]];                                                                                        
            if (dphij>3.1415926) dphij -= 2.0*3.1415926;                                                                                                           
            if (dphij<-3.1415926) dphij += 2.0*3.1415926;                                                                                                
            float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[2]];                                                                                         
            deltaR13 = sqrt( dphij*dphij + detaj*detaj );  
          }
        



        // determine jet flavor from MC info
        int theFlav[nJetsCut];
       

        for (int iJ=0; iJ<nJetsCut; ++iJ) {
          theFlav[iJ] = -1;
        }
        for (int iJ=0; iJ<nJetsCut; ++iJ) {
          theFlav[iJ] = jetFlavorCodeNewCat(leadingJets[iJ]);
          if(theFlav[iJ] != -1 ) {


            const double fraction0 = 0.2;
            const double fraction1= 1.0-fraction0;
            if(theFlav[iJ] == 3) {
              const double cc_eff = bTagEffOffline->eff(3,"CSVT",ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
              //c eff at half pt
              const double c_eff_0 = bTagEffOffline->eff(1,"CSVT",ptJet[leadingJets[iJ]]*fraction0,etaJet[leadingJets[iJ]]);
              const double c_eff_1 = bTagEffOffline->eff(1,"CSVT",ptJet[leadingJets[iJ]]*fraction1,etaJet[leadingJets[iJ]]);
              
              const double filtered_cc_eff = 1.0 - (1.0 - c_eff_0) * (1.0 - c_eff_1);

              const double eff_ratio = filtered_cc_eff > 0.0 ? cc_eff / filtered_cc_eff : 0.0;
              h_fat_cc_jet_efficiencies[iJetCuts]->Fill(eff_ratio,weight);
            } else if(theFlav[iJ] == 4) {
              const double bb_eff = bTagEffOffline->eff(4,"CSVT",ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
              //c eff at half pt
              const double b_eff_0 = bTagEffOffline->eff(2,"CSVT",ptJet[leadingJets[iJ]]*fraction0,etaJet[leadingJets[iJ]]);
              const double b_eff_1 = bTagEffOffline->eff(2,"CSVT",ptJet[leadingJets[iJ]]*fraction1,etaJet[leadingJets[iJ]]);
              
              const double filtered_bb_eff = 1.0 - (1.0 - b_eff_0) * (1.0 - b_eff_1);
              const double eff_ratio = filtered_bb_eff > 0.0 ? bb_eff / filtered_bb_eff : 0.0;
              h_fat_bb_jet_efficiencies[iJetCuts]->Fill(eff_ratio,weight);
            }
            


          }
        }
        if(nSelJet >= 1)
          {
            //if((theFlav[0] != -1)  
            if((theFlav[0] != -1 && theFlav[1] != -1 && theFlav[2] != -1))
              {
                for (int ibtag=0; ibtag<nbtag; ++ibtag)
                  {
                    // require ONE btag for the  leading jets
                    //if( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) )
                    if(1==1)
                      {
                        for (int iJ=0; iJ<nSelJet; ++iJ) {
                          if( (theBJetTag[ibtag][leadingJets[iJ]]>btcut[ibtag])  && isJetMatchedPartonExist[leadingJets[iJ]] ) {

                            //double CalculateDeltaR (double phi1, double phi2, double eta1, double eta2)
                            double mindeltaR = 9999.0;
                            for (int iJet=0; iJet<numberOfJets; ++iJet) {
                              if(iJet != leadingJets[iJ] ) {
                                const double dR = CalculateDeltaR(phiJet[leadingJets[iJ]], phiJet[iJet], etaJet[leadingJets[iJ]], etaJet[iJet]);
                                if(dR < mindeltaR) {
                                  mindeltaR = dR;
                                }
                              }
                            }

                            



                            SingleTag_OffBtagptJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                            SingleTag_OffBtagetaJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                            SingleTag_OffBtagphiJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                            
                            if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0)
                              {
                                SingleTag_OffBtagptJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                                SingleTag_OffBtagetaJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                                SingleTag_OffBtagphiJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                                
                                if(deltaEta < deltaEtaCut[iJetCuts])
                                  {
                                    SingleTag_OffBtagptJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                                    SingleTag_OffBtagetaJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                                    SingleTag_OffBtagphiJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                                    
                                  
                                    
                                    SingleTag_OffBtagdeltaRJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,mindeltaR,weight);

                                  }
                              }
                          }
                        }
                      }
                    if (_doMC) {
              
                      double wbt[nJetsCut];
                      double wbterror[nJetsCut];
              
                      for (int iJ=0; iJ<nJetsCut; ++iJ) {
                        if(theFlav[iJ] != -1) {
                          wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                          wbterror[iJ] = bTagEffOffline->erreff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                        }
                        else {
                          std::cout << "error in single tag closure test" << std::endl;
                        }
                      }

                     
                      


                      for(int istat = 0; istat < 1; istat++)
                        {
                        
                          for (int iJ=0; iJ<nSelJet; ++iJ) {
                            if ( isJetMatchedPartonExist[leadingJets[iJ]]){


                              double mindeltaR = 9999.0;
                              for (int iJet=0; iJet<numberOfJets; ++iJet) {
                                if(iJet != leadingJets[iJ] ) {
                                  const double dR = CalculateDeltaR(phiJet[leadingJets[iJ]], phiJet[iJet], etaJet[leadingJets[iJ]], etaJet[iJet]);
                                  if(dR < mindeltaR) {
                                    mindeltaR = dR;
                                  }
                                }
                              }




                            const double error = ( istat == 1 ? +1.0 : -1.0) * wbterror[iJ];
                            double wtotDoubleTag = weight * 
                              (
                               // wbt[iJ] + (istat == 0 ? 0.0 : error)
                               wbt[iJ] + (istat == 0 ? 0.0 : error)
                               ); //single tag weight?
                            if(wtotDoubleTag <= 0.0)
                              {
                                std::cout << "error in single tag closure test" << std::endl;
                                throw std::exception();
                              wtotDoubleTag = 1.0;
                              }



                            SingleTag_OffBtagptJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                            SingleTag_OffBtagetaJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                            SingleTag_OffBtagphiJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);

                            if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0)
                              {
                                SingleTag_OffBtagetaJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                                SingleTag_OffBtagptJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                                SingleTag_OffBtagphiJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);

                                if(deltaEta < deltaEtaCut[iJetCuts])
                                  {
                                    SingleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                                    SingleTag_OffBtagptJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                                    SingleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);

                                    SingleTag_OffBtagdeltaRJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,mindeltaR,wtotDoubleTag);



                                    
                                  }
                              }
                              }
                          }
                        }
                    }
          
                  }
              }
          }
        //
        

        if(nSelJet >= 2)
          {

            if((theFlav[0] != -1 && theFlav[1] != -1 && theFlav[2] != -1))
              {
                //calculate mjj
                const float energyTot = energyJet[leadingJets[0]] + energyJet[leadingJets[1]];
                const float pxTot = pxJet[leadingJets[0]] + pxJet[leadingJets[1]];
                const float pyTot = pyJet[leadingJets[0]] + pyJet[leadingJets[1]];
                const float pzTot = pzJet[leadingJets[0]] + pzJet[leadingJets[1]];
                
                const float dijet_mass_sq = energyTot * energyTot - pxTot * pxTot - pyTot * pyTot - pzTot * pzTot;
                
                const float dijet_mass = sqrt( dijet_mass_sq );
                for (int ibtag=0; ibtag<nbtag; ++ibtag)
                  {
                    // require ONE btag for the two leading jets
                    if( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag])
                         && isJetMatchedPartonExist[leadingJets[0]] && isJetMatchedPartonExist[leadingJets[1]]
                        )
                      {
                        for (int iJ=0; iJ<nSelJet; ++iJ) {
                          DoubleTag_OffBtagptJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                          DoubleTag_OffBtagetaJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                          DoubleTag_OffBtagphiJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);

                          if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0)
                            {
                              DoubleTag_OffBtagptJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                              DoubleTag_OffBtagetaJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                              DoubleTag_OffBtagphiJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                          
                              if(deltaEta < deltaEtaCut[iJetCuts])
                                {
                                  DoubleTag_OffBtagptJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                                  DoubleTag_OffBtagetaJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                                  DoubleTag_OffBtagphiJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                                  DoubleTag_OffBtagptJetBt_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][theFlav[iJ]]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                                  DoubleTag_OffBtag_nJetConstituents_Bt_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][theFlav[iJ]]->fill(trgSelect,numberOfConstituentsInJet[leadingJets[iJ]],weight);
                                  if(iJ == 0)
                                    {
                                      DoubleTag_OffBtagMjjBt_deltaRcut_deltaEta[ibtag][iJetCuts]->fill(trgSelect,dijet_mass,weight);
                                    }
                                  if(sbtag[ibtag] == "CSVT" && iJetCuts == 0) {
                                    h_csv_discr->Fill(theBJetTag[ibtag][leadingJets[iJ]],weight);
                                  }
                                }
                            }
                        }
                      }
                    if (_doMC) {
                      double wbt[nSelJet];
                      double wbterror[nJetsCut];
              
                      for (int iJ=0; iJ<nSelJet; ++iJ) {
                        if(theFlav[iJ] != -1) {
                         //  wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                          //wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);

                          int local_iflav = theFlav[iJ];
                          if(theFlav[iJ] == 4)
                            local_iflav = 2;
                          if(theFlav[iJ] == 3)
                            local_iflav = 1;
                          

                          //wbt[iJ] = flavmix.eff(local_iflav,sbtag[ibtag].c_str(), true,true,false, iJ, pxJet[leadingJets[iJ]],pyJet[leadingJets[iJ]],pzJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                          wbt[iJ] =  bTagEffOffline->eff_merged(local_iflav,sbtag[ibtag].c_str(), sc[iJetCuts], ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                          
                          if(wbt[iJ] <= 0.0) {
                            throw std::exception();
                          }



                          wbterror[iJ] = bTagEffOffline->erreff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                        }
                      }
                      for(int istat = 0; istat < nstat; istat++)
                        {

                          if( isJetMatchedPartonExist[leadingJets[0]] && isJetMatchedPartonExist[leadingJets[1]] ) {

                          const double error0 = ( istat == 1 ? +1.0 : -1.0) * wbterror[0];
                          const double error1 = ( istat == 1 ? +1.0 : -1.0) * wbterror[1];
                  
                          double wtotDoubleTag = weight * 
                            (
                             wbt[0]*wbt[1]
                             //(wbt[0]+ (istat == 0 ? 0.0 : error0) ) * (wbt[1]+ (istat == 0 ? 0.0 : error1))
                             ); //double tag weight?
                          if(wtotDoubleTag <= 0.0)
                            {
                              std::cout << "error in double tag closure test. weight <=0: " << wtotDoubleTag << std::endl;
                              wtotDoubleTag = 1.0;
                              throw std::exception();
                            }
                          for (int iJ=0; iJ<nSelJet; ++iJ) {
                    
                            DoubleTag_OffBtagptJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                            DoubleTag_OffBtagetaJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                            DoubleTag_OffBtagphiJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);


                            if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0)
                              {
                                DoubleTag_OffBtagptJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag); 
                                DoubleTag_OffBtagetaJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag); 
                                DoubleTag_OffBtagphiJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag); 
                        
                                if(deltaEta<deltaEtaCut[iJetCuts])
                                  {
                                    DoubleTag_OffBtagptJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag); 
                                    DoubleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag); 
                                    DoubleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag); 
                                    
                                    DoubleTag_OffBtagptJetBtw_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat][theFlav[iJ]]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                                    DoubleTag_OffBtag_nJetConstituents_Btw_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat][theFlav[iJ]]->fill(trgSelect,numberOfConstituentsInJet[leadingJets[iJ]],wtotDoubleTag);
                                    if(iJ == 0 && istat == 0)
                                      {
                                        DoubleTag_OffBtagMjjBtw_deltaRcut_deltaEta[ibtag][iJetCuts]->fill(trgSelect,dijet_mass,wtotDoubleTag);
                                        
                                      }
                                   
                                  }

                              }
                          }
 
                          }
                        }
                    }
                  }
              }
          }



        if(nSelJet >= 3)
          {
            if((theFlav[0] != -1 && theFlav[1] != -1 && theFlav[2] != -1))
              {
                for (int ibtag=0; ibtag<nbtag; ++ibtag)
                  {
                    // require ONE btag for the three leading jets
                    if( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]))
                      {
                        for (int iJ=0; iJ<nSelJet; ++iJ) {
                          TripleTag_OffBtagptJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                          TripleTag_OffBtagetaJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                          TripleTag_OffBtagphiJetBt[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                      
                          if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0)
                            {
                              TripleTag_OffBtagptJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                              TripleTag_OffBtagetaJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                              TripleTag_OffBtagphiJetBt_deltaRcut[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                          
                              if(deltaEta<deltaEtaCut[iJetCuts])
                                {
                                  TripleTag_OffBtagptJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                                  TripleTag_OffBtagetaJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
                                  TripleTag_OffBtagphiJetBt_deltaRcut_deltaEta[ibtag][iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
                                  
                                  TripleTag_OffBtagptJetBt_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][theFlav[iJ]]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                                }
                            }
                        }
                      }
                    if (_doMC) {
                      double wbt[nSelJet];
                      double wbterror[nJetsCut];
              
                      for (int iJ=0; iJ<nSelJet; ++iJ) {
                        if(theFlav[iJ] != -1) {
                          wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                          wbterror[iJ] = bTagEffOffline->erreff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
                        }
                      }
                      for(int istat = 0; istat < nstat; istat++)
                        {
                          if( isJetMatchedPartonExist[leadingJets[0]] && isJetMatchedPartonExist[leadingJets[1]] && isJetMatchedPartonExist[leadingJets[2]]) {
                            const double error0 = ( istat == 1 ? +1.0 : -1.0) * wbterror[0];
                            const double error1 = ( istat == 1 ? +1.0 : -1.0) * wbterror[1];
                            const double error2 = ( istat == 1 ? +1.0 : -1.0) * wbterror[2];
                            double wtotDoubleTag = weight * 
                              (
                               (wbt[0] + (istat == 0 ? 0.0 : error0)) * (wbt[1] + (istat == 0 ? 0.0 : error1)) * (wbt[2] + (istat == 0 ? 0.0 : error2))
                               ); //tag weight?
                            if(wtotDoubleTag <= 0.0)
                              wtotDoubleTag = 1.0;
                            for (int iJ=0; iJ<nSelJet; ++iJ) {
                    
                              TripleTag_OffBtagptJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                              TripleTag_OffBtagetaJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                              TripleTag_OffBtagphiJetBtw[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);

                              if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0)
                                {
                                  TripleTag_OffBtagptJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                                  TripleTag_OffBtagetaJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                                  TripleTag_OffBtagphiJetBtw_deltaRcut[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);
                      
                                  if(deltaEta < deltaEtaCut[iJetCuts])
                                    {
                                      TripleTag_OffBtagptJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                                      TripleTag_OffBtagetaJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,etaJet[leadingJets[iJ]],wtotDoubleTag);
                                      TripleTag_OffBtagphiJetBtw_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][istat]->fill(trgSelect,phiJet[leadingJets[iJ]],wtotDoubleTag);

                                      if(istat == 0) {
                                        TripleTag_OffBtagptJetBtw_Fc_deltaRcut_deltaEta[ibtag][iJ][iJetCuts][theFlav[iJ]]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
                                      }

                                    }
                                }
                            }

                          }
                        }
                    }
                  }
              }
          }
      }



  } // loop over the tree
  // termination
 
  delete bTagEffOffline;
  

  hout->Write();
  hout->Close();
}
