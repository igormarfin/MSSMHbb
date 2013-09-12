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
#include <boost/filesystem.hpp>

//#include "Analysis/Utilities/interface/HBBTo4B.h"
//#include "Analysis/Utilities/src/HBBTo4B.cc"
#include "Analysis/Utilities/interface/svMass.h"
#include "Analysis/Utilities/interface/svMass2012.h"

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


#include "Analysis/HbbMSSMAnalysis/interface/SystReadInput.h"
#include "Analysis/Utilities/interface/jetFlavorCode.h"


TCanvas* canvas;

//... A snip of code has to be after all nonexecutable statements
#include "Analysis/HbbMSSMAnalysis/interface/HbbPuritySet.cc"
//..............................................................

//void TripleBtagAnalysis_SF() {
int main(int narg,char** varg) {


  // macro for analysis in channel with three btagged jets

  canvas = new TCanvas ("cg1","mycanvas",10,10,800,600);
  // open an ntuple file

  bool _doMC = false;  // set to true for MC, false for real data

  //  bool _doOnlineRelBtag= false; ///true; ///to do OnlineRelBtag
  bool _doOnlineRelBtag= true; ///to do OnlineRelBtag
  bool _doBkgPred= false; ///true; ///to do Bkg prediction

  bool isMCWithTrigSim = true;
  bool doLumiWeighting = true;  // weight MC events according to intgrated lumi
  bool forceTripleOnlineBtagMatch = false;//true;  // require triple online btag match on triple btag signature 

  const int nSelJet = 3;   // number of jets considered in b(bb) final state

  std::string onlineBtagVersion("V4"); // choices: smpf or V4

  bool doMVASel = false;
  const std::string MVAselFilename = "MVAselTree.root";

  if (boost::filesystem::exists(MVAselFilename)) {
    doMVASel = true;
  }
 
  std::cout << "MVA BG selection activation is " << doMVASel << std::endl;

 
  bool doTrigSelect = true;  // perform trigger selection according to scenario
  std::string Scenario("MediumMass2012");  // Higgs mass scenario // available: MediumMass2012, HighMass2012, VeryHighMass2012

  //  bool doBtagSF = true;   // apply btag scaling factors on MC btag efficiencies
  bool doBtagSF = false;   // apply btag scaling factors on MC btag efficiencies
  bool doJESVar = false;  // apply JES scaling variation (for evaluation of systematics)
  bool doBtagSFVar_bc = false; // apply BtagSF variation for bc (for evaluation of systematics)
  bool doBtagSFVar_q = false; // apply BtagSF variation for q (for evaluation of systematics)
  double JESVar = 0;
  double BtagSFVar = 0;
  
  //  bool doJERSF = true;  // apply jet energy resolution scaling factors 
  bool doJERSF = false;  // apply jet energy resolution scaling factors 
  bool doJERSFVar = false;  // apply jet energy resolution scaling variation (for evaluation of systematics)
  double JERVar = 2;

  // create systematics control object
  HbbSystControl theSystControl;  // constructor does nothing

  int systype=1;

  SystReadInput(systype,BtagSFVar);

  if (systype==2) doBtagSFVar_bc=true;
  if (systype==3) doBtagSFVar_q=true;

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
    std::cout << "Assume this is RealData" << std::endl;
    _doMC = false;
  } else {
    std::cout << "Assume this is MonteCarlo" << std::endl;
    _doMC = true;
  }  

  // create the trigger filter list. Only histograms for these triggers will be created
  std::vector<std::string> genericTriggerList;
  std::vector<std::string> HLTBtagTriggerObjectList;
  std::vector<std::string> triggerFilterList;
  std::vector<std::string> HLTBtagTriggerObjectFilterList;

  std::vector<std::string> L25BtagTriggerObjectFilterList;
  std::vector<std::string> L25BtagTriggerObjectList;

  if (_doMC) {
    std::cout << "warning: trigger selection for MC needs to be defined in the code! not yet done!"<< std::endl;
    //  if (isMCWithTrigSim) {
    //       //triggerFilterList.push_back("HLT_bbPhi_CentralJet46_CentralJet38_DiBTagIP3D_L25MinTag4_v1");
    //       //triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D"); 
    //       triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D");
    //       triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D");
    //     } else {
    //       triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D");
    //     }


    triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v");
    triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV_v");
    triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose_v");
    
    //    triggerFilterList.push_back("HLT_L1DoubleJet36Central_v");
    //    triggerFilterList.push_back("HLT_DiJet40Eta2p6_BTagIP3DFastPV_v");
    //    triggerFilterList.push_back("HLT_DiJet80Eta2p6_BTagIP3DFastPVLoose_v");



    //there are two HLT Btag trigger objects. The first two triggers specified above use the same HLT Btag objects
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");
    
    //for monitor trigger. first entry meaningless ...
    //     HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");
    //     HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BL1FastJetFastPV");
    //     HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BLooseL1FastJetFastPV");

    //Now the same for L25 objects (first element corresponds to first element of triggerFilterList and so on)
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhi1stTrackL1FastJetFastPV");

  } else {
    //  triggerFilterList.push_back("HLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D");
    //     triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D");
    //     triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D");
    //     triggerFilterList.push_back("HLT_CentralJet60_CentralJet53_DiBTagIP3D");
    triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v");
    triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV_v");
    triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose_v");
    
    //    triggerFilterList.push_back("HLT_L1DoubleJet36Central_v");
    //    triggerFilterList.push_back("HLT_DiJet40Eta2p6_BTagIP3DFastPV_v");
    //    triggerFilterList.push_back("HLT_DiJet80Eta2p6_BTagIP3DFastPVLoose_v");



    //there are two HLT Btag trigger objects. The first two triggers specified above use the same HLT Btag objects
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");
    
    //for monitor trigger. first entry meaningless ...
    //    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");
    //    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BL1FastJetFastPV");
    //    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhi1BLooseL1FastJetFastPV");

    //Now the same for L25 objects (first element corresponds to first element of triggerFilterList and so on)
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
    L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhi1stTrackL1FastJetFastPV");

  }

  // extract generic trigger list and number of input events from ntuple files
  float lumiScaleFac = 1;
  float xsect=0e0;
  int myentries=getHbbMetadata(theFileColl,genericTriggerList,HLTBtagTriggerObjectList,1000.,lumiScaleFac,_doMC,xsect);
  getHbbL25Metadata(theFileColl,L25BtagTriggerObjectList);

  std::cout<<"# of entries given by etHbbMetadata"<<myentries<<std::endl;  

  //  for(unsigned int i = 0; i < genericTriggerList.size(); i++)
  //     {
  //       std::cout << "generic trigger list " << genericTriggerList.at(i) << std::endl;
  //     }
  //   std::cout << "triggerfilterlist.size " << triggerFilterList.size() << std::endl;
  //   for(unsigned int i = 0; i <triggerFilterList.size(); i++)
  //     {
  //       std::cout << " triggerFilterList " << triggerFilterList.at(i) << std::endl;
  //     }

  

  // this trigger weight array will be used later in each event
  float* triggerWeight = new float[genericTriggerList.size()];
  // create trigger weight helper object
  HbbTrigWeight* theHbbTrigWeight = new HbbTrigWeight(&genericTriggerList,&triggerFilterList);

  f.AddFileInfoList((TCollection*) theFileColl.GetList());
  TTree* hbbtree = &f;


  // create trigger run selection object
  TriggerRunSelection trigSelector( Scenario, &genericTriggerList );

  TriggerRunSelection * ptrTrigSelector = & trigSelector;

  if (_doMC) {
    doTrigSelect = false;
    ptrTrigSelector=NULL;
  }

  hbbtree->SetBranchStatus("*L3*",1);

  checkHLTBtagMatch HLTBtagMatchObj(
                                    triggerFilterList,
                                    HLTBtagTriggerObjectFilterList,
                                    HLTBtagTriggerObjectList,
                                    doTrigSelect,
                                    ptrTrigSelector
                                    );


  //add the information for the L25 objects
  HLTBtagMatchObj.SetL25(L25BtagTriggerObjectFilterList,L25BtagTriggerObjectList);



  //HLTBtagMatchObj.check(,run)

  // try HBBTo4B class
  //   HBBTo4B * hbbTo4b = new HBBTo4B();
  //   hbbTo4b->Init(hbbtree);


#include "Analysis/HbbMSSMAnalysis/interface/HbbPurityInit.cc"

  // open seltree if needed

  TTree* seltree = NULL;
  TFile* seltreefile = NULL;
  if ( doMVASel ) {
    seltreefile = new TFile(TString(MVAselFilename.c_str()));
    if (seltreefile == NULL) {
      std::cout << "error: Opening of " << MVAselFilename.c_str() << " failed" << std::endl;
      return 0;
    }
    seltree = (TTree*) seltreefile->Get("selTree");
#include "Analysis/Utilities/interface/HbbSelect.cc"
    std::cout << "opening MVA selection tree" << std::endl;
  }
  


  //switch off TTree branches that are currently not needed?
  hbbtree->SetBranchStatus("*Muon*",0);
  hbbtree->SetBranchStatus("neutralHadronEnergyFraction",0);
  hbbtree->SetBranchStatus("photonEnergyFraction",0);
  hbbtree->SetBranchStatus("chargedHadronEnergyFraction",0);
  hbbtree->SetBranchStatus("electronEnergyFraction",0);
  hbbtree->SetBranchStatus("DeltaRWithHltBtag",0);
  hbbtree->SetBranchStatus("DeltaPtWithHltBtag",0);
  hbbtree->SetBranchStatus("*BJetTag*",0);
  hbbtree->SetBranchStatus("*CombSV*",1);
  hbbtree->SetBranchStatus("*L2*",0);
  hbbtree->SetBranchStatus("L2NumberOfJets",1);
  hbbtree->SetBranchStatus("L2JetPhi",1);
  hbbtree->SetBranchStatus("L3*",1);
  hbbtree->SetBranchStatus("*L3*",1);
  hbbtree->SetBranchStatus("*L25*",1);

  //end of switch off TTree branches that are currently not needed?


  //kinematic analysis cuts
  double deltaEtaCut12 = 1.7;//1.7;
  double maxEta = 1.7;
  double jetPtMin[nSelJet] = { 60.0, 53.0, 20};  // leading jet lower pt thresholds
  // double jetPtMin[nSelJet] = { 60.0, 60.0, 20};  // leading jet lower pt thresholds
  if (Scenario == "HighMass2012") {
    jetPtMin[0] = 80.0;
    jetPtMin[1] = 70.0;
    maxEta = 1.7;
  }
  if (Scenario == "VeryHighMass2012") {
    jetPtMin[0] = 160.0;
    jetPtMin[1] = 120.0;
    maxEta = 2.4;
    deltaEtaCut12 = 999999.9;//no delta eta cut
  }

  std::cout << "applied Jet pt cuts: " << jetPtMin[0] << "  " 
	    << jetPtMin[1] << "  " << jetPtMin[2] << std::endl;
  double jetPtMax[nSelJet] = {3500, 3500, 3500}; // leading jet upper pt thresholds
  
  // here we define dimensions and labels
  //was increased to 5! 'cc' has position '3' and 'bb' -- '4'
  const int nflav=5;  // number of flavor groups to be distinguished 
  //  const int nflav=3;  // number of flavor groups to be distinguished 
  const int nbtag = 4;  // number of offline btag working points
  const int ncateg = 3; // ranks of untagged jet
  const int ncorr = 2; // correction levels
  const int ntpat = 4; // trigger pattern (online btag). (0,1,2)=rank of non-btagged jet, 3=all combined

  //  const std::string flavlabel[nflav] = {"udsg","c","b"};  // labels for flavor groups
  const std::string flavlabel[] = {"udsg","c","b","cc","bb"};  // labels for flavor groups
  // define names of btag variants
  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };
  // corresponding discriminant cuts
  double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags
  //  bool bbPurityCorr[nbtag] = { true, true, true, true };  // specify if bbPurity collection should be applied
  bool bbPurityCorr[nbtag] = { false, false, false, false };  // specify if bbPurity collection should be applied

  // define flavor pair classes
  const int nfcDijet =  15;	///6;
  string sfcDijet[] = { "b1b1", "b1c1", "b1q", "c1c1", "c1q", "qq","c2b1","c2c1","c2q","c2c2","b2b1","b2c1","b2q","b2b2","b2c2" };
  const int nfc3rd = 3;
  string sfc3rd[nfc3rd] = { "q", "c", "b" };
  // flavor triplet classes
  const int nfcTrip = 6;
  string sfcTrip[nfcTrip] = { "bbb", "bbc", "bbq", "bcb", "bqb", "non-bb" };
  // the following are used for templates
  //  const int nfc = 3;  
  const int nfc = 5;  
  string sfc[] = { "q", "c", "b" ,"cc","bb"};

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

  ///2012 Medium
  // create the btag efficiency objects
  // offline btag object


  bTagEff* bTagEffOffline = NULL;
  bTagEff* bTagReleffOnline=NULL;


  // open the SVMass template files
  TString SVMassFileName("/data/user/rasp/Run/Hbb/svMass.root");
  TString SVMassOnlineFileName("/data/user/rasp/Run/Hbb/svMass1BL1FastJetFastPVHltMatch.root");
  





  ///2011 Medium

  //  bTagEff* bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt-N14/btagMatrixOffline-csv/plotmistag-b.root","offline",flavlabel,sbtag,nflav,nbtag);
  //  bTagEff* bTagReleffOnline  = new bTagEff("/afs/naf.desy.de/user/r/rmankel/public/btagMatrix/btagMatrixOnline-V4.root","online",flavlabel,sbtag,nflav,nbtag);

  ///  TString SVMassFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_svMass.root");
  //  TString SVMassOnlineFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_BTagHLTmatched_svMass.root");

  
  if (Scenario == "MediumMass2012") {
    bTagEffOffline= new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_mediummass.root","offline",flavlabel,sbtag,5,nbtag);
    // online btag object (with smp smoothing)
    bTagReleffOnline =  new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_normal_mediummass.root","online",flavlabel,sbtag,5,nbtag);
  } else if(Scenario == "VeryHighMass2012") {
    bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_veryhighmass.root","offline",flavlabel,sbtag,5,nbtag);
    bTagReleffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_loose_veryhighmass.root","online",flavlabel,sbtag,5,nbtag); 
    SVMassOnlineFileName="/data/user/rasp/Run/Hbb/svMass1BLooseL1FastJetFastPVHltMatch.root";
  } else if(Scenario == "HighMass2012") {
    bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_highmass.root","offline",flavlabel,sbtag,5,nbtag);
    bTagReleffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_normal_highmass.root","online",flavlabel,sbtag,5,nbtag); 
    SVMassOnlineFileName="/data/user/rasp/Run/Hbb/svMass1BL1FastJetFastPVHltMatch.root" ;
  } else {
    std::cout << "Scenario '" << Scenario << "' not known! exit." << std::endl;
    throw std::exception();
  }


  ///2012	
  svMass2012 * sv = new svMass2012(SVMassFileName);
  svMass2012 * svOnline = new svMass2012(SVMassOnlineFileName);

  ///2011

  //  svMass * sv = new svMass(SVMassFileName);
  //  svMass * svOnline = new svMass(SVMassOnlineFileName);


  std::cout << "Using online btag version: " << onlineBtagVersion << std::endl;

  // open the bbPurity correction functions
  TF1* fbbfrac[nbtag][ncateg];
  TFile* bbPur = new TFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/calib/bbPurity-csv-online.root");
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      fbbfrac[ibtag][icateg] = (TF1*) bbPur->Get( Form("fbbfrac-ww-%s-Cat%d",sbtag[ibtag].c_str(),icateg) );
      if ( fbbfrac[ibtag][icateg] == NULL ) {
	std::cout << "bbPur correction function not found for" << sbtag[ibtag].c_str() << " categ " << icateg << std::endl;
	if (bbPurityCorr[ibtag]) return 0;
      }
    }
  }



  // create the Root output file
  TFile* hout = new TFile("TripleBtagAnalysis.root","recreate");

  TH1::AddDirectory(true);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();  // not really needed (we inherit from TH1)

  Int_t nentries = (Int_t) hbbtree->GetEntries();
  std::cout << "Number of events in ntuple: " << nentries << std::endl;

#include "Analysis/Utilities/interface/HbbNtuple.cc"
#include "Analysis/HbbMSSMAnalysis/interface/HbbPurityBook.cc"

  // book the general histograms (before btag selection)
  hout->mkdir("general");
  hout->cd("general");
  std::cout << "book general histograms" << std::endl;
  const int ngeneralhistoscuts = 8; //number of pt cut scenarios used for the 'general' histograms
  const double generalHistosJetCuts[3*ngeneralhistoscuts] = {
    60.0, 53.0, 20.0,
    80.0, 70.0, 20.0,
    160.0, 120.0, 20.0,
    46.0, 38.0, 20.0, //this one is just for tests

    60.0, 53.0, 20.0,//with loose jet pu id
    80.0, 70.0, 20.0,//with loose jet pu id
    160.0, 120.0, 20.0,//with loose jet pu id
    46.0, 38.0, 20.0//with loose jet pu id

   
  };

  const bool applyLooseJetPUID[ngeneralhistoscuts] = {
    false,
    false,
    false,
    false,
    true,
    true,
    true,
    true
  };


  TrigHistArray* nJetA;
  TrigHistArray* nJetPostselA;
  TrigHistArray* ptJetA[3][ngeneralhistoscuts];
  TrigHistArray* etaJetA[3][ngeneralhistoscuts];
  TrigHistArray* phiJetA[3][ngeneralhistoscuts];
  TrigHistArray* phiL2JetA[3][ngeneralhistoscuts];
  TrigHistArray* nL2JetA[ngeneralhistoscuts];
  
  TrigHistArray* tcheJetA[3][ngeneralhistoscuts];
  TrigHistArray* tchpJetA[3][ngeneralhistoscuts];
  TrigHistArray* isJetWithBtagA[3][ngeneralhistoscuts];

  TrigHistArray* deltaEtaHLTBTagMatch2outof3[ngeneralhistoscuts]; //delta eta between two leading jets
  TrigHistArray* deltaEtaNoHLTBTagMatch[ngeneralhistoscuts]; //delta eta between two leading jets
  
  TrigHistArray* inclusivePtHLTBTagMatch2outof3[ngeneralhistoscuts]; //inclusive jet pt distribution
  TrigHistArray* inclusivePtNoHLTBTagMatch[ngeneralhistoscuts];//inclusive jet pt distribution

  TrigHistArray* CSV_discr_inclusive_HLTBTagMatch2outof3[ngeneralhistoscuts];
  TrigHistArray* CSV_discr_inclusive_NoHLTBTagMatch[ngeneralhistoscuts];

  TrigHistArray* CSV_discr_jet_HLTBTagMatch2outof3[3][ngeneralhistoscuts];
  TrigHistArray* CSV_discr_jet_NoHLTBTagMatch[3][ngeneralhistoscuts];

  TrigHistArray* CSV_discr_jet_HLTBTagMatch2outof3_isJetWithOnlineBtag[3][ngeneralhistoscuts];
  TrigHistArray* CSV_discr_jet_NoHLTBTagMatch_isJetWithOnlineBtag[3][ngeneralhistoscuts];

  nJetA = new TrigHistArray(&genericTriggerList,&triggerFilterList,"nJetA","nJets",41,-0.5,40.5);
  nJetPostselA = new TrigHistArray(&genericTriggerList,&triggerFilterList,"nJetPostselA","nJets Post Selection",41,-0.5,40.5);

  for(int iJetCuts = 0; iJetCuts < ngeneralhistoscuts; iJetCuts++)
    {
      deltaEtaHLTBTagMatch2outof3[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("deltaEtaHLTBTagMatch2outof3_JetCuts%d",iJetCuts),"#Delta#eta",60,-5.0,5.0);
      deltaEtaNoHLTBTagMatch[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("deltaEtaNoHLTBTagMatch_JetCuts%d",iJetCuts),"#Delta#eta",60,-5.0,5.0);

 
  
      inclusivePtHLTBTagMatch2outof3[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("inclusivePtHLTBTagMatch2outof3_JetCuts%d",iJetCuts),"inclusive jet p_{T} (GeV)",75,0.0,600.0);
      inclusivePtNoHLTBTagMatch[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("inclusivePtNoHLTBTagMatch_JetCuts%d",iJetCuts),"inclusive jet p_{T} (GeV)",75,0.0,600.0);
  
 
  
      CSV_discr_inclusive_HLTBTagMatch2outof3[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("CSV_discr_inclusive_HLTBTagMatch2outof3_JetCuts%d",iJetCuts),"CSV btag discriminator", 50, -0.05, 1.05);
      CSV_discr_inclusive_NoHLTBTagMatch[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("CSV_discr_inclusive_NoHLTBTagMatch_JetCuts%d",iJetCuts),"CSV btag discriminator", 50, -0.05, 1.05);

      //the same per jet
      for (int ii=0; ii<3; ++ii) {
        CSV_discr_jet_HLTBTagMatch2outof3[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("CSV_discr_jet%d_HLTBTagMatch2outof3_JetCuts%d",ii,iJetCuts),Form("CSV discr. of %dth leading jet",ii+1), 50, -0.05, 1.05);
        CSV_discr_jet_NoHLTBTagMatch[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("CSV_discr_jet%d_NoHLTBTagMatch_JetCuts%d",ii,iJetCuts),Form("CSV discr. of %dth leading jet",ii+1), 50, -0.05, 1.05);
 
        CSV_discr_jet_HLTBTagMatch2outof3_isJetWithOnlineBtag[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("CSV_discr_jet%d_HLTBTagMatch2outof3_isJetWithOnlineBtag_JetCuts%d",ii,iJetCuts),Form("CSV discr. of %dth leading matched jet",ii+1), 50, -0.05, 1.05);
        CSV_discr_jet_NoHLTBTagMatch_isJetWithOnlineBtag[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("CSV_discr_jet%d_NoHLTBTagMatch_isJetWithOnlineBtag_JetCuts%d",ii,iJetCuts),Form("CSV discr. of %dth leading matched jet",ii+1), 50, -0.05, 1.05);
      }

      nL2JetA[iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("nL2Jets_JetCuts%d",iJetCuts),"number of L2 jets", 100, 0.0, 100.0);
    
      for (int ii=0; ii<3; ++ii) {
        ptJetA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("ptj%d_JetCuts%d",ii,iJetCuts),Form("pt of %dth leading jet",ii+1),75,0,600);
        etaJetA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("aetaj%d_JetCuts%d",ii,iJetCuts),Form("eta of %dth leading jet",ii+1), 100, -3.5, 3.5);
        phiJetA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("aphij%d_JetCuts%d",ii,iJetCuts),Form("phi of %dth leading jet",ii+1), 100, -3.15, 3.15);
        phiL2JetA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("aphiL2j%d_JetCuts%d",ii,iJetCuts),Form("phi of %dth leading L2jet",ii+1), 100, -3.15, 3.15);


        tcheJetA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("tchej%d_JetCuts%d",ii,iJetCuts),Form("TCHE of %dth leading jet, ptjet>50",ii+1), 60, -10, 20 );
        tchpJetA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("tchpj%d_JetCuts%d",ii,iJetCuts),Form("TCHP of %dth leading jet, ptjet>50",ii+1), 60, -10, 20 );
        isJetWithBtagA[ii][iJetCuts] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("isWithBtagj%d_JetCuts%d",ii,iJetCuts),Form("IsJetWithBtag for %dth leading jet",ii+1),10,-3.5,6.5);
      }
    }


  TrigHistArray* matchPatternA = new TrigHistArray(&genericTriggerList,&triggerFilterList,"matchPattern","matchPattern",20,-0.5,19.5);  
  TrigHistArray* dPhiJet1Jet2A = new TrigHistArray(&genericTriggerList,&triggerFilterList,"adphij","delta-phi of two leading jets", 100, -3.15, 3.15); 

  // for monitoring the trigger weights
  TrigHistArray* triggerWeightA = new TrigHistArray(&genericTriggerList,&triggerFilterList,
						    "triggerWeight","triggerWeight",100, 0, 1);
  hout->cd();

  // book histogram arrays for dijet mass and pt
  //  Btcut    = btag selection by explicit cuts on discriminant
  //  Btweight = btag emulated by applying btag weights (MC only)
  //  Fc       = split into flavor pattern classes (MC only)
  hout->mkdir("mjjBtcut");
  hout->mkdir("mjjBtweight");
  hout->mkdir("mjjBtcutFc");
  hout->mkdir("mjjBtweightFc");
  
  hout->mkdir("DoubleOffBtagPtJetBtcut");
  hout->mkdir("DoubleOffBtagPtJetBtweight");

  hout->mkdir("DoubleOffBtagPtJetBtcutFc");
  hout->mkdir("DoubleOffBtagPtJetBtweightFc");


  hout->mkdir("ptjBtcut");
  hout->mkdir("ptjBtcutFc");
  hout->mkdir("ptjBtweight");
  hout->mkdir("ptjBtweightFc");

  hout->mkdir("dPhiBtweight");

  TrigHistArray* DoubleOffBtagptJetBtA[nbtag][nSelJet];
  TrigHistArray* DoubleOffBtagptJetBtwA[nbtag][nSelJet];
 
  TrigHistArray* DoubleOffBtagptJetBtA_Fc[nbtag][nfcTrip][nSelJet];
  TrigHistArray* DoubleOffBtagptJetBtwA_Fc[nbtag][nfcTrip][nSelJet];
 
 
  TrigHistArray* mDijetBtA[nbtag];
  TrigHistArray* ptJetBtA[nbtag][nSelJet];
  TrigHistArray* mDijetFcBtA[nbtag][nfcTrip];
  TrigHistArray* ptJetFcBtA[nbtag][nfcTrip][nSelJet];

  TrigHistArray* mDijetBtwA[nbtag];
  TrigHistArray* ptJetBtwA[nbtag][nSelJet];
  TrigHistArray* mDijetFcBtwA[nbtag][nfcTrip];
  TrigHistArray* ptJetFcBtwA[nbtag][nfcTrip][nSelJet];

  TrigHistArray* dPhiJet1Jet3BtwA[nbtag];
  TrigHistArray* dPhiJet2Jet3BtwA[nbtag];

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hout->cd("mjjBtcut");
    mDijetBtA[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                         Form("mjjbt%s",sbtag[ibtag].c_str()),
                                         Form("m(Jet1Jet2) 3*%s",sbtag[ibtag].c_str()),50,0,500 );
    hout->cd("mjjBtcutFc");
    for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
      mDijetFcBtA[ibtag][ifcTrip] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                      Form("mjjbt%sfc%s",sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),
                                                      Form("m(Jet1Jet2) 3*%s fc %s",sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),50,0,500 );
    }

    if (_doMC) {
      hout->cd("mjjBtweight");
      mDijetBtwA[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                            Form("mjjww%s",sbtag[ibtag].c_str()),
                                            Form("m(Jet1Jet2) weighted as 3*%s",sbtag[ibtag].c_str()),50,0,500 );
      hout->cd("mjjBtweightFc");
      for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
	mDijetFcBtwA[ibtag][ifcTrip] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                         Form("mjjww%sfc%s",sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),
                                                         Form("m(Jet1Jet2) weighted as 3*%s fc %s",sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),50,0,500 );
      }
    }
    hout->cd();

    hout->cd("dPhiBtweight");
    dPhiJet1Jet3BtwA[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
						Form("dPhiJ1J3%s",sbtag[ibtag].c_str()),Form("dPhiJ1J3 weighted as 3* %s",sbtag[ibtag].c_str()),50,-3.14,3.14);
    dPhiJet2Jet3BtwA[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
						Form("dPhiJ2J3%s",sbtag[ibtag].c_str()),Form("dPhiJ2J3 weighted as 3* %s",sbtag[ibtag].c_str()),50,-3.14,3.14);
    hout->cd();

    for (int iJ=0; iJ<nSelJet; ++iJ) {
      hout->cd("ptjBtcut");
      ptJetBtA[ibtag][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                              Form("ptj%dbt%s",iJ,sbtag[ibtag].c_str()),
                                              Form("pt leading jet %d 3*%s",iJ,sbtag[ibtag].c_str()),30,0,300);
      hout->cd("ptjBtcutFc");
      for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
	ptJetFcBtA[ibtag][ifcTrip][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                           Form("ptj%dbt%sfc%s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),
                                                           Form("pt leading jet %d 3*%s fc %s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),30,0,300);
      }
      if (_doMC) {
	hout->cd("ptjBtweight");
	ptJetBtwA[ibtag][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                 Form("ptj%dww%s",iJ,sbtag[ibtag].c_str()),
                                                 Form("pt leading jet %d weighted as 3*%s",iJ,sbtag[ibtag].c_str()),30,0,300);
	hout->cd("ptjBtweightFc");
	for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
	  ptJetFcBtwA[ibtag][ifcTrip][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                              Form("ptj%dww%sfc%s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),
                                                              Form("pt leading jet %d weighted as 3*%s fc %s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),30,0,300);
	}
      }
    }
    for (int iJ=0; iJ<nSelJet; ++iJ) {
      hout->cd("DoubleOffBtagPtJetBtcut");
      DoubleOffBtagptJetBtA[ibtag][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                           Form("DoubleOffBtag_ptj%dbt%s",iJ,sbtag[ibtag].c_str()),
                                                           Form("pt leading jet %d 3*%s",iJ,sbtag[ibtag].c_str()),30,0,300);

      for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
	DoubleOffBtagptJetBtA_Fc[ibtag][ifcTrip][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                         Form("DoubleOffBtag_ptj%dbt%s_Fc_%s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),
                                                                         Form("pt leading jet %d 3*%s fc %s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),30,0,300);
      }
      
      if (_doMC) {
	hout->cd("DoubleOffBtagPtJetBtweight");
	DoubleOffBtagptJetBtwA[ibtag][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                              Form("DoubleOffBtag_ptj%dww%s",iJ,sbtag[ibtag].c_str()),
                                                              Form("pt leading jet %d weighted as 3*%s",iJ,sbtag[ibtag].c_str()),30,0,300);
        for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
          DoubleOffBtagptJetBtwA_Fc[ibtag][ifcTrip][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                            Form("DoubleOffBtag_ptj%dww%s_Fc_%s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),
                                                                            Form("pt leading jet %d 3*%s fc %s",iJ,sbtag[ibtag].c_str(),sfcTrip[ifcTrip].c_str()),30,0,300);
        }
      }
    }
  }

  hout->cd();

  // svMass histograms
  hout->mkdir("svMass");
  hout->cd("svMass");
  TrigHistArray* svMassA[nSelJet];
  TrigHistArray* svMassBtA[nbtag][nSelJet];
  for (int iJ=0; iJ<nSelJet; ++iJ) {
    svMassA[iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
				    Form("svMassj%d",iJ),Form("svMass jet %d",iJ),20,0.,10.);
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      svMassBtA[ibtag][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
					       Form("svMassj%d %s",iJ,sbtag[ibtag].c_str()),
					       Form("svMass jet %d %s",iJ,sbtag[ibtag].c_str()),20,0.,10.);
    }
  }
  hout->cd();

  // evBtag and dijetMassEvBTag histograms
  hout->mkdir("evBtag");
  hout->mkdir("massEvBtag");
  TrigHistArray* evBtagBtA[nbtag];
  TrigHistArray2D* massEvBtagBtA[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hout->cd("evBtag");
    evBtagBtA[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
					 Form("evBTag_%s",sbtag[ibtag].c_str()),
					 Form("evBTag %s",sbtag[ibtag].c_str()),6,-0.5,5.5);
    hout->cd("massEvBtag");
    massEvBtagBtA[ibtag] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
					       Form("mjjEvBTag_%s",sbtag[ibtag].c_str()),
					       Form("mjjEvBTag %s",sbtag[ibtag].c_str()),25,0.,500.,6,-0.5,5.5);
  }
  hout->cd();

  // dijetMassEvBTag trigger-weighted histograms (MC only)
  hout->mkdir("massEvBtagTW");
  TrigHistArray2D* massEvBtagBtTWA[nbtag];
  if (_doMC) {
    hout->cd("massEvBtagTW");
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      massEvBtagBtTWA[ibtag] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
						   Form("mjjEvBTagTW_%s",sbtag[ibtag].c_str()),
						   Form("mjjEvBTag TW %s",sbtag[ibtag].c_str()),25,0.,500.,6,-0.5,5.5);
    }
  }
  hout->cd();

  // create the counter histograms for flavor pair classes (purity statistics). Meaningful for MC only
  hout->mkdir("fcCount");
  hout->cd("fcCount");
  TH1F* hfc[nbtag];
  TH1F* hfcm[nbtag];
  TH1F* hfcww[nbtag];
  TH1F* hfcmww[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hfc[ibtag] = new TH1F(Form("hfc%s",sbtag[ibtag].c_str()),Form("Flavor triple code 3* %s",sbtag[ibtag].c_str()),46,-1.5,44.5);
    hfcm[ibtag] = new TH1F(Form("hfcm%s",sbtag[ibtag].c_str()),Form("Flavor triple code 3* %s, mass 100-140",sbtag[ibtag].c_str()),46,-1.5,44.5);
    hfcww[ibtag] = new TH1F(Form("hfc%sww",sbtag[ibtag].c_str()),Form("Flavor triple code weighted as 3* %s",sbtag[ibtag].c_str()),46,-1.5,44.5);
    hfcmww[ibtag] = new TH1F(Form("hfcm%sww",sbtag[ibtag].c_str()),Form("Flavor triple code weighted as 3* %s, mass 100-140",sbtag[ibtag].c_str()),46,-1.5,44.5);
    for (int ifcDijet=0; ifcDijet<nfcDijet; ++ifcDijet) {
      for (int ifc3rd=0; ifc3rd<nfc3rd; ++ifc3rd) {
	hfc[ibtag]->GetXaxis()->SetBinLabel(nfcDijet* ifc3rd + ifcDijet +2,Form("%s%s",sfcDijet[ifcDijet].c_str(),sfc3rd[ifc3rd].c_str()));
	hfcm[ibtag]->GetXaxis()->SetBinLabel(nfcDijet* ifc3rd + ifcDijet +2,Form("%s%s",sfcDijet[ifcDijet].c_str(),sfc3rd[ifc3rd].c_str()));
	hfcww[ibtag]->GetXaxis()->SetBinLabel(nfcDijet* ifc3rd + ifcDijet +2,Form("%s%s",sfcDijet[ifcDijet].c_str(),sfc3rd[ifc3rd].c_str()));
	hfcmww[ibtag]->GetXaxis()->SetBinLabel(nfcDijet* ifc3rd + ifcDijet +2,Form("%s%s",sfcDijet[ifcDijet].c_str(),sfc3rd[ifc3rd].c_str()));
      }
    }
  }
  TH1F* hfctrip = new TH1F("hfctrip","Flavor condensed triple code",7,-1.5,6.5);
  for (int ifcTrip=0; ifcTrip<nfcTrip; ++ifcTrip) {
    hfctrip->GetXaxis()->SetBinLabel(ifcTrip+2,sfcTrip[ifcTrip].c_str());
  }
  hout->cd();

  // Mass templates derived from double btag sample
  hout->mkdir("mjjTemplates");
  hout->cd("mjjTemplates");
  TrigHistArray* massTemplateA[nfc][nbtag][ncateg][ncorr][ntpat];
  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	for (int icorr=0; icorr<ncorr; ++icorr) {
	  for (int itpat=0; itpat<ntpat; ++itpat) {
	    massTemplateA[ifc][ibtag][icateg][icorr][itpat]
	      = new TrigHistArray(&genericTriggerList,&triggerFilterList,
				  Form("mjjTemp_%s_%s_Cat%dCorr%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,icorr,itpat),
				  Form("mJet1Jet2 Template, %s %s Cat%d Corr%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,icorr,itpat),
				  50,0,500);
	  }
	}
      }
    }
  }
  hout->cd();
  
  TrigHistArray* bTagTemplateA[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* massBTagTemplateA[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* errorMassBTagTemplateA[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* massBTagTemplateUncldA[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* massBTagTemplateCldA[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* massBTagTemplateNonbbRA[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* massBTagTemplateCldRA[nfc][nbtag][ncateg][ntpat];
  // event bTag Templates
  hout->mkdir("bTagTemplates");
  hout->cd();
  hout->mkdir("massBTagTemplates");
  hout->cd();
  // error for mass btag templates
  hout->mkdir("errorMassBTagTemplates");
  hout->cd();
  hout->mkdir("massBTagTemplatesCld");
  hout->cd();
  hout->mkdir("massBTagTemplatesNonbbR");
  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	for (int itpat=0; itpat<ntpat; ++itpat) {
	  hout->cd("bTagTemplates");
	  bTagTemplateA[ifc][ibtag][icateg][itpat] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                                       Form("bTagTemplate_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                       Form("bTagTemplate, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                       6,-0.5,5.5);
	  hout->cd("massBTagTemplates");
	  massBTagTemplateA[ifc][ibtag][icateg][itpat] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
                                                                             Form("MassBTagTemplate_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                             Form("MassBTagTemplate, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
									     25,0.,500.,6,-0.5,5.5);
	  hout->cd("errorMassBTagTemplates");
	  errorMassBTagTemplateA[ifc][ibtag][icateg][itpat] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
                                                                                  Form("ErrorMassBTagTemplate_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                  Form("ErrorMassBTagTemplate, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                  25,0.,500.,6,-0.5,5.5);
	  hout->cd("massBTagTemplatesCld");  // templates cleaned with variable threshold
	  massBTagTemplateCldA[ifc][ibtag][icateg][itpat] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
                                                                                Form("MassBTagTemplateCld_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                Form("MassBTagTemplate cleaned, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                25,0.,500.,6,-0.5,5.5);
	  massBTagTemplateUncldA[ifc][ibtag][icateg][itpat] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,  // for comparison, uncleaned
                                                                                  Form("MassBTagTemplateUncld_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                  Form("MassBTagTemplate non-cleaned, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                  25,0.,500.,6,-0.5,5.5);
	  hout->cd("massBTagTemplatesNonbbR");  // templates cleaned with fixed threshold method (R method)
	  massBTagTemplateNonbbRA[ifc][ibtag][icateg][itpat] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
                                                                                   Form("MassBTagTemplateNonbbR_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                   Form("MassBTagTemplate nonbb R, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                   25,0.,500.,6,-0.5,5.5);	
	  massBTagTemplateCldRA[ifc][ibtag][icateg][itpat] = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
                                                                                 Form("MassBTagTemplateCldR_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                 Form("MassBTagTemplate cleaned R, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
                                                                                 25,0.,500.,6,-0.5,5.5);
        }
      }
    }
  }
  hout->cd();

  // statistics for online btag pattern (double btag)
  hout->mkdir("onlineBtag");
  hout->cd("onlineBtag");
  TrigHistArray* tpatA[nbtag][ncateg];
  TrigHistArray* tpatAllA[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      tpatA[ibtag][icateg] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                               Form("tpat_%s_Cat%d",sbtag[ibtag].c_str(),icateg),
                                               Form("Online btag trig pattern, 2*%s, Cat%d",sbtag[ibtag].c_str(),icateg),9,-0.5,8.5);
    }
    tpatAllA[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                        Form("tpat_%s_all",sbtag[ibtag].c_str()),
                                        Form("Online btag trig pattern, 2*%s, all",sbtag[ibtag].c_str()),9,-0.5,8.5);
  }

  // statistics for online btag pattern (triple btag)
  TrigHistArray* atpattripall[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    atpattripall[ibtag] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                            Form("tpattrip_%s_all",sbtag[ibtag].c_str()),
                                            Form("Online btag trig pattern, 3*%s, all",sbtag[ibtag].c_str()),9,-0.5,8.5);
  }
  hout->cd();

  // double btag kinematic histograms
  hout->mkdir("doubleBtag");
  hout->cd("doubleBtag");
  TrigHistArray* mDibBtcutA[nbtag][ncateg];
  TrigHistArray* mDibBtweightA[nbtag][ncateg];
  TrigHistArray* ptDibBtcutA[nbtag][ncateg][nSelJet];
  TrigHistArray* ptDibBtweightA[nbtag][ncateg][nSelJet];

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      mDibBtcutA[ibtag][icateg] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
						    Form("mdib_%s_Cat%d_bt",sbtag[ibtag].c_str(),icateg),
						    Form("Jet1Jet2 2* %s Cat%d",sbtag[ibtag].c_str(),icateg),
						    25,0,500);
      mDibBtweightA[ibtag][icateg] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
						       Form("mdib_%s_Cat%d_ww",sbtag[ibtag].c_str(),icateg),
						       Form("Jet1Jet2 weighted as 2* %s Cat%d",sbtag[ibtag].c_str(),icateg),
						       25,0,500);
      for (int iJ=0; iJ<nSelJet; ++iJ) {
	ptDibBtcutA[ibtag][icateg][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
							   Form("ptdib_%s_Cat%d_j%d_bt",sbtag[ibtag].c_str(),icateg,iJ),
							   Form("pt dibtag %s Cat%d jet %d",sbtag[ibtag].c_str(),icateg,iJ),
							   30,0,300);
	ptDibBtweightA[ibtag][icateg][iJ] = new TrigHistArray(&genericTriggerList,&triggerFilterList,
                                                              Form("ptdib_%s_Cat%d_j%d_ww",sbtag[ibtag].c_str(),icateg,iJ),
                                                              Form("pt dibtagweight %s Cat%d jet %d",sbtag[ibtag].c_str(),icateg,iJ),
                                                              30,0,300);
      }
    }
  }


  hout->cd();

  // double btag purity histograms. These are produced to create bbPurity corrections.
  //  Btcut    = btag selection by explicit cuts on discriminant
  //  Btweight = btag emulated by applying btag weights (MC only)
  hout->mkdir("bbPurityBtcut");
  hout->mkdir("bbPurityBtweight");
  TH1F* mDibBtcutFcH[nfcDijet][nbtag][ncateg];
  TH1F* mDibBtweightFcH[nfcDijet][nbtag][ncateg];
  if (_doMC) {
    for (int ifcDijet=0; ifcDijet<nfcDijet; ++ifcDijet) {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	for (int icateg=0; icateg<ncateg; ++icateg) {
	  hout->cd("bbPurityBtcut");
	  mDibBtcutFcH[ifcDijet][ibtag][icateg] = new TH1F(Form("mdib_%s_%s_Cat%dbt",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
                                                           Form("Jet1Jet2 %s 2* %s Cat%d",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
                                                           50,0,500);
	  hout->cd("bbPurityBtweight");
	  mDibBtweightFcH[ifcDijet][ibtag][icateg] = new TH1F(Form("mdib_%s_%s_Cat%dww",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
                                                              Form("Jet1Jet2 %s weighted as 2* %s Cat%d",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
                                                              50,0,500);
	}
      }
    }
  } /// _doMC


  hout->cd();

  //cout<<"Book igor's histo"<<endl;


  hout->mkdir("diJetMassFlavor");
  TH1F* mDiJetMassFlavor[nfcDijet][nbtag];
  if (_doMC) {

    for (int ifcDijet=0; ifcDijet<nfcDijet; ++ifcDijet) {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {

        hout->cd("diJetMassFlavor");
        mDiJetMassFlavor[ifcDijet][ibtag]=new TH1F(Form("mDiJetMassFlavor_%s_%s",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str()),
                                                   Form("Jet1Jet2 %s 2* %s ",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str()), 50,0,500);

      }
    }

  } ///_doMC 

  hout->cd(); 

  // for the double btag signature, histogram the pt of the btagged jets
  const int nAB = 2;
  string sAB[nAB] = {"jetA","jetB"};
  TH1F* hptdibbt[nfc][nAB][nbtag][ncateg];
  hout->cd("bbPurityBtcut");
  if (_doMC) {
    for (int ifc=0; ifc<nfc; ++ifc) {
      for (int iAB=0; iAB<nAB; ++iAB) {
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	  for (int icateg=0; icateg<ncateg; ++icateg) {
	    hptdibbt[ifc][iAB][ibtag][icateg] = new TH1F(Form("ptdibbt_%s_%s_%s_Cat%d",sfc[ifc].c_str(),sAB[iAB].c_str(),
                                                              sbtag[ibtag].c_str(),icateg),
                                                         Form("ptdibbt_%s_%s_%s_Cat%d",sfc[ifc].c_str(),sAB[iAB].c_str(),
                                                              sbtag[ibtag].c_str(),icateg),30,0,300);
	  }
	}
      }
    }
  }
  hout->cd();

  // background prediction histograms
  hout->mkdir("bgPredict");
  hout->cd("bgPredict");
  TrigHistArray* massPred[nfc][nbtag][ncateg][ntpat];
  TrigHistArray2D* massBTagPred[nfc][nbtag][ncateg][ntpat];
  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	for (int itpat=0; itpat<ntpat; ++itpat) {
	  massPred[ifc][ibtag][icateg][itpat]
	    = new TrigHistArray(&genericTriggerList,&triggerFilterList,
				Form("massPred_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
				Form("mJet1Jet2 Prediction, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
				50,0,500);
	  massBTagPred[ifc][ibtag][icateg][itpat]
	    = new TrigHistArray2D(&genericTriggerList,&triggerFilterList,
				  Form("MassBTagPred_%s_%s_Cat%dTpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
				  Form("MassBTag prediction, %s %s Cat%d Tpat%d",sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,itpat),
				  25,0.,500.,6,-0.5,5.5);
	}
      }
    }
  }
  hout->cd();
  hout->mkdir("debug_L3");
  hout->cd("debug_L3");
  TH1D* hptL3Objects_btagtrigger = new TH1D("hptL3Objects_btagtrigger","hptL3Objects_btagtrigger", 300, 0.0, 300.0);
  TH1D* hptL3Objects_matchedtoL2_btagtrigger = new TH1D("hptL3Objects_matchedtoL2_btagtrigger","hptL3Objects_matchedtoL2_btagtrigger", 300, 0.0, 300.0);
  TH1D* h_distance_L3Objects_vs_L2Objects_btagtrigger = new TH1D("h_distance_L3Objects_vs_L2Objects_btagtrigger","h_distance_L3Objects_vs_L2Objects_btagtrigger", 80.0, 0.0, 2.0);
  
  hout->cd();


  //check MVA tree
  if ( doMVASel ) {
    if(seltree->GetEntries() != nentries) {
      std::cout << "MVASel: number of entries mismatch: " << seltree->GetEntries() << " " << nentries << std::endl;
      return 0;
    }
  }

  std::cout<<"before loop, the # of entries is "<<nentries<<std::endl;
  // loop over the tree
  for (Int_t iE=0; iE<nentries; iE++) {
    //    if (iE>10000) {
    //      std::cout << "Event limit reached !" << std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
    //      break;
    //    }
    hbbtree->GetEntry(iE);
    //hbbTo4b->GetEntry(iE);

    if(iE % 200000 == 0)
      {
        const double percent = 100.0*(double)iE/((double)nentries);
        std::cout << "process event " << iE << " - " << percent << " % done." << std::endl;
      }

    //	std::cout << "process event " << iE << std::endl;


    if ( doMVASel ) {
      seltree->GetEntry(iE);
      if ( (selRun != run) || (selEvent != event)) {
	std::cout << " MVASel: synchronization mismatch "
		  << " run/event= " << run << " / " << event
		  << " selRun/selEvent= " << selRun << " / " << selEvent
		  << std::endl;
	return 0;
      }
      if (! (selPass)) continue;
    }
   
 



    nJetA->fill(trgAccept,float(numberOfJets),1);

    // perform JES variation if switched on
    double theBtagSF[nbtag] = {1.0, 1.0, 1.0, 1.0};
    if (doJESVar) theSystControl(0, JESVar, -1, -1, 0., 0., theBtagSF[0]);  // no SF returned
    
    // perform JER scaling and variation if switched on
    double theJERSF = 1.0;
    if ( doJERSF && ! doJERSFVar )    theSystControl(4, 0, -1, -1, 0., 0., theJERSF);  // no SF returned
    if ( doJERSFVar ) theSystControl(4, JERVar, -1, -1, 0., 0., theJERSF);  // no SF returned




    // trigger selection
    unsigned trgSelect = trgAccept;
    if (doTrigSelect) trgSelect = trgAccept & trigSelector.mask( run );


    // compute the weight for this event
    float weight = 1;
    if ( _doMC && doLumiWeighting )  weight = lumiScaleFac;






    //check the distributions of L3 objects
    for (int iJet=0; iJet < l3NumberOfBJets; ++iJet) {

      std::string theTrigNameBtag;
      theTrigNameBtag.assign("HLT_DiJet40Eta2p6_BTagIP3DFastPV_v");


      std::vector<std::string>::iterator tSlotBtag = std::find(genericTriggerList.begin(), genericTriggerList.end(), theTrigNameBtag);
      if (tSlotBtag != genericTriggerList.end()) {
        //std::cout << "Btag trigger found at slot " << tSlotBtag - genericTriggerList.begin() << std::endl;
      } else {
        std::cout << "Btag trigger not found in any slot" << std::endl;
        return 0;
      }
      unsigned int tNumberBtag = tSlotBtag - genericTriggerList.begin();


      //ask for dijet monitoring trigger
      if(trgAccept & (1<<tNumberBtag)) {

        hptL3Objects_btagtrigger->Fill(l3PtBJet[iJet],weight);

         
        double mindeltaR = 1000.0;

        for (int iJetL2=0; iJetL2 < l2NumberOfJets; ++iJetL2) {
            
          float dphij =  l2PhiJet[iJetL2] - l3PhiBJet[iJet];
          if (dphij>3.1415926) dphij -= 2*3.1415926;
          if (dphij<-3.1415926) dphij += 2*3.1415926;
            
          float detaj = l2EtaJet[iJetL2] - l3EtaBJet[iJet];
          double deltaR = sqrt( dphij*dphij + detaj*detaj );
          
          if(deltaR < mindeltaR) {
            mindeltaR = deltaR;
          }
        }
        if(mindeltaR < 0.5)
          hptL3Objects_matchedtoL2_btagtrigger->Fill(l3PtBJet[iJet],weight);
        h_distance_L3Objects_vs_L2Objects_btagtrigger->Fill(mindeltaR,weight);
      }    
    }
       
    //end of check L3 object distributions


    //general jet histograms
    // fill general jet histograms 
    // do jet selection with various pt cut scenarios according to trigger


    for(int iJetCuts = 0; iJetCuts < ngeneralhistoscuts; iJetCuts++)
      {
        std::vector<int> leadingJets;

        // find set of leading jets
        int nJet = 0;
        // loop over the jets
        for (int iJet=0; iJet<numberOfJets; ++iJet) {
          if (nJet >= nSelJet) break;
          if (! (fabs(etaJet[iJet])<maxEta) ) continue;

          if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
          if(applyLooseJetPUID[iJetCuts] && !puJetIDLoose[iJet]) continue;

          if ( (ptJet[iJet] > generalHistosJetCuts[3*iJetCuts+nJet]) && (ptJet[iJet] < jetPtMax[nJet]) ) {
            leadingJets.push_back(iJet);
            ++nJet;
          }
        }



        if (nJet < nSelJet) continue;


        float deltaRj = -1.0f;

        float deltaR12 = -1.0f;
        float deltaR23 = -1.0f;
        float deltaR13 = -1.0f;

        float deltaEta = 9999.9; //between leading jets

        if(nJet >= 2)
          {
            const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[1]];
            float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
            if (dphij>3.1415926) dphij -= 2*3.1415926;
            if (dphij<-3.1415926) dphij += 2*3.1415926;      
            deltaR12 = sqrt( dphij*dphij + detaj*detaj );
            deltaEta = TMath::Abs(etaJet[leadingJets[1]] - etaJet[leadingJets[0]]);
          }
        if(nJet >= 3)
          {
            const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[2]];
            float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[0]];
            if (dphij>3.1415926) dphij -= 2*3.1415926;
            if (dphij<-3.1415926) dphij += 2*3.1415926;      
            deltaR13 = sqrt( dphij*dphij + detaj*detaj );
          }
        if(nJet >= 3)
          {
            const float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[2]];
            float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[1]];
            if (dphij>3.1415926) dphij -= 2*3.1415926;
            if (dphij<-3.1415926) dphij += 2*3.1415926;      
            deltaR23 = sqrt( dphij*dphij + detaj*detaj );
          }
                



        // check the matching flags
        int mJet = 0;
        int matchPat = 0;
        for (int iJet=0; iJet<numberOfJets; ++iJet) {
          if (! (ptJet[iJet] > 15) ) continue;
          if (! (fabs(etaJet[iJet])<maxEta) ) continue;
          if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
          {


            if (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[iJet],isJetWithL25JetBitPattern[iJet],run)) {
	
              matchPat = matchPat | (1<<mJet);


            }

          }
          ++mJet;
        }
  

        bool tripleOnlineBtagMatchOK = true;
        if (forceTripleOnlineBtagMatch) {
          tripleOnlineBtagMatchOK = false;
          if (nJet>=nSelJet ) {
            tripleOnlineBtagMatchOK = (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run))
              || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run))
              || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run));
          }
        }
        
        bool BtagMatch2outof3 = (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run))
          || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run))
          || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run));
    

        if ( (nJet>=nSelJet) && BtagMatch2outof3 && deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12)
          {
            const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[1]];
            deltaEtaHLTBTagMatch2outof3[iJetCuts]->fill(trgSelect,detaj,weight);

            for (int iJ=0; iJ<nSelJet; ++iJ) {
              inclusivePtHLTBTagMatch2outof3[iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
              CSV_discr_inclusive_HLTBTagMatch2outof3[iJetCuts]->fill(trgSelect,combSVBJetTag[leadingJets[iJ]],weight);
              CSV_discr_jet_HLTBTagMatch2outof3[iJ][iJetCuts]->fill(trgSelect,combSVBJetTag[leadingJets[iJ]],weight);
              if(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[iJ]],isJetWithL25JetBitPattern[leadingJets[iJ]],run))
                {
                  CSV_discr_jet_HLTBTagMatch2outof3_isJetWithOnlineBtag[iJ][iJetCuts]->fill(trgSelect,combSVBJetTag[leadingJets[iJ]],weight);
                }
            }
          }

    
        if ( (nJet>=nSelJet) && deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) //no HLT btag matching
          {
            const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[1]];
            deltaEtaNoHLTBTagMatch[iJetCuts]->fill(trgSelect,detaj,weight);

            for (int iJ=0; iJ<nSelJet; ++iJ) {
              inclusivePtNoHLTBTagMatch[iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
              CSV_discr_inclusive_NoHLTBTagMatch[iJetCuts]->fill(trgSelect,combSVBJetTag[leadingJets[iJ]],weight);
              CSV_discr_jet_NoHLTBTagMatch[iJ][iJetCuts]->fill(trgSelect,combSVBJetTag[leadingJets[iJ]],weight);
              if(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[iJ]],isJetWithL25JetBitPattern[leadingJets[iJ]],run))
                {
                  CSV_discr_jet_NoHLTBTagMatch_isJetWithOnlineBtag[iJ][iJetCuts]->fill(trgSelect,combSVBJetTag[leadingJets[iJ]],weight);
                }
            }

          }

      
        if ( (nJet>=nSelJet) && deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12 ) {
          nL2JetA[iJetCuts]->fill(trgSelect,l2NumberOfJets,weight);
          for (int iJ=0; iJ<l2NumberOfJets; ++iJ) {
            if(iJ<3)
              phiL2JetA[iJ][iJetCuts]->fill(trgSelect,l2PhiJet[iJ],weight);
          }
        }
        if ( (nJet>=nSelJet) && tripleOnlineBtagMatchOK && deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) {
          for (int iJ=0; iJ<nSelJet; ++iJ) {
            // jet kinematics
            ptJetA[iJ][iJetCuts]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
            etaJetA[iJ][iJetCuts]->fill(trgSelect,etaJet[leadingJets[iJ]],weight);
            phiJetA[iJ][iJetCuts]->fill(trgSelect,phiJet[leadingJets[iJ]],weight);
            isJetWithBtagA[iJ][iJetCuts]->fill(trgSelect,HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[iJ]],isJetWithL25JetBitPattern[leadingJets[iJ]],run),weight);
          }
        }
      }
    //end of general jet histograms

    float deltaR12 = -1.0f;
    float deltaR23 = -1.0f;
    float deltaR13 = -1.0f;

    float deltaEta = 9999.0; //between leading jets
    
    std::vector<int> leadingJets;

    // find set of leading jets
    int nJet = 0;
    // loop over the jets
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      if (nJet >= nSelJet) break;
      if (! (fabs(etaJet[iJet])<maxEta) ) continue;

      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
      if(!puJetIDLoose[iJet]) continue;
      if(!jetIDLoose[iJet] ) continue;
      if ( (ptJet[iJet] > jetPtMin[nJet]) && (ptJet[iJet] < jetPtMax[nJet]) ) {
	leadingJets.push_back(iJet);
	++nJet;
      }
    }

    if (nJet < nSelJet) continue;




    ///Added,  to select using new steps deltaR13, deltaR12, deltaR23 and deltaEtaCut12

    float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[1]];
    float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
    if (dphij>3.1415926) dphij -= 2*3.1415926;
    if (dphij<-3.1415926) dphij += 2*3.1415926;      
    deltaR12 = sqrt( dphij*dphij + detaj*detaj );
    deltaEta = TMath::Abs(etaJet[leadingJets[1]] - etaJet[leadingJets[0]]);
        
    detaj=etaJet[leadingJets[0]] - etaJet[leadingJets[2]];
    dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[0]];
    if (dphij>3.1415926) dphij -= 2*3.1415926;
    if (dphij<-3.1415926) dphij += 2*3.1415926;      
    deltaR13 = sqrt( dphij*dphij + detaj*detaj );

    detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[2]];
    dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[1]];
    if (dphij>3.1415926) dphij -= 2*3.1415926;
    if (dphij<-3.1415926) dphij += 2*3.1415926;      
    deltaR23 = sqrt( dphij*dphij + detaj*detaj );
        

    if ( !(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) ) continue;
                

#include "Analysis/HbbMSSMAnalysis/interface/HbbPurityFill.cc"

   
    // check the matching flags
    int mJet = 0;
    int matchPat = 0;
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      if (! (ptJet[iJet] > 15) ) continue;
      if (! (fabs(etaJet[iJet])<maxEta) ) continue;
      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
      if (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[iJet],isJetWithL25JetBitPattern[iJet],run)) {
	matchPat = matchPat | (1<<mJet);
      }
      ++mJet;
    }
    matchPatternA->fill(trgSelect,matchPat);



    bool tripleOnlineBtagMatchOK = true;
    if (forceTripleOnlineBtagMatch) {
      tripleOnlineBtagMatchOK = false;
      if (nJet>=nSelJet ) {
	tripleOnlineBtagMatchOK = (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run))
	  || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run))
	  || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run));
      }
    }

    
    if ( (nJet>=nSelJet) && tripleOnlineBtagMatchOK ) {
      float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
      if (dphij>3.1415926) dphij -= 2*3.1415926;
      if (dphij<-3.1415926) dphij += 2*3.1415926;
      float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[0]];
      deltaR12 = sqrt( dphij*dphij + detaj*detaj );
      dPhiJet1Jet2A->fill(trgSelect,dphij,weight);
      deltaEta = TMath::Abs(etaJet[leadingJets[1]] - etaJet[leadingJets[0]]);
    }

    if(nJet >= 3)
      {
        const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[2]];
        float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[0]];
        if (dphij>3.1415926) dphij -= 2*3.1415926;
        if (dphij<-3.1415926) dphij += 2*3.1415926;      
        deltaR13 = sqrt( dphij*dphij + detaj*detaj );
      }
    if(nJet >= 3)
      {
        const float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[2]];
        float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[1]];
        if (dphij>3.1415926) dphij -= 2*3.1415926;
        if (dphij<-3.1415926) dphij += 2*3.1415926;      
        deltaR23 = sqrt( dphij*dphij + detaj*detaj );
      }


    // determine jet flavor from MC info
    int theFlav[nSelJet];
    for (int iJ=0; iJ<nSelJet; ++iJ) {
      theFlav[iJ] = jetFlavorCodeNewCat(leadingJets[iJ]);
    }

    // determine the dijet flavor code for the first two leading jets
    //    cout<<"Before first diJetCode"<<endl;
    int theFcDijet = diJetFlavorCode( leadingJets[0], leadingJets[1] );

    // now combine with flavor of third jet
    int theFc3rd = theFlav[2];
    int theFc = nfcDijet * theFc3rd + theFcDijet;

    // now make the triplet code
    int theFcTrip = -1;
    switch (theFcDijet) {
    case 0:  // (bb)
      switch (theFc3rd) {
      case 2: // bbb
	theFcTrip = 0;
	break;
      case 1: // bbc
	theFcTrip = 1;
	break;
      case 0: // bbq
	theFcTrip = 2;
	break;
      default:
        //	std::cout << "Bad 3rd flavor code " << theFcDijet << " " << theFc3rd << std::endl;
        break;

      }
      break;
    case 1:  // (bc)
      switch (theFc3rd) {
      case 2: // bcb
	theFcTrip = 3;
	break;
      case 1:
      case 0:
	theFcTrip = 5; // non-bb
	break;
      default:
        //	std::cout << "Bad 3rd flavor code " << theFcDijet << " " << theFc3rd << std::endl;
        break;

      }
      break;
    case 2:  // (bq)
      switch (theFc3rd) {
      case 2: // bqb
	theFcTrip = 4;
	break;
      case 1:
      case 0:
	theFcTrip = 5; // non-bb
	break;
      }
      break;
    case 3:
    case 4:
    case 5:
      theFcTrip = 5; // non-bb
      break;
    default:
      //       std::cout<<"theFcTrip="<<theFcTrip<<std::endl;
      //      std::cout << "Bad flavor codes " << theFcDijet << " " << theFc3rd << std::endl;
      break;
    }

    hfctrip->Fill(theFcTrip,weight);

    // determine the svMass and svMassIndex

    ///FIX is svMassIndex  correct?!


    float svMassIndex[nSelJet];
    for (int iJ=0; iJ<nSelJet; ++iJ) {
      svMassA[iJ]->fill(trgSelect,svMassJet[leadingJets[iJ]],weight);
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	// svMass under offline btag
	if (theBJetTag[ibtag][leadingJets[iJ]]>btcut[ibtag]) {
	  svMassBtA[ibtag][iJ]->fill(trgSelect,svMassJet[leadingJets[iJ]],weight);
	}
      }
      svMassIndex[iJ] = sv->getBinFromSvMass( svMassJet[leadingJets[iJ]] );
    }
    
    //closure test for offline btag efficiencies in the double btag sample
    if ( (nJet>=nSelJet) && deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && (deltaEta < deltaEtaCut12)) {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
        // require TWO btags for the two leading jets
        if( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag]) )
          {
            for (int iJ=0; iJ<nSelJet; ++iJ) {
              DoubleOffBtagptJetBtA[ibtag][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
              if (theFcTrip>=0)
                {
                  // ptJetFcBtA[ibtag][theFcTrip][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                  DoubleOffBtagptJetBtA_Fc[ibtag][theFcTrip][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
                }
            }
          }
        if (_doMC) {
          double wbt[nSelJet];
          for (int iJ=0; iJ<nSelJet; ++iJ) {
            wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
          }
          const double wtotDoubleTag = weight * 
            (
             wbt[0] * wbt[1]
             ); //double tag weight?
          for (int iJ=0; iJ<nSelJet; ++iJ) {
            DoubleOffBtagptJetBtwA[ibtag][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
            if (theFcTrip>=0)
              {
                DoubleOffBtagptJetBtwA_Fc[ibtag][theFcTrip][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],wtotDoubleTag);
              }
          }
        }
      }
    }

    // compute invariant mass of two leading jets, and event btag
    if ( (nJet>=nSelJet) && deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && (deltaEta < deltaEtaCut12) && tripleOnlineBtagMatchOK ) {
      float energyTot = energyJet[leadingJets[0]] + energyJet[leadingJets[1]];
      float pxTot = pxJet[leadingJets[0]] + pxJet[leadingJets[1]];
      float pyTot = pyJet[leadingJets[0]] + pyJet[leadingJets[1]];
      float pzTot = pzJet[leadingJets[0]] + pzJet[leadingJets[1]];
      
      float dijet_mass_sq = energyTot * energyTot - pxTot * pxTot - pyTot * pyTot - pzTot * pzTot;

      if (dijet_mass_sq >= 0) {
	float dijet_mass = sqrt( dijet_mass_sq );

	// determine the deltaPhi
	float dPhiJet1Jet3 = phiJet[leadingJets[0]] - phiJet[leadingJets[2]];
	if (dPhiJet1Jet3>3.1415926) dPhiJet1Jet3 -= 2*3.1415926;
	if (dPhiJet1Jet3<-3.1415926) dPhiJet1Jet3 += 2*3.1415926;
	float dPhiJet2Jet3 = phiJet[leadingJets[1]] - phiJet[leadingJets[2]];
	if (dPhiJet2Jet3>3.1415926) dPhiJet2Jet3 -= 2*3.1415926;
	if (dPhiJet2Jet3<-3.1415926) dPhiJet2Jet3 += 2*3.1415926;

	// fill the triple-btag histograms
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	  if ( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) &&  (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag])
	       && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {
	    if (ibtag == 0) nJetPostselA->fill(trgAccept,float(numberOfJets),1);
	    mDijetBtA[ibtag]->fill(trgSelect,dijet_mass,weight);
	    // triple-jet specific flavor code histograms
	    if (theFcTrip>=0) mDijetFcBtA[ibtag][theFcTrip]->fill(trgSelect,dijet_mass,weight);
	    // fc counting
	    hfc[ibtag]->Fill( theFc,weight );
	    if ( (dijet_mass > 100) && (dijet_mass < 140 ) ) hfcm[ibtag]->Fill( theFc,weight );
	    for (int iJ=0; iJ<nSelJet; ++iJ) {
	      ptJetBtA[ibtag][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
	      if (theFcTrip>=0) ptJetFcBtA[ibtag][theFcTrip][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
	    }
	    // tpat counting
	    int theTpat = 4 * int(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run)) + 2 * int(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run))
	      + int(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run));
	    atpattripall[ibtag]->fill(trgSelect,float(theTpat),weight);
	  }


          //std::cout<<"after dijetmass"<<std::endl;

	  if (_doMC) {
	    // now do the same without btag cut and with weighting
	    double wbt[nSelJet];
	    for (int iJ=0; iJ<nSelJet; ++iJ) {
	      wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
	    }
	    double wtot = weight * wbt[0] * wbt[1] * wbt[2];
	    mDijetBtwA[ibtag]->fill(trgSelect,dijet_mass,wtot);
	    // triple-jet specific flavor code histograms
	    if (theFcTrip>=0) mDijetFcBtwA[ibtag][theFcTrip]->fill(trgSelect,dijet_mass,wtot);
	    // fc counting
	    hfcww[ibtag]->Fill( theFc,wtot );
	    if ( (dijet_mass > 100) && (dijet_mass < 140 ) ) hfcmww[ibtag]->Fill( theFc,wtot );
	    for (int iJ=0; iJ<nSelJet; ++iJ) {
	      ptJetBtwA[ibtag][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],wtot);
	      if (theFcTrip>=0) ptJetFcBtwA[ibtag][theFcTrip][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],wtot);
	    }
	    // fill also the deltaPhi histograms
	    dPhiJet1Jet3BtwA[ibtag]->fill(trgSelect,dPhiJet1Jet3,wtot);
	    dPhiJet2Jet3BtwA[ibtag]->fill(trgSelect,dPhiJet2Jet3,wtot);
	  } ///_domC
	} ///ibtag

	// compute the event btag
	float svMassJets[3];
	for (int iJ=0; iJ<nSelJet; ++iJ) {
	  svMassJets[iJ] = svMassJet[leadingJets[iJ]];
	  if (svMassJets[iJ]<=0) svMassJets[iJ] = 0.01;
	  if (svMassJets[iJ]>=6) svMassJets[iJ] = 5.99;
	}
	int nEvBTag = sv->eventXBTag(svMassJets);

	// fill the event btag and 2D histograms (triple btag)
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	  if ( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) &&  (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag])
	       && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {
	    evBtagBtA[ibtag]->fill(trgSelect,nEvBTag,weight);
	    massEvBtagBtA[ibtag]->fill(trgSelect,dijet_mass,float(nEvBTag),weight);
	  }
	}

	// in case of MC, compute trigger-weighted 2D histograms
	if (_doMC) {
	  float thePt[nSelJet];
	  float theEta[nSelJet];
	  for (int iJ=0; iJ<nSelJet; ++iJ) {
	    thePt[iJ] = ptJet[leadingJets[iJ]];
	    theEta[iJ] = etaJet[leadingJets[iJ]];
	  }

	  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	    if ( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) 
		 &&  (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag])
		 && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {
	      for (unsigned int ib=0; ib<genericTriggerList.size(); ++ib) {
		triggerWeight[ib] = 0;
	      }
	      theHbbTrigWeight->getTrigWeight(bTagReleffOnline,sbtag[ibtag].c_str(),
					      theFlav,thePt,theEta,triggerWeight);
	      
	      // fill with trigger weights
	      massEvBtagBtTWA[ibtag]->fillTW(dijet_mass,float(nEvBTag),triggerWeight,weight);

	      // monitor L1L2-only efficiency
	      theHbbTrigWeight->getTrigWeight(bTagReleffOnline,sbtag[ibtag].c_str(),
					      theFlav,thePt,theEta,triggerWeight,true);
	      triggerWeightA->fillMonitor(triggerWeight,1);
	    }
	  }
	}

        //			std::cout<<"bkg templ"<<std::endl;

	// background templates
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	  // set category membership for this btag selection
	  bool lcateg[3] = {false, false, false};
	  int firstB[ncateg] = { 1, 0, 0 };  // to address b jets in double btag in pt ordered way
	  int secondB[ncateg] = { 2, 2, 1 };
	  int* theB[2];
	  theB[0] = firstB;
	  theB[1] = secondB;
	  if ( (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {
	    lcateg[0] = true;
	  }
	  if ((theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {
	    lcateg[1] = true;
	  }
	  if ((theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) && (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag]) ) {
	    lcateg[2] = true;
	  };
	  // determine the online btag pattern
	  int theTpat = 4 * int(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run)) + 2 * int(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run))
	    + int(HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run));
	  for (int icateg=0; icateg<ncateg; ++icateg) {
	    if (lcateg[icateg]) tpatA[ibtag][icateg]->fill(trgSelect,float(theTpat),weight);
	  }
	  // for cases with double online btag, we define the itpat
	  int thisTpat = -1;
	  switch (theTpat) {
	  case 3:
	    thisTpat = 2; // TTx
	    break;
	  case 5:
	    thisTpat = 1; // TxT
	    break;
	  case 6:
	    thisTpat = 0; // xTT
	    break;
	  case 7:
	    thisTpat = 3; // TTT
	    break;
	  }

	  tpatAllA[ibtag]->fill(trgSelect,float(theTpat),weight);


          ///Here I (Igor) added a few lines to fill flavor dependent plots of M12 
          //cout<<"Igor code"<<endl;
          if (_doMC) {

            int mistagFlav1 = theFlav[0];
            int mistagFlav2 = theFlav[1];

            // btag scaling factor for the untagged jet
            double _theSF1 = 1;
            double _theSF2 = 1;

            if (doBtagSF) {
              theSystControl(1,0,ibtag,mistagFlav1,ptJet[leadingJets[0]],etaJet[leadingJets[0]],_theSF1);
              theSystControl(1,0,ibtag,mistagFlav2,ptJet[leadingJets[1]],etaJet[leadingJets[1]],_theSF2);
            }

            float _mistagWeight1 = bTagEffOffline->eff(mistagFlav1,sbtag[ibtag].c_str(),ptJet[leadingJets[0]],etaJet[leadingJets[0]]) * _theSF1;
            //          float _errorMistagWeight1 = bTagEffOffline->erreff(mistagFlav1,sbtag[ibtag].c_str(),ptJet[leadingJets[0]],etaJet[leadingJets[0]]) * _theSF1;
            float _mistagWeight2 = bTagEffOffline->eff(mistagFlav2,sbtag[ibtag].c_str(),ptJet[leadingJets[1]],etaJet[leadingJets[1]]) * _theSF2;
            //   float _errorMistagWeight2 = bTagEffOffline->erreff(mistagFlav2,sbtag[ibtag].c_str(),ptJet[leadingJets[1]],etaJet[leadingJets[1]]) * _theSF2;

            float _relmistagWeight1= bTagReleffOnline->eff(mistagFlav1,sbtag[ibtag].c_str(),ptJet[leadingJets[0]],etaJet[leadingJets[0]]);
            float _relmistagWeight2= bTagReleffOnline->eff(mistagFlav2,sbtag[ibtag].c_str(),ptJet[leadingJets[1]],etaJet[leadingJets[1]]);

            int dijet = diJetFlavorCode( leadingJets[0],leadingJets[1]);
            mDiJetMassFlavor[dijet][ibtag]->Fill(dijet_mass, weight*_mistagWeight1*_mistagWeight2*_relmistagWeight1*_relmistagWeight2);

          } ///_doMC



          //std::cout<<"Before ifc"<<endl;
	  // fill the templates
	  for (int ifc=0; ifc<nfc; ++ifc) {

            double mistagFlav=ifc;

            ///replace 'c' by 'cc', 'b' by 'bb'
            //	if (ifc>0) mistagFlav+=2;		


	    for (int icateg=0; icateg<ncateg; ++icateg) {
	      if (thisTpat == icateg || thisTpat == 3) {  // online btag pattern must match offline
		if (lcateg[icateg]) {
		  // btag scaling factor for the untagged jet
		  double theSF = 1;
		  if (doBtagSFVar_bc && (ifc >= 1)) {
                    theSystControl(2,BtagSFVar,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],theSF);
		  } else if (doBtagSFVar_q && (ifc ==0)) {
		    theSystControl(3,BtagSFVar,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],theSF);
		  } else if (doBtagSF) {
		    theSystControl(1,        0,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],theSF);
		  }
                  //		  float mistagWeight = bTagEffOffline->eff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]) * theSF;
                  //		  float errorMistagWeight = bTagEffOffline->erreff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]) * theSF;
		  float mistagWeight = bTagEffOffline->eff(mistagFlav,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]) * theSF;
		  float errorMistagWeight = bTagEffOffline->erreff(mistagFlav,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]) * theSF;

		  // mass templates
		  for (int icorr=0; icorr<ncorr; ++icorr) {
		    double wbbpur = 1;
		    if (bbPurityCorr[ibtag] && (icorr == 1)) wbbpur =  fbbfrac[ibtag][icateg]->Eval( dijet_mass );
		    for (int itpat=0; itpat<3; ++itpat) {
		      double wtpat = 1;
			
			
                      //			std::cout<<"_doOnlineRelBtag"<<std::endl;
                      //			std::cout<<"icor "<<icorr<<std::endl;
                      //			std::cout<<"ifc "<<ifc<<std::endl;

                      if (_doOnlineRelBtag)
                        if ( (itpat != icateg) && thisTpat != 3) {  // weight correction if offline & online btag pattern are different (except TTT)
                          if (bTagReleffOnline->eff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]) == 0) {       
                            std::cout << "call eff icateg=" << icateg << " ptJet=" << ptJet[leadingJets[icateg]] 
                                      << " etaJet=" << etaJet[leadingJets[icateg]] << std::endl;
                          }
                          if (bTagReleffOnline->eff(2,sbtag[ibtag].c_str(),ptJet[leadingJets[itpat]],etaJet[leadingJets[itpat]]) == 0) {
                            std::cout << "call eff  itpat=" << itpat <<  " ptJet=" << ptJet[leadingJets[itpat]] 
                                      << " etaJet=" << etaJet[leadingJets[itpat]] << std::endl;
                          }
		      
                          //			wtpat = bTagReleffOnline->eff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]])
                          //			  / bTagReleffOnline->eff(2,sbtag[ibtag].c_str(),ptJet[leadingJets[itpat]],etaJet[leadingJets[itpat]]);
                          wtpat = bTagReleffOnline->eff(mistagFlav,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]])
                            / bTagReleffOnline->eff(2,sbtag[ibtag].c_str(),ptJet[leadingJets[itpat]],etaJet[leadingJets[itpat]]);
                        }
                      //			std::cout<<"wtpat "<<wtpat<<std::endl;
                      //			std::cout<<"wbbpur "<<wbbpur<<std::endl;

                      //			 if (ifc>2)	cout<<"ifc="<<ifc<<endl;
		      massTemplateA[ifc][ibtag][icateg][icorr][itpat]->fill(trgSelect,dijet_mass,weight * mistagWeight * wbbpur * wtpat);

		      // fill the predictions
                      if (_doBkgPred)
                        if (icorr == 1) {
                          float theFlavFrac[3] = { .95, .03, .02 };  // replace by actual numbers
                          massPred[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,weight * mistagWeight * theFlavFrac[ifc] * wtpat);
                        }


                      ///FIX EvtBtag!!!
                      double mistagSV=ifc;
                      //	if (mistagSV>2) mistagSV-=2;

		      // event btag templates (no purity correction yet)
		      if (icorr == 1) {
			// get the svMassIndex probability vector for jet "categ", flavor=ifc

			float probSvMassIndexOfflineIfc[3]; /// FIX change to 5?!
			sv->getSVbins(mistagSV,ibtag,ptJet[leadingJets[icateg]],fabs(etaJet[leadingJets[icateg]]),probSvMassIndexOfflineIfc);
			float probSvMassIndexOnlineIfc[3]; /// FIX change to 5?!
			svOnline->getSVbins(mistagSV,ibtag,ptJet[leadingJets[icateg]],fabs(etaJet[leadingJets[icateg]]),probSvMassIndexOnlineIfc);
			    
			// same for jet "itpat", flavor=b
			float probSvMassIndexOfflineB[3]; 
			sv->getSVbins(2,ibtag,ptJet[leadingJets[itpat]],fabs(etaJet[leadingJets[itpat]]),probSvMassIndexOfflineB);
			float probSvMassIndexOnlineB[3]; 
			svOnline->getSVbins(2,ibtag,ptJet[leadingJets[itpat]],fabs(etaJet[leadingJets[itpat]]),probSvMassIndexOnlineB);

                        //cout<<"svMassIndexOfTemplate"<<endl;

			for (int iSvM=0; iSvM<3; ++iSvM) {  // loop over possible svMassIndex values
			  int svMassIndexOfTemplate[nSelJet];
			  for (int iJ=0; iJ<nSelJet; ++iJ) {
			    if (iJ == icateg) {
			      svMassIndexOfTemplate[iJ] = iSvM;
			    } else {
			      svMassIndexOfTemplate[iJ] = svMassIndex[iJ]; // copy 
			    }
			  }


			  // make event btag
			  int evtBTag = sv->eventBTag( svMassIndexOfTemplate );

                          //			std::cout<<"Get evtBTag"<<std::endl;
			  // extra weight according to online btag pattern
			  float wtpatEvtBtag = 1;
			  float wtpatEvtBtagSum = 0;
			  if ( (itpat != icateg) && thisTpat != 3) {  // weight correction if offline & online btag pattern are different (except TTT)
			    if ( probSvMassIndexOnlineB[iSvM] > 0) {
			      // to model this background, the untagged jet must take the role of an online btagged jet
			      // plus, correct for the fact that one b jet is no longer online btagged
			      wtpatEvtBtag = probSvMassIndexOnlineIfc[iSvM] *  probSvMassIndexOfflineB[iSvM] / probSvMassIndexOnlineB[iSvM];
			    } else {
			      wtpatEvtBtag = probSvMassIndexOnlineIfc[iSvM]; // omit second correction if we would divide by zero
			    }


			  } else {
			    wtpatEvtBtag = probSvMassIndexOfflineIfc[iSvM]; // simplest case: jet is not biased by online btag
                                                                            // this also includes the borderline case of three online btags
			  }
			  wtpatEvtBtagSum += wtpatEvtBtag;

			  bTagTemplateA[ifc][ibtag][icateg][itpat]->fill(trgSelect,float(evtBTag),weight * mistagWeight * wbbpur * wtpat * wtpatEvtBtag);
			  massBTagTemplateA[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,float(evtBTag),weight * mistagWeight * wbbpur * wtpat * wtpatEvtBtag);
			  errorMassBTagTemplateA[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,float(evtBTag),weight * errorMistagWeight * wbbpur * wtpat * wtpatEvtBtag);

			  // fill the predictions
                          //			 std::cout<<"pred beg"<<std::endl;
			  float theFlavFrac[3] = { .95, .03, .02 };  // replace by actual numbers
                          if (_doBkgPred)
                            massBTagPred[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,float(evtBTag),weight * mistagWeight * theFlavFrac[ifc] * wtpat * wtpatEvtBtag);

                          //			 std::cout<<"pred end"<<std::endl;


			  // data-driven bb purity correction
			  bool hasNegBTag[2];
			  float Rmistag[2];
			  const int bbPurity_Version_VarThresh = 3;
			  for (int iiJ=0; iiJ<2; ++iiJ) {
			    hasNegBTag[iiJ] = ifNegBTag( ntcHEBJetTag[leadingJets[theB[iiJ][icateg]]],
                                                         ptJet[leadingJets[theB[iiJ][icateg]]],
                                                         sbtag[ibtag],icateg,iiJ,bbPurity_Version_VarThresh,
                                                         Rmistag[iiJ] );
			  }
			  massBTagTemplateUncldA[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,float(evtBTag),weight * mistagWeight * wtpat * wtpatEvtBtag);
			  if (! (hasNegBTag[0] || hasNegBTag[1] )) {
			    massBTagTemplateCldA[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,float(evtBTag),weight * mistagWeight * wtpat * wtpatEvtBtag);
			  }

			  // fixed-threshold method
			  const int bbPurity_Version_FixedThresh = 4;
			  for (int iiJ=0; iiJ<2; ++iiJ) {
			    hasNegBTag[iiJ] = ifNegBTag( ntcHEBJetTag[leadingJets[theB[iiJ][icateg]]],
                                                         ptJet[leadingJets[theB[iiJ][icateg]]],
                                                         sbtag[ibtag],icateg,iiJ,bbPurity_Version_FixedThresh,
                                                         Rmistag[iiJ] );
			  }
			  float weightR = 0;
			  if ( hasNegBTag[0] && (! hasNegBTag[1] )) {
			    weightR = Rmistag[0];
			  } else if ( hasNegBTag[1] && (! hasNegBTag[0] )) {
			    weightR = Rmistag[1];
			  } else if ( hasNegBTag[0] && hasNegBTag[1] ) {
			    weightR = Rmistag[0] + Rmistag[1] - Rmistag[0]*Rmistag[1];
			  }
			  if ( hasNegBTag[0] || hasNegBTag[1] ) {
			    massBTagTemplateNonbbRA[ifc][ibtag][icateg][itpat]->fill(trgSelect,dijet_mass,float(evtBTag),weight * mistagWeight * wtpat * wtpatEvtBtag * weightR);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }	  
	}  // end of templates

        //std::cout<<"end of templ"<<std::endl;

	// here we fill the double btag histograms
	// auxiliary array to perform permutations among first three jets
	int shift[6] = { 0, 1, 2, 0, 1, 2 };
	
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {


	  for (int icateg=0; icateg<ncateg; ++icateg) {
            //cout<<"Before icateg"<<endl;

	    int jetA = leadingJets[shift[icateg+1]];
	    int jetB = leadingJets[shift[icateg+2]];

	    // first the cut-based hist independent of flavor code
	    if ( (theBJetTag[ibtag][jetA]>btcut[ibtag]) &&  (theBJetTag[ibtag][jetB]>btcut[ibtag]) ) {
	      mDibBtcutA[ibtag][icateg]->fill(trgSelect,dijet_mass, weight );
	      for (int iJ=0; iJ<nSelJet; ++iJ) {
		ptDibBtcutA[ibtag][icateg][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight);
	      }
	    }

	    if (_doMC) {
	      float wDoubleBtag = bTagEffOffline->eff(theFlav[shift[icateg+1]],sbtag[ibtag].c_str(),ptJet[jetA],etaJet[jetA])
		* bTagEffOffline->eff(theFlav[shift[icateg+2]],sbtag[ibtag].c_str(),ptJet[jetB],etaJet[jetB]);
	      mDibBtweightA[ibtag][icateg]->fill(trgSelect,dijet_mass, weight * wDoubleBtag);
	      for (int iJ=0; iJ<nSelJet; ++iJ) {
		ptDibBtweightA[ibtag][icateg][iJ]->fill(trgSelect,ptJet[leadingJets[iJ]],weight * wDoubleBtag );
	      }
		
              //cout<<"Before second diJetCode"<<endl;
	      int dibFcDijet = diJetFlavorCode( jetA, jetB );
	      if ( (theBJetTag[ibtag][jetA]>btcut[ibtag]) &&  (theBJetTag[ibtag][jetB]>btcut[ibtag]) ) {
                //cout<<"Before fill dijet"<<endl;
                //cout<<dibFcDijet<<endl;
                //cout<<"jetA="<<jetA<<endl;
                //cout<<"jetB="<<jetB<<endl;
                //cout<<"codeA="<<jetFlavorCodeNewCat( jetA)<<endl;
                //cout<<"2 codeA="<<theFlav[jetA]<<endl;
                //cout<<"codeB="<<jetFlavorCodeNewCat( jetB)<<endl;
                //cout<<"2 codeB="<<theFlav[jetB]<<endl;

		mDibBtcutFcH[dibFcDijet][ibtag][icateg]->Fill( dijet_mass, weight );
		// also fill the pt histograms under double btag selection
                //cout<<jetFlavorCode(jetA)<<endl;
                //cout<<"jetA"<<endl;
                //cout<<jetFlavorCode(jetB)<<endl;
                //cout<<"jetB"<<endl;
		hptdibbt[jetFlavorCode(jetA)][0][ibtag][icateg]->Fill( ptJet[jetA], weight );
		hptdibbt[jetFlavorCode(jetB)][1][ibtag][icateg]->Fill( ptJet[jetB], weight );
	      }
              //cout<<"Before fill diJetCode"<<endl;
              //cout<<"icateg="<<icateg<<endl;
              //cout<<ptJet[jetA]<<endl;
              //cout<<ptJet[jetB]<<endl;
              //cout<<theFlav[shift[icateg+1]]<<endl;
              //cout<<theFlav[shift[icateg+2]]<<endl;
              //cout<<"dibFcDijet="<<dibFcDijet<<endl;
              //cout<<"Edn"<<endl;

	      // same with btag weighting. New: include online btag weight
	      mDibBtweightFcH[dibFcDijet][ibtag][icateg]->Fill( dijet_mass, weight
                                                                * bTagEffOffline->eff(theFlav[shift[icateg+1]],sbtag[ibtag].c_str(),ptJet[jetA],etaJet[jetA])
                                                                * bTagEffOffline->eff(theFlav[shift[icateg+2]],sbtag[ibtag].c_str(),ptJet[jetB],etaJet[jetB])
                                                                // 	       * bTagReleffOnline->eff(theFlav[shift[icateg+1]],sbtag[ibtag].c_str(),ptJet[jetA],etaJet[jetA])
                                                                // 	       * bTagReleffOnline->eff(theFlav[shift[icateg+2]],sbtag[ibtag].c_str(),ptJet[jetB],etaJet[jetB])  );
								);
              //cout<<"After fill diJetCode"<<endl;

	    }
            //std::cout<<"end of icateg"<<std::endl;
	  } // loop over categ
          //std::cout<<"end of ibtag"<<std::endl;
	} // loop over ibtag
      } else {
        std::cout << "Unphysical dijet_mass_sq = " << dijet_mass_sq << std::endl;
        //dijet_mass = -1;
      } // if (dijet_mass...)
    } // if (njet...)
  } // loop over the tree
  // termination
  
  std::cout << "Here we are in the termination" << std::endl;
  // add the different tpat contributions, taking correlations into account
  for (int ifc=0; ifc<nfc; ++ifc) {
    for (int ibtag=0; ibtag<nbtag; ++ibtag) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	for (int icorr=0; icorr<ncorr; ++icorr) {
	  // for mass template
	  TrigHistArray::mergeHistos( 3, massTemplateA[ifc][ibtag][icateg][icorr] );
	}
	// for event btag template
	TrigHistArray::mergeHistos( 3, bTagTemplateA[ifc][ibtag][icateg] );
	// for 2D mass / event btag template
	TrigHistArray2D::mergeHistos( 3, massBTagTemplateA[ifc][ibtag][icateg] );
	TrigHistArray2D::mergeHistos( 3, errorMassBTagTemplateA[ifc][ibtag][icateg] );
	TrigHistArray2D::mergeHistos( 3, massBTagTemplateUncldA[ifc][ibtag][icateg] );
	TrigHistArray2D::mergeHistos( 3, massBTagTemplateCldA[ifc][ibtag][icateg] );
	TrigHistArray2D::mergeHistos( 3, massBTagTemplateNonbbRA[ifc][ibtag][icateg] );
	// for prediction histograms
        // 	TrigHistArray::mergeHistos( 3, massPred[ifc][ibtag][icateg] );
        // 	TrigHistArray2D::mergeHistos( 3, massBTagPred[ifc][ibtag][icateg] );
	// data-driven subtraction of non-bb (fixed threshold method). Only for tpat-corrected hist.
	int mtpat = 3;
	std::cout << "Call the adder for ifc=" << ifc << " ibtag=" << ibtag
		  << " icateg=" << icateg 
		  << " mtpat=" << mtpat << std::endl;
	massBTagTemplateCldRA[ifc][ibtag][icateg][mtpat]->add( massBTagTemplateUncldA[ifc][ibtag][icateg][mtpat],
							       massBTagTemplateNonbbRA[ifc][ibtag][icateg][mtpat],
							       1, -1);
      }
    }
  }

  hout->Write();
  hout->Close();
}
