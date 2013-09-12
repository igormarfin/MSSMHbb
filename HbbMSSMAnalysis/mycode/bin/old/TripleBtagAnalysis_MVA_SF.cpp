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
#include <TFormula.h>

//#include "Analysis/Utilities/interface/HBBTo4B.h"
//#include "Analysis/Utilities/src/HBBTo4B.cc"

#include "Analysis/Utilities/interface/svMass.h"
#include "Analysis/Utilities/interface/bTagEff.h"
#include "Analysis/Utilities/interface/TrigHistArray.h"
#include "Analysis/Utilities/interface/TrigHistArray2D.h"
#include "Analysis/Utilities/interface/getHbbMetadata.h"
#include "Analysis/Utilities/interface/jetFlavorCode.h"
#include "Analysis/Utilities/interface/HbbNtuple.h"
#include "Analysis/Utilities/interface/HbbTrigWeight.h"





///_MVA_BEGIN_ headers

#include "Analysis/HbbMSSMAnalysis/interface/MVAReadInput.h"


#if defined(CMSSW_CPP)

#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

#endif




#include <fstream>
#include <TString.h>
#include <TSystem.h>
#include <TROOT.h>
#include <map>
#include <TObjArray.h>
#include <TLorentzVector.h>


using namespace std;



#include "Analysis/HbbMSSMAnalysis/test/NegativeBTag.C"



///_MVA_END_


///_PU_BEGIN_

#if defined(PU)

#include "Analysis/Utilities/interface/PileUpProducer.h"

#endif 

#if defined(PU_ALEX)

#include "Analysis/Utilities/interface/PUWeight.h"

#endif 

///_PU_END_


///_SF_BEGIN_

#if defined(SF)  || defined(SYST) || defined(JER)

#include "Analysis/Utilities/interface/HbbSystControl.h"
#include "Analysis/HbbMSSMAnalysis/interface/SystReadInput.h"

#endif
///_SF_END_




///_OVERLAP_BEGIN

#if defined(TRUNCATED)

#include "Analysis/HbbMSSMAnalysis/interface/Overlap.h"

#endif

///_OVERLAP_END


TCanvas* canvas;




void TripleBtagAnalysis_MVA_SF(TString _mvasettings="Analysis/HbbMSSMAnalysis/data/mvasettings.lst") {

#if !(defined(CMSSW_CPP))

TString __base=gSystem->GetFromPipe("echo $CMSSW_BASE");
__base+="/src/";
_mvasettings=__base+_mvasettings;


#endif

  // macro for analysis in channel with three btagged jets

  canvas = new TCanvas ("cg1","mycanvas",10,10,800,600);
  // open an ntuple file

  bool _doMC = false;  // set to true for MC, false for real data

  bool isMCWithTrigSim = true;
  bool doLumiWeighting = true;  // weight MC events according to intgrated lumi


//set it to false for data
  bool forceTripleOnlineBtagMatch = true;  // require triple online btag match on triple btag signature 

  const int nSelJet = 3;   // number of jets considered in b(bb) final state

  // chain mode
  TChain f("hbbanalysis/HBBTo4B");
  TFileCollection theFileColl("myNtFileCollection","","theMergeList.txt");
  if (theFileColl.GetNFiles() <=0) {
    std::cout << "Problem opening list of input files " << std::endl;
    return;
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
  std::vector<std::string> triggerFilterList;
  if (_doMC) {
    if (isMCWithTrigSim) {
      //triggerFilterList.push_back("HLT_bbPhi_CentralJet46_CentralJet38_DiBTagIP3D_L25MinTag4_v1");
      //triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D"); 

///Low mass scenario (2011)
//      triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D");
//      triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D");


///Medium mass scenario (2011)
//      triggerFilterList.push_back("HLT_CentralJet60_CentralJet53_DiBTagIP3D");


	/// 2012 analysis

	 triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV");
	 triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV");
	 triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose");


    } else {
      triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D");
    }
  } else {

///2011 data
/*
    triggerFilterList.push_back("HLT_CentralJet46_BTagIP3D_CentralJet38_BTagIP3D");
    triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_DiBTagIP3D");
    triggerFilterList.push_back("HLT_CentralJet46_CentralJet38_CentralJet20_DiBTagIP3D");
     triggerFilterList.push_back("HLT_CentralJet60_CentralJet53_DiBTagIP3D");
*/
    /// 2012 analysis



     triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV");
     triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV");
     triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose");

//     triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3D");
//     triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3D");
//     triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DLoose");
  }

  // extract generic trigger list and number of input events from ntuple files
  float lumiScaleFac = 1;
  int nInputEvents= getHbbMetadata(theFileColl,genericTriggerList,1000.,lumiScaleFac,_doMC);
//  getHbbMetadata(theFileColl,genericTriggerList,1000.,lumiScaleFac,_doMC);

  // this trigger weight array will be used later in each event
  float* triggerWeight = new float[genericTriggerList.size()];
  // create trigger weight helper object
  HbbTrigWeight* theHbbTrigWeight = new HbbTrigWeight(&genericTriggerList,&triggerFilterList);


  f.AddFileInfoList((TCollection*) theFileColl.GetList());
  TTree* hbbtree = &f;

// try HBBTo4B class
//   HBBTo4B * hbbTo4b = new HBBTo4B();
//   hbbTo4b->Init(hbbtree);

//raise the analysis level pt cuts to
//(60,53,20) GeV.



  // kinematic analysis cuts

#if !defined(MEDIUM_MASS)

///low-mass-scenario
  double jetPtMin[nSelJet] = { 46, 38, 20};  // leading jet lower pt thresholds
#else
///medium-mass-scenario
  double jetPtMin[nSelJet] = { 60, 53, 20};  // leading jet lower pt thresholds

#endif

  double jetPtMax[nSelJet] = {3500, 3500, 3500}; // leading jet upper pt thresholds
  double maxEta = 2.2;

  // here we define dimensions and labels
  const int nflav=3;  // number of flavor groups to be distinguished
  const int nbtag = 4;  // number of offline btag working points
  const int ncateg = 3; // ranks of untagged jet
  const int ncorr = 2; // correction levels
  const int ntpat = 4; // trigger pattern (online btag). (0,1,2)=rank of non-btagged jet, 3=all combined

 std::string onlineBtagVersion("V4"); // choices: smpf or V4


  const std::string flavlabel[nflav] = {"udsg","c","b"};  // labels for flavor groups
  // define names of btag variants
  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };
  // corresponding discriminant cuts
  double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags
  bool bbPurityCorr[nbtag] = { true, true, true, true };  // specify if bbPurity collection should be applied

  // define flavor pair classes
  const int nfcDijet = 6;
  string sfcDijet[nfcDijet] = { "bb", "bc", "bq", "cc", "cq", "qq" };
  const int nfc3rd = 3;
  string sfc3rd[nfc3rd] = { "q", "c", "b" };
  // flavor triplet classes
  const int nfcTrip = 6;
  string sfcTrip[nfcTrip] = { "bbb", "bbc", "bbq", "bcb", "bqb", "non-bb" };
  // the following are used for templates
  const int nfc = 3;  
  string sfc[nfc] = { "q", "c", "b" };

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
  bTagEff* bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v3.root","offline",flavlabel,sbtag,nflav,nbtag); 
///afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt-N14/btagMatrixOffline-csv/plotmistag-b.root
//","offline",flavlabel,sbtag,nflav,nbtag); 

  // online btag object (with smp smoothing)
 //   bTagEff* bTagReleffOnline = new bTagEff("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt_bEnriched-T14b/btagMatrixOnline-csv/smooth/btagMatrixOnline-csv-smpf.root","online",flavlabel,sbtag,nflav,nbtag);

  // online btag object (with smp smoothing)
  bTagEff* bTagReleffOnline = NULL;
  if (onlineBtagVersion == "smpf") {
    bTagReleffOnline  = new bTagEff("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/results/job014/QCD_Pt_bEnriched-T14b/btagMatrixOnline-csv/smooth/btagMatrixOnline-csv-smpf.root","online",flavlabel,sbtag,nflav,nbtag);
  } else if (onlineBtagVersion == "V4") {
    bTagReleffOnline  = new bTagEff("/afs/naf.desy.de/user/r/rmankel/public/btagMatrix/btagMatrixOnline-V4.root","online",flavlabel,sbtag,nflav,nbtag);
  } else {
    std::cout << "No such online btag version available, " << onlineBtagVersion
              << std::endl;
    return ;
  }
  std::cout << "Using online btag version: " << onlineBtagVersion << std::endl;


  // open the bbPurity correction functions
  TF1* fbbfrac[nbtag][ncateg];
  TFile* bbPur = new TFile("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_3_patch2/src/DesyHiggsAnalyses/HBBAnalyses/test/calib/bbPurity-csv-online.root");
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      fbbfrac[ibtag][icateg] = (TF1*) bbPur->Get( Form("fbbfrac-ww-%s-Cat%d",sbtag[ibtag].c_str(),icateg) );
      if ( fbbfrac[ibtag][icateg] == NULL ) {
	std::cout << "bbPur correction function not found for" << sbtag[ibtag].c_str() << " categ " << icateg << std::endl;
	if (bbPurityCorr[ibtag]) return;
      }
    }
  }

  // open the SVMass template files
  //TString SVMassFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_svMass.root");
  //TString SVMassOnlineFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_BTagHLTmatched_svMass.root");
  TString SVMassFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_svMass.root");
  TString SVMassOnlineFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_BTagHLTmatched_svMass.root");

  svMass * sv = new svMass(SVMassFileName);
  svMass * svOnline = new svMass(SVMassOnlineFileName);


  // create the Root output file
  TFile* hout = new TFile("TripleBtagAnalysis.root","recreate");

  TH1::AddDirectory(true);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();  // not really needed (we inherit from TH1)

  Int_t nentries = (Int_t) hbbtree->GetEntries();
  cout << "Number of events in ntuple: " << nentries << endl;

#include "Analysis/Utilities/interface/HbbNtuple.cc"

  // book the general histograms (before btag selection)
  hout->mkdir("general");
  hout->cd("general");
  TrigHistArray* ptJetA[3];
  TrigHistArray* etaJetA[3];
  TrigHistArray* phiJetA[3];
  TrigHistArray* tcheJetA[3];
  TrigHistArray* tchpJetA[3];
  TrigHistArray* isJetWithBtagA[3];

  for (int ii=0; ii<3; ++ii) {
    ptJetA[ii] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("ptj%d",ii),Form("pt of %dth leading jet",ii+1),25,0,200);
    etaJetA[ii] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("aetaj%d",ii),Form("eta of %dth leading jet",ii+1), 100, -3.5, 3.5);
    phiJetA[ii] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("aphij%d",ii),Form("phi of %dth leading jet",ii+1), 100, -3.15, 3.15); 
    tcheJetA[ii] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("tchej%d",ii),Form("TCHE of %dth leading jet, ptjet>50",ii+1), 60, -10, 20 );
    tchpJetA[ii] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("tchpj%d",ii),Form("TCHP of %dth leading jet, ptjet>50",ii+1), 60, -10, 20 );
    isJetWithBtagA[ii] = new TrigHistArray(&genericTriggerList,&triggerFilterList,Form("isWithBtagj%d",ii),Form("IsJetWithBtag for %dth leading jet",ii+1),10,-3.5,6.5);
  }
  TrigHistArray* matchPatternA = new TrigHistArray(&genericTriggerList,&triggerFilterList,"matchPattern","matchPattern",20,-0.5,19.5);  
  TrigHistArray* dPhiJet1Jet2A = new TrigHistArray(&genericTriggerList,&triggerFilterList,"adphij","delta-phi of two leading jets", 100, -3.15, 3.15); 
  hout->cd();

  // book histogram arrays for dijet mass and pt
  //  Btcut    = btag selection by explicit cuts on discriminant
  //  Btweight = btag emulated by applying btag weights (MC only)
  //  Fc       = split into flavor pattern classes (MC only)
  hout->mkdir("mjjBtcut");
  hout->mkdir("mjjBtweight");
  hout->mkdir("mjjBtcutFc");
  hout->mkdir("mjjBtweightFc");

  hout->mkdir("ptjBtcut");
  hout->mkdir("ptjBtcutFc");
  hout->mkdir("ptjBtweight");
  hout->mkdir("ptjBtweightFc");

  hout->mkdir("dPhiBtweight");

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
    hfc[ibtag] = new TH1F(Form("hfc%s",sbtag[ibtag].c_str()),Form("Flavor triple code 3* %s",sbtag[ibtag].c_str()),19,-1.5,17.5);
    hfcm[ibtag] = new TH1F(Form("hfcm%s",sbtag[ibtag].c_str()),Form("Flavor triple code 3* %s, mass 100-140",sbtag[ibtag].c_str()),19,-1.5,17.5);
    hfcww[ibtag] = new TH1F(Form("hfc%sww",sbtag[ibtag].c_str()),Form("Flavor triple code weighted as 3* %s",sbtag[ibtag].c_str()),19,-1.5,17.5);
    hfcmww[ibtag] = new TH1F(Form("hfcm%sww",sbtag[ibtag].c_str()),Form("Flavor triple code weighted as 3* %s, mass 100-140",sbtag[ibtag].c_str()),19,-1.5,17.5);
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

  // double btag purity histograms. These are produced to create bbPurity corrections.
  //  Btcut    = btag selection by explicit cuts on discriminant
  //  Btweight = btag emulated by applying btag weights (MC only)
  hout->mkdir("bbPurityBtcut");
  hout->mkdir("bbPurityBtweight");
  TH1F* hmdibbt[nfcDijet][nbtag][ncateg];
  TH1F* hmdibww[nfcDijet][nbtag][ncateg];
  if (_doMC) {
    for (int ifcDijet=0; ifcDijet<nfcDijet; ++ifcDijet) {
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	for (int icateg=0; icateg<ncateg; ++icateg) {
	  hout->cd("bbPurityBtcut");
	  hmdibbt[ifcDijet][ibtag][icateg] = new TH1F(Form("mdib_%s_%s_Cat%dbt",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
						     Form("Jet1Jet2 %s 2* %s Cat%d",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
						     50,0,500);
	  hout->cd("bbPurityBtweight");
	  hmdibww[ifcDijet][ibtag][icateg] = new TH1F(Form("mdib_%s_%s_Cat%dww",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
						     Form("Jet1Jet2 %s weighted as 2* %s Cat%d",sfcDijet[ifcDijet].c_str(),sbtag[ibtag].c_str(),icateg),
						     50,0,500);
	}
      }
    }
  }
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



///_MVA_BEGIN_ declaration && definition && ini

hout->cd();
hout->mkdir("countSample");
hout->cd("countSample");
TH1 * _count = new TH1F("count","count",1,0,2);
hout->cd();

hout->mkdir("countPUSample");
hout->cd("countPUSample");
TH1 * _countPU = new TH1F("countPUSample","countPUSample",1,0,2);
hout->cd();

hout->mkdir("lumiweight");
hout->cd("lumiweight");
TH1 * _lumiweight = new TH1F("lumiweight","lumiweight",1,0,2);
hout->cd();



#ifdef EXCLUDE_TRAINER_FROM_SELECTION


ExcludeVariableFromMVA();

#endif


ReadMVASettings( _mvasettings);

int _goodEvents=0;
int _afterMVA=0;
///_MVA_END_ ini


///_PU_BEGIN_ ini 


#if defined(PU)

PileUpProducer * _pileup=0;


if (_doMC) {
////Compute pileup reweighting
///version 1

//_pileup = new PileUpProducer("PUProducer","PUProducer","pudist.root","pileup");

/// version 2, take TH1 from data and MC root files
_pileup = new PileUpProducer("PUProducer","PUProducer","MyDataPileupHistogram.root","PileUp_Fall11.root","pileup");

}



#endif


#if defined(PU_ALEX)

  TString baseDir("/afs/naf.desy.de/user/r/rasp/public/HTauTauAnalysis/RooT");
  PUWeight * pileupWeight = new PUWeight(baseDir);
  pileupWeight->InitPUWeightsS6("PileUp_2011.root","PileUp_Fall11.root");

#endif
 

///_PU_END_ ini




///_SF_BEGIN_ ini




#if defined(SF)  || defined(SYST) || defined (JER)


HbbSystControl  _sf_sys;

 



/**
	systype defines what we want to do with HbbSystControl.
	HbbSystControl can do JES and SF systematics and can return SF central value
		systype :
				0 -- make JES systematic variation
				1 -- return SF for all flavors 
				2 -- return SF b/c  central value or SF b/c btag syst. variation
				3 -- return SF udsg central value or SF udsg btag syst. variation
				4 -- make JER systematic variation


		variation in sigmas:
				-2sigma or + 2sigma ; -2.5sigma 2.5sigma; etc
**/


/// by default, it does SF reweighting 
int systype=1;
double variation=0;


SystReadInput(systype,variation);



#endif


///_SF_END_ ini


///_OVERLAP_BEGIN ini 

#if defined(TRUNCATED)

///file 'padova_events' must exist.
Overlap_Ini("padova_events");

#endif

///_OVERLAP_END



///main loop!

_count->SetBinContent(1,nInputEvents);

double _PU_nentries=0.;
double _eff_presel =(nInputEvents>0 && nentries>=0)?(double) nentries/(double) nInputEvents:1.0;
//cout<<"_eff_presel="<<_eff_presel<<endl;

_lumiweight->SetBinContent(1,lumiScaleFac);


double _non3bb=0;
bool _lightIsOk=true;

double _passKinematics=0;
double _passKinematicsMVA=0;
double _passKinematicsMVAAfter=0;
double _passTrig0=0;
double _beforeSys=0;
double _beforeJet=0;
double _afterJet=0;

  // loop over the tree
  for (Int_t iE=0; iE<nentries; iE++) {

//  if (iE>3) return ;

    hbbtree->GetEntry(iE);
    //hbbTo4b->GetEntry(iE);
    //if (iE<100) std::cout << numberOfJets << std::endl;
//	if  (iE>5000) break;
/*     if (run > 5000) {
       std::cout << "Run number above " << 167151 << " --> terminate" << std::endl;
       break;
     }

*/


///Make cuts very loose if we reqire no selection
#if defined(NO_OFFLINESEL)

double jetPtMinLoose[nSelJet] = { 5, 5, 5};  
for(int iii = 0; iii < nSelJet; iii++) jetPtMin[iii]=jetPtMinLoose[iii];
maxEta = 5.2;

double btcutLoose[nbtag] = { -50., -50.,-50., -50. }; 
for(int iii = 0; iii < nbtag; iii++) btcut[iii]=btcutLoose[iii];

#endif



    std::vector<int> leadingJets;

// cout<<"Event : "<<iE<<endl;

#if !defined(RUNEVENTLUMI_LIST)
if (iE%5000==0) { cout<<"Event : "<<iE; 	cout << '\xd';}
#endif



///_OVERLAP_BEGIN

#if defined(TRUNCATED)

if (Overlap_Check(run,event,lumi)) {
//std::cout<<"Overlap Event: "<<run<<" "<<event<<" "<<lumi<<std::endl; 
continue;}

#endif

///_OVERLAP_END


///_MVA_BEGIN_ declaration of objects to feed MVA!


vector<math::XYZTLorentzVector> _vecJets ; /// objects to feed mvacomputer and mvatrainer
vector<double> * res=0; /// response of mva computer
double minPtMVA=5e0; /// some cut on object
double maxEtaMVA = 2.2;
std::vector<int> _theFlav; /// to catch flavor of objects
bool _isplaced=false;
double _offlineBtagWeight=1e0; /// needed by MVATrainer
double _onlineBtagWeight=1e0; /// needed by MVATrainer

#if defined(NO_OFFLINESEL)
maxEtaMVA = 5.2;
#endif

///_MVA_END_



///_SF_BEGIN_ 


double _SF_MC[nbtag]={1.0,1.0,1.0,1.0};

double _JERF=1.0;

#if defined(JER)
if (_doMC) {
	if (systype!=4){
//cout<<"JER is defined"<<endl;
	_sf_sys(4, 0, -1, -1, 0., 0.,_JERF);  // JER
	}
}
#endif



#if defined(SYST) 

//cout<<"Syst defiend"<<endl; 

_sf_sys(systype,variation,-1,-1,0.,0.,_SF_MC[0]); ///possible JES systematic


#endif

_beforeSys++;




///fill _PU_nentries
double _puwgt=1.0;
#if defined(PU)
	_puwgt=(_pileup)?_pileup->GetWeight(nPUInTime):1.0;
#endif

 _PU_nentries+=_puwgt;



///_SF_END_
 
_beforeJet++;

    // find set of leading jets
    int nJet = 0;
    // loop over the jets
    for (int iJet=0; iJet<numberOfJets; ++iJet) {

 _isplaced=false;



      if (nJet<nSelJet)
      if ( (fabs(etaJet[iJet])<maxEta) ) 
      if ( (numberOfConstituentsInJet[iJet] > 1) )      
      if ( (ptJet[iJet] > jetPtMin[nJet]) && (ptJet[iJet] < jetPtMax[nJet]) ) {
	
//	cout<<"Jet["<<nJet<<"] :"<< "Pt="<<ptJet[iJet]<<"  Px= "<<pxJet[iJet]<<endl;

	leadingJets.push_back(iJet);
	++nJet;

///_MVA_BEGIN_ getting of the first nSelJet-1 th objects' info 

	 _theFlav.push_back(jetFlavorCode(iJet)); 
	  _vecJets.push_back(math::XYZTLorentzVector(pxJet[iJet], pyJet[iJet], pzJet[iJet], energyJet[iJet]));
	  _isplaced=true;

//_MVA_END_ 	
      }

///_MVA_BEGIN_ getting of  nSelJet th object info

  	         if (nJet == nSelJet && !_isplaced)
		  if ( (fabs(etaJet[iJet])<maxEtaMVA) )
		  if ( (numberOfConstituentsInJet[iJet] > 1) )
		   if ( ptJet[iJet]  > minPtMVA ) {
		_theFlav.push_back(jetFlavorCode(iJet));
		 _vecJets.push_back(math::XYZTLorentzVector(pxJet[iJet], pyJet[iJet], pzJet[iJet], energyJet[iJet]));

		}	


	
//_MVA_END_ 


    } /// end of loop over jets



_afterJet++;

    if (nJet < nSelJet) continue;

///End of preselection?

_passKinematics++;


    // check the matching flags
    int mJet = 0;
    int matchPat = 0;
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      if (! (ptJet[iJet] > 15) ) continue;
      if (! (fabs(etaJet[iJet])<maxEta) ) continue;
      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
      if (isJetWithHltBtag[iJet]) {
	matchPat = matchPat | (1<<mJet);
      }
      ++mJet;
    }
    matchPatternA->fill(trgAccept,matchPat);

    // compute the weight for this event
    float weight = 1;
    if ( _doMC && doLumiWeighting )  weight = lumiScaleFac;


///_PU_BEGIN_ weight


#if defined(PU)

/// _pileup is non-zero only for the case of data

if (_pileup) {

///v.2, interface is taken from Alexei code

weight*=_pileup->GetWeight(nPUInTime);
//cout<<"PU weight is "<<_pileup->GetWeight(nPUInTime)<<endl;

}


#endif 


/// used only for validation

#if defined(PU_ALEX)

if (pileupWeight)
{
//cout<<"PU from Alexei weight is "<< pileupWeight->PUWeightS6(nPU)<<endl;
cout<<"PU from Alexei weight is "<< pileupWeight->PUWeightS6(nPUInTime)<<endl;
}

#endif



///_PU_END_ weight






///_MVA_BEGIN_ mva selection



//  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
//  string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };
  // corresponding discriminant cuts
//  double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags
//  bool bbPurityCorr[nbtag] = { true, true, true, true };  // specify if bbPurity collection should be applied

///We are working with "TCHPT"


///      int iiovl=1; /// "TCHP6" 
        int iiovl=0; /// "TCHPT" 

///We are working with CSVT
	iiovl=2;

 float _dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
 if (_dphij>3.1415926) _dphij -= 2*3.1415926;
 if (_dphij<-3.1415926) _dphij += 2*3.1415926;
 float _detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[0]];
 float  _deltaR = sqrt( _dphij*_dphij + _detaj*_detaj );

#if defined(NO_OFFLINESEL)
	_deltaR=1e2;
#endif

#if defined(NO_TRIGSEL)
	trgAccept=1;
#endif



bool _isOk=false;


///data
if (!_doMC) {


bool _btag=false;

if ( (theBJetTag[iiovl][leadingJets[0]]>btcut[iiovl]) && (theBJetTag[iiovl][leadingJets[1]]>btcut[iiovl]) &&  (theBJetTag[iiovl][leadingJets[2]]>btcut[iiovl]))
_btag=true;

if ( _deltaR>1  && _btag )  
if (trgAccept  & (1<<0) || trgAccept  & (1<<1) || trgAccept  & (1<<2)) ///triggered data by Trigger0 or Trigger1 or Trigger2
//if (trgAccept  & (1<<9) || trgAccept  & (1<<14) || trgAccept  & (1<<19) ) ///triggered data by Trigger9 or Trigger14 or Trigger19
_isOk=true;



}





///MC
if (_doMC) {
if (  _deltaR>1 )  _isOk=true;

///SF 
///_SF_BEGIN_

#if defined(SF) || defined(SYST)

for (int ibtag=0; ibtag<nbtag; ++ibtag) {

double _SF_2=1; ///multiplicative factor
bool _dosf=false;
double _sfvariation=variation;

#if defined(SYST)  
if (systype==2) _dosf=true;
if (systype==3) _dosf=true;
#endif

if (systype==0) {_dosf=true; _sfvariation=0e0;}
if (systype==1) {_dosf=true; _sfvariation=0e0;}
if (systype==4) {_dosf=true; _sfvariation=0e0;}

if (_dosf) /// if it's SF reweighting or SF systematics 
for (int iJ=0; iJ<nSelJet; ++iJ)
{

if (systype==0) 
{
#if defined(SF)
_sf_sys(1,_sfvariation,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics
#endif
}
else if (systype==4) {
#if defined(SF)
_sf_sys(1,_sfvariation,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics
#endif
}
else if (systype==1)
_sf_sys(systype,_sfvariation,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics
else if (systype==2) {
#if defined(SF)
	if (_theFlav[iJ]==0) _sf_sys(1,0.0,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics
#endif
      if (_theFlav[iJ]>=1) _sf_sys(systype,_sfvariation,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics

}
else if (systype==3) {
#if defined(SF)
        if (_theFlav[iJ]>=1) _sf_sys(1,0.0,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics
#endif
      if (_theFlav[iJ]==0) {
_sf_sys(systype,_sfvariation,ibtag,_theFlav[iJ],ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]],_SF_2); ///possible SF systematics
//	_non3bb+=1;
}

}


_SF_MC[ibtag]*=_SF_2;


}


}
#endif



///_SF_END_





///OfflineBtag (triple-tag) weigting

          for (int iJ=0; iJ<nSelJet; ++iJ) 
          _offlineBtagWeight *=   bTagEffOffline->eff(_theFlav[iJ],sbtag[iiovl].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);        
                
///OnlineBtag (double-tag)

double sumdb=0e0;
for (int iJ=0; iJ<nSelJet; ++iJ) 
for (int iiJ=0; iiJ<nSelJet; ++iiJ) 
{
        if (iiJ==iJ) continue;
        int kkk=nSelJet-iJ - iiJ;
///     if (kkk<0) continue;

        double _wgt = bTagReleffOnline->eff(_theFlav[iJ],sbtag[iiovl].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
        _wgt *=  bTagReleffOnline->eff(_theFlav[iiJ],sbtag[iiovl].c_str(),ptJet[leadingJets[iiJ]],etaJet[leadingJets[iiJ]]);
        _wgt *= (1e0-bTagReleffOnline->eff(_theFlav[kkk],sbtag[iiovl].c_str(),ptJet[leadingJets[kkk]],etaJet[leadingJets[kkk]]));

        sumdb+=_wgt;    
}

	_onlineBtagWeight=sumdb/2e0;


///add pattern TTT
	double _T1 =  bTagReleffOnline->eff(_theFlav[0],sbtag[iiovl].c_str(),ptJet[leadingJets[0]],etaJet[leadingJets[0]]);
	double _T2 =  bTagReleffOnline->eff(_theFlav[1],sbtag[iiovl].c_str(),ptJet[leadingJets[1]],etaJet[leadingJets[1]]);
	double _T3 =  bTagReleffOnline->eff(_theFlav[2],sbtag[iiovl].c_str(),ptJet[leadingJets[2]],etaJet[leadingJets[2]]);

	_onlineBtagWeight += _T1*_T2*_T3;


} /// if _doMC


if (_isOk) _goodEvents++;




bool _exclude=false;
#ifdef EXCLUDE_TRAINER_FROM_SELECTION
	_exclude=true;
#endif



bool _passMVA=MVASelectAndFillEvent(_isOk,_exclude,_vecJets, _theFlav,weight*_onlineBtagWeight*_offlineBtagWeight*_SF_MC[iiovl]);

_passKinematicsMVA++;

if (!_passMVA) continue; /// mva selection, here
//cout<<"Event "<<iE<<" has been passed  "<< _lastComp <<" mva selection"<<endl;

_passKinematicsMVAAfter++;

///_MVA_END 

if (_isOk) _afterMVA++;




    bool tripleOnlineBtagMatchOK = true;
    bool tripleOnlineBtagMatchOK2=true;
    if (forceTripleOnlineBtagMatch) {
      tripleOnlineBtagMatchOK2 = false;
      if (nJet>=nSelJet ) {
	tripleOnlineBtagMatchOK2 = (isJetWithHltBtag[leadingJets[0]] && isJetWithHltBtag[leadingJets[1]])
	  || (isJetWithHltBtag[leadingJets[0]] && isJetWithHltBtag[leadingJets[2]])
	  || (isJetWithHltBtag[leadingJets[1]] && isJetWithHltBtag[leadingJets[2]]);
      }
    }

    // fill general jet histograms
    float deltaRj = -1;
    if ( (nJet>=nSelJet) && tripleOnlineBtagMatchOK ) {

      for (int iJ=0; iJ<nSelJet; ++iJ) {
	// jet kinematics
    if (tripleOnlineBtagMatchOK2) { ///cut-based analysis
	ptJetA[iJ]->fill(trgAccept,ptJet[leadingJets[iJ]],weight);
	etaJetA[iJ]->fill(trgAccept,etaJet[leadingJets[iJ]],weight);
	phiJetA[iJ]->fill(trgAccept,phiJet[leadingJets[iJ]],weight);
	isJetWithBtagA[iJ]->fill(trgAccept,isJetWithHltBtag[leadingJets[iJ]],weight);
}
      }
      float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
      if (dphij>3.1415926) dphij -= 2*3.1415926;
      if (dphij<-3.1415926) dphij += 2*3.1415926;
      dPhiJet1Jet2A->fill(trgAccept,dphij,weight);
      float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[0]];
      deltaRj = sqrt( dphij*dphij + detaj*detaj );
    }

    // determine jet flavor from MC info
    int theFlav[nSelJet];
    for (int iJ=0; iJ<nSelJet; ++iJ) {
      theFlav[iJ] = jetFlavorCode(leadingJets[iJ]);
    }

    // determine the dijet flavor code for the first two leading jets
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
	std::cout << "Bad 3rd flavor code " << theFcDijet << " " << theFc3rd << std::endl;
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
	std::cout << "Bad 3rd flavor code " << theFcDijet << " " << theFc3rd << std::endl;
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
      std::cout << "Bad flavor codes " << theFcDijet << " " << theFc3rd << std::endl;
    }
    if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
    hfctrip->Fill(theFcTrip,weight);

    // determine the svMass and svMassIndex
    float svMassIndex[nSelJet];
    for (int iJ=0; iJ<nSelJet; ++iJ) {
    if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
      svMassA[iJ]->fill(trgAccept,svMassJet[leadingJets[iJ]],weight);
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	// svMass under offline btag
	if (theBJetTag[ibtag][leadingJets[iJ]]>btcut[ibtag]) {
    if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
	  svMassBtA[ibtag][iJ]->fill(trgAccept,svMassJet[leadingJets[iJ]],weight);
	}
      }
      svMassIndex[iJ] = sv->getBinFromSvMass( svMassJet[leadingJets[iJ]] );
    }

    // compute invariant mass of two leading jets, and event btag
    if ( (nJet>=nSelJet) && (deltaRj>1) && tripleOnlineBtagMatchOK ) {
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
     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
	    mDijetBtA[ibtag]->fill(trgAccept,dijet_mass,weight);

     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
	    // triple-jet specific flavor code histograms
	    if (theFcTrip>=0) mDijetFcBtA[ibtag][theFcTrip]->fill(trgAccept,dijet_mass,weight);
	    // fc counting
     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
	    hfc[ibtag]->Fill( theFc,weight );
     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
	    if ( (dijet_mass > 100) && (dijet_mass < 140 ) ) hfcm[ibtag]->Fill( theFc,weight );
	    for (int iJ=0; iJ<nSelJet; ++iJ) {
     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
{
	      ptJetBtA[ibtag][iJ]->fill(trgAccept,ptJet[leadingJets[iJ]],weight);
	      if (theFcTrip>=0) ptJetFcBtA[ibtag][theFcTrip][iJ]->fill(trgAccept,ptJet[leadingJets[iJ]],weight);
}
	    }
	    // tpat counting
	    int theTpat = 4 * int(isJetWithHltBtag[leadingJets[2]]) + 2 * int(isJetWithHltBtag[leadingJets[1]])
	      + int(isJetWithHltBtag[leadingJets[0]]);
     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
	    atpattripall[ibtag]->fill(trgAccept,float(theTpat),weight);
	  }

	  if (_doMC) {
	    // now do the same without btag cut and with weighting
	    double wbt[nSelJet];
	    for (int iJ=0; iJ<nSelJet; ++iJ) {
	      wbt[iJ] = bTagEffOffline->eff(theFlav[iJ],sbtag[ibtag].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);
	    }
	    double wtot = weight * wbt[0] * wbt[1] * wbt[2]*_SF_MC[ibtag];
     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
{
	    mDijetBtwA[ibtag]->fill(trgAccept,dijet_mass,wtot);
	    // triple-jet specific flavor code histograms
	    if (theFcTrip>=0) mDijetFcBtwA[ibtag][theFcTrip]->fill(trgAccept,dijet_mass,wtot);
	    // fc counting
	    hfcww[ibtag]->Fill( theFc,wtot );
	    if ( (dijet_mass > 100) && (dijet_mass < 140 ) ) hfcmww[ibtag]->Fill( theFc,wtot );
	    for (int iJ=0; iJ<nSelJet; ++iJ) {
	      ptJetBtwA[ibtag][iJ]->fill(trgAccept,ptJet[leadingJets[iJ]],wtot);
	      if (theFcTrip>=0) ptJetFcBtwA[ibtag][theFcTrip][iJ]->fill(trgAccept,ptJet[leadingJets[iJ]],wtot);
	    }
	    // fill also the deltaPhi histograms
	    dPhiJet1Jet3BtwA[ibtag]->fill(trgAccept,dPhiJet1Jet3,wtot);
	    dPhiJet2Jet3BtwA[ibtag]->fill(trgAccept,dPhiJet2Jet3,wtot);
}
	  }
	}

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


     if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
 {
	  if ( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) &&  (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag])
	       && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {

	  if (ibtag==3)  _passTrig0++;
	    evBtagBtA[ibtag]->fill(trgAccept,nEvBTag,weight*_SF_MC[ibtag]);
	    massEvBtagBtA[ibtag]->fill(trgAccept,dijet_mass,float(nEvBTag),weight*_SF_MC[ibtag]);
		#if defined(RUNEVENTLUMI_LIST)
		if (ibtag==2 && trgAccept!=0) std::cout<<run<<" "<<event<<" "<<lumi<<std::endl;
		#endif


	  }
}
	}


        // in case of MC, compute trigger-weighted 2D histograms
        if (_doMC) {
          float thePt[nSelJet];
          float theEta[nSelJet];
          _lightIsOk=false;
          for (int iJ=0; iJ<nSelJet; ++iJ) {
            thePt[iJ] = ptJet[leadingJets[iJ]];
            theEta[iJ] = etaJet[leadingJets[iJ]];
            if (theFlav[iJ]==0) _lightIsOk=true;
          }

          for (int ibtag=0; ibtag<nbtag; ++ibtag) {
            if ( (theBJetTag[ibtag][leadingJets[0]]>btcut[ibtag]) 
                 &&  (theBJetTag[ibtag][leadingJets[1]]>btcut[ibtag])
                 && (theBJetTag[ibtag][leadingJets[2]]>btcut[ibtag]) ) {
              for (unsigned int ib=0; ib<genericTriggerList.size(); ++ib) {
                triggerWeight[ib] = 0;
              }

//  cout<<"SF_MC: "<< _SF_MC[ibtag]<<endl;

	if ( ibtag==2  && _lightIsOk ) _non3bb+=1;

              theHbbTrigWeight->getTrigWeight(bTagReleffOnline,sbtag[ibtag].c_str(),
                                              theFlav,thePt,theEta,triggerWeight);
              
              // fill with trigger weights
              massEvBtagBtTWA[ibtag]->fillTW(dijet_mass,float(nEvBTag),triggerWeight,weight*_SF_MC[ibtag]);
            }
          }
        }

    
 if (tripleOnlineBtagMatchOK2)  ///cut-based analysis
 {

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
	  int theTpat = 4 * int(isJetWithHltBtag[leadingJets[2]]) + 2 * int(isJetWithHltBtag[leadingJets[1]])
	    + int(isJetWithHltBtag[leadingJets[0]]);
	  for (int icateg=0; icateg<ncateg; ++icateg) {
	    if (lcateg[icateg]) tpatA[ibtag][icateg]->fill(trgAccept,float(theTpat),weight);
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

	  tpatAllA[ibtag]->fill(trgAccept,float(theTpat),weight);
	  // fill the templates
	  for (int ifc=0; ifc<nfc; ++ifc) {

	    for (int icateg=0; icateg<ncateg; ++icateg) {
	      if (thisTpat == icateg || thisTpat == 3) {  // online btag pattern must match offline
		if (lcateg[icateg]) {
		  float mistagWeight = bTagEffOffline->eff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]);
                  float errorMistagWeight = bTagEffOffline->erreff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]);


///SF 
///_SF_BEGIN_

double _SF=1.0;

#if defined(SF) || defined(SYST)

bool _dosf=false;
double _sfvariation=variation;

#if defined(SYST)  
if (systype==2) _dosf=true;
if (systype==3) _dosf=true;
#endif

if (systype==0) {_dosf=true; _sfvariation=0e0;}
if (systype==1) {_dosf=true; _sfvariation=0e0;}
if (systype==4) {_dosf=true; _sfvariation=0e0;}


if (_dosf) /// if it's SF reweighting or SF systematics
for (int iJ=0; iJ<nSelJet; ++iJ)
{

if (systype==0) {
#if defined(SF)
_sf_sys(1,_sfvariation,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics
#endif
}
else if (systype==4) {
#if defined(SF)
_sf_sys(1,_sfvariation,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics
#endif
}
else if (systype==1)
_sf_sys(systype,_sfvariation,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics
else if (systype==2) {
#if defined(SF)
        if (ifc==0) 
        _sf_sys(1,0e0,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics
#endif
      if (ifc>=1) 
     _sf_sys(systype,_sfvariation,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics

}

else if (systype==3) {
#if defined(SF)
        if (ifc>=1) 
        _sf_sys(1,0e0,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics
#endif
      if (ifc==0)
     _sf_sys(systype,_sfvariation,ibtag,ifc,ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]],_SF); ///possible SF systematics

}



}



#endif

///_SF_END_




		  // mass templates
		  for (int icorr=0; icorr<ncorr; ++icorr) {
		    double wbbpur = 1;
		    if (bbPurityCorr[ibtag] && (icorr == 1)) wbbpur =  fbbfrac[ibtag][icateg]->Eval( dijet_mass );
		    for (int itpat=0; itpat<3; ++itpat) {
		      double wtpat = 1;

		      if ( (itpat != icateg) && thisTpat != 3) {  // weight correction if offline & online btag pattern are different (except TTT)
			if (bTagReleffOnline->eff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]]) == 0) {       
			  std::cout << "call eff icateg=" << icateg << " ptJet=" << ptJet[leadingJets[icateg]] 
				    << " etaJet=" << etaJet[leadingJets[icateg]] << std::endl;
			}
			if (bTagReleffOnline->eff(2,sbtag[ibtag].c_str(),ptJet[leadingJets[itpat]],etaJet[leadingJets[itpat]]) == 0) {
			  std::cout << "call eff  itpat=" << itpat <<  " ptJet=" << ptJet[leadingJets[itpat]] 
				    << " etaJet=" << etaJet[leadingJets[itpat]] << std::endl;
			}
		      
			wtpat = bTagReleffOnline->eff(ifc,sbtag[ibtag].c_str(),ptJet[leadingJets[icateg]],etaJet[leadingJets[icateg]])
			  / bTagReleffOnline->eff(2,sbtag[ibtag].c_str(),ptJet[leadingJets[itpat]],etaJet[leadingJets[itpat]]);
		      }

		      massTemplateA[ifc][ibtag][icateg][icorr][itpat]->fill(trgAccept,dijet_mass,weight * mistagWeight * wbbpur * wtpat * _SF);

		      // event btag templates (no purity correction yet)
		      if (icorr == 1) {
			// get the svMassIndex probability vector for jet "categ", flavor=ifc
			float probSvMassIndexOfflineIfc[3];
			sv->getSVbins(ifc,ibtag,ptJet[leadingJets[icateg]],fabs(etaJet[leadingJets[icateg]]),probSvMassIndexOfflineIfc);
			float probSvMassIndexOnlineIfc[3];
			svOnline->getSVbins(ifc,ibtag,ptJet[leadingJets[icateg]],fabs(etaJet[leadingJets[icateg]]),probSvMassIndexOnlineIfc);
			    
			// same for jet "itpat", flavor=b
			float probSvMassIndexOfflineB[3];
			sv->getSVbins(2,ibtag,ptJet[leadingJets[itpat]],fabs(etaJet[leadingJets[itpat]]),probSvMassIndexOfflineB);
			float probSvMassIndexOnlineB[3];
			svOnline->getSVbins(2,ibtag,ptJet[leadingJets[itpat]],fabs(etaJet[leadingJets[itpat]]),probSvMassIndexOnlineB);

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

			  bTagTemplateA[ifc][ibtag][icateg][itpat]->fill(trgAccept,float(evtBTag),weight * mistagWeight * wbbpur * wtpat * wtpatEvtBtag * _SF);
			  massBTagTemplateA[ifc][ibtag][icateg][itpat]->fill(trgAccept,dijet_mass,float(evtBTag),weight * mistagWeight * wbbpur * wtpat * wtpatEvtBtag * _SF);
                          errorMassBTagTemplateA[ifc][ibtag][icateg][itpat]->fill(trgAccept,dijet_mass,float(evtBTag),weight * errorMistagWeight * wbbpur * wtpat * wtpatEvtBtag*_SF);

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
                          massBTagTemplateUncldA[ifc][ibtag][icateg][itpat]->fill(trgAccept,dijet_mass,float(evtBTag),weight * mistagWeight * wtpat * wtpatEvtBtag*_SF);
                          if (! (hasNegBTag[0] || hasNegBTag[1] )) {
                            massBTagTemplateCldA[ifc][ibtag][icateg][itpat]->fill(trgAccept,dijet_mass,float(evtBTag),weight * mistagWeight * wtpat * wtpatEvtBtag*_SF);
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
                            massBTagTemplateNonbbRA[ifc][ibtag][icateg][itpat]->fill(trgAccept,dijet_mass,float(evtBTag),weight * mistagWeight * wtpat * wtpatEvtBtag * weightR*_SF);
                          }
                        
       

			  // still need to fill template with combined weight
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }	  
	}

	if (_doMC) {
	  // here we fill the double btag purity histograms
	  // auxiliary array to perform permutations among first three jets
	  int shift[6] = { 0, 1, 2, 0, 1, 2 };
	
	  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	    for (int icateg=0; icateg<ncateg; ++icateg) {
	      int jetA = leadingJets[shift[icateg+1]];
	      int jetB = leadingJets[shift[icateg+2]];
	      int dibFcDijet = diJetFlavorCode( jetA, jetB );
	      if ( (theBJetTag[ibtag][jetA]>btcut[ibtag]) &&  (theBJetTag[ibtag][jetB]>btcut[ibtag]) ) {
		hmdibbt[dibFcDijet][ibtag][icateg]->Fill( dijet_mass, weight *  _SF_MC[ibtag]);
		// also fill the pt histograms under bouble btag selection
		hptdibbt[jetFlavorCode(jetA)][0][ibtag][icateg]->Fill( ptJet[jetA], weight  * _SF_MC[ibtag]);
		hptdibbt[jetFlavorCode(jetB)][1][ibtag][icateg]->Fill( ptJet[jetB], weight  *  _SF_MC[ibtag]);
	      }
	      // same with btag weighting. New: include online btag weight
	      hmdibww[dibFcDijet][ibtag][icateg]->Fill( dijet_mass, weight* _SF_MC[ibtag]
						       * bTagEffOffline->eff(theFlav[shift[icateg+1]],sbtag[ibtag].c_str(),ptJet[jetA],etaJet[jetA])
						       * bTagEffOffline->eff(theFlav[shift[icateg+2]],sbtag[ibtag].c_str(),ptJet[jetB],etaJet[jetB])
	      // 						     * bTagReleffOnline->eff(theFlav[shift[icateg+1]],sbtag[ibtag].c_str(),ptJet[jetA],etaJet[jetA])
	      // 						     * bTagReleffOnline->eff(theFlav[shift[icateg+2]],sbtag[ibtag].c_str(),ptJet[jetB],etaJet[jetB])  );
						       );
	    }
	  }
	}
  } // tripleOnlineBtagMatchOK2
      } else {
	std::cout << "Unphysical dijet_mass_sq = " << dijet_mass_sq << std::endl;
	//dijet_mass = -1;
      } // if (dijet_mass...)

    } // if (njet...)

_sf_sys.Reset();


  }

  // termination

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


  // data-driven subtraction of non-bb (fixed threshold method). Only for tpat-corrected hist.
  int mtpat = 3;
  massBTagTemplateCldRA[ifc][ibtag][icateg][mtpat]->add( massBTagTemplateUncldA[ifc][ibtag][icateg][mtpat],
                                                               massBTagTemplateNonbbRA[ifc][ibtag][icateg][mtpat],
                                                               1, -1);


      }
    }


  } 


///prepare  _PU_nentries for writting
 _PU_nentries/=_eff_presel;
_countPU->SetBinContent(1, _PU_nentries);



  hout->Write();
  hout->Close();


///_MVA_BEGIN_ clean up memory and calculation of events went trought mva computers

MVAClean();


cout<<"Number of events passed kinematics preselection: MC && (+trigger, btagging cuts):DATA  ==>" <<_goodEvents<<endl;
cout<<"Number of events passed kinematics preselection & MVA: MC && (+trigger, btagging cuts):DATA  ==>" <<_afterMVA<<endl;

cout<<"Number of events passed selection non 3 b-jets "<<_non3bb<<endl;

cout<<"Number of events  before Syst"<<_beforeSys<<endl;
cout<<"Number of events  before Jet"<<_beforeJet<<endl;
cout<<"Number of events  after Jet"<<_afterJet<<endl;
cout<<"Number of events  passed Kin.Off.Selection"<<_passKinematics<<endl;
cout<<"Number of events  passed Kin.Off.Selection+MVA"<<_passKinematicsMVA<<endl;
cout<<"Number of events  passed Kin.Off.Selection+MVA+AFter"<<_passKinematicsMVAAfter<<endl;
cout<<"Number of events  passed CSVT to trigger 0 "<<_passTrig0<<endl;



///_MVA_END_


///_PU_BEGIN_ clean up memmory

#if defined(PU)

if (_pileup) delete _pileup;

#endif
///_PU_END_ 

}

///To test only: create binary


#if defined(CMSSW_CPP)

int main(int argc, char * argv[])
{

  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();



  // initialize command line parser
optutl::CommandLineParser parser ("mva parser");

parser.addOption("mvasettings",optutl::VariableMapCont::kString,"file having mvasettings");

  // set defaults
  parser.stringValue ("mvasettings") =  "/Analysis/HbbMSSMAnalysis/data/mvasettings.lst";

  // parse arguments
  parser.parseArguments (argc, argv);

  std::string mvaset_ = parser.stringValue("mvasettings");

TString _fileIn;
TString __str=gSystem->GetFromPipe("echo $CMSSW_BASE");
__str+="/src/";
_fileIn = __str+ TString(mvaset_.c_str());

  cout<<"file full path: "<<_fileIn<<endl;
  TString _mvasettings = _fileIn;


TripleBtagAnalysis_MVA_SF(_mvasettings);
}

#endif

