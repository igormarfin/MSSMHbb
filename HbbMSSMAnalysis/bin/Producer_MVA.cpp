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

#include "Analysis/Utilities/interface/svMass2012.h"


#include "Analysis/Utilities/interface/checkHLTBTagMatch.h"
#include "Analysis/Utilities/interface/TriggerRunSelection.h"



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


///MVA selection tree support
#include "Analysis/Utilities/interface/HbbSelect.h"


///_MVA_END_




///_PU_END_


///_SF_BEGIN_

#if defined(SF)  || defined(SYST) || defined(JER)

#include "Analysis/Utilities/interface/HbbSystControl.h"
#include "Analysis/HbbMSSMAnalysis/interface/SystReadInput.h"

#endif
///_SF_END_





TCanvas* canvas;




void Producer_MVA(TString _mvasettings="Analysis/HbbMSSMAnalysis/data/mvasettings.lst") {

#if !(defined(CMSSW_CPP))

TString __base=gSystem->GetFromPipe("echo $CMSSW_BASE");
__base+="/src/";
_mvasettings=__base+_mvasettings;


#endif

  // macro for analysis in channel with three btagged jets

  canvas = new TCanvas ("cg1","mycanvas",10,10,800,600);
  // open an ntuple file

  bool _doMC = false;  // set to true for MC, false for real data

//  bool doLumiWeighting = true;  // weight MC events according to intgrated lumi



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

  std::vector<std::string> HLTBtagTriggerObjectList;
  std::vector<std::string> HLTBtagTriggerObjectFilterList;

std::vector<std::string> L25BtagTriggerObjectFilterList;
std::vector<std::string> L25BtagTriggerObjectList;


  if (_doMC) {


	/// 2012 analysis

	 triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v");
	 triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV_v");
	 triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose_v");

    //there are two HLT Btag trigger objects. The first two triggers specified above use the same HLT Btag objects
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");


//Now the same for L25 objects (first element corresponds to first element of triggerFilterList and so on)
L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhi1stTrackL1FastJetFastPV");
	

  } else {

    /// 2012 analysis



     triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v");
     triggerFilterList.push_back("HLT_Jet80Eta1p7_Jet70Eta1p7_DiBTagIP3DFastPV_v");
     triggerFilterList.push_back("HLT_Jet160Eta2p4_Jet120Eta2p4_DiBTagIP3DFastPVLoose_v");


    //there are two HLT Btag trigger objects. The first two triggers specified above use the same HLT Btag objects
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
    HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");



//Now the same for L25 objects (first element corresponds to first element of triggerFilterList and so on)
L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhiL1FastJetFastPV");
L25BtagTriggerObjectFilterList.push_back("hltBLifetimeL25FilterbbPhi1stTrackL1FastJetFastPV");


  }

  // extract generic trigger list and number of input events from ntuple files
  float lumiScaleFac = 1;
  float xsect=0e0;
  //int nInputEvents= 
  getHbbMetadata(theFileColl,genericTriggerList,HLTBtagTriggerObjectList,1000.,lumiScaleFac,_doMC,xsect);
//  getHbbMetadata(theFileColl,genericTriggerList,1000.,lumiScaleFac,_doMC);
getHbbL25Metadata(theFileColl,L25BtagTriggerObjectList);

  // this trigger weight array will be used later in each event
//  float* triggerWeight = new float[genericTriggerList.size()];
  // create trigger weight helper object
//  HbbTrigWeight* theHbbTrigWeight = new HbbTrigWeight(&genericTriggerList,&triggerFilterList);


  f.AddFileInfoList((TCollection*) theFileColl.GetList());
  TTree* hbbtree = &f;


  //kinematic analysis cuts
  double maxEta = 1.7;
  double deltaEtaCut12 = 1.7;//1.7;

  double jetPtMin[nSelJet] = { 60.0, 53.0, 20};  // leading jet lower pt thresholds
  std::string Scenario("MediumMass2012");  // Higgs mass scenario // available: LowMass2012, MediumMass2012, HighMass2012


#if defined (PU)
	TString PU_mc="MC_Summer12_PU_S10-600bins.root";
	TString PU_data="MyDataPileupHistogram_Jet60Eta1p7_Jet53Eta1p7_600bins_true_190389-208686.root";

#endif

#if defined(HIGH_MASS)
  Scenario= std::string("HighMass2012");
  jetPtMin[0] = 80;
  jetPtMin[1] = 70.0;

  #if defined (PU)
		 PU_data="MyDataPileupHistogram_Jet80Eta1p7_Jet70Eta1p7_600bins_true_190389-208686.root";
  #endif
#endif

#if defined(VERY_HIGH_MASS)
    Scenario= std::string("VeryHighMass2012");
    jetPtMin[0] = 160.0;
    jetPtMin[1] = 120.0;
    maxEta = 2.4;

  #if defined (PU)
         PU_data="MyDataPileupHistogram_Jet160Eta2p4_Jet120Eta2p4_600bins_true_190389-208686.root";
  #endif

#endif

  double jetPtMax[nSelJet] = {3500, 3500, 3500}; // leading jet upper pt thresholds

//set it to false for data
 bool forceTripleOnlineBtagMatch =false;  // require triple online btag match on triple btag signature 
 bool doTrigSelect = true;  // perform trigger selection according to scenario


#if defined(ONLINE_BTAG_MATCH)
forceTripleOnlineBtagMatch =true;
#endif





  // create trigger run selection object
  TriggerRunSelection trigSelector( Scenario, &genericTriggerList );

   TriggerRunSelection * ptrTrigSelector = & trigSelector;

 if (_doMC) {
	doTrigSelect = false;
	ptrTrigSelector=NULL;
	}

  checkHLTBtagMatch HLTBtagMatchObj(
                                    triggerFilterList,
                                    HLTBtagTriggerObjectFilterList,
                                    HLTBtagTriggerObjectList,
                                    doTrigSelect,
                                    ptrTrigSelector
                                    );


//add the information for the L25 objects
HLTBtagMatchObj.SetL25(L25BtagTriggerObjectFilterList,L25BtagTriggerObjectList);


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
  hbbtree->SetBranchStatus("*L25*",1);

  //end of switch off TTree branches that are currently not needed?




  // here we define dimensions and labels
  const int nflav=5;  // number of flavor groups to be distinguished
  const int nbtag = 4;  // number of offline btag working points
//  const int ncateg = 3; // ranks of untagged jet
//  const int ncorr = 2; // correction levels
//  const int ntpat = 4; // trigger pattern (online btag). (0,1,2)=rank of non-btagged jet, 3=all combined

 std::string onlineBtagVersion("V4"); // choices: smpf or V4




  const std::string flavlabel[nflav] = {"udsg","c","b","cc","bb"};  // labels for flavor groups
  // define names of btag variants
  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };
  // corresponding discriminant cuts
  double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags
//  bool bbPurityCorr[nbtag] = { true, true, true, true };  // specify if bbPurity collection should be applied

  // define flavor pair classes
//  const int nfcDijet = 6;
//  string sfcDijet[nfcDijet] = { "bb", "bc", "bq", "cc", "cq", "qq" };
//  const int nfc3rd = 3;
//  string sfc3rd[nfc3rd] = { "q", "c", "b" };
  // flavor triplet classes
//  const int nfcTrip = 6;
//  string sfcTrip[nfcTrip] = { "bbb", "bbc", "bbq", "bcb", "bqb", "non-bb" };
  // the following are used for templates
//  const int nfc = 3;  
//  string sfc[nfc] = { "q", "c", "b" };

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

bTagEff* bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_mediummass.root","offline",flavlabel,sbtag,5,nbtag);
bTagEff* bTagReleffOnline =  new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_normal_mediummass.root","online",flavlabel,sbtag,5,nbtag);

  // open the SVMass template files
TString SVMassFileName("/data/user/rasp/Run/Hbb/svMass.root");
TString SVMassOnlineFileName("/data/user/rasp/Run/Hbb/svMass1BL1FastJetFastPVHltMatch.root");


  if(Scenario == "HighMass2012") {
    //FIXME: use the same relative online btag efficiencies as for the 'normal' triggers -> has to be changed!
  bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_highmass.root","offline",flavlabel,sbtag,5,nbtag);
  bTagReleffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_normal_highmass.root","online",flavlabel,sbtag,5,nbtag);
  SVMassOnlineFileName="/data/user/rasp/Run/Hbb/svMass1BL1FastJetFastPVHltMatch.root" ;
  }


  if(Scenario == "VeryHighMass2012") {
    //FIXME: use the same relative online btag efficiencies as for the 'normal' triggers -> has to be changed!
  bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_veryhighmass.root","offline",flavlabel,sbtag,5,nbtag);
  bTagReleffOnline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagMatrixOnline-V61_53X_loose_veryhighmass.root","online",flavlabel,sbtag,5,nbtag);
  SVMassOnlineFileName="/data/user/rasp/Run/Hbb/svMass1BLooseL1FastJetFastPVHltMatch.root";
  }


///2012 
//  svMass2012 * sv = new svMass2012(SVMassFileName);
//  svMass2012 * svOnline = new svMass2012(SVMassOnlineFileName);



  // open the bbPurity correction functions
/*
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
*/
  // open the SVMass template files

//  TString SVMassFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_svMass.root");
//  TString SVMassOnlineFileName("/afs/naf.desy.de/user/r/rasp/public/Hbb4b/TTJets_BTagHLTmatched_svMass.root");
//  svMass * sv = new svMass(SVMassFileName);
//  svMass * svOnline = new svMass(SVMassOnlineFileName);


  // create the Root output file
  TFile* hout = new TFile("SelTree.root","recreate");


  TH1::AddDirectory(true);
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();  // not really needed (we inherit from TH1)

  Int_t nentries = (Int_t) hbbtree->GetEntries();
  cout << "Number of events in ntuple: " << nentries << endl;

#include "Analysis/Utilities/interface/HbbNtuple.cc"


TTree *seltree = new TTree("selTree","selTree");
#include "Analysis/Utilities/interface/HbbSelectBranch.cc"


ReadMVASettings( _mvasettings);

int _goodEvents=0;
//int _afterMVA=0;
///_MVA_END_ ini





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




double _passKinematics=0;
double _passKinematicsMVA=0;


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



    float deltaR12 = -1.0f;
    float deltaR23 = -1.0f;
    float deltaR13 = -1.0f;

    float deltaEta = 9999.0; //between leading jets



std::vector<int> leadingJets;

// cout<<"Event : "<<iE<<endl;

#if !defined(RUNEVENTLUMI_LIST)
if (iE%5000==0) { cout<<"Event : "<<iE; 	cout << '\xd';}
// cout<<"Event : "<<iE<<endl;
#endif




///_MVA_BEGIN_ declaration of objects to feed MVA!


vector<math::XYZTLorentzVector> _vecJets ; /// objects to feed mvacomputer and mvatrainer


//vector<double> * res=0; /// response of mva computer
double minPtMVA=5e0; /// some cut on object
double maxEtaMVA = 5.2;
std::vector<int> _theFlav; /// to catch flavor of objects
std::vector<int> _leading; /// to store jet index for MVA
bool _isplaced=false;
double _offlineBtagWeight=1e0; /// needed by MVATrainer
double _onlineBtagWeight=1e0; /// needed by MVATrainer


///_MVA_END_



///_SF_BEGIN_ 


double _SF_MC[nbtag]={1.0,1.0,1.0,1.0};





 // trigger selection

    if (doTrigSelect) {
//    cout<<"Trigger before: "<<trgAccept<<endl;
    unsigned int trgSelect = trgAccept & trigSelector.mask( run );
//     cout<<"Mask "<<trigSelector.mask( run )<<endl;
	 trgAccept=trgSelect;
//	cout<<"Trigger after: "<<trgAccept<<endl;

     }

 
    // find set of leading jets
    int nJet = 0;
    // loop over the jets
    for (int iJet=0; iJet<numberOfJets; ++iJet) {

#if defined(MVA_WITH_KINSELECTION)
    _isplaced=false;

      if (nJet<nSelJet)
      if ( (fabs(etaJet[iJet])<maxEta) ) 
      if ( (numberOfConstituentsInJet[iJet] > 1) )      
      if(puJetIDLoose[iJet])
      if(jetIDLoose[iJet] )

      if ( (ptJet[iJet] > jetPtMin[nJet]) && (ptJet[iJet] < jetPtMax[nJet]) ) {
	
//	cout<<"Jet["<<nJet<<"] :"<< "Pt="<<ptJet[iJet]<<"  Px= "<<pxJet[iJet]<<endl;

	leadingJets.push_back(iJet);
	++nJet;

///_MVA_BEGIN_ getting of the first nSelJet-1 th objects' info 

 _theFlav.push_back(jetFlavorCodeNewCat(iJet));
 _leading.push_back(iJet);
	  _vecJets.push_back(math::XYZTLorentzVector(pxJet[iJet], pyJet[iJet], pzJet[iJet], energyJet[iJet]));
	  _isplaced=true;

//_MVA_END_ 	
      }

///_MVA_BEGIN_ getting of  nSelJet th object info

 	      if (nJet == nSelJet && !_isplaced)
           if(puJetIDLoose[iJet])
          if(jetIDLoose[iJet] )
		  if ( (numberOfConstituentsInJet[iJet] > 1) )

#endif
		  if ( (fabs(etaJet[iJet])<maxEtaMVA) )
		   if ( ptJet[iJet]  > minPtMVA ) {

#if !defined(MVA_WITH_KINSELECTION)
 		leadingJets.push_back(iJet);
#endif		
         _theFlav.push_back(jetFlavorCodeNewCat(iJet));
		 _leading.push_back(iJet);
		 _vecJets.push_back(math::XYZTLorentzVector(pxJet[iJet], pyJet[iJet], pzJet[iJet], energyJet[iJet]));

		}	


	
//_MVA_END_ 


    } /// end of loop over jets









selRun=run;
selEvent= event;
selPass= 0;
mjj=-1e13;

#if defined(MVA_WITH_KINSELECTION)


if (nJet < nSelJet) {

seltree->Fill();

continue;
}

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

if ( !(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) )
{
seltree->Fill();
continue;
}




#endif 

if ((int) _vecJets.size()<nSelJet) 
{

#if !defined(MVA_WITH_KINSELECTION) &&  !defined(MVA_WITH_TRIGSELECTION)
if (_vecJets.size()>1) mjj=(_vecJets[0]+_vecJets[1]).M();
else mjj=-1e13;
#else
mjj=-1e13;

#endif


seltree->Fill();
continue;

}




///End of preselection?
_passKinematics++;



  
      // check the matching flags
/*
        int mJet = 0;
        int matchPat = 0;
        for (int iJet=0; iJet<numberOfJets; ++iJet) {
          if (! (ptJet[iJet] > 15) ) continue;
          if (! (fabs(etaJet[iJet])<maxEta) ) continue;
          if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
          if (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[iJet],run)) {
            matchPat = matchPat | (1<<mJet);
          }
          ++mJet;
        }
  
*/

    // compute the weight for this event
    float weight = 1;
//    if ( _doMC && doLumiWeighting )  weight = lumiScaleFac;






///_MVA_BEGIN_ mva selection



//  const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
//  string btdiscr[nbtag] = { "TCHP", "TCHP", "CSV", "SSVHP" };
  // corresponding discriminant cuts
//  double btcut[nbtag] = { 3.41, 6 , 0.898, 2 };  // working points for btags
//  bool bbPurityCorr[nbtag] = { true, true, true, true };  // specify if bbPurity collection should be applied


///We are working with CSVT
   	int iiovl=2;

/*
 float _dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
 if (_dphij>3.1415926) _dphij -= 2*3.1415926;
 if (_dphij<-3.1415926) _dphij += 2*3.1415926;
 float _detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[0]];
 float  _deltaR = sqrt( _dphij*_dphij + _detaj*_detaj );
*/


bool _isOk=false;


///data
if (!_doMC) {


bool _btag=false;
bool _tripleOnlineBtagMatchOK = true;

if (forceTripleOnlineBtagMatch) {
      _tripleOnlineBtagMatchOK = false;
      if (nJet>=nSelJet ) {
    _tripleOnlineBtagMatchOK = (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) &&
     HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run))
      || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[0]],isJetWithL25JetBitPattern[leadingJets[0]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run))
      || (HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[1]],isJetWithL25JetBitPattern[leadingJets[1]],run) && HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[leadingJets[2]],isJetWithL25JetBitPattern[leadingJets[2]],run));
      }
    }




if ( (theBJetTag[iiovl][leadingJets[0]]>btcut[iiovl]) && (theBJetTag[iiovl][leadingJets[1]]>btcut[iiovl]) &&  (theBJetTag[iiovl][leadingJets[2]]>btcut[iiovl]))
_btag=true;

if (_tripleOnlineBtagMatchOK)
if ( (deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) )
if (  _btag )  
if (trgAccept)
//if (trgAccept  & (1<<0) || trgAccept  & (1<<1) || trgAccept  & (1<<2)) ///triggered data by Trigger0 or Trigger1 or Trigger2
//if (trgAccept  & (1<<9) || trgAccept  & (1<<14) || trgAccept  & (1<<19) ) ///triggered data by Trigger9 or Trigger14 or Trigger19
_isOk=true;


}





///MC
if (_doMC) {
if ( (deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) )
 _isOk=true;

///SF 
///_SF_BEGIN_






///OfflineBtag (triple-tag) weigting

          for (int iJ=0; iJ<nSelJet; ++iJ) 
          _offlineBtagWeight *=   bTagEffOffline->eff(_theFlav[iJ],sbtag[iiovl].c_str(),ptJet[leadingJets[iJ]],etaJet[leadingJets[iJ]]);        
                
///OnlineBtag (double-tag)



#if !defined(ONLINE_BTAG_OFF)

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

#endif 

} /// if _doMC


if (_isOk) _goodEvents++;






bool _exclude=false;
#ifdef EXCLUDE_TRAINER_FROM_SELECTION
	_exclude=true;
#endif


bool _passMVA=MVASelectAndFillEvent(_isOk,_exclude,_vecJets, _theFlav,_leading,weight*_onlineBtagWeight*_offlineBtagWeight*_SF_MC[iiovl]);


#if defined(MVA_WITH_TRIGSELECTION)

selPass=(_passMVA && _isOk ) ? 1 : 0;

#else
selPass=(_passMVA) ? 1 : 0;

#endif



#if defined(MVA_WITH_TRIGSELECTION)
if (_vecJets.size()>1 && _isOk) mjj=(_vecJets[0]+_vecJets[1]).M();
else mjj=-1e13;
#else
if (_vecJets.size()>1) mjj=(_vecJets[0]+_vecJets[1]).M();
else mjj=-1e13;
#endif



seltree->Fill();




if (!_passMVA) continue; /// mva selection, here

//cout<<"Event "<<iE<<" has been passed  "<< _lastComp <<" mva selection"<<endl;


_passKinematicsMVA++;


///_MVA_END 






} //nevt 


seltree->Write();


  hout->Write();
  hout->Close();


///_MVA_BEGIN_ clean up memory and calculation of events went trought mva computers

MVAClean();



cout<<"Number of events  passed Kin.Off.Selection"<<_passKinematics<<endl;
cout<<"Number of events  passed Kin.Off.Selection+MVA"<<_passKinematicsMVA<<endl;



///_MVA_END_



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


Producer_MVA(_mvasettings);
}

#endif

