// -*- C++ -*-
//
// Package:    HBBAnalysis
// Class:      HBBAnalysis
// 
/**\class HBBAnalysis HBBAnalysis.cc Jets/HBBAnalysis/src/HBBAnalysis.cc
 
 Description: [one line class summary]
 
 Implementation:
 [Notes on implementation]
 */
//
// Original Author:  Aliaksei Raspiareza,,,DESY
//         Created:  Fri Apr 15 10:25:07 CEST 2011
// $Id: HbbMSSMAnalysis.cc,v 1.5 2011/08/26 15:17:18 rmankel Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <DataFormats/Candidate/interface/Candidate.h>


#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticleFwd.h"

#include "CommonTools/Utils/interface/PtComparator.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "TTree.h"

#include "TVector3.h"

//
// class declaration
//

class HbbMSSMAnalysis : public edm::EDAnalyzer {
public:
    explicit HbbMSSMAnalysis(const edm::ParameterSet&);
    ~HbbMSSMAnalysis();
    
    
private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;
    bool IsJetMatchedToHLTBtag( const pat::JetRef & jet, std::vector<reco::Particle> HLTBtag , const double DR,  
                               float & theDeltaRWithHltBtag, float & theDeltaPtWithHltBtag );
    int jetHflContent( const pat::JetRef & jet, const edm::Handle<reco::GenParticleCollection> & theGenParticles );
    
    
    bool IsJetMatchedToL1Jet( const pat::JetRef & , const l1extra::L1JetParticleCollection & , const double ,  
                             float & , float &  );
    
    // ----------member data ---------------------------
    
    
    edm::InputTag _triggerEvent;
    edm::InputTag _jetCollection;
    edm::InputTag _muonCollection;
    edm::InputTag _primaryVtxCollection;
    edm::InputTag _genParticleCollection; 
    edm::InputTag _metCollection;
    edm::InputTag _triggerPathSrc;
    edm::InputTag _triggerEventTag;
    std::vector<edm::InputTag> _l1ExtraJets;
    edm::InputTag _genJetCollection;
    edm::InputTag _l1Muons;
    
  //std::string _hltBTagJetFilterName;
  std::vector<std::string> _hltBTagJetFilterName;

  //std::string _hltL2JetFilterName;
   std::vector<std::string> _hltL2JetFilterName;
    std::string _hltMuonFilterName;
    
    
    
    // Cuts --->
    float _jetPtCut;
    float _muonPtCut;
    int _minNumberOfJets;
    
    // number of processed events -->
    int _events;
    // number of selected events  -->
    int _selectedEvents;
    
    bool _doMC;
    bool _printOut;
    bool _doNegativeBtags;
    
    // Tree variables
    
    // run, event, lumi
    int run;
    int event;
    int lumi;
    
    
    // Monte Carlo information -->
    
    // Higgs Boson ;
    // PDG code of neutral MSSM Higgs ( PDG(h)=25, PDG(H)=35, PDG(A)=36 ) 
    int pdgHiggsBoson; 
    float etaHiggsBoson;
    float phiHiggsBoson;
    float ptHiggsBoson;
    float massHiggsBoson;
    float pxHiggsBoson;
    float pyHiggsBoson;
    float pzHiggsBoson;
    float energyHiggsBoson;
    
    // partons (quarks) from Higgs boson decay 
    int numberOfPartons; // should be two : b b~ ... but nonetheless... keep option for larger number of partons :)
    int pdgParton[10];
    float pxParton[10];
    float pyParton[10];
    float pzParton[10];
    float energyParton[10];
    float etaParton[10];
    float phiParton[10];
    float ptParton[10];
    
    // ** generator info (x-sections)
    float _internalXsec;
    float _internalXsecE;
    float _externalXsecLO;
    float _externalXsecLOE;
    float _externalXsecNLO;
    float _externalXsecNLOE;
    float _filterEfficiency;
    
    // pileup information from MC
    int nPU;
    int nPUInTime;
    
    // bit pattern for triggers defining PD
    int trgAccept;
    
    // reconstructed primary vertices -->
    int nPV; // number of primary vertices
    float probPV[100]; // prob associated with chi2 of vertex fit
    float ndofPV[100]; // ndof of vertex fit
    int   nTrkPV[100]; // number of tracks in primary vertex
    float chi2PV[100]; // chi2 of vertex fit
    float xPV[100]; // x coordinate of primary vertex
    float yPV[100]; // y coordinate of primary vertex
    float zPV[100]; // z coordinate of primary vertex
    float sumPtPV[100]; // scalar sum of tracks' pt
    
    // reconstructed jets -->  
    // general information
    int numberOfJets;
    float etaJet[100];
    float phiJet[100];
    float ptJet[100];
    float pxJet[100];
    float pyJet[100];
    float pzJet[100];
    float energyJet[100];
    int numberOfTracksJet[100];
    float trackPxJet[100][3];
    float trackPyJet[100][3];
    float trackPzJet[100][3];
    
    int numberOfConstituentsInJet[100];
    int numberOfChargedConstituentsInJet[100];
    float neutralHadronEnergyFraction[100];
    float photonEnergyFraction[100];
    float electronEnergyFraction[100];
    float chargedHadronEnergyFraction[100];
    
    bool jetIDLoose[100];
    bool jetIDMedium[100];
    bool jetIDTight[100];
    
    int nSvTracksJet[100];
    float svMassJet[100];
    float svFDsigJet[100];
    
    // gen jet
    int genNumberOfJets;
    float genEtaJet[100];
    float genPhiJet[100];
    float genPtJet[100];
    float genPxJet[100];
    float genPyJet[100];
    float genPzJet[100];
    float genEnergyJet[100];
    
    // jet - L1 Jet
    int l1NumberOfJets;
    float l1EtaJet[100];
    float l1PhiJet[100];
    float l1PtJet[100];
    float l1PxJet[100];
    float l1PyJet[100];
    float l1PzJet[100];
    float l1EnergyJet[100];
    
    // jet - L2 Jet
    int l2NumberOfJets;
    float l2EtaJet[100];
    float l2PhiJet[100];
    float l2PtJet[100];
    float l2PxJet[100];
    float l2PyJet[100];
    float l2PzJet[100];
    float l2EnergyJet[100];
    
    // jet - HLT BTag Jet
    int l3NumberOfBJets;
    float l3EtaBJet[100];
    float l3PhiBJet[100];
    float l3PtBJet[100];
    float l3PxBJet[100];
    float l3PyBJet[100];
    float l3PzBJet[100];
    float l3EnergyBJet[100];
    
    
    // jet - parton match 
    bool isJetMatchedPartonExist[100];
    bool isJetMatchedPartonFromHiggsBosonDecay[100]; 
    float pxJetMatchedParton[100];
    float pyJetMatchedParton[100];
    float pzJetMatchedParton[100];
    float energyJetMatchedParton[100];
    int flavorJetMatchedParton[100];
    int partonFlavorJet[100];         // jet->partonFlavour()
    int hflContentJet[100];           // from deltaR matching of bc partons
    // BTag information ---->
    float jetBProbBJetTag[100]; // track combined B-probability tag
    float jetProbBJetTag[100];  // track combined probability tag
    float tcHPBJetTag[100]; // track counting high purity b-tag
    float tcHEBJetTag[100]; // track counting high efficiency b-tag
    float ntcHPBJetTag[100]; // negative track counting high purity b-tag
    float ntcHEBJetTag[100]; // negative track counting high efficiency b-tag
    float svHEBJetTag[100]; // simple secondary vertex high efficiency b-tag
    float svHPBJetTag[100]; // simple secondary vertex high purity b-tag
    float combSVBJetTag[100]; // combined secondary vertex tag
    float combSVMVABJetTag[100]; // combined MVA secondary vertex tag
    float tcMinus2ndBJetTag[100]; // ipsig3d of second track from negative end
    float tcMinus3rdBJetTag[100]; // ipsig3d of third  track from negative end
  int tcNumberOfSelectedTracks[100]; // number of selected tracks for TC btag
    
  // HLT btag match information
  //bool isJetWithHltBtag[100]; // true if matched to HLT btag object
  unsigned int isJetWithHltBtagBitPattern[100]; //true if matched to HLT btag object
  float deltaRWithHltBtag[100]; //deltaR for btag match
  float deltaPtWithHltBtag[100];// pt (jet) - pt(btag match)
  


    // L1 match information
    bool isJetWithL1Jet[100]; // true if matched to L1 object
    float deltaRWithL1Jet[100]; // deltaR for L1 match
    float deltaPtWithL1Jet[100]; // pt (jet) - pt(L2 match)
    
    // L2 match information
    //bool isJetWithL2Jet[100]; // true if matched to L2 object
  unsigned int isJetWithL2JetBitPattern[100];
    float deltaRWithL2Jet[100]; // deltaR for L1 match
    float deltaPtWithL2Jet[100]; // pt (jet) - pt(L2 match)
    
    // reconstructed muons -->
    int numberOfMuons;
    int chargeMuon[20];
    bool isTrackerMuon[20];
    bool isGlobalMuon[20];
    // absolute muon PFlow isolation variables -->
    float pfChargedHadronIsoMuon[20]; 
    float pfNeutralHadronIsoMuon[20];
    float pfGammaIsoMuon[20];
    bool isMuonInJet[20]; // whether muon is assigned to jet;
    int associatedJetIndex[20]; // index of jet containing muon (index in the created array of jets)
    float etaMuon[20];
    float phiMuon[20];
    float ptMuon[20];
    float pxMuon[20];
    float pyMuon[20];
    float pzMuon[20];
    // matching generator muon -->
    bool isMatchedGeneratorMuonExist[20];
    int chargeGenMuon[20];
    float pxGenMuon[20];
    float pyGenMuon[20];
    float pzGenMuon[20];
    int pdgMotherGenMuon[20];
 
  // L1 muons
  int l1NumberOfMuons;
  float l1EtaMuon[100];
  float l1PhiMuon[100];
  float l1PtMuon[100];
  float l1PxMuon[100];
  float l1PyMuon[100];
  float l1PzMuon[100];
  float l1EnergyMuon[100];
  int l1ChargeMuon[100];

  // HLT Muons
  int hltNumberOfMuons;
  float hltEtaMuon[100];
  float hltPhiMuon[100];
  float hltPtMuon[100];
  float hltPxMuon[100];
  float hltPyMuon[100];
  float hltPzMuon[100];
  float hltEnergyMuon[100];
  int hltChargeMuon[100];

    
       
    // Missing transverse momentum;
    float met;
    float missingPx;
    float missingPy;
    
    // trigger prescale factors
    float trigPrescale[25];
    
    // TFile Service
    edm::Service<TFileService> fs;
    
    // Trees
    TTree * tree;
    
    // Tree with info on cross-sections
    TTree * treeGenInfo;
    
    // histogram with trigger dictionary
    TH1F* gtlHist;
   // histogram with hltBTagJetFilterName dictionary
  TH1F* tlHist_hltBTagJetFilterName;
  TH1F* tlHist_hltL2JetFilterName;

    bool firstEvent;
    
    std::vector<std::string> _genericTriggerList;
    
    // trigger counters
    std::map<std::string, int> triggerCounter;
    std::map<std::string, int> triggerObjectCounter;
    
};

//
// constructors and destructor
//
HbbMSSMAnalysis::HbbMSSMAnalysis(const edm::ParameterSet& iConfig)

{
    //now do what ever initialization is needed
    
    // collection names --->
    _muonCollection        = iConfig.getParameter<edm::InputTag>                ( "MuonSource" );
    _metCollection         = iConfig.getParameter<edm::InputTag>                ( "METSource" );
    _jetCollection         = iConfig.getParameter<edm::InputTag>                ( "JetSource" );
    _primaryVtxCollection  = iConfig.getParameter<edm::InputTag>                ( "PrimaryVertexSource" ); 
    _genParticleCollection = iConfig.getParameter<edm::InputTag>                ( "GenParticleSource" );
    _triggerEvent          = iConfig.getParameter<edm::InputTag>                ( "triggerEvent" );
    _triggerPathSrc        = iConfig.getParameter<edm::InputTag>                ( "triggerPathCollection" );
    _l1ExtraJets           = iConfig.getParameter< std::vector<edm::InputTag> > ( "L1ExtraJets" );  // could be tau, central or forward
    _genJetCollection      = iConfig.getParameter<edm::InputTag>                ( "GenJetSource" );
  _l1Muons               = iConfig.getParameter<edm::InputTag>                ( "L1Muons" );
    
    
    // trigger event tag (first declare default)
    const edm::InputTag dTriggerEventTag("hltTriggerSummaryAOD","","HLT");
    _triggerEventTag = iConfig.getParameter<edm::InputTag>("triggerEventTag");


    //_hltBTagJetFilterName = iConfig.getParameter<std::string>("hltBTagJetFilterName");
    _hltBTagJetFilterName = iConfig.getParameter< std::vector<std::string> > ( "hltBTagJetFilterName");
 // package the list of trigger names into a histogram, so that it is included in the ntuple as an object that can be cleanly merged
    tlHist_hltBTagJetFilterName = fs->make<TH1F>("tlHist_hltBTagJetFilterName","tlHist_hltBTagJetFilterName directory",_hltBTagJetFilterName.size(),0,_hltBTagJetFilterName.size());
    for (unsigned int ibin=0; ibin<_hltBTagJetFilterName.size(); ++ibin) {
        tlHist_hltBTagJetFilterName->GetXaxis()->SetBinLabel(ibin+1,_hltBTagJetFilterName[ibin].c_str());
        tlHist_hltBTagJetFilterName->Fill(ibin+0.5);  // just to avoid histogram being empty
    }
   
    //_hltL2JetFilterName = iConfig.getParameter<std::string>("hltL2JetFilterName");
    _hltL2JetFilterName = iConfig.getParameter< std::vector<std::string> > ("hltL2JetFilterName");
    tlHist_hltL2JetFilterName = fs->make<TH1F>("tlHist_hltL2JetFilterName","tlHist_hltL2JetFilterName directory",_hltL2JetFilterName.size(),0,_hltL2JetFilterName.size());
    for (unsigned int ibin=0; ibin<_hltL2JetFilterName.size(); ++ibin) {
        tlHist_hltL2JetFilterName->GetXaxis()->SetBinLabel(ibin+1,_hltL2JetFilterName[ibin].c_str());
        tlHist_hltL2JetFilterName->Fill(ibin+0.5);  // just to avoid histogram being empty
    }
    

    _hltMuonFilterName    = iConfig.getParameter<std::string>("hltMuonFilterName");
    
    // cuts --->
    _jetPtCut  = float(iConfig.getParameter<double>( "JetPtCut" ));
    _muonPtCut = float(iConfig.getParameter<double>( "MuonPtCut" )); 
    _minNumberOfJets = iConfig.getParameter<int>( "MinimumNumberOfJets" );
    
    _doMC = iConfig.getParameter<bool> ( "DoMonteCarlo" );
    
    _printOut = iConfig.getParameter<bool> ( "PrintOut" );
    
    _doNegativeBtags = iConfig.getParameter<bool> ( "DoNegativeBtags" );

    // initialize the generic trigger list
    _genericTriggerList = iConfig.getParameter< std::vector<std::string> > ( "triggerList" );
    
    
    // package the list of trigger names into a histogram, so that it is included in the ntuple as an object that can be cleanly merged
    gtlHist = fs->make<TH1F>("gtlHist","Generic trigger bit directory",_genericTriggerList.size(),0,_genericTriggerList.size());
    for (unsigned int ibin=0; ibin<_genericTriggerList.size(); ++ibin) {
        gtlHist->GetXaxis()->SetBinLabel(ibin+1,_genericTriggerList[ibin].c_str());
        gtlHist->Fill(ibin+0.5);  // just to avoid histogram being empty
    }
    
}


HbbMSSMAnalysis::~HbbMSSMAnalysis()
{
    // nothing to do... 
    
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
HbbMSSMAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    // increment event counter;
    if (_printOut) {
        std::cout << "we are in HbbMSSMAnalysis::analyze() " << std::endl;
    }
    _events++;
    
    // Accessing run, event and lumi-block
    run = iEvent.id().run();
    event = iEvent.id().event();
    lumi = iEvent.getLuminosityBlock().luminosityBlock();
    
    // Accessing trigger information ---> ALL TRIGGER PATHS
    
    edm::Handle<pat::TriggerPathCollection> triggerPathCollection;
    iEvent.getByLabel(_triggerPathSrc, triggerPathCollection);
    int ntrigs = int(triggerPathCollection->size());
    
    // clear precale array
    for (int trgSlot=0; trgSlot<25; ++trgSlot) {
        trigPrescale[trgSlot] = -1;
    }
    
    // loop over generic trigger list and set trigger accept word
    trgAccept = 0;
    for (std::vector<std::string>::iterator gtlIt = _genericTriggerList.begin();
         gtlIt != _genericTriggerList.end(); ++gtlIt)
    {
        
        // find appropriate trigger in list of all paths
        int trgSlot = -1;
        for (int it=0;it<ntrigs;++it) {
            pat::TriggerPath triggerPath = triggerPathCollection->at(it);
            bool accept = triggerPath.wasAccept();
            TString trigName(triggerPath.name());
            int iprescale = int(triggerPath.prescale());	
            std::cout << it << "   " << trigName << "  accepted : " << accept << "    prescale = " << iprescale << std::endl;
            if (trigName.Contains(*gtlIt)) {
                 	 std::cout << std::setw(3) << it << " "
                 	   << std::setw(30) << *gtlIt << "  " << std::setw(30) << trigName << "  accepted : " << accept << "    prescale = " << iprescale << std::endl;
                trgSlot = gtlIt - _genericTriggerList.begin();
                if (accept) {
                    if (trgSlot>=0 && trgSlot < 32) {
                        trgAccept = trgAccept | (1<<trgSlot);
                    }
                    else {
                        std::cout << "Bad trigger slot " << trgSlot << std::endl;
                    }
                }
                if (trgSlot<25) {
                    trigPrescale[trgSlot] = iprescale;
                }
            }
        }
    }
    
    // count triggers
    for (int it=0;it<ntrigs;++it)
    {
        pat::TriggerPath triggerPath = triggerPathCollection->at(it);
        bool accept = triggerPath.wasAccept();
        TString trigName(triggerPath.name());
        if (accept) {
            if (triggerCounter.find(triggerPath.name()) == triggerCounter.end()) {
                triggerCounter[triggerPath.name()] = 1;
            }
            else {
                ++triggerCounter[triggerPath.name()];
            }
        }
    }
    
    // look for trigger summary information
    
    // get trigger event
    edm::Handle< trigger::TriggerEvent > handleTriggerEvent_;
    
    if ( ! iEvent.getByLabel("hltTriggerSummaryAOD", handleTriggerEvent_ ) ) {
      edm::LogError( "errorTriggerEventValid" ) << "trigger::TriggerEvent product with InputTag " << "hltTriggerSummaryAOD::HLT" << " not in event";
      return;
    }
    iEvent.getByLabel( _triggerEventTag, handleTriggerEvent_ ); 
    
    //std::vector<reco::Particle>  HLTBtagMatched;
    std::vector<std::vector<reco::Particle> >  HLTBtagMatched;
    for(unsigned int ihltBTagJetFilterName = 0;
        ihltBTagJetFilterName < _hltBTagJetFilterName.size();
        ihltBTagJetFilterName++) {
      HLTBtagMatched.push_back(std::vector<reco::Particle> () ); //initiate with empty vector
    }
   //  if(HLTBtagMatched.size() != _hltBTagJetFilterName.size())
//       {
//         edm::LogError( "errorsizemismatch" ) << " " << HLTBtagMatched.size() << " "<< _hltBTagJetFilterName.size();
//         return;
//       }
//     for(unsigned int ihltBTagJetFilterName = 0;
//         ihltBTagJetFilterName < _hltBTagJetFilterName.size();
//         ihltBTagJetFilterName++) {
//       if(HLTBtagMatched.at(ihltBTagJetFilterName).size() != 0)
//         {
//           edm::LogError( "errorsizemismatch2" ) << " " << HLTBtagMatched.at(ihltBTagJetFilterName).size() << " "<< _hltBTagJetFilterName.size(); 
//           return; 
//         }
//     }

    //    std::vector<reco::Particle>  L2Matched;
std::vector<std::vector<reco::Particle> > L2Matched;
     for(unsigned int ihltL2JetFilterName = 0;
        ihltL2JetFilterName < _hltL2JetFilterName.size();
        ihltL2JetFilterName++) {
      L2Matched.push_back(std::vector<reco::Particle> () ); //initiate with empty vector
    }
   //  if(L2Matched.size() != _hltL2JetFilterName.size())
//       {
//         edm::LogError( "errorsizemismatch" ) << " " << L2Matched.size() << " "<< _hltL2JetFilterName.size();
//         return;
//       }
//     for(unsigned int ihltL2JetFilterName = 0;
//         ihltL2JetFilterName < _hltL2JetFilterName.size();
//         ihltL2JetFilterName++) {
//       if(L2Matched.at(ihltL2JetFilterName).size() != 0)
//         {
//           edm::LogError( "errorsizemismatch2" ) << " " << L2Matched.at(ihltL2JetFilterName).size() << " "<< _hltL2JetFilterName.size(); 
//           return; 
//         }
//     }



    std::vector<reco::Particle>  HltMuonMatched;
    
    const trigger::TriggerObjectCollection & toc(handleTriggerEvent_->getObjects());
    for ( unsigned int ia = 0; ia < handleTriggerEvent_->sizeFilters(); ++ ia)
    {
        std::string fullname = handleTriggerEvent_->filterTag(ia).encode();
        
        // filter statistics
        if (triggerObjectCounter.find(fullname) == triggerObjectCounter.end()) {
            triggerObjectCounter[fullname] = 1;
        }
        else {
            ++triggerObjectCounter[fullname];
        }
        
        // strip the process name
        std::string name;
        size_t p = fullname.find_first_of(':');
        if ( p != std::string::npos) {
            name = fullname.substr(0, p);
        }
        else {
            name = fullname;
        }
        
        for(unsigned int ihltBTagJetFilterName = 0;
            ihltBTagJetFilterName < _hltBTagJetFilterName.size();
            ihltBTagJetFilterName++) {
          if ( name == _hltBTagJetFilterName.at(ihltBTagJetFilterName) ) {
            int nkeys = -1;
            if ( &toc !=0 ) {
              const trigger::Keys & k = handleTriggerEvent_->filterKeys(ia);
              nkeys = k.size();
              for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) { 
                HLTBtagMatched.at(ihltBTagJetFilterName).push_back(toc[*ki].particle());
              }
            } else {
              std::cout << "Problem with toc" << std::endl;
            }
          }
        }

        for(unsigned int ihltL2JetFilterName = 0;
            ihltL2JetFilterName < _hltL2JetFilterName.size();
            ihltL2JetFilterName++) {
          if ( name ==  _hltL2JetFilterName.at(ihltL2JetFilterName)) {
            int nkeys = -1;
            if ( &toc !=0 ) {
              const trigger::Keys & k = handleTriggerEvent_->filterKeys(ia);
              nkeys = k.size();
              for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) { 
                L2Matched.at(ihltL2JetFilterName).push_back(toc[*ki].particle());
              }
            } else {
              std::cout << "Problem with toc" << std::endl;
            }
          }
        }
        if ( name == _hltMuonFilterName ) {
          int nkeys = -1;
          if ( &toc !=0 ) {
             const trigger::Keys & k = handleTriggerEvent_->filterKeys(ia);
             nkeys = k.size();
             for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) { 
               HltMuonMatched.push_back(toc[*ki].particle());
             }
          } else {
             std::cout << "Problem with toc" << std::endl;
          }
        }
    }	
    // 
    //   // jet - L1 Jet
    
    // L1 Jets
    l1extra::L1JetParticleCollection l1AllJets;
    for ( size_t j = 0 ; j < _l1ExtraJets.size(); ++j )  // loop over all input l1extrajets collections
    {
        edm::Handle<l1extra::L1JetParticleCollection> l1ExtraJetsHandler;
        iEvent.getByLabel(_l1ExtraJets.at(j), l1ExtraJetsHandler);
        const l1extra::L1JetParticleCollection & l1Jets = *(l1ExtraJetsHandler.product());
        for ( size_t i = 0 ; i < l1Jets.size() ; ++i ) l1AllJets.push_back(l1Jets[i]);
    }
    
    NumericSafeGreaterByPt<l1extra::L1JetParticle> compL1Jets;
    std::sort (l1AllJets.begin (), l1AllJets.end (), compL1Jets);
    
    l1NumberOfJets = 0;
    for ( size_t j = 0 ; j < l1AllJets.size(); ++j )
    {
        if ( l1NumberOfJets < 100 )
        {
            l1PtJet[l1NumberOfJets]     = l1AllJets[j].pt();
            l1EtaJet[l1NumberOfJets]    = l1AllJets[j].eta();
            l1PhiJet[l1NumberOfJets]    = l1AllJets[j].phi();
            l1PxJet[l1NumberOfJets]     = l1AllJets[j].px();
            l1PyJet[l1NumberOfJets]     = l1AllJets[j].py();
            l1PzJet[l1NumberOfJets]     = l1AllJets[j].pz();
            l1EnergyJet[l1NumberOfJets] = l1AllJets[j].energy();
            l1NumberOfJets++;
        }
    }
    
    // L2 Jets
    NumericSafeGreaterByPt<reco::Particle> compL2Jets;
    //std::sort (L2Matched.begin (), L2Matched.end (), compL2Jets);
    l2NumberOfJets = 0;
    if(L2Matched.size() > 0) {
      std::sort (L2Matched.front().begin (), L2Matched.front().end (), compL2Jets);
      for ( size_t j = 0 ; j < L2Matched.front().size(); ++j )
        {
          if ( l2NumberOfJets < 100 )
            {
              l2PtJet[l2NumberOfJets]     = L2Matched.front().at(j).pt();
              l2EtaJet[l2NumberOfJets]    = L2Matched.front().at(j).eta();
              l2PhiJet[l2NumberOfJets]    = L2Matched.front().at(j).phi();
              l2PxJet[l2NumberOfJets]     = L2Matched.front().at(j).px();
              l2PyJet[l2NumberOfJets]     = L2Matched.front().at(j).py();
              l2PzJet[l2NumberOfJets]     = L2Matched.front().at(j).pz();
              l2EnergyJet[l2NumberOfJets] = L2Matched.front().at(j).energy();
              l2NumberOfJets++;
            }
        }
    }
    else {
      edm::LogWarning("emptyL2Matched") << "L2Matched is empty";
    }
    // HLT BTag Jets
    NumericSafeGreaterByPt<reco::Particle> compL3BJets;
    //use only the first row
    //FIXME: will .front() throw an exception if for some reason HLTBtagMatched is empty? -> check for the size
    l3NumberOfBJets = 0;
    if(HLTBtagMatched.size() > 0) {
      std::sort (HLTBtagMatched.front().begin(), HLTBtagMatched.front().end(), compL3BJets);
      for ( size_t j = 0 ; j < HLTBtagMatched.front().size(); ++j )
        {
          if ( l3NumberOfBJets < 100 )
            {
              l3PtBJet[l3NumberOfBJets]     = HLTBtagMatched.front().at(j).pt();
              l3EtaBJet[l3NumberOfBJets]    = HLTBtagMatched.front().at(j).eta();
              l3PhiBJet[l3NumberOfBJets]    = HLTBtagMatched.front().at(j).phi();
              l3PxBJet[l3NumberOfBJets]     = HLTBtagMatched.front().at(j).px();
              l3PyBJet[l3NumberOfBJets]     = HLTBtagMatched.front().at(j).py();
              l3PzBJet[l3NumberOfBJets]     = HLTBtagMatched.front().at(j).pz();
              l3EnergyBJet[l3NumberOfBJets] = HLTBtagMatched.front().at(j).energy();
              l3NumberOfBJets++;
            }
        }
    }
    else {
      edm::LogWarning("emptyHLTBtagMatched") << "HLTBtagMatched is empty";
    }
    
    // Accessing generator level information
    if (_doMC)
    {
        
        if (firstEvent) {
            try {
                
                // accessing information on x-section, MC filter efficiency ======>
                
                edm::Handle<GenRunInfoProduct>  gi;
                iEvent.getRun().getByType(gi);
                
                _internalXsec = float(gi->internalXSec().value());
                _internalXsecE = float(gi->internalXSec().error());
                
                _externalXsecLO = float(gi->externalXSecLO().value());
                _externalXsecLOE = float(gi->externalXSecLO().error());
                
                _externalXsecNLO = float(gi->externalXSecNLO().value());
                _externalXsecNLOE = float(gi->externalXSecNLO().error());
                
                _filterEfficiency = float(gi->filterEfficiency());
                
                
                std::cout << "XSec              : " << _internalXsec << "+/-" << _internalXsecE << " pb" << std::endl;
                std::cout << "XSec LO           : " << _externalXsecLO << "+/-" << _externalXsecLOE << " pb" << std::endl;
                std::cout << "XSec NLO          : " << _externalXsecNLO << "+/-" << _externalXsecNLOE << " pb" << std::endl;
                std::cout << "Filter efficiency : " << _filterEfficiency << std::endl;
                std::cout << std::endl;
                
                treeGenInfo->Fill();
                
                firstEvent = false;
                
            }
            catch(...){
                std::cout << "Information on cross sections is not available in Event record..." << std::endl;
            }
        }
        
    }
    
    nPU = 0;
    nPUInTime = 0;
    if (_doMC)
    {
        // access pileup information from simulation
        // this will not work generally with releases before 42X
        edm::InputTag PileupSrc_("addPileupInfo");
        Handle<std::vector< PileupSummaryInfo > >  PupInfo;
        iEvent.getByLabel(PileupSrc_, PupInfo);
        
        std::vector<PileupSummaryInfo>::const_iterator PVI;
        for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
            int BX = PVI->getBunchCrossing();
            if ( _printOut ) {
                std::cout << " Pileup Information: bunchXing, nvtx: " << BX << " " << PVI->getPU_NumInteractions() << std::endl;
            }
            nPU += PVI->getPU_NumInteractions();
            if(BX == 0) { 
                nPUInTime = PVI->getPU_NumInteractions();
            }
        }
    }
    
    // accessing jet information;
    
    numberOfJets = 0;
    Handle<pat::JetCollection> jets;
    
    // vector of jet constituents;
    std::vector<std::vector<reco::PFCandidatePtr> > jetsConstituents;
    jetsConstituents.clear();
    
    edm::ParameterSet paramsLoose;
    edm::ParameterSet paramsMedium;
    edm::ParameterSet paramsTight;
  
    paramsLoose.addParameter(std::string("version"),std::string("FIRSTDATA"));
    paramsLoose.addParameter(std::string("quality"),std::string("LOOSE"));
  
    paramsMedium.addParameter(std::string("version"),std::string("FIRSTDATA"));
    paramsMedium.addParameter(std::string("quality"),std::string("MEDIUM"));
  
    paramsTight.addParameter(std::string("version"),std::string("FIRSTDATA"));
    paramsTight.addParameter(std::string("quality"),std::string("TIGHT"));
  
    PFJetIDSelectionFunctor jetIDLooseFunctor( paramsLoose );
    PFJetIDSelectionFunctor jetIDMediumFunctor( paramsMedium );
    PFJetIDSelectionFunctor jetIDTightFunctor( paramsTight );
  
    pat::strbitset ret_loose = jetIDLooseFunctor.getBitTemplate();
    pat::strbitset ret_medium = jetIDMediumFunctor.getBitTemplate();
    pat::strbitset ret_tight = jetIDTightFunctor.getBitTemplate();

    try {
        iEvent.getByLabel(_jetCollection,jets);
        int numberOfAllJets = jets->size();
        for (int iJet=0; iJet<numberOfAllJets; ++iJet) {
            pat::JetRef jet(jets,iJet);
            float jetPt = jet->pt();
            if (jetPt>_jetPtCut&&numberOfJets<100) { // cut on jet pt and protection against jet array size overflow
                
                etaJet[numberOfJets] = jet->eta();
                phiJet[numberOfJets] = jet->phi();
                ptJet[numberOfJets]  = jet->pt();
                pxJet[numberOfJets]  = jet->px();
                pyJet[numberOfJets]  = jet->py();
                pzJet[numberOfJets]  = jet->pz();
                energyJet[numberOfJets] = jet->energy();
                
                // jetId
                ret_loose.set(false);
                jetIDLoose[numberOfJets] = jetIDLooseFunctor(*jet, ret_loose);
                ret_medium.set(false);
                jetIDMedium[numberOfJets] = jetIDMediumFunctor(*jet, ret_medium);
                ret_tight.set(false);
                jetIDTight[numberOfJets] = jetIDTightFunctor(*jet, ret_tight);

                // HLT btag matching
                isJetWithHltBtagBitPattern[numberOfJets] = 0;//reset bit pattern
                deltaRWithHltBtag[numberOfJets] = -10000.0;
                deltaPtWithHltBtag[numberOfJets] = -10000.0;
                //isJetWithHltBtag[numberOfJets] = 0; //reset 'old' variable
                float theDeltaRWithHltBtag = -1.;
                float theDeltaPtWithHltBtag = -1.;
                //loop over all trigger objects
                for(unsigned int ihltBTagJetFilterName = 0;
                    ihltBTagJetFilterName < _hltBTagJetFilterName.size();
                    ihltBTagJetFilterName++) {
                  float tmptheDeltaRWithHltBtag = -1.;
                  float tmptheDeltaPtWithHltBtag = -1.;
                  bool trgMatch = IsJetMatchedToHLTBtag( jet, HLTBtagMatched.at(ihltBTagJetFilterName), 0.5, tmptheDeltaRWithHltBtag, tmptheDeltaPtWithHltBtag);
                  
                  if (trgMatch) {
                    if ( _printOut ) {
                      std::cout << "Jet match successful for event " << event << " jet " << numberOfJets 
                                << " deltaR=" << tmptheDeltaRWithHltBtag << " deltaPt=" << tmptheDeltaPtWithHltBtag << std::endl;
                    }
                    //isJetWithHltBtag[numberOfJets]  = 1;
                    if(ihltBTagJetFilterName < 32) {
                      isJetWithHltBtagBitPattern[numberOfJets] = isJetWithHltBtagBitPattern[numberOfJets] | (1<<ihltBTagJetFilterName);
                    }
                    else {
                      std::cout << "Bad trigger object slot " << ihltBTagJetFilterName << std::endl;
                    }
                    if(theDeltaRWithHltBtag < 0.0 || tmptheDeltaRWithHltBtag < theDeltaRWithHltBtag) {
                      theDeltaRWithHltBtag = tmptheDeltaRWithHltBtag;
                      theDeltaPtWithHltBtag = tmptheDeltaPtWithHltBtag;
                    }
                  }
                }
                deltaRWithHltBtag[numberOfJets]  = theDeltaRWithHltBtag;
                deltaPtWithHltBtag[numberOfJets]  = theDeltaPtWithHltBtag;
                // L1 matching
                isJetWithL1Jet[numberOfJets] = 0;
                deltaRWithL1Jet[numberOfJets] = -10000.0;
                deltaPtWithL1Jet[numberOfJets] = -10000.0;
                
                float theDeltaRWithL1Jet = -1;
                float theDeltaPtWithL1Jet = -1;
                bool l1trgMatch = IsJetMatchedToL1Jet( jet, l1AllJets, 0.5, theDeltaRWithL1Jet, theDeltaPtWithL1Jet);
                
                if (l1trgMatch) {
                    if ( _printOut ) {
                        std::cout << "Jet match with L1 successful for event " << event << " jet " << numberOfJets 
                        << " deltaR=" << theDeltaRWithL1Jet << " deltaPt=" << theDeltaPtWithL1Jet << std::endl;
                    }
                    isJetWithL1Jet[numberOfJets]  = 1;
                    deltaRWithL1Jet[numberOfJets]  = theDeltaRWithL1Jet;
                    deltaPtWithL1Jet[numberOfJets]  = theDeltaPtWithL1Jet;
                }
                
                // L2 matching
                isJetWithL2JetBitPattern[numberOfJets] = 0;
                deltaRWithL2Jet[numberOfJets] = -10000.0;
                deltaPtWithL2Jet[numberOfJets] = -10000.0;
                float theDeltaRWithL2Jet = -1;
                float theDeltaPtWithL2Jet = -1;
                for(unsigned int ihltL2JetFilterName = 0;
                    ihltL2JetFilterName < _hltL2JetFilterName.size();
                    ihltL2JetFilterName++) {
                 
                  float tmptheDeltaRWithL2Jet = -1.;
                  float tmptheDeltaPtWithL2Jet = -1.;
                  bool l2trgMatch = IsJetMatchedToHLTBtag( jet, L2Matched.at(ihltL2JetFilterName), 0.5, tmptheDeltaRWithL2Jet, tmptheDeltaPtWithL2Jet);
                
                  if (l2trgMatch) {
                    if ( _printOut ) {
                      std::cout << "Jet match with L2 successful for event " << event << " jet " << numberOfJets 
                                << " deltaR=" << tmptheDeltaRWithL2Jet << " deltaPt=" << tmptheDeltaPtWithL2Jet << std::endl;
                    }
                    if(ihltL2JetFilterName < 32) {
                      //  isJetWithL2Jet[numberOfJets]  = 1;
                      isJetWithL2JetBitPattern[numberOfJets] = isJetWithL2JetBitPattern[numberOfJets] | (1<<ihltL2JetFilterName);
                    }
                    else  {
                      std::cout << "Bad trigger object slot " << ihltL2JetFilterName << std::endl;
                    }
                    if(theDeltaRWithL2Jet < 0.0 || tmptheDeltaRWithL2Jet < theDeltaRWithL2Jet) {
                      theDeltaRWithL2Jet = tmptheDeltaRWithL2Jet;
                      theDeltaPtWithL2Jet = tmptheDeltaPtWithL2Jet;
                    }
                  }
                }
                deltaRWithL2Jet[numberOfJets]  = theDeltaRWithL2Jet;
                deltaPtWithL2Jet[numberOfJets]  = theDeltaPtWithL2Jet;
                // accessing information on jet constituents
                std::vector<reco::PFCandidatePtr> jetConstituents = jet->getPFConstituents(); 
                
                numberOfConstituentsInJet[numberOfJets] = jetConstituents.size();
                neutralHadronEnergyFraction[numberOfJets]=jet->neutralHadronEnergyFraction();
                photonEnergyFraction[numberOfJets]=jet->photonEnergyFraction();
                numberOfChargedConstituentsInJet[numberOfJets]=int(jet->chargedMultiplicity());
                electronEnergyFraction[numberOfJets]=jet->chargedEmEnergyFraction();
                chargedHadronEnergyFraction[numberOfJets]=jet->chargedHadronEnergyFraction();
                
                jetsConstituents.push_back(jetConstituents);
                
                
                // accessing information on matched parton/quark
                if (_doMC) {

		    Handle<reco::GenParticleCollection> genParticles;
		    iEvent.getByLabel(_genParticleCollection,genParticles);

		    partonFlavorJet[numberOfJets] = jet->partonFlavour();
		    hflContentJet[numberOfJets] = jetHflContent( jet, genParticles );

                    const reco::GenParticle * matchedGenParton = jet->genParton();
                    
                    isJetMatchedPartonExist[numberOfJets] = false;
                    isJetMatchedPartonFromHiggsBosonDecay[numberOfJets] = false;
                    
                    if (matchedGenParton!=NULL) {
                        
                        isJetMatchedPartonExist[numberOfJets] = true;
                        flavorJetMatchedParton[numberOfJets]  = matchedGenParton->pdgId();
                        pxJetMatchedParton[numberOfJets]      = matchedGenParton->px();
                        pyJetMatchedParton[numberOfJets]      = matchedGenParton->py();
                        pzJetMatchedParton[numberOfJets]      = matchedGenParton->pz();	   
                        energyJetMatchedParton[numberOfJets]  = matchedGenParton->pz();	   
                        
                        isJetMatchedPartonFromHiggsBosonDecay[numberOfJets] = false;
                        
                        bool traceBackGenParticleRecord = true;
                        
                        const reco::Candidate * intermediateParton = NULL;
                        bool isFirstIteration = true;
                        
                        // this is trick to check if parton stems from Higgs decay
                        while (traceBackGenParticleRecord) {
                            
                            int numberOfMothers = 0;
                            if (isFirstIteration)
                                numberOfMothers = matchedGenParton->numberOfMothers();
                            else 
                                numberOfMothers = intermediateParton->numberOfMothers();
                            
                            if (numberOfMothers==0||numberOfMothers>1) 
                                traceBackGenParticleRecord = false;
                            else {
                                
                                const reco::Candidate * mother = NULL;
                                
                                if (isFirstIteration) 
                                    mother = matchedGenParton->mother(0);
                                else 
                                    mother = intermediateParton->mother(0);
                                
                                int pdgMother = mother->pdgId();
                                
                                if (pdgMother==25||pdgMother==35||pdgMother==36) { // mother is Higgs boson ! 
                                    isJetMatchedPartonFromHiggsBosonDecay[numberOfJets] = true;
                                    traceBackGenParticleRecord = false;
                                }
                                else if (pdgMother==flavorJetMatchedParton[numberOfJets]) { // seems like it is the same parton before gluon radiation 
                                    isFirstIteration = false;
                                    intermediateParton = mother;
                                }
                                else  // mother is NOT a Higgs boson : stop tracing back generator particle list!
                                    traceBackGenParticleRecord = false;
                                
                            }
                            
                        }
                        
                    }
                    
                }
                
                // BTag information ---->
                jetBProbBJetTag[numberOfJets] = -10000.0;
                jetProbBJetTag[numberOfJets] = -10000.0;
                tcHPBJetTag[numberOfJets] = -10000.0;
                tcHEBJetTag[numberOfJets] = -10000.0;
		ntcHPBJetTag[numberOfJets] = -10000.0;
		ntcHEBJetTag[numberOfJets] = -10000.0;
                svHEBJetTag[numberOfJets] = -10000.0;
                svHPBJetTag[numberOfJets] = -10000.0;
                combSVBJetTag[numberOfJets] = -10000.0;
                combSVMVABJetTag[numberOfJets] = -10000.0;
                tcMinus2ndBJetTag[numberOfJets] = +10000.0;
                tcMinus3rdBJetTag[numberOfJets] = +10000.0;
                tcNumberOfSelectedTracks[numberOfJets] = -1;
                
                const std::vector< std::pair< std::string, float > > pairDiscriVector = jet->getPairDiscri();
                int nDiscri = pairDiscriVector.size();
                if ( _printOut ) { // printOut of Jet info
                    
                    std::cout << " Jet " << numberOfJets << " : Et = " << jet->et()
                    << " ; Eta = " << jet->eta() 
                    << " ; Phi = " << jet->phi(); 
                    
                    if (_doMC)
                        std::cout << " ; is matched parton exist :  " << isJetMatchedPartonExist[numberOfJets]
                        << " ; is parton from Higgs : " << isJetMatchedPartonFromHiggsBosonDecay[numberOfJets] 
                        << " ; parton flavor = " << flavorJetMatchedParton[numberOfJets];
                    
                    std::cout << std::endl;
                }
                
                for (int iD=0;iD<nDiscri;++iD) {
                    std::pair<std::string, float> pairDiscri = pairDiscriVector[iD];
                    TString nameOfDiscri(pairDiscri.first);
                    float tag = pairDiscri.second;
                    
                    if (_printOut) 
                        std::cout << "     " << nameOfDiscri << " : " << tag << std::endl;    
                    
                    if (nameOfDiscri==TString("jetBProbabilityBJetTags")) {
                        jetBProbBJetTag[numberOfJets] = tag;
                        //	     std::cout << " jetBProbabilityBJetTags OK " << std::endl;
                    }
                    if (nameOfDiscri==TString("jetProbabilityBJetTags")) {
                        jetProbBJetTag[numberOfJets] = tag;
                        //	     std::cout << " jetProbabilityBJetTags OK " << std::endl;
                    }
                    if (nameOfDiscri==TString("trackCountingHighPurBJetTags")) {
                        tcHPBJetTag[numberOfJets] = tag;
                        //	     std::cout << " trackCountingHighPurBJetTags OK" << std::endl;
                    }
                    if (nameOfDiscri==TString("trackCountingHighEffBJetTags")) {
                        tcHEBJetTag[numberOfJets] = tag;
                        //	     std::cout << " trackCountingHighEffBJetTags OK " << std::endl;
                    }
		    if (_doNegativeBtags) {
		      if (nameOfDiscri==TString("negativeTrackCountingHighPurAODPFlow")) {
			ntcHPBJetTag[numberOfJets] = tag;
			//	     std::cout << " negativeTrackCountingHighPurBJetTags OK" << std::endl;
		      }
		      if (nameOfDiscri==TString("negativeTrackCountingHighEffJetTags")) {
			ntcHEBJetTag[numberOfJets] = tag;
			//	     std::cout << " negativeTrackCountingHighEffBJetTags OK " << std::endl;
		      }
		    }
                    if (nameOfDiscri==TString("simpleSecondaryVertexHighEffBJetTags")) {
                        svHEBJetTag[numberOfJets] = tag;
                        //	     std::cout << " simpleSecondaryVertexHighEffBJetTags OK " << std::endl;
                    }
                    if (nameOfDiscri==TString("simpleSecondaryVertexHighPurBJetTags")) {
                        svHPBJetTag[numberOfJets] = tag;
                        //	     std::cout << " simpleSecondaryVertexHighPurBJetTags OK " << std::endl;
                    }
                    if (nameOfDiscri==TString("combinedSecondaryVertexBJetTags")) {
                        combSVBJetTag[numberOfJets] = tag;
                        //	     std::cout << " combinedSecondaryVertexBJetTags OK " << std::endl;
                    }
                    if (nameOfDiscri==TString("combinedSecondaryVertexMVABJetTags")) {
                        combSVMVABJetTag[numberOfJets] = tag;	  
                        //	     std::cout << " combinedSecondaryVertexBJetTags OK " << std::endl;
                    }
                }
                
                // try to access trackIPTagInfo (needed for "negative" IP tags)
                std::string s3("impactParameter");
                if (_printOut) {
                    std::cout << "Jet hasTagInfo for label " << s3 << " = " << jet->hasTagInfo(s3) << std::endl;
                }
                if (jet->hasTagInfo(s3)) {
                    const reco::TrackIPTagInfo* ipti = jet->tagInfoTrackIP(s3);
                    if (ipti != NULL) {
                        const std::vector<reco::TrackIPTagInfo::TrackIPData> & theIpData = ipti->impactParameterData(); 
                        if (_printOut) std::cout << "ipData size=" << theIpData.size() << std::endl;
                        
                        reco::TrackRefVector selTracks=ipti->selectedTracks();
                        for(unsigned int j=0;j<selTracks.size(); j++) {
                            double ipsig = ipti->impactParameterData()[j].ip3d.significance();
                            if (_printOut) std::cout << "  selTrack " << j << " ipsid3d=" << ipsig << std::endl;
                        }
                        
                        // test sortedIndices
                        if (_printOut) {
                            for (unsigned int k=0; k< (ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )).size(); ++k) {
                                std::cout << " sorted index " << ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )[k] 
                                << " ipsig3d=" << ipti->impactParameterData()[ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )[k]].ip3d.significance()
                                << std::endl;
                            }
                        }
                        reco::TrackRefVector sortedTracks  =  ipti -> sortedTracks( ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig ) );
                        numberOfTracksJet[numberOfJets] = 0;
//                        std::cout << "Jet index = " << numberOfJets << std::endl;
                        for ( unsigned int  k = 0 ; k < 3 ; ++k ) 
                        {
                              trackPxJet[numberOfJets][k] = -99999;
                              trackPyJet[numberOfJets][k] = -99999;
                              trackPzJet[numberOfJets][k] = -99999;
                        }
                        
                        for (unsigned int k=0; k < sortedTracks.size(); ++k)
                        {
                           reco::Track track( *(sortedTracks.at(k)) );
                           if ( k < 3 ) 
                           {
                              trackPxJet[numberOfJets][k] = track.px();
                              trackPyJet[numberOfJets][k] = track.py();
                              trackPzJet[numberOfJets][k] = track.pz();
                              ++numberOfTracksJet[numberOfJets];
                           }
                        }
                        
                        
                        // find positive and negative tags
                        float ipsig3d_2nd = -100;
                        float ipsig3d_3rd = -100;
                        
                        float ipsig3d_minus_2nd = +100;
                        float ipsig3d_minus_3rd = +100;
                        
                        int nSel = ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig ).size();
                        
                        int nth = 2;
                        if (nSel>=nth) {
                            ipsig3d_2nd = ipti->impactParameterData()[ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )[nth-1]].ip3d.significance();
                        }
                        nth = 3;
                        if (nSel>=nth) {
                            ipsig3d_3rd = ipti->impactParameterData()[ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )[nth-1]].ip3d.significance();
                        }
                        nth = nSel -1;
                        if (nth>0) {
                            ipsig3d_minus_2nd = ipti->impactParameterData()[ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )[nth-1]].ip3d.significance();
                        }
                        nth = nSel -2;
                        if (nth>0) {
                            ipsig3d_minus_3rd = ipti->impactParameterData()[ipti->sortedIndexes( reco::TrackIPTagInfo::IP3DSig )[nth-1]].ip3d.significance();
                        }
                        
                        if (_printOut) {
                            std::cout << " ipsig3d_2nd = " << ipsig3d_2nd << "   TCHE=" << tcHEBJetTag[numberOfJets] << std::endl;
                            std::cout << " ipsig3d_3rd = " << ipsig3d_3rd << "   TCHP=" << tcHPBJetTag[numberOfJets] << std::endl;
                            std::cout << " ipsig3d_minus_2nd = " << ipsig3d_minus_2nd  << std::endl;
                            std::cout << " ipsig3d_minus_3rd = " << ipsig3d_minus_3rd  << std::endl;
                        }
                        tcMinus2ndBJetTag[numberOfJets] = ipsig3d_minus_2nd;
                        tcMinus3rdBJetTag[numberOfJets] = ipsig3d_minus_3rd;
                        tcNumberOfSelectedTracks[numberOfJets] = nSel;
                        
                        // 	     Handle<reco::VertexCollection> pVtxs2;
                        // 	     iEvent.getByLabel(_primaryVtxCollection,pVtxs2);
                        // 	     reco::VertexCollection::const_iterator vtxItr = find(pVtxs2->begin(),pVtxs2->end(),ipti->primaryVertex() );
                        
                    }
                }
                
                // try to access trackIPTagInfo (needed for "negative" IP tags)
                nSvTracksJet[numberOfJets] = 0;
                svMassJet[numberOfJets] = -1;
                svFDsigJet[numberOfJets] = -10000;
                
                
                std::string sv3("secondaryVertex");
                if (_printOut) {
                    std::cout << "Jet hasTagInfo for label " << sv3 << " = " << jet->hasTagInfo(sv3) << std::endl;
                }
                if (jet->hasTagInfo(sv3)) {
                    const reco::SecondaryVertexTagInfo* svti = jet->tagInfoSecondaryVertex(sv3);
                    if (svti != NULL) {
                       // Vertices in a jet are not ordered in flight distance significance!!!
                       // Find the vertex with highest flight distance significance
                       double fdsig;
                       double highFDsig = -10000;
                       int highFDsigIndx = -1;
                       
                       for ( unsigned int j = 0; j < svti -> nVertices() ; ++j )
                       {
                          fdsig =  svti->reco::SecondaryVertexTagInfo::flightDistance(j).significance();
                          if ( fdsig >= highFDsig ) 
                          {
                             highFDsig = fdsig;
                             highFDsigIndx = j;
                          }
                       }
                       if ( highFDsigIndx >=0 ) // jet has at least one vertex
                       {
                          const reco::Vertex & secondaryVertex = svti -> secondaryVertex(highFDsigIndx);
                          nSvTracksJet[numberOfJets] = secondaryVertex.nTracks();
                          svMassJet[numberOfJets] = secondaryVertex.p4().mass();
                          svFDsigJet[numberOfJets] = svti->reco::SecondaryVertexTagInfo::flightDistance(highFDsigIndx).significance();
                       }
                    }
               }
                        
                
                
                // incrementing jet counter
                numberOfJets++;
                
            }
            
        }
        
    }
    catch(...) { 
        std::cout << "pat::JetCollection " << _jetCollection << " is not found in Event record..." << std::endl;
    }
    
    
    
    // cut on number of jets --->
    if (numberOfJets<_minNumberOfJets) {
        if (_printOut) {
            std::cout << std::endl;
            std::cout << "Number of jets passing pt cut is less than " << _jetPtCut << " < " << _minNumberOfJets << std::endl;
            std::cout << "Event is skipped..." << std::endl;    
            std::cout << std::endl;
            std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
            std::cout << std::endl;
        }
        return;
    }
    
    // all primary vertices ---->
    
    nPV = 0;
    
    Handle<reco::VertexCollection> pVtxs;
    
    try {
        
        iEvent.getByLabel(_primaryVtxCollection,pVtxs);
        int numberOfPV = int(pVtxs->size());
        
        if (_printOut) {
            std::cout << std::endl;
            std::cout << "Initial number of primary vertices = " << numberOfPV << std::endl;
            std::cout << std::endl;
        }
        for (int iPV=0;iPV<numberOfPV;++iPV) {
            if (nPV<100) {
                reco::Vertex vertex = pVtxs->at(iPV);
                chi2PV[nPV] = vertex.chi2();
                ndofPV[nPV] = vertex.ndof();
                if (_printOut) {
                    std::cout << "vertex " << iPV << " ndof() = " <<  vertex.ndof() << std::endl;
                }
                probPV[nPV] = float(TMath::Prob(double(chi2PV[nPV]),int(ndofPV[nPV])));
                nTrkPV[nPV] = vertex.tracksSize();
                float ptsum = 0.0;
                reco::Vertex::trackRef_iterator iTrack;
                for(iTrack  = vertex.tracks_begin();
                    iTrack != vertex.tracks_end();++iTrack){
                    ptsum += (*iTrack)->pt();
                }
                sumPtPV[nPV] = ptsum;
                xPV[nPV] = vertex.x();
                yPV[nPV] = vertex.y();
                zPV[nPV] = vertex.z();
                nPV++;
            }
        }
    }
    catch(...) {
        std::cout << std::endl;
        std::cout << "Vertex Collection " << _primaryVtxCollection << " is not found in Event record..." << std::endl;
        std::cout << std::endl;
    }
    
    
    if (nPV<1) {
        if (_printOut) {
            std::cout << std::endl;
            std::cout << "There are no good primary vertices in an Event..." << std::endl;
            std::cout << "Event is skipped..." << std::endl;
            std::cout << std::endl;
            std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
            std::cout << std::endl;
        }
        return; 
    }
    
    // acessing MET collection
    Handle<pat::METCollection> patMet;
    // initializing MET --->
    met = 0;
    missingPx = 0;
    missingPy = 0;
    try {
        iEvent.getByLabel(_metCollection,patMet);
        int nPatMET = int(patMet->size());
        met = 0;
        missingPx = 0;
        missingPy = 0;
        
        if (nPatMET>0) {
            pat::METRef PATMet(patMet,0);
            met = PATMet->et();
            missingPx = PATMet->px();
            missingPy = PATMet->py();
        }
        
        if (_printOut) {
            std::cout << std::endl;
            std::cout << "PAT MET  :   (Px,Py)=(" << missingPx << "," << missingPy << ") ;  MET = " << met << std::endl;
            std::cout << std::endl;
        }
        
    }
    catch(...) {
        std::cout << "MET Collection " << _metCollection << " is not found in Event record..." << std::endl;
    }
    
    // accessing muons 
    numberOfMuons = 0;
    Handle<pat::MuonCollection> patMuons;
    
    try {
        iEvent.getByLabel(_muonCollection,patMuons);
        int nAllMuons = int(patMuons->size());
        
        for (int iMuon=0; iMuon<nAllMuons; ++iMuon) {
            pat::MuonRef patMuon(patMuons,iMuon);
            float muonPt = patMuon->pt();
            
            if (muonPt>_muonPtCut&&numberOfMuons<20) {
                float charge = patMuon->charge();
                if (charge>0.1) 
                    chargeMuon[numberOfMuons] = 1;
                else 
                    chargeMuon[numberOfMuons] = -1;
                
                isTrackerMuon[numberOfMuons] = patMuon->isTrackerMuon();
                isGlobalMuon[numberOfMuons] = patMuon->isGlobalMuon();
                
                etaMuon[numberOfMuons] = patMuon->eta();
                phiMuon[numberOfMuons] = patMuon->phi();
                ptMuon[numberOfMuons]  = patMuon->pt();
                
                pxMuon[numberOfMuons] = patMuon->px();
                pyMuon[numberOfMuons] = patMuon->py();
                pzMuon[numberOfMuons] = patMuon->pz();
                
                TVector3 muon3P(pxMuon[numberOfMuons],pyMuon[numberOfMuons],pzMuon[numberOfMuons]);
                
                isMuonInJet[numberOfMuons] = false;
                associatedJetIndex[numberOfMuons] = -1;
                
                bool jetIsFound = false;
                
                for (int iSelectedJet=0; iSelectedJet<numberOfJets; ++iSelectedJet) {
                    std::vector<reco::PFCandidatePtr> jetConstituents = jetsConstituents.at(iSelectedJet);
                    int numberOfConstituents = int(jetConstituents.size());
                    for (int iConstituent=0; iConstituent<numberOfConstituents;++iConstituent) {
                        float constituentPx = jetConstituents.at(iConstituent)->p4().Px();
                        float constituentPy = jetConstituents.at(iConstituent)->p4().Py();
                        float constituentPz = jetConstituents.at(iConstituent)->p4().Pz();	     
                        TVector3 constituent3P(constituentPx,constituentPy,constituentPz);
                        
                        TVector3 difference3P = muon3P - constituent3P;
                        float deltaP = float(difference3P.Mag());
                        
                        if (deltaP<1e-5) {
                            isMuonInJet[numberOfMuons] = true;
                            associatedJetIndex[numberOfMuons] = iSelectedJet;
                            jetIsFound = true;
                            break;
                        }
                        
                        
                    }
                    
                    if (jetIsFound)
                        break; 
                    
                }
                
                // PFlow based absolute isolation variables ---->
                pfChargedHadronIsoMuon[numberOfMuons] = patMuon->chargedHadronIso();
                pfNeutralHadronIsoMuon[numberOfMuons] = patMuon->neutralHadronIso();
                pfGammaIsoMuon[numberOfMuons] = patMuon->photonIso();
                
                // accessing matching generator muon
                if (_doMC) {
                    const reco::GenParticle * genMuon = patMuon->genLepton();
                    
                    isMatchedGeneratorMuonExist[numberOfMuons] = false;
                    if (genMuon!=NULL) {
                        isMatchedGeneratorMuonExist[numberOfMuons] = true;
                        float qGenMuon = genMuon->charge();
                        if (qGenMuon>0.1)
                            chargeGenMuon[numberOfMuons] = 1;
                        else 
                            chargeGenMuon[numberOfMuons] = -1;
                        
                        pxGenMuon[numberOfMuons] = genMuon->px();
                        pyGenMuon[numberOfMuons] = genMuon->py();
                        pzGenMuon[numberOfMuons] = genMuon->pz();
                        
                        
                        // accessing PDG code of muon's mother
                        
                        pdgMotherGenMuon[numberOfMuons] = 0;
                        
                        int pdgGeneratedMuon = genMuon->pdgId();
                        
                        bool traceBackGenParticleRecord = true;
                        
                        const reco::Candidate * intermediateMuon = NULL;
                        bool isFirstIteration = true;
                        
                        // this is a trick to get muons' mother PDG code 
                        while (traceBackGenParticleRecord) {
                            
                            int numberOfMothers = 0;
                            if (isFirstIteration)
                                numberOfMothers = genMuon->numberOfMothers();
                            else 
                                numberOfMothers = intermediateMuon->numberOfMothers();
                            
                            if (numberOfMothers==0||numberOfMothers>1) 
                                traceBackGenParticleRecord = false;
                            else {
                                
                                const reco::Candidate * mother = NULL;
                                
                                if (isFirstIteration) 
                                    mother = genMuon->mother(0);
                                else 
                                    mother = intermediateMuon->mother(0);
                                
                                int pdgMother = mother->pdgId();
                                
                                if (pdgMother==pdgGeneratedMuon) { // seems like it is the same muon before final state radiation 
                                    isFirstIteration = false;
                                    intermediateMuon = mother;
                                }
                                else  { // this is genuine muon's mother  
                                    traceBackGenParticleRecord = false;
                                    pdgMotherGenMuon[numberOfMuons] = pdgMother;
                                }
                            }
                        }
                    }
                }
                
                if (_printOut) {
                    std::cout << "Muon " << numberOfMuons 
                    << " ; (Px,Py,Pz)=(" << pxMuon[numberOfMuons] << "," << pyMuon[numberOfMuons] << "," << pzMuon[numberOfMuons] << ")";
                    if (_doMC)
                        std::cout << " ; GenMuon (Px,Py,Pz)=(" << pxMuon[numberOfMuons] << "," << pyMuon[numberOfMuons] << "," << pzMuon[numberOfMuons] << ")"
                        << "  PDG(Mother)=" <<  pdgMotherGenMuon[numberOfMuons];
                    std::cout << std::endl;
                }
                
                // incrementing muon counter
                numberOfMuons++;
                
            }
            
        }
        
    }
    catch(...) {
        std::cout << "Muon Collection " << _muonCollection << " is not found in Event record..." << std::endl;
    }
    
   // L1 muons
   edm::Handle<l1extra::L1MuonParticleCollection> l1MuonsHandler;
   iEvent.getByLabel(_l1Muons, l1MuonsHandler);
   const l1extra::L1MuonParticleCollection & l1Muons = *(l1MuonsHandler.product());
   
   l1NumberOfMuons = 0;
   for ( size_t j = 0 ; j < l1Muons.size(); ++j )
   {
      if ( l1NumberOfMuons < 100 )
      {
         l1PtMuon[l1NumberOfMuons]     = l1Muons[j].pt();
         l1EtaMuon[l1NumberOfMuons]    = l1Muons[j].eta();
         l1PhiMuon[l1NumberOfMuons]    = l1Muons[j].phi();
         l1PxMuon[l1NumberOfMuons]     = l1Muons[j].px();
         l1PyMuon[l1NumberOfMuons]     = l1Muons[j].py();
         l1PzMuon[l1NumberOfMuons]     = l1Muons[j].pz();
         l1EnergyMuon[l1NumberOfMuons] = l1Muons[j].energy();
         l1ChargeMuon[l1NumberOfMuons] = l1Muons[j].charge();
         l1NumberOfMuons++;
      }
   }

   // HLT Muons
   NumericSafeGreaterByPt<reco::Particle> compHltMuons;
   std::sort (HltMuonMatched.begin (), HltMuonMatched.end (), compHltMuons);
   
   hltNumberOfMuons = 0;
   for ( size_t j = 0 ; j < HltMuonMatched.size(); ++j )
   {
      if ( hltNumberOfMuons < 100 )
      {
         hltPtMuon[hltNumberOfMuons]     = HltMuonMatched[j].pt();
         hltEtaMuon[hltNumberOfMuons]    = HltMuonMatched[j].eta();
         hltPhiMuon[hltNumberOfMuons]    = HltMuonMatched[j].phi();
         hltPxMuon[hltNumberOfMuons]     = HltMuonMatched[j].px();
         hltPyMuon[hltNumberOfMuons]     = HltMuonMatched[j].py();
         hltPzMuon[hltNumberOfMuons]     = HltMuonMatched[j].pz();
         hltEnergyMuon[hltNumberOfMuons] = HltMuonMatched[j].energy();
         hltChargeMuon[hltNumberOfMuons] = HltMuonMatched[j].charge();
         hltNumberOfMuons++;
      }
   }
   
    if (_doMC)
    {
        // accessing generator level information (Higgs Boson and its decay products)
        
        // initialization of tree variables
        pdgHiggsBoson = 0;
        etaHiggsBoson = 0;
        phiHiggsBoson = 0;
        ptHiggsBoson = 0;
        massHiggsBoson = 0;
        pxHiggsBoson = 0;
        pyHiggsBoson = 0;
        pzHiggsBoson = 0;
        energyHiggsBoson = 0;
        numberOfPartons = 0;   // added to avoid non-initialization  RM
        
        Handle<reco::GenParticleCollection> genParticles;
        
        bool bosonFound = false; 
        try {
            
            iEvent.getByLabel(_genParticleCollection,genParticles);
            
            int numberOfGenParticles = genParticles->size();
            
            for (int iP=0; iP<numberOfGenParticles; iP++) {
                reco::GenParticleRef particle(genParticles,iP);
                int partPdg = particle->pdgId();
                if ( partPdg==25 || partPdg==35 || partPdg==36 ) { // one of the Higgs bosons 
                    // h (pdg=25) , H (pdg=35) or A (pdg=36)
                    pdgHiggsBoson = partPdg;
                    etaHiggsBoson = particle->eta();
                    phiHiggsBoson = particle->phi();
                    ptHiggsBoson  = particle->pt();
                    massHiggsBoson = particle->mass();
                    pxHiggsBoson = particle->px();
                    pyHiggsBoson = particle->py();
                    pzHiggsBoson = particle->pz();
                    energyHiggsBoson = particle->energy();
                    numberOfPartons = 0;
                    
                    // accessing daughters --->
                    // this is a trick to get partons/quarks from Higgs decay
                    bool partonsNotYetFound = true;
                    bool isFirstIteration = true;
                    const reco::Candidate * intermediateDaugher = NULL;
                    numberOfPartons = 0; // number of partons is initially set to 0
                    while (partonsNotYetFound) {
                        
                        int numberOfDaughters = 0;
                        
                        if (isFirstIteration) 
                            numberOfDaughters = particle->numberOfDaughters();
                        else 
                            numberOfDaughters = intermediateDaugher->numberOfDaughters();
                        
                        if (numberOfDaughters==0)	       
                            partonsNotYetFound = false;
                        else if (numberOfDaughters==1) { 
                            const reco::Candidate * daughter = NULL;
                            if (isFirstIteration) 
                                daughter = particle->daughter(0);
                            else 
                                daughter = intermediateDaugher->daughter(0);
                            intermediateDaugher = daughter;
                            isFirstIteration = false;
                        }
                        else {
                            int nOfPartons = numberOfDaughters ; 
                            for (int iDaughter=0; iDaughter<nOfPartons; ++iDaughter) {
                                const reco::Candidate * daughter = NULL;
                                if (isFirstIteration)
                                    daughter = particle->daughter(iDaughter);
                                else 
                                    daughter = intermediateDaugher->daughter(iDaughter);
                                int pdgDaughter = daughter->pdgId();
                                if ( (pdgDaughter!=pdgHiggsBoson) && (numberOfPartons<10) ) {
                                    pdgParton[numberOfPartons] = daughter->pdgId();
                                    pxParton[numberOfPartons]  = daughter->px();
                                    pyParton[numberOfPartons]  = daughter->py();
                                    pzParton[numberOfPartons]  = daughter->pz();
                                    etaParton[numberOfPartons]  = daughter->eta();
                                    phiParton[numberOfPartons]  = daughter->phi();
                                    ptParton[numberOfPartons]  = daughter->pt();		   
                                    energyParton[numberOfPartons]  = daughter->energy();
                                    numberOfPartons++;
                                }
                            }
                            partonsNotYetFound = false;
                        }
                        
                    }
                    
                    
                    bosonFound = true;
                    break;
                }
                
            }
            if ( _printOut ) {
                std::cout << std::endl;
                if (bosonFound) {
                    std::cout << "Higgs Boson PDG Code = " << pdgHiggsBoson
                    << " ; Px = " << pxHiggsBoson 
                    << " ; Py = " << pyHiggsBoson 
                    << " ; Pz = " << pzHiggsBoson 
                    << " ; Energy = " << energyHiggsBoson
                    << " ; Mass = " << massHiggsBoson << std::endl;
                    for (int iParton=0; iParton<numberOfPartons; ++iParton) {
                        std::cout << "Parton  " << iParton << " : PDG Code = " << pdgParton[iParton]
                        << " ; Px = " << pxParton[iParton] 
                        << " ; Py = " << pyParton[iParton] 
                        << " ; Pz = " << pzParton[iParton]
                        << " ; Energy = " << energyParton[iParton]
                        << " ; Pt = " << ptParton[iParton] 
                        << " ; Eta = " << etaParton[iParton] 
                        << " ; Phi = " << phiParton[iParton] << std::endl;
                    }
                }
                else {
                    std::cout << "Higgs boson is not found in the generated particle list..." << std::endl;
                }
                std::cout << std::endl;
            }
            
        }
        catch(...){
            std::cout << "reco::GenParticle collection " << _genParticleCollection << " is not found in Event record..." << std::endl;
        }
        
        genNumberOfJets = 0;
        Handle<reco::GenJetCollection> genJets;
        iEvent.getByLabel(_genJetCollection,genJets);
        
        for ( unsigned int j = 0 ; j < genJets->size(); ++j )
        {
            const reco::GenJet genJet = (*genJets)[j];
            if ( genNumberOfJets < 100 )
            {
                genPtJet[genNumberOfJets]     = genJet.pt();
                genEtaJet[genNumberOfJets]    = genJet.eta();
                genPhiJet[genNumberOfJets]    = genJet.phi();
                genPxJet[genNumberOfJets]     = genJet.px();
                genPyJet[genNumberOfJets]     = genJet.py();
                genPzJet[genNumberOfJets]     = genJet.pz();
                genEnergyJet[genNumberOfJets] = genJet.energy();
                genNumberOfJets++;
            }
            
        }
        
        
    }
    
    // filling tree --->
    tree->Fill();
    
    // event is selected --->
    _selectedEvents++;
    
}


// ------------ method called once each job just before starting event loop  ------------
void 
HbbMSSMAnalysis::beginJob()
{
    
    std::cout << std::endl;
    std::cout << " HbbMSSMAnalysis::beginJob() --> " << std::endl;
    std::cout << " Creating Trees..." << std::endl;
    
    // Creating dynamical TTree object with TFileService
    tree = fs->make<TTree>("HBBTo4B","HBBTo4B");
    
    tree->Branch("Run",&run,"Run/I"); // Run number
    tree->Branch("Event",&event,"Event/I"); // Event number
    tree->Branch("Lumi",&lumi,"Lumi/I"); // Lumi block number
    
    // pileup information
    tree->Branch("NumberOfPU",&nPU,"NumberOfPU/I"); // number of PU events in simulation
    tree->Branch("NumberOfPUInTime",&nPUInTime,"NumberOfPUInTime/I"); // number of PU events in simulation
    
    // bit pattern for triggers defining PD
    tree->Branch("trgAccept",&trgAccept,"trgAccept/I");  // trigger bit pattern
    
    // primary vertices
    tree->Branch("NumberOfPV",&nPV,"NumberOfPV/I"); // number of PVs
    tree->Branch("ProbPV",probPV,"ProbPV[NumberOfPV]/F"); // probability of PV fit
    tree->Branch("NdofPV",ndofPV,"ProbPV[NumberOfPV]/F"); // ndof of PV fit (float not integer!)
    tree->Branch("NTrkPV",nTrkPV,"NTrkPV[NumberOfPV]/I"); // number of tracks in PV
    tree->Branch("Chi2PV",chi2PV,"Chi2PV[NumberOfPV]/F"); // chi2 of PV fit
    tree->Branch("XPV",xPV,"XPV[NumberOfPV]/F"); // x coordinate of PV
    tree->Branch("YPV",yPV,"YPV[NumberOfPV]/F"); // y coordinate of PV
    tree->Branch("ZPV",zPV,"ZPV[NumberOfPV]/F"); // z coordinate of PV
    tree->Branch("SumTrkPtPV",sumPtPV,"SumTrkPtPV[NumberOfPV]/F"); // scalar sum of tracks' pt 
    
    // l1 jets
    tree->Branch("L1NumberOfJets",&l1NumberOfJets,"L1NumberOfJets/I");
    tree->Branch("L1JetEta",l1EtaJet,"L1JetEta[L1NumberOfJets]/F");
    tree->Branch("L1JetPhi",l1PhiJet,"L1JetPhi[L1NumberOfJets]/F");
    tree->Branch("L1JetPt",l1PtJet,"L1JetPt[L1NumberOfJets]/F");
    tree->Branch("L1JetPx",l1PxJet,"L1JetPx[L1NumberOfJets]/F");
    tree->Branch("L1JetPy",l1PyJet,"L1JetPy[L1NumberOfJets]/F");
    tree->Branch("L1JetPz",l1PzJet,"L1JetPz[L1NumberOfJets]/F");
    tree->Branch("L1JetEnergy",l1EnergyJet,"L1JetEnergy[L1NumberOfJets]/F");
    
    // l2 jets
    tree->Branch("L2NumberOfJets",&l2NumberOfJets,"L2NumberOfJets/I");
    tree->Branch("L2JetEta",l2EtaJet,"L2JetEta[L2NumberOfJets]/F");
    tree->Branch("L2JetPhi",l2PhiJet,"L2JetPhi[L2NumberOfJets]/F");
    tree->Branch("L2JetPt",l2PtJet,"L2JetPt[L2NumberOfJets]/F");
    tree->Branch("L2JetPx",l2PxJet,"L2JetPx[L2NumberOfJets]/F");
    tree->Branch("L2JetPy",l2PyJet,"L2JetPy[L2NumberOfJets]/F");
    tree->Branch("L2JetPz",l2PzJet,"L2JetPz[L2NumberOfJets]/F");
    tree->Branch("L2JetEnergy",l2EnergyJet,"L2JetEnergy[L2NumberOfJets]/F");
    
    // l3 B jets
    tree->Branch("L3NumberOfBJets",&l3NumberOfBJets,"L3NumberOfBJets/I");
    tree->Branch("L3BJetEta",l3EtaBJet,"L3BJetEta[L3NumberOfBJets]/F");
    tree->Branch("L3BJetPhi",l3PhiBJet,"L3BJetPhi[L3NumberOfBJets]/F");
    tree->Branch("L3BJetPt",l3PtBJet,"L3BJetPt[L3NumberOfBJets]/F");
    tree->Branch("L3BJetPx",l3PxBJet,"L3BJetPx[L3NumberOfBJets]/F");
    tree->Branch("L3BJetPy",l3PyBJet,"L3BJetPy[L3NumberOfBJets]/F");
    tree->Branch("L3BJetPz",l3PzBJet,"L3BJetPz[L3NumberOfBJets]/F");
    tree->Branch("L3BJetEnergy",l3EnergyBJet,"L3BJetEnergy[L3NumberOfBJets]/F");
    
    
    // reconstructed jets
    tree->Branch("NumberOfJets",&numberOfJets,"NumberOfJets/I");
    tree->Branch("JetEta",etaJet,"JetEta[NumberOfJets]/F");
    tree->Branch("JetPhi",phiJet,"JetPhi[NumberOfJets]/F");
    tree->Branch("JetPt",ptJet,"JetPt[NumberOfJets]/F");
    tree->Branch("JetPx",pxJet,"JetPx[NumberOfJets]/F");
    tree->Branch("JetPy",pyJet,"JetPy[NumberOfJets]/F");
    tree->Branch("JetPz",pzJet,"JetPz[NumberOfJets]/F");
    tree->Branch("JetEnergy",energyJet,"JetEnergy[NumberOfJets]/F");
    
    tree->Branch("NumberOfTracksJet",numberOfTracksJet,"NumberOfTracksJet[NumberOfJets]/I");
    tree->Branch("JetTrackPx",trackPxJet,"JetTrackPx[NumberOfJets][3]/F");
    tree->Branch("JetTrackPy",trackPyJet,"JetTrackPy[NumberOfJets][3]/F");
    tree->Branch("JetTrackPz",trackPzJet,"JetTrackPz[NumberOfJets][3]/F");
    
    tree->Branch("nJetConstituents",numberOfConstituentsInJet,"nJetConstituents[NumberOfJets]/I");
    tree->Branch("nJetChargedConstituents",numberOfChargedConstituentsInJet,"nJetChargedConstituents[NumberOfJets]/I");
    tree->Branch("neutralHadronEnergyFraction",neutralHadronEnergyFraction,"neutralHadronEnergyFraction[NumberOfJets]/F");
    tree->Branch("photonEnergyFraction",photonEnergyFraction,"photonEnergyFraction[NumberOfJets]/F");
    tree->Branch("chargedHadronEnergyFraction",chargedHadronEnergyFraction,"chargedHadronEnergyFraction[NumberOfJets]/F");
    tree->Branch("electronEnergyFraction",electronEnergyFraction,"electronEnergyFraction[NumberOfJets]/F");
    
    tree->Branch("JetIDLoose",jetIDLoose,"JetIDLoose[NumberOfJets]/O");
    tree->Branch("JetIDMedium",jetIDMedium,"JetIDMedium[NumberOfJets]/O");
    tree->Branch("JetIDTight",jetIDTight,"JetIDTight[NumberOfJets]/O");
    
    tree->Branch("JetSvNumTracks",nSvTracksJet,"JetSvNumTracks[NumberOfJets]/I");
    tree->Branch("JetSvMass",svMassJet,"JetSvMass[NumberOfJets]/F");
    tree->Branch("JetSvFDSig",svFDsigJet,"JetSvFDSig[NumberOfJets]/F");
    
    if (_doMC) {
        tree->Branch("IsJetMatchedPartonExist",isJetMatchedPartonExist,"IsJetMatchedPartonExist[NumberOfJets]/O");
        tree->Branch("IsJetFromHiggs",isJetMatchedPartonFromHiggsBosonDecay,"IsJetFromHiggs[NumberOfJets]/O");
        tree->Branch("MatchedPartonFlavor",flavorJetMatchedParton,"MatchedPartonFlavor[NumberOfJets]/I");
	tree->Branch("PartonFlavorJet",partonFlavorJet,"PartonFlavorJet[NumberOfJets]/I");
	tree->Branch("HflContentJet", hflContentJet,"HflContentJet[NumberOfJets]/I");
    }
    // B-Tag Discriminants -->
    tree->Branch("JetBProbBJetTag",jetBProbBJetTag,"jetBProbBJetTag[NumberOfJets]/F");
    tree->Branch("JetProbBJetTag",jetProbBJetTag,"jetProbBJetTag[NumberOfJets]/F");
    tree->Branch("TCHPBJetTag",tcHPBJetTag,"TCHPBJetTag[NumberOfJets]/F");
    tree->Branch("TCHEBJetTag",tcHEBJetTag,"TCHEBJetTag[NumberOfJets]/F");
    if (_doNegativeBtags) {
      std::cout << "Branches for negative btags added" << std::endl;
      tree->Branch("nTCHPBJetTag",ntcHPBJetTag,"nTCHPBJetTag[NumberOfJets]/F");
      tree->Branch("nTCHEBJetTag",ntcHEBJetTag,"nTCHEBJetTag[NumberOfJets]/F");
    }
    tree->Branch("SVHPBJetTag",svHPBJetTag,"SVHPBJetTag[NumberOfJets]/F");
    tree->Branch("SVHEBJetTag",svHEBJetTag,"SVHEBJetTag[NumberOfJets]/F");
    tree->Branch("CombSVBJetTag",combSVBJetTag,"CombSVBJetTag[NumberOfJets]/F");
    tree->Branch("CombSVMVABJetTag",combSVMVABJetTag,"CombSVMVABJetTag[NumberOfJets]/F");
    tree->Branch("TCMinus2ndBJetTag",tcMinus2ndBJetTag,"TCMinus2ndBJetTag[NumberOfJets]/F");
    tree->Branch("TCMinus3rdBJetTag",tcMinus3rdBJetTag,"TCMinus3rdBJetTag[NumberOfJets]/F");
    tree->Branch("TCNumberOfSelectedTracks",tcNumberOfSelectedTracks,"TCNumberOfSelectedTracks[NumberOfJets]/I");
    
    // HLT btag match information
    //tree->Branch("IsJetWithHltBtag",isJetWithHltBtag,"isJetWithHltBtag[NumberOfJets]/O");
    tree->Branch("IsJetWithHltBtagBitPattern",isJetWithHltBtagBitPattern,"isJetWithHltBtagBitPattern[NumberOfJets]/i");
    tree->Branch("DeltaRWithHltBtag",deltaRWithHltBtag,"deltaRWithHltBtag[NumberOfJets]/F");
    tree->Branch("DeltaPtWithHltBtag",deltaPtWithHltBtag,"deltaPtWithHltBtag[NumberOfJets]/F");
    
    // L1 match information
    tree->Branch("IsJetWithL1Jet",isJetWithL1Jet,"isJetWithL1Jet[NumberOfJets]/O");
    tree->Branch("DeltaRWithL1Jet",deltaRWithL1Jet,"deltaRWithL1Jet[NumberOfJets]/F");
    tree->Branch("DeltaPtWithL1Jet",deltaPtWithL1Jet,"deltaPtWithL1Jet[NumberOfJets]/F");
    
    // L2 match information
    //tree->Branch("IsJetWithL2Jet",isJetWithL2Jet,"isJetWithL2Jet[NumberOfJets]/O");
    tree->Branch("IsJetWithL2JetBitPattern",isJetWithL2JetBitPattern,"isJetWithL2JetBitPattern[NumberOfJets]/i");
    tree->Branch("DeltaRWithL2Jet",deltaRWithL2Jet,"deltaRWithL2Jet[NumberOfJets]/F");
    tree->Branch("DeltaPtWithL2Jet",deltaPtWithL2Jet,"deltaPtWithL2Jet[NumberOfJets]/F");
    
    // reconstructed muons --->
    
    tree->Branch("NumberOfMuons",&numberOfMuons,"NumberOfMuons/I");
    tree->Branch("MuonCharge",chargeMuon,"MuonCharge[NumberOfMuons]/I");
    tree->Branch("IsTrackerMuon",isTrackerMuon,"IsTrackerMuon[NumberOfMuons]/O");
    tree->Branch("IsGlobalMuon",isGlobalMuon,"IsGlobalMuon[NumberOfMuons]/O");
    tree->Branch("PFChHadronIsoMuon",pfChargedHadronIsoMuon,"PFChHadronIsoMuon[NumberOfMuons]/F");
    tree->Branch("PFNeHadronIsoMuon",pfNeutralHadronIsoMuon,"PFNeHadronIsoMuon[NumberOfMuons]/F");
    tree->Branch("PFGammaIsoMuon",pfGammaIsoMuon,"PFGammaIsoMuon[NumberOfMuons]/F");
    tree->Branch("IsMuonInJet",isMuonInJet,"IsMuonInJet[NumberOfMuons]/O");
    tree->Branch("MuonAssociatedJetIndex",associatedJetIndex,"MuonAssociatedJetIndex[NumberOfMuons]/I");
    tree->Branch("MuonEta",etaMuon,"MuonEta[NumberOfMuons]/F");
    tree->Branch("MuonPhi",phiMuon,"MuonPhi[NumberOfMuons]/F");
    tree->Branch("MuonPt",ptMuon,"MuonPt[NumberOfMuons]/F");
    tree->Branch("MuonPx",pxMuon,"MuonPx[NumberOfMuons]/F");
    tree->Branch("MuonPy",pyMuon,"MuonPy[NumberOfMuons]/F");
    tree->Branch("MuonPz",pzMuon,"MuonPz[NumberOfMuons]/F");
    tree->Branch("TrigPrescale",trigPrescale,"trigPrescale[25]/F");
    
  // l1 muons
  tree->Branch("L1NumberOfMuons",&l1NumberOfMuons,"L1NumberOfMuons/I");
  tree->Branch("L1MuonEta",l1EtaMuon,"L1MuonEta[L1NumberOfMuons]/F");
  tree->Branch("L1MuonPhi",l1PhiMuon,"L1MuonPhi[L1NumberOfMuons]/F");
  tree->Branch("L1MuonPt",l1PtMuon,"L1MuonPt[L1NumberOfMuons]/F");
  tree->Branch("L1MuonPx",l1PxMuon,"L1MuonPx[L1NumberOfMuons]/F");
  tree->Branch("L1MuonPy",l1PyMuon,"L1MuonPy[L1NumberOfMuons]/F");
  tree->Branch("L1MuonPz",l1PzMuon,"L1MuonPz[L1NumberOfMuons]/F");
  tree->Branch("L1MuonEnergy",l1EnergyMuon,"L1MuonEnergy[L1NumberOfMuons]/F");
  tree->Branch("L1MuonCharge",l1ChargeMuon,"L1MuonCharge[L1NumberOfMuons]/I");

  // hlt muons
  tree->Branch("HltNumberOfMuons",&hltNumberOfMuons,"HltNumberOfMuons/I");
  tree->Branch("HltMuonEta",hltEtaMuon,"HltMuonEta[HltNumberOfMuons]/F");
  tree->Branch("HltMuonPhi",hltPhiMuon,"HltMuonPhi[HltNumberOfMuons]/F");
  tree->Branch("HltMuonPt",hltPtMuon,"HltMuonPt[HltNumberOfMuons]/F");
  tree->Branch("HltMuonPx",hltPxMuon,"HltMuonPx[HltNumberOfMuons]/F");
  tree->Branch("HltMuonPy",hltPyMuon,"HltMuonPy[HltNumberOfMuons]/F");
  tree->Branch("HltMuonPz",hltPzMuon,"HltMuonPz[HltNumberOfMuons]/F");
  tree->Branch("HltMuonEnergy",hltEnergyMuon,"HltMuonEnergy[HltNumberOfMuons]/F");
  tree->Branch("HltMuonCharge",hltChargeMuon,"HltMuonCharge[HltNumberOfMuons]/I");
    
    if (_doMC) {
        tree->Branch("IsMatchedGenMuonExist",isMatchedGeneratorMuonExist,"isMatchedGenMuonExist[NumberOfMuons]/O");
        tree->Branch("ChargeGenMuon",chargeGenMuon,"ChargeGenMuon[NumberOfMuons]/I");
        tree->Branch("GenMuonPx",pxGenMuon,"GenMuonPx[NumberOfMuons]/F");
        tree->Branch("GenMuonPy",pyGenMuon,"GenMuonPy[NumberOfMuons]/F");
        tree->Branch("GenMuonPz",pzGenMuon,"GenMuonPz[NumberOfMuons]/F");
        tree->Branch("GenMuonMotherPDG",pdgMotherGenMuon,"GenMuonMotherPDG[NumberOfMuons]/I");
        // Gen jets
        tree->Branch("GenNumberOfJets",&genNumberOfJets,"GenNumberOfJets/I");
        tree->Branch("GenJetEta",genEtaJet,"GenJetEta[GenNumberOfJets]/F");
        tree->Branch("GenJetPhi",genPhiJet,"GenJetPhi[GenNumberOfJets]/F");
        tree->Branch("GenJetPt",genPtJet,"GenJetPt[GenNumberOfJets]/F");
        tree->Branch("GenJetPx",genPxJet,"GenJetPx[GenNumberOfJets]/F");
        tree->Branch("GenJetPy",genPyJet,"GenJetPy[GenNumberOfJets]/F");
        tree->Branch("GenJetPz",genPzJet,"GenJetPz[GenNumberOfJets]/F");
        tree->Branch("GenJetEnergy",genEnergyJet,"GenJetEnergy[GenNumberOfJets]/F");
    }
    
    if (_doMC) {
        // generator information (cross-sections)
        
        treeGenInfo=fs->make<TTree>("GenInfo","GenInfo");
        
        treeGenInfo->Branch("internalXsec",&_internalXsec,"internalXsec/F");
        treeGenInfo->Branch("internalXsecE",&_internalXsecE,"internalXsecE/F");
        
        treeGenInfo->Branch("externalXsecLO",&_externalXsecLO,"externalXsecLO/F");
        treeGenInfo->Branch("externalXsecLOE",&_externalXsecLOE,"externalXsecLOE/F");
        
        treeGenInfo->Branch("externalXsecNLO",&_externalXsecNLO,"externalXsecNLO/F");
        treeGenInfo->Branch("externalXsecNLOE",&_externalXsecNLOE,"externalXsecNLOE/F");
        
        treeGenInfo->Branch("filterEfficiency",&_filterEfficiency,"filterEfficiency/F");
        
        // information on Higgs boson kinematics and its decay products
        
        tree->Branch("PdgHiggsBoson",&pdgHiggsBoson,"PdgHiggsBoson/I");
        tree->Branch("EtaHiggsBoson",&etaHiggsBoson,"EtaHiggsBoson/F");
        tree->Branch("PhiHiggsBoson",&phiHiggsBoson,"PhiHiggsBoson/F");
        tree->Branch("PtHiggsBoson",&ptHiggsBoson,"PtHiggsBoson/F");
        tree->Branch("PxHiggsBoson",&pxHiggsBoson,"PxHiggsBoson/F");
        tree->Branch("PyHiggsBoson",&pyHiggsBoson,"PyHiggsBoson/F");
        tree->Branch("PzHiggsBoson",&pzHiggsBoson,"PzHiggsBoson/F");
        tree->Branch("MassHiggsBoson",&massHiggsBoson,"MassHiggsBoson/F");
        tree->Branch("EnergyHiggsBoson",&energyHiggsBoson,"EnergyHiggsBoson/F");
        
        tree->Branch("NumberOfPartons",&numberOfPartons,"NumberOfPartons/I");
        tree->Branch("PartonPdg",pdgParton,"PartonPdg[NumberOfPartons]/I");
        tree->Branch("PartonPx",pxParton,"PartonPx[NumberOfPartons]/F");
        tree->Branch("PartonPy",pyParton,"PartonPy[NumberOfPartons]/F");
        tree->Branch("PartonPz",pzParton,"PartonPz[NumberOfPartons]/F");
        tree->Branch("PartonEnergy",energyParton,"PartonEnergy[NumberOfPartons]/F");
        tree->Branch("PartonEta",etaParton,"PartonEta[NumberOfPartons]/F");
        tree->Branch("PartonPhi",phiParton,"PartonPhi[NumberOfPartons]/F");
        tree->Branch("PartonPt",ptParton,"PartonPt[NumberOfPartons]/F");
        
        
    }
  
    // missing energy and momentum
    tree -> Branch( "met", &met, "met/F" );
    tree -> Branch( "missingPx", &missingPx, "missingPx/F" );
    tree -> Branch( "missingPy", &missingPy, "missingPy/F" );
    
    
    std::cout << "Trees are created..." << std::endl;
    std::cout << std::endl;
    
    std::cout << "Trigger names are stored in gtlHist" << std::endl;
    
    _events = 0;
    _selectedEvents = 0;
    
    firstEvent = true;
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HbbMSSMAnalysis::endJob() {
    
    std::cout << std::endl;
    std::cout << " HbbMSSMAnalysis::endJob() --> " << std::endl;
    std::cout << " Number of processed events = " << _events << std::endl;
    std::cout << " Number of selected events  = " << _selectedEvents << std::endl;
    std::cout << std::endl;
    
    tree->GetDirectory()->cd();
    tree->Write();
    
    if (_doMC)  { 
        treeGenInfo->GetDirectory()->cd();
        treeGenInfo->Write();
    }
    
    // trigger count statistics
    std::cout << "=======================================================" << std::endl;
    std::cout << "Statistics of all trigger paths encountered through all events:" << std::endl << std::endl;
    for (std::map<std::string,int>::iterator mit=triggerCounter.begin(); mit != triggerCounter.end(); ++mit) {
        std::cout << "Trigger " << std::setw(40) << mit->first << " count " << std::setw(8) << mit->second << std::endl;
    }
    
    // trigger object count statistics
    std::cout << "=======================================================" << std::endl;
    std::cout << "Statistics of all trigger objects encountered through all events:" << std::endl << std::endl;
    for (std::map<std::string,int>::iterator mit=triggerObjectCounter.begin(); mit != triggerObjectCounter.end(); ++mit) {
        std::cout << "Trigger object  " << std::setw(15) << mit->first << " count " << std::setw(8) << mit->second << std::endl;
    }
    
    
    std::cout << "=======================================================" << std::endl;
    // print generic trigger list 
    std::cout << "List of generic triggers for this analysis:" << std::endl << std::endl;
    for (std::vector<std::string>::iterator gtlIt = _genericTriggerList.begin();
         gtlIt != _genericTriggerList.end(); ++gtlIt) {
        std::cout << std::setw(3) << gtlIt - _genericTriggerList.begin()
        << "  " << *gtlIt << std::endl;
    }
    std::cout << "=======================================================" << std::endl;
}

bool HbbMSSMAnalysis::IsJetMatchedToHLTBtag( const pat::JetRef & jet, std::vector<reco::Particle> HLTBtag , const double DR , 
                                            float & theDeltaRWithHltBtag, float & theDeltaPtWithHltBtag) {
    unsigned int dim =  HLTBtag.size();
    unsigned int nPass=0;
    if (dim==0) return false;
    for (unsigned int k =0; k< dim; k++ ) {
        if (  (deltaR(HLTBtag[k], *jet) < DR) ){     
            nPass++ ;
            theDeltaRWithHltBtag = deltaR(HLTBtag[k], *jet);
            theDeltaPtWithHltBtag = jet->pt() - HLTBtag[k].pt();
        }
    }
    return (nPass>0);
}

bool HbbMSSMAnalysis::IsJetMatchedToL1Jet( const pat::JetRef & jet, const l1extra::L1JetParticleCollection & l1jets, const double DR , 
                                          float & theDeltaR, float & theDeltaPt) {
    unsigned int dim =  l1jets.size();
    unsigned int nPass=0;
    if (dim==0) return false;
    if (dim>100) dim = 100;
    for (unsigned int k =0; k< dim; k++ ) {
        if (  (deltaR(l1jets[k], *jet) < DR) ){     
            nPass++ ;
            theDeltaR = deltaR(l1jets[k], *jet);
            theDeltaPt = jet->pt() - l1jets[k].pt();
        }
    }
    return (nPass>0);
}

int HbbMSSMAnalysis::jetHflContent( const pat::JetRef & jet, const edm::Handle<reco::GenParticleCollection> & theGenParticles ) {
  // return -4, -5, +4 or +5 if if b or c parton matched to jet
  // return zero otherwise
  int numberOfGenParticles = theGenParticles->size();
  for (int iP=0; iP<numberOfGenParticles; iP++) {
    reco::GenParticleRef particle(theGenParticles,iP);
    int partPdg = particle->pdgId();
    if ( (abs(partPdg) == 4) || (abs(partPdg) == 5) ) {
      if (deltaR( *particle, *jet) < 0.3) {
	return partPdg;
      }
    }
  }
  return 0;
}
  


//define this as a plug-in
DEFINE_FWK_MODULE(HbbMSSMAnalysis);
