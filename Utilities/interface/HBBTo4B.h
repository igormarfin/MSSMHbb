//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 13 23:01:08 2012 by ROOT version 5.27/06b
// from TTree HBBTo4B/HBBTo4B
// found on file: /scratch/hh/lustre/cms/user/walsh/HbbMSSMAnalysis/CMSSW424p1/v1/data/Run2011A/PromptReco-v4/HbbMSSMAnalysis_10_1_luY.root
//////////////////////////////////////////////////////////

#ifndef HBBTo4B_h
#define HBBTo4B_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class HBBTo4B {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Event;
   Int_t           Lumi;
   Int_t           NumberOfPU;
   Int_t           NumberOfPUInTime;
   Int_t           trgAccept;
   Int_t           NumberOfPV;
   Float_t         ProbPV[20];   //[NumberOfPV]
   Float_t         NdofPV[20];   //[NumberOfPV]
   Int_t           NTrkPV[20];   //[NumberOfPV]
   Float_t         Chi2PV[20];   //[NumberOfPV]
   Float_t         XPV[20];   //[NumberOfPV]
   Float_t         YPV[20];   //[NumberOfPV]
   Float_t         ZPV[20];   //[NumberOfPV]
   Float_t         SumTrkPtPV[20];   //[NumberOfPV]
   Int_t           L1NumberOfJets;
   Float_t         L1JetEta[8];   //[L1NumberOfJets]
   Float_t         L1JetPhi[8];   //[L1NumberOfJets]
   Float_t         L1JetPt[8];   //[L1NumberOfJets]
   Float_t         L1JetPx[8];   //[L1NumberOfJets]
   Float_t         L1JetPy[8];   //[L1NumberOfJets]
   Float_t         L1JetPz[8];   //[L1NumberOfJets]
   Float_t         L1JetEnergy[8];   //[L1NumberOfJets]
   Int_t           L2NumberOfJets;
   Float_t         L2JetEta[1];   //[L2NumberOfJets]
   Float_t         L2JetPhi[1];   //[L2NumberOfJets]
   Float_t         L2JetPt[1];   //[L2NumberOfJets]
   Float_t         L2JetPx[1];   //[L2NumberOfJets]
   Float_t         L2JetPy[1];   //[L2NumberOfJets]
   Float_t         L2JetPz[1];   //[L2NumberOfJets]
   Float_t         L2JetEnergy[1];   //[L2NumberOfJets]
   Int_t           L3NumberOfBJets;
   Float_t         L3BJetEta[1];   //[L3NumberOfBJets]
   Float_t         L3BJetPhi[1];   //[L3NumberOfBJets]
   Float_t         L3BJetPt[1];   //[L3NumberOfBJets]
   Float_t         L3BJetPx[1];   //[L3NumberOfBJets]
   Float_t         L3BJetPy[1];   //[L3NumberOfBJets]
   Float_t         L3BJetPz[1];   //[L3NumberOfBJets]
   Float_t         L3BJetEnergy[1];   //[L3NumberOfBJets]
   Int_t           NumberOfJets;
   Float_t         JetEta[12];   //[NumberOfJets]
   Float_t         JetPhi[12];   //[NumberOfJets]
   Float_t         JetPt[12];   //[NumberOfJets]
   Float_t         JetPx[12];   //[NumberOfJets]
   Float_t         JetPy[12];   //[NumberOfJets]
   Float_t         JetPz[12];   //[NumberOfJets]
   Float_t         JetEnergy[12];   //[NumberOfJets]
   Int_t           NumberOfTracksJet[12];   //[NumberOfJets]
   Float_t         JetTrackPx[12][3];   //[NumberOfJets]
   Float_t         JetTrackPy[12][3];   //[NumberOfJets]
   Float_t         JetTrackPz[12][3];   //[NumberOfJets]
   Int_t           nJetConstituents[12];   //[NumberOfJets]
   Int_t           nJetChargedConstituents[12];   //[NumberOfJets]
   Float_t         neutralHadronEnergyFraction[12];   //[NumberOfJets]
   Float_t         photonEnergyFraction[12];   //[NumberOfJets]
   Float_t         chargedHadronEnergyFraction[12];   //[NumberOfJets]
   Float_t         electronEnergyFraction[12];   //[NumberOfJets]
   Bool_t          JetIDLoose[12];   //[NumberOfJets]
   Bool_t          JetIDMedium[12];   //[NumberOfJets]
   Bool_t          JetIDTight[12];   //[NumberOfJets]
   Int_t           JetSvNumTracks[12];   //[NumberOfJets]
   Float_t         JetSvMass[12];   //[NumberOfJets]
   Float_t         JetSvFDSig[12];   //[NumberOfJets]
   Float_t         JetBProbBJetTag[12];   //[NumberOfJets]
   Float_t         JetProbBJetTag[12];   //[NumberOfJets]
   Float_t         TCHPBJetTag[12];   //[NumberOfJets]
   Float_t         TCHEBJetTag[12];   //[NumberOfJets]
   Float_t         nTCHPBJetTag[12];   //[NumberOfJets]
   Float_t         nTCHEBJetTag[12];   //[NumberOfJets]
   Float_t         SVHPBJetTag[12];   //[NumberOfJets]
   Float_t         SVHEBJetTag[12];   //[NumberOfJets]
   Float_t         CombSVBJetTag[12];   //[NumberOfJets]
   Float_t         CombSVMVABJetTag[12];   //[NumberOfJets]
   Float_t         TCMinus2ndBJetTag[12];   //[NumberOfJets]
   Float_t         TCMinus3rdBJetTag[12];   //[NumberOfJets]
   Int_t           TCNumberOfSelectedTracks[12];   //[NumberOfJets]
   Bool_t          IsJetWithHltBtag[12];   //[NumberOfJets]
   Float_t         DeltaRWithHltBtag[12];   //[NumberOfJets]
   Float_t         DeltaPtWithHltBtag[12];   //[NumberOfJets]
   Bool_t          IsJetWithL1Jet[12];   //[NumberOfJets]
   Float_t         DeltaRWithL1Jet[12];   //[NumberOfJets]
   Float_t         DeltaPtWithL1Jet[12];   //[NumberOfJets]
   Bool_t          IsJetWithL2Jet[12];   //[NumberOfJets]
   Float_t         DeltaRWithL2Jet[12];   //[NumberOfJets]
   Float_t         DeltaPtWithL2Jet[12];   //[NumberOfJets]
   Int_t           NumberOfMuons;
   Int_t           MuonCharge[3];   //[NumberOfMuons]
   Bool_t          IsTrackerMuon[3];   //[NumberOfMuons]
   Bool_t          IsGlobalMuon[3];   //[NumberOfMuons]
   Float_t         PFChHadronIsoMuon[3];   //[NumberOfMuons]
   Float_t         PFNeHadronIsoMuon[3];   //[NumberOfMuons]
   Float_t         PFGammaIsoMuon[3];   //[NumberOfMuons]
   Bool_t          IsMuonInJet[3];   //[NumberOfMuons]
   Int_t           MuonAssociatedJetIndex[3];   //[NumberOfMuons]
   Float_t         MuonEta[3];   //[NumberOfMuons]
   Float_t         MuonPhi[3];   //[NumberOfMuons]
   Float_t         MuonPt[3];   //[NumberOfMuons]
   Float_t         MuonPx[3];   //[NumberOfMuons]
   Float_t         MuonPy[3];   //[NumberOfMuons]
   Float_t         MuonPz[3];   //[NumberOfMuons]
   Float_t         TrigPrescale[25];
   Int_t           L1NumberOfMuons;
   Float_t         L1MuonEta[4];   //[L1NumberOfMuons]
   Float_t         L1MuonPhi[4];   //[L1NumberOfMuons]
   Float_t         L1MuonPt[4];   //[L1NumberOfMuons]
   Float_t         L1MuonPx[4];   //[L1NumberOfMuons]
   Float_t         L1MuonPy[4];   //[L1NumberOfMuons]
   Float_t         L1MuonPz[4];   //[L1NumberOfMuons]
   Float_t         L1MuonEnergy[4];   //[L1NumberOfMuons]
   Int_t           L1MuonCharge[4];   //[L1NumberOfMuons]
   Int_t           HltNumberOfMuons;
   Float_t         HltMuonEta[1];   //[HltNumberOfMuons]
   Float_t         HltMuonPhi[1];   //[HltNumberOfMuons]
   Float_t         HltMuonPt[1];   //[HltNumberOfMuons]
   Float_t         HltMuonPx[1];   //[HltNumberOfMuons]
   Float_t         HltMuonPy[1];   //[HltNumberOfMuons]
   Float_t         HltMuonPz[1];   //[HltNumberOfMuons]
   Float_t         HltMuonEnergy[1];   //[HltNumberOfMuons]
   Int_t           HltMuonCharge[1];   //[HltNumberOfMuons]

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Lumi;   //!
   TBranch        *b_NumberOfPU;   //!
   TBranch        *b_NumberOfPUInTime;   //!
   TBranch        *b_trgAccept;   //!
   TBranch        *b_NumberOfPV;   //!
   TBranch        *b_ProbPV;   //!
   TBranch        *b_NdofPV;   //!
   TBranch        *b_NTrkPV;   //!
   TBranch        *b_Chi2PV;   //!
   TBranch        *b_XPV;   //!
   TBranch        *b_YPV;   //!
   TBranch        *b_ZPV;   //!
   TBranch        *b_SumTrkPtPV;   //!
   TBranch        *b_L1NumberOfJets;   //!
   TBranch        *b_L1JetEta;   //!
   TBranch        *b_L1JetPhi;   //!
   TBranch        *b_L1JetPt;   //!
   TBranch        *b_L1JetPx;   //!
   TBranch        *b_L1JetPy;   //!
   TBranch        *b_L1JetPz;   //!
   TBranch        *b_L1JetEnergy;   //!
   TBranch        *b_L2NumberOfJets;   //!
   TBranch        *b_L2JetEta;   //!
   TBranch        *b_L2JetPhi;   //!
   TBranch        *b_L2JetPt;   //!
   TBranch        *b_L2JetPx;   //!
   TBranch        *b_L2JetPy;   //!
   TBranch        *b_L2JetPz;   //!
   TBranch        *b_L2JetEnergy;   //!
   TBranch        *b_L3NumberOfBJets;   //!
   TBranch        *b_L3BJetEta;   //!
   TBranch        *b_L3BJetPhi;   //!
   TBranch        *b_L3BJetPt;   //!
   TBranch        *b_L3BJetPx;   //!
   TBranch        *b_L3BJetPy;   //!
   TBranch        *b_L3BJetPz;   //!
   TBranch        *b_L3BJetEnergy;   //!
   TBranch        *b_NumberOfJets;   //!
   TBranch        *b_JetEta;   //!
   TBranch        *b_JetPhi;   //!
   TBranch        *b_JetPt;   //!
   TBranch        *b_JetPx;   //!
   TBranch        *b_JetPy;   //!
   TBranch        *b_JetPz;   //!
   TBranch        *b_JetEnergy;   //!
   TBranch        *b_NumberOfTracksJet;   //!
   TBranch        *b_JetTrackPx;   //!
   TBranch        *b_JetTrackPy;   //!
   TBranch        *b_JetTrackPz;   //!
   TBranch        *b_nJetConstituents;   //!
   TBranch        *b_nJetChargedConstituents;   //!
   TBranch        *b_neutralHadronEnergyFraction;   //!
   TBranch        *b_photonEnergyFraction;   //!
   TBranch        *b_chargedHadronEnergyFraction;   //!
   TBranch        *b_electronEnergyFraction;   //!
   TBranch        *b_JetIDLoose;   //!
   TBranch        *b_JetIDMedium;   //!
   TBranch        *b_JetIDTight;   //!
   TBranch        *b_JetSvNumTracks;   //!
   TBranch        *b_JetSvMass;   //!
   TBranch        *b_JetSvFDSig;   //!
   TBranch        *b_JetBProbBJetTag;   //!
   TBranch        *b_JetProbBJetTag;   //!
   TBranch        *b_TCHPBJetTag;   //!
   TBranch        *b_TCHEBJetTag;   //!
   TBranch        *b_nTCHPBJetTag;   //!
   TBranch        *b_nTCHEBJetTag;   //!
   TBranch        *b_SVHPBJetTag;   //!
   TBranch        *b_SVHEBJetTag;   //!
   TBranch        *b_CombSVBJetTag;   //!
   TBranch        *b_CombSVMVABJetTag;   //!
   TBranch        *b_TCMinus2ndBJetTag;   //!
   TBranch        *b_TCMinus3rdBJetTag;   //!
   TBranch        *b_TCNumberOfSelectedTracks;   //!
   TBranch        *b_IsJetWithHltBtag;   //!
   TBranch        *b_DeltaRWithHltBtag;   //!
   TBranch        *b_DeltaPtWithHltBtag;   //!
   TBranch        *b_IsJetWithL1Jet;   //!
   TBranch        *b_DeltaRWithL1Jet;   //!
   TBranch        *b_DeltaPtWithL1Jet;   //!
   TBranch        *b_IsJetWithL2Jet;   //!
   TBranch        *b_DeltaRWithL2Jet;   //!
   TBranch        *b_DeltaPtWithL2Jet;   //!
   TBranch        *b_NumberOfMuons;   //!
   TBranch        *b_MuonCharge;   //!
   TBranch        *b_IsTrackerMuon;   //!
   TBranch        *b_IsGlobalMuon;   //!
   TBranch        *b_PFChHadronIsoMuon;   //!
   TBranch        *b_PFNeHadronIsoMuon;   //!
   TBranch        *b_PFGammaIsoMuon;   //!
   TBranch        *b_IsMuonInJet;   //!
   TBranch        *b_MuonAssociatedJetIndex;   //!
   TBranch        *b_MuonEta;   //!
   TBranch        *b_MuonPhi;   //!
   TBranch        *b_MuonPt;   //!
   TBranch        *b_MuonPx;   //!
   TBranch        *b_MuonPy;   //!
   TBranch        *b_MuonPz;   //!
   TBranch        *b_trigPrescale;   //!
   TBranch        *b_L1NumberOfMuons;   //!
   TBranch        *b_L1MuonEta;   //!
   TBranch        *b_L1MuonPhi;   //!
   TBranch        *b_L1MuonPt;   //!
   TBranch        *b_L1MuonPx;   //!
   TBranch        *b_L1MuonPy;   //!
   TBranch        *b_L1MuonPz;   //!
   TBranch        *b_L1MuonEnergy;   //!
   TBranch        *b_L1MuonCharge;   //!
   TBranch        *b_HltNumberOfMuons;   //!
   TBranch        *b_HltMuonEta;   //!
   TBranch        *b_HltMuonPhi;   //!
   TBranch        *b_HltMuonPt;   //!
   TBranch        *b_HltMuonPx;   //!
   TBranch        *b_HltMuonPy;   //!
   TBranch        *b_HltMuonPz;   //!
   TBranch        *b_HltMuonEnergy;   //!
   TBranch        *b_HltMuonCharge;   //!

   HBBTo4B(TTree *tree=0);
   virtual ~HBBTo4B();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif // #ifdef 
