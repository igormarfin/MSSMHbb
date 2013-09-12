#include "Analysis/Utilities/interface/HBBTo4B_MC.h"
HBBTo4B_MC::HBBTo4B_MC(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/hh/current/cms/user/walsh/HbbMSSMAnalysis/CMSSW424p1/v1/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/HbbMSSMAnalysis_10_1_FPu.root");
      if (!f) {
         f = new TFile("/scratch/hh/current/cms/user/walsh/HbbMSSMAnalysis/CMSSW424p1/v1/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/HbbMSSMAnalysis_10_1_FPu.root");
         f->cd("/scratch/hh/current/cms/user/walsh/HbbMSSMAnalysis/CMSSW424p1/v1/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/HbbMSSMAnalysis_10_1_FPu.root:/hbbanalysis");
      }
      tree = (TTree*)gDirectory->Get("HBBTo4B");

   }
   Init(tree);
}

HBBTo4B_MC::~HBBTo4B_MC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t HBBTo4B_MC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t HBBTo4B_MC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void HBBTo4B_MC::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_Lumi);
   fChain->SetBranchAddress("NumberOfPU", &NumberOfPU, &b_NumberOfPU);
   fChain->SetBranchAddress("NumberOfPUInTime", &NumberOfPUInTime, &b_NumberOfPUInTime);
   fChain->SetBranchAddress("trgAccept", &trgAccept, &b_trgAccept);
   fChain->SetBranchAddress("NumberOfPV", &NumberOfPV, &b_NumberOfPV);
   fChain->SetBranchAddress("ProbPV", ProbPV, &b_ProbPV);
   fChain->SetBranchAddress("NdofPV", NdofPV, &b_NdofPV);
   fChain->SetBranchAddress("NTrkPV", NTrkPV, &b_NTrkPV);
   fChain->SetBranchAddress("Chi2PV", Chi2PV, &b_Chi2PV);
   fChain->SetBranchAddress("XPV", XPV, &b_XPV);
   fChain->SetBranchAddress("YPV", YPV, &b_YPV);
   fChain->SetBranchAddress("ZPV", ZPV, &b_ZPV);
   fChain->SetBranchAddress("SumTrkPtPV", SumTrkPtPV, &b_SumTrkPtPV);
   fChain->SetBranchAddress("L1NumberOfJets", &L1NumberOfJets, &b_L1NumberOfJets);
   fChain->SetBranchAddress("L1JetEta", L1JetEta, &b_L1JetEta);
   fChain->SetBranchAddress("L1JetPhi", L1JetPhi, &b_L1JetPhi);
   fChain->SetBranchAddress("L1JetPt", L1JetPt, &b_L1JetPt);
   fChain->SetBranchAddress("L1JetPx", L1JetPx, &b_L1JetPx);
   fChain->SetBranchAddress("L1JetPy", L1JetPy, &b_L1JetPy);
   fChain->SetBranchAddress("L1JetPz", L1JetPz, &b_L1JetPz);
   fChain->SetBranchAddress("L1JetEnergy", L1JetEnergy, &b_L1JetEnergy);
   fChain->SetBranchAddress("L2NumberOfJets", &L2NumberOfJets, &b_L2NumberOfJets);
   fChain->SetBranchAddress("L2JetEta", L2JetEta, &b_L2JetEta);
   fChain->SetBranchAddress("L2JetPhi", L2JetPhi, &b_L2JetPhi);
   fChain->SetBranchAddress("L2JetPt", L2JetPt, &b_L2JetPt);
   fChain->SetBranchAddress("L2JetPx", L2JetPx, &b_L2JetPx);
   fChain->SetBranchAddress("L2JetPy", L2JetPy, &b_L2JetPy);
   fChain->SetBranchAddress("L2JetPz", L2JetPz, &b_L2JetPz);
   fChain->SetBranchAddress("L2JetEnergy", L2JetEnergy, &b_L2JetEnergy);
   fChain->SetBranchAddress("L3NumberOfBJets", &L3NumberOfBJets, &b_L3NumberOfBJets);
   fChain->SetBranchAddress("L3BJetEta", L3BJetEta, &b_L3BJetEta);
   fChain->SetBranchAddress("L3BJetPhi", L3BJetPhi, &b_L3BJetPhi);
   fChain->SetBranchAddress("L3BJetPt", L3BJetPt, &b_L3BJetPt);
   fChain->SetBranchAddress("L3BJetPx", L3BJetPx, &b_L3BJetPx);
   fChain->SetBranchAddress("L3BJetPy", L3BJetPy, &b_L3BJetPy);
   fChain->SetBranchAddress("L3BJetPz", L3BJetPz, &b_L3BJetPz);
   fChain->SetBranchAddress("L3BJetEnergy", L3BJetEnergy, &b_L3BJetEnergy);
   fChain->SetBranchAddress("NumberOfJets", &NumberOfJets, &b_NumberOfJets);
   fChain->SetBranchAddress("JetEta", JetEta, &b_JetEta);
   fChain->SetBranchAddress("JetPhi", JetPhi, &b_JetPhi);
   fChain->SetBranchAddress("JetPt", JetPt, &b_JetPt);
   fChain->SetBranchAddress("JetPx", JetPx, &b_JetPx);
   fChain->SetBranchAddress("JetPy", JetPy, &b_JetPy);
   fChain->SetBranchAddress("JetPz", JetPz, &b_JetPz);
   fChain->SetBranchAddress("JetEnergy", JetEnergy, &b_JetEnergy);
   fChain->SetBranchAddress("NumberOfTracksJet", NumberOfTracksJet, &b_NumberOfTracksJet);
   fChain->SetBranchAddress("JetTrackPx", JetTrackPx, &b_JetTrackPx);
   fChain->SetBranchAddress("JetTrackPy", JetTrackPy, &b_JetTrackPy);
   fChain->SetBranchAddress("JetTrackPz", JetTrackPz, &b_JetTrackPz);
   fChain->SetBranchAddress("nJetConstituents", nJetConstituents, &b_nJetConstituents);
   fChain->SetBranchAddress("nJetChargedConstituents", nJetChargedConstituents, &b_nJetChargedConstituents);
   fChain->SetBranchAddress("neutralHadronEnergyFraction", neutralHadronEnergyFraction, &b_neutralHadronEnergyFraction);
   fChain->SetBranchAddress("photonEnergyFraction", photonEnergyFraction, &b_photonEnergyFraction);
   fChain->SetBranchAddress("chargedHadronEnergyFraction", chargedHadronEnergyFraction, &b_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("electronEnergyFraction", electronEnergyFraction, &b_electronEnergyFraction);
   fChain->SetBranchAddress("JetIDLoose", JetIDLoose, &b_JetIDLoose);
   fChain->SetBranchAddress("JetIDMedium", JetIDMedium, &b_JetIDMedium);
   fChain->SetBranchAddress("JetIDTight", JetIDTight, &b_JetIDTight);
   fChain->SetBranchAddress("JetSvNumTracks", JetSvNumTracks, &b_JetSvNumTracks);
   fChain->SetBranchAddress("JetSvMass", JetSvMass, &b_JetSvMass);
   fChain->SetBranchAddress("JetSvFDSig", JetSvFDSig, &b_JetSvFDSig);
   fChain->SetBranchAddress("IsJetMatchedPartonExist", IsJetMatchedPartonExist, &b_IsJetMatchedPartonExist);
   fChain->SetBranchAddress("IsJetFromHiggs", IsJetFromHiggs, &b_IsJetFromHiggs);
   fChain->SetBranchAddress("MatchedPartonFlavor", MatchedPartonFlavor, &b_MatchedPartonFlavor);
   fChain->SetBranchAddress("PartonFlavorJet", PartonFlavorJet, &b_PartonFlavorJet);
   fChain->SetBranchAddress("HflContentJet", HflContentJet, &b_HflContentJet);
   fChain->SetBranchAddress("JetBProbBJetTag", JetBProbBJetTag, &b_JetBProbBJetTag);
   fChain->SetBranchAddress("JetProbBJetTag", JetProbBJetTag, &b_JetProbBJetTag);
   fChain->SetBranchAddress("TCHPBJetTag", TCHPBJetTag, &b_TCHPBJetTag);
   fChain->SetBranchAddress("TCHEBJetTag", TCHEBJetTag, &b_TCHEBJetTag);
   fChain->SetBranchAddress("nTCHPBJetTag", nTCHPBJetTag, &b_nTCHPBJetTag);
   fChain->SetBranchAddress("nTCHEBJetTag", nTCHEBJetTag, &b_nTCHEBJetTag);
   fChain->SetBranchAddress("SVHPBJetTag", SVHPBJetTag, &b_SVHPBJetTag);
   fChain->SetBranchAddress("SVHEBJetTag", SVHEBJetTag, &b_SVHEBJetTag);
   fChain->SetBranchAddress("CombSVBJetTag", CombSVBJetTag, &b_CombSVBJetTag);
   fChain->SetBranchAddress("CombSVMVABJetTag", CombSVMVABJetTag, &b_CombSVMVABJetTag);
   fChain->SetBranchAddress("TCMinus2ndBJetTag", TCMinus2ndBJetTag, &b_TCMinus2ndBJetTag);
   fChain->SetBranchAddress("TCMinus3rdBJetTag", TCMinus3rdBJetTag, &b_TCMinus3rdBJetTag);
   fChain->SetBranchAddress("TCNumberOfSelectedTracks", TCNumberOfSelectedTracks, &b_TCNumberOfSelectedTracks);
   fChain->SetBranchAddress("IsJetWithHltBtag", IsJetWithHltBtag, &b_IsJetWithHltBtag);
   fChain->SetBranchAddress("DeltaRWithHltBtag", DeltaRWithHltBtag, &b_DeltaRWithHltBtag);
   fChain->SetBranchAddress("DeltaPtWithHltBtag", DeltaPtWithHltBtag, &b_DeltaPtWithHltBtag);
   fChain->SetBranchAddress("IsJetWithL1Jet", IsJetWithL1Jet, &b_IsJetWithL1Jet);
   fChain->SetBranchAddress("DeltaRWithL1Jet", DeltaRWithL1Jet, &b_DeltaRWithL1Jet);
   fChain->SetBranchAddress("DeltaPtWithL1Jet", DeltaPtWithL1Jet, &b_DeltaPtWithL1Jet);
   fChain->SetBranchAddress("IsJetWithL2Jet", IsJetWithL2Jet, &b_IsJetWithL2Jet);
   fChain->SetBranchAddress("DeltaRWithL2Jet", DeltaRWithL2Jet, &b_DeltaRWithL2Jet);
   fChain->SetBranchAddress("DeltaPtWithL2Jet", DeltaPtWithL2Jet, &b_DeltaPtWithL2Jet);
   fChain->SetBranchAddress("NumberOfMuons", &NumberOfMuons, &b_NumberOfMuons);
   fChain->SetBranchAddress("MuonCharge", MuonCharge, &b_MuonCharge);
   fChain->SetBranchAddress("IsTrackerMuon", IsTrackerMuon, &b_IsTrackerMuon);
   fChain->SetBranchAddress("IsGlobalMuon", IsGlobalMuon, &b_IsGlobalMuon);
   fChain->SetBranchAddress("PFChHadronIsoMuon", PFChHadronIsoMuon, &b_PFChHadronIsoMuon);
   fChain->SetBranchAddress("PFNeHadronIsoMuon", PFNeHadronIsoMuon, &b_PFNeHadronIsoMuon);
   fChain->SetBranchAddress("PFGammaIsoMuon", PFGammaIsoMuon, &b_PFGammaIsoMuon);
   fChain->SetBranchAddress("IsMuonInJet", IsMuonInJet, &b_IsMuonInJet);
   fChain->SetBranchAddress("MuonAssociatedJetIndex", MuonAssociatedJetIndex, &b_MuonAssociatedJetIndex);
   fChain->SetBranchAddress("MuonEta", MuonEta, &b_MuonEta);
   fChain->SetBranchAddress("MuonPhi", MuonPhi, &b_MuonPhi);
   fChain->SetBranchAddress("MuonPt", MuonPt, &b_MuonPt);
   fChain->SetBranchAddress("MuonPx", MuonPx, &b_MuonPx);
   fChain->SetBranchAddress("MuonPy", MuonPy, &b_MuonPy);
   fChain->SetBranchAddress("MuonPz", MuonPz, &b_MuonPz);
   fChain->SetBranchAddress("TrigPrescale", TrigPrescale, &b_trigPrescale);
   fChain->SetBranchAddress("L1NumberOfMuons", &L1NumberOfMuons, &b_L1NumberOfMuons);
   fChain->SetBranchAddress("L1MuonEta", L1MuonEta, &b_L1MuonEta);
   fChain->SetBranchAddress("L1MuonPhi", L1MuonPhi, &b_L1MuonPhi);
   fChain->SetBranchAddress("L1MuonPt", L1MuonPt, &b_L1MuonPt);
   fChain->SetBranchAddress("L1MuonPx", L1MuonPx, &b_L1MuonPx);
   fChain->SetBranchAddress("L1MuonPy", L1MuonPy, &b_L1MuonPy);
   fChain->SetBranchAddress("L1MuonPz", L1MuonPz, &b_L1MuonPz);
   fChain->SetBranchAddress("L1MuonEnergy", L1MuonEnergy, &b_L1MuonEnergy);
   fChain->SetBranchAddress("L1MuonCharge", L1MuonCharge, &b_L1MuonCharge);
   fChain->SetBranchAddress("HltNumberOfMuons", &HltNumberOfMuons, &b_HltNumberOfMuons);
   fChain->SetBranchAddress("HltMuonEta", HltMuonEta, &b_HltMuonEta);
   fChain->SetBranchAddress("HltMuonPhi", HltMuonPhi, &b_HltMuonPhi);
   fChain->SetBranchAddress("HltMuonPt", HltMuonPt, &b_HltMuonPt);
   fChain->SetBranchAddress("HltMuonPx", HltMuonPx, &b_HltMuonPx);
   fChain->SetBranchAddress("HltMuonPy", HltMuonPy, &b_HltMuonPy);
   fChain->SetBranchAddress("HltMuonPz", HltMuonPz, &b_HltMuonPz);
   fChain->SetBranchAddress("HltMuonEnergy", HltMuonEnergy, &b_HltMuonEnergy);
   fChain->SetBranchAddress("HltMuonCharge", HltMuonCharge, &b_HltMuonCharge);
   fChain->SetBranchAddress("IsMatchedGenMuonExist", IsMatchedGenMuonExist, &b_IsMatchedGenMuonExist);
   fChain->SetBranchAddress("ChargeGenMuon", ChargeGenMuon, &b_ChargeGenMuon);
   fChain->SetBranchAddress("GenMuonPx", GenMuonPx, &b_GenMuonPx);
   fChain->SetBranchAddress("GenMuonPy", GenMuonPy, &b_GenMuonPy);
   fChain->SetBranchAddress("GenMuonPz", GenMuonPz, &b_GenMuonPz);
   fChain->SetBranchAddress("GenMuonMotherPDG", GenMuonMotherPDG, &b_GenMuonMotherPDG);
   fChain->SetBranchAddress("GenNumberOfJets", &GenNumberOfJets, &b_GenNumberOfJets);
   fChain->SetBranchAddress("GenJetEta", GenJetEta, &b_GenJetEta);
   fChain->SetBranchAddress("GenJetPhi", GenJetPhi, &b_GenJetPhi);
   fChain->SetBranchAddress("GenJetPt", GenJetPt, &b_GenJetPt);
   fChain->SetBranchAddress("GenJetPx", GenJetPx, &b_GenJetPx);
   fChain->SetBranchAddress("GenJetPy", GenJetPy, &b_GenJetPy);
   fChain->SetBranchAddress("GenJetPz", GenJetPz, &b_GenJetPz);
   fChain->SetBranchAddress("GenJetEnergy", GenJetEnergy, &b_GenJetEnergy);
   fChain->SetBranchAddress("PdgHiggsBoson", &PdgHiggsBoson, &b_PdgHiggsBoson);
   fChain->SetBranchAddress("EtaHiggsBoson", &EtaHiggsBoson, &b_EtaHiggsBoson);
   fChain->SetBranchAddress("PhiHiggsBoson", &PhiHiggsBoson, &b_PhiHiggsBoson);
   fChain->SetBranchAddress("PtHiggsBoson", &PtHiggsBoson, &b_PtHiggsBoson);
   fChain->SetBranchAddress("PxHiggsBoson", &PxHiggsBoson, &b_PxHiggsBoson);
   fChain->SetBranchAddress("PyHiggsBoson", &PyHiggsBoson, &b_PyHiggsBoson);
   fChain->SetBranchAddress("PzHiggsBoson", &PzHiggsBoson, &b_PzHiggsBoson);
   fChain->SetBranchAddress("MassHiggsBoson", &MassHiggsBoson, &b_MassHiggsBoson);
   fChain->SetBranchAddress("EnergyHiggsBoson", &EnergyHiggsBoson, &b_EnergyHiggsBoson);
   fChain->SetBranchAddress("NumberOfPartons", &NumberOfPartons, &b_NumberOfPartons);
   fChain->SetBranchAddress("PartonPdg", &PartonPdg, &b_PartonPdg);
   fChain->SetBranchAddress("PartonPx", &PartonPx, &b_PartonPx);
   fChain->SetBranchAddress("PartonPy", &PartonPy, &b_PartonPy);
   fChain->SetBranchAddress("PartonPz", &PartonPz, &b_PartonPz);
   fChain->SetBranchAddress("PartonEnergy", &PartonEnergy, &b_PartonEnergy);
   fChain->SetBranchAddress("PartonEta", &PartonEta, &b_PartonEta);
   fChain->SetBranchAddress("PartonPhi", &PartonPhi, &b_PartonPhi);
   fChain->SetBranchAddress("PartonPt", &PartonPt, &b_PartonPt);
   Notify();
}

Bool_t HBBTo4B_MC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void HBBTo4B_MC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t HBBTo4B_MC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
