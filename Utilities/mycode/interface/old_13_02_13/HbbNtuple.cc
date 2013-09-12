
hbbtree->SetBranchAddress("Run",&run);
hbbtree->SetBranchAddress("Event",&event);
hbbtree->SetBranchAddress("Lumi",&lumi);

  // pileup information
hbbtree->SetBranchAddress("NumberOfPU",&nPU);
hbbtree->SetBranchAddress("NumberOfPUInTime",&nPUInTime);

//#if defined(CMSSW535)
hbbtree->SetBranchAddress("nPUI",&nPUI);
hbbtree->SetBranchAddress("nPUIM1",&nPUIM1);
hbbtree->SetBranchAddress("nPUIP1",&nPUIP1);
hbbtree->SetBranchAddress("nPUTruth",&nPUTruth);

hbbtree->SetBranchAddress("PUJetMVA",puJetMVA);
hbbtree->SetBranchAddress("PUJetIDLoose",puJetIDLoose);
hbbtree->SetBranchAddress("PUJetIDMedium",puJetIDMedium);
hbbtree->SetBranchAddress("PUJetIDTight",puJetIDTight);
//#endif


  // bit pattern for triggers defining PD
hbbtree->SetBranchAddress("trgAccept",&trgAccept);

  // primary vertices
hbbtree->SetBranchAddress("NumberOfPV",&nPV);
hbbtree->SetBranchAddress("ProbPV",probPV);
hbbtree->SetBranchAddress("NdofPV",ndofPV);
hbbtree->SetBranchAddress("NTrkPV",nTrkPV);
hbbtree->SetBranchAddress("Chi2PV",chi2PV);
hbbtree->SetBranchAddress("XPV",xPV);
hbbtree->SetBranchAddress("YPV",yPV);
hbbtree->SetBranchAddress("ZPV",zPV);
hbbtree->SetBranchAddress("SumTrkPtPV",sumPtPV);

// reconstructed jets
hbbtree->SetBranchAddress("NumberOfJets",&numberOfJets);
hbbtree->SetBranchAddress("JetEta",etaJet);
hbbtree->SetBranchAddress("JetPhi",phiJet);
hbbtree->SetBranchAddress("JetPt",ptJet);
hbbtree->SetBranchAddress("JetPx",pxJet);
hbbtree->SetBranchAddress("JetPy",pyJet);
hbbtree->SetBranchAddress("JetPz",pzJet);
hbbtree->SetBranchAddress("JetEnergy",energyJet);
hbbtree->SetBranchAddress("nJetConstituents",numberOfConstituentsInJet);
hbbtree->SetBranchAddress("nJetChargedConstituents",numberOfChargedConstituentsInJet);
hbbtree->SetBranchAddress("neutralHadronEnergyFraction",neutralHadronEnergyFraction);
hbbtree->SetBranchAddress("photonEnergyFraction",photonEnergyFraction);
hbbtree->SetBranchAddress("chargedHadronEnergyFraction",chargedHadronEnergyFraction);
hbbtree->SetBranchAddress("electronEnergyFraction",electronEnergyFraction);

hbbtree->SetBranchAddress("JetIDLoose",jetIDLoose);
hbbtree->SetBranchAddress("JetIDMedium",jetIDMedium);
hbbtree->SetBranchAddress("JetIDTight",jetIDTight);


hbbtree->SetBranchAddress("JetSvNumTracks",nSvTracksJet);
hbbtree->SetBranchAddress("JetSvMass",svMassJet);
hbbtree->SetBranchAddress("JetSvFDSig",svFDsigJet);

if (_doMC) {
  hbbtree->SetBranchAddress("IsJetMatchedPartonExist",isJetMatchedPartonExist);
  hbbtree->SetBranchAddress("IsJetFromHiggs",isJetMatchedPartonFromHiggsBosonDecay);
  hbbtree->SetBranchAddress("MatchedPartonFlavor",flavorJetMatchedParton);
  hbbtree->SetBranchAddress("PartonFlavorJet",partonFlavorJet);
  hbbtree->SetBranchAddress("HflContentJet", hflContentJet);
 }
// B-Tag Discriminants -->
hbbtree->SetBranchAddress("JetBProbBJetTag",jetBProbBJetTag);
hbbtree->SetBranchAddress("JetProbBJetTag",jetProbBJetTag);
hbbtree->SetBranchAddress("TCHPBJetTag",tcHPBJetTag);
hbbtree->SetBranchAddress("TCHEBJetTag",tcHEBJetTag);
hbbtree->SetBranchAddress("nTCHPBJetTag",ntcHPBJetTag);
hbbtree->SetBranchAddress("nTCHEBJetTag",ntcHEBJetTag);
hbbtree->SetBranchAddress("SVHPBJetTag",svHPBJetTag);
hbbtree->SetBranchAddress("SVHEBJetTag",svHEBJetTag);
hbbtree->SetBranchAddress("CombSVBJetTag",combSVBJetTag);
hbbtree->SetBranchAddress("CombSVMVABJetTag",combSVMVABJetTag);
hbbtree->SetBranchAddress("TCMinus2ndBJetTag",tcMinus2ndBJetTag);
hbbtree->SetBranchAddress("TCMinus3rdBJetTag",tcMinus3rdBJetTag);
hbbtree->SetBranchAddress("TCNumberOfSelectedTracks",tcNumberOfSelectedTracks);

// HLT btag match information
//hbbtree->SetBranchAddress("IsJetWithHltBtag",isJetWithHltBtag);
hbbtree->SetBranchAddress("IsJetWithHltBtagBitPattern",isJetWithHltBtagBitPattern);
hbbtree->SetBranchAddress("DeltaRWithHltBtag",deltaRWithHltBtag);
hbbtree->SetBranchAddress("DeltaPtWithHltBtag",deltaPtWithHltBtag);

// reconstructed muons --->

hbbtree->SetBranchAddress("NumberOfMuons",&numberOfMuons);
hbbtree->SetBranchAddress("MuonCharge",chargeMuon);
hbbtree->SetBranchAddress("IsTrackerMuon",isTrackerMuon);
hbbtree->SetBranchAddress("IsGlobalMuon",isGlobalMuon);
hbbtree->SetBranchAddress("PFChHadronIsoMuon",pfChargedHadronIsoMuon);
hbbtree->SetBranchAddress("PFNeHadronIsoMuon",pfNeutralHadronIsoMuon);
hbbtree->SetBranchAddress("PFGammaIsoMuon",pfGammaIsoMuon);
hbbtree->SetBranchAddress("IsMuonInJet",isMuonInJet);
hbbtree->SetBranchAddress("MuonAssociatedJetIndex",associatedJetIndex);
hbbtree->SetBranchAddress("MuonEta",etaMuon);
hbbtree->SetBranchAddress("MuonPhi",phiMuon);
hbbtree->SetBranchAddress("MuonPt",ptMuon);
hbbtree->SetBranchAddress("MuonPx",pxMuon);
hbbtree->SetBranchAddress("MuonPy",pyMuon);
hbbtree->SetBranchAddress("MuonPz",pzMuon);
hbbtree->SetBranchAddress("TrigPrescale",trigPrescale);
if (_doMC) {
  hbbtree->SetBranchAddress("IsMatchedGenMuonExist",isMatchedGeneratorMuonExist);
  hbbtree->SetBranchAddress("ChargeGenMuon",chargeGenMuon);
  hbbtree->SetBranchAddress("GenMuonPx",pxGenMuon);
  hbbtree->SetBranchAddress("GenMuonPy",pyGenMuon);
  hbbtree->SetBranchAddress("GenMuonPz",pzGenMuon);
  hbbtree->SetBranchAddress("GenMuonMotherPDG",pdgMotherGenMuon);

  hbbtree->SetBranchAddress("PdgHiggsBoson",&pdgHiggsBoson);
  hbbtree->SetBranchAddress("EtaHiggsBoson",&etaHiggsBoson);
  hbbtree->SetBranchAddress("PhiHiggsBoson",&phiHiggsBoson);
  hbbtree->SetBranchAddress("PtHiggsBoson",&ptHiggsBoson);
  hbbtree->SetBranchAddress("PxHiggsBoson",&pxHiggsBoson);
  hbbtree->SetBranchAddress("PyHiggsBoson",&pyHiggsBoson);
  hbbtree->SetBranchAddress("PzHiggsBoson",&pzHiggsBoson);
  hbbtree->SetBranchAddress("MassHiggsBoson",&massHiggsBoson);
  hbbtree->SetBranchAddress("EnergyHiggsBoson",&energyHiggsBoson);

  hbbtree->SetBranchAddress("NumberOfPartons",&numberOfPartons);
  hbbtree->SetBranchAddress("PartonPdg",pdgParton);
  hbbtree->SetBranchAddress("PartonPx",pxParton);
  hbbtree->SetBranchAddress("PartonPy",pyParton);
  hbbtree->SetBranchAddress("PartonPz",pzParton);
  hbbtree->SetBranchAddress("PartonEnergy",energyParton);
  hbbtree->SetBranchAddress("PartonEta",etaParton);
  hbbtree->SetBranchAddress("PartonPhi",phiParton);
  hbbtree->SetBranchAddress("PartonPt",ptParton);

  hbbtree->SetBranchAddress("GenNumberOfJets",&numberOfGenJets);
  hbbtree->SetBranchAddress("GenJetPx",pxGenJet);
  hbbtree->SetBranchAddress("GenJetPy",pyGenJet);
  hbbtree->SetBranchAddress("GenJetPz",pzGenJet);
  hbbtree->SetBranchAddress("GenJetEnergy",energyGenJet);
  hbbtree->SetBranchAddress("GenJetEta",etaGenJet);
  hbbtree->SetBranchAddress("GenJetPhi",phiGenJet);
  hbbtree->SetBranchAddress("GenJetPt",ptGenJet);
 }

// missing transverse energy and momentum
hbbtree->SetBranchAddress("met",&met);
hbbtree->SetBranchAddress("missingPx",&missingPx);
hbbtree->SetBranchAddress("missingPy",&missingPy);




hbbtree->SetBranchAddress("L2NumberOfJets",&l2NumberOfJets);
hbbtree->SetBranchAddress("L2JetEta",l2EtaJet);
hbbtree->SetBranchAddress("L2JetPhi",l2PhiJet);
hbbtree->SetBranchAddress("L2JetPt",l2PtJet);
hbbtree->SetBranchAddress("L2JetPx",l2PxJet);
hbbtree->SetBranchAddress("L2JetPy",l2PyJet);
hbbtree->SetBranchAddress("L2JetPz",l2PzJet);
hbbtree->SetBranchAddress("L2JetEnergy",l2EnergyJet);
