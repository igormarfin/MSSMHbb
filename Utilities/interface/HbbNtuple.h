#ifndef HbbNtuple_h
#define HbbNtuple_h

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
float nPUTruth;
int nPUI;
int nPUIM1; // bx = -1
int nPUIP1; // bx = +1


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
float areaJet[100];
float ptdJet[100];

int numberOfConstituentsInJet[100];
int numberOfChargedConstituentsInJet[100];
float neutralHadronEnergyFraction[100];
float photonEnergyFraction[100];
float electronEnergyFraction[100];
float chargedHadronEnergyFraction[100];

bool jetIDLoose[100];
bool jetIDMedium[100];
bool jetIDTight[100];

float puJetMVA[100]; 
bool puJetIDLoose[100];
bool puJetIDMedium[100];
bool puJetIDTight[100];

unsigned char BContentJet3[100]; 
unsigned char CContentJet3[100]; 
unsigned char BContentJet4[100]; 
unsigned char CContentJet4[100]; 
unsigned char BContentJet5[100]; 
unsigned char CContentJet5[100]; 
unsigned char BContentJet6[100]; 
unsigned char CContentJet6[100]; 
unsigned char BContentJet7[100]; 
unsigned char CContentJet7[100]; 
unsigned char BContentJet10[100]; 
unsigned char CContentJet10[100]; 

int nSvTracksJet[100];
float svMassJet[100];
float svFDsigJet[100];

// jet - parton match 
bool isJetMatchedPartonExist[100];
int statusJetMatchedParton[100];
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
unsigned int isJetWithHltBtagBitPattern[100];
unsigned int isJetWithL25JetBitPattern[100];

float deltaRWithHltBtag[100]; // deltaR for btag match
float deltaPtWithHltBtag[100]; // pt (jet) - pt(btag match)

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

// Missing transverse momentum;
float met;
float missingPx;
float missingPy;

// trigger prescale factors
float trigPrescale[25];

  
// Generated Jets
int numberOfGenJets;
float pxGenJet[100];
float pyGenJet[100];
float pzGenJet[100];
float energyGenJet[100];
float etaGenJet[100];
float phiGenJet[100];
float ptGenJet[100];
  
// jet - L2 Jet
int l2NumberOfJets;
float l2EtaJet[100];
float l2PhiJet[100];
float l2PtJet[100];
float l2PxJet[100];
float l2PyJet[100];
float l2PzJet[100];
float l2EnergyJet[100];

// jet - L2 Jet
int l25NumberOfJets;
float l25EtaJet[100];
float l25PhiJet[100];
float l25PtJet[100];
float l25PxJet[100];
float l25PyJet[100];
float l25PzJet[100];
float l25EnergyJet[100];



// jet - HLT BTag Jet
int l3NumberOfBJets;
float l3EtaBJet[100];
float l3PhiBJet[100];
float l3PtBJet[100];
float l3PxBJet[100];
float l3PyBJet[100];
float l3PzBJet[100];
float l3EnergyBJet[100];

float BminDeltaR[100]; //min. delta R between reconstructed jet and b quark
float CminDeltaR[100]; //min. delta R between reconstructed jet and c quark





int numberOfTracksJet[100];
float trackPxJet[100][3];
float trackPyJet[100][3];
float trackPzJet[100][3];



int numberOfPartonsFull;
int stablePartonStatusFull;
float ptPartonFull[500];
float etaPartonFull[500];
float phiPartonFull[500];
float energyPartonFull[500];
int   pdgIdPartonFull[500];
int   statusPartonFull[500];


#endif // #ifdef 
