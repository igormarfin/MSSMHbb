#ifndef JETFINDER_H
#define JETFINDER_H

// Virtual JetFinder base class
//
// V0.0 Mar 01/99 : R. Shanks  Derived from Java routines written by G.Bower. 
// V1.0 May 15/00 : M.Iwasaki  Make necessary modifications, and change classes
//                             Merge JadeEJetFinder, JadeJetFinder, 
//                             DurhamJetfinder, and JetFinder
//                             into one JetFinder.
// V1.1 Aug 15 00 : M.Iwasaki  Change calcinvmass. Remove DurhamJetFinder, 
//                             JadeJetFinder, and JadeEJetFinder classes.
//                             Add setDURHAM, setJADE, or setJADEE. 

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TClass.h"

//_____________________________________________________________________
//  ----------------------
//   JetFinder Base Class
//  ----------------------
//
//class JetFinder : public TObject{  
class JetFinder {
public:
  JetFinder(Double_t ycut = 0.); 
  virtual ~JetFinder();             
  
  void setPartList(TObjArray* e);    // Input the particle 4(3)-vector list
                                     // e: 4-vector  TLorentzVector ..(px,py,pz,E) or
                                     //    3-vector  TVector3       ..(px,py,pz) 
                                     // If you input TVector3, the energy of particle
                                     // will be E = sqrt(px**2 + py**2 + pz**2) 
  void setYCut(Double_t ycut);       // Set the YCut value

  void doFindJets();                 // Clustering the particles into Jets  
    
  Int_t njets(){ return m_njets; };     // The number of jets found.  
  TLorentzVector* jet4vec(Int_t index); // Return the 4 vector of a jet (particle sum).
                                        // index: The index of the jet of interest  
  Int_t nParticlesPerJet(Int_t index);  // Number of particles in a particular jet
                                        // index: The index of the jet of interest  
  TArrayI getPartIndex(){ return m_ipart_jet_assoc; };	// Return the particle index array.
                                        // m_ipart_jet_assoc[i] = j means 
                                        // particle i is placed in jet j.  
  Int_t fewestParts(){ return m_ifewest_parts; }; // minimum number of particles to make a jet.  
  Double_t getYCut() { return m_dycut; }; // Obtain the current ycut
  
  void setDURHAM(); // Select DURHAM algorithm
  void setJADE();   //        JADE   algorithm
  void setJADEE();  //        JADE E algorithm

  Double_t calcinvmass(const TLorentzVector &jet1, 
		       const TLorentzVector &jet2);

private:
  Int_t m_njets;     // Number of jets found  
  TObjArray m_jet;   // m_jet[i] is the 4 vector sum of all the particles in jet i.  
  TArrayI m_ipart_jet_assoc; // m_ipart_jet_assoc[i] = j means particle i was placed in jet j.  
  TArrayI m_inparts_per_jet; // m_inparts_per_jet[i] = j means jet i has j particles in it.  
  Int_t m_ifewest_parts; // m_ifewest_parts is number of particles in the jet 
                         // with the fewest particles.  
  Double_t m_evis;
  Double_t m_dycut;  
  Int_t m_algorithm; // Algorithm used in Jet clustering

  const static Int_t UNASSOC;
  const static Int_t DURHAM;
  const static Int_t JADE;  
  const static Int_t JADEE;

protected:  
  TObjArray* m_4vec;

  ClassDef(JetFinder,1) // Jetfinder base class	    
};       


#endif
    

