// Jet finder base class routines.
//
// V0.0 Mar 01 99 : R. Shanks  Derived from Java routines written by G.Bower. 
// V0.1 Jul 02 99 : M.Iwasaki  Make some Modification on masscut and ymass 
//                             calculation (in JadeEJetfonder.cxx 
//                             or DurhamJetFinder.cxx) on doFindJets().
//                             Adding Jade cluster (yclus).
// V0.2 Jul 07 99 : M.Iwasaki  Change Mod to Mag function in TVector3
//                             so as to use root 2.22.x
// V0.3 Jul 21 99 : M.Iwasaki  Make temporary 4-vectors copied from input 
//                             particle 4-vectors in doFindJets() 
//                             to protect the initial 4-vectors.
// V0.4 Sep 23 99 : M.Iwasaki  Fix setEvent memory leak.
//
// V1.0 May 15 00 : M.Iwasaki  Make necessary modifications, and change classes
//                             Merge JadeEJetFinder, JadeJetFinder, 
//                             DurhamJetfinder, and JetFinder into one JetFinder.
// V1.1 Aug 15 00 : M.Iwasaki  Change calcinvmass. Remove DurhamJetFinder, 
//                             JadeJetFinder, and JadeEJetFinder classes.
//                             Add setDURHAM, setJADE, or setJADEE. 
//                             Now we can cange algorithms by one of them.
//

#include "JetFinder.h"

//_____________________________________________________________________
//  ------------------
//   JetFinder Class
//  ------------------
//
ClassImp(JetFinder)

const Int_t JetFinder::UNASSOC = -999;

const Int_t JetFinder::DURHAM = 1;
const Int_t JetFinder::JADE   = 2;
const Int_t JetFinder::JADEE  = 3;

JetFinder::JetFinder(Double_t ycut): 
    m_evis(0.),m_dycut(ycut),
    m_algorithm(JetFinder::DURHAM){ // Default constructor  uses Durham
    m_4vec = new TObjArray();
}
//______________________________________________________

JetFinder::~JetFinder() {
    m_jet.Delete();
    m_4vec->Delete(); delete m_4vec;
}
//______________________________________________________

TLorentzVector* JetFinder::jet4vec(Int_t index) {
  return (TLorentzVector*)m_jet.At(index);
}
//_______________________________________________________

Int_t JetFinder::nParticlesPerJet(Int_t index) {
  return m_inparts_per_jet[index];
}
//_______________________________________________________

void JetFinder::setYCut(Double_t ycut) {
  m_dycut = ycut;
}
//_______________________________________________________

// Input the particle 4(3)-vector list
// e: 4-vector  TLorentzVector ..(px,py,pz,E) or
//    3-vector  TVector3       ..(px,py,pz) 
// If you input TVector3, the energy of particle
// will be E = sqrt(px**2 + py**2 + pz**2)
void JetFinder::setPartList(TObjArray* e) {

  m_evis = 0;
  m_4vec->Delete();

  Int_t ne = e->GetEntries();
  for(Int_t i=0;i<ne;i++) {
    TObject* o = e->At(i);
    TString nam(o->IsA()->GetName());
    if (nam.Contains("TLorentzVector")) {
      TVector3 vec3(((TLorentzVector*) o)->X(),
		    ((TLorentzVector*) o)->Y(),
		    ((TLorentzVector*) o)->Z());   	
      TLorentzVector* in = 
	new TLorentzVector(vec3,((TLorentzVector *) o)->T());
      m_evis += in->T();
      m_4vec->Add(in);
    }    
    else if (nam.Contains("TVector3")) {
      TVector3 vec3(((TVector3 *) o)->X(),
		    ((TVector3 *) o)->Y(),
		    ((TVector3 *) o)->Z());
      TLorentzVector* in = 
	new TLorentzVector(vec3,((TVector3 *) o)->Mag());
      m_evis += in->T();
      m_4vec->Add(in);
    }
    else {
      printf("JetFinder::setEvent input is not a TVector3 or a TLorentzVector\n");
    }
  }
}
//_______________________________________________________

void JetFinder::doFindJets(){

  Int_t np = m_4vec->GetEntries(); 	
  if (np<2) return;

  TObjArray* part = new TObjArray();
  for (Int_t Ipart=0; Ipart <np; Ipart++) {
    TVector3 vec3(((TLorentzVector*)m_4vec->At(Ipart))->X(),
		  ((TLorentzVector*)m_4vec->At(Ipart))->Y(),
		  ((TLorentzVector*)m_4vec->At(Ipart))->Z());
    TLorentzVector* vec4 = 
      new TLorentzVector(vec3,
			 ((TLorentzVector *)(m_4vec->At(Ipart)))->T());
    part->Add(vec4);
  }

  m_ipart_jet_assoc.Reset();
  m_ipart_jet_assoc.Set(np);
  for (Int_t m=0; m<np; m++) m_ipart_jet_assoc[m] = UNASSOC;

  m_njets = 0;

  //
  // create invariant mass pair array.
  //
  TMatrix ymass = TMatrix(np,np);
  
  for (Int_t i1 = 0; i1 < np - 1; i1++ ) {
    for (Int_t i2 = i1 + 1 ; i2 < np ; i2++ ) {
      TLorentzVector &jeti1 = *(TLorentzVector*)part->At(i1);
      TLorentzVector &jeti2 = *(TLorentzVector*)part->At(i2);

      ymass(i1,i2) = calcinvmass(jeti1, jeti2);

    }
  }
  
  Double_t masscut = m_dycut * m_evis * m_evis ;

  for (;;)
    {
      Int_t im = -1;
      Int_t jm = -1;
      Double_t minmass = 100000000.;
      if ( minmass <= masscut ) minmass = masscut * 10.;
      //
      // find least invariant mass pair.
      //
      for(Int_t i = 0 ; i < np - 1 ; i++ ) {
	for(Int_t j = i + 1 ; j < np ; j++ ) {
	  if (m_ipart_jet_assoc[i] != JetFinder::UNASSOC) continue;
	  if (m_ipart_jet_assoc[j] != JetFinder::UNASSOC) continue;
	  if (ymass(i,j) > minmass) continue;
	  
	  minmass = ymass(i,j);
	  im = i;  jm = j;
	}
      }

      if (minmass > masscut) break;
      //
      // combine particles im and jm.
      //
      *(TLorentzVector*)part->At(im) += *(TLorentzVector*)part->At(jm);

      for(Int_t ipart = 0; ipart < np; ipart++ ){
	if(m_ipart_jet_assoc[ipart] == jm) m_ipart_jet_assoc[ipart] = im;
      }
      //
      // Recalculate invariant masses for newly combined particle
      //
      m_ipart_jet_assoc[jm] = im;
      for (Int_t ipar = 0; ipar < np ; ipar++) {
	if (ipar == im) continue;
	if (m_ipart_jet_assoc[ipar] != UNASSOC ) continue;
	
	Int_t imin = TMath::Min(ipar,im);
	Int_t imax = TMath::Max(ipar,im);
	
	TLorentzVector &jetimin = *(TLorentzVector*)part->At(imin);
	TLorentzVector &jetimax = *(TLorentzVector*)part->At(imax);

	ymass(imin,imax) = calcinvmass(jetimin, jetimax);

      }

    }
  //
  // finish up by filling jet array.
  //
  
  for(Int_t ip = 0 ; ip < np ; ip++) {
    if (m_ipart_jet_assoc[ip] == UNASSOC) m_njets++;			
  }
  
  m_jet.Delete();
  m_inparts_per_jet.Reset();
  m_inparts_per_jet.Set(m_njets);  

  Int_t nj = 0;
  Int_t npart;
  m_ifewest_parts = 100000; // Starting min value
  for(Int_t i = 0 ; i < np ; i++ ){
    if (m_ipart_jet_assoc[i] != UNASSOC) continue;

    TVector3 vec3(((TLorentzVector*)part->At(i))->X(),
		  ((TLorentzVector*)part->At(i))->Y(),
		  ((TLorentzVector*)part->At(i))->Z());   	
    TLorentzVector* jet = 
      new TLorentzVector(vec3,((TLorentzVector*)part->At(i))->T());
    m_jet.Add(jet);

    npart = 1;
    for (Int_t j = 0 ; j < np ; j++) {
      if(m_ipart_jet_assoc[j] == i) {
	m_ipart_jet_assoc[j] = nj;
	npart++;
      }
    }
    m_ipart_jet_assoc[i] = nj;
    m_inparts_per_jet[nj] = npart;
    if( npart < m_ifewest_parts) m_ifewest_parts = npart;
    nj++;
  }  
  part->Delete(); delete part;  
}

void JetFinder::setDURHAM(){
  m_algorithm = JetFinder::DURHAM; 
}
void JetFinder::setJADE(){
  m_algorithm = JetFinder::JADE;
}
void JetFinder::setJADEE(){
  m_algorithm = JetFinder::JADEE;
}

Double_t JetFinder::calcinvmass(const TLorentzVector& jet1,
				   const TLorentzVector& jet2){
  TVector3 P_jet1 = jet1.Vect();
  TVector3 P_jet2 = jet2.Vect();
  Double_t costh = (P_jet1 * P_jet2)/P_jet1.Mag()/P_jet2.Mag();

  if     (m_algorithm == JetFinder::DURHAM) {  // DURHAM
    Double_t minT = TMath::Min(jet1.E(),jet2.E());
    return 2. * minT*minT * (1.-costh);
  }  
  else if (m_algorithm == JetFinder::JADE)   {  // JADE    
    return 2. * (jet1.E())*(jet2.E()) * (1.-costh) ; 
  }  
  else if (m_algorithm == JetFinder::JADEE)  {  // JADE E
    return (jet1 + jet2).M2();
  }

  printf(" Strange Algorithm!!!! \n");
  return 0.;

};
