#ifndef EVENTSHAPE
#define EVENTSHAPE

#include "TMath.h"
#include "TMatrix.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TClass.h"

class EventShape : public TObject{
	

public:
	EventShape();
	~EventShape();

	void setPartList(TObjArray* e);
	// Input the particle 3(4)-vector list
	// e: 3-vector  TVector3       ..(px,py,pz) or
	//    4-vector  TLorentzVector ..(px,py,pz,E) 
	// Even input the TLorentzVector, we don't use Energy 
	
	void     setThMomPower(Double_t tp);
	Double_t getThMomPower();
	void     setFast(Int_t nf);
	Int_t    getFast();

	TVector3 thrustAxis();
	TVector3 majorAxis();
	TVector3 minorAxis();
	
	TVector3 thrust(); 
	// thrust :: Corresponding thrust, major, and minor value.

	Double_t oblateness();

private:	
	Double_t ulAngle(Double_t x, Double_t y);
	Double_t sign(Double_t a, Double_t b);
	void     ludbrb(TMatrix *mom, 
			Double_t the, 
			Double_t phi, 
			Double_t bx, 
			Double_t by,
			Double_t bz);
	
	Int_t iPow(Int_t man, Int_t exp);
	
	Double_t m_dSphMomPower; 
	// PARU(41): Power of momentum dependence in sphericity finder.

	Double_t m_dDeltaThPower;
	// PARU(42): Power of momentum dependence in thrust finder.	

	Int_t m_iFast; 
	// MSTU(44): # of initial fastest particles choosen to start search.

	Double_t m_dConv;
	// PARU(48): Convergence criteria for axis maximization.

	Int_t m_iGood;
	// MSTU(45): # different starting configurations that must
	// converge before axis is accepted as correct.	

	TMatrix m_dAxes;
	// data: results
	// m_dAxes[1] is the Thrust axis.
	// m_dAxes[2] is the Major axis.
	// m_dAxes[3] is the Minor axis.

	TVector3 m_ThrustAxis;
	TVector3 m_MajorAxis;
	TVector3 m_MinorAxis;
	TVector3 m_Thrust;

	TRandom m_random;

	Double_t m_dThrust[4];
	Double_t m_dOblateness;
	
	static Int_t m_maxpart;
public:
	ClassDef(EventShape,0)
};

#endif


