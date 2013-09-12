//  Sample macro using the event shape analyzers in libEventShape.so.
//  This macro reads the events generated by a test program of Event.
//
//                                          May/13/2000  M. Iwasaki


#include <stdlib.h>

int sample_thrust(int nEvent=10)
{
  gROOT->Reset();

//   Load libJetFind library 
  gSystem->Load("libEventShape.so");
//   Load libEvent library in $ROOTSYS/test
  gSystem->Load("libEvent.so");

//   Use the events generated in $ROOTSYS/test
  TFile f("Event.root");
//   Read Tree named "T" in memory. Tree pointer is assigned the same name
  TTree *T = (TTree*)f.Get("T");

//   Set up the EventShapeFinder
  EventShape* eshape = new EventShape();

//   Start main  analysis
  Event *event = new Event();   

  TBranch *branch  = T->GetBranch("event");
  branch->SetAddress(&event);
  Int_t nevent = T->GetEntries();
  if (nevent < nEvent) nEvent = nevent;
  
// Event loop
  for (int iEvent = 0; iEvent < nEvent; iEvent++) {
      
    T->GetEvent(iEvent);      

//  Number of the tracks
    Int_t nTracks = event->GetNtrack(); 
    
//  Event cut # track >= 2
    if (nTracks < 2) continue;
    
    Float_t Pvec3[3];

    TObjArray* veclist = new TObjArray();
//  Track loop    
    TClonesArray *tracks = event->GetTracks();
    for (Int_t ntrk=0; ntrk<nTracks; ntrk++) {
      Track* trk = (Track*)tracks->UncheckedAt(ntrk);
	
      Pvec3[0] = trk->GetPx();
      Pvec3[1] = trk->GetPy();
      Pvec3[2] = trk->GetPz();
      
      TVector3* vec3 = new TVector3(Pvec3);

      veclist->Add(vec3); // Make particle 3-vector list
    }

//  Do Thrust finding
    eshape->setPartList(veclist); // Input particle 3-vector list
    
    double thrust = eshape->thrust()->X();
    double major  = eshape->thrust()->Y();
    double minor  = eshape->thrust()->Z();
    
    printf (" #Tracks %i  Thrust %f  Major %f  Minor %f\n", 
	    nTracks, thrust, major, minor);


    veclist->Delete(); delete veclist;    
    event->Clear();
  }

  f.Close();
}