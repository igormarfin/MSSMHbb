#ifndef CalculateDeltaR_h
#define CalculateDeltaR_h

#include <TMath.h>
#include <math.h>
#include <exception>

double CalculateDeltaR (double phi1, double phi2, double eta1, double eta2)
{
  if(
     phi1 > TMath::Pi() || phi1 < -1.0*TMath::Pi()
     || phi2 > TMath::Pi() || phi2 < -1.0*TMath::Pi()
     ) {
    std::cout << "error: call of CalculateDeltaR(double phi1, double phi2, double eta1, double eta2) with phi values out of range (maybe phi and eta are swapped?): " << phi1 << " " << phi2 << std::endl;
    throw std::exception();
  }
  double delta_phi = phi1 - phi2;
  
  if(delta_phi > TMath::Pi()) delta_phi -= 2.0*TMath::Pi();
  if (delta_phi< -1.0*TMath::Pi()) delta_phi += 2.0*TMath::Pi();
  
  const double delta_eta = eta1 - eta2;
  
  const double delta_R = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
  
  return delta_R;
}

#endif
