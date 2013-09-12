#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TH1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TBranch.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <TFileCollection.h>
#include <TLorentzVector.h>
#include <TGraph2D.h>
#include <TLatex.h>
#include <TMarker.h>
#include <math.h>
#include <TVector3.h>
#include <TMath.h>
#include <TFileInfo.h>
#include "Analysis/Utilities/interface/PileUpWeight.h"
#include "Analysis/Utilities/interface/deltaR.h"
#include "Analysis/Utilities/interface/getHbbMetadata.h"
#include "Analysis/Utilities/interface/checkHLTBTagMatch.h"
#include "Analysis/Utilities/interface/HbbNtuple.h"
#include "Analysis/Utilities/interface/bTagEff.h"
//#include "Analysis/Utilities/interface/bTagEff.h_backup"
//#include "../interface/tiltt.h"
#include <map>
#include "Analysis/Utilities/interface/RelativeOnlineSFandUncertainties.h"
#define readmistag_C
#include <exception>


bool testJetIsolation(
                      const float phireferencejet, 
                      const float etareferencejet, 
                      const int njets,
                      const float *jetphi,
                      const float *jeteta,
                      const float *jetpt
                      )
{
  int match = 0;
  for(int i = 0; i < njets; i++)
    {
      if(jetpt[i] > 20.0)
        {
          float delta_phi = phireferencejet - jetphi[i];
      
          if(delta_phi > TMath::Pi()) delta_phi -= 2.0*TMath::Pi();
          if (delta_phi< -1.0*TMath::Pi()) delta_phi += 2.0*TMath::Pi();
      
          const float delta_eta = etareferencejet - jeteta[i];
      
          const float delta_R = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
          if(delta_R <= 1.0)
            match++;
          if(match >= 2)
            break;
        }
    }
  //  return match > 1 ? false : true; //true if isolated, false if not. match > 1 because the reference jet is included in the list.
  return true;
  
}

bool JetIsIsolated(
                   const float phireferencejet, 
                   const float etareferencejet, 
                   const float btcut,
                   const int njets,
                   const float *discr,
                   const float *jetphi,
                   const float *jeteta,
                   const float *jetpt
                   ) 
{
  int match = 0;
  for(int i = 0; i < njets; i++)
    {
      if(discr[i] > btcut && jetpt[i] > 20.0)
        {
          float delta_phi = phireferencejet - jetphi[i];
      
          if(delta_phi > TMath::Pi()) delta_phi -= 2.0*TMath::Pi();
          if (delta_phi< -1.0*TMath::Pi()) delta_phi += 2.0*TMath::Pi();
      
          const float delta_eta = etareferencejet - jeteta[i];
      
          const float delta_R = TMath::Sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
          if(delta_R <= 1.0)
            match++;
          if(match >= 2)
            break;
        }
    }
  //return match > 1 ? false : true; //true if isolated, false if not. match > 1 because the reference jet is included in the list.
  return true; //is always isolated
}








class JetConstituents {
public:
  JetConstituents(TFile *hout,
                  const unsigned int nbt_tmp,
                  const unsigned int nflav_tmp,
                
                  const std::string *btlabel,
                  std::string *flavlabel_tmp
                  // const std::string *modelabel
                  
                  ) : nbt(nbt_tmp), nflav_external(nflav_tmp), nflav(2)
  {

    flavlabel.push_back("slim");
    flavlabel.push_back("fat");
    flavlabel_external.assign(&flavlabel_tmp[0], &flavlabel_tmp[0]+nflav_external);

    //  const unsigned int nflav = 2;
    //     const std::string flavlabel[nflav] = {"slim","fat"};

    const unsigned int nmode = 2;
    const std::string modelabel[nmode] = {"all", "btag"};

    hout->mkdir("jet_constituents");
    hout->cd("jet_constituents");
    
    std::vector<std::vector<std::vector<TH1D*> > > initvec
      = std::vector<std::vector<std::vector<TH1D*> > > (nbt, std::vector<std::vector<TH1D*> > (nflav, std::vector<TH1D*> (nmode, NULL)) );
    h_nchargedConstituents = initvec;
    h_nConstituents = initvec;
    h_nchargedConstituents_ptWindow = initvec;
    h_nConstituents_ptWindow = initvec;
    h_prel_ptWindow = initvec;
    h_prel = initvec;
    h_p_ratio_tracks_jet = initvec;
    h_jetarea = initvec;
    h_jetarea_over_pt = initvec;
    h_ptd = initvec;
    h_csv = initvec;


    std::vector<std::vector<std::vector<TH2D*> > > initvec2d
      = std::vector<std::vector<std::vector<TH2D*> > > (nbt, std::vector<std::vector<TH2D*> > (nflav, std::vector<TH2D*> (nmode, NULL)) );
    
    h_nchargedConstituents_vs_pt =  initvec2d;
    h_nConstituents_vs_pt=  initvec2d;
    h_prel_vs_pt=  initvec2d;
    h_jetarea_vs_pt=  initvec2d;
    h_jetarea_over_pt_vs_pt=  initvec2d;
    h_ptd_vs_pt =  initvec2d;
    h_csv_vs_pt =  initvec2d;
 
    const int nptbins = 15;
    const double ptbinning[nptbins+1] = {
      2.00000e+01, 
      2.55761e+01, 
      3.27068e+01, 
      4.18256e+01, 
      5.34867e+01, 
      6.83990e+01, 
      8.74690e+01, 
      1.11856e+02, 
      1.43041e+02, 
      1.82922e+02, 
      2.33921e+02, 
      2.99140e+02, 
      3.82541e+02, 
      4.89195e+02, 
      6.25585e+02, 
      800.0
    };

    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_jetarea_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_jetarea_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 20, 0.0, 1.5);
        }
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_jetarea_over_pt_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_jetarea_over_pt_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 20, 0.0, 0.07);
        }
      }
    }
    
    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_nchargedConstituents_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nchargedConstituents_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 30, 0.0, 60.0);
        }
      }
    }
    
    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_nConstituents_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nConstituents_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 25, 0.0, 50.0);
        }
      }
    }
    


    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_csv_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_csv_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 30, 0.0, 1.0);
        }
      }
    }
    
  








    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_ptd_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_ptd_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 20, 0.0, 1.0);
        }
      }
    }
    

    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_prel_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_prel_vs_pt.at(ibt).at(iflav).at(imode) = new TH2D(scharged.c_str(),scharged.c_str(), nptbins, ptbinning, 20, 0.0, 1.0);
        }
      }
    }
    




    




    //now 1d


    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_jetarea_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_jetarea.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 100, 0.0, 1.5);
          
        
        }

        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_jetarea_over_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_jetarea_over_pt.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 100, 0.0, 0.07);
        }
      }
    }


    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_ptd_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_ptd.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 60, 0.0, 1.0);
          
        
        }
      }
    }



    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("h_csv_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_csv.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 60, 0.0, 1.0);
          
        
        }
      }
    }









    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("prel_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_prel.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 60, 0.0, 1.0);
          
        
        }
      }
    }


    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("p_ratio_tracks_jet_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_p_ratio_tracks_jet.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 400, 0.0, 4.0);
          
        
        }
      }
    }






    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("prel_ptWindow_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_prel_ptWindow.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 60, 0.0, 1.0);
        }
      }
    }




    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("nchargedConstituents_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nchargedConstituents.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 60, 0.0, 60.0);
          
          std::string sallconsti(Form("nConstituents_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nConstituents.at(ibt).at(iflav).at(imode) = new TH1D(sallconsti.c_str(),sallconsti.c_str(), 80, 0.0, 80.0);
        }
      }
    }


    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string scharged(Form("nchargedConstituents_ptWindow_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nchargedConstituents_ptWindow.at(ibt).at(iflav).at(imode) = new TH1D(scharged.c_str(),scharged.c_str(), 60, 0.0, 60.0);
          
          std::string sallconsti(Form("nConstituents_ptWindow_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nConstituents_ptWindow.at(ibt).at(iflav).at(imode) = new TH1D(sallconsti.c_str(),sallconsti.c_str(), 80, 0.0, 80.0);
        }
      }
    }
  }
  
  void Fill(const unsigned int ibt, const int iflav_external, const unsigned int imode_out, const double jetpt, const double ncharged, const double n,  const double weight, const int iJet)
  {
    return;
    //turned off!

    if(imode_out != 1 && imode_out != 0 )
      return;
    if(iflav_external == -1)
      return;
    const unsigned int imode = imode_out == 1 ? 0 : 1;

    unsigned int iflav = -1;

    //mapping jets to fat and slim categories
    if(flavlabel_external.at(iflav_external) == "c" || flavlabel_external.at(iflav_external) == "b") {
      iflav = 0;
    } else  if(flavlabel_external.at(iflav_external) == "cc" || flavlabel_external.at(iflav_external) == "bb") {
      iflav = 1;
    } else if(abs(partonFlavorJet[iJet]) == 21) {
      iflav = 1;
    } else {
      iflav = 0;
    }
    
    h_nchargedConstituents.at(ibt).at(iflav).at(imode)->Fill(ncharged, weight);
    h_nConstituents.at(ibt).at(iflav).at(imode)->Fill(n, weight);

    h_jetarea.at(ibt).at(iflav).at(imode)->Fill(areaJet[iJet], weight);
    h_jetarea_over_pt.at(ibt).at(iflav).at(imode)->Fill(areaJet[iJet]/jetpt, weight);
    h_ptd.at(ibt).at(iflav).at(imode)->Fill(ptdJet[iJet], weight);
  
    h_csv.at(ibt).at(iflav).at(imode)->Fill(combSVBJetTag[iJet], weight);
    h_csv_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,combSVBJetTag[iJet], weight);
  
    h_nchargedConstituents_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,ncharged, weight);
    h_nConstituents_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,n, weight);

    h_jetarea_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,areaJet[iJet], weight);
    h_jetarea_over_pt_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,areaJet[iJet]/jetpt, weight);
    h_ptd_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,ptdJet[iJet], weight);
  
  
    if(jetpt >= 53. && jetpt <= 120.) {
      h_nchargedConstituents_ptWindow.at(ibt).at(iflav).at(imode)->Fill(ncharged, weight);
      h_nConstituents_ptWindow.at(ibt).at(iflav).at(imode)->Fill(n, weight);
    }

    




    //calculate relative p of tracks w.r.t. the jet

    //const double pjet = TMath::Sqrt(pxJet[iJet] * pxJet[iJet] + pyJet[iJet] * pyJet[iJet] + pzJet[iJet] * pzJet[iJet]);
    TVector3 vec_jet(pxJet[iJet], pyJet[iJet], pzJet[iJet]);

    double sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0;
    for(int iTrack = 0; iTrack < numberOfTracksJet[iJet]; iTrack++) {
      sum_px += trackPxJet[iJet][iTrack];
      sum_py += trackPyJet[iJet][iTrack];
      sum_pz += trackPzJet[iJet][iTrack];
    }
    TVector3 vec_tracks(sum_px, sum_py, sum_pz);
 
    
    const double p_rel_tracks = vec_tracks.Mag() * TMath::Cos(vec_jet.Angle(vec_tracks));
    const double p_rel = p_rel_tracks / vec_jet.Mag();

    h_prel.at(ibt).at(iflav).at(imode)->Fill(p_rel, weight);
    h_p_ratio_tracks_jet.at(ibt).at(iflav).at(imode)->Fill( vec_tracks.Mag() / vec_jet.Mag(), weight);

    h_prel_vs_pt.at(ibt).at(iflav).at(imode)->Fill(jetpt,p_rel, weight);

    if(jetpt >= 53. && jetpt <= 120.) {
      h_prel_ptWindow.at(ibt).at(iflav).at(imode)->Fill(p_rel, weight);
    }





  }
 
  

  const unsigned int nbt;
  const unsigned int nflav_external;
  std::vector<std::string> flavlabel_external;
  const unsigned int nflav;
  std::vector<std::string> flavlabel;
  

  //  const unsigned int nmode;
  std::vector<std::vector<std::vector<TH1D*> > > h_nchargedConstituents;
  std::vector<std::vector<std::vector<TH1D*> > > h_nConstituents;
  
  std::vector<std::vector<std::vector<TH1D*> > > h_nchargedConstituents_ptWindow;
  std::vector<std::vector<std::vector<TH1D*> > > h_nConstituents_ptWindow;

  std::vector<std::vector<std::vector<TH1D*> > > h_prel_ptWindow;
  std::vector<std::vector<std::vector<TH1D*> > > h_prel;

  std::vector<std::vector<std::vector<TH1D*> > > h_p_ratio_tracks_jet;

  

  std::vector<std::vector<std::vector<TH1D*> > > h_jetarea;
  std::vector<std::vector<std::vector<TH1D*> > > h_jetarea_over_pt;
  std::vector<std::vector<std::vector<TH1D*> > > h_ptd;
  std::vector<std::vector<std::vector<TH1D*> > > h_csv;


  //2d histograms
  std::vector<std::vector<std::vector<TH2D*> > > h_nchargedConstituents_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_nConstituents_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_prel_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_jetarea_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_jetarea_over_pt_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_ptd_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_csv_vs_pt;

  


  
  
};


class JetFlavourFractions {
public:
  JetFlavourFractions(
                      std::string scenario,
                      TFile *hout,
                      const unsigned int nbt_tmp,
                      const unsigned int nflav_tmp,
                
                      const std::string *flavlabel,

                      const std::string *btlabel_tmp,
                      const std::string *btdiscr_tmp,
                      const float *btcut_tmp,
                      const bool *online_tmp,
                      const bool *btagMatch_tmp
                 
                  
                      ) : nbt(nbt_tmp), nflav(nflav_tmp), btlabel(btlabel_tmp),btdiscr(btdiscr_tmp),btcut(btcut_tmp),online(online_tmp),btagMatch(btagMatch_tmp)
  {
    // const std::string flavlabel[] = {"udsg","c","b","cc","bb"};
    const std::string sbtag[1] = { "CSVT" };
    bTagEffOffline = new bTagEff("/afs/naf.desy.de/user/j/jbehr/public/cms/btag/btagefficiencies_v110_53X_mediummass.root","offline",flavlabel,sbtag,5,nbt);

    hout->mkdir(Form("jet_flavour_fractions_debug_%s",scenario.c_str()));
    hout->cd(Form("jet_flavour_fractions_debug_%s",scenario.c_str()));
    
    h_debug_pt_vs_p = new TH2D("h_debug_pt_vs_p","h_debug_pt_vs_p", 20, 0.0, 5000.0, 20, 0.0, 800.0);
    

    hout->mkdir(Form("jet_flavour_fractions_%s",scenario.c_str()));
    hout->cd(Form("jet_flavour_fractions_%s",scenario.c_str()));

    std::vector<std::vector<std::vector<TH1D*> > > initvec
      = std::vector<std::vector<std::vector<TH1D*> > > (nbt, std::vector<std::vector<TH1D*> > (nflav, std::vector<TH1D*> (nrank, NULL)) );
    
    h_pt_pat0 = initvec;
    h_p_pat0 = initvec;
    h_eta_pat0 = initvec;


    h_pt_pat1 = initvec;
    h_p_pat1 = initvec;
    h_eta_pat1 = initvec;


    h_pt_pat2 = initvec;
    h_p_pat2 = initvec;
    h_eta_pat2 = initvec;


    h_pt_pat3 = initvec;
    h_p_pat3 = initvec;
    h_eta_pat3 = initvec;

    h_pt_pat4 = initvec;
    h_p_pat4 = initvec;
    h_eta_pat4 = initvec;




    h_pt_pat0_wtx = initvec;
    h_p_pat0_wtx = initvec;
    h_eta_pat0_wtx = initvec;


    h_pt_pat1_wtx = initvec;
    h_p_pat1_wtx = initvec;
    h_eta_pat1_wtx = initvec;


    h_pt_pat2_wtx = initvec;
    h_p_pat2_wtx = initvec;
    h_eta_pat2_wtx = initvec;


    h_pt_pat3_wtx = initvec;
    h_p_pat3_wtx = initvec;
    h_eta_pat3_wtx = initvec;

    h_pt_pat4_wtx = initvec;
    h_p_pat4_wtx= initvec;
    h_eta_pat4_wtx = initvec;







  
    //  const double pt_binning_1d[nch+1] = {2.00000e+01, 2.46667e+01, 2.93333e+01, 3.40000e+01, 3.86667e+01, 4.33333e+01, 4.80000e+01, 5.26667e+01, 5.73333e+01, 6.20000e+01, 6.66667e+01, 7.13333e+01, 7.60000e+01, 8.06667e+01, 8.53333e+01, 9.00000e+01, 9.46667e+01, 9.93333e+01, 1.04000e+02, 1.08667e+02, 1.13333e+02, 1.18000e+02, 1.22667e+02, 1.27333e+02, 1.32000e+02, 1.36667e+02, 1.41333e+02, 1.46000e+02, 1.50667e+02, 1.55333e+02, 1.60000e+02, 1.64667e+02, 1.69333e+02, 1.74000e+02, 1.78667e+02, 1.83333e+02, 1.88000e+02, 1.92667e+02, 1.97333e+02, 2.02000e+02, 2.06667e+02, 2.11333e+02, 2.16000e+02, 2.20667e+02, 2.25333e+02, 2.30000e+02, 2.34667e+02, 2.39333e+02, 2.44000e+02, 2.48667e+02, 2.53333e+02, 2.58000e+02, 2.62667e+02, 2.67333e+02, 2.72000e+02, 2.76667e+02, 2.81333e+02, 2.86000e+02, 2.90667e+02, 2.95333e+02, 3.00000e+2, 400.0, 550.0, 800.0};
    //  const int nch=6;
    //     const double pt_binning_1d[nch+1] = {
    //       2.00000e+01, 
    //       80.0,
    //       150.0,
    //       250.0,
    //       400.0,
    //       550.0,
    //       800.0};

    // const int nch=3;
    //     const double pt_binning_1d[nch+1] = {
    //       2.00000e+01, 
    //       150.0,
    //       300.0,
    //       800.0};


    // const int nch=4;
    //     const double pt_binning_1d[nch+1] = {
    //       2.00000e+01, 
    //       100.0,
    //       250.0,
    //       500.0,
    //       800.0};


    const int nch=9;
    const double pt_binning_1d[nch+1] = {
      2.00000e+01, 
      3.17167e+01, 
      5.02973e+01, 
      7.97632e+01, 
      1.26491e+02, 
      2.00594e+02, 
      3.18108e+02, 
      5.04467e+02, 
      800.0,
      1200.0
    };


    const int nch_p=10;
    const double p_binning_1d[nch_p+1] = {
      2.00000e+01, 
      3.35223e+01, 
      5.61872e+01, 
      9.41763e+01, 
      1.57850e+02, 
      2.64575e+02, 
      4.43458e+02, 
      7.43287e+02, 
      1.24583e+03, 
      2.08816e+03, 
      3500.0
    };


    

    //  const int nch=5;
    //     const double pt_binning_1d[nch+1] = {
    //       2.00000e+01, 
    //       150.0,
    //       300.0,
    //       400.0,
    //       600.0,
    //       800.0};


    //the following is really ugly because of copy/paste....

    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat0_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat0.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat0_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat0.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat0_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat0.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--



        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat1_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat1.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat1_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat1.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat1_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat1.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--



        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat2_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat2.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat2_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat2.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat2_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat2.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--


        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat3_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat3.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat3_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat3.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat3_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat3.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--



        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat4_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat4.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat4_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat4.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat4_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat4.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--










      }
    }








    hout->mkdir(Form("jet_flavour_fractions_wtx_%s",scenario.c_str()));
    hout->cd(Form("jet_flavour_fractions_wtx_%s",scenario.c_str()));



    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat0_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat0_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat0_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat0_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat0_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat0_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--



        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat1_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat1_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat1_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat1_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat1_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat1_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--



        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat2_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat2_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat2_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat2_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat2_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat2_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--


        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat3_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat3_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat3_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat3_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat3_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat3_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--



        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_pt_pat4_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_pt_pat4_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch,pt_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_p_pat4_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_p_pat4_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), nch_p,p_binning_1d);
        }
        for (unsigned int irank=0; irank<nrank; ++irank) {
          std::string s(Form("h_eta_pat4_wtx_%s_%s_%i_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),irank,scenario.c_str()));
          h_eta_pat4_wtx.at(ibt).at(iflav).at(irank) = new TH1D(s.c_str(),s.c_str(), 3, 0.0, 2.6);
        }

        //--


      }
    }





    //hout->mkdir(Form("jet_flavour_fractions_debug_%s",scenario.c_str()));
    hout->cd(Form("jet_flavour_fractions_debug_%s",scenario.c_str()));
    h_debug_njets = new TH1D(Form("h_debug_njets_%s",scenario.c_str()), Form("h_debug_njets_%s",scenario.c_str()),10, 0.0, 10.0);
    h_debug_njets_twobtags = new TH1D(Form("h_debug_njets_twobtags_%s",scenario.c_str()), Form("h_debug_njets_twobtags_%s",scenario.c_str()),10, 0.0, 10.0);





  }
  bool isOfflineBtagged(unsigned int ibt, unsigned int iJet)
  {
    bool tagged = false;
    bool csv = false;
    bool csvtagged = false;
    if (btdiscr[ibt] == "TCHE") {
      if (tcHEBJetTag[iJet]>btcut[ibt]) {
        tagged = true; 
      }
    } else if (btdiscr[ibt] == "TCHP") {
      if (tcHPBJetTag[iJet]>btcut[ibt]) {
        tagged = true; 
      }
    } else if (btdiscr[ibt] == "nTCHE") {
      if (ntcHEBJetTag[iJet]>btcut[ibt]) {
        tagged = true;  
      }
    }
    else if (btdiscr[ibt] == "nTCHP") {
      if (ntcHPBJetTag[iJet]>btcut[ibt]) {
        tagged = true;     
      }
    } else if (btdiscr[ibt] == "CSV") {
    
      csv = true;
      if (combSVBJetTag[iJet]>btcut[ibt]) {
        tagged = true;   
        csvtagged = true;
      }
    } else if (btdiscr[ibt] == "SSVHP") {
      if (svHPBJetTag[iJet]>btcut[ibt]) {
              
        tagged = true; 
      }
    }
    if(tagged != csvtagged) {
      std::cout << "error: bug in fraction code: tagging does not make sense" << std::endl;
      throw std::exception();
    }
    if(tagged && !csv) {
      std::cout << "error: bug in fraction code2: tagging does not make sense" << std::endl;
      throw std::exception();
    }
    return tagged;
  }
 
  void DetermineFlavFractions(const std::vector<int> &leadingJets, const double weight) {
    h_debug_njets->Fill(leadingJets.size(),weight);
  


    if(leadingJets.size() < 3)
      return;
    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      std::vector<bool> leadingJetsBtag;
      unsigned int nbtag = 0;
      for(unsigned int i =0; i < leadingJets.size() && i <= 2; i++) {
         
        const bool btag = isOfflineBtagged(ibt,leadingJets.at(i));
        leadingJetsBtag.push_back(btag);
        if(btag) {
          nbtag++;
        }
      }      
      //the 'weighting method'
      // dont fill for pat0 because weight too complicated for the time being?
        
      double jetflaveff[3];
      //int flav[3];
      JetFlavor::Codes flav(3,JetFlavor::UNDEFINED);
      for(unsigned int i =0; i <= 2; i++) {
        flav[i] = JetFlavor::code(leadingJets.at(i));
        if(flav[i] == JetFlavor::UNDEFINED) {
          std::cout << "error " << flav[i] << " " << leadingJets.at(i) << " " << isJetMatchedPartonExist[leadingJets.at(i)] << std::endl;
          throw std::exception();
        }
        //only for CSVT
        jetflaveff[i] = bTagEffOffline->eff(flav[i],"CSVT",ptJet[leadingJets.at(i)],etaJet[leadingJets.at(i)]);
      }
        


        
      //not very generic but it should work ...
      const double weight_pat1 = weight * jetflaveff[0] * jetflaveff[1] * jetflaveff[2];
      const double weight_pat2 = weight *    1.0        * jetflaveff[1] * jetflaveff[2];
      const double weight_pat3 = weight * jetflaveff[0] *   1.0         * jetflaveff[2];
      const double weight_pat4 = weight * jetflaveff[0] * jetflaveff[1] * 1.0;

      for(unsigned int i =0; i < leadingJets.size() && i <= 2; i++) {
        const int iflav = flav[i];
          
        const double pJet = TMath::Sqrt(pxJet[leadingJets.at(i)] * pxJet[leadingJets.at(i)] + pyJet[leadingJets.at(i)] * pyJet[leadingJets.at(i)] + pzJet[leadingJets.at(i)] * pzJet[leadingJets.at(i)]);

        h_debug_pt_vs_p->Fill(pJet,ptJet[leadingJets.at(i)]);

        const double weight_pat0 = weight * jetflaveff[i];
        if(isJetMatchedPartonExist[leadingJets.at(i)]) {
          h_pt_pat0_wtx.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight_pat0);
          h_p_pat0_wtx.at(ibt).at(iflav).at(i)->Fill(pJet,weight_pat0);
          h_eta_pat0_wtx.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight_pat0);
        }

        if(isJetMatchedPartonExist[leadingJets.at(0)] && isJetMatchedPartonExist[leadingJets.at(1)] && isJetMatchedPartonExist[leadingJets.at(2)]) {
          h_pt_pat1_wtx.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight_pat1);
          h_p_pat1_wtx.at(ibt).at(iflav).at(i)->Fill(pJet,weight_pat1);
          h_eta_pat1_wtx.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight_pat1);
        }

        if(isJetMatchedPartonExist[leadingJets.at(1)] && isJetMatchedPartonExist[leadingJets.at(2)]) {
          h_pt_pat2_wtx.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight_pat2);
          h_p_pat2_wtx.at(ibt).at(iflav).at(i)->Fill(pJet,weight_pat2);
          h_eta_pat2_wtx.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight_pat2);
        }

        if(isJetMatchedPartonExist[leadingJets.at(0)]  && isJetMatchedPartonExist[leadingJets.at(2)]) {
          h_pt_pat3_wtx.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight_pat3);
          h_p_pat3_wtx.at(ibt).at(iflav).at(i)->Fill(pJet,weight_pat3);
          h_eta_pat3_wtx.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight_pat3);
        }

        if(isJetMatchedPartonExist[leadingJets.at(0)] && isJetMatchedPartonExist[leadingJets.at(1)]) {
          h_pt_pat4_wtx.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight_pat4);
          h_p_pat4_wtx.at(ibt).at(iflav).at(i)->Fill(pJet,weight_pat4);
          h_eta_pat4_wtx.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight_pat4);
        }
      }
        
      //the 'hard cut method'
      if(nbtag >=2) {
        //pattern
        //pat0 -> all
        //pat1 -> 111
        //pat2 -> 011
        //pat3 -> 101
        //pat4 -> 110
        if(ibt == 0) {
          h_debug_njets_twobtags->Fill(leadingJets.size(),weight);
        }
        int pat = -1;
        if(leadingJetsBtag.at(0) && leadingJetsBtag.at(1) && leadingJetsBtag.at(2)) {
          pat = 1;
        } else if(!leadingJetsBtag.at(0) && leadingJetsBtag.at(1) && leadingJetsBtag.at(2)) {
          pat = 2;
        } else if(leadingJetsBtag.at(0) && !leadingJetsBtag.at(1) && leadingJetsBtag.at(2)) {
          pat = 3;
        } else if(leadingJetsBtag.at(0) && leadingJetsBtag.at(1) && !leadingJetsBtag.at(2)) {
          pat = 4;
        } 
          
        //std::cout << "pat " << pat << " " << leadingJets.size() << " " << leadingJetsBtag.at(0) << " " <<  leadingJetsBtag.at(1) << " " << leadingJetsBtag.at(2) << std::endl;
        for(unsigned int i =0; i < leadingJets.size(); i++) {
          //if(!(isJetMatchedPartonExist[leadingJets.at(i)]) ) continue; //not needed to increase statistics?
          const int iflav = flav[i];//jetFlavorCodeNew(leadingJets.at(i));
          const double pJet = TMath::Sqrt(pxJet[leadingJets.at(i)] * pxJet[leadingJets.at(i)] + pyJet[leadingJets.at(i)] * pyJet[leadingJets.at(i)] + pzJet[leadingJets.at(i)] * pzJet[leadingJets.at(i)]);
            
          if(iflav >= 0) {
            if(pat >= 0) {
              h_pt_pat0.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight);
              h_p_pat0.at(ibt).at(iflav).at(i)->Fill(pJet,weight);
              h_eta_pat0.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight);
            }
            if(pat == 1) {
              h_pt_pat1.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight);
              h_p_pat1.at(ibt).at(iflav).at(i)->Fill(pJet,weight);
              h_eta_pat1.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight);
            } else if(pat == 2) {
              h_pt_pat2.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight);
              h_p_pat2.at(ibt).at(iflav).at(i)->Fill(pJet,weight);
              h_eta_pat2.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight);
            } else if(pat == 3) {
              h_pt_pat3.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight);
              h_p_pat3.at(ibt).at(iflav).at(i)->Fill(pJet,weight);
              h_eta_pat3.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight);
            } else if(pat == 4) {
              h_pt_pat4.at(ibt).at(iflav).at(i)->Fill(ptJet[leadingJets.at(i)],weight);
              h_p_pat4.at(ibt).at(iflav).at(i)->Fill(pJet,weight);
              h_eta_pat4.at(ibt).at(iflav).at(i)->Fill(TMath::Abs(etaJet[leadingJets.at(i)]),weight);
            }
          } 
        }
      }
    }
  }
  ~JetFlavourFractions()
  {
    delete bTagEffOffline;
  }
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat0, h_p_pat0, h_eta_pat0;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat1, h_p_pat1, h_eta_pat1;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat2, h_p_pat2, h_eta_pat2;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat3, h_p_pat3, h_eta_pat3;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat4, h_p_pat4, h_eta_pat4;

  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat0_wtx, h_p_pat0_wtx, h_eta_pat0_wtx;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat1_wtx, h_p_pat1_wtx, h_eta_pat1_wtx;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat2_wtx, h_p_pat2_wtx, h_eta_pat2_wtx;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat3_wtx, h_p_pat3_wtx, h_eta_pat3_wtx;
  std::vector<std::vector<std::vector<TH1D*> > > h_pt_pat4_wtx, h_p_pat4_wtx, h_eta_pat4_wtx;

  

  bTagEff *bTagEffOffline;
  
  static const unsigned int nrank = 3;
  static const unsigned int npattern = 5;
  const unsigned int nbt;
  const unsigned int nflav;

  const std::string *btlabel,*btdiscr;
  const float *btcut;
  const bool *online;
  const bool *btagMatch;
  
  TH2D *h_debug_pt_vs_p;
  TH1D *h_debug_njets;
  TH1D *h_debug_njets_twobtags;
};


class quarkgluon_separator {
public:


  void Norm( std::vector<std::vector<std::vector<TH2D*> > > &h) {

   
    //    for(unsigned int x = 0; x < h.size(); x++) {
    //      for(unsigned int y = 0; y < h.at(x).size(); y++) {
    //        for(unsigned int z = 0; z < h.at(x).at(y).size(); z++) {
    //          h.at(x).at(y).at(z)->Scale(1.0/h.at(x).at(y).at(z)->Integral());
    //        }
    //      }
    //    }

    for(unsigned int x = 0; x < h.size(); x++) {
      for(unsigned int y = 0; y < h.at(x).size(); y++) {
        for(unsigned int z = 0; z < h.at(x).at(y).size(); z++) {
          for(int iBinX = 1; iBinX <= h.at(x).at(y).at(z)->GetXaxis()->GetLast(); iBinX++) {
            
            double sumy = 0.0;
            for(int iBinY = 1; iBinY <= h.at(x).at(y).at(z)->GetYaxis()->GetLast(); iBinY++) {
              sumy += h.at(x).at(y).at(z)->GetBinContent(iBinX,iBinY);
              
            }
            if(sumy > 0.0) {
              for(int iBinY = 1; iBinY <= h.at(x).at(y).at(z)->GetYaxis()->GetLast(); iBinY++) {
                const double content = h.at(x).at(y).at(z)->GetBinContent(iBinX,iBinY);
                const double newcontent = content / sumy;
                if(newcontent > 1.0) {
                  std::cout << "error: bin content larger 1.0" << std::endl;
                  throw std::exception();
                }
                const double error = h.at(x).at(y).at(z)->GetBinError(iBinX,iBinY);
                h.at(x).at(y).at(z)->SetBinContent(iBinX,iBinY, newcontent);
                h.at(x).at(y).at(z)->SetBinError(iBinX,iBinY, error / sumy);
              }
            }
          }
        }
      }
    }
  }
  void Normalise()
  {
    Norm(h_nchargedConstituents_vs_pt);
    Norm(h_nConstituents_vs_pt);
    Norm(h_prel_vs_pt);
    Norm(h_jetarea_vs_pt);
    Norm(h_jetarea_over_pt_vs_pt);
    Norm(h_ptd_vs_pt);
    Norm(h_csv_vs_pt);
  }













  quarkgluon_separator() {
    file = new TFile("/data/user/jbehr/quarkgluon_sep.root");
    
    std::vector<std::string> flavlabel;
    flavlabel.push_back("slim");
    flavlabel.push_back("fat");

    const unsigned int nmode = 2;
    const std::string modelabel[nmode] = {"all", "btag"};


    const unsigned int nbt = 1;
    std::string btlabel[nbt] = {"CSVT"};

    const unsigned int nflav = flavlabel.size();
  

    std::vector<std::vector<std::vector<TH2D*> > > initvec2d
      = std::vector<std::vector<std::vector<TH2D*> > > (nbt, std::vector<std::vector<TH2D*> > (nflav, std::vector<TH2D*> (nmode, NULL)) );
    
    h_nchargedConstituents_vs_pt =  initvec2d;
    h_nConstituents_vs_pt=  initvec2d;
    h_prel_vs_pt=  initvec2d;
    h_jetarea_vs_pt=  initvec2d;
    h_jetarea_over_pt_vs_pt=  initvec2d;
    
    h_ptd_vs_pt =  initvec2d;
    h_csv_vs_pt =  initvec2d;

    file->cd();




    for(unsigned int ibt = 0; ibt < nbt; ibt++) {
      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_jetarea_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_jetarea_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //   h_jetarea_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_jetarea_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        

        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_jetarea_over_pt_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_jetarea_over_pt_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //  h_jetarea_over_pt_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_jetarea_over_pt_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        

        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_nchargedConstituents_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nchargedConstituents_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //  h_nchargedConstituents_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_nchargedConstituents_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        



        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_nConstituents_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_nConstituents_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //   h_nConstituents_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_nConstituents_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        



        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_prel_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_prel_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //  h_prel_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_prel_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        


        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_ptd_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_ptd_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //   h_ptd_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_ptd_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        
        for (unsigned int imode=0; imode<nmode; ++imode) {
          std::string s(Form("jet_constituents/h_csv_vs_pt_%s_%s_%s",btlabel[ibt].c_str(),flavlabel[iflav].c_str(),modelabel[imode].c_str()));
          h_csv_vs_pt.at(ibt).at(iflav).at(imode)  = (TH2D*)(file->Get(TString(s.c_str())))->Clone();
          //   h_ptd_vs_pt.at(ibt).at(iflav).at(imode)->Scale(1.0/h_ptd_vs_pt.at(ibt).at(iflav).at(imode)->Integral());
        }
        



      }
    }


    this->Normalise();



  }

  double prob(TH2D *h, double x, double y) {
    double r = 0.0;
    if(h->GetXaxis()->GetBinLowEdge(1) < x && h->GetXaxis()->GetBinUpEdge(h->GetXaxis()->GetLast()) > x )
      {
        if(h->GetYaxis()->GetBinLowEdge(1) < y && h->GetYaxis()->GetBinUpEdge(h->GetYaxis()->GetLast()) > y )
          {
            r = h->Interpolate(x,y);
            
          }
      }
    return r;
  }
  unsigned int GetSlimFat(unsigned int iJet, std::string *flavlabel, unsigned int iflav_external) {

    int iflav = -1;

    //mapping jets to fat and slim categories
    if(flavlabel[iflav_external] == "c" || flavlabel[iflav_external] == "b") {
      iflav = 0;
    } else  if(flavlabel[iflav_external] == "cc" || flavlabel[iflav_external] == "bb") {
      iflav = 1;
    } else if(abs(partonFlavorJet[iJet]) == 21) {
      iflav = 1;
    } else {
      iflav = 0;
    }

    if(iflav == -1) {
      std::cout << "error: unknown slim/fat type found" << std::endl;
      throw std::exception();
    }

    return iflav;

  }


  double GetDiscriminant(unsigned int iJet) {
    
    const double jetpt = ptJet[iJet];
    // numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet];
    //areaJet[iJet] areaJet[iJet]/jetpt ptdJet[iJet]
    TVector3 vec_jet(pxJet[iJet], pyJet[iJet], pzJet[iJet]);

    double sum_px = 0.0, sum_py = 0.0, sum_pz = 0.0;
    for(int iTrack = 0; iTrack < numberOfTracksJet[iJet]; iTrack++) {
      sum_px += trackPxJet[iJet][iTrack];
      sum_py += trackPyJet[iJet][iTrack];
      sum_pz += trackPzJet[iJet][iTrack];
    }
    TVector3 vec_tracks(sum_px, sum_py, sum_pz);
 
    
    const double p_rel_tracks = vec_tracks.Mag() * TMath::Cos(vec_jet.Angle(vec_tracks));
    const double p_rel = p_rel_tracks / vec_jet.Mag();



    double gluon_prob = 1.0;
    double quark_prob = 1.0;



    const unsigned int imode = 0;
    const unsigned int ibt = 0;
    //gluon_prob *= prob(h_jetarea_vs_pt.at(ibt).at(1).at(imode),jetpt, areaJet[iJet]);
    //gluon_prob *=  prob(h_jetarea_over_pt_vs_pt.at(ibt).at(1).at(imode),jetpt, areaJet[iJet] / jetpt);
    gluon_prob *=  prob(h_nchargedConstituents_vs_pt.at(ibt).at(1).at(imode),jetpt, numberOfChargedConstituentsInJet[iJet]);
    gluon_prob *=  prob(h_nConstituents_vs_pt.at(ibt).at(1).at(imode),jetpt, numberOfConstituentsInJet[iJet]);
    gluon_prob *=  prob(h_prel_vs_pt.at(ibt).at(1).at(imode),jetpt, p_rel);
    gluon_prob *=  prob(h_ptd_vs_pt.at(ibt).at(1).at(imode),jetpt, ptdJet[iJet]);
    gluon_prob *=  prob(h_csv_vs_pt.at(ibt).at(1).at(imode),jetpt,combSVBJetTag[iJet]);
    
    
    //quark_prob *=  prob(h_jetarea_vs_pt.at(ibt).at(0).at(imode),jetpt, areaJet[iJet]);
    //quark_prob *=  prob(h_jetarea_over_pt_vs_pt.at(ibt).at(0).at(imode),jetpt, areaJet[iJet] / jetpt);
    quark_prob *=  prob(h_nchargedConstituents_vs_pt.at(ibt).at(0).at(imode),jetpt, numberOfChargedConstituentsInJet[iJet]);
    quark_prob *=  prob(h_nConstituents_vs_pt.at(ibt).at(0).at(imode),jetpt, numberOfConstituentsInJet[iJet]);
    quark_prob *=  prob(h_prel_vs_pt.at(ibt).at(0).at(imode),jetpt, p_rel);
    quark_prob *=  prob(h_ptd_vs_pt.at(ibt).at(0).at(imode),jetpt, ptdJet[iJet]);
    quark_prob *=  prob(h_csv_vs_pt.at(ibt).at(1).at(imode),jetpt,combSVBJetTag[iJet]);
    
   


    const double discr = (quark_prob + gluon_prob) > 0.0 ? quark_prob / (quark_prob + gluon_prob) : 0.0;
    

    if(discr > 1.0 || discr < 0.0)  {
      std::cout << "error: discriminator > 1 or < 0 : " << discr << std::endl;
      throw std::exception();

    }

    return discr;

  }



  ~quarkgluon_separator() {
    delete file;
  }
  
  TFile *file;

  std::vector<std::vector<std::vector<TH2D*> > > h_nchargedConstituents_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_nConstituents_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_prel_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_jetarea_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_jetarea_over_pt_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_ptd_vs_pt;
  std::vector<std::vector<std::vector<TH2D*> > > h_csv_vs_pt;
  


};


class jetpartonmatchingefficiency {
public:
  jetpartonmatchingefficiency(TFile *hout) 
  {

    const int nch=63;
    //double low= 0, high=1000;
    
    const double pt_binning_1d[nch+1] = {2.00000e+01, 2.46667e+01, 2.93333e+01, 3.40000e+01, 3.86667e+01, 4.33333e+01, 4.80000e+01, 5.26667e+01, 5.73333e+01, 6.20000e+01, 6.66667e+01, 7.13333e+01, 7.60000e+01, 8.06667e+01, 8.53333e+01, 9.00000e+01, 9.46667e+01, 9.93333e+01, 1.04000e+02, 1.08667e+02, 1.13333e+02, 1.18000e+02, 1.22667e+02, 1.27333e+02, 1.32000e+02, 1.36667e+02, 1.41333e+02, 1.46000e+02, 1.50667e+02, 1.55333e+02, 1.60000e+02, 1.64667e+02, 1.69333e+02, 1.74000e+02, 1.78667e+02, 1.83333e+02, 1.88000e+02, 1.92667e+02, 1.97333e+02, 2.02000e+02, 2.06667e+02, 2.11333e+02, 2.16000e+02, 2.20667e+02, 2.25333e+02, 2.30000e+02, 2.34667e+02, 2.39333e+02, 2.44000e+02, 2.48667e+02, 2.53333e+02, 2.58000e+02, 2.62667e+02, 2.67333e+02, 2.72000e+02, 2.76667e+02, 2.81333e+02, 2.86000e+02, 2.90667e+02, 2.95333e+02, 3.00000e+2, 400.0, 550.0, 800.0};
    
    

    const std::string foldername = "jet_parton_matching_efficiency";

    hout->mkdir(TString(foldername.c_str()));
    hout->cd(TString(foldername.c_str()));
    
    
    std::string name1d = "eff_pt";
    profile_1d["pt"] = new TProfile(name1d.c_str(), name1d.c_str(), nch, pt_binning_1d);

    for(int iJet = 0; iJet < nj; iJet++) {
      std::stringstream name;
      name << "eff_pt_j" << iJet;
      profile_1d[name.str()] = new TProfile(name.str().c_str(), name.str().c_str(), nch, pt_binning_1d);

    }
    
    name1d = "avg_npartons_vs_pt_matched";
    profile_1d[name1d] = new TProfile(name1d.c_str(), name1d.c_str(), nch, pt_binning_1d);
    
    name1d = "avg_npartons_vs_pt_unmatched";
    profile_1d[name1d] = new TProfile(name1d.c_str(), name1d.c_str(), nch, pt_binning_1d);
    

    name1d = "avg_npartons_RMS_vs_pt_matched";
    profile_1d[name1d] = new TProfile(name1d.c_str(), name1d.c_str(), nch, pt_binning_1d);
    profile_1d[name1d]->SetErrorOption("s");

    name1d = "avg_npartons_RMS_vs_pt_unmatched";
    profile_1d[name1d] = new TProfile(name1d.c_str(), name1d.c_str(), nch, pt_binning_1d);
    profile_1d[name1d]->SetErrorOption("s");






    name1d = "eff_eta";
    profile_1d["eta"] = new TProfile(name1d.c_str(), name1d.c_str(), 30, -2.6, 2.6);

    name1d = "eff_nPUTruth";
    profile_1d["nPUTruth"] = new TProfile(name1d.c_str(), name1d.c_str(), 80, 0.0, 80);

    std::string name2d = "eff_eta_vs_pt";
    profile_2d["eta_vs_pt"] = new TProfile2D(name2d.c_str(), name2d.c_str(), nch, pt_binning_1d, 10, -2.6, 2.6);


    name2d = "eff_nPUTruth_vs_eta";
    profile_2d["nPUTruth_vs_eta"] = new TProfile2D(name2d.c_str(), name2d.c_str(), 40, 0.0, 80.0, 10, -2.6, 2.6);


    name2d = "eff_nPUTruth_vs_pt";
    profile_2d["nPUTruth_vs_pt"] = new TProfile2D(name2d.c_str(), name2d.c_str(), nch, pt_binning_1d, 40, 0.0, 80.0);
 

    hout->cd();
  }
  void Fill(double weight)
  {
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      if ( (fabs(etaJet[iJet]) > 2.6) ) continue;
      if ( ptJet[iJet] < 20. ) continue;
      if(!puJetIDLoose[iJet])
        continue;
      if(!jetIDLoose[iJet] ) continue;
      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
      
      const double matched = isJetMatchedPartonExist[iJet] ? 1.0 : 0.0;
      
      profile_1d["pt"]->Fill(ptJet[iJet], matched, weight);
      profile_1d["eta"]->Fill(etaJet[iJet], matched, weight);
      profile_1d["nPUTruth"]->Fill(nPUTruth, matched, weight);
      profile_2d["eta_vs_pt"]->Fill(ptJet[iJet],etaJet[iJet],matched,weight);
      profile_2d["nPUTruth_vs_eta"]->Fill(etaJet[iJet],nPUTruth,matched,weight);
      profile_2d["nPUTruth_vs_pt"]->Fill(ptJet[iJet],nPUTruth,matched,weight);
      
      if(iJet < nj) {
        std::stringstream name;
        name << "eff_pt_j" << iJet;
        profile_1d[name.str()]->Fill(ptJet[iJet], matched, weight);
      }
      const TLorentzVector vec_jet(pxJet[iJet], pyJet[iJet], pzJet[iJet], energyJet[iJet]);
        

      int npartons = 0;
      for(int iParton = 0; iParton < numberOfPartonsFull && iParton < 500; iParton++) {
        if(statusPartonFull[iParton] != 2) continue;
        if(!(ptPartonFull[iParton] > 0.0)) continue;
        TLorentzVector vec_parton(1.0,1.0,1.0,1.0);
        if(ptPartonFull[iParton] > 0.0) {
        } else {
          std::cout << "error " << numberOfPartonsFull << " " << iParton << " " << ptPartonFull[iParton] << " " << etaPartonFull[iParton]
                    << " " <<  phiPartonFull[iParton] << " " <<  energyPartonFull[iParton]
                    << " " << stablePartonStatusFull
                    << " " << statusPartonFull[iParton]
                    << " " << pdgIdPartonFull[iParton]
            
                    << std::endl;
          throw std::exception();
        }

        vec_parton.SetPtEtaPhiE(ptPartonFull[iParton], etaPartonFull[iParton], phiPartonFull[iParton], energyPartonFull[iParton]);
        if(vec_parton.DrEtaPhi(vec_jet) < 0.3) {
          npartons++;
        }
      }

      if(isJetMatchedPartonExist[iJet]) {
        profile_1d["avg_npartons_vs_pt_matched"]->Fill(ptJet[iJet], npartons, weight);
        profile_1d["avg_npartons_RMS_vs_pt_matched"]->Fill(ptJet[iJet], npartons, weight);
        
      } else {
        profile_1d["avg_npartons_vs_pt_unmatched"]->Fill(ptJet[iJet], npartons, weight);
        profile_1d["avg_npartons_RMS_vs_pt_unmatched"]->Fill(ptJet[iJet], npartons, weight);
      }
      
    }
  }
private:
  static const int nj = 5;
  std::map<std::string, TProfile*> profile_1d;
  std::map<std::string, TProfile2D*> profile_2d;
  
};



std::vector<int> trijetselection(const double *jetPtMin, const double etacut, const double deltaEtaCut12) {
  std::vector<int> leadingJets;
  std::vector<int> returnleadingjets;
  // find set of leading jets
  int nJet = 0;
  // loop over the jets
  for (int iJet=0; iJet<numberOfJets; ++iJet) {
    if ( (fabs(etaJet[iJet]) > etacut) ) continue;
    if ( ptJet[iJet] < 20. ) continue;
    if(!puJetIDLoose[iJet])
      continue;
    if(!jetIDLoose[iJet] ) continue;
    if (nJet >= 3) break;
    if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
    if ( (ptJet[iJet] > jetPtMin[nJet]) && (ptJet[iJet] < 3500.0) ) {
      leadingJets.push_back(iJet);
      ++nJet;
    }
  }
  
  bool flav_match_ok = true;
  for(unsigned int i =0; i < leadingJets.size(); i++) {
    if(i <= 2) {
      if(flav_match_ok && JetFlavor::code(leadingJets.at(i)) == JetFlavor::UNDEFINED) 
        flav_match_ok = false;
    }
  }

  if(flav_match_ok && leadingJets.size() >= 3 ) {
    double deltaR12 = -1.0;
    double deltaR23 = -1.0;
    double deltaR13 = -1.0;
    
    double deltaEta = 9999.9; //between leading jets
    if(nJet >= 2)
      {
        const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[1]];
        float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
        if (dphij>3.1415926) dphij -= 2*3.1415926;
        if (dphij<-3.1415926) dphij += 2*3.1415926;      
        deltaR12 = sqrt( dphij*dphij + detaj*detaj );
        deltaEta = TMath::Abs(etaJet[leadingJets[1]] - etaJet[leadingJets[0]]);
      }
    if(nJet >= 3)
      {
        const float detaj = etaJet[leadingJets[0]] - etaJet[leadingJets[2]];
        float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[0]];
        if (dphij>3.1415926) dphij -= 2*3.1415926;
        if (dphij<-3.1415926) dphij += 2*3.1415926;      
        deltaR13 = sqrt( dphij*dphij + detaj*detaj );
      }
    if(nJet >= 3)
      {
        const float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[2]];
        float dphij =  phiJet[leadingJets[2]] - phiJet[leadingJets[1]];
        if (dphij>3.1415926) dphij -= 2*3.1415926;
        if (dphij<-3.1415926) dphij += 2*3.1415926;      
        deltaR23 = sqrt( dphij*dphij + detaj*detaj );
      }

    if(deltaR12 > 1.0 && deltaR23 > 1.0 && deltaR13 > 1.0 && deltaEta < deltaEtaCut12) {
      returnleadingjets = leadingJets;
    }

  }
  return returnleadingjets;
}

void readmistag(const int pileupwtx = 0) {

  TH1::SetDefaultSumw2(kTRUE);
  TH2::SetDefaultSumw2(kTRUE);
  const bool TRIGSELECT = false;

  std::cout << "TRIGSELECT " << TRIGSELECT << std::endl;
  std::cout << "pileupwtx " << pileupwtx << std::endl;


  // here we require three btagged jets

  // open an ntuple file
  std::cout << " starting..." << std::endl;

 
  std::vector<std::string> genericTriggerList;

  
  const bool _doMC = true;


  TTree* hbbtree;
  float lumiScaleFac = 1.0;
 

  // chain mode
  TChain f("hbbanalysis/HBBTo4B");
  TFileCollection fc("dum","","theMergeList.txt");

  // extract generic trigger list and number of input events from ntuple files
  lumiScaleFac = 1.0;
  



  std::vector<std::string> HLTBtagTriggerObjectList;

  getHbbMetadata(fc,genericTriggerList,HLTBtagTriggerObjectList,1000.,lumiScaleFac,_doMC);
  


  f.AddFileInfoList((TCollection*) fc.GetList());
  hbbtree = &f;

  TIter next( (TCollection*) fc.GetList() );
  TFileInfo* tfi = (TFileInfo*) next();
  TString firstFile( tfi->GetFirstUrl()->GetFile() );

  //  quarkgluon_separator qgdiscriminant;


  TFile* hout = new TFile("BtagMatrix.root","recreate");
  PileUpWeight PileUpWeightObj(_doMC,pileupwtx, hout);


  jetpartonmatchingefficiency jpm(hout);

  

  //const bool oldAddDir = TH1::AddDirectoryStatus();
  TH1::AddDirectory(true);

  TH1::SetDefaultSumw2(kTRUE);

  Int_t nentries = (Int_t) hbbtree->GetEntries();
  
  std::cout << "Number of events in ntuple: " << nentries << std::endl;

  //#include "readhbbset.C"
#include  "Analysis/Utilities/interface/HbbNtuple.cc"

  //disable everything
  hbbtree->SetBranchStatus("*",0);
  hbbtree->SetBranchStatus("MatchedPartonFlavor",1);
  //enable certain needed branches
  hbbtree->SetBranchStatus("nJetConstituents",1);
  hbbtree->SetBranchStatus("nJetChargedConstituents",1);
  hbbtree->SetBranchStatus("NumberOfJets",1);
  hbbtree->SetBranchStatus("JetEta",1);
  hbbtree->SetBranchStatus("JetPhi",1);
  hbbtree->SetBranchStatus("JetPt",1);
  hbbtree->SetBranchStatus("JetPx",1);
  hbbtree->SetBranchStatus("JetEnergy",1);
  hbbtree->SetBranchStatus("JetPy",1);
  hbbtree->SetBranchStatus("JetPz",1);
  hbbtree->SetBranchStatus("JetArea",0);
  hbbtree->SetBranchStatus("ptdArea",0);

  hbbtree->SetBranchStatus("IsJetMatchedPartonExist",1);
  hbbtree->SetBranchStatus("statusJetMatchedParton",1);
  hbbtree->SetBranchStatus("ptJetMatchedParton",1);
  hbbtree->SetBranchStatus("etaJetMatchedParton",1);



  hbbtree->SetBranchStatus("HflContentJet",0);
  hbbtree->SetBranchStatus("PartonFlavorJet",1);

  hbbtree->SetBranchStatus("NumberOfPU",1);
  hbbtree->SetBranchStatus("NumberOfPUInTime",1);

  hbbtree->SetBranchStatus("CombSVBJetTag",1);
  hbbtree->SetBranchStatus("TCHEBJetTag",0);
  hbbtree->SetBranchStatus("TCHPBJetTag",0);
  hbbtree->SetBranchStatus("nTCHEBJetTag",0);
  hbbtree->SetBranchStatus("nTCHPBJetTag",0);
  hbbtree->SetBranchStatus("SVHPBJetTag",0);
 
  hbbtree->SetBranchStatus("nPUTruth",1);



  //  hbbtree->SetBranchStatus("",1);
  //  hbbtree->SetBranchStatus("",1);
  //   hbbtree->SetBranchStatus("",1);
  //  hbbtree->SetBranchStatus("",1);
  //   hbbtree->SetBranchStatus("",1);




  hbbtree->SetBranchStatus("NumberOfTracksJet",1);
  hbbtree->SetBranchStatus("JetTrackPx",1);
  hbbtree->SetBranchStatus("JetTrackPy",1);
  hbbtree->SetBranchStatus("JetTrackPz",1);
 
  if(hbbtree->GetBranch("NumberOfPartonsFull") != NULL){
    hbbtree->SetBranchStatus("*Full",1);
  } 





  std::cout << "first file " << firstFile << std::endl;

  
 
    

  //   hbbtree->SetBranchAddress("BContentJet4", BContentJet4);
  //       hbbtree->SetBranchAddress("CContentJet4", CContentJet4);
  //   hbbtree->SetBranchAddress("BContentJet5", BContentJet5);
  //       hbbtree->SetBranchAddress("CContentJet5", CContentJet5);
  //       hbbtree->SetBranchAddress("BContentJet6", BContentJet6);
  //       hbbtree->SetBranchAddress("CContentJet6", CContentJet6);
  //       hbbtree->SetBranchAddress("BContentJet7", BContentJet7);
  //       hbbtree->SetBranchAddress("CContentJet7", CContentJet7);
  //       hbbtree->SetBranchAddress("BContentJet10", BContentJet10);
  //       hbbtree->SetBranchAddress("CContentJet10", CContentJet10);


  hbbtree->SetBranchStatus("PUJetMVA",0);
  hbbtree->SetBranchStatus("PUJetIDLoose",1);
  hbbtree->SetBranchStatus("PUJetIDMedium",0);
  hbbtree->SetBranchStatus("PUJetIDTight",0);


  hbbtree->SetBranchStatus("JetIDLoose",1);


  hbbtree->SetBranchStatus("BContentJet3", 1);
  hbbtree->SetBranchStatus("CContentJet3", 1);

  hbbtree->SetBranchStatus("BContentJet5", 0);
  hbbtree->SetBranchStatus("CContentJet5", 0);
   
  
  for(int i =0;i < 100; i++)
    {
      //  puJetMVA[i] = 0.0f;
      //           puJetIDLoose[i] = false;
      //           puJetIDMedium[i] = false;
      //           puJetIDTight[i] = false;
      //           BContentJet3[i] = 0;
      //           CContentJet3[i] = 0;
      BContentJet4[i] = 0;
      CContentJet4[i] = 0;
      //    BContentJet5[i] = 0;
      //       CContentJet5[i] = 0;
      BContentJet6[i] = 0;
      CContentJet6[i] = 0;
      BContentJet7[i] = 0;
      CContentJet7[i] = 0;
      BContentJet10[i] = 0;
      CContentJet10[i] = 0;
    }
   
  




  //int nch=200;
  const int nch=63;
  //double low= 0, high=1000;
 
  const double pt_binning_1d[nch+1] = {2.00000e+01, 2.46667e+01, 2.93333e+01, 3.40000e+01, 3.86667e+01, 4.33333e+01, 4.80000e+01, 5.26667e+01, 5.73333e+01, 6.20000e+01, 6.66667e+01, 7.13333e+01, 7.60000e+01, 8.06667e+01, 8.53333e+01, 9.00000e+01, 9.46667e+01, 9.93333e+01, 1.04000e+02, 1.08667e+02, 1.13333e+02, 1.18000e+02, 1.22667e+02, 1.27333e+02, 1.32000e+02, 1.36667e+02, 1.41333e+02, 1.46000e+02, 1.50667e+02, 1.55333e+02, 1.60000e+02, 1.64667e+02, 1.69333e+02, 1.74000e+02, 1.78667e+02, 1.83333e+02, 1.88000e+02, 1.92667e+02, 1.97333e+02, 2.02000e+02, 2.06667e+02, 2.11333e+02, 2.16000e+02, 2.20667e+02, 2.25333e+02, 2.30000e+02, 2.34667e+02, 2.39333e+02, 2.44000e+02, 2.48667e+02, 2.53333e+02, 2.58000e+02, 2.62667e+02, 2.67333e+02, 2.72000e+02, 2.76667e+02, 2.81333e+02, 2.86000e+02, 2.90667e+02, 2.95333e+02, 3.00000e+2, 400.0, 550.0, 800.0};











  const int nmode=3;
  std::string modelabel[nmode] = {"btag","all","eff"};



  //   const int nbt = 8;
  //   string btlabel[nbt] = { "TCHEL","TCHEM","TCHE6","TCHE10", "TCHPM", "TCHPT", "TCHP6", "TCHP10" };
  //   string btdiscr[nbt] = { "TCHE", "TCHE", "TCHE", "TCHE", "TCHP", "TCHP", "TCHP", "TCHP" };
  //   float btcut[nbt] = { 1.7, 3.3, 6, 10, 1.93, 3.41, 6, 10};
  const int nbt = 10;
  //string btlabel[nbt] = { "TCHPT", "TCHP6", "tTCHPT", "btTCHPT", "tTCHP6", "btTCHP6"};
  std::string btlabel[nbt] = { "TCHPT", "TCHP6", "tTCHPT", "btTCHPT", "tTCHP6", "btTCHP6", "CSVT", "tCSVT", "btCSVT", "SSVHPT" };
  std::string btdiscr[nbt] = { "TCHP", "TCHP", "TCHP", "TCHP", "TCHP", "TCHP", "CSV", "CSV", "CSV", "SSVHP"};
  float btcut[nbt] = { 3.41, 6, 3.41, 3.41, 6, 6, 0.898, 0.898, 0.898, 2};
  bool  online[nbt] = { false, false, true, true, true, true, false, true, true, false};
  bool  btagMatch[nbt] = { false, false, false, true, false, true, false, false, true, false};


  //  const int nflav=5;
  //     std::string flavlabel[nflav] = {
  //       "udsg","c","b", "cc", "bb"
   
  //     };
    
  const std::vector<JetFlavor::Code> jetFlavorCodes  = JetFlavor::getList();
  const unsigned int nflav = jetFlavorCodes.size();
  std::vector<std::string> flavlabel;          // flavor codes
  for(std::vector<JetFlavor::Code>::const_iterator it = jetFlavorCodes.begin();
      it != jetFlavorCodes.end(); ++it) {
    flavlabel.push_back( JetFlavor::toString( *it ) );
  }
    


  TH1D* h_nvtx[nbt][nflav][nmode];
  TH1D* h_nPU[nbt][nflav][nmode];

  TH1D* h_nvtx_JetIDLoose[nbt][nflav][nmode];
  TH1D* h_nPU_JetIDLoose[nbt][nflav][nmode];

  TH1D* h_nvtx_trijet[nbt][nflav][nmode];
  TH1D* h_nPU_trijet[nbt][nflav][nmode];
  
  TH1D* hpt_fine[nbt][nflav][nmode];


  //number of flavor definitions
  const int n_flav_def = 14*2;
 
  TH1D* hpt[nbt][nflav][nmode][n_flav_def];
 



  TH1D* h_bcontent[nbt][nflav][nmode];
  TH1D* h_ccontent[nbt][nflav][nmode];
  


  TH1D* hpt_trijet_doubletag[nbt][nflav][nmode];
 

  TH1D* heta[nbt][nflav][nmode];
  // TH2D* hpteta[nbt][nflav][nmode];//used for default efficiency functions
  const int nchargedconstituents = 5;
  //const int nconstjet[nchargedconstituents+1] = {0,6,12,18,24,999999};
  

  TH2D* hpteta[nbt][nflav][nmode][nchargedconstituents];
  


  TH1D* h_btagdiscriminator[nbt][nflav][nmode];

  TH1D *h_csv_discr = new TH1D("h_csv_discr","h_csv_discr", 100, -20.0, 1.1);
  TH1D *h_csv_discr_jetisolation = new TH1D("h_csv_discr_jetisolation","h_csv_discr_jetisolation", 100, -20.0, 1.1);

 
  TProfile* hpteta_pt_bincenter[nbt][nflav][nmode];
  TProfile* hpteta_eta_bincenter[nbt][nflav][nmode];

  //  const double pt_binning_b[35] = {0,9.67742,19.3548,29.0323,38.7097,48.3871,58.0645,67.7419,77.4194,87.0968,96.7742,106.452,116.129,125.806,135.484,145.161,154.839,164.516,174.194,183.871,193.548,203.226,212.903,222.581,232.258,241.935,251.613,261.29,270.968,280.645,290.323,300,400.0,550.0,800.0};


//   //
//   const double pt_binning_b[22] = {
//     0.0,
//     2.00000e+01, 
//     2.40510e+01, 
//     2.89225e+01, 
//     3.47808e+01, 
//     4.18256e+01, 
//     5.02973e+01, 
//     6.04850e+01, 
//     7.27363e+01, 
//     8.74690e+01, 
//     1.05186e+02, 
//     1.26491e+02, 
//     1.52112e+02, 
//     1.82922e+02, 
//     2.19973e+02, 
//     2.64528e+02, 
//     3.18108e+02, 
//     3.82541e+02, 
//     4.60025e+02, 
//     5.53202e+02, 
//     6.65253e+02, 
//     800.0
//   };

    //used for default efficiency functions
  const double pt_binning_b[40] = {
    0
    , 20.0
    ,29.0323
    ,38.7097
    ,48.3871
    ,58.0645
    ,67.7419
    ,77.4194
    ,87.0968
    ,96.7742
    ,106.452
    ,116.129
    ,125.806
    ,135.484
    ,145.161
    ,154.839
    ,164.516
    ,174.194
    ,183.871
    ,193.548
    ,203.226
    ,212.903
    ,222.581
    ,232.258
    ,241.935
    ,251.613
    ,261.29
    ,270.968
    ,280.645
    ,290.323,

    3.00000e+02, 
    3.34542e+02, 
    3.73062e+02, 
    4.16017e+02, 
    4.63917e+02, 
    5.17333e+02, 
    5.76900e+02, 
    6.43325e+02, 
    7.17398e+02,
    800.0


  };

  


  //  const double pt_binning[18] = { 0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300, 500.0, 800.0};


  //   //used for default efficiency functions
  //   const double pt_binning[25] = { 0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,

  //  3.00000e+02, 
  // 3.34542e+02, 
  // 3.73062e+02, 
  // 4.16017e+02, 
  // 4.63917e+02, 
  // 5.17333e+02, 
  // 5.76900e+02, 
  // 6.43325e+02, 
  // 7.17398e+02,
  // 800.0


  // };



  const double pt_binning[12] = {
    0.0,
    2.00000e+01, 
    2.89225e+01, 
    4.18256e+01, 
    6.04850e+01, 
    8.74690e+01, 
    1.26491e+02, 
    1.82922e+02, 
    2.64528e+02, 
    3.82541e+02, 
    5.53202e+02,

    800.0

  };





  const int netabins = 5;
  //    const int netabins_fine = 10;


  //const double eta_binning[7] = {0.00000e+00, 4.00000e-01, 8.00000e-01, 1.20000e+00, 1.60000e+00, 2.2, 2.6};
  //const double eta_binning_default[netabins+1] = {0.00000e+00, 4.00000e-01, 8.00000e-01, 1.20000e+00, 1.60000e+00, 2.2, 2.6};
  const double eta_binning_default[6] = {0, 0.52, 1.04, 1.56, 2.08, 2.6};
  //  const double eta_binning_fine[netabins_fine+1] = {
  //       0.00000e+00, 
  //       2.40000e-01, 
  //       4.80000e-01, 
  //       7.20000e-01, 
  //       9.60000e-01, 
  //       1.20000e+00, 
  //       1.44000e+00, 
  //       1.68000e+00, 
  //       1.92000e+00, 
  //       2.16000e+00, 
  //       2.4
  //     };
  
  // const double eta_binning[6] = {0, 0.52, 1.04, 1.56, 2.08, 2.6};
  //   const double eta_binning_coarselargeeta[6] = {
  //       0.00000e+00, 
  //       4.50000e-01, 
  //       9.00000e-01, 
  //       1.35000e+00, 
  //       1.8,
  //       2.6
  //     };
  const double *eta_binning = eta_binning_default;//eta_binning_coarselargeeta;//eta_binning_default;

  //  const double eta_binning[6] = _{0.00000e+00, 3.00000e-01, 8.00000e-01, 1.30000e+00, 2.00000e+00, 2.6};
  //   //const double eta_binning[6] = {0.00000e+00, 5.20000e-01, 1.04000e+00, 1.56000e+00, 2.08000e+00, 2.6};

  //    const double eta_binning_udscg[4] = {0.00000e+00, 8.0e-01, 1.8e+00, 2.6};

  TH1D *h_all[nbt][nflav][nmode];

  //create (charged) constituents histograms
  JetConstituents jetconsti(hout,nbt,nflav,btlabel,&flavlabel.front());

  //create histogram object which handles the determination of the flavour fractions in the trijet sample

  const unsigned int nbt_jetfrac = 1;
  std::string btlabel_jetfrac[nbt_jetfrac] = { "CSVT" };
  std::string btdiscr_jetfrac[nbt_jetfrac] = { "CSV" };
  float btcut_jetfrac[nbt_jetfrac] = { 0.898 };
  bool  online_jetfrac[nbt_jetfrac] = { false };
  bool  btagMatch_jetfrac[nbt_jetfrac] = {false};
  


  JetFlavourFractions jetfrac_mediummass("mediummass",hout,nbt_jetfrac,nflav,&flavlabel.front(),btlabel_jetfrac,btdiscr_jetfrac,btcut_jetfrac,online_jetfrac,btagMatch_jetfrac);
  JetFlavourFractions jetfrac_highmass("highmass",hout,nbt_jetfrac,nflav,&flavlabel.front(),btlabel_jetfrac,btdiscr_jetfrac,btcut_jetfrac,online_jetfrac,btagMatch_jetfrac);
  JetFlavourFractions jetfrac_veryhighmass("veryhighmass",hout,nbt_jetfrac,nflav,&flavlabel.front(),btlabel_jetfrac,btdiscr_jetfrac,btcut_jetfrac,online_jetfrac,btagMatch_jetfrac);
  

  // create histograms
  
  for (int bt=0; bt<nbt; ++bt) {
    for (int ich=0; ich<nchargedconstituents; ++ich) {
      hout->mkdir(Form("btag_%s_%i",btlabel[bt].c_str(),ich));
      hout->cd(Form("btag_%s_%i",btlabel[bt].c_str(),ich));

  

      for (unsigned int iflav=0; iflav<nflav; ++iflav) {
        for (int mode=0; mode<nmode; ++mode) {
        







          std::string sall(Form("integral_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_all[bt][iflav][mode]  = new TH1D(sall.c_str(), sall.c_str(), 1, 0.0, 1);

          std::string snvtx(Form("nvtx_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_nvtx[bt][iflav][mode] = new TH1D(snvtx.c_str(),snvtx.c_str(), 50, 0.0, 50.0);

          std::string snPu(Form("nPu_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_nPU[bt][iflav][mode] = new TH1D(snPu.c_str(),snPu.c_str(), 100, 0.0, 200.0);
        
          std::string snvtx_JetIDLoose(Form("nvtx_%s_%s_%s_JetIDLoose",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_nvtx_JetIDLoose[bt][iflav][mode] = new TH1D(snvtx_JetIDLoose.c_str(),snvtx_JetIDLoose.c_str(), 50, 0.0, 50.0);

          std::string snPu_JetIDLoose(Form("nPu_%s_%s_%s_JetIDLoose",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_nPU_JetIDLoose[bt][iflav][mode] = new TH1D(snPu_JetIDLoose.c_str(),snPu_JetIDLoose.c_str(), 100, 0.0, 200.0);
        




          std::string snvtx_trijet(Form("nvtx_%s_%s_%s_trijet",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_nvtx_trijet[bt][iflav][mode] = new TH1D(snvtx_trijet.c_str(),snvtx_trijet.c_str(), 50, 0.0, 50.0);

          std::string snPu_trijet(Form("nPu_%s_%s_%s_trijet",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_nPU_trijet[bt][iflav][mode] = new TH1D(snPu_trijet.c_str(),snPu_trijet.c_str(), 100, 0.0, 200.0);
        



          //h_btagdiscriminator
          std::string sbtagdiscriminator(Form("btagdiscr_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          h_btagdiscriminator[bt][iflav][mode] = new TH1D(sbtagdiscriminator.c_str(),sbtagdiscriminator.c_str(), 50, -0.05, 1.05);
        
          for(int i_flav_def =0; i_flav_def < n_flav_def; i_flav_def++)
            {
              std::string s(Form("%s_%s_%s_%d",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str(),i_flav_def));
              hpt[bt][iflav][mode][i_flav_def] = new TH1D(s.c_str(),s.c_str(), nch, pt_binning_1d);
              //	hpt[bt][iflav][mode] = new TH1D(s.c_str(),s.c_str(), 1, 0.0, 800.0);
            }



          std::string s_trijet_doubletag(Form("%s_%s_%s_trijet_doubletag",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          hpt_trijet_doubletag[bt][iflav][mode] = new TH1D(s_trijet_doubletag.c_str(),s_trijet_doubletag.c_str(), nch, pt_binning_1d);


          std::string sfine(Form("fine_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          hpt_fine[bt][iflav][mode] = new TH1D(sfine.c_str(),sfine.c_str(), 1200, 0.0, 600.0);
 

          std::string seta(Form("eta_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          heta[bt][iflav][mode] = new TH1D(seta.c_str(),seta.c_str(), 60, -3, 3);
          std::string spteta(Form("pteta_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          std::string spteta_eta_bincenter(Form("pteta_eta_bincenter_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
          std::string spteta_pt_bincenter(Form("pteta_pt_bincenter_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));


        




        

          // double number of bins for b only
          // if (iflav == 2) {
          // 	  hpteta[bt][iflav][mode] = new TH2D(spteta.c_str(),spteta.c_str(), 31, 0, 300, 5, 0, 2.6);
          // 	} else {
          // 	  hpteta[bt][iflav][mode] = new TH2D(spteta.c_str(),spteta.c_str(), 15, 0, 300, 5, 0, 2.6);
          // 	}
          // hpteta_pt_bincenter[bt][iflav][mode]->Sumw2();
          //         hpteta_eta_bincenter[bt][iflav][mode]->Sumw2();
        
          if(iflav == 2 || iflav == 4) {
         
            //   if( firstFile.Contains("TTJets") && iflav == 2 ) {
            //                 netabins = netabins_fine;
            //               }
            //               const double *eta_bin = firstFile.Contains("TTJets") && iflav == 2 ? eta_binning_fine : eta_binning;
            const double *eta_bin = eta_binning;
            hpteta[bt][iflav][mode][ich] = new TH2D(spteta.c_str(),spteta.c_str(), 39,pt_binning_b, netabins,eta_bin);
            hpteta_pt_bincenter[bt][iflav][mode] = new TProfile(spteta_pt_bincenter.c_str(),spteta_pt_bincenter.c_str(), 39,pt_binning_b);
            hpteta_eta_bincenter[bt][iflav][mode] = new TProfile(spteta_eta_bincenter.c_str(),spteta_eta_bincenter.c_str(), netabins,eta_bin);
          } else {
          
            hpteta[bt][iflav][mode][ich] = new TH2D(spteta.c_str(),spteta.c_str(), 11, pt_binning, netabins,eta_binning);
            hpteta_pt_bincenter[bt][iflav][mode] = new TProfile(spteta_pt_bincenter.c_str(),spteta_pt_bincenter.c_str(), 11,pt_binning);
            hpteta_eta_bincenter[bt][iflav][mode] = new TProfile(spteta_eta_bincenter.c_str(),spteta_eta_bincenter.c_str(), netabins,eta_binning);
          }


        }
      }
    }
  }
  // create histograms
  hout->mkdir("bc_content");
  hout->cd("bc_content");
  for (int bt=0; bt<nbt; ++bt) {
    for (unsigned int iflav=0; iflav<nflav; ++iflav) {
      for (int mode=0; mode<nmode; ++mode) {

        std::string s_bcontent(Form("bcontent_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
        h_bcontent[bt][iflav][mode] = new TH1D(s_bcontent.c_str(),s_bcontent.c_str(), 10, 0.0, 10.0);

        std::string s_ccontent(Form("ccontent_%s_%s_%s",btlabel[bt].c_str(),flavlabel[iflav].c_str(),modelabel[mode].c_str()));
        h_ccontent[bt][iflav][mode] = new TH1D(s_ccontent.c_str(),s_ccontent.c_str(), 10, 0.0, 10.0);
      }
    }
  }

 

  hout->mkdir("debug");
  hout->cd("debug");

  //  TH1D *h_test_bcontent_b = new TH1D("h_test_bcontent_b","h_test_bcontent_b", 10, 0.0, 10.0);
  //   TH1D *h_test_bcontent_b_gluonsplit = new TH1D("h_test_bcontent_b_gluonsplit","h_test_bcontent_b_gluonsplit", 10, 0.0, 10.0);

  //   TH1D *h_test = new TH1D("h_test", "h_test", 10.0, 0.0, 10.0);
  

  //jet flavour content histograms
  hout->mkdir("JetFlavContent");
  hout->cd("JetFlavContent");
    

  //   const int nflavevent = 2;
  //  const std::string JetFlavEventSelection[nflavevent] = {"incjets","trijets"};

  //  TH1D *h_partonflavorJet[nflavevent];
  //   TH1D *h_hflContentJet[nflavevent];
  //   TH1D *h_iflav[nflavevent];
  //  TH1D *h_theFlavor[nflavevent];
  //  TH1D *h_jetflav_bcontent[nflavevent];
  //   TH1D *h_jetflav_ccontent[nflavevent];
    
  //  const int nbccontent = 2;
  //   const std::string str_bcontent[nbccontent] = {"bcontent","ccontent"};
    
  //   TH1D *h_partonflavorJet_bc_Content[nflavevent][22][nbccontent];
  //   TH1D *h_hflContentJet_bc_Content[nflavevent][6][nbccontent];
  //   TH1D *h_iflav_bc_Content[nflavevent][3][nbccontent];


  //  const int nflavarray = 7;
  //   const std::string flavarray[nflavarray] = {"hard_g","hard_uds", "hard_c", "hard_b", "gsplit_udsg", "gsplit_c", "gsplit_b"};

    
  //   TH1D *h_jet_flav_BC_content_array[nflavarray][nflavevent][nbccontent];
    
  //   for(int iflavarray; iflavarray < nflavevent; iflavarray++)
  //     {
  //       for(int iflavevent; iflavevent < nflavevent; iflavevent++)
  //         {
  //           for(int ibccontent; ibccontent<nbccontent; ibccontent++)
  //             {
  //               h_jet_flav_BC_content_array[iflavarray][iflavevent][ibccontent] = new TH1D( Form("h_jet_flav_BC_content_array_%s_%s_%s", flavarray[iflavarray].c_str(), JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str()),
  //                                                                                           Form("h_jet_flav_BC_content_array_%s_%s_%s", flavarray[iflavarray].c_str(), JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str()),
  //                                                                                           10.0, 0.0, 10.0
  //                                                                                           );
  //             }
  //         }
  //     }


  // for(int iflavevent = 0; iflavevent < nflavevent; iflavevent++)
  //     {

  //       h_jetflav_bcontent[iflavevent] = new TH1D(
  //                                                 Form("h_jetflav_bcontent_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                 Form("h_jetflav_bcontent_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                 10,0.0,10.0);  
        
  //       h_jetflav_ccontent[iflavevent] = new TH1D(
  //                                                 Form("h_jetflav_ccontent_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                 Form("h_jetflav_ccontent_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                 10,0.0,10.0);  
        
  //   h_partonflavorJet[iflavevent] = new TH1D(
  //                                                  Form("h_partonflavorJet_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                  Form("h_partonflavorJet_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                  23, 0.0, 23);  
  //         h_hflContentJet[iflavevent] = new TH1D(
  //                                                Form("h_hflContentJet_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                Form("h_hflContentJet_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                                6, 0.0, 6.0); 
  //         h_iflav[iflavevent] = new TH1D(
  //                                            Form("h_iflav_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                            Form("h_iflav_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                            4, 0.0, 4.0); 

        
  //           h_theFlavor[iflavevent] = new TH1D(
  //                                              Form("h_theFlavor_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                              Form("h_theFlavor_%s", JetFlavEventSelection[iflavevent].c_str()),
  //                                              22, 0.0, 22.0); 



  // for(int ibccontent = 0; ibccontent < nbccontent; ibccontent++)
  //           {
  //             for(int i = 0; i < 22; i++)
  //               {
  //                 h_partonflavorJet_bc_Content[iflavevent][i][ibccontent] = new TH1D(
  //                                                                                    Form("h_partonflavorJet_bc_Content_%s_%s_%d", JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str(),i),
  //                                                                                    Form("h_partonflavorJet_bc_Content_%s_%s_%d", JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str(),i),
  //                                                                                    10.0, 0.0, 10.0);
               


  //                 if(i < 6)
  //                   {
  //                     h_hflContentJet_bc_Content[iflavevent][i][ibccontent] = new TH1D(
  //                                                                                      Form("h_hflContentJet_bc_Content_%s_%s_%d", JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str(),i),
  //                                                                                      Form("h_hflContentJet_bc_Content_%s_%s_%d", JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str(),i),
  //                                                                                      10.0, 0.0, 10.0);
                    
  //                   }

  //                 if(i < 3)
  //                   {

  //                     h_iflav_bc_Content[iflavevent][i][ibccontent] = new TH1D(
  //                                                                                  Form("h_iflav_bc_Content_%s_%s_%d", JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str(),i),
  //                                                                                  Form("h_iflav_bc_Content_%s_%s_%d", JetFlavEventSelection[iflavevent].c_str(),str_bcontent[ibccontent].c_str(),i),
  //                                                                                  10.0, 0.0, 10.0);
  //                   }
  //               }
  //           }
  //   }

  //   // search for the single btag trigger
  //   //string theTrigName("HLT_CentralJet46_CentralJet38_DiBTagIP3D_v1");
  //   //std::string theTrigNameBtag("HLT_L1DoubleJet36Central_BTagIP3D_v1");
  //   std::string theTrigNameBtag("HLT_bbPhi_L1DoubleJet36Central_BTagIP3D_L25MinTag3_v1");
  //   //string theTrigNameBtag("HLT_L1DoubleJet36Central_BTagIP3D_v2");
  //   //std::string theTrigNameRef("HLT_L1DoubleJet36Central");
  //   std::string theTrigNameRef("HLT_bbPhi_L1DoubleJet36Central_v1");



  hout->mkdir("quark_gluon_discrimination");
  hout->cd("quark_gluon_discrimination");
 

  const unsigned int nflav_qg = 2;
  const std::string flavlabel_qg[nflav_qg] = {"slim", "fat"};
  TH1D *h_quark_gluon_discrimination[nflav_qg];
  const unsigned int npt_quark_gluon_discrimination = 15;
  //const double ptbinning_quark_gluon_discrimination[npt_quark_gluon_discrimination+1] = {
  //   2.00000e+01, 
  //       2.55761e+01, 
  //       3.27068e+01, 
  //       4.18256e+01, 
  //       5.34867e+01, 
  //       6.83990e+01, 
  //       8.74690e+01, 
  //       1.11856e+02, 
  //       1.43041e+02, 
  //       1.82922e+02, 
  //       2.33921e+02, 
  //       2.99140e+02, 
  //       3.82541e+02, 
  //       4.89195e+02, 
  //       6.25585e+02, 
  //    800.0
  //  };
  TH1D *h_quark_gluon_discrimination_ptbin[nflav_qg][npt_quark_gluon_discrimination];


 

  for(unsigned int iflav = 0; iflav < nflav_qg; iflav++) {
    std::string s(Form("h_quark_gluon_discrimination_%s", flavlabel_qg[iflav].c_str()));
    h_quark_gluon_discrimination[iflav] = new TH1D(s.c_str(), s.c_str(), 30, 0.0, 1.0);
   
    for(unsigned int ipt = 0; ipt < npt_quark_gluon_discrimination; ipt++) {
      std::string sm(Form("h_quark_gluon_discrimination_ptbin_%d_%s", ipt,flavlabel_qg[iflav].c_str()));
      h_quark_gluon_discrimination_ptbin[iflav][ipt] = new TH1D(sm.c_str(), sm.c_str(), 30, 0.0, 1.0);
    }
  }




  //test
  // search for the single btag trigger
 
  std::string theTrigNameBtag("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v");
  std::string theTrigNameRef("HLT_L1DoubleJet36Central_v");

  std::vector<std::string> triggerFilterList;
  triggerFilterList.push_back("HLT_Jet60Eta1p7_Jet53Eta1p7_DiBTagIP3DFastPV_v");

  std::vector<std::string> HLTBtagTriggerObjectFilterList;
  HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
  //   HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiL1FastJetFastPV");
  //     HLTBtagTriggerObjectFilterList.push_back("hltBLifetimeL3FilterbbPhiLooseL1FastJetFastPV");

  checkHLTBtagMatch HLTBtagMatchObj(
                                    triggerFilterList, 
                                    HLTBtagTriggerObjectFilterList, 
                                    HLTBtagTriggerObjectList,
                                    false,
                                    NULL
                                    );

  //end of test








  std::vector<std::string>::iterator tSlotBtag = std::find(genericTriggerList.begin(), genericTriggerList.end(), theTrigNameBtag);
  if (tSlotBtag != genericTriggerList.end()) {
    std::cout << "Btag trigger found at slot " << tSlotBtag - genericTriggerList.begin() << std::endl;
  } else {
    std::cout << "Btag trigger not found in any slot" << std::endl;
    // return;

    
  }
  unsigned int tNumberBtag = tSlotBtag - genericTriggerList.begin();
 

  

  // loop over the tree
  for (Int_t i=0; i<nentries; i++) {
    hbbtree->GetEntry(i);
    if(i % 200000 == 0)
      {
        const double percent = 100.0*(double)i/((double)nentries);
        std::cout << "process event " << i << " - " << percent << " % done." << std::endl;
      }
   

    //test relative online sf
    //    RelativeOnlineSFandUncertainties sf;
    
    const double pileup_weight = PileUpWeightObj.GetWeight(nPUTruth);
    
    const double weight = pileup_weight * lumiScaleFac;


    jpm.Fill(weight);
    // loop over the jets

    

    //determine flavour fractions
    //medium mass
    const double jetPtMin_mediummass[3] = { 60.0, 53.0, 20.};
    const double jetPtMin_highmass[3] = {80.0, 70.0, 20.0};
    const double jetPtMin_veryhighmass[3] = {160.0, 120.0, 20.0};
    const std::vector<int> leadingjets_mediummass = trijetselection(jetPtMin_mediummass, 1.7, 2.4);
    const std::vector<int> leadingjets_highmass = trijetselection(jetPtMin_highmass, 1.7, 2.4);
    const std::vector<int> leadingjets_veryhighmass = trijetselection(jetPtMin_veryhighmass, 2.4, 9999999.0);
    
    jetfrac_mediummass.DetermineFlavFractions(leadingjets_mediummass,weight);
    jetfrac_highmass.DetermineFlavFractions(leadingjets_highmass,weight);
    jetfrac_veryhighmass.DetermineFlavFractions(leadingjets_veryhighmass,weight);
    //end of determination of flavour fractions
    



    //trijet selection
    bool trijet_event = false;
    const double jetPtMin[3] = { 60.0, 53.0, 20};  

    std::vector<int> leadingJets;
    // find set of leading jets
    int nJet = 0;
    // loop over the jets
    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      if(!puJetIDLoose[iJet])
        continue;
      if(!jetIDLoose[iJet] ) continue;
      if (nJet >= 3) break;
      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;
      if ( (ptJet[iJet] > jetPtMin[nJet]) && (ptJet[iJet] < 3500.0) ) {
       
        leadingJets.push_back(iJet);
        ++nJet;
      }
    }
    
    float deltaRj = -1;
    
    if(nJet>=2)
      {
        float dphij =  phiJet[leadingJets[1]] - phiJet[leadingJets[0]];
        if (dphij>3.1415926) dphij -= 2*3.1415926;
        if (dphij<-3.1415926) dphij += 2*3.1415926;
      
        float detaj = etaJet[leadingJets[1]] - etaJet[leadingJets[0]];
        deltaRj = sqrt( dphij*dphij + detaj*detaj );
      
      
        

      }
    if(nJet>=3 && deltaRj > 1.0)
      trijet_event = true;

    for (int iJet=0; iJet<numberOfJets; ++iJet) {
      
      if(!puJetIDLoose[iJet])
        continue;
      if(!jetIDLoose[iJet] ) continue;
      //if ( ptJet[iJet] < 70.0 || ptJet[iJet] > 80. ) continue;
      if ( ptJet[iJet] < 20. ) continue;
      //   if(! testJetIsolation(
      //                            phiJet[iJet],
      //                            etaJet[iJet],
      //                            numberOfJets,
      //                            phiJet,
      //                            etaJet,
      //                            ptJet    
      //                            )) continue;
      //      if (! (fabs(etaJet[iJet])<2.6) ) continue;
      if (! (fabs(etaJet[iJet])<2.4) ) continue; //2.6 eta cut is DEFAULT! but only the very high mass trigger goes to 2.4, the others to 1.7
      //if (! (fabs(etaJet[iJet])<1.5) ) continue;

      if (! (numberOfConstituentsInJet[iJet] > 1) ) continue;

      if (! (isJetMatchedPartonExist[iJet]) ) continue;

      //       // HERE: explicit pt cut!
      //       if (ptJet[iJet]<80) continue;
      //       if ( (numberOfChargedConstituentsInJet[iJet] < 8)
      // 	   || (numberOfChargedConstituentsInJet[iJet] >9) ) continue;
      



      //  const int bcarray[nbccontent] = {(int)BContentJet3[iJet],(int)CContentJet3[iJet]};

      // for(int iflavevent; iflavevent < nflavevent; iflavevent++)
      //         {
      //           if(JetFlavEventSelection[iflavevent] == "trijets" && !trijet_event)
      //             continue;

      //           for(int ibccontent; ibccontent<nbccontent; ibccontent++)
      //             {
      //               if(
      //                  abs(partonFlavorJet[iJet]) == 1 || 
      //                  abs(partonFlavorJet[iJet]) == 2 ||
      //                  abs(partonFlavorJet[iJet]) == 3
      //                  )
      //                 {
      //                   //"hard" uds
      //                   h_jet_flav_BC_content_array[0][iflavevent][ibccontent]->Fill(bcarray[ibccontent],weight);
      //                 }

      //               if(
      //                  abs(partonFlavorJet[iJet]) == 4
      //                  )
      //                 {
      //                   //"hard" c
      //                   h_jet_flav_BC_content_array[1][iflavevent][ibccontent]->Fill(bcarray[ibccontent],weight);
      //                 }

      //               if(abs(partonFlavorJet[iJet]) == 5)
      //                 {
      //                   //"hard" b
      //                   h_jet_flav_BC_content_array[2][iflavevent][ibccontent]->Fill(bcarray[ibccontent],weight);
      //                 }
          
      //               if(abs(partonFlavorJet[iJet]) == 21)
      //                 {
      //                   if(abs(hflContentJet[iJet]) == 0)
      //                     {
      //                       //"g split" udsg
      //                       h_jet_flav_BC_content_array[3][iflavevent][ibccontent]->Fill(bcarray[ibccontent],weight);
      //                     }
      //                   if(abs(hflContentJet[iJet]) == 4)
      //                     {
      //                       //"g split" c
      //                       h_jet_flav_BC_content_array[4][iflavevent][ibccontent]->Fill(bcarray[ibccontent],weight);
      //                     }
      //                   if(abs(hflContentJet[iJet]) == 5)
      //                     {
      //                       //"g split" b
      //                       h_jet_flav_BC_content_array[5][iflavevent][ibccontent]->Fill(bcarray[ibccontent],weight);
      //                     }
               
              
      //                 }
          

      //             }
      //         }



















      //0 -> old
      //1 -> old plus multi quarks separately // not sure whether this makes sense ...
      //2 -> new (cone 0.X)
      //3 -> new (cone 0.X) plus multi quarks separately

      //const int parton_matching_mode = 3;
      
      JetFlavor::Code iflav = JetFlavor::code(iJet);
      

      //    if(parton_matching_mode == 0 || parton_matching_mode == 1)
      //         {
      //           //switch (abs(flavorJetMatchedParton[iJet])) {

      //           // flavor fix
      //           int theFlavor = partonFlavorJet[iJet]; 
      //           if (hflContentJet[iJet] != 0) theFlavor = hflContentJet[iJet];

      //           switch (abs(theFlavor)) {
      //           case 1:
      //           case 2:
      //           case 3:
      //             iflav = 0;
      //             break;
      //           case 4:
      //             iflav = 1;
      //             break;
      //           case 5:
      //             iflav = 2;
      //             break;
      //           case 21:
      //             iflav = 0;
      //             break;
      //           }
      //           if( abs(partonFlavorJet[iJet]) == 4 && (int)CContentJet3[iJet] > 0)
      //             h_test->Fill(0.0 ,weight);
      //           else if( abs(partonFlavorJet[iJet]) == 4 && (int)CContentJet3[iJet] == 0)
      //             h_test->Fill(1.0 ,weight);//inconsistent!!!
          
      //           if(parton_matching_mode == 1)
      //             {
      //               if((int)BContentJet3[iJet] == 1)
      //                 iflav = 2;
      //               else if((int)BContentJet3[iJet] > 1)
      //                 iflav = 4;
      //               else if((int)CContentJet3[iJet] == 1)
      //                 iflav = 1;
      //               else if((int)CContentJet3[iJet] > 1)
      //                 iflav = 3;
      //             }
      //         }
      //       else if(parton_matching_mode == 2 || parton_matching_mode == 3)
      //         {
      //           if( (int)BContentJet3[iJet] > 0 )
      //             {
      //               if(parton_matching_mode == 3)
      //                 {
      //                   if((int)BContentJet3[iJet] == 1)
      //                     iflav = 2;
      //                   else 
      //                     iflav = 4;  
      //                 }
      //               else
      //                 iflav = 2;
      //             }
      //           else if( (int)CContentJet3[iJet] > 0 )
      //             {
      //               if(parton_matching_mode == 3)
      //                 {
      //                   if( (int)CContentJet3[iJet] == 1 )
      //                     iflav = 1;
      //                   else
      //                     iflav = 3; 
      //                 }
      //               else
      //                 iflav = 1;
      //             }
      //           else
      //             {
      //               //everything else is assumed to be udsg ... 
      //               switch (abs(partonFlavorJet[iJet])) {
      //               case 1:
      //               case 2:
      //               case 3:
      //               case 4:
      //               case 5:
      //               case 21:
      //                 iflav = 0;
      //                 break;
      //               }
      //             }
      //         }
      //       else
      //         {
      //           std::cout << "wrong parton matching mode?" << std::endl;
      //           return;
      //         }



      //now p_t plotting of different flavor definitions for the efficiencies

      JetFlavor::Codes vec_iflav;
      
      //      vec_iflav.push_back(iflav);//0
     
      //         if(isJetMatchedPartonExist[iJet])
      //           {
      //             int iflav_tmp = -1;
      //             switch (abs(flavorJetMatchedParton[iJet])) {
      //             case 1:
      //             case 2:
      //             case 3:
      //               iflav_tmp = 0;
      //               break;
      //             case 4:
      //               iflav_tmp = 1;
      //               break;
      //             case 5:
      //               iflav_tmp = 2;
      //               break;
      //             case 21:
      //               iflav_tmp = 0;
      //               break;
      //             }
      //             vec_iflav.push_back(iflav_tmp);//1
      //           }
      //         else
      //           {
      //             vec_iflav.push_back(-1);//1
      //           }
      
      //         const int b_array_def[6] = {(int)BContentJet3[iJet],(int)BContentJet4[iJet],(int)BContentJet5[iJet],
      //                                     (int)BContentJet6[iJet],(int)BContentJet7[iJet],(int)BContentJet10[iJet]};
      //         const int c_array_def[6] = {(int)CContentJet3[iJet],(int)CContentJet4[iJet],(int)CContentJet5[iJet],
      //                                     (int)CContentJet6[iJet],(int)CContentJet7[iJet],(int)CContentJet10[iJet]};
      //         for(int ibcarray = 0; ibcarray < 6; ibcarray++)
      //           {
      //             if( b_array_def[ibcarray] > 0 )
      //               {
      //                 vec_iflav.push_back(2);
      //               }
      //             else if(c_array_def[ibcarray] > 0 )
      //               {
      //                 vec_iflav.push_back(1);
      //               }
      //             else
      //               {
      //                 int iflav_tmp = -1;
      //                 //everything else is assumed to be udsg ... 
      //                 switch (abs(partonFlavorJet[iJet])) {
      //                 case 1:
      //                 case 2:
      //                 case 3:
      //                 case 4:
      //                 case 5:
      //                 case 21:
      //                   iflav_tmp = 0;
      //                   break;
      //                 }
      //                 vec_iflav.push_back(iflav_tmp);
      //               }
      //           }


      //         //    
      //         for(int ibcarray = 0; ibcarray < 6; ibcarray++)
      //           {
      //             if( b_array_def[ibcarray] > 0 )
      //               {
      //                 if(b_array_def[ibcarray] == 1)
      //                   vec_iflav.push_back(2);
      //                 else 
      //                   vec_iflav.push_back(4);  
      //               }
      //             else if(c_array_def[ibcarray] > 0 )
      //               {
      //                 if(c_array_def[ibcarray] == 1 )
      //                   vec_iflav.push_back(1);
      //                 else
      //                   vec_iflav.push_back(3); 
      //               }
      //             else
      //               {
      //                 int iflav_tmp = -1;
      //                 //everything else is assumed to be udsg ... 
      //                 switch (abs(partonFlavorJet[iJet])) {
      //                 case 1:
      //                 case 2:
      //                 case 3:
      //                 case 4:
      //                 case 5:
      //                 case 21:
      //                   iflav_tmp = 0;
      //                   break;
      //                 }
      //                 vec_iflav.push_back(iflav_tmp);
      //               }
      //           }
      //         //end of multi b,c quarks in cone



      









      //         if(vec_iflav.size() != (int)(n_flav_def*0.5))
      //           {
      //             std::cout << "error: n_flav_def != vec_iflav.size(): " << n_flav_def << " " << vec_iflav.size() << std::endl;
      //             //          return;
      //           }


      //fill here hpt histo
      for(size_t i_flav_def = 0; i_flav_def < vec_iflav.size(); i_flav_def++) 
        {
          if(vec_iflav.at(i_flav_def) != JetFlavor::UNDEFINED)
            continue;
          
          for (int bt=0; bt<nbt; ++bt) {
            hpt[bt][vec_iflav.at(i_flav_def)][1][i_flav_def]->Fill(ptJet[iJet],weight);
        
            hpt[bt][vec_iflav.at(i_flav_def)][1][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);

            if (btdiscr[bt] == "TCHE") {
              if (tcHEBJetTag[iJet]>btcut[bt]) {
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def]->Fill(ptJet[iJet],weight);
              
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);
              }
            } else if (btdiscr[bt] == "TCHP") {
              if (tcHPBJetTag[iJet]>btcut[bt]) {
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def]->Fill(ptJet[iJet],weight);
               
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);

              }
            } else if (btdiscr[bt] == "nTCHE") {
              if (ntcHEBJetTag[iJet]>btcut[bt]) {
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def]->Fill(ptJet[iJet],weight);
              
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);

              }
            }
            else if (btdiscr[bt] == "nTCHP") {
              if (ntcHPBJetTag[iJet]>btcut[bt]) {
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def]->Fill(ptJet[iJet],weight);

              
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);
              }
            } else if (btdiscr[bt] == "CSV") {
              if (combSVBJetTag[iJet]>btcut[bt]) {
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def]->Fill(ptJet[iJet],weight);

                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);


              }
            } else if (btdiscr[bt] == "SSVHP") {
              if (svHPBJetTag[iJet]>btcut[bt]) {
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def]->Fill(ptJet[iJet],weight);

     
                hpt[bt][vec_iflav.at(i_flav_def)][0][i_flav_def+vec_iflav.size()]->Fill(ptJet[iJet],weight);


              }
            }
          }
        }
      //end of fill hpt histo

      //end of p_t plotting of different flavor definitions for the efficiencies



      //       if(abs(partonFlavorJet[iJet]) == 21)
      //         {
      //           if(iflav == 2)
      //             {
      //               h_test_bcontent_b_gluonsplit->Fill(bcarray[0],weight);
      //             }
      //         }
      //       if(abs(partonFlavorJet[iJet]) == 5)
      //         {
        
      //           {
      //             h_test_bcontent_b->Fill(bcarray[0],weight);
      //           }
      //         }









      




      // //       if (iflav == -1)
      // //         {
      // //           continue;
      // //         }
    
      
      //       //       h_partonflavorJet[0]->Fill(abs(partonFlavorJet[iJet]),weight);
      //       //       h_hflContentJet[0]->Fill(abs(hflContentJet[iJet]),weight);
      //       //       h_iflav[0]->Fill(iflav,weight);
      //       //       h_theFlavor[0]->Fill(abs(theFlavor),weight);
      
      //       h_jetflav_bcontent[0]->Fill((int)BContentJet3[iJet],weight);
      //       h_jetflav_ccontent[0]->Fill((int)CContentJet3[iJet],weight);
      
      //       //       for(int ibccontent = 0; ibccontent < nbccontent; ibccontent++)
      //       //         {
      //       //           if(abs(partonFlavorJet[iJet]) < 22)
      //       //             h_partonflavorJet_bc_Content[0][abs(partonFlavorJet[iJet])][ibccontent]->Fill(bcarray[ibccontent],weight);
          
      //       //           if(iflav < 3)
      //       //             h_iflav_bc_Content[0][iflav][ibccontent]->Fill(bcarray[ibccontent],weight);
          
      //       //           if(abs(hflContentJet[iJet]) < 6)
      //       //             h_hflContentJet_bc_Content[0][abs(hflContentJet[iJet])][ibccontent]->Fill(bcarray[ibccontent],weight);
      //       //         }
      
      //       if(trijet_event)
      //         {
      //           //    h_partonflavorJet[1]->Fill(abs(partonFlavorJet[iJet]),weight);
      //           //           h_hflContentJet[1]->Fill(abs(hflContentJet[iJet]),weight);
      //           //           h_iflav[1]->Fill(iflav,weight);
      //           //           h_theFlavor[1]->Fill(abs(theFlavor),weight);

      //           h_jetflav_bcontent[1]->Fill((int)BContentJet3[iJet],weight);
      //           h_jetflav_ccontent[1]->Fill((int)CContentJet3[iJet],weight);
     

      //           //           for(int ibccontent = 0; ibccontent < nbccontent; ibccontent++)
      //           //             {
      //           //               if(abs(partonFlavorJet[iJet]) < 22)
      //           //                 h_partonflavorJet_bc_Content[1][abs(partonFlavorJet[iJet])][ibccontent]->Fill(bcarray[ibccontent],weight);
          
      //           //               if(iflav < 3)
      //           //                 h_iflav_bc_Content[1][iflav][ibccontent]->Fill(bcarray[ibccontent],weight);
          
      //           //               if(abs(hflContentJet[iJet]) < 6)
      //           //                 h_hflContentJet_bc_Content[1][abs(hflContentJet[iJet])][ibccontent]->Fill(bcarray[ibccontent],weight);
      //           //             }
   
      //         }
     
      
      
      //    if(iflav != -1 && ptJet[iJet] > 0.0 && ptJet[iJet] < 800.0) {
      //         const double discr = qgdiscriminant.GetDiscriminant(iJet);
      //         const unsigned int iflav_qg = qgdiscriminant.GetSlimFat(iJet, flavlabel, iflav);

      //         h_quark_gluon_discrimination[iflav_qg]->Fill(discr > 0.0 ? discr : -1.0, weight);
        
      //         for(unsigned int ibin =0; ibin < npt_quark_gluon_discrimination; ibin++) {
      //           if(ptJet[iJet] >= ptbinning_quark_gluon_discrimination[ibin] && ptJet[iJet] < ptbinning_quark_gluon_discrimination[ibin+1]) {
      //             h_quark_gluon_discrimination_ptbin[iflav_qg][ibin]->Fill(discr > 0.0 ? discr : -1.0,weight);
      //             break;
      //           }
      //         }
      //       }


      

      for (int bt=0; bt<nbt; ++bt)  {
        if (iflav == JetFlavor::UNDEFINED)
          {
            continue;
          }












        const int ich = 0;
        //        int ich = -1;

        //  for(int ichbin = 0; ichbin < nchargedconstituents; ichbin++) {
        //           if(
        //              numberOfChargedConstituentsInJet[iJet] >= nconstjet[ichbin]
        //              && numberOfChargedConstituentsInJet[iJet] < nconstjet[ichbin+1]
        //              ) {
        //             ich = ichbin;
        //           }
        //         }
        //         if(iflav == 0) {
        //           ich = 0;
        //         }
        //         if(ich == -1) {
        //           continue;
        //         }
        //         //numberOfChargedConstituentsInJet
        

        //         //merge fat and slim jets (not used for default efficiency functions)
        //         if(iflav == 4) {
        //           iflav = 2;
        //         }
        //         if(iflav == 3) {
        //           iflav = 1;
        //         }
        
        







        if (btdiscr[bt] == "CSV") {
          h_csv_discr->Fill( combSVBJetTag[iJet], weight);
          if(testJetIsolation(
                              phiJet[iJet],
                              etaJet[iJet],
                              numberOfJets,
                              phiJet,
                              etaJet,
                              ptJet    
                              )) {
            h_csv_discr_jetisolation->Fill( combSVBJetTag[iJet], weight);
          }
        }
            
        
        // if ((btdiscr[bt] == "CSV" && combSVBJetTag[iJet] > 0.0 && combSVBJetTag[iJet] < 1.0) || btdiscr[bt] != "CSV")
        
      
        jetconsti.Fill(bt,iflav,1,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);

       



          
        hpteta[bt][iflav][1][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);


        hpteta_eta_bincenter[bt][iflav][1]->Fill(fabs(etaJet[iJet]),fabs(etaJet[iJet]),weight);
        hpteta_pt_bincenter[bt][iflav][1]->Fill(ptJet[iJet],ptJet[iJet],weight);

        hpteta_eta_bincenter[bt][iflav][0]->Fill(fabs(etaJet[iJet]),fabs(etaJet[iJet]),weight);
        hpteta_pt_bincenter[bt][iflav][0]->Fill(ptJet[iJet],ptJet[iJet],weight);
       

        hpt_fine[bt][iflav][1]->Fill(ptJet[iJet],weight);
        
        h_bcontent[bt][iflav][1]->Fill((int)BContentJet3[iJet],weight);
        h_ccontent[bt][iflav][1]->Fill((int)CContentJet3[iJet],weight);

        heta[bt][iflav][1]->Fill(etaJet[iJet],weight);

        if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
          h_all[bt][iflav][1]->Fill(0.5, weight);
        }
        
        h_nvtx[bt][iflav][1]->Fill((double)nPV,weight);
        h_nPU[bt][iflav][1]->Fill((double)nPU,weight);
        
        if(puJetIDLoose[iJet])
          {
            h_nvtx_JetIDLoose[bt][iflav][1]->Fill((double)nPV,weight);
            h_nPU_JetIDLoose[bt][iflav][1]->Fill((double)nPU,weight); 
          }

        if(trijet_event)
          {
            h_nvtx_trijet[bt][iflav][1]->Fill((double)nPV,weight);
            h_nPU_trijet[bt][iflav][1]->Fill((double)nPU,weight);
          }

        if(trijet_event && (combSVBJetTag[leadingJets[0]]>btcut[bt]) && leadingJets[0] != iJet)
          {
            hpt_trijet_doubletag[bt][iflav][1]->Fill(ptJet[iJet],weight);  
          }
        
        

        

        if ( online[bt] && (! (trgAccept & (1<<tNumberBtag)))) continue;
        if ( btagMatch[bt] && ! HLTBtagMatchObj.check(isJetWithHltBtagBitPattern[iJet],isJetWithL25JetBitPattern[iJet]) ) continue;





        if (btdiscr[bt] == "TCHE") {
          h_btagdiscriminator[bt][iflav][0]->Fill(tcHEBJetTag[iJet],weight);
          if (tcHEBJetTag[iJet]>btcut[bt]) {
          
            if(trijet_event)
              {
                h_nvtx_trijet[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_trijet[bt][iflav][0]->Fill((double)nPU,weight);

              }
            h_nvtx[bt][iflav][0]->Fill((double)nPV,weight);
            h_nPU[bt][iflav][0]->Fill((double)nPU,weight);

            if(puJetIDLoose[iJet])
              {
                h_nvtx_JetIDLoose[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_JetIDLoose[bt][iflav][0]->Fill((double)nPU,weight); 
              }
            if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
              h_all[bt][iflav][0]->Fill(0.5, weight);
            }
          
            hpt_fine[bt][iflav][0]->Fill(ptJet[iJet],weight);
                
            h_bcontent[bt][iflav][0]->Fill((int)BContentJet3[iJet],weight);
            h_ccontent[bt][iflav][0]->Fill((int)CContentJet3[iJet],weight);

            jetconsti.Fill(bt,iflav,0,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);

            heta[bt][iflav][0]->Fill(etaJet[iJet],weight);
            hpteta[bt][iflav][0][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);
          }
        
        } else if (btdiscr[bt] == "TCHP") {
          h_btagdiscriminator[bt][iflav][0]->Fill(tcHPBJetTag[iJet],weight);
          if (tcHPBJetTag[iJet]>btcut[bt]) {
          
            if(trijet_event)
              {
                h_nvtx_trijet[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_trijet[bt][iflav][0]->Fill((double)nPU,weight);

              }

            if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
              h_all[bt][iflav][0]->Fill(0.5, weight);
            }


            h_nvtx[bt][iflav][0]->Fill((double)nPV,weight);
            h_nPU[bt][iflav][0]->Fill((double)nPU,weight);

            if(puJetIDLoose[iJet])
              {
                h_nvtx_JetIDLoose[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_JetIDLoose[bt][iflav][0]->Fill((double)nPU,weight); 
              }

               
            hpt_fine[bt][iflav][0]->Fill(ptJet[iJet],weight);
            
            h_bcontent[bt][iflav][0]->Fill((int)BContentJet3[iJet],weight);
            h_ccontent[bt][iflav][0]->Fill((int)CContentJet3[iJet],weight);

            jetconsti.Fill(bt,iflav,0,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);
            heta[bt][iflav][0]->Fill(etaJet[iJet],weight);
            hpteta[bt][iflav][0][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);
          }
      
        } else if (btdiscr[bt] == "nTCHE") {
          h_btagdiscriminator[bt][iflav][0]->Fill(ntcHEBJetTag[iJet],weight);
          if (ntcHEBJetTag[iJet]>btcut[bt]) {
          
            if(trijet_event)
              {
                h_nvtx_trijet[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_trijet[bt][iflav][0]->Fill((double)nPU,weight);

              }

            if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
              h_all[bt][iflav][0]->Fill(0.5, weight);
            }


            h_nvtx[bt][iflav][0]->Fill((double)nPV,weight);
            h_nPU[bt][iflav][0]->Fill((double)nPU,weight);

            if(puJetIDLoose[iJet])
              {
                h_nvtx_JetIDLoose[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_JetIDLoose[bt][iflav][0]->Fill((double)nPU,weight); 
              }

            
            hpt_fine[bt][iflav][0]->Fill(ptJet[iJet],weight);
            
            h_bcontent[bt][iflav][0]->Fill((int)BContentJet3[iJet],weight);
            h_ccontent[bt][iflav][0]->Fill((int)CContentJet3[iJet],weight);
            jetconsti.Fill(bt,iflav,0,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);

            heta[bt][iflav][0]->Fill(etaJet[iJet],weight);
            hpteta[bt][iflav][0][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);
          }
    
        } else if (btdiscr[bt] == "nTCHP") {
          h_btagdiscriminator[bt][iflav][0]->Fill(ntcHPBJetTag[iJet],weight);
          if (ntcHPBJetTag[iJet]>btcut[bt]) {
        
            if(trijet_event)
              {
                h_nvtx_trijet[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_trijet[bt][iflav][0]->Fill((double)nPU,weight);

              }

            if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
              h_all[bt][iflav][0]->Fill(0.5, weight);
            }

            h_nvtx[bt][iflav][0]->Fill((double)nPV,weight);
            h_nPU[bt][iflav][0]->Fill((double)nPU,weight);

            if(puJetIDLoose[iJet])
              {
                h_nvtx_JetIDLoose[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_JetIDLoose[bt][iflav][0]->Fill((double)nPU,weight); 
              }
  
          
            hpt_fine[bt][iflav][0]->Fill(ptJet[iJet],weight);
            
            h_bcontent[bt][iflav][0]->Fill((int)BContentJet3[iJet],weight);
            h_ccontent[bt][iflav][0]->Fill((int)CContentJet3[iJet],weight);
            jetconsti.Fill(bt,iflav,0,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);

            heta[bt][iflav][0]->Fill(etaJet[iJet],weight);
            hpteta[bt][iflav][0][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);
          }
  
        } else if (btdiscr[bt] == "CSV") {
          h_btagdiscriminator[bt][iflav][0]->Fill(combSVBJetTag[iJet],weight);
          if (combSVBJetTag[iJet]>btcut[bt]) {
          
            if(trijet_event)
              {
                h_nvtx_trijet[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_trijet[bt][iflav][0]->Fill((double)nPU,weight);
              }
            if(trijet_event && (combSVBJetTag[leadingJets[0]]>btcut[bt]) && leadingJets[0] != iJet)
              {
                hpt_trijet_doubletag[bt][iflav][0]->Fill(ptJet[iJet],weight);  
              }

            if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
              h_all[bt][iflav][0]->Fill(0.5, weight);
            }

            jetconsti.Fill(bt,iflav,0,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);
            h_nvtx[bt][iflav][0]->Fill((double)nPV,weight);
            h_nPU[bt][iflav][0]->Fill((double)nPU,weight);

            if(puJetIDLoose[iJet])
              {
                h_nvtx_JetIDLoose[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_JetIDLoose[bt][iflav][0]->Fill((double)nPU,weight); 
              }

            //     if(trijet_event && (combSVBJetTag[leadingJets[0]]>btcut[bt]) && leadingJets[0] != iJet)
            {
             
              hpt_fine[bt][iflav][0]->Fill(ptJet[iJet],weight);

              h_bcontent[bt][iflav][0]->Fill((int)BContentJet3[iJet],weight);
              h_ccontent[bt][iflav][0]->Fill((int)CContentJet3[iJet],weight);


              heta[bt][iflav][0]->Fill(etaJet[iJet],weight);
              hpteta[bt][iflav][0][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);
            }
          }
 
        } else if (btdiscr[bt] == "SSVHP") {
          h_btagdiscriminator[bt][iflav][0]->Fill(svHPBJetTag[iJet],weight);
          if (svHPBJetTag[iJet]>btcut[bt]) {
        
            if(trijet_event)
              {
                h_nvtx_trijet[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_trijet[bt][iflav][0]->Fill((double)nPU,weight);

              }

            if(ptJet[iJet] > 30.0 && ptJet[iJet] < 670.0) {
              h_all[bt][iflav][0]->Fill(0.5, weight);
            }

            jetconsti.Fill(bt,iflav,0,ptJet[iJet],numberOfChargedConstituentsInJet[iJet],numberOfConstituentsInJet[iJet],weight, iJet);
            h_nvtx[bt][iflav][0]->Fill((double)nPV,weight);
            h_nPU[bt][iflav][0]->Fill((double)nPU,weight);

            if(puJetIDLoose[iJet])
              {
                h_nvtx_JetIDLoose[bt][iflav][0]->Fill((double)nPV,weight);
                h_nPU_JetIDLoose[bt][iflav][0]->Fill((double)nPU,weight); 
              }

            
            hpt_fine[bt][iflav][0]->Fill(ptJet[iJet],weight);

            h_bcontent[bt][iflav][0]->Fill((int)BContentJet3[iJet],weight);
            h_ccontent[bt][iflav][0]->Fill((int)CContentJet3[iJet],weight);


            heta[bt][iflav][0]->Fill(etaJet[iJet],weight);
            hpteta[bt][iflav][0][ich]->Fill(ptJet[iJet],fabs(etaJet[iJet]),weight);
          }
        
        } else {
          std::cout << "error: Invalid TC discriminant " << btdiscr[bt] << std::endl;
          throw std::exception();
          return;
        }
    
      }
    }
  
  }

  // termination

  // for (int bt=0; bt<nbt; ++bt) {
  //     for (int iflav=0; iflav<nflav; ++iflav) {
  //       std::string sn = hpt[bt][iflav][2]->GetName();
  //       std::string st = hpt[bt][iflav][2]->GetTitle();
  //       //      delete hpt[bt][iflav][2];
  //       // hpt[bt][iflav][2] = (TH1D*) hpt[bt][iflav][0]->Clone();
  //       //       hpt[bt][iflav][2]->SetName(sn.c_str());
  //       //       hpt[bt][iflav][2]->SetTitle(st.c_str());
  //       //       hpt[bt][iflav][2]->Divide( hpt[bt][iflav][1] );
  //     }
  //   }

  //  jetconsti.Normalise();

  hout->Write();
  hout->Close();
}
                                                                             

int main(int narg,char** varg) {   
  const std::string option = varg[1];
               
  int pileupwtx = 0;
  if(option == "0")
    pileupwtx = 0;
  else if(option == "1")
    pileupwtx = 1;
  else if(option == "2")
    pileupwtx = 2;
  else if(option == "3")
    pileupwtx = 3;
  else {
    std::cout << "error: unknown PU option : " << option << std::endl;
    throw std::exception();
  }
                                                   
  readmistag(pileupwtx);
}
