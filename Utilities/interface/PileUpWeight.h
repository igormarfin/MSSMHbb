#ifndef _PileUpWeight_
#define _PileUpWeight_

#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <exception>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TH1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TChain.h>
#include <TFileCollection.h>
#include <TGraph2D.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TProfile.h>
#include <math.h>
#include <TMath.h>
#include <TFileInfo.h>

class PileUpWeight {
public:
  PileUpWeight(const bool mc,const int mode, TFile *hout) 
    : _doMC(mc), pileupwtx(mode)
  {
    std::cout << "Pileup weights mode: " << pileupwtx << " " << _doMC << std::endl;

    PU_rootfile_MC = "/data/user/jbehr/PU_reweighting/MC_Summer12_PU_S10-600bins.root";
   
    PU_scenario.push_back("PU_mediummass");
    PU_scenario.push_back("PU_highmass");
    PU_scenario.push_back("PU_veryhighmass");
    
    //now control trigger
    PU_scenario.push_back("PU_dijet40");
    PU_scenario.push_back("PU_dijet80");
    PU_scenario.push_back("PU_L1DoubleJet36C");
    
    

    hout->mkdir("pileup_weights");
    hout->cd("pileup_weights");

//     PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_Jet60Eta1p7_Jet53Eta1p7_600bins_true.root");
//     PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_Jet80Eta1p7_Jet70Eta1p7_600bins_true.root");
//     PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_Jet160Eta2p4_Jet120Eta2p4_600bins_true.root");
    
    PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_Jet60Eta1p7_Jet53Eta1p7_600bins_true_190389-208686.root");
    PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_Jet80Eta1p7_Jet70Eta1p7_600bins_true_190389-208686.root");
    PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_Jet160Eta2p4_Jet120Eta2p4_600bins_true_190389-208686.root");
    PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_DiJet40Eta2p6_BTagIP3DFastPV_600bins_true_190389-208686.root");
    PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_DiJet80Eta2p6_BTagIP3DFastPVLoose_600bins_true_190389-208686.root");
    PU_rootfile_Data.push_back("/data/user/jbehr/PU_reweighting/MyDataPileupHistogram_L1DoubleJet36Central_600bins_true_190389-208686.root");
   


    if(pileupwtx < 0 || (unsigned int)pileupwtx > PU_rootfile_Data.size() ||  PU_rootfile_Data.size() != PU_scenario.size()) {
      std::cout << "error: misconfiguration of PU reweighting object" << std::endl;
      throw std::exception();
    }


  
    TH1::AddDirectory(kFALSE);
    TFile *file_mc = new TFile(PU_rootfile_MC.c_str());
    if(file_mc) {
    pu_histo_mc = (TH1D*)(file_mc->Get("pileup"))->Clone();
 if(!pu_histo_mc) {
          std::cout << "pileup histogram not found in root file " << PU_rootfile_MC.c_str() << std::endl;
          throw std::exception();
        }

    pu_histo_mc->SetName(Form("%s_MC",pu_histo_mc->GetName()));
    pu_histo_mc->Scale(1.0/pu_histo_mc->Integral());
    hout->cd("pileup_weights");
    pu_histo_mc->Write();
    delete file_mc;

    }
else {
        std::cout << "root file could not be opened: " << PU_rootfile_MC.c_str() << std::endl;
          throw std::exception();
        }

  

    for(unsigned int i = 0; i < PU_rootfile_Data.size(); i++) {


   
      TH1::AddDirectory(kFALSE);
      TFile *file_data = new TFile(PU_rootfile_Data.at(i).c_str());
      if(file_data) {
        TH1D *h_tmp = (TH1D*)(file_data->Get("pileup"));//->Clone();
        if(h_tmp) {
        pu_histos_data.push_back(h_tmp);
        pu_histos_data.back()->SetName(Form("%s_%i",pu_histos_data.back()->GetName(),i));
        pu_histos_data.back()->Scale(1.0/pu_histos_data.back()->Integral());
        hout->cd("pileup_weights");
        pu_histos_data.back()->Write();
        } else {
          std::cout << "pileup histogram not found in root file " << PU_rootfile_Data.at(i).c_str() << std::endl;
          throw std::exception();
        }
        delete file_data;
   
      }
      else {
        std::cout << "root file could not be opened: " << PU_rootfile_Data.at(i).c_str() << std::endl;
          throw std::exception();
        }


      TH1::AddDirectory(kTRUE);





    }

   

    for(unsigned int i = 0; i < PU_scenario.size(); i++) {
      std::stringstream name;
      name << "weights_" << PU_scenario.at(i);
      h_weights.push_back(new TH1D(name.str().c_str(), name.str().c_str(), 50, 0.0, 5.0));
      std::stringstream name2d;
      name2d << "weights_vs_nPUTruth_" << PU_scenario.at(i);

      h_weights_vs_nPUTruth.push_back(new TH2D(name2d.str().c_str(), name2d.str().c_str(), 100, 0.0, 5.0, 80,0.0, 80.0));

    }
    //
    for(unsigned int i = 0; i < PU_scenario.size(); i++) {
      TH1D *h_ratio = (TH1D*) pu_histos_data.at(i)->Clone();
      h_ratio->SetName(Form("%s_ratio_data_over_mc_%i",h_ratio->GetName(),i));
      h_ratio->Divide(pu_histos_data.at(i), pu_histo_mc);
      h_ratio->Write();
    }

    hout->cd();

  }



  double GetWeight(double x) {
    double pileup_weight =1.0;
    if(_doMC && pileupwtx > 0) {
   
      int i = pileupwtx - 1;
     
      if(x >= pu_histo_mc->GetXaxis()->GetBinLowEdge(1) && x <= pu_histo_mc->GetXaxis()->GetBinUpEdge(pu_histo_mc->GetXaxis()->GetLast())) {
        const double pu_data = pu_histos_data.at(i)->GetBinContent(pu_histos_data.at(i)->FindBin(x));
        const double pu_mc = pu_histo_mc->GetBinContent(pu_histo_mc->FindBin(x));
        pileup_weight = pu_data / pu_mc;
         
        
      }
     
     
    }
    if(pileupwtx > 0) {
 h_weights.at(pileupwtx-1)->Fill(pileup_weight);
    h_weights_vs_nPUTruth.at(pileupwtx-1)->Fill(pileup_weight,x);
    }
 return pileup_weight;

  }
  const bool _doMC;
  const int pileupwtx;
  std::vector<TH1D*> h_weights;
  std::vector<TH2D*> h_weights_vs_nPUTruth;
  
  std::vector<TH1D*> pu_histos_data;
  std::vector<std::string> PU_rootfile_Data;
  std::string PU_rootfile_MC;
  std::vector<std::string> PU_scenario;
  TH1D *pu_histo_mc;
};




#endif
