#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLorentzVector.h"
#include "TRFIOFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"

#include "Analysis/Utilities/interface/HBBTo4B_MC.h"
#include "Analysis/Utilities/interface/PUWeight.h"

int main(int argc, char * argv[])
{
  // usage                                                                                                                                       
  // > PUWeight [filelist]                                                                                                               
  // example filelist is Analysis/HbbMSSMAnalysis/test/SUSYBBHToBB_M-120                                                                              

  std::ifstream inputfile(argv[1]);

  int numberOfFiles = 0;
  inputfile >> numberOfFiles;

  std::string filelist(argv[1]);
  TString FileName(filelist);

  TString ChainName("hbbanalysis/HBBTo4B");
  TString baseDir("/afs/naf.desy.de/user/r/rasp/public/HTauTauAnalysis/RooT");

  PUWeight * pileupWeight = new PUWeight(baseDir);

  pileupWeight->InitPUWeightsS6("PileUp_2011.root","PileUp_Fall11.root");

  for (int iFile=0; iFile<numberOfFiles; ++iFile) {

    TString filename;
    inputfile >> filename;
    TFile * file = new TFile(filename);
    TTree * tree = (TTree*)file->Get(ChainName);
    int numberOfEntries = tree->GetEntries();
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    std::cout << std::endl;
    std::cout << filename << " : " << numberOfEntries << std::endl;

    HBBTo4B_MC * hbbTo4b = new HBBTo4B_MC(tree);
    hbbTo4b->Init(tree);

    for (int iE=0; iE<numberOfEntries; ++iE) {

      hbbTo4b->GetEntry(iE);

      int nPU = hbbTo4b->NumberOfPUInTime;

      float puWeight = pileupWeight->PUWeightS6(nPU);

      std::cout << "number of PU : " << nPU << "  PU weight = " << puWeight << std::endl;

    }
  }
}

