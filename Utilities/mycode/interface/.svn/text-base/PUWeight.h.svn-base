#include <TFile.h>
#include <TH1.h>
#include <TF1.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>
#include <assert.h>


class PUWeight {

 public:
  PUWeight(TString baseDir);
  ~PUWeight();

  void InitPUWeightsS6(TString puFall11DataFileName,
		       TString puFall11MCFileName);

  float PUWeightS6(int nPU);

 private:

  TString _baseDir;

  TH1F * _puFall11MC;
  TH1F * _puFall11Data;


};
