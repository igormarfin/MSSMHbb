#include "Analysis/Utilities/interface/PUWeight.h"

PUWeight::PUWeight(TString baseDir) {

  _baseDir = baseDir;

}

PUWeight::~PUWeight() {

}

void PUWeight::InitPUWeightsS6(TString puFall11DataFileName,
			       TString puFall11MCFileName) {


  TFile * puFall11MCFile = new TFile(_baseDir+"/"+puFall11MCFileName);
  TFile * puFall11DataFile = new TFile(_baseDir+"/"+puFall11DataFileName);

  if (puFall11MCFile && puFall11DataFile)  {
    _puFall11MC = (TH1F*) puFall11MCFile->Get("pileup");
    _puFall11Data = (TH1F*) puFall11DataFile->Get("pileup");

    /*
      The applied rebinning is not really necessary as the
      reweighting function uses FindBin(). The same applies to the
      range (0,25) for the MC histogram. However, the histograms for
      
      both data and MC have the same binning to avoid edge/binning
      effects. A coarse binning should also limit the effect of the
      missing random seed during the pile-up creation process.
    */
    _puFall11Data->Rebin(10.);
    _puFall11MC->Scale(1./_puFall11MC->Integral());
    _puFall11Data->Scale(1./_puFall11Data->Integral());
    for (int idx1 = 1; idx1 < _puFall11MC->GetNbinsX(); idx1++)
      std::cout << idx1 << "\t" << _puFall11MC->GetBinContent(idx1) << "\t" << _puFall11Data->GetBinContent(idx1) << "\n";
  }
  else {
    _puFall11MC = 0;
    _puFall11Data = 0;
  }
  
}

float PUWeight::PUWeightS6(int nPU) {

  if (!_puFall11MC || !_puFall11Data)
    return 1.;

  float mean = float(nPU);
  float weightMC = _puFall11MC->GetBinContent(_puFall11MC->FindBin(mean));
  float weightData = _puFall11Data->GetBinContent(_puFall11Data->FindBin(mean));

  if (weightMC==0.)
    return 1.;
  else
    return weightData/weightMC;


}
