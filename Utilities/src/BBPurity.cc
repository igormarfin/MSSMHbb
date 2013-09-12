#include "Analysis/Utilities/interface/BBPurity.h"

BBPurity::BBPurity(TString fileName) {

  TFile * file = new TFile(fileName);

  taggers[0] = "TCHPT"; 
  taggers[1] = "TCHP6"; 
  taggers[2] = "CSVT";
  taggers[3] = "SSVHPT"; 

  categories[0] = "Xbb";
  categories[1] = "bXb";
  categories[2] = "bbX";
  
  btagConfigurations[0][0] = "00";
  btagConfigurations[0][1] = "10";
  btagConfigurations[0][2] = "20";
  btagConfigurations[0][3] = "01";
  btagConfigurations[0][4] = "11";
  btagConfigurations[0][5] = "21";

  btagConfigurations[1][0] = "00";
  btagConfigurations[1][1] = "10";
  btagConfigurations[1][2] = "20";
  btagConfigurations[1][3] = "01";
  btagConfigurations[1][4] = "11";
  btagConfigurations[1][5] = "21";

  btagConfigurations[2][0] = "Sum0";
  btagConfigurations[2][1] = "Sum1";
  btagConfigurations[2][2] = "Sum2";
  btagConfigurations[2][3] = "Sum0";
  btagConfigurations[2][4] = "Sum1";
  btagConfigurations[2][5] = "Sum2";

  for (int itagger=0; itagger<4; ++itagger) {
    for (int icat=0; icat<3; ++icat) {
      for (int ibconf=0; ibconf<6; ++ibconf) {
	TString funcName = "bb_" + taggers[itagger] + "_" + categories[icat] + "_" + btagConfigurations[icat][ibconf];
	func[itagger][icat][ibconf] = (TF1*)file->Get(funcName);
      }
    }
  }

}

BBPurity::~BBPurity() {

}

float BBPurity::getBBPurity(int tagger, int categ, float dijetMass, int n1, int n2) {


  int iconf = 0;
  if (categ==2) { // bbX
    int sum = n1 + n2;
    if (sum<2) 
      iconf = 0;
    else if (sum<3) 
      iconf = 1;
    else 
      iconf = 2;
  }
  else { // bXb or Xbb
    int T2 = 0;
    if (n2>1)
      T2 = 3;
    iconf = n1 + T2;
  }

  double xminD = 0;
  double xmaxD = 500;
  func[tagger][categ][iconf]->GetRange(xminD,xmaxD);

  float xmin = float(xminD);
  float xmax = float(xmaxD);

  float non2bb = 0;

  if (dijetMass<xmin)
    non2bb = func[tagger][categ][iconf]->Eval(xmin);
  else if (dijetMass>xmax)
    non2bb = func[tagger][categ][iconf]->Eval(xmax);
  else
    non2bb = func[tagger][categ][iconf]->Eval(dijetMass);

  if (non2bb<0) non2bb = 0;
 
  float purity = 1/(1+non2bb);

  return purity;

}

float BBPurity::getBBPurity(int tagger, int categ, float dijetMass, float svMass1, float svMass2) {

  int n1 = getSVMassBin(svMass1);
  int n2 = getSVMassBin(svMass2);

  float purity = getBBPurity(tagger, categ, dijetMass, n1, n2);

  return purity;

}
