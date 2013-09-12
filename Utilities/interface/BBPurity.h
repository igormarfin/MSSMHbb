#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <TString.h>

class BBPurity {

public:

  BBPurity(TString fileName);
  ~BBPurity();

  float getBBPurity(int tagger, int categ, float dijetMass, float svMass1, float svMass2);
  float getBBPurity(int tagger, int categ, float dijetMass, int iMass1, int iMass2);

private:

  int getSVMassBin(float svMass) {

    int bin = 0;
    
    if (svMass>1&&svMass<=2)
      bin = 1;

    if (svMass>2)
      bin = 2;

    return bin;

  }

  TString taggers[4];
  TString categories[3];
  TString btagConfigurations[3][6];

  TF1 * func[4][3][6];

};

