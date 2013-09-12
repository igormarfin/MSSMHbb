#include <vector>
#include <TString.h>
#include <TH1F.h>
#include <TFile.h>

class Trigger {

public: 

  int filterPattern;
  std::vector<TString> genericTriggerList;
  
  Trigger(TFile * file, std::vector<TString> triggerList) {
    
    genericTriggerList.clear();

    TH1F* gtlHist = (TH1F*) file->Get("hbbanalysis/gtlHist");
    if (gtlHist != NULL) {
      for (int ibin=1; ibin<= gtlHist->GetXaxis()->GetNbins(); ++ibin) {
	genericTriggerList.push_back(gtlHist->GetXaxis()->GetBinLabel( ibin ));
	std::cout << gtlHist->GetXaxis()->GetBinLabel( ibin ) << std::endl;
      }
    }

    int numberOfGenericTriggers = genericTriggerList.size();
    int numberOfTriggers = triggerList.size();
    filterPattern = 0;

    for (int iT=0; iT<numberOfTriggers; ++iT) {
      bool triggerFound = false; 
      TString currentTriggerName = triggerList.at(iT);
      for (int ibin = 0; ibin<numberOfGenericTriggers; ++ibin) {
	TString triggerName = genericTriggerList.at( ibin );
	if (triggerName.Contains(currentTriggerName)) {
	  std::cout << "Trigger " << triggerName << " found in the slot " << ibin << " of generic trigger list" << std::endl;
	  triggerFound = true;
	  filterPattern = filterPattern | (1<<ibin);
	}
      }
      if (!triggerFound)
	std::cout << "Trigger " << currentTriggerName << " is not found in the generic trigger list " << std::endl;

    }

    std::cout << std::endl;
    std::cout << "Filter pattern = " << filterPattern << std::endl;
    std::cout << std::endl;

  }

  bool accept(int trgAccept) {
    return trgAccept & filterPattern; 
  }

  int getFilterPattern() {
    return filterPattern;
  }

  void PrintOut() {
    std::cout << "Printing all available triggers" << std::endl;
    int numberOfTriggers = genericTriggerList.size();
    for (int iT=0; iT<numberOfTriggers; ++iT) {
      std::cout << genericTriggerList.at(iT) << std::endl;
    }
    std::cout << std::endl;
  }


};
