#ifndef Overlap_h
#define Overlap_h


#include <fstream>
#include <vector>
#include <algorithm>

void Overlap(char * nameLep, char * nameHad) {

  std::ifstream fileLep(nameLep);
  std::ifstream fileHad(nameHad);

  std::vector<long int> LumiLep;
  std::vector<long int> RunLep;
  std::vector<long int> EventLep;

 long  int lumi, run, event;
  while (fileLep >> run >> event >> lumi) {
    LumiLep.push_back(lumi);
    RunLep.push_back(run);
    EventLep.push_back(event);
  }

  std::cout << "Vector size = " << LumiLep.size() << std::endl;

  int nEventsLep = LumiLep.size();
  int nEventsHad = 0;
  int nEventsOverlap = 0;

  while (fileHad >> run >> event >> lumi) {

    for (int iE=0; iE<nEventsLep; ++iE) {
      if (lumi==LumiLep.at(iE) && run==RunLep.at(iE) && event==EventLep.at(iE)) {
	nEventsOverlap++;
//	std::cout<<"run= "<<run<<" event= "<<event<<" lumi= "<<lumi<<std::endl;
	break;
      }

//	std::vector<long int>::iterator it_lum=std::find(LumiLep.begin(),LumiLep.end(),lumi);
//	std::vector<long int>::iterator it_run=std::find(RunLep.begin(),RunLep.end(),run);
//	std::vector<long int>::iterator it_evt=std::find(EventLep.begin(),EventLep.end(),event);

//	if (it_lum!=LumiLep.end() && it_run!=RunLep.end() && it_evt!=EventLep.end()) nEventsOverlap++;
    }
    nEventsHad++;
  }

  int nEventsTot = nEventsLep + nEventsHad - nEventsOverlap;

  float fractLep = float(nEventsOverlap)/float(nEventsLep);
  float fractHad = float(nEventsOverlap)/float(nEventsHad);
  float fractTot = float(nEventsOverlap)/float(nEventsTot);

  std::cout << "Overlaping events : " << nEventsOverlap << std::endl;
  std::cout << "Leptonic : " << nEventsLep << "  ;  " << fractLep << std::endl;
  std::cout << "Hadronic : " << nEventsHad << "  ;  " << fractHad << std::endl;
  std::cout << "Total    : " << nEventsTot << "  ;  " << fractTot << std::endl;

}


  std::vector<long int> LumiLep_global;
  std::vector<long int> RunLep_global;
  std::vector<long int> EventLep_global;
bool _wasIni=false;


void Overlap_Ini(const char * nameLep)
{

std::cout<<"Overlap file "<<nameLep<<std::endl;
  std::ifstream fileLep_global(nameLep);


 long  int _lumi, _run, _event;
  while (fileLep_global >> _run >> _event >> _lumi) {
    LumiLep_global.push_back(_lumi);
    RunLep_global.push_back(_run);
    EventLep_global.push_back(_event);
_wasIni=true;

//std::cout<<RunLep_global.back()<<" "<<EventLep_global.back()<<" "<< LumiLep_global.back()<<std::endl;

  }


//std::cout<<"size 1 "<< LumiLep_global.size()<<std::endl;
//std::cout<<"size 2 "<< RunLep_global.size()<<std::endl;
//std::cout<<"size 3 "<< EventLep_global.size()<<std::endl;
return;

}

///return true if overlap occurs
///return false if no overlap occurs

bool Overlap_Check(long int run,long int event, long int lumi )
{

if (_wasIni && LumiLep_global.empty()) {std::cout<<"Something wrong. Exit"<<std::endl; return false;}
if (!_wasIni)  return false;

//        std::vector<long int>::iterator it_lum=std::find(LumiLep_global.begin(),LumiLep_global.end(),lumi);
//        std::vector<long int>::iterator it_run=std::find(RunLep_global.begin(),RunLep_global.end(),run);
//        std::vector<long int>::iterator it_evt=std::find(EventLep_global.begin(),EventLep_global.end(),event);

    for (unsigned int iE=0; iE<LumiLep_global.size(); ++iE) {
      if (lumi==LumiLep_global.at(iE) && run==RunLep_global.at(iE) && event==EventLep_global.at(iE)) 
     return true;
      }

return false;

 //      return  (it_lum!=LumiLep_global.end() && it_run!=RunLep_global.end() && it_evt!=EventLep_global.end());


}

#endif
