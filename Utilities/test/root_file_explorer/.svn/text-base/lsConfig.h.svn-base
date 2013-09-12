#ifndef LSCONFIG_H
#define LSCONFIG_H

//This Class is used to hold all config options for lsRoot

#include <set>
#include <vector>
#include <map>
#include <fstream>

#include <TString.h>
#include <TPostScript.h>
#include <TObject.h>
#include <TCanvas.h>

#include "lsUtils.h"

class lsConfig {

 public:
  lsConfig();
  ~lsConfig();

  int set(int, const char**,  std::vector<TString> &);
  void setInitialPrefs();
  void Initialize();

  void setColors();

  void getColor(TString &color, const TString& className);
  double getColor(const TString &color);

  void sortKeys( std::vector<JKey*> &objs, std::vector<JKey*> &dirs );
  
  // private:

  //generic globals
  int rootOpens;
  int rootCloses;
  bool quiet; // = 0; //If true suppress some outputs
  bool doNewLine; // = 1; //print extra new line after end of each new directory
  bool readObj; // = 0;//Do not unpack TObject from TKey if possible
  bool longFormat; // = 0;
  bool showRootInfo; // = 0;
  bool isTerminal; // = 0;//Initialize()
  bool isPipe; // = 0;//Initialize()
  bool doRecursion; // = 0;//-Z
  bool makeHumanReadable; // = 0;//-h
  bool doRootRecursion; // = 0;//-R
  bool displayClassic; // = 1;//??
  bool m_multiCol; // = 1;//let lsRoot default be -C
  bool m_reverse; // = 0;//-r
  bool m_ave; // = 0;//-A
  bool m_stdDev; // = 0;//-D
  bool m_N;// = 0 ;//-N
  bool m_integral; // = 0;//-I
  int m_N_binN; // = -1;//-N[]
  bool m_NBins; // = 0;//-B
  bool m_doColor; // = 1;//-G
  int m_nCols;
  int m_nRows;
  int m_precision; // = 3;
  //for -S
  TString m_scriptName;
  fstream m_scriptFile;
  TPostScript* m_script_ps;
  TCanvas* m_can;
  std::vector<TObject*> m_obj;
  TString m_can_paths;
  TString m_histName;
  std::map<TString, histDisplayOpt> m_histNames;//map of regex w/ displayOpts
  std::vector<histDisplayOpt> m_histDisplayOpt;//redundat info for histDisplayOpt but seems best since map has regex not pure string
  int m_histNumber;
  int m_canRow;
  int m_canCol;
  unsigned m_canPerPage;
  //for -m
  bool map_directories; // = 0;
  //for --cycle
  TString m_cycle; // = "last";
  //for -T[TT]/--tree
  bool m_printTree; // = 0;
  bool m_printTree_topOnly; // = 0;
  bool m_printTree_verbose; // = 0;
  TString m_tree_code;
  std::vector<TString> m_tree_code_v;
  //for --custom
  TString m_customFormat;
  std::vector<TString> m_customFormat_v;
  //for -d
  bool m_dir; // = 0;
  bool m_longDir; // = 0;
  //for sorting
  bool m_sortByAlpha; // = 1;
  bool m_sortByDate; // = 0;
  bool m_sortBySize; // = 0;
  bool m_sortByNBins; // = 0;
  bool m_sortByN; // = 0;
  //still used??
  //for -p/-P
  TString m_initialRootPath; // = "";
  bool m_initialPath_regex; // = false;
  //for --draw/--interactive/-i
  bool m_drawAsci; // = false;
  bool m_interactive; // = false;
  long nBytes; // = 0;
  //for --types
  bool m_types; // = 0;
  std::set<TString> m_types_s;

  //Totals for --stats
  bool m_stats; // = 0;
  long m_totKeys;
  long m_totBins;
  long m_totEntries;
  long m_totBytes;
  long double m_totKeys_glob;
  long double m_totBins_glob;
  long double m_totEntries_glob;
  long double m_totBytes_glob;

  //--draw
  asciDrawOptions* m_drawOpt;// = 0;
  asciDrawOptions* m_drawOpt2;// = 0;


  std::vector<TString> m_colors;

  std::map< TString, TString> colorMap;
  std::map< TString, TString> utilMap;

  TString colorNorm;// = "\033[0m";

  
};

#endif 
