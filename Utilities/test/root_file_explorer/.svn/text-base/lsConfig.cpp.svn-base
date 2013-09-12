#include "lsConfig.h"
#include "lsPrint.h"
#include "stringUtils.h"

#include <TRegexp.h>
#include <TPRegexp.h>
#include <TList.h>
#include <TError.h>
#include <TSystem.h>

#include <iostream>
#include <algorithm>

//C-Style includes
#ifndef WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/ioctl.h>
#endif

using namespace std;

lsConfig::lsConfig() 
:
//generic globals
rootOpens(0),
rootCloses(0),
quiet(0), //If true suppress some outputs
doNewLine(1), //prextra new line after end of each new directory
readObj(0),//Do not unpack TObject from TKey if possible
longFormat(0),
showRootInfo(0),
isTerminal(1),//Initialize()
isPipe(0),//Initialize()
doRecursion(0),//-Z
makeHumanReadable(0),//-h
doRootRecursion(0),//-R
displayClassic(1),//??
m_multiCol(1),//let lsRoot default be -C
m_reverse(0),//-r
m_ave(0),//-A
m_stdDev(0),//-D
m_N(0),//-N
m_integral(0),//-I
m_N_binN(-1),//-N[]
m_NBins(0),//-B
m_doColor(1),//-G
m_precision(3),
//for -S
m_script_ps(0),
m_can(0),
//for -m
map_directories(0),
//for --cycle
m_cycle("last"),
//for -T[TT]/--tree
m_printTree(0),
m_printTree_topOnly(0),
m_printTree_verbose(0),
//for -d
m_dir(0),
m_longDir(0),
//for sorting
m_sortByAlpha(1),
m_sortByDate(0),
m_sortBySize(0),
m_sortByNBins(0),
m_sortByN(0),
//for -p/-P
m_initialPath_regex(0),
//for --draw/--interactive/-i
m_drawAsci(0),
m_interactive(0),
nBytes(0),
//for --types
m_types(0),
//Totals for --stats
m_stats(0),
//--draw
m_drawOpt(0),
m_drawOpt2(0)
{

}


lsConfig::~lsConfig(){}

int lsConfig::set(int argc, const char* argv[], vector<TString>& allArgs){

  //check if pipe
  #ifndef WIN32
  isPipe = (isatty(0)==0);
  #else
  m_doColor = false;
  #endif


  setInitialPrefs();
  
  TString arg0 = argv[1];
  int beginArgC = 1;
  if(arg0(0,1)=="-"){
    bool arg = true;
    while(arg){
      beginArgC++;
      if(arg0.BeginsWith("-") && !arg0.BeginsWith("--")){
	if(arg0.Contains("A")) m_ave = true;
	if(arg0.Contains("B")) m_NBins = true;
	if(arg0.Contains("N")) { 
	  m_N = true; 
	  if(arg0.Contains("N[")) {
	    TString binN = arg0( TRegexp("N\\[[0-9]*\\]") );
	    binN.ReplaceAll("N[","");
	    binN.ReplaceAll("]","");
	    m_N_binN = binN.Atoi();
	  }	     
	}
	if(arg0.Contains("D")) m_stdDev = true;
	if(arg0.Contains("I")) m_integral = true;
	if(arg0.Contains("d")) m_dir = true;
	if(arg0.Contains("dd")) m_longDir = true;
	if(arg0.Contains("e")) showRootInfo = true;
	if(arg0.Contains("C")) { displayClassic = true; m_multiCol = true; }
	if(arg0.Contains("c") && argc > beginArgC) { m_customFormat = argv[beginArgC]; ++beginArgC; }
	if(arg0.Contains("G")) m_doColor = false;
	if(arg0.Contains("i")) { m_interactive = true; m_drawAsci = true; }
	if(arg0.Contains("t")) m_sortByDate = true;
	if(arg0.Contains("l")) longFormat = true;
	if(arg0.Contains("r")) m_reverse = true;
	if(arg0.Contains("S") && argc > beginArgC+1) { m_scriptName = argv[ beginArgC ]; m_histName = argv[ beginArgC +1]; ++beginArgC; ++beginArgC; }
	if(arg0.Contains("m")) map_directories = true;
	if(arg0.Contains("Z")) doRecursion = true; 
	if(arg0.Contains("h")) makeHumanReadable = true;
	if(arg0.Contains("R")) doRootRecursion = true; 
	if(arg0.Contains("T")) { m_printTree_topOnly = true; m_printTree = true;}
	if(arg0.Contains("TT")) { m_printTree = true; m_printTree_topOnly = false;}
	if(arg0.Contains("TTT")) { m_printTree = true; m_printTree_topOnly = false; m_printTree_verbose = true;}
	if(arg0.Contains("P") && argc > beginArgC) { m_initialRootPath = argv[beginArgC]; m_initialPath_regex = true; ++ beginArgC; } 
	if(arg0.Contains("p")) {
	  if(argc > beginArgC) { 
	    m_initialRootPath = argv[beginArgC]; 
	    if(m_initialRootPath == "regex"){
	      m_initialPath_regex = true;
	      if( argc >beginArgC+1) {
		m_initialRootPath = argv[beginArgC+1];
		++beginArgC;
	      }
	      else cerr << "ERROR    did not supply regex w/ -p regex";
	    }
	    ++beginArgC; 
	  }
	}
      }
      if(arg0.BeginsWith("--")){
	if(arg0.Contains("--draw")) m_drawAsci = true;
	if(arg0.Contains("--interactive")) {m_interactive = true; m_drawAsci = true;}
	if(arg0.Contains("--noColor")) m_doColor = false;
	if(arg0.Contains("--size")) m_sortBySize = true;
	if(arg0.Contains("--alpha")) m_sortByAlpha = true;
	if(arg0.Contains("--cycle") && argc > beginArgC) { m_cycle = argv[beginArgC]; ++beginArgC; }
	if(arg0.Contains("--precision") && argc > beginArgC) { m_precision = ((TString) argv[beginArgC]).Atoi(); ++beginArgC; }
	if(arg0.Contains("--custom") && argc > beginArgC) { m_customFormat = argv[beginArgC]; ++beginArgC; }
	if(arg0.Contains("--tree") && argc > beginArgC) { m_tree_code = argv[beginArgC]; ++beginArgC; }
	if(arg0.Contains("--sortN")) m_sortByN = true;
	if(arg0.Contains("--stats")) m_stats = true;
	if(arg0.Contains("--sortNBins")) { m_sortByN = false; m_sortByNBins = true;}
	if(arg0.Contains("--types") && argc > beginArgC) { m_types = true; TString arg1 = argv[beginArgC]; m_types_s = parseSet( arg1, ',' ); ++beginArgC;}
	if(arg0.Contains("help")){
	  TString next;
	  if(argc > beginArgC) {
	    next = argv[beginArgC];
	    if(next=="script") PrintScriptHelp();
	    else if(next=="path") PrintPathHelp();
	    else if(next=="color") showColors();
	    else if(next=="draw") PrintDrawHelp();
	    else if(next=="custom") PrintCustomHelp();
	    else PrintHelp();
	    return 0;
	  }
	  else {
	    PrintHelp();
	    return 0;
	  }
	}
      }
      arg0 = argv[beginArgC];
      if(!(arg0(0,1)=="-")) arg = false;
    }
  }  

  //support lsRoot -p 'fname:/path/blah'
  TString initialRootPath_fName, res;
  if( splitString(m_initialRootPath, initialRootPath_fName, res, ':') && m_initialRootPath.Contains(":/") ){
    initialRootPath_fName = initialRootPath_fName(0, initialRootPath_fName.Length());
    m_initialRootPath(initialRootPath_fName+":/") = "";
  }

  //loop over remaing args or cin if receiving pipe
  if(initialRootPath_fName.Length()) allArgs.push_back(initialRootPath_fName);
  if(!isPipe){ for(int i = beginArgC; i < argc; ++i) allArgs.push_back(argv[i]); }
  else { //read from cin if from pipe
    char arg1[256];
    while(cin.getline(arg1,256)) allArgs.push_back(arg1);
  }


  // so that 'lsRoot -opt'  will assume entire directory
  if( allArgs.size() == 0 ) {    
    char buffer[512];
    getcwd(buffer,512);
    TString buffer_t = buffer;
    if(!buffer_t.EndsWith("/")) buffer_t+="/";  
    allArgs.push_back(buffer_t);
  }
  

  Initialize();
  return beginArgC;
}

void lsConfig::Initialize(){
  
  //   gErrorIgnoreLevel = 1001;//suppress root info messages
  //   gErrorIgnoreLevel = 2001;//suppress root warning messages
  if(showRootInfo) gErrorIgnoreLevel = 1001;
  else gErrorIgnoreLevel = 2001;

  if( m_cycle != "all" && m_cycle != "last" && !m_cycle.IsDigit() ) {
    cout << "WARNING   " << "--cycle " << m_cycle << " invalid option.  Use 'all','last', or 'n=1,2,3,...'" << endl;
    m_cycle = "last";
  }
      
  #ifndef WIN32
  //Determine if writing to terminal or not  
  if(isatty(1)) isTerminal=true;
  else isTerminal=false;
  //Determine window size
  struct winsize ws;
  ioctl(1, TIOCGWINSZ, &ws);
  m_nCols = ws.ws_col;
  m_nRows = ws.ws_row;
  #else
  m_nCols=60;
  m_nRows=80;
  #endif

  m_totKeys_glob = 0;
  m_totBins_glob = 0;
  m_totEntries_glob = 0;
  m_totBytes_glob = 0;

  //If m_terminal=0 We must turn off color and multi column
  if(!isTerminal){
    m_multiCol = false;//turn off multiCol
    m_nRows = 58;
    m_nCols = 200; 
    m_doColor = false;
  }  

  setColors();

  //Revert to long format if T TT B A N D I or some flavor of --draw ( -i, --interactive )
  if(m_N || m_NBins|| m_ave || m_stdDev || m_integral|| m_printTree_topOnly || m_printTree || m_drawAsci || !m_tree_code.IsNull() ) { longFormat = true; readObj = true; }

  if(!m_customFormat.IsNull()) {
    std::map<TString, TString>::const_iterator m_itr = utilMap.find(m_customFormat);
    if(m_itr != utilMap.end()) m_customFormat = (*m_itr).second;
    longFormat = true;
    m_customFormat_v = parseList(m_customFormat, "\\n");
    if(m_customFormat.Contains("\\N") || m_customFormat.Contains("\\B") || m_customFormat.Contains("\\A") || m_customFormat.Contains("\\D") || m_customFormat.Contains("\\I") || m_customFormat.Contains("\\F") ) 
      readObj = true;
  }

  if(!m_customFormat.IsNull()) {
    m_customFormat_v = parseList(m_customFormat, "\\n");
  }

  if(longFormat==true) {m_multiCol = false; displayClassic = false;}
  //Turn of sortByAlpha if another sort invoked
  if(m_sortByDate || m_sortBySize || m_sortByN || m_sortByNBins) m_sortByAlpha = false;
  //Determine if readObj is necessary
  if( m_sortByNBins || m_sortByN || m_sortBySize || m_stats || m_scriptName!="") readObj = true;

  if( !m_tree_code.IsNull() ) m_tree_code_v = parseList(m_tree_code, "\\n");
}

void lsConfig::setInitialPrefs(){
  TString colorFile = ((TString) gSystem->Getenv ("HOME") )+"/.lsRoot_colors";
  fstream colors(colorFile,ios::in);
  while ( !colors.eof() && colors.is_open() ){
    string line_;
    getline( colors, line_);
    TString line = line_;
    if( line.BeginsWith("#") || line.Sizeof() < 5) continue;
    if( line.Contains("#") ) {
      TString newline, res;
      splitString( line, newline, res, '#' );
      line = newline;
    }
    if( !line.BeginsWith("m_") ) continue;
    TString parm, val;
    splitString(line, parm, val, ' ');
    parm.ReplaceAll(" ","");
    if(parm=="m_precision") m_precision = val.Atoi();
    if(parm=="m_quiet") quiet = val.Atoi();
    if(parm=="m_doNewLine") doNewLine = val.Atoi();
  }
}

void lsConfig::setColors(){
  m_colors.clear();

  TString colorFile = ((TString) gSystem->Getenv("HOME"))+"/.lsRoot_colors";
  fstream colors(colorFile,ios::in);
  while ( !colors.eof() && colors.is_open() ){
    string line_;
    getline( colors, line_);
    TString line = line_;
    if( line.BeginsWith("#") || line.Sizeof() < 5) continue;
    if( line.Contains("#") ) {
      TString newline, res;
      splitString( line, newline, res, '#' );
      line = newline;
    }
    TString type, color;
    splitString(line, type, color, ' ');
    type.ReplaceAll(" ","");
    TString original_color = color;
    color.ReplaceAll(" ","");

    if(type=="nRows") {m_nRows = color.Atoi(); continue; }
    if(type=="nCols") {m_nCols = color.Atoi(); continue; }

    if( type == "m_customFormat" ) {
      if( m_customFormat.IsNull() && longFormat ) m_customFormat = original_color;	
    }
    if( type == "customFormat" ){
      TString custom_name, custom_string;
      splitString( original_color, custom_name, custom_string , '=');
      custom_name.ReplaceAll(" ","");
      std::map<TString,TString>::iterator itr = utilMap.find(custom_name);
      if(itr==utilMap.end()) utilMap.insert(make_pair(custom_name, custom_string));
      else (*itr).second = custom_string;      
      continue;
    }
    if( type == "SetPalette" && m_doColor) {
      vector<TString> theColors = parseList(color, ',');
      for(vector<TString>::const_iterator itr = theColors.begin(); itr != theColors.end(); ++itr){
	TString color_tmp = "\033["+(*itr)(1,1)+";"+(*itr)(3,2)+"m";    
	m_colors.push_back(color_tmp);
      }
      continue;
    }
    if( type == "SetPaletteAsci" && !m_doColor) {
      vector<TString> theColors = parseList(color, ',');
      for(vector<TString>::const_iterator itr = theColors.begin(); itr != theColors.end(); ++itr){
	TString color_tmp = *itr;
	if(color_tmp!="null") m_colors.push_back(color_tmp);
	else m_colors.push_back(" ");
      }
      continue;
    }
    if(type!="FillSymbol" && type!="SetOptStat" && type!="LogyThreshold") color = "\033["+color(1,1)+";"+color(3,2)+"m";    
    std::map<TString,TString>::iterator itr = colorMap.find(type);
    if(itr==colorMap.end()) colorMap.insert(make_pair(type, color));
    else (*itr).second = color;
    if(type.Contains("default") && m_doColor) colorNorm = color;//quick method to get default color
    else if (!m_doColor) colorNorm = "";
  }  
  if(m_colors.size()==0 && m_doColor){
    m_colors.push_back("\033[0;31m");  m_colors.push_back("\033[0;33m");  m_colors.push_back("\033[0;32m");
    m_colors.push_back("\033[0;36m");  m_colors.push_back("\033[0;34m");  m_colors.push_back("\033[1;35m"); m_colors.push_back("");      
  }
  else if(m_colors.size()==0) {
    m_colors.push_back("*");
    m_colors.push_back("+"); 
    m_colors.push_back("O"); 
    m_colors.push_back("o"); 
    m_colors.push_back("."); 
    m_colors.push_back(" "); 
  }
  if(!m_doColor) { 
    colorMap.clear(); 
  }
}

void lsConfig::getColor(TString &color, const TString& className){

  if(this->colorMap.size()==0 || className=="default") color = this->colorNorm;
  else {
    std::map<TString,TString>::const_iterator itr = this->colorMap.find(className);
    if(itr==this->colorMap.end()) { color = this->colorNorm; return;}
    else color = (*itr).second; 
  }    
  if(this->colorMap.size()==0 && className=="FillSymbol") color = "XO";
}

double lsConfig::getColor(const TString &className){

  TString color;
  if(this->colorMap.size()==0 || className=="default") color = this->colorNorm;
  else {
    std::map<TString,TString>::const_iterator itr = this->colorMap.find(className);
    if(itr==this->colorMap.end()) { color = this->colorNorm; return 100;}//this default value is assuming only LogyThreshold
    else color = (*itr).second; 
  }    
  return color.Atof();
}


void lsConfig::sortKeys( vector<JKey*> &objs, vector<JKey*> &dirs ){

  if(m_sortByAlpha){ //Sorting keys alphabetically (are sorted by time by root by default)
    if(!m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyComp);
      std::sort(dirs.begin(), dirs.end(), jkeyComp);
    }
    if(m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompR);
      std::sort(dirs.begin(), dirs.end(), jkeyCompR);
    }      
  }
  else if(m_sortBySize){
    if(!m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompSize);
      std::sort(dirs.begin(), dirs.end(), jkeyCompSize);
    }
    if(m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompSizeR);
      std::sort(dirs.begin(), dirs.end(), jkeyCompSizeR);
    }      
  }
  else if(m_sortByN){
    if(!m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompN);
      std::sort(dirs.begin(), dirs.end(), jkeyCompN);
    }
    if(m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompNR);
      std::sort(dirs.begin(), dirs.end(), jkeyCompNR);
    }      
  }
  else if(m_sortByNBins){
    if(!m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompNBins);
      std::sort(dirs.begin(), dirs.end(), jkeyCompNBins);
    }
    if(m_reverse){
      std::sort(objs.begin(), objs.end(), jkeyCompNBinsR);
      std::sort(dirs.begin(), dirs.end(), jkeyCompNBinsR);
    }      
  }

  if(!m_reverse && m_sortByDate){//by default if we do not sort keys, keys are in reverse time order so to get normal time order just reverse
    std::reverse(objs.begin(), objs.end());
    std::reverse(dirs.begin(), dirs.end());
  }

}
