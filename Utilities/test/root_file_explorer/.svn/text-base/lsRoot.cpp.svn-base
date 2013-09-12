#include "lsRoot.h"
#include "lsPrint.h"
#include "lsUtils.h"

#include <TClass.h>
#include <TTree.h>
#include <TPRegexp.h>
#include <TRegexp.h>
#include <TH1.h>
#include <TH3.h>
#include <THStack.h>
#include <TError.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TFormula.h>
#include <TGraphErrors.h>
#include <TSystem.h>

#include <vector>
#include <algorithm>

using namespace std;

lsConfig conf;

int main(int argc, const char* argv[]){

  //Here we set all global variables in this order:
  //   1. .lsRoot_colors 
  //   2. command-line args
  //   3. .lsRoot_colors again
  //We also extract all non '-/--' args into allArgs including those from standard input
  vector<TString> allArgs;
  conf.set( argc, argv, allArgs );

  if(argc==1 && !conf.isPipe) {
    PrintHelp();
    return 0;
  }

  //Loop over all args including those from cin
  for(vector<TString>::iterator itr = allArgs.begin(); itr != allArgs.end(); ++itr){
    TString file = *itr;
    //truly I should figure out how not to follow symbolic links
    if(file.Contains("afshome")) continue;
    
    void *dp = 0;
    if(!file.Contains(":/")) dp = gSystem->OpenDirectory(file);
    
    if(dp){
      //If arg is a bash directory
      TString bashPath = file+"/";
      listDirectory(dp, bashPath, conf.m_initialRootPath);
      gSystem->FreeDirectory(dp);
      continue;
    }
    
    openFile(file);
    if(!f) continue;
    
    TString bashPath = f->GetName();
    
    if(conf.m_initialRootPath!="") { 
      listDirectoryFromPath(conf.m_initialRootPath, bashPath);
      closeFile();
      continue;
    } 
    else if(conf.map_directories==true) mapDirectory(f,0);
    else if(conf.m_initialRootPath=="") listDirectory(f, gSystem->Getenv("CWD") );
    closeFile();
  }
  
  if( conf.m_stats){
    cout << "All " << conf.rootOpens << " Files " << conf.m_initialRootPath << " : " << "Total Keys:  ";
    cout << conf.m_totKeys+conf.m_totKeys_glob << "  Total Bytes:  "  << conf.m_totBytes+conf.m_totBytes_glob ;
    cout << "  Total Bins:  " << conf.m_totBins+conf.m_totBins_glob << " Total Entries: " << conf.m_totEntries+conf.m_totEntries_glob << endl;  
  }
  
  return 0;
}

void listDirectory(void* dp, TString bashPath, TString initialRootPath){
  const char* dir_entry;
  while((dir_entry = gSystem->GetDirEntry(dp)) != NULL){
    TString file = dir_entry;
    if(file=="." || file == "..")
      continue;
    void *dp2 = 0;
    
    //Fix me I don't want to follow symbolic links
    if(file.Contains("afshome")) continue;
    TString fullFile = bashPath+file;
    fullFile("//") = "/";
    dp2 = gSystem->OpenDirectory(fullFile);
    
    if(dp2 != NULL && conf.doRecursion){
      TString bashPath2=bashPath+file;
      bashPath2+="/";
      //      cout << "Opening " << bashPath2 << endl;
      listDirectory(dp2, bashPath2, initialRootPath);
      gSystem->FreeDirectory(dp2);
    }
    else {
      TString fullFile = bashPath+file;
      fullFile("//") = "/";
      openFile(fullFile);
      if(!f) continue;

      if(initialRootPath!="") { 
	listDirectoryFromPath(initialRootPath, bashPath);
	closeFile();
	continue;
      } 
      else if(conf.map_directories==true) mapDirectory(f,0);
      else if(initialRootPath=="") listDirectory(f, fullFile);

      closeFile();
    }
  }
}

bool listDirectory(TDirectory* dir, TString bashPath, TString pattern){

  bool hasPattern = (pattern!="");
  TString path = dir->GetPath();
  if( dir!=f ) path+="/";

  TString shortPath = path(path.Index(":/")+2,path.Sizeof()-path.Index(":/")-2);
  while(shortPath.Contains(":/")) shortPath = shortPath(shortPath.Index(":/")+2,shortPath.Sizeof()-shortPath.Index(":/")-2);

  if(!dir) {std::cerr << "ERROR" << "Directory Not Found" << std::endl; return false; }

  TIter itr(dir->GetListOfKeys());

  TKey* key = 0;
  vector<JKey*> dirs; // store keys of TDirectories of mother TDirectory
  vector<JKey*> objs; // store all objects of TDirectory mother
  int nKeys_loc = 0;
  int nBytes_loc = 0;
  int nBins_loc = 0;
  int nEntries_loc = 0;
  unsigned maxName = 10;

  TKey* oldkey = 0;
  while(( key = dynamic_cast<TKey*> (itr()))){
    TString keyName = key->GetName();
    if( oldkey && conf.m_cycle != "all" ){
      TString oldKeyName = oldkey->GetName();
      if( conf.m_cycle == "last" && oldKeyName == keyName ) continue;
      if (conf.m_cycle.IsDigit()){
	stringstream ss;
	ss << conf.m_cycle;
	int cycle;
	ss >> cycle;
	if( cycle != key->GetCycle() ) continue;
      }
    }
    oldkey = key;
    TString className = key->GetClassName();
    bool isDirectory = (className=="TDirectory"||className=="TDirectoryFile");
    //Complicated Classes that bareroot can't understand
    if(className.Contains("SMatrix") || className.Contains("TASImage")) continue;
    //continue if pattern
    if(hasPattern){
      TPRegexp reg = (TPRegexp) pattern;
      if(keyName(reg)!=keyName) continue;
    }
    //continue if not requested type
    bool requested_type = true;
    if(conf.m_types) requested_type = !(conf.m_types_s.find(className) == conf.m_types_s.end());
    if(!isDirectory && !requested_type) continue;
    TObject* obj = 0;
    if(conf.readObj || isDirectory || className=="TTree") obj = dynamic_cast<TObject*> (key->ReadObj());

    if(conf.m_stats && requested_type) nBytes_loc+=key->GetKeylen();
    if(conf.m_stats && requested_type){
      TH1* h = dynamic_cast<TH1*> (obj);
      TTree* t = dynamic_cast<TTree*> (obj);
      nKeys_loc++;      
      if(h) { nBins_loc+=h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ(); nEntries_loc+=(long) h->GetEntries(); }
      if(t) { nBytes_loc+=t->GetTotBytes(); nEntries_loc+=(long) t->GetEntries(); }
      else { nBytes_loc+=key->GetNbytes(); }
    }

    //add to keyName if -d or -dd
    if(conf.m_longDir && path!="") {
      keyName = path+keyName;
    }
    else if(conf.m_dir && shortPath!="") keyName = shortPath+keyName;

    if(keyName.Sizeof() > (int) maxName) maxName = keyName.Sizeof();
    JKey* jkey = new JKey(obj, key, dir);
    objs.push_back(jkey);
    if(isDirectory){
      JKey* jkey_dir = new JKey(obj, key, dir);
      dirs.push_back(jkey_dir);
    }

  }
  
  conf.sortKeys( objs, dirs);

  //Display Current TDirectory w nKeys
  if(objs.size()) {  
    TString color;
    conf.getColor( color, "TDirectoryFile");
    cout << color << path << " " << conf.colorNorm << dir->GetNkeys() << " Keys" << endl;
  }

  if(conf.m_stats) {
    cout << "Total Keys:  " << nKeys_loc << "  Total Bytes:  " << nBytes_loc << "  Total Bins:  " << nBins_loc << " Total Entries: " << nEntries_loc << endl;  
    conf.m_totBytes+=nBytes_loc;
    conf.m_totEntries+=nEntries_loc;
    conf.m_totBins+=nBins_loc;
    conf.m_totKeys+=nKeys_loc;
  }

  //for -S
  if(conf.m_scriptFile.is_open() || !conf.m_scriptName.IsNull()){
    for( vector<JKey*>::const_iterator itr = objs.begin(); itr != objs.end(); ++itr){
      JKey* jkey = *itr;
      if(isInHistList(jkey->key->GetName())) addHistToScript(shortPath,jkey->obj);
    }
  }


  //Display multi-column
  if(conf.displayClassic && conf.m_multiCol && maxName+2 < (unsigned) conf.m_nCols/2){
    maxName+=2;
    int nCols = conf.m_nCols/maxName;
    int nRows = objs.size()/nCols;
    if(objs.size()%nCols!=0) nRows++;
    for(int row = 0; row != nRows; ++row){
      for(int col = 0; col != nCols; ++col){
	unsigned element = nRows*col+row;	
	if ( element >= objs.size()) break;
	JKey* jkey = objs.at(element);
	TString name = jkey->key->GetName();
	if(conf.m_longDir && path != "") {
	  name = path+name;
	}
	else if (conf.m_dir && shortPath != "") name = shortPath+name;
	TString className = jkey->key->GetClassName();
	addPad(name, maxName-name.Sizeof()); 
	printKeyClassic(name, className);
	bool isDirectory = (className=="TDirectory"||className=="TDirectoryFile");
	if(!isDirectory) delete jkey;
      }
      cout << endl;
    }   
  }
  else if (conf.displayClassic) {
    for( vector<JKey*>::const_iterator itr = objs.begin(); itr != objs.end(); ++itr){
      JKey* jkey = *itr;
      TString className = jkey->key->GetClassName();
      TString name = jkey->key->GetName();
      if(conf.m_longDir && path != "") {
	name = path+name;
      }
      else if (conf.m_dir && shortPath != "") name = shortPath+name;
      bool isDirectory = (className=="TDirectory"||className=="TDirectoryFile");
      printKeyClassic(name, className); 
      cout << endl;
      if(!isDirectory) delete jkey;
    }
  }    
  else {
    for( vector<JKey*>::const_iterator itr = objs.begin(); itr != objs.end(); ++itr){
      JKey* j = *itr;
      TDirectory* dir = j->dir;
      TKey* key = j->key;
      TObject* obj = j->obj;
      TString keyName = key->GetName();
      TString className = j->key->GetClassName();
      bool isDirectory = (className=="TDirectory"||className=="TDirectoryFile");

      if(conf.m_drawAsci) {
	drawAsci((TH1*) obj);
	if( ((TString) key->GetClassName()).BeginsWith("TGraph"))  drawGraph((TGraph*) obj);
      }      

      if(!conf.m_customFormat.IsNull()) printKeyCustom(dir, key, obj, "");
      else if(conf.longFormat && !conf.m_dir ) printKey(dir, key, obj, "");
      else if(conf.longFormat && conf.m_dir && !conf.m_longDir) printKey(dir, key, obj, shortPath);
      else if(conf.longFormat && conf.m_longDir) printKey(dir, key, obj, path);
      else printKeyBasic(dir, key, obj);
      if(!isDirectory) delete j;
    }
  }
  if(objs.size() && conf.doNewLine) cout << endl;
  //perform recurse
  for( vector<JKey*>::const_iterator itr = dirs.begin(); itr != dirs.end(); ++itr){
    JKey* j = *itr;
    if(conf.doRootRecursion || (hasPattern && !conf.m_dir)) listDirectory(dynamic_cast<TDirectory*> (j->obj), bashPath, "");
    delete j;
  }

  return (objs.size() > 0);
}

void printKeyCustom(TDirectory* dir, TKey* key, TObject* obj, TString path){

  int pad = 7 + conf.m_precision;

  for(vector<TString>::const_iterator itr = conf.m_customFormat_v.begin(); itr != conf.m_customFormat_v.end(); ++itr){
    TString line = (*itr);
    int level = 1;
    TH1* h = dynamic_cast<TH1*> (obj);
    TTree* t = dynamic_cast<TTree*> (obj);
    if(line == "\\F"){
      if(!h) continue;
      for( int ibin = 1; ibin != h->GetNbinsX()+1; ++ibin){
        TString label = h->GetXaxis()->GetBinLabel(ibin);
	if(label.IsNull()) label = returnString( ibin );
        label += ":";
        label+=makePad(level, std::max( 30-label.Length(), 1) );
        cout << label << " " << returnString(h->GetBinContent(ibin), conf.m_precision) << endl;
      }
      continue;
    }
    if(line.Contains("\\p")) { level = 0; line.ReplaceAll("\\p",""); }
    TString color;
    conf.getColor(color, key->GetClassName());
    if(line.Contains("\\N")) {    
      int count = 0;
      while( line.Contains("\\N") && count < 1000){
	count++; //safeguard agains bad config
	TString val;
	TRegexp bincon_nx_y("\\\\N\\[[ ]*[0-9]*[ ]*-[ ]*[0-9]*[ ]*\\]");	
	TRegexp bincon_n("\\\\N\\[[ ]*[0-9]*[ ]*\\]");
	TRegexp bincon_all("\\\\N\\[[ ]*[Aa][Ll][Ll][ ]*\\]");
	if(line.Contains(bincon_all) && h){
	  TString res = line(bincon_all);
	  for( int ibin = 0; ibin !=  h->GetNbinsX() +1; ++ibin){ 
	      TString val_local = returnString(h->GetBinContent( ibin ), conf.m_precision);
	      val_local+=makePad(level, std::max(pad-val_local.Length(), 1));
	      val += val_local;
	  }
	  line.ReplaceAll( res, val );
	}
	else if(line.Contains(bincon_nx_y)){
	  TString res = line(bincon_nx_y);
	  TString bin_str = res; bin_str.ReplaceAll(" ",""); bin_str.ReplaceAll("[",""); bin_str.ReplaceAll("]",""); bin_str.ReplaceAll("\\N","");
	  TString bin1, bin2;
	  splitString(bin_str, bin1, bin2, '-');
	  int ibin1 = bin1.Atoi(); int ibin2 = bin2.Atoi();
	  if(h)
	    for( int ibin = std::max(ibin1, 0); ibin != std::min(ibin2+1, h->GetNbinsX() +1); ++ibin){ 
	      TString val_local = returnString(h->GetBinContent( ibin ), conf.m_precision);
	      val_local+=makePad(level, std::max(pad-val_local.Length(), 1));
	      val += val_local;
	    }
	  line.ReplaceAll( res, val );
	}
	else if(line.Contains(bincon_n)){
	  TString res = line(bincon_n);
	  TString bin_str = res; bin_str.ReplaceAll(" ",""); bin_str.ReplaceAll("[",""); bin_str.ReplaceAll("]",""); bin_str.ReplaceAll("\\N","");
	  if(h) val = returnString(h->GetBinContent( std::min( bin_str.Atoi(), h->GetNbinsX() + 1 ) ), conf.m_precision);
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll( res, val );
	}
	else {
	  if( h ) val = returnString( (double) h->GetEntries(), conf.m_precision );
	  else if( t ) val = returnString( (double) t->GetEntries(), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  //line.ReplaceAll("\\N", val );
	  line.Replace(line.Index("\\N"), 2, val);
	}
	if(count == 999) cerr << "ERROR   " << "bad config in --custom " << line << endl;
      }
    }
    if(line.Contains("\\c")) insertString( line, key->GetClassName(), "\\c", 1, level );
    if(line.Contains("\\I")) {    
      int count = 0;
      while( line.Contains("\\I") && count < 1000) {
	count++;
	TString val;
	TRegexp integral_x("\\\\I\\[[ ]*[0-9]*[ ]*-[ ]*[0-9]*[ ]*\\]");
	TRegexp integral_y("\\\\I\\[[ ]*[0-9]*[ ]*-[ ]*[0-9]*[ ]*,[ ]*[0-9]*[ ]*-[ ]*[0-9]*[ ]*\\]");			   
	if(line.Contains(integral_x)){
	  TString res = line(integral_x);
	  TString bin_str = res; bin_str.ReplaceAll(" ",""); bin_str.ReplaceAll("\\I",""); bin_str.ReplaceAll("[",""); bin_str.ReplaceAll("]","");
	  TString bin1, bin2;
	  splitString(bin_str, bin1, bin2, '-');
	  if(h) val = returnString(h->Integral( std::max(0, bin1.Atoi()), std::min( h->GetNbinsX() +1, bin2.Atoi()) ), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll(res,val);
	}
	else if (line.Contains(integral_y)){
	  TString res = line(integral_y);
	  if(!res.BeginsWith("\\I")) res = "\\"+res;
	  TString bin_str = res; bin_str.ReplaceAll(" ",""); bin_str.ReplaceAll("\\I",""); bin_str.ReplaceAll("[",""); bin_str.ReplaceAll("]","");
	  TString binsx, binsy;
	  splitString(bin_str, binsx, binsy, ',');
	  TString binx1, binx2, biny1, biny2;
	  splitString(binsx, binx1, binx2, '-');
	  splitString(binsy, biny1, biny2, '-');
	  TH2* h2 = dynamic_cast <TH2*> (h);
	  if(h2) val = returnString(h2->Integral( std::max(0, binx1.Atoi()), std::min( h2->GetNbinsX()+1, binx2.Atoi()), 
						  std::max(0, binx1.Atoi()), std::min( h2->GetNbinsX()+1, binx2.Atoi()) ), conf.m_precision); 
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll(res, val);
	}
	else {
	  if( h ) val = returnString( h->Integral(), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.Replace(line.Index("\\I"), 2, val);
	}
	if(count == 999) cerr << "ERROR   " << "bad config in --custom " << line << endl;
      }
    }
    if(line.Contains("\\D")) {  
      while(line.Contains("\\D")){
	if(line.Contains("\\DX")){
	    TString val;
	    if( h ) val = returnString( h->GetRMS(1), conf.m_precision );
	    val+=makePad(level, std::max(pad-val.Length(), 1));
	    line.ReplaceAll("\\DX", val );
	}
	else if(line.Contains("\\DY")){
	  TString val;
	  TH2* h2 = dynamic_cast<TH2*> (h);
	  TH3* h3 = dynamic_cast<TH3*> (h);
	  if( h2 || h3 ) val = returnString( h->GetRMS(2), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\DY", val );
	}
	else if(line.Contains("\\DZ")){
	  TString val;
	  TH3* h3 = dynamic_cast<TH3*> (h);
	  if( h3 ) val = returnString( h->GetRMS(3), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\DZ", val );
	}
	else {
	  TString val;
	  if( h ) val = returnString( h->GetRMS(), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\D", val );
	}
      }
    }
    if(line.Contains("\\A")) {    
      while(line.Contains("\\A")){
	if(line.Contains("\\AX")){
	    TString val;
	    if( h ) val = returnString( h->GetMean(1), conf.m_precision );
	    val+=makePad(level, std::max(pad-val.Length(), 1));
	    line.ReplaceAll("\\AX", val );
	}
	else if(line.Contains("\\AY")){
	  TString val;
	  TH2* h2 = dynamic_cast<TH2*> (h);
	  TH3* h3 = dynamic_cast<TH3*> (h);
	  if( h2 || h3 ) val = returnString( h->GetMean(2), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\AY", val );
	}
	else if(line.Contains("\\AZ")){
	  TString val;
	  TH3* h3 = dynamic_cast<TH3*> (h);
	  if( h3 ) val = returnString( h->GetMean(3), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\AZ", val );
	}
	else {
	  TString val;    
	  if( h ) val = returnString( h->GetMean(), conf.m_precision );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\A", val );
	}
      }
    }
    if(line.Contains("\\B")) {    
      while(line.Contains("\\B")){
	if(line.Contains("\\BX")){
	    TString val;
	    if( h ) val = returnString( h->GetNbinsX() );
	    val+=makePad(level, std::max(pad-val.Length(), 1));
	    line.ReplaceAll("\\BX", val );
	}
	else if(line.Contains("\\BY")){
	  TString val;
	  TH2* h2 = dynamic_cast<TH2*> (h);
	  TH3* h3 = dynamic_cast<TH3*> (h);
	  if( h2 || h3 ) val = returnString( h->GetNbinsY() );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\BY", val );
	}
	else if(line.Contains("\\BZ")){
	  TString val;
	  TH3* h3 = dynamic_cast<TH3*> (h);
	  if( h3 ) val = returnString( h->GetNbinsZ() );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\BZ", val );
	}
	else {
	  TString val;    
	  if( h ) val = returnString( h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ() );
	  val+=makePad(level, std::max(pad-val.Length(), 1));
	  line.ReplaceAll("\\B", val );
	}
      }
    }
    if(line.Contains("\\t")){
      TString val;
      val = key->GetDatime().AsString();
      val+=makePad(level, std::max(30-val.Length(), 1));
      line.ReplaceAll("\\t", val );
    }
    if(line.Contains("\\s")){
      TString val;
      long int nBytes = key->GetNbytes();
      TTree* t = dynamic_cast<TTree*> (obj);
      if( t ) {
	nBytes = t->GetTotBytes();
      }      
      if(conf.makeHumanReadable) {
	TString str = returnString( key->GetNbytes(), conf.m_precision);
	if(nBytes>1024 && conf.makeHumanReadable)  str = MakeReadable(nBytes);
	else if(conf.makeHumanReadable) str ="B";
	val = returnString(nBytes, std::max(conf.m_precision, 4));
	val += str;
      }
      else val = returnString(nBytes, conf.m_precision);
      val+=makePad(level, std::max(pad-val.Length(), 1));
      line.ReplaceAll("\\s", val );
    }
    if(line.Contains("\\k")) {
      TDirectory* dir = dynamic_cast<TDirectory*> (obj);
      TString val;
      if( dir ) {
	val = returnString(dir->GetNkeys(), 5);
      }
      val+=makePad(level, std::max(pad-val.Length() - 2, 1));
      line.ReplaceAll("\\k", val );
    }

    if(line.Contains("\\C")){
      TString val = returnString(key->GetCycle());
      val+=makePad(level, std::max(pad-val.Length() - 2, 1));
      line.ReplaceAll("\\C", val );      
    }

    //Anything that requires color last
    if(line.Contains("\\d")) {    
      TString path = dir->GetPath();
      TString shortPath = path(path.Index(":/")+2,path.Sizeof()-path.Index(":/")-2);
      while(shortPath.Contains(":/")) shortPath = shortPath(shortPath.Index(":/")+2,shortPath.Sizeof()-shortPath.Index(":/")-2);
      if( !shortPath.IsNull() ) shortPath+="/";
      shortPath+=makePad(level, std::max(35-shortPath.Length(), 1));
      shortPath.ReplaceAll(" ","");
      line.ReplaceAll("\\d", color+shortPath+conf.colorNorm );
    }
    if(line.Contains("\\h")) {    
      TString name = key->GetName();
      name+=makePad(level, std::max(35-name.Length(), 1));
      line.ReplaceAll("\\h", color+name+conf.colorNorm );
    }

    while(line.Contains("expr<")){
      TString expr = line(line.Index("expr<")+5, line.Index(">") - line.Index("expr<") - 5);
      TString val;
      
      if(expr.Length() && !(((TString) key->GetClassName())=="TDirectoryFile") ) {
	TFormula form("dummy",expr);
	val = returnString( form.Eval(-1), conf.m_precision ); 
      }
      val+=makePad(level, std::max(pad-val.Length(), 1) );
      line.ReplaceAll(expr, val);
      line.Replace(line.Index("expr<"), 5, "");
      line.Replace(line.Index(">"),1,"");
    }
    //  cout << color << line << conf.colorNorm << endl;
    while( line.EndsWith(" ") ) line = line(0, line.Length() - 1);
    cout << line << endl;
  }
}

void printKey(TDirectory* dir, TKey* key, TObject* obj, TString path){

  TString className = key->GetClassName();
  long nBytes = key->GetNbytes();
  if(className == "TTree" ) nBytes = ((TTree*) obj)->GetTotBytes();
  TString color;
  conf.getColor( color, className);
  TDirectory* dir2 = 0;
  if( (TString) key->GetClassName() == "TDirectoryFile" ) dir2 = dir->GetDirectory(key->GetName());
  TString str = "";
  TString classPad = makePad(1, 15-className.Sizeof());
  if(!(conf.m_dir||conf.m_longDir) && !dir2) {
    path = "";
  }

  int pad = 6 + conf.m_precision;
  if(pad < 9) pad = 9;
  int extra_pad = 0;//if user requests N[bin]
  if( conf.m_N_binN >=0 ) extra_pad = ((TString) "[" + returnString(conf.m_N_binN) + "]").Length();

  if(nBytes>1024 && conf.makeHumanReadable)  str = MakeReadable(nBytes);
  else if(conf.makeHumanReadable) str ="B";
  cout << key->GetDatime().AsString() << "  " << right << setprecision(conf.m_precision) << setw(pad) << nBytes << str << " " ;      
  cout << className << classPad << " " ;

  if ( conf.m_N || conf.m_ave || conf.m_stdDev || conf.m_NBins || conf.m_integral){
    TH1* h = dynamic_cast<TH1*> (obj);
    TTree* t = dynamic_cast<TTree*> (obj);
    if( h && conf.m_N ) { 
      if( conf.m_N_binN >= 0 && h->GetNbinsX() >= conf.m_N_binN ) cout << "N["<< conf.m_N_binN << "]=" << left << setprecision(conf.m_precision) << setw(pad) << h->GetBinContent(conf.m_N_binN) << " " ;
      else cout << "N=" << left << setprecision(conf.m_precision) << setw(pad+extra_pad) << h->GetEntries() << " " ;
    }
    else if (conf.m_N && !t) cout << "   " << setw(pad+extra_pad) << " " ;
    if( t && conf.m_N ) cout << "N=" << left << setprecision(conf.m_precision) << setw(pad+extra_pad) << t->GetEntries() << " " ;
    if( h && conf.m_NBins ) cout << "NB=" << left << setprecision(conf.m_precision) << setw(pad) << h->GetNbinsX()*h->GetNbinsY()*h->GetNbinsZ() << " " ;
    else if (conf.m_NBins ) cout << "    " << left << setw(pad) << " " ;
    if( h && conf.m_ave ) cout << "ave=" << left << setprecision(conf.m_precision) << setw(pad) << h->GetMean() << " " ;
    else if (conf.m_ave) cout << "     " << left << setw(pad) << " " ;
    if( h && conf.m_stdDev ) cout << "rms=" << left << setprecision(conf.m_precision) << setw(pad) << h->GetRMS() << " " ;
    else if (conf.m_stdDev) cout << "     " << left << setw(pad) << " " ;
    if( h && conf.m_integral ) cout << "I=" << left << setprecision(conf.m_precision) << setw(pad) << h->Integral() << " " ;
    else if (conf.m_integral) cout << "     " << left << setw(pad) << " " ;
  }
  cout << color << path << key->GetName() << conf.colorNorm << " " << key->GetCycle() ;

  if(dir2) 
    cout <<"; " << dir2->GetNkeys() << " Keys" << endl;
  else if ((conf.m_printTree || !conf.m_tree_code.IsNull())&& className=="TTree"){    
    TTree* t = dynamic_cast<TTree*> (obj);
    cout << "\n" << endl;
    if(conf.m_printTree){
      if(conf.m_printTree_topOnly) t->Print("toponly");
      else if( !conf.m_printTree_verbose) t->Print();
      else printTreeVerbose( t );    
    }
    if( !conf.m_tree_code.IsNull() ) printTreeCode( t );
    cout << "\n\n" << endl;
  }
  else cout << endl;
}

void printKeyBasic(TDirectory* dir, TKey* key, TObject* obj){
  TString className =  key->GetClassName();
  TString color;
  TString path = dir->GetPath();
  if(dir!=f) path+="/";
  conf.getColor( color, className);
  TDirectory* dir2 = dir->GetDirectory(key->GetName());
  if(dir2) cout << color << path << key->GetName() << conf.colorNorm << " " << dir2->GetNkeys() << " Keys" << endl;
  else     cout << color << path << key->GetName() << conf.colorNorm << endl;
  if(className == "TTree" && conf.m_printTree){
    TTree* t = dynamic_cast<TTree*> (obj);
    if(conf.m_printTree_topOnly) t->Print("toponly");
    else t->Print();
    cout << "\n\n" << endl;
  }
}

void printKeyClassic(const TString &name, const TString &className){
  TString color;
  conf.getColor( color, className);
  cout << color << name << conf.colorNorm ;
}

void printTreeVerbose( TTree* t , TBranch* parent, int level) {
  if( level ) cout << endl;
  if( parent ) {
    TIter itr(parent->GetListOfBranches());
    while( TBranch* br = dynamic_cast<TBranch*> (itr()) ) {
      printBranch( t, br, level );
      if(br->GetListOfBranches()->GetEntries()) printTreeVerbose( t, br, level+1 );
    }
  }
  else {
    TIter itr(t->GetListOfBranches());
    while( TBranch* br = dynamic_cast<TBranch*> (itr()) ) {
      printBranch( t, br, level );
      if(br->GetListOfBranches()->GetEntries()) {
	printTreeVerbose( t, br, level+1 );
      }
    }
  }
}

void printBranch( TTree* t, TBranch* b, int level){
  cout << "branch: " << left << setw(50) << b->GetName() << " Type " << getBranchType(b) << endl;
}

void printTreeCode( TTree* t, TBranch* parent, int level){
  if( level ) cout << endl;
  if( parent ) {
    TIter itr(parent->GetListOfBranches());
    while( TBranch* br = dynamic_cast<TBranch*> (itr()) ) {
      printBranchCode( t, br, level );
      if(br->GetListOfBranches()->GetEntries()) printTreeCode( t, br, level+1 );
    }
  }
  else {
    TIter itr(t->GetListOfBranches());
    while( TBranch* br = dynamic_cast<TBranch*> (itr()) ) {
      printBranchCode( t, br, level );
      if(br->GetListOfBranches()->GetEntries()) {
	printTreeCode( t, br, level+1 );
      }
    }
  }
}

void printBranchCode( TTree* t, TBranch* b, int level){
  TString name = b->GetName();
  TString type = getBranchType(b);
  type.ReplaceAll("(","");
  type.ReplaceAll(")","");
  if( type == "a character string terminated by the 0 character" ) type = "const char*";
  
  for(vector<TString>::const_iterator itr = conf.m_tree_code_v.begin(); itr != conf.m_tree_code_v.end(); ++itr){
    TString line = (*itr);
    line.ReplaceAll("\\t",type);
    line.ReplaceAll("\\b",name);
    cout << line << endl;
  }
}

TString getBranchType( TBranch* b ){
  TString c_name = b->GetClassName() ;
  TString c_title = b->GetTitle();
  TString p1, p2;
  splitString(c_title,p1,p2,'/');
  c_title = p2;
  //From class TTree
  //   - C : a character string terminated by the 0 character
  //     - B : an 8 bit signed integer (Char_t)
  //     - b : an 8 bit unsigned integer (UChar_t)
  //     - S : a 16 bit signed integer (Short_t)
  //     - s : a 16 bit unsigned integer (UShort_t)
  //     - I : a 32 bit signed integer (Int_t)
  //     - i : a 32 bit unsigned integer (UInt_t)
  //     - F : a 32 bit floating point (Float_t)
  //     - D : a 64 bit floating point (Double_t)
  //     - L : a 64 bit signed integer (Long64_t)
  //     - l : a 64 bit unsigned integer (ULong64_t)
  //     - O : a boolean (Bool_t)
  if( c_title == "C" ) c_title = "a character string terminated by the 0 character";
  if( c_title == "B" ) c_title = "(Char_t)"; 
  if( c_title == "b" ) c_title = "(UChar_t)"; 
  if( c_title == "S" ) c_title = "(Short_t)"; 
  if( c_title == "s" ) c_title = "(UShort_t)"; 
  if( c_title == "I" ) c_title = "(Int_t)"; 
  if( c_title == "i" ) c_title = "(UInt_t)"; 
  if( c_title == "F" ) c_title = "(Float_t)"; 
  if( c_title == "D" ) c_title = "(Double_t)"; 
  if( c_title == "L" ) c_title = "(Long64_t)"; 
  if( c_title == "l" ) c_title = "(ULong64_t)"; 
  if( c_title == "O" ) c_title = "(Bool_t)"; 

  if(c_title.Length()) return c_title;
  else {
    //    c_name.ReplaceAll("vector<","std::vector<");
    if (c_name.Contains("vector<")) c_name+="*";
    return c_name;
  }

}


TString MakeReadable(long &nBytes){
  TString str = "B";
  if(nBytes/1024 < 1) { return str; }
  str = "K";
  nBytes/=1024;
  if(nBytes/1024 < 1) { return str; }
  str = "M";
  nBytes/=1024;
  if(nBytes/1024 < 1) { return str; }
  str = "G";
  nBytes/=1024;
  return str;
}

void mapDirectory(TDirectory* dir, int level){
  TIter itr(dir->GetListOfKeys());
  
  while( TKey* key = dynamic_cast<TKey*> (itr())) {
    TDirectory* dir1 = 0;
    if(((TString) key->GetClassName())=="TDirectoryFile") dir1 = dir->GetDirectory(key->GetName());
    if(dir1){
      printPad(level);
      cout << key->GetName() << "  -----> " << endl;
      mapDirectory( dir1, level+1 );
      delete dir1;
    }
    else continue;
  }
}

TDirectory* listDirectoryFromPath(TString& str, TString &bashPath){

  if(str.BeginsWith("'") || str.EndsWith("'")) str.ReplaceAll("'","");

  TDirectory* dirBase = f;
  listDirectoryFromPathWild(dirBase, str, bashPath);

  return 0;
}

void listDirectoryFromPathWild(TDirectory* dirMother, TString& str, TString &bashPath, int level){

  stringResiduals strRes(str, conf.m_initialPath_regex);
  //strRes.print();
  TDirectory* dir = dirMother;

  if(strRes.firstPart!=""){//Jump To Directory if given 
    dir = dirMother->GetDirectory(strRes.firstPart);
    if(!dir && level==0) { cerr << "ERROR" << f->GetName() << " Has No Directory " << strRes.firstPart << endl; return; }
    else if (!dir) return;
  }
  else if(dir->GetMotherDir()) dirMother = dir->GetMotherDir();

  //if(strRes.secondPart=="" && ( (strRes.isDirectory && !strRes.hasWildCard) || (!strRes.isDirectory && strRes.hasWildCard) )){
  if(strRes.secondPart=="" && ( (strRes.isDirectory && !strRes.hasWildCard) || (!strRes.isDirectory && strRes.hasWildCard) || (!strRes.wildcard.IsNull() && !strRes.hasWildCard) )){
    if(!strRes.isDirectory ) //Begin main program w/ pattern on TObjects
      listDirectory(dir, bashPath, strRes.wildcard);
    else if (strRes.isDirectory && !strRes.hasWildCard) //Begin main program w/o pattern on TObjects
      listDirectory(dir, bashPath);
    return;
  }
  
  TIter nextcd(dir->GetListOfKeys());
  while(TKey* key = dynamic_cast<TKey*> (nextcd()) ){
    TString className = key->GetClassName();
    if(className!="TDirectory" && className!="TDirectoryFile") continue;
    TDirectory* dir1 = dir->GetDirectory(key->GetName());
    if(!dir1) continue;
    if(!dir1->InheritsFrom("TDirectory") ) continue;
    TString dirName = dir1->GetName();
    TPRegexp pattern = (TPRegexp) strRes.wildcard;
    if(dirName(pattern) == dirName) {
      if (strRes.secondPart=="") listDirectory(dir1, bashPath); 
      else listDirectoryFromPathWild(dir1, strRes.secondPart, bashPath, 1); 
    }    
    delete dir1;
  }
}

void openFile(TString & file){
  if(f){ closeFile(); }

  //strip :/initialRootPath if there
  //could cause problems??
  if(file.Contains(":/") && !file.Contains("rfio:/") && !file.Contains("://")){
    TString fName;
    splitString(file,fName,conf.m_initialRootPath,':');
    conf.m_initialRootPath = conf.m_initialRootPath(1,conf.m_initialRootPath.Length() -1);
    file = fName;
  }
  //Determine if file is a root file
  fstream fbin(file,ios::in | ios::binary);
  //Check header of file
  string line1;
  // line1 = (string) gSystem->GetFromPipe( (((TString)"head -c 4 ")+"\""+file+"\"").Data());
  char buff[4];
  fbin.read(buff,4);
  line1=buff;
  //  getline(fbin, line1);
  //See if file is likely on castor or webserver
  bool localfile = fbin.is_open();
  //Determine if file is castor file
  bool castorFile = false;
  if(!localfile && (file.Contains("/castor") || file.Contains("http://") || file.Contains("https://")) ) castorFile = true;

  //  if(line1.find("root")==string::npos) return;
  if(strncmp(line1.c_str(), "root", 4) && !castorFile) { fbin.close(); return; }

  bool binary = false;
  if(localfile) {
    int byte;
    while( (byte = fbin.get()) != EOF ){
      if(byte > 127) { binary = true; break; }
    }
  }
  if(!binary && localfile) { fbin.close(); return; }
  fbin.close();
  f = TFile::Open(file); 
  if(f) { 
    if(conf.m_scriptName!="")createScript();
    conf.rootOpens++;// 

    if(!conf.quiet) cout << "Opening : " << file << " Total Opens: " << conf.rootOpens << endl;
    conf.m_totKeys_glob += conf.m_totKeys;
    conf.m_totBins_glob += conf.m_totBins;
    conf.m_totEntries_glob += conf.m_totEntries;
    conf.m_totBytes_glob += conf.m_totBytes;
    conf.m_totKeys = 0;
    conf.m_totBins = 0;
    conf.m_totEntries = 0;
    conf.m_totBytes = 0;
  }
}

void closeFile(){
  //  cout << "Closing File " << endl;
  if(conf.m_stats) cout << f->GetPath() << conf.m_initialRootPath << " : " << "Total Keys:  " << conf.m_totKeys << "  Total Bytes:  " 
		   << conf.m_totBytes << "  Total Bins:  " << conf.m_totBins << " Total Entries: " << conf.m_totEntries << endl;  
  if(f){
    f->Close(); conf.rootCloses++;
    if(conf.m_scriptName != "") closeScript();
    if(!conf.quiet) cout << conf.rootOpens << " And Closes " << conf.rootCloses << endl;
    delete f;
    f = 0;
  }
}

void createScript(){
  if(!conf.quiet) cout << "Opening script file" << endl;
  setStyle();
  TString name=conf.m_scriptName+".C";
  conf.m_histNumber = 1;
  conf.m_can = new TCanvas();
  //Remove any spaces from histName
  conf.m_histName.ReplaceAll(" ","");
  parseHistName(conf.m_histName);
  //Now we know how many hists to per page figure out conf.m_canCol/conf.m_canRow
  conf.m_canPerPage = conf.m_histNames.size();
  if(conf.m_canPerPage==1) { conf.m_canCol = 1; conf.m_canRow =1;}
  else if(conf.m_canPerPage==2) { conf.m_canCol = 1; conf.m_canRow =2;}
  else if(conf.m_canPerPage==3 || conf.m_canPerPage==4) { conf.m_canCol = 2; conf.m_canRow =2;}
  else if(conf.m_canPerPage < 7) { conf.m_canCol = 2; conf.m_canRow = 3;}
  else  { conf.m_canCol = 2; conf.m_canRow = 4;}

  conf.m_script_ps = new TPostScript( conf.m_scriptName+".ps" , 112);

}

void closeScript(){
  if(conf.m_obj.size()) fillCanvas(true); 
  else { 
    conf.m_script_ps->On();
    conf.m_script_ps->Close();
    delete conf.m_can;
  }
  #ifndef WIN32
  if(conf.m_scriptName.EndsWith(".pdf")){
    TString fileName = conf.m_scriptName+".ps";
    TString command;
    if( gSystem->Exec("which ps2pdf") == 0) command = "ps2pdf ";
    else if( gSystem->Exec("which pstopdf") ==0) command = "pstopdf ";
    else {
      cerr << "ERROR    cannot convert ps to pdf cannot find command 'pstopdf' or 'convert' if such another command exists send email to griffith@cern.ch" << endl;
      return;
    }
    command += conf.m_scriptName+".ps "+conf.m_scriptName;
    gSystem->Exec(command.Data());
    gSystem->Exec( ((TString) "rm " + conf.m_scriptName+".ps").Data() );
    if(command.BeginsWith("pstopdf") ) gSystem->Exec( ((TString) "mv " + conf.m_scriptName+".pdf " + conf.m_scriptName).Data() );
    //    else gSystem->Exec( ((TString) "mv " + conf.m_scriptName+".pdf " + conf.m_scriptName).Data() );
  }
  #endif 

}

void parseHistName(const TString &name){
  TString hName1;
  TString res;
  splitString(name, hName1, res, ',');

  int statLevel = 1111;
  if(hName1.Contains("<") && hName1.Contains(">")){
    statLevel = ((TString) hName1(hName1.Index("<") + 1, hName1.Sizeof() - 1 - hName1.Index("<")) ).Atoi();
    hName1 = hName1(0, hName1.Index("<") );
  }

  int logLevel = 0;
  if(hName1.EndsWith("Logy")||hName1.EndsWith("LogY")||hName1.EndsWith("logy")){
    logLevel = 1;
    hName1("Logy") = "";
    hName1("LogY") = "";
    hName1("logy") = "";
  }
  if(hName1.EndsWith("Logz")||hName1.EndsWith("LogZ")||hName1.EndsWith("logz")){
    logLevel = 2;
    hName1("Logz") = "";
    hName1("LogZ") = "";
    hName1("logz") = "";
  }

  bool foundname = true;
  while(foundname) { //Allow histName = '.*,.*,.*,.*' to work by changing 2nd .* to .*.*
    std::map<TString, histDisplayOpt>::const_iterator itr = conf.m_histNames.find(hName1);
    if (itr == conf.m_histNames.end()) foundname = false;
    else hName1 += ".*";//ensure unique entries to the map
  }
  conf.m_histNames.insert(make_pair(hName1,histDisplayOpt(logLevel, statLevel)));

  if(res=="" || conf.m_histNames.size()==8) return;
  else parseHistName(res);
}

bool isInHistList(const TString &name){
  for(map<TString, histDisplayOpt>::const_iterator itr = conf.m_histNames.begin(); itr != conf.m_histNames.end(); ++itr){
    TPRegexp reg = (TPRegexp) (*itr).first;
    TString name_ = name(reg);
    if(name_==name) {
      conf.m_histDisplayOpt.push_back((*itr).second);
      return true;
    }
  }
  return false;
}


void addHistToScript(TString path, TObject* obj){

  conf.m_histNumber++;

  //Add pdf page for hist
  if(conf.m_can_paths.IsNull()) conf.m_can_paths = path;
  else conf.m_can_paths = conf.m_can_paths + "," + path;
  TH1* obj_c = dynamic_cast<TH1*> (obj);
  TCanvas* obj_c_can = dynamic_cast<TCanvas*> (obj);
  if(!obj_c && !obj_c_can) return;
  if(obj_c){
    TH1* h = (TH1*) obj_c->Clone();
    h->SetDirectory(0);//Separate this object from rest of program delete in fillCanvas
    conf.m_obj.push_back( h );
  }
  if(obj_c_can){
    TCanvas* c = (TCanvas*) obj_c_can->Clone();
    conf.m_obj.push_back(c);
  }
  if(conf.m_obj.size() == conf.m_canPerPage) {fillCanvas();}

}

void fillCanvas(bool lastCan /*=false*/){
  conf.m_can->Clear();
  conf.m_script_ps->On();
  vector<TString> paths = parseList( conf.m_can_paths, "," );
  conf.m_can->Divide(conf.m_canCol, conf.m_canRow);
  TString hName1;
  if(conf.m_obj.size()) hName1 = conf.m_obj.at(0)->GetName();
  conf.m_can->SetTitle(paths[0]+conf.m_histName);

  for(unsigned i = 0; i != conf.m_obj.size(); ++i){

    conf.m_can->cd(i+1);
    TH1* h = dynamic_cast<TH1*> (conf.m_obj.at(i));
    TCanvas* c = dynamic_cast<TCanvas*> (conf.m_obj.at(i));
    int logLevel = 0;
    if(i<conf.m_histDisplayOpt.size()) logLevel=conf.m_histDisplayOpt.at(i).loglevel;
    if(logLevel==1) gPad->SetLogy();
    if(logLevel==2) gPad->SetLogz();

    if(!h && !c) continue;
    int statLevel = 1111;
    statLevel = conf.m_histDisplayOpt.at(i).statOpt;
    gStyle->SetOptStat(statLevel);
    if(h){
      h->SetTitle(paths[i]+"/"+h->GetTitle());
      //    h->Draw();
      gPad->GetListOfPrimitives()->Add(h);
    }
    else {
      //   cout << "Adding to Pad " << gPad << " Can " << conf.m_can->cd(i+1) << endl;
      //      gPad->GetListOfPrimitives()->Add(c);      


      gPad->SetTopMargin( c->GetTopMargin() );
      gPad->SetBottomMargin( c->GetBottomMargin() );
      gPad->SetRightMargin( c->GetRightMargin() );
      gPad->SetLeftMargin( c->GetLeftMargin() );

      TIter iter(c->GetListOfPrimitives());
      int hists = 0;
      while( TObject* obj = dynamic_cast<TObject*> (iter()) ){
	TH1* h = dynamic_cast<TH1*> (obj);
	if(h){
	  if(hists==0) {
	    h->SetTitle(paths[i]+"/"+h->GetTitle());
	    h->Draw();
	  }
	  else h->Draw("sames");
	  hists++;
	}
	else {
	  THStack* hs = dynamic_cast<THStack*> (obj);
	  TGraph* g = dynamic_cast<TGraph*> (obj);
	  TGraphErrors* ge = dynamic_cast<TGraphErrors*> (obj);
	  if(hs) {
	    if(hists==0) {
	      hs->SetTitle(paths[i]+"/"+hs->GetTitle());
	      hs->Draw();
	    }
	    else hs->Draw("sames");
	    hists++;
	  }	  
	  else if(g) {
	    if(hists==0) {
	      g->SetTitle(paths[i]+"/"+g->GetTitle());
	      g->Draw("AP");
	    }
	    else g->Draw("Psame");
	    hists++;
	  }
	  else if(ge) {
	    if(hists==0) {
	      ge->SetTitle(paths[i]+"/"+ge->GetTitle());
	      ge->Draw("AP");
	    }
	    else ge->Draw("Psame");
	    hists++;
	  }
	  else obj->Draw("same");
	}

	  
      }//objs of canvas
    }//not h
  }//loop over conf.m_obj

  conf.m_can->Update();
  if(lastCan) conf.m_script_ps->Close();

  conf.m_script_ps->NewPage();
  //Turn off file so that objs can be safely deleted
  conf.m_script_ps->Off();

  for(unsigned i = 0; i != conf.m_obj.size(); ++i) delete conf.m_obj.at(i);
  conf.m_obj.clear();
  conf.m_histDisplayOpt.clear();
  conf.m_can_paths = "";
}



void setStyle(){

  gStyle->SetPalette(1);
  //Use AtlasStyle
  Int_t icol=0;
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetStatColor(icol);

  // set the paper & margin sizes
  gStyle->SetPaperSize(20,26);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.12);

  // use large fonts
  //Int_t font=72;
  Int_t font=42;
  Double_t tsize=0.05;
  gStyle->SetTextFont(font);

  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font,"x");
  gStyle->SetTitleFont(font,"x");
  gStyle->SetLabelFont(font,"y");
  gStyle->SetTitleFont(font,"y");
  gStyle->SetLabelFont(font,"z");
  gStyle->SetTitleFont(font,"z");

  gStyle->SetLabelSize(tsize,"x");
  gStyle->SetTitleSize(tsize,"x");
  gStyle->SetLabelSize(tsize,"y");
  gStyle->SetTitleSize(tsize,"y");
  gStyle->SetLabelSize(tsize,"z");
  gStyle->SetTitleSize(tsize,"z");


  //use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth((Width_t) 2.);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  //do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
//   //gStyle->SetOptStat(1111);
//   if(!stats)  gStyle->SetOptStat(0);
//   else gStyle->SetOptStat(1111); 
 //gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gROOT->SetStyle("Plain");
}

