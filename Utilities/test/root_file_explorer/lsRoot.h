#ifndef LSROOT_H
#define LSROOT_H

#include "stringUtils.h"
#include "lsConfig.h"
#include "lsGlobals.h"

#include <TFile.h>
#include <TKey.h>
#include <TDatime.h>
#include <TString.h>
#include <TStyle.h>
#include <TPRegexp.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TPostScript.h>

#include <iomanip>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>

//My File 
TFile* f = 0;

void mergeFiles();
bool listDirectory(TDirectory* dir, TString bashPath, TString pattern = "");
void listDirectory(void* dp, TString bashPath, TString initialRootPath = "");
void mapDirectory(TDirectory* dir, int level);
TString MakeReadable(long &nBytes);

void printKeyCustom(TDirectory* dir, TKey* key, TObject* obj, TString path);
void printKey(TDirectory*, TKey* key, TObject*, TString);
void printKeyBasic(TDirectory*, TKey* key, TObject*);
void printKeyClassic(const TString &, const TString &);
void printTreeVerbose( TTree* t , TBranch* parent = 0, int level = 0);
void printBranch(TTree* t, TBranch* b, int level=0);
void printTreeCode( TTree* t, TBranch* parent = 0, int level = 0);

TString getBranchType( TBranch* b );
void printBranchCode( TTree* t, TBranch* b, int level=0);

void openFile(TString & file);
void closeFile();

TDirectory* listDirectoryFromPath(TString& str, TString &bashPath);
void listDirectoryFromPathWild(TDirectory* dir, TString& str, TString &bashPath, int level = 0);

void setStyle();
void createScript();
void parseHistName(const TString &name);
bool isInHistList(const TString &name);
void addHistToScript(TString path, TObject* obj);
void fillCanvas(bool lastCan = false);


void drawAsci(TH1* h, asciDrawOptions* opt=conf.m_drawOpt, bool originalCall = true);
void drawAsci2(TH2* h, asciDrawOptions* opt=conf.m_drawOpt2, bool originalCall = true);
void drawGraph(TGraph* g);

void closeScript();

void Initialize();

/* void setColors(); */
/* void setInitialPrefs(); */

#endif
