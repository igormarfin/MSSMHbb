#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TMarker.h>
#include <TColor.h>
#include <TFractionFitter.h>
#include <THStack.h>
#include <TMinuit.h>
#include <TMath.h>

/// systematics read-out
#include <map>
#include <TString.h>
#include <TObjArray.h>
#include <TSystem.h>

#include "Analysis/Utilities/interface/getHbbCfg.h"

TCanvas* canvas;

TH2F* mergeSignal(TFile* fa,char* hname,const char* newname) {
  TH2F* theData[2];
  if (fa == NULL) {
    std::cout << "mergeSignal: no hfile opened" << std::endl;
    return NULL;
  }

//   // weights
//   double tWeight[2];
//   tWeight[0] = (524.90+265.75)*.94*0.94 + 193.58+251.12+453.33;
//   tWeight[1] = 246.53+732.73;
//   double sumWeight = tWeight[0]+tWeight[1];
//   for (int ii=0; ii<2; ++ii) {
//     tWeight[ii] = tWeight[ii] / sumWeight;
//   }
//   std::cout << "Trig weights: " << tWeight[0] << " " << tWeight[1] << std::endl;
  



  // only Trig0 is needed
  for (int ii=0; ii<1; ++ii) {
    std::cout << "mergeSignal: pick histogram " << Form("%sTrig%d",hname,ii) << std::endl;
    theData[ii] = (TH2F*) fa->Get(Form("%sTrig%d",hname,ii));
    if (theData[ii] == NULL) {
      std::cout << "mergeSignal: histogram " << Form("%sTrig%d",hname,ii) << " not found" << std::endl;
      return 0;
    }
    //theData[ii]->Draw();
  }
  std::string theName( hname );
  // strip off any path in histogram name
  long unsigned int nsl = theName.find_last_of('/');
  if (nsl != theName.npos) {
    theName = theName.substr(nsl+1,theName.size()-nsl-1);
  }
  std::string theTitle( theData[0]->GetTitle() );
  TH2F* mergedData = new TH2F( *theData[0] );

  //mergedData->Add(theData[0],theData[1],tWeight[0],tWeight[1]);

  if (newname == 0) {
    mergedData->SetName( Form("%s%s",theName.c_str(),"TrigMerged") );
  } else {
    mergedData->SetName( newname );
  }
  mergedData->SetTitle( mergedData->GetName() );
  std::cout << "mergeSignal: hist has name " << mergedData->GetName() << std::endl;
  //mergedData->Draw();
  //canvas->Print(Form("%s.png",newname));
  return mergedData;
}


TH2F* getTrigsAndMerge(TFile* fa,char* hname,const int nTCombData,string tCombData[]) {
  TH2F* theData[nTCombData];
  if (fa == NULL) {
    std::cout << "getTrigsAndMerge: no hfile opened" << std::endl;
    return NULL;
  }
  std::cout << "getTrigsAndMerge: hname=" << hname << std::endl;
  for (int iTCombData=0; iTCombData<nTCombData; ++iTCombData) {
//     std::cout << "getTrigsAndMerge: iTCombData= " << iTCombData << " tCombData=" << tCombData[iTCombData] << std::endl;

    theData[iTCombData] =  (TH2F*) fa->Get(Form("%s%s",hname,tCombData[iTCombData].c_str()));
    if (theData[iTCombData] == NULL) {
      std::cout << "getTrigsAndMerge: histogram " << Form("%s%s",hname,tCombData[iTCombData].c_str()) << " not found" << std::endl;
      return 0;
    }
  }
  string theName( hname );
  // strip off any path in histogram name
  long unsigned int nsl = theName.find_last_of('/');
  if (nsl != theName.npos) {
    theName = theName.substr(nsl+1,theName.size()-nsl-1);
  }
  string theTitle( theData[0]->GetTitle() );
  TH2F* mergedData = new TH2F( *theData[0] );
  for (int iTCombData=0; iTCombData<nTCombData; ++iTCombData) {
    theName += tCombData[iTCombData];
    if (iTCombData != 0) theTitle += tCombData[iTCombData];
  }
  mergedData->SetName( theName.c_str() );
  for (int iTCombData=1; iTCombData<nTCombData; ++iTCombData) {
    mergedData->Add( theData[iTCombData] );
  }
  std::cout << "getTrigsAndMerge: new name: " << mergedData->GetName() << std::endl;
  return mergedData;
}





TH2F* mergeBackground(TFile* fa,const char* pattern,const char * theSbtag,const char* newname,const string scenario) {
  std::cout << "mergeBackground: newname = " << newname << std::endl;
  const int nTrig=3;
  const int nfc=3;
  const int ncateg=3;
  string sfc[nfc] = { "q", "c", "b" };
  TH2F* theData[nTrig][nfc][ncateg];
  if (fa == NULL) {
    std::cout << "mergeBackground: no hfile opened" << std::endl;
    return NULL;
  }

  // no weights
  for (int iTrig=0; iTrig<nTrig; ++iTrig) {
    for (int ifc=0; ifc<nfc; ++ifc) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	//string sname( Form(Form("%sTrig%d",pattern,iTrig),sfc[ifc].c_str(),theSbtag,icateg,icateg) );
	string sname( Form("bgPredict/MassBTagPred_%s_%s_Cat%dTpat%dTrig%d",
			   sfc[ifc].c_str(),theSbtag,icateg,icateg,iTrig) );
	theData[iTrig][ifc][icateg] = (TH2F*) fa->Get( sname.c_str());
	if (theData[iTrig][ifc][icateg] == NULL) {
	  std::cout << "mergeBackground: histogram " << sname << " not found" << std::endl;
	  return 0;
	}
    //theData[ii]->Draw();
      }
    }
  }

  TH2F* mergedData = new TH2F( *theData[0][0][0] );
  mergedData->SetName( newname );
  mergedData->SetTitle( newname );
  for (int iTrig=0; iTrig<nTrig; ++iTrig) {
    for (int ifc=0; ifc<nfc; ++ifc) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	if (! ((iTrig == 0) && (ifc == 0) && (icateg == 0) )) {
	  mergedData->Add( theData[iTrig][ifc][icateg] );
	  std::cout << "mergeBackground: adding " << theData[iTrig][ifc][icateg]->GetName() << std::endl;
	}
      }
    }
  }  
  std::cout << "mergeBackground: newname(4) = " << newname << std::endl;

  // cross check total weight
  double tw = 0;
  for (int iTrig=0; iTrig<nTrig; ++iTrig) {
    for (int ifc=0; ifc<nfc; ++ifc) {
      for (int icateg=0; icateg<ncateg; ++icateg) {
	tw += theData[iTrig][ifc][icateg]->GetSumOfWeights();
      }
    }
  }
  std::cout << newname << " Sum = " << tw << "  mergedData = " << mergedData->GetSumOfWeights() << std::endl;
  std::cout << mergedData->GetName() << std::endl;
  
  // final normalization (from control region)
  double wScale = 1;
  string sb( theSbtag );
  if (scenario == "LowMass2011") {
    {
      if ( sb == "CSVT" ) {
	wScale = 1.07817;
      } else if ( sb == "TCHPT" ) {
	wScale = 0.970401;
      }
    }
  } else if (scenario == "MediumMass2011") {
    {
      if ( sb == "CSVT" ) {
	wScale = 0.9505;
      } else if ( sb == "TCHPT" ) {
	wScale = 0.874897;
      }
    }
  } else {
    std::cout << "merge background: bad scenario " << scenario << std::endl;
    return 0;
  }
  std::cout << "Use expected background scale factor of " << wScale << " for " << theSbtag 
	    << " scenario " << scenario << std::endl;
  mergedData->Scale( wScale );
  std::cout << newname << " After final normalization: " << mergedData->GetSumOfWeights() << std::endl;

  return mergedData;
}


class templateId {
public:
  int flav;
  int categ;
  templateId(int theFlav,int theCateg) : flav(theFlav), categ(theCateg) {}
  string name() {
    string sFlav[3] = {"Q", "C", "B"};
    string sLegend("bbb");
    sLegend.replace(categ,1,sFlav[flav]);
    return sLegend;
  }
};




void packTemplatesMass(const int iMass) {
  // This Root macro is for the purpose of providing input for expected limit computation
  // It packs predicted background and a set of signal samples into a root file per Higgs mass
  // 
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  canvas = new TCanvas ("cg1","PadX",10,10,800,600);
  gStyle->SetPadColor(0);
  canvas->SetFillColor(0);

//   const int nbtag = 4;
  const int nbtag = 1;
  //const std::string sbtag[nbtag] = { "TCHPT", "TCHP6", "CSVT", "SSVHPT" };
  const std::string sbtag[nbtag] = { "CSVT" };

  const int nfc=3;
  const int ncateg=3;
  string sfc[nfc] = { "q", "c", "b" };

  // this is for the combination of triggers in real data
  const int nTCombData = 4;
  std::string tCombData[nTCombData] = {"Trig0", "Trig1", "Trig2", "Trig3"};

  std::string L1L2Mode("Weight");
  std::string signalMode("PU_WEIGHTED-NEW");
//   std::string L1L2Mode("Cut");
//   std::string signalMode("CUT_BASED");

  //std::string scenario("LowMass2011");
  //std::string scenario("MediumMass2011");

  std::string scenario;
  bool useTemplateError;
  bool useNP;
  if (  getHbbCfg(scenario,useTemplateError,useNP) != 0 ) return;

  string IgorVersion("V6");
#include "Analysis/Utilities/interface/HbbMass.h"

  if (iMass >= nSignal) {
    std::cout << "Bad iMass=" << iMass << std::endl;
    return;
  }
//   const int nSignal=7;
//   int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };
  // int signalMass[nSignal] = { 90, 100, 120, 130, 140, 160, 180, 200, 250, 350 }
  //  there are also : 450, 500, 600, 700, 800, 900, 1000
//   double efficiency[nSignal] = { 0.0022081, 0.00324694, 0.00600146, 0.00918135,
// 				0.0138382, 0.0189684, 0.0206572 };
  double efficiency[nSignal][nbtag];
  double intLumi = 0;
  string IgorScen("");
  string spacer("");
  string SashaPath("");
  string IgorPath("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/SignalTemplates-Production");

  if (IgorVersion == "V4") {
    IgorPath.assign("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-4");
  } else if (IgorVersion == "V6") {
    IgorPath.assign("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/SignalTemplates-Production2");
  }

  if (scenario == "LowMass2011") {
    //intLumi = 2.66794; // in fb-1
    intLumi = 2.692643;  // with new method
    IgorScen.assign("low");
    spacer.assign("");
    //SashaPath.assign("Data-Run2011AB");
    SashaPath.assign("Data-Run2011AB/TripleBtagAnalysis_CR3_SF7");
  } else if (scenario == "MediumMass2011") {
    //intLumi = 3.99983; // in fb-1
    intLumi = 4.040802;
    IgorScen.assign("medium");
    spacer.assign("/MEDIUM");
    //SashaPath.assign("Data-Run2011AB-Medium");
    SashaPath.assign("Data-Run2011AB/TripleBtagAnalysis_CR3_SF7_med");
  } else {
    std::cout << "Bad scenario" << std::endl;
    return;
  }
  string signalHistPattern("massEvBtag/mjjEvBTag_%s");
  if (L1L2Mode == "Weight") {
    signalHistPattern.assign("massEvBtagTW/mjjEvBTagTW_%s");
  }

  double fScal[nSignal][nbtag];
  
  // systematics
//   const int nSyst = 3;
//   std::string systName[nSyst] = { "JES", "SFbc", "SFudsg" };
//   bool realDataNuisance[nSyst] = { false, true, true }; // indicate if relevant for real data

  string bbPurity("DataDriven"); //  "MC", "DataDrivenR", "None"
  //string bbPurity("None"); //  "MC", "DataDrivenR", "None"

  bool onlineBtagCorr = true;

  const int nSyst = 4;
  std::string systName[nSyst] = { "JES", "SFbc", "SFudsg", "JER" };
  bool realDataNuisance[nSyst] = { false, true, true, false }; // indicate if relevant for real data

  const int nUpDown = 2;
  std::string signalFile = "";


  // signal templates
  if (IgorVersion != "V3") {
    std::cout << "Using signal files " << IgorVersion << std::endl;
    signalFile.assign( Form("%s/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola%s/SF/job_1/TripleBtagAnalysisM-%d_%s.root",IgorPath.c_str(),signalMass[iMass],spacer.c_str(),signalMass[iMass],IgorScen.c_str() ) );
  } else {
    std::cout << "Using V3 signal files" << std::endl;
    signalFile.assign( Form("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-3/%s/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola%s/SF/job_1/TripleBtagAnalysisM-%d_%s.root",signalMode.c_str(),signalMass[iMass],spacer.c_str(),signalMass[iMass],IgorScen.c_str() ) );
  }
  std::string signalSystFiles[nSyst][nUpDown];
  //bool activeNuisance[nSyst] = { true, true, true };
  std::string upDownName[nUpDown] = { "Up", "Down" };
  TH2F* hSignalSyst[nSyst][nUpDown][nbtag];

  // output file
  TFile* hout = new TFile(Form("packedTemplates-M-%d.root",signalMass[iMass]),"recreate");
  hout->cd();
  TH2::AddDirectory(true);

  TFile* fSig = new TFile( signalFile.c_str() );
  if ( fSig == NULL ) {
    std::cout << "Could not open signal central file " << signalFile.c_str() << std::endl;
    return;
  } else {
    std::cout << "Open signal file " << signalFile.c_str() << std::endl;
  }

  TH2F* hSignalCentral[nbtag];

  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hSignalCentral[ibtag] = mergeSignal(fSig,Form(signalHistPattern.c_str(),sbtag[ibtag].c_str()),
					Form("bbH_%s",sbtag[ibtag].c_str()));
    // read the efficiency
    TH1F* histEffMerged = (TH1F*) fSig->Get(Form("TrigEff/EffMerged%s",sbtag[ibtag].c_str()));
    if ( histEffMerged == NULL) {
      std::cout << "Efficiency histo not found" << std::endl;
      return;
    }
    double newEff = histEffMerged->GetBinContent(1);
    std::cout << "Mass= " << signalMass[iMass] 
	      << " btag= " << sbtag[ibtag]
	      << " Efficiency = " << newEff << std::endl;
    efficiency[iMass][ibtag] = newEff;

    double normShould = 1000 * intLumi * efficiency[iMass][ibtag];
    double normIs = hSignalCentral[ibtag]->GetSumOfWeights();
    std::cout << hSignalCentral[ibtag]->GetName() << " TotalContents=" << hSignalCentral[ibtag]->GetSumOfWeights()
	      << std::endl;
    fScal[iMass][ibtag] = normShould / normIs;
    std::cout << "normShould = " << normShould << " normIs " << normIs
	      << " rescale by " << fScal[iMass][ibtag] << std::endl;
    hSignalCentral[ibtag]->Scale( fScal[iMass][ibtag] );
    hout->cd();
    hSignalCentral[ibtag]->Write();
    histEffMerged->Write();

    // create empty file just as marker
    ofstream markerFile;
    markerFile.open(Form("pack-%s-%s.txt",sbtag[ibtag].c_str(),scenario.c_str()),ios::app);
    markerFile << "Template for mass " << signalMass[iMass] << std::endl;
    markerFile.close();
  }

  // read the nominal cross section
  TH1F* histXSect = (TH1F*) fSig->Get("xsection/xsect");
  if ( histXSect == NULL) {
    std::cout << "xsection/xsect" << " not found" << std::endl;
    return;
  }
  histXSect->Write();

  for (int iSyst=0; iSyst<nSyst; ++iSyst) {
    for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
      if (IgorVersion != "V3") {
	signalSystFiles[iSyst][iUpDown] = Form( "%s/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola%s/%s_Sys%s/job_1/TripleBtagAnalysisM-%d_%s.root",
						IgorPath.c_str(),signalMass[iMass],spacer.c_str(),systName[iSyst].c_str(),
						upDownName[iUpDown].c_str(),signalMass[iMass],IgorScen.c_str());
      } else {	
	signalSystFiles[iSyst][iUpDown] = Form( "/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-3/%s/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola%s/%s_Sys%s/job_1/TripleBtagAnalysisM-%d_%s.root",
						signalMode.c_str(),signalMass[iMass],spacer.c_str(),systName[iSyst].c_str(),
						upDownName[iUpDown].c_str(),signalMass[iMass],IgorScen.c_str());
      }
      std::cout << "Signal systematics file " << signalSystFiles[iSyst][iUpDown] << std::endl;
      TFile* fSigSys = new TFile( signalSystFiles[iSyst][iUpDown].c_str() );
      if ( fSigSys == NULL ) {
	std::cout << "Could not open signal syst file " << signalSystFiles[iSyst][iUpDown].c_str() << std::endl;
	return;
      }
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	hSignalSyst[iSyst][iUpDown][ibtag] 
	  = mergeSignal(fSigSys,Form(signalHistPattern.c_str(),sbtag[ibtag].c_str()),
			Form("bbH_%s_%s_%s",systName[iSyst].c_str(),
			     upDownName[iUpDown].c_str(),sbtag[ibtag].c_str()));
	std::cout << "The merged hist has name " << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << std::endl;
	std::cout << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << " TotalContents=" << hSignalSyst[iSyst][iUpDown][ibtag]->GetSumOfWeights()
		  << std::endl;
	
	hSignalSyst[iSyst][iUpDown][ibtag]->Scale( fScal[iMass][ibtag] );
	hout->cd();
	hSignalSyst[iSyst][iUpDown][ibtag]->Write();
      }
      fSigSys->Close();
    }
  }

  // real data
  std::string backgroundFile( Form("/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/v1/%s/TripleBtagAnalysis_SF/TripleBtagAnalysis.root",SashaPath.c_str()) );
  std::cout << "Background central file : " << backgroundFile << std::endl;
  TFile* fBac = new TFile( backgroundFile.c_str() );
  if ( fBac == NULL ) {
    std::cout << "Could not open background central file " << signalFile.c_str() << std::endl;
    return;
  }

  // hist-to-be-fitted
  TH2F* mjjEbtdata[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    mjjEbtdata[ibtag] = getTrigsAndMerge(fBac,Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()),nTCombData,tCombData);

    if (mjjEbtdata[ibtag] == NULL) {
      std::cout << "Histogram not found: " << Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()) << std::endl;
      return;
    }
    // rename
    mjjEbtdata[ibtag]->SetName( Form("Data_%s",sbtag[ibtag].c_str() ) );
    hout->cd();
    mjjEbtdata[ibtag]->Write();
  }

  TH2F* hBackgroundCentral[nbtag][ncateg][nfc];
  TH2F* hBackgroundCentralError[nbtag][ncateg][nfc];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      int theTpat;
      if (onlineBtagCorr) {
	theTpat = 3;
      } else {
	theTpat = icateg;
      }

      for (int ifc=0; ifc<nfc; ++ifc) {
	string templateCore("massBTagTemplatesCld/MassBTagTemplateCld");
	if (bbPurity == "None") {
	   templateCore.assign("massBTagTemplatesCld/MassBTagTemplateUncld");
	}
	string hbSystName( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
				sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,theTpat) );
	hBackgroundCentral[ibtag][icateg][ifc] = getTrigsAndMerge(fBac,Form("%s_%s_%s_Cat%dTpat%d",
									   templateCore.c_str(),
									   sfc[ifc].c_str(),
									   sbtag[ibtag].c_str(),icateg,theTpat),nTCombData,tCombData);

	if ( hBackgroundCentral[ibtag][icateg][ifc] == NULL ) {
	  std::cout << "Hist not found: " << hbSystName << std::endl;
	  return;
	}
	// rename
	templateId tName(ifc,icateg);
	hBackgroundCentral[ibtag][icateg][ifc]->SetName( Form("%s_%s",tName.name().c_str(),sbtag[ibtag].c_str()) );
	// read the template errors
	templateCore.assign("errorMassBTagTemplates/ErrorMassBTagTemplate");
	string hbSystNameError( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
				sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,theTpat) );
	hBackgroundCentralError[ibtag][icateg][ifc] 
	  = getTrigsAndMerge(fBac,Form("%s_%s_%s_Cat%dTpat%d",
				       templateCore.c_str(),
				       sfc[ifc].c_str(),
				       sbtag[ibtag].c_str(),icateg,theTpat),nTCombData,tCombData);
	if ( hBackgroundCentralError[ibtag][icateg][ifc] == NULL ) {
	  std::cout << "Hist not found: " << hbSystNameError << std::endl;
	  return;
	}
	if (useTemplateError) {
	  // add the template error
	  std::cout << " ==== Adding Btag Errors ==== " << hBackgroundCentral[ibtag][icateg][ifc]->GetName() << std::endl;
	  for (int ibinx=1; ibinx<= (hBackgroundCentral[ibtag][icateg][ifc]->GetXaxis()->GetNbins()); ++ibinx) {
	    for (int ibiny=1; ibiny<= (hBackgroundCentral[ibtag][icateg][ifc]->GetYaxis()->GetNbins()); ++ibiny) {
	      float oldError = hBackgroundCentral[ibtag][icateg][ifc]->GetBinError(ibinx,ibiny);
	      float addError = hBackgroundCentralError[ibtag][icateg][ifc]->GetBinContent(ibinx,ibiny);
	      float newError = sqrt( oldError * oldError + addError * addError );
	      hBackgroundCentral[ibtag][icateg][ifc]->SetBinError(ibinx,ibiny,newError);
	    }
	  }
	}
	hout->cd();
	hBackgroundCentral[ibtag][icateg][ifc]->Write();
      }
    }
  }

  std::string backgroundSystFiles[nSyst][nUpDown];
  std::string systNameSasha[nSyst] = { "JES", "SFbc", "SFq" };
  std::string upDownNameSasha[nUpDown] = { "plus2", "minus2" };
  TH2F* hBackgroundSyst[nSyst][nUpDown][nbtag][ncateg][nfc];
  TH2F* hBackgroundSystError[nSyst][nUpDown][nbtag][ncateg][nfc];

  // for nuisances like JEC, the up/down templates are just copies of the central templates
  for (int iSyst=0; iSyst<nSyst; ++iSyst) {
    if (! realDataNuisance[iSyst]) {
      std::cout << "Non-real data relevant nuisance: " << systName[iSyst].c_str() << std::endl;
      for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	  for (int icateg=0; icateg<ncateg; ++icateg) {
	    for (int ifc=0; ifc<nfc; ++ifc) {
	      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]
		= new TH2F( *hBackgroundCentral[ibtag][icateg][ifc] );
	      // rename
	      templateId tName(ifc,icateg);
	      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]
		->SetName( Form("%s_%s_%s_%s",
				tName.name().c_str(),systName[iSyst].c_str(),
				upDownName[iUpDown].c_str(),sbtag[ibtag].c_str()) );
	      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->Write();
	    }
	  }
	}
      }
    }
  }
  fBac->Close();


  for (int iSyst=0; iSyst<nSyst; ++iSyst) {
    if (realDataNuisance[iSyst]) {
      std::cout << "Real data relevant nuisance: " << systName[iSyst].c_str() << std::endl;
      for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	backgroundSystFiles[iSyst][iUpDown] = Form( "/afs/naf.desy.de/user/r/rmankel/scratch/HbbPat/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/v1/%s/TripleBtagAnalysis_%s%s/TripleBtagAnalysis.root",SashaPath.c_str(),systNameSasha[iSyst].c_str(),upDownNameSasha[iUpDown].c_str() );
	TFile* fBacSys = new TFile( backgroundSystFiles[iSyst][iUpDown].c_str() );
	std::cout << "Background systematics file " << backgroundSystFiles[iSyst][iUpDown] << std::endl;
	if ( fBacSys == NULL ) {
	  std::cout << "Could not open background syst file " << backgroundSystFiles[iSyst][iUpDown] << std::endl;
	  return;
	}
	for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	  for (int icateg=0; icateg<ncateg; ++icateg) {
	    int theTpat;
	    if (onlineBtagCorr) {
	      theTpat = 3;
	    } else {
	      theTpat = icateg;
	    }

	    for (int ifc=0; ifc<nfc; ++ifc) {
	      string templateCore("massBTagTemplatesCld/MassBTagTemplateCld");
	      if (bbPurity == "None") {
		templateCore.assign("massBTagTemplatesCld/MassBTagTemplateUncld");
	      }
	      string hbSystName( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
				      sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,theTpat) );
	      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] 
		= getTrigsAndMerge(fBacSys,Form("%s_%s_%s_Cat%dTpat%d",
						templateCore.c_str(),
						sfc[ifc].c_str(),
						sbtag[ibtag].c_str(),icateg,theTpat),nTCombData,tCombData);
	      if ( hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc] == NULL ) {
		std::cout << "Hist not found: " << hbSystName << std::endl;
		return;
	      }
	      // rename
	      templateId tName(ifc,icateg);
	      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]
		->SetName( Form("%s_%s_%s_%s",
				tName.name().c_str(),systName[iSyst].c_str(),
				upDownName[iUpDown].c_str(),sbtag[ibtag].c_str()) );
	      
	      // read template errors
	      templateCore.assign("errorMassBTagTemplates/ErrorMassBTagTemplate");
	      string hbSystNameError( Form("%s_%s_%s_Cat%dTpat%d",templateCore.c_str(),
					   sfc[ifc].c_str(),sbtag[ibtag].c_str(),icateg,theTpat) );
	      hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc] 
		= getTrigsAndMerge(fBacSys,Form("%s_%s_%s_Cat%dTpat%d",
						templateCore.c_str(),
						sfc[ifc].c_str(),
						sbtag[ibtag].c_str(),icateg,theTpat),nTCombData,tCombData);
	      if ( hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc] == NULL ) {
		std::cout << "Hist not found: " << hbSystNameError << std::endl;
		return;
	      }
	      if (useTemplateError) {
		// add the template error
		std::cout << " ==== ErrorAdd ==== " << hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName() << std::endl;
		for (int ibinx=1; ibinx<= (hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetXaxis()->GetNbins()); ++ibinx) {
		  for (int ibiny=1; ibiny<= (hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetYaxis()->GetNbins()); ++ibiny) {
		    float oldError = hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetBinError(ibinx,ibiny);
		    float addError = hBackgroundSystError[iSyst][iUpDown][ibtag][icateg][ifc]->GetBinContent(ibinx,ibiny);
		    float newError = sqrt( oldError * oldError + addError * addError );
		    hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->SetBinError(ibinx,ibiny,newError);
		  }
		}
	      }
	      hout->cd();
	      hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->Write();
	    }
	  }
	}
	fBacSys->Close();
      }
    }
  }

  std::cout << "Everything done " << std::endl;

#ifdef PROJECTIONS
  // loop over background templates
  //TH1D* bgProX[nbtag][ncateg][nfc];
  TH1D* bgSystProX[nSyst][nUpDown][nbtag][ncateg][nfc];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    for (int icateg=0; icateg<ncateg; ++icateg) {
      for (int ifc=0; ifc<nfc; ++ifc) {
	TH2F* theTemp = hBackgroundCentral[ibtag][icateg][ifc];
	TH1D* theTempProX = theTemp->ProjectionX(Form("%sProX",theTemp->GetName()),0,-1,"e");
	theTempProX->SetName( Form("%sProX",theTemp->GetName()) );
	TH1D* theTempProY = theTemp->ProjectionY(Form("%sProY",theTemp->GetName()),0,-1,"e");
	std::cout << "Made projection " << theTempProX->GetName() << std::endl;

	if ( (icateg == 0) && (ifc == 0) ) {
	  for (int ibinx=1; ibinx<= theTempProX->GetXaxis()->GetNbins(); ++ibinx) {
	    std::cout << ibinx << " content " << theTempProX->GetBinContent(ibinx)
		      << " error " << theTempProX->GetBinError(ibinx) << std::endl;
	  }
	}
	theTempProX->SetMarkerStyle(20);
	theTempProX->SetMarkerColor(1);
	theTempProX->SetLineColor(1);
	theTempProX->SetMarkerSize(1);
	std::cout << "Draw" << std::endl;
	theTempProX->Draw("EP");
	//theTempProX->Draw("LP,SAME");
	// draw the SFbc systematics
	int colSyst[3] = {1, 2, 4};
	int lstyleUpDown[2] = {1, 1};
	for (int iSyst=1; iSyst<nSyst; ++iSyst) {
	  for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
	    bgSystProX[iSyst][iUpDown][ibtag][icateg][ifc] =  hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->ProjectionX(Form("%sProXX",hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName()),0,-1,"e");
	    bgSystProX[iSyst][iUpDown][ibtag][icateg][ifc]->SetName( Form("%sProXX",hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc]->GetName()) );
	    std::cout << "Made projection " << bgSystProX[iSyst][iUpDown][ibtag][icateg][ifc]->GetName() << std::endl;
	     bgSystProX[iSyst][iUpDown][ibtag][icateg][ifc]->SetLineColor( colSyst[iSyst] );
	     bgSystProX[iSyst][iUpDown][ibtag][icateg][ifc]->Draw("HIST,SAME");


// 	    TH2F* theTemp = hBackgroundSyst[iSyst][iUpDown][ibtag][icateg][ifc];
// 	    TH1D* theTempProX = theTemp->ProjectionX(Form("%sProX",theTemp->GetName()),0,-1,"e");
// 	    theTempProX->SetLineColor( colSyst[iSyst] );
// 	    TH1D* theTempProY = theTemp->ProjectionY(Form("%sProY",theTemp->GetName()),0,-1,"e");
// 	    theTempProY->SetLineColor( colSyst[iSyst] );
// 	    //theTempProX->Draw("HIST");
	  }
	}

	canvas->Print(Form("Template_%s_%s_Cat%d_ProX.png",sbtag[ibtag].c_str(),sfc[ifc].c_str(),icateg));
      }
    }
  }
#endif

  hout->Write();
  hout->Close();

  // close the signal central file
  fSig->Close();

  return;
}

void packTemplatesX() {
  //int iMass = 2;
#include "Analysis/Utilities/interface/HbbMass.h"
  for (int iMass=0; iMass<nSignal; ++iMass) {
    packTemplatesMass(iMass);
  }
}
