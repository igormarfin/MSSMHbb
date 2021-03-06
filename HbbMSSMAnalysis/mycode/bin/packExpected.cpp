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

TCanvas* canvas;

TH2F* mergeSignal(TFile* fa,char* hname,const char* newname) {
  TH2F* theData[2];
  if (fa == NULL) {
    std::cout << "mergeSignal: no hfile opened" << std::endl;
    return NULL;
  }

  // weights
  double tWeight[2];
  tWeight[0] = (524.90+265.75)*.94*0.94 + 193.58+251.12+453.33;
  tWeight[1] = 246.53+732.73;
  double sumWeight = tWeight[0]+tWeight[1];
  for (int ii=0; ii<2; ++ii) {
    tWeight[ii] = tWeight[ii] / sumWeight;
  }
  std::cout << "Trig weights: " << tWeight[0] << " " << tWeight[1] << std::endl;
  
  for (int ii=0; ii<2; ++ii) {
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
  mergedData->Add(theData[0],theData[1],tWeight[0],tWeight[1]);
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


TH2F* mergeBackground(TFile* fa,const char* pattern,const char * theSbtag,const char* newname) {
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
  if ( sb == "CSVT" ) {
    wScale = 1.07817;
  } else if ( sb == "TCHPT" ) {
    wScale = 0.970401;
  }
  mergedData->Scale( wScale );
  std::cout << newname << " After final normalization: " << mergedData->GetSumOfWeights() << std::endl;

  return mergedData;
}




void packExpected() {
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

  const int nSignal=7;
  int signalMass[nSignal] = { 90, 100, 120, 140, 180, 250, 350 };
  int iMass = 6;
  double efficiency[nSignal] = { 0.0022081, 0.00324694, 0.00600146, 0.00918135,
				0.0138382, 0.0189684, 0.0206572 };
  double intLumi = 2.66794; // in fb-1
  double fScal[nSignal][nbtag];
  
  // signal templates
  const int nSyst = 3;
  const int nUpDown = 2;
  std::string signalFile( Form("/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-2/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola/SF/job_1/TripleBtagAnalysis.root",signalMass[iMass]) );
  std::string signalSystFiles[nSyst][nUpDown];
  std::string systName[nSyst] = { "JES", "SFbc", "SFudsg" };
  std::string upDownName[nUpDown] = { "Up", "Down" };
  TH2F* hSignalSyst[nSyst][nUpDown][nbtag];

  

  // output file
  TFile* hout = new TFile(Form("ExpectedLimit-M-%d.root",signalMass[iMass]),"recreate");
  hout->cd();
  TH2::AddDirectory(true);

  TFile* fSig = new TFile( signalFile.c_str() );
  if ( fSig == NULL ) {
    std::cout << "Could not open signal central file " << signalFile.c_str() << std::endl;
    return;
  }
  TH2F* hSignalCentral[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hSignalCentral[ibtag] = mergeSignal(fSig,Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()),
				       "bbH");
    double normShould = 1000 * intLumi * efficiency[iMass];
    double normIs = hSignalCentral[ibtag]->GetSumOfWeights();
    std::cout << hSignalCentral[ibtag]->GetName() << " TotalContents=" << hSignalCentral[ibtag]->GetSumOfWeights()
	      << std::endl;
    fScal[iMass][ibtag] = normShould / normIs;
    std::cout << "normShould = " << normShould << " normIs " << normIs
	      << " rescale by " << fScal[iMass][ibtag] << std::endl;
    hSignalCentral[ibtag]->Scale( fScal[iMass][ibtag] );
    hout->cd();
    hSignalCentral[ibtag]->Write();
  }


  for (int iSyst=0; iSyst<nSyst; ++iSyst) {
    for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
      signalSystFiles[iSyst][iUpDown] = Form( "/data/user/marfin/CMSSW_5_0_1/src/Analysis/HbbMSSMAnalysis/test/Systematics-test-2/theMergeList-SUSYBBHToBB_M-%d_7TeV-pythia6-tauola/%s_Sys%s/job_1/TripleBtagAnalysis.root",signalMass[iMass],systName[iSyst].c_str(),upDownName[iUpDown].c_str() );
      TFile* fSigSys = new TFile( signalSystFiles[iSyst][iUpDown].c_str() );
      if ( fSigSys == NULL ) {
	std::cout << "Could not open signal syst file " << signalSystFiles[iSyst][iUpDown].c_str() << std::endl;
	return;
      }
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	hSignalSyst[iSyst][iUpDown][ibtag] = mergeSignal(fSigSys,Form("massEvBtag/mjjEvBTag_%s",sbtag[ibtag].c_str()),
							 Form("bbH_%s_%s",systName[iSyst].c_str(),upDownName[iUpDown].c_str()));
	std::cout << "The merged hist has name " << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << std::endl;
	std::cout << hSignalSyst[iSyst][iUpDown][ibtag]->GetName() << " TotalContents=" << hSignalSyst[iSyst][iUpDown][ibtag]->GetSumOfWeights()
	      << std::endl;

	hSignalSyst[iSyst][iUpDown][ibtag]->Scale( fScal[iMass][ibtag] );
	hout->cd();
	hSignalSyst[iSyst][iUpDown][ibtag]->Write();
      }
      //fSigSys->Close();
    }
  }

  // backgrounds
  std::string backgroundFile( "/afs/naf.desy.de/user/s/spiridon/scratch/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/v1/Data-Run2011AB/TripleBtagAnalysis_SFzero/TripleBtagAnalysis.root" );
  TFile* fBac = new TFile( backgroundFile.c_str() );
  if ( fBac == NULL ) {
    std::cout << "Could not open background central file " << signalFile.c_str() << std::endl;
    return;
  }
  TH2F* hBackgroundCentral[nbtag];
  for (int ibtag=0; ibtag<nbtag; ++ibtag) {
    hBackgroundCentral[ibtag] = mergeBackground(fBac,"bgPredict/MassBTagPred_%s_%s_Cat%dTpat%d",sbtag[ibtag].c_str(),"BBB");
    hout->cd();
    hBackgroundCentral[ibtag]->Write();
  }

  std::string backgroundSystFiles[nSyst][nUpDown];
  std::string systNameSasha[nSyst] = { "JES", "SFbc", "SFq" };
  std::string upDownNameSasha[nUpDown] = { "plus2", "minus2" };
  TH2F* hBackgroundSyst[nSyst][nUpDown][nbtag];

  for (int iSyst=1; iSyst<nSyst; ++iSyst) {
    for (int iUpDown=0; iUpDown<nUpDown; ++iUpDown) {
      backgroundSystFiles[iSyst][iUpDown] = Form( "/afs/naf.desy.de/user/s/spiridon/scratch/CMSSW_4_2_4_patch1/src/Analysis/HbbMSSMAnalysis/test/results/v1/Data-Run2011AB/TripleBtagAnalysis_%s%s/TripleBtagAnalysis.root",systNameSasha[iSyst].c_str(),upDownNameSasha[iUpDown].c_str() );
      TFile* fBacSys = new TFile( backgroundSystFiles[iSyst][iUpDown].c_str() );
      if ( fBacSys == NULL ) {
	std::cout << "Could not open background syst file " << backgroundSystFiles[iSyst][iUpDown] << std::endl;
	return;
      }
      for (int ibtag=0; ibtag<nbtag; ++ibtag) {
	hBackgroundSyst[iSyst][iUpDown][ibtag] = mergeBackground(fBacSys,
			      "bgPredict/MassBTagPred_%s_%s_Cat%dTpat%d",
			     sbtag[ibtag].c_str(),Form("BBB_%s_%s",systName[iSyst].c_str(),upDownName[iUpDown].c_str()));
	hout->cd();
	hBackgroundSyst[iSyst][iUpDown][ibtag]->Write();
      }
    }
  }

  hout->Write();
  hout->Close();

  return;
}
