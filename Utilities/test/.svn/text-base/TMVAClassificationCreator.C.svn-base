///usage:  help() {  sed  -ne 's/^\/\/\/usage://p' <(grep usage $1 ) } 
///usage: 
///usage:
///usage:  new help
///usage: sed  's/^\/\/\/usage://p' <(grep usage TMVAClassificationCreator.C )
///usage: 
///usage: root -l `cat tmva_list | awk '{print $1}'` 'TMVAClassificationCreator.C("tmva_list","Likelihood")'
///usage:  where tmva_list must have 7 fields: 
///usage: /data/user/marfin/MVAPythia/theMergeList-SUSYBBHToBB_M-100_7TeV-pythia6-tauola.root  SUSYBBHToBB_M-100_7TeV-pythia6-tauola SGN 1.0  1.0  1.0 1
///usage: 
///usage: 
///usage: 
///usage: 
///usage:  ls ../*root | grep -v "readtribn" | xargs -I {} echo  " echo -n \`pwd\`/{} ;echo -n \" \" ; aa=\`basename {} | sed -ne 's/theMergeList-|//p; s/.root//p;' \`;  echo -n \" \$aa \" ; bb=\`expr index \$aa SUSY\`; [[ \$bb < 1 ]] &&  echo -n  \" BKG \";  [[ \$bb > 0 ]] &&  echo -n  \" SGN \";  echo -n \" 1.0 \"; echo -n \" 1.0 \";echo -n \" 1.0 \"; echo  \" 1 \"; " | sh > mylist;
///usage: root -l `cat mylist | awk '{print $1}'` 'TMVAClassificationCreator.C("mylist","Likelihood")'
///usage:  rm mylist;
///usage: 
///usage: 
///usage: 
///usage:  if we want to restrict our samples to some concrete files:
///usage:  mm="*";mm="M-120|50_150";ls ../*root | grep -v "readtribn" | egrep  $mm | xargs -I {} echo  " echo -n \`pwd\`/{} ;echo -n \" \" ; aa=\`basename {} | sed -ne 's/theMergeList-|//p; s/.root//p;' \`;  echo -n \" \$aa \" ; bb=\`expr index \$aa SUSY\`; [[ \$bb < 1 ]] &&  echo -n  \" BKG \";  [[ \$bb > 0 ]] &&  echo -n  \" SGN \";  echo -n \" 1.0 \"; echo -n \" 1.0 \";echo -n \" 1.0 \"; echo  \" 1 \"; " | sh > mylist;
///usage:  root -l `cat mylist | awk '{print $1}'` 'TMVAClassificationCreator.C("mylist","Likelihood")'
///usage:  rm mylist;
///usage: 
///usage: 
///usage: how to create exclusion list of varables:
///usage:  source utils.sh
///usage:
///usage: cccc=( $(ReadXMLParameter ../tree0.xml Variable VarIndex | tr ' ' '\n'))
///usage: dddd=( $(ReadXMLParameter ../tree0.xml Variable Expression | tr ' ' '\n'))
///usage: [[ -f exclude.lst ]] && rm exclude.lst ;touch exclude.lst; for (( i=1;i<`echo ${#cccc}`;i++ )); do [[ $i -gt 25 ]] && echo ${dddd[$i]} >> exclude.lst; done
///usage:
///usage: last example, how to use exclusion list and new TTree 'KinVars4'
///usage:  root -l `cat mylist | awk '{print $1}'` 'TMVAClassificationCreator.C("mylist","Likelihood","KinVars4","exclude.lst")'
///usage:
///usage: example how to use option_book list
///usage: root -l `cat mylist_tmp_dir2 | awk '{print $1}'` 'TMVAClassificationCreator.C("mylist_tmp_dir2","LikelihoodD","KinVars4","exclude.lst","options_book.lst")'
///usage:
///usage:
///usage: example of likelihood with smaller number of evetns for trainig:
///usage: root -l `cat mylist_tmp_dir3 | awk '{print $1}'` 'TMVAClassificationCreator.C("mylist_tmp_dir3","Likelihood","KinVars4","exclude.lst","options_book_likelihood.lst","options_train_test.lst")'
///usage:
///usage: example of CFMlpANN on 50K events for training
///usage: root -l `cat mylist_tmp_dir3 | awk '{print $1}'` 'TMVAClassificationCreator.C("mylist_tmp_dir3","CFMlpANN","KinVars4","exclude.lst","options_book_CFMlpANN.lst","options_train_test_50k.lst")'
///usage:



// @(#)root/tmva $Id: TMVAClassification.C,v 1.36 2009-04-14 13:08:13 andreas.hoecker Exp $
/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l TMVAClassification.C\(\"Fisher,Likelihood\"\)                       *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set is used.                                     *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 **********************************************************************************/

#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <algorithm>
#include <vector>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"

#include "TMVAGui.C"

#include "functions.h"

#include "TMVA/Config.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif


  
void TMVAClassificationCreator( TString list_of_files="list_of_files",TString myMethodList = "", TString _treeName="", TString _exclusion="", TString _options_training="", TString _train_test_options = "" ) ///usage:
{

//(TMVA::gConfig().GetVariablePlotting()).fNbins1D = 50.0;
//(TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve=50;

///Set search path properly:
TString  env =  gEnv->GetValue("Unix.*.Root.MacroPath","A");
FILE * myPipe = gSystem->OpenPipe("pwd","r");
TString pwd;
pwd.Gets(myPipe);
gSystem->ClosePipe(myPipe);

env=env + pwd  + ":";

//gEnv->SetValue("Unix.*.Root.MacroPath",env.Data());
gROOT->SetMacroPath(env.Data());


///Load function to check how many fields we have inside list_of_files
///HowManyFieldsInFile.C was created by calling: 'createCPPCode HowManyFieldsInFile'
/// HowManyFieldsInFile lives in utils.sh
gROOT->LoadMacro("HowManyFieldsInFile.C");

TString res=HowManyFieldsInFile(list_of_files);
Int_t nField=res.Atoi();

//cout<<nField<<endl;

///Part to read list_of_files

///All fields that might be are listed below
/// Add a new field and create a new version of Readinput if you add the field in list_of_file 
std::vector<std::string>  filesList;
std::vector<std::string> processList;
std::vector<std::string>  TypeList;
std::vector<double>  xsectsList;
std::vector<double>  effsList;
std::vector<double>  wgtList;
std::vector<int>  mergeList;
int Num;
/**
	Only filesList, TypeList, wgtList are relevant for this analysis!

**/


///Here we are reading
if (nField==(Int_t) 7) 
Readinput_v2(list_of_files,filesList, processList, TypeList, xsectsList,  effsList ,wgtList,mergeList,Num);
else {std::cout<<"There's no Readinput_v function to read the list_of_files"<<std::endl;return;}

// define what sample we use
   std::map<std::string,double> SampleWeight;
   std::map<std::string,double> SampleType; ///>0 -sig ,<0 -bkg , 0 -- data!, data must be skiped!


///Fill SampleWeight && SampleType
for (int k = 0; k<Num; k++) {





TString type=TypeList[k];

if (type.Contains("SGN")) SampleType[filesList[k]]=1.0;
if (type.Contains("BKG")) SampleType[filesList[k]]=-1.0;
if (type.Contains("DATA")) SampleType[filesList[k]]=0.0;

//cout<<"Weight = "<<wgtList[k]<<endl;

//if  (SampleType[filesList[k]]<0)
//SampleWeight[filesList[k]]=0.0510891/113.373;
//else
SampleWeight[filesList[k]]=wgtList[k];

cout<<"Weight = "<<wgtList[k]<<endl;

}


///needed files with TTrees

//cout<<"Check"<<endl;

	if (gROOT->GetListOfFiles()->GetEntries()==0) return; /// no input files 

//cout<<"Check"<<endl;	




   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the 
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   // 
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   // this loads the library
   TMVA::Tools::Instance();

   //---------------------------------------------------------------
   // default MVA methods to be trained + tested
   std::map<std::string,int> Use;

TString _nameMVA="";

   // --- Allowed MVA discriminators
   Use["Likelihood"]      = 1; ///usage: 
   Use["LikelihoodD"]     = 1; ///usage: the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 1; ///usage: the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 1;  ///usage: 
   Use["LikelihoodMIX"]   = 1; ///usage:    
Use["LikelihoodBoost"]   = 1;  ///usage: 
   Use["PDERS"]   = 1;  ///usage: 
   Use["MLP"]   = 1; ///usage: 
   Use["BDT"]   = 1; ///usage: 
   Use["BDTSmallMaxDepth"]=1; ///usage: 
 Use["BDTBagging"] =1; ///usage: 
 Use["BDTGrad"]=1; ///usage: 
 Use["MLPBNN"]=1; ///usage: 
 Use["CFMlpANN"]=1; ///usage: 
 Use["TMlpANN"]=1; ///usage: 
   // ---
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }





   // Create a new root output file.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory will
   // then run the performance analysis for you.
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/ 
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in 
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile, 
//                                               "!V:!Silent:Color:DrawProgressBar" );
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;G,D;N" );
//                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D" );

   // If you wish to modify default settings 
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]


///Suppoused name of tree
TString treeName;
if (_treeName.Length()>0) treeName=_treeName;
else treeName="KinVars";

///Our tree
TTree * tree=0;

///Here is all varibles used for learning && variable skipped
/// set to 0 if you don't want to use them!
std::map<std::string,int> variables_to_use;
	


///Let's extrack TTree from first input file!
for (Int_t j=1;j<gROOT->GetListOfFiles()->GetEntries();j++)
{
         ((TFile*) gROOT->GetListOfFiles()->At(j))->cd();
	if (!TString(gFile->GetName()).Contains("TMVA")) break;
}
// ((TFile*) gROOT->GetListOfFiles()->At(1))->cd();

   tree = (TTree *)gFile->Get(treeName);
if (!tree) {std::cout<<"Can't extract tree from "<<gFile->GetName()<<std::endl; return ;}
	TObjArray *branches = tree->GetListOfBranches();

///Let's fill variables_to_use
	for (Int_t jjj=0;jjj<tree->GetListOfBranches()->GetSize();jjj++)
	{
		TObject * obj = tree->GetListOfBranches()->At(jjj);
		if (obj) {
		std::string str(obj->GetName());
  	        variables_to_use[str]=1;
}
	}


//return;

///Variables to skip

/*

variables_to_use["sumEt"]=0;
variables_to_use["Et1"]=0;
variables_to_use["Et2"]=0;
variables_to_use["Et3"]=0;
//variables_to_use["maxdijetmass"]=0;
//variables_to_use["detajet3jet4"]=0;
//variables_to_use["detajet2jet3"]=0;

//variables_to_use["dthetajet1jet2_boost"]=0;
//variables_to_use["dphijet1jet2_boost"]=0;
//variables_to_use["dphijet1jet2"]=0;
variables_to_use["maxdijetmass"]=0;
variables_to_use["mindijetmass"]=0;
variables_to_use["sphericity_boost"]=0;
//variables_to_use["sphericity"]=0;
variables_to_use["aplanarity"]=0;
variables_to_use["nEvt"]=0;



///Correlated

//variables_to_use["dptjet1jet2"]=0;


//variables_to_use["Et2byEt1"]=0; // 1 -- for all case, 0 -- 4 heavy flavour 
//variables_to_use["Et3byEt1"]=0; // 1 -- for all case, 0 -- 4 heavy flavour 
//variables_to_use["Et3byEt1"]=1; ///was 0
//variables_to_use["dptjet1jet3"]=0; ///was 0 
variables_to_use["djet1jet2pt"]=0; ///was 0

//variables_to_use["Eta3"]=0;  ///was 0
//variables_to_use["Eta2"]=0;  ///was 0



//variables_to_use["Et3byEt2"]=0;
//variables_to_use["detajet1jet2"]=0;
//variables_to_use["detajet2jet3"]=0;
variables_to_use["detajet3jet4"]=0;
variables_to_use["D"]=0;
variables_to_use["isotropy"]=0;
*/


/*
        <Variable VarIndex="0" Expression="detajet2jet3"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="1" Expression="detajet3jet4" Type="F" BINS="100"  MIN="0" MAX="3.5"/> 
        <Variable VarIndex="2" Expression="maxdijetmass" Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="3" Expression="mindijetmass"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="4" Expression="aplanarity"  Type="F" BINS="100"  MIN="0" MAX="0.5"/>
        <Variable VarIndex="5" Expression="sphericity"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="6" Expression="Et1"  Type="F" BINS="100"  MIN="0" MAX="1000"/>
        <Variable VarIndex="7" Expression="Et2"  Type="F" BINS="100"  MIN="0" MAX="1000"/>
        <Variable VarIndex="8" Expression="Et3"  Type="F" BINS="100"  MIN="0" MAX="1000"/>
        <Variable VarIndex="9" Expression="Eta1"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="10" Expression="Eta2"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="11" Expression="Eta3"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="12" Expression="dptjet1jet2"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="13" Expression="dptjet1jet3"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="14" Expression="djet1jet2pt"  Type="F" BINS="100"  MIN="0" MAX="1"/>
        <Variable VarIndex="15" Expression="detajet1jet2"  Type="F" BINS="100"  MIN="0" MAX="3.5"/>
        <Variable VarIndex="16" Expression="dthetajet1jet2_boost"  Type="F" BINS="100"  MIN="-7.0" MAX=7.0"/>
        <Variable VarIndex="17" Expression="dphijet1jet2_boost"  Type="F" BINS="100"  MIN="-4.0" MAX=4.0"/>
        <Variable VarIndex="18" Expression="dphijet1jet2"  Type="F" BINS="100"  MIN="-4.0" MAX=4.0"/>
        <Variable VarIndex="19" Expression="dphijet2jet3"  Type="F" BINS="100"  MIN="-4.0" MAX=4.0"/>
        <Variable VarIndex="20" Expression="dphijet1jet3"  Type="F" BINS="100"  MIN="-4.0" MAX=4.0"/>
        <Variable VarIndex="21" Expression="Et2byEt1"  Type="F" BINS="100"  MIN="0.0" MAX=1.0"/>
        <Variable VarIndex="22" Expression="Et3byEt1"  Type="F" BINS="100"  MIN="0.0" MAX=1.0"/>
        <Variable VarIndex="23" Expression="Et3byEt2"  Type="F" BINS="100"  MIN="0.0" MAX=1.0"/>
        <Variable VarIndex="24" Expression="sphericity_boost"  Type="F" BINS="100"  MIN="0" MAX="1"/>

*/




///weight as well
variables_to_use["weight"]=0;

///if _exclusion list exists
if (_exclusion.Length()>0)
{

ifstream f(_exclusion);

///Parsing of file:
TString  str222;

///Read first line
str222.ReadLine(f);




while (str222.Length()>0 ) {

//cout<<"str222 = "<<str222<<endl;

if (!str222.Contains("/")  && str222.Length()>0 ) variables_to_use[str222.Data()]=0;
if (str222.Contains("/") && str222.Contains("NameMVA") ) {
TString _nameMVA2 = str222; _nameMVA2.ReplaceAll("/","");
 _nameMVA2.ReplaceAll("NameMVA=","");

_nameMVA+=_nameMVA2;

}



str222.ReadLine(f);

}

}


 for (std::map<std::string,int>::iterator it2 = variables_to_use.begin(); it2 != variables_to_use.end(); it2++)
{	
//	cout<<it2->second<<endl;
	if (it2->second	>0) {
        factory->AddVariable(TString(it2->first), TString(it2->first), 'F' );
	 cout<<"I'm adding the variable "<<it2->first<<endl;
//	cout<<"2 "<<it2->first<<endl;

	}
	else 
{
	cout<<"I'm skipping the variable "<<it2->first<<endl;

	}	
}

///Add weight as spectacor.  From Roberval's code:
//factory->AddSpectator("weight");




///Used to calculate Optimal cut
Double_t Nsgn=0;
Double_t Nbkg=0;

   // You can add so-called "Spectator variables", which are not used in the MVA training, 
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the 
   // input variables, the response values of all trained MVAs, and the spectator variables


for (Int_t j=0;j<gROOT->GetListOfFiles()->GetEntries();j++)
{
	 ((TFile*) gROOT->GetListOfFiles()->At(j))->cd();

	 if (TString(gFile->GetName()).Contains("TMVA")) continue;

	cout<<"Current file is "<<gFile->GetName()<<endl;
	if (SampleType.find(gFile->GetName())==SampleType.end() || SampleWeight.find(gFile->GetName()) == SampleWeight.end()) continue;

	cout<<"SampleType is "<<SampleType[gFile->GetName()]<<endl;
	cout<<"SampleWeight is "<<SampleWeight[gFile->GetName()]<<endl;

	 tree = (TTree *)gFile->Get(treeName);
	if (!tree) continue;

// ====== register trees ====================================================
	if ( SampleType[gFile->GetName()]>0) { 
	factory->AddSignalTree    ( tree,    SampleWeight[gFile->GetName()]      ); 
	cout<<"Signal "<<gFile->GetName()<<" is added"<<endl;
	tree->Draw("weight>>tmpth1");
	cout<<endl<<endl<<endl; 	
	TH1F * tmpth1 =  (TH1F*)gDirectory->Get("tmpth1"); 
	for (Int_t k=1;k<=tmpth1->GetNbinsX();k++) {
//		cout<<"I'm adding Nevents  as "<< SampleWeight[gFile->GetName()]*tmpth1->GetBinContent(k)*tmpth1->GetBinCenter(k)<<endl;	
	Nsgn+= SampleWeight[gFile->GetName()]*tmpth1->GetBinContent(k)*tmpth1->GetBinCenter(k);
	}
	cout<<"The current number of signal events for optimal Cut finding is "<<Nsgn<<endl;
	cout<<endl<<endl<<endl; 	

	}
	if ( SampleType[gFile->GetName()]<0) {
	factory->AddBackgroundTree    ( tree,    SampleWeight[gFile->GetName()]      ); 
	cout<<"Bkg "<<gFile->GetName()<<" is added"<<endl;
	cout<<endl<<endl<<endl; 	
	tree->Draw("weight>>tmpth1"); 	TH1F * tmpth1 =  (TH1F*)gDirectory->Get("tmpth1"); 
	 for (Int_t k=1;k<=tmpth1->GetNbinsX();k++) {
//   cout<<"I'm adding Nevents  as "<< SampleWeight[gFile->GetName()]*tmpth1->GetBinContent(k)*tmpth1->GetBinCenter(k)<<endl;
	Nbkg+= SampleWeight[gFile->GetName()]*tmpth1->GetBinContent(k)*tmpth1->GetBinCenter(k);
	}
	cout<<"The current number of bkg events for optimal Cut finding is "<<Nbkg<<endl;
	cout<<endl<<endl<<endl; 	


	}



}	

   // This would set individual event weights (the variables defined in the
   // expression need to exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
  
   factory->SetSignalWeightExpression("weight");
   factory->SetBackgroundWeightExpression("weight");

		


TString test_options = "";


///if _exclusion list exists
if (_train_test_options.Length()>0)
{

ifstream f(_train_test_options);

///Parsing of file:
TString  str222;

///Read first line
str222.ReadLine(f);


while (str222.Length()>0 ) {

//cout<<"str222 = "<<str222<<endl;

if (!str222.Contains("/")  && str222.Length()>1 ) test_options+=str222+TString(":");

if (str222.Contains("/") && str222.Contains("NameMVA") ) {
TString _nameMVA2 = str222; _nameMVA2.ReplaceAll("/","");
 _nameMVA2.ReplaceAll("NameMVA=","");

_nameMVA+=_nameMVA2;

}


str222.ReadLine(f);

}

}


if (test_options.Length()>0) test_options.Chop();
else test_options =  "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:SplitSeed=100:NormMode=EqualNumEvents:!V";

cout<<"test_options = "<< test_options<<endl;
      
   

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // tell the factory to use all remaining events in the trees after training for testing:
   factory->PrepareTrainingAndTestTree( mycuts, mycutb, test_options.Data());
//                                        "nTrain_Signal=0:nTrain_Background=0:nTest_Signal=30000:nTest_Background=30000:SplitMode=Random:SplitSeed=0:NormMode=NumEvents:!V" );
//                                        "nTest_Signal=0:nTest_Background=15000:SplitMode=Block:NormMode=NumEvents:!V" );
//                                       "nTest_Signal=300000:nTest_Background=300000:SplitMode=Block:NormMode=NumEvents:!V" );
//                                       "nTest_Signal=700000:nTest_Background=700000:SplitMode=Block:NormMode=NumEvents:!V" );

///I'm using NormMode=None to skip changing weights!
//                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:SplitSeed=100:NormMode=None:!V" );
//                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:SplitSeed=100:NormMode=EqualNumEvents:!V" );
//                                        "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:SplitSeed=100:NormMode=EqualNumEvents:!V" );
//                                        "nTrain_Signal=100000:nTrain_Background=100000:SplitMode=Random:SplitSeed=0:NormMode=NumEvents:!V" );
//                                        "nTrain_Signal=100000:nTrain_Background=100000:SplitMode=Random:SplitSeed=10:NormMode=NumEvents:!V" );
//                                        "nTrain_Signal=100000:nTrain_Background=300000:SplitMode=Block:SplitSeed=0:NormMode=NumEvents:!V" );
//                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:SplitSeed=0:NormMode=NumEvents:!V" );
//                                        "nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:SplitSeed=0:NormMode=EqualNumEvents:!V" );
//				         "nTest_Signal=3000:nTest_Background=3000:SplitMode=Random:SplitSeed=0:NormMode=NumEvents:!V" );

   // If no numbers of events are given, half of the events in the tree are used for training, and 
   // the other half for testing:
   //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );  
   // To also specify the number of testing events, use:
   //    factory->PrepareTrainingAndTestTree( mycut, 
   //                                            "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );  

   // ---- Book MVA methods
   //
   // please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

//return;


///From Roberval's code + my modification


/// if we need to book several options among many variables!

TString _option_book ="";

/// _options_training may contain the file name with booking options for the methods

if (_options_training.Length()>0)
{

ifstream f(_options_training);

///Parsing of file:
TString  str222;

///Read first line
str222.ReadLine(f);


while (str222.Length()>0 ) {

cout<<"1212 str222 = "<<str222<<endl;
if (!str222.Contains("/") && str222.Length()>1) {
_option_book = _option_book +  str222 + TString(":");
cout<<"!!!!!! _option_book = "<<_option_book<<endl;
}

if (str222.Contains("/") && str222.Contains("NameMVA") ) {
TString _nameMVA2 = str222; _nameMVA2.ReplaceAll("/","");
 _nameMVA2.ReplaceAll("NameMVA=","");

_nameMVA+=_nameMVA2;

}
str222.ReadLine(f);

}

}

if (_option_book.Length()>0) _option_book.Chop();

cout<<" _option_book = "<<_option_book<<endl;

TString _main_option_book;
_main_option_book="H:!V";

if (_option_book.Length()>0)
_main_option_book+=TString(":")+_option_book;



 // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"])
   {


///_main_option_book="H:!V:!TransformOutput:PDFInterpol=Spline3:Nbins=50:NSmooth=15";



if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kLikelihood, _nameMVA, _main_option_book.Data());
else           factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood", _main_option_book.Data());

//	  "H:!V:!TransformOutput:CreateMVAPdfs=True:PDFInterpol=Spline5:Nbins=100:NSmooth=5" );
//	  "H:!V:!TransformOutput:PDFInterpol=Spline5:Nbins=100:NSmooth=5" );
//	  "H:!V:!TransformOutput:PDFInterpol=Spline1:Nbins=100:NSmooth=5:NbinsMVAPdf=50:NSmoothMVAPdf=10");
//	  "H:!V:!TransformOutput:PDFInterpol=Spline2:Nbins=40:MinNSmooth=150:MaxNSmooth=2050");
//	  "H:!V:!TransformOutput:PDFInterpol=Spline3:Nbins=50:MinNSmooth=150:MaxNSmooth=2050");
//	  "H:!V:!TransformOutput:PDFInterpol=Spline3:Nbins=50");
//	  "H:!V:!TransformOutput:PDFInterpol=Spline3:Nbins=50:NSmooth=15");
//	  "H:!V:!TransformOutput:PDFInterpol=Spline2:Nbins=100:NbinsBkg[0]=50:NbinsBkg[1]=50:NbinsBkg[2]=50:NbinsBkg[3]=50:NbinsBkg[4]=50:NbinsBkg[5]=50:NbinsBkg[6]=50:NbinsBkg[7]=50:NbinsBkg[8]=50" );
//        "H:!V:!TransformOutput:CreateMVAPdfs=True:PDFInterpol=Spline2:NAvEvtPerBin=60:NSmooth=1" );
//       "H:!V:!TransformOutput:CreateMVAPdfs=True:PDFInterpol=Spline0:Nbins=100:NSmooth=3" ); ///It's better of previous
//        "H:!V:!TransformOutput:PDFInterpol=Spline5:Nbins=100:NSmooth=5" );
   }

 // Likelihood ("naive Bayes estimator")
   if (Use["LikelihoodBoost"])
   {


	_main_option_book+=TString(":")+TString("Boost_Num=30:Boost_Type=AdaBoost");

if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kLikelihood, _nameMVA, _main_option_book.Data());
else                   factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodBoost", _main_option_book.Data());

   }


   // Decorrelated likelihood
   if (Use["LikelihoodD"])
   {

	_main_option_book+=TString(":")+TString("VarTransform=G,D");

if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kLikelihood, _nameMVA, _main_option_book.Data());
else        factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD", _main_option_book.Data());

   }

   if (Use["PDERS"])
   {


	_main_option_book+=TString(":")+TString("VolumeRangeMode=Adaptive:KernelEstimator=Gauss");
	  if (_option_book.Length()>0) _main_option_book+=TString(":")+_option_book;	
if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kPDERS, _nameMVA, _main_option_book.Data());
else                   factory->BookMethod( TMVA::Types::kPDERS, "PDERS", _main_option_book.Data());

   }


///
   if (Use["MLP"])
   {


/*        
         _main_option_book+=TString(":")+TString("NCycles=500");
        _main_option_book+=TString(":")+TString("HiddenLayers=N,N-1");
        _main_option_book+=TString(":")+TString("SamplingTraining=True");
        _main_option_book+=TString(":")+TString("Sampling=1");
        _main_option_book+=TString(":")+TString("SamplingEpoch=1");
        _main_option_book+=TString(":")+TString("SamplingImportance=1");
        _main_option_book+=TString(":")+TString("VarTransform=Norm");
        _main_option_book+=TString(":")+TString("NeuronInputType=sqsum");
        _main_option_book+=TString(":")+TString("TrainingMethod=BP");
        _main_option_book+=TString(":")+TString("NeuronType=tanh");
*/
	
	//  if (_option_book.Length()>0) _main_option_book+=TString(":")+_option_book;

if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kMLP, _nameMVA, _main_option_book.Data());
else           factory->BookMethod( TMVA::Types::kMLP, "MLP", _main_option_book.Data());

///too long
//   factory->BookMethod( TMVA::Types::kMLP, "MLPBNN4", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=N,N-1:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators
// factory->BookMethod( TMVA::Types::kMLP, "MLPBNN2", "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators
// factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

//      factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

// Tmlp(Root)ANN
//      factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...



   }


if ( Use["MLPBNN"]) 
{
//	  if (_option_book.Length()>0) _main_option_book+=TString(":")+_option_book;
	 _main_option_book+=TString(":UseRegulator");

if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kMLP, _nameMVA, _main_option_book.Data());
else	 factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", _main_option_book.Data());


}

if ( Use["CFMlpANN"]) 
{
	if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kCFMlpANN, _nameMVA, _main_option_book.Data());
else	 factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", _main_option_book.Data()  );
}


if ( Use["TMlpANN"])
{
	if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kTMlpANN, _nameMVA, _main_option_book.Data());
else         factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", _main_option_book.Data()  );
}



///
   if (Use["BDT"])
   {

/*
        _main_option_book+=TString(":")+TString("NTrees=200");
//        _main_option_book+=TString(":")+TString("NTrees=1");
        _main_option_book+=TString(":")+TString("BoostType=AdaBoost");
        _main_option_book+=TString(":")+TString("UseWeightedTrees=True");
        _main_option_book+=TString(":")+TString("NodePurityLimit=0.5");
//        _main_option_book+=TString(":")+TString("SeparationType=GiniIndex");
        _main_option_book+=TString(":")+TString("SeparationType=SDivSqrtSPlusB");
        _main_option_book+=TString(":")+TString("PruneStrength=-1");
//        _main_option_book+=TString(":")+TString("PruneStrength=1");
        _main_option_book+=TString(":")+TString("NCuts=20");
        _main_option_book+=TString(":")+TString("UseYesNoLeaf=True");
//        _main_option_book+=TString(":")+TString("PruneMethod=CostComplexity");

	 _main_option_book+=TString(":")+TString("PruneMethod=NoPruning");

*/

cout<<"Main options " <<_main_option_book.Data()<<endl;

if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kBDT, _nameMVA, _main_option_book.Data());
else           factory->BookMethod( TMVA::Types::kBDT, "BDT", _main_option_book.Data());

//      factory->BookMethod( TMVA::Types::kBDT, "BDT","!H:!V:NTrees=450:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.1:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning:NodePurityLimit=0.1" );

///works
//     factory->BookMethod( TMVA::Types::kBDT, "BDT","!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );

///works
//factory->BookMethod( TMVA::Types::kBDT, "BDT","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

   }

if (Use["BDTSmallMaxDepth"]) 
{
if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kBDT, _nameMVA, _main_option_book.Data());
else  factory->BookMethod( TMVA::Types::kBDT, "BDTSmallMaxDepth",_main_option_book.Data());

}

if (Use["BDTBagging"])
{

///works
//     factory->BookMethod( TMVA::Types::kBDT, "BDTBagging","!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );



if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kBDT, _nameMVA, _main_option_book.Data());
else  factory->BookMethod( TMVA::Types::kBDT, "BDTBagging",_main_option_book.Data());

//cout<<"Nmae" <<_nameMVA<<endl;
//cout<<"Main option"<<_main_option_book<<endl;

}

if (Use["BDTGrad"])
{


cout<<"options : "<< _main_option_book.Data()<<endl;

///works
//factory->BookMethod( TMVA::Types::kBDT, "BDTGrad","!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=20:NNodesMax=5" );

if (_nameMVA.Length()>0)  factory->BookMethod( TMVA::Types::kBDT, _nameMVA, _main_option_book.Data());
else  factory->BookMethod( TMVA::Types::kBDT, "BDTGrad",_main_option_book.Data());

}



   // PCA-transformed likelihood

   if (Use["LikelihoodPCA"])
/**
Here my modification:
        *remove: PDFInterpolSig[4]=Spline0:PDFInterpolBkg[4]=Spline0
        *change: NAvEvtPerBin to 30 for all variables!
        *add CreateMVAPdfs


**/

      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
   "H:!V:!TransformOutput:CreateMVAPdfs=True::PDFInterpol=Spline2:NAvEvtPerBin=30:VarTransform=PCA" );


   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"])

/**
Here my modification:
        *change: NAvEvtPerBin to 30 for all variables!
        *add CreateMVAPdfs                  

**/


      factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
       "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:Nbins=70" ); 

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"])

/**
Here my modification:
        *change: NAvEvtPerBin to 30 for all variables!
        *add CreateMVAPdfs


**/



     factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
           "!H:!V:CreateMVAPdfs=True:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=30" ); 




   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory->TrainAllMethodsForClassification();

cout<<"Afterr train"<<endl;

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

cout<<"Afterr test"<<endl;

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();    

cout<<"Afterr eval"<<endl;


   // --------------------------------------------------------------
   
   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;      


std::cout << "==>  Nsgn = "<<Nsgn<<std::endl;
std::cout << "==>  Nbkg = "<<Nbkg<<std::endl;



   delete factory;

   // Launch the GUI for the root macros

///(TMVA::gConfig().GetVariablePlotting()).fNbins1D = 40.0;

   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
