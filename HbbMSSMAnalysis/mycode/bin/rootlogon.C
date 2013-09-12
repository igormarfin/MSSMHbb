{

     gSystem->Load("libFWCoreFWLite");
     gSystem->Load("libAnalysisUtilities");
     AutoLibraryLoader::enable();

     cout<<"Hi"<<endl;
   gSystem->AddIncludePath("-DMEDIUM2012=1");
//   gSystem->AddIncludePath("-DSYST=1");
//    gSystem->AddIncludePath("-DPU=1");
//    gSystem->AddIncludePath("-DEXCLUDE_TRAINER_FROM_SELECTION=1");
//    gSystem->AddIncludePath("-DPU_ALEX=1");


 }  
