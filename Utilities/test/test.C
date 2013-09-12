void test()
{
cout<<gROOT->GetClass("MVAComputer")->GetName()<<endl;
cout<<gROOT->GetClass("MVATrainer")->GetName()<<endl;
 gROOT->GetClass("MVAComputer")->Dump();
 gROOT->GetClass("MVATrainer")->Dump();


}
