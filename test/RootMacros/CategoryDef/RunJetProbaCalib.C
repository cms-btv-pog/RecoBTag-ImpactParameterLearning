{


 gROOT->Reset() ;
 
  
 gROOT->ProcessLine(".L CategoryDef.C++");
 gSystem->Load("CategoryDef_C.so");

 gROOT->ProcessLine(".L CategoryDefCollection.C++");
 gSystem->Load("CategoryDefCollection_C.so");

 gROOT->ProcessLine(".L CategoriesDefinition.C++");
 gSystem->Load("CategoriesDefinition_C.so");


 
 cout << "compile JetProbaCalib " << endl;
 gROOT->ProcessLine(".L JetProbaCalib.C+g") ;


 //if (gROOT->GetClass("btagana/ttree")==0) return;
 
 cout << "initialize ttree " << endl;
 TChain c("btagana/ttree");


 cout << "add root files to the ttree" << endl;
 

 c.Add("../../../../PerformanceMeasurements/test/8bdaf0383a3c233e3e7fde7c7f1c868d/*.root");


 cout << "construct  JetProbaCalib" << endl;
 JetProbaCalib* theanalyzer = new JetProbaCalib(&c);
 
 cout << "loop on events " << endl;
 
 
 theanalyzer->Loop();
 
 
 
 
}
