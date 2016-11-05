{
  cout<<"Start!\n";
  #include <string>
  std::string base(".L ~/Dropbox/mmtrigger/cxx_root/");///n/atlascode/backedup/stchan/mmtrigger/cxx_root/");
  gROOT->ProcessLine((base+"MMT_struct.cxx+").c_str());
  gROOT->ProcessLine((base+"MMT_Loader.cxx+").c_str());
  gROOT->ProcessLine((base+"MMT_Fitter.cxx+").c_str());
  gROOT->ProcessLine((base+"MMT_Finder.cxx+").c_str());
  gROOT->ProcessLine((base+"MM_Trigger.cxx+").c_str());
  gROOT->ProcessLine(".q");
  MM_Trigger lark("../ntuple/NSWPRDValAlg_allDet_xxuv_1000GeV15_vxp_001_001_70.root");
  lark.analysis(2);
//   delete lark;
}
//   MM_Trigger *lark=new MM_Trigger("/n/atlascode/backedup/stchan/mmtrigger/ntuple/blc_home/NSWPRDValAlg_allDet_xxuv_1000GeV15_vxp_001_001_70.root");
