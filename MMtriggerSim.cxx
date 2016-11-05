#include "MM_Trigger.h"
const bool lxplus=false;

int main(int argc, char **argv){
  int nev=-1;
  string mmt_base=(lxplus?"/afs/cern.ch/work/s/stchan/mmtrigger/":"/n/atlasfs/atlascode/backedup/stchan/mmtrigger/");
  if(argc<4){
    par_par ok=dlm;//ok.genbg=true;//ok.set_mal_par(5,1.e-3);
    //ok.set_mal_par(2,5.);ok.ctx=2;ok.ctuv=1;//ok.set_correct_to_neg_misal();
    MM_Trigger lark = MM_Trigger(mmt_base+"ntuple/dlm_new/NSWPRDValAlgPt100GeV.digi.ntuple.root",ok);
    lark.analysis(nev);
    lark.save_file("/n/atlasfs/atlasdata/atlasdata1/stchan/mmtrigger/raw_hist/test/","lark");
  }
  else{
    par_par ok;
    for(int i=4;i<argc;i++)ok.set_parameter(argv[i]);
    cout<<"the par_par looks like: "<<ok.print_pars()<<"; file: "<<argv[1]<<endl;
    MM_Trigger lark = MM_Trigger(argv[1],ok);
    lark.save_analysis(argv[2],nev,argv[3]);
  }
  return 0;
}
