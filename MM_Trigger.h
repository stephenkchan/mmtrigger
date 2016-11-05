#ifndef MM_TRIGGER_H
#define MM_TRIGGER_H

//C++ language libararies
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

//ROOT libraries
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "MMT_struct.h" 
#include "MMT_Loader.h"
#include "MMT_Loader.h"
#include "MMT_Finder.h" 
#include "MMT_Fitter.h" 
using namespace std;

class MM_Trigger : public MMT_Loader{
 public:
  MM_Trigger(string fname, par_par var=par_par(), const string& tname="");
  ~MM_Trigger(){}
  void analysis(int nevents=-1);
  void save_file(const string& filedir,const string& opttag="",bool pcrep=false);
  void save_analysis(const string& filedir,int nevents=-1,const string& opttag="",bool pcrep=false);

/*   void Run_File(const string& background="bgon"); */
  
 protected:
  string filename(const string& filedir,const string& opttag)const;
  bool file_exists(const string&name)const;
  bool process_event(int event,ofstream&hit_printer);
  void print_header(ofstream& print)const;
  void print_ev_header(ofstream& print,int ev)const;
  void print_hit(ofstream& print,const hdst_entry&hit)const;
  void print_event_summary(ofstream& print,int ev=-999,double tru_th=-999,double tru_ph=-999,double tru_dth=-999,double mx=-999,double my=-999,double mxl=-999,double theta=-999,double phi=-999,double dtheta=-999)const;
  void print_rd_summary(const map<pair<int,int>,finder_entry>&road_tracks,int ev);

  //TH1's
  TH1D *m_fit_the,*m_fit_phi,*m_fit_dth;
  TH1D *m_theta_tru_entpos;
  TH2D *m_hcode;
  vector<TH1D*> m_res_the,m_res_phi,m_res_dth,m_res_mxl,m_res_mx,m_res_my,m_res_mu,m_res_mv,m_res_mxuv,m_res_mxn,m_all_charge,m_fit_charge;
  vector<vector<vector<TH1D*> > > m_crep_the,m_crep_phi,m_crep_dth;
  vector<TH2D*> m_bcid_eff;
  TH1D *m_store_no;

  TH1D *nev_ent_wedge,*nev_ct,*nev_fit_clean,*nev_fit_bg,*nev_bcid_pass;
  TH3D *nev_tru_xuv,*nev_fit_xuv,*nev_rec_xuv,*nev_rft_xuv,*nev_rpl_xuv;
  int nev,nana;
  //progress bar:
  void progress_bar(int current, int total);

 private:
  int n_etabins,n_phibins,correct_bcid,bad_bcid;
  vector<double> m_etabins,m_phibins;
  double etalo,etahi;
  MMT_Finder m_find;//finds tracks 
  bool redo;
  string nom;
  int plotcounter;
  //ClassDef(MM_Trigger,0);
};
#endif
