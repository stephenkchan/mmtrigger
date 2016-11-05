#ifndef MMT_LOADER_H
#define MMT_LOADER_H

//C++ language libararies
#include <dirent.h>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iostream>

//ROOT libraries
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TH2D.h"
#include "Rtypes.h"
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TF1.h"
#include "MMT_struct.h" 
//#include "../plots/mmt_plot/atlasstyle-00-03-05/AtlasStyle.h"
using namespace std;

class MMT_Loader{// : public TObject{
 public:
  MMT_Loader(string fname, par_par var=par_par(), const string& tname="");
  virtual ~MMT_Loader(){ }
  int get_nevents() const {return tr->GetEntries();}

 protected:
  //debug stuff
  map<int,int> zplanes;
  vector<TH1D*> m_diff_xuv;
  //parameters for analysis
  MMT_Parameters *m_par;
  vector<vector<TH2D*> > gposx_strippos_l,glposx_strippos_l;
  vector<vector<TH1D*> > ybases_solved,ybases_pos,ybases_ent,strip_widths;

  //# track tot
  int n_bar,n_ent;

  //load event into containers
  bool load_event(int event);
  vector<hdst_entry> event_hdsts(int find_event) const;
  vector<hdst_key> event_hdst_keys(int find_event) const;

  //pretty plot time
  void print_event(const map<int,hdst_entry>&fit_hits,const TString& evtag)const;
  TLatex plot_latex() const;
  TStyle* set_style() const;

  //bg events
  bool gen_bg;
  vector<hdst_entry> Generate_Incoherent_Background(const vector<hdst_key>& keys);
/*   vector<hdst_entry> blc_Generate_Incoherent_Background(const vector<hdst_key>& keys); */
/*   void Build_Occupancy_Table(); */

  //containers and maintenance
  void clear_buffers(int event);
  map<int,evFit_entry> Event_Fit;//key is event no.
  map<int,evInf_entry> Event_Info;//key is event no.
  map<int,evAna_entry> Event_Analysis;//key is event no.
  map<hdst_key,hdst_entry> Hits_Data_Set_Time;//key is hit_index? <BC_time,time>?

  //VMM info
  vector<vector<vector<bool> > > VMM_chip_status;
  vector<vector<vector<double> > > VMM__chip_last_hit_time;
  void reset_VMMs();
  bool Mimic_VMM_Chip_Deadtime(hdst_entry& candy);
  double VMM_deadtime;
  int num_VMM_per_plane;
  string file_name;

 private:
  //incoherent bg stuff
  TRandom3 *m_rand;
  TF1 *hit_rate,*hit_ct_rand;
  double hit_rate_GHz,tlo,thi, wedgex,wedgey;
  vector<vector<double> > Occupancy_Table;

  //x <---> u/v switches
  bool uvxxmod;
  void xxuv_to_uvxx(TVector3& hit,int plane)const;
  void hit_rot_stereo_fwd(TVector3& hit)const;//x to v, u to x
  void hit_rot_stereo_bck(TVector3& hit)const;//x to u, v to x

  //calculate misalignment hit coordinates; takes a truth hit and translates
  //to misaligned wedge-centric coordinates
  hdst_info misaligned_hit(double x, double y, double z, int plane) const{}

  //Import_Athena..._.m stuff
  double phi_shift(double athena_phi) const;
  int Get_VMM_chip(int strip) const;  //*** Not Finished... Rough
  int strip_number(int station,int plane,int spos)const;
  int Get_Strip_ID(double X,double Y,int plane,int&station) const;

  //NSWHitsTree stuff
  void Init(TTree *t);
  void GetEntry(int i);

  TTree *tr;
  TFile *f;

  //branch variables
  UInt_t          runNumber,eventNumber,TruthVertex_n,TruthParticle_n;
  vector<int>     *TruthParticle_Pdg;
  vector<double>  *TruthVertex_X,*TruthVertex_Y,*TruthVertex_Z;
  vector<double>  *TruthParticle_Pt,*TruthParticle_Eta,*TruthParticle_Phi,*TruthParticle_E;
  vector<int>     *TruthParticle_Barcode;
  UInt_t          MuEntry_Particle_n;
  vector<double>  *MuEntry_Particle_Pt,*MuEntry_Particle_Eta,*MuEntry_Particle_Phi;
  vector<double>  *MuEntry_Position_Eta,*MuEntry_Position_Phi;
  UInt_t          Hits_MM_n;
  vector<double>  *Hits_MM_globalTime;
  vector<double>   *Hits_MM_hitGlobalPositionX,*Hits_MM_hitGlobalPositionY,*Hits_MM_hitGlobalPositionZ;
  UInt_t          Digits_MM;
  vector<int>     *Digits_MM_multiplet,*Digits_MM_gas_gap;
  vector<vector<float> > *Digits_MM_time, *Digits_MM_charge;
  vector<vector<int> > *Digits_MM_stripPosition;
  vector<int> *Digits_MM_stationEta,*Digits_MM_stationPhi;
  vector<string> *Digits_MM_stationName;
  vector<vector<double> > *Digits_MM_stripLposX,*Digits_MM_stripLposY;
  vector<vector<double> > *Digits_MM_stripGposX,*Digits_MM_stripGposY,*Digits_MM_stripGposZ;
  vector<int>     *Digits_MM_truth_barcode;
  vector<double>  *Digits_MM_truth_localPosX,*Digits_MM_truth_localPosY;
  vector<float>   *Digits_MM_truth_XZ_angle;
  vector<int>     *Digits_MM_stripForTrigger;
  vector<float>   *Digits_MM_stripTimeForTrigger;
  //the actual branches
  TBranch *b_runNumber,*b_eventNumber,*b_TruthVertex_n,*b_TruthParticle_n,*b_TruthVertex_X,*b_TruthVertex_Y,*b_TruthVertex_Z;
  TBranch *b_TruthParticle_Pt,*b_TruthParticle_Eta,*b_TruthParticle_Phi,*b_TruthParticle_E,*b_TruthParticle_Pdg;
  TBranch *b_TruthParticle_Barcode,*b_MuEntry_Particle_n,*b_MuEntry_Particle_Pt,*b_MuEntry_Particle_Eta,*b_MuEntry_Particle_Phi;
  TBranch *b_MuEntry_Position_Eta,*b_MuEntry_Position_Phi,*b_Hits_MM_n,*b_Hits_MM_globalTime,*b_Digits_MM,*b_Digits_MM_multiplet,*b_Digits_MM_gas_gap;
  TBranch *b_Hits_MM_hitGlobalPositionX,*b_Hits_MM_hitGlobalPositionY,*b_Hits_MM_hitGlobalPositionZ;
  TBranch *b_Digits_MM_time,*b_Digits_MM_charge,*b_Digits_MM_stripPosition,*b_Digits_MM_stripLposX,*b_Digits_MM_stripLposY;
  TBranch *b_Digits_MM_stationEta,*b_Digits_MM_stationPhi,*b_Digits_MM_stationName;
  TBranch *b_Digits_MM_stripGposX,*b_Digits_MM_stripGposY,*b_Digits_MM_stripGposZ,*b_Digits_MM_truth_barcode,*b_Digits_MM_truth_localPosX,*b_Digits_MM_truth_localPosY;
  TBranch *b_Digits_MM_truth_XZ_angle,*b_Digits_MM_stripForTrigger,*b_Digits_MM_stripTimeForTrigger;

  //ClassDef(MMT_Loader,0); 
};

#endif
