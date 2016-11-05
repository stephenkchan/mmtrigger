#ifndef MMT_FINDER_H
#define MMT_FINDER_H

//C++ language libararies
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

//ROOT libraries
#include "TTree.h"
#include "TH1F.h"
#include "MMT_struct.h"
#include "MMT_Fitter.h"
using namespace std;

class MMT_Finder : public MMT_Fitter{
 public:
  MMT_Finder(MMT_Parameters *par);
  ~MMT_Finder(){}
  void Check_Coincidence_Gates();
  int Coincidence_Gate(const vector<bool>& plane_hits) const;
  int finder_rd_cgate(const map<pair<int,int>,finder_entry>& finder,int road,vector<Hit>&track)const;
  void set_roads(int _roads) { roads=_roads; }

  int get_roads()const{return roads;}
  map<pair<int,int>,finder_entry> finder_event(int event,vector<hdst_entry>& hdsts,vector<evFit_entry>& fit,vector<pair<double,double> >&mfits,double&mxl) const;
  int Finder_Control(map<hdst_key,hdst_entry>& data_hdst,int BC_window, map<int,evFit_entry>& Event_Fit, map<hdst_key,hdst_entry>& Hits_Data_Set_Time);

 protected:
  //event-wise version function exclusives
  void record_hit_evFinder(map<pair<int,int>,finder_entry>& evFinder, const Hit& hit) const;
  void sort_hdsts(vector<hdst_entry>& hdsts) const;

  //functions translated
  Hit Translate_Hit(const hdst_entry& data) const;
  int Check_Road_for_Coincidence(int road) const;
  vector<Hit> Receive_Hits(map<hdst_key,hdst_entry>& data_hdst) const;
  void Record_Hits_in_Finder(const vector<Hit>& hits_to_write);
  vector<Hit> Read_Finder_to_Track(int road) const;
  void Fit_Gate(int road, map<int,evFit_entry>& Event_Fit, map<hdst_key,hdst_entry>& Hits_Data_Set_Time);
  void Clear_Finder(int road, string type);
  
 private:
  vector<int> q_planes(const string& type) const;
  //Finder components
  double clock,max_age;
  int roads;
  double slope_min,slope_max;
  vector<vector<double> > Gate_Flags; 
  vector<vector<finder_entry> > Finder;
};
#endif
