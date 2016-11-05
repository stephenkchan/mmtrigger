#ifndef MMT_FITTER_H
#define MMT_FITTER_H

//C++ language libararies
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <cstdio>

//ROOT libraries
#include "TTree.h"
#include "TH1F.h"
#include "MMT_struct.h"
using namespace std;

class MMT_Fitter{
 public:
  MMT_Fitter(MMT_Parameters *par, int nlg=256, double LG_min=0., double LG_max=0.5);
  ~MMT_Fitter(){}
  void Get_Fit(vector<Hit>& track, map<int,evFit_entry>& Event_Fit, map<hdst_key,hdst_entry>& Hits_Data_Set_Time);
  evFit_entry fit_event(int event, vector<Hit>& track, vector<hdst_entry>& hdsts, int& nfit,vector<pair<double,double> >&mfits) const;
  int get_last() const {return last;}
  int SC_ROI_n_x() const {return m_par->n_x;}
  int SC_ROI_n_y() const {return m_par->n_y;}
/*   vector<int> xent,yent; */
  int find_hdst(const vector<hdst_entry>& hdsts, const hdst_key& key) const;

 protected:
  int last;
  //some variables
  MMT_Parameters *m_par;

  //functions translated
  int Filter_UV(vector<Hit>& track) const;
  float32fixed<2> Get_Global_Slope(const vector<Hit>& track, const string& type) const;
  float32fixed<2> Get_Local_Slope(const vector<Hit>& Track,double theta=-999.,double phi=-999.) const;
  ROI Get_ROI(float32fixed<2> M_x,float32fixed<2> M_u,float32fixed<2> M_v,const vector<Hit>&track,vector<pair<double,double> >&mfits) const;
  double phi_correct_factor(const vector<Hit>&track)const;
  float32fixed<2> Get_Delta_Theta(float32fixed<2> M_local,float32fixed<2> M_global) const;
  float32fixed<2> Get_Delta_Theta_division(float32fixed<2> M_local,float32fixed<2> M_global, float32fixed<4> a=1.) const;
  int Rough_ROI_temp(float32fixed<4> theta, float32fixed<4> phi) const;

  //sim hit code stuff
  int track_to_index(const vector<Hit>&track)const;

  //ideal local slope for debugging purposes
  double ideal_local_slope(const vector<Hit>& Track)const;
  double ideal_z(const Hit& hit)const;
  double ideal_ak(const vector<Hit>& Track)const;
  double ideal_zbar(const vector<Hit>& Track)const;

  //translated from Table_Generators.m
  float32fixed<2> LG_lgr(int ilgr, double a, int number_LG_regions, float32fixed<2> _min, float32fixed<2> _max) const;
  float32fixed<2> mult_factor_lgr(int ilgr, double a, int number_LG_regions, float32fixed<2> _min, float32fixed<2> _max) const;
  float32fixed<4> Slope_Components_ROI_val(int jy, int ix, int thetaphi) const;
  float32fixed<4> Slope_Components_ROI_theta(int jy, int ix) const;
  float32fixed<4> Slope_Components_ROI_phi(int jy, int ix) const;
  float32fixed<2> DT_Factors_val(int i, int j) const;

 private:
  //some functions
  vector<Hit> q_hits(const string& type,const vector<Hit>& hits) const;
  double degtorad(double degree_value) const;
  double radtodeg(double radian_value) const;

  //Fitter components
  int number_LG_regions,n_fit;
  float32fixed<2> LG_min,LG_max;
  vector<int> q_planes(char type) const;//return position of what planes are where
};
#endif
