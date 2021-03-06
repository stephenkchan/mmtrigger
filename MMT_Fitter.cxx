#include "MMT_Fitter.h"
const bool squack=debug&&false;
const bool div_hack=false;
const bool circumvent=false;

MMT_Fitter::MMT_Fitter(MMT_Parameters *par, int nlg, double lgmin, double lgmax): /*m_par(par),*/ number_LG_regions(nlg), LG_min(lgmin), LG_max(lgmax){
  if(squack) cout << "MMT_F::building fitter"<<endl;
  m_par=par;
  last=0;
  n_fit=0;
//   if(m_par->correct.type==2){m_x_min*=cos(m_par->correct.rotate.Y());m_x_min*=cos(m_par->correct.rotate.Y());}
  //use floats, because we are initializing
  if(squack) cout << "MMT_F::built fitter"<<endl;
}

void MMT_Fitter::Get_Fit(vector<Hit>& track, map<int,evFit_entry>& Event_Fit, map<hdst_key,hdst_entry>& Hits_Data_Set_Time){
  bool fit=false, Generate_FPGA_Sample=false; //the latter sometimes used when debugging

  //@@@@@@@@@@ Begin Track Fitting @@@@@@@@@@@@
  //----- UV filter --------------
  Filter_UV(track);

  //---- Calc global slopes and local X slope -----
  float32fixed<2> M_x_global=Get_Global_Slope(track,"x"),M_u_global=Get_Global_Slope(track,"u"),M_v_global=Get_Global_Slope(track,"v"),M_x_local=Get_Local_Slope(track);

  //---------Debugging------------
  /*
  double LG=M_x_local*M_x_global;
  if(LG>LG_max) LG_max=LG;
  if(LG<LG_min) LG_min=LG;
  */

  //----  Calc delta theta ----------
  float32fixed<2> Delta_Theta_division = Get_Delta_Theta_division(M_x_local,M_x_global,1.);//_division?
  float32fixed<2> Delta_Theta = Get_Delta_Theta(M_x_local,M_x_global);
  //----- Calc ROI ----------
  vector<pair<double,double> >dummy;
  ROI ROI = Get_ROI(M_x_global,M_u_global,M_v_global,track,dummy);

 //----- Abandon fit if ROI comes back as out of bounds ------
  if(ROI.theta==-999){
    cerr<<"SOMETHING IS OFF!  Get_Fit\n";
    exit(-999);
  }

  //@@@@@@@@ Begin Info Storage for Later Analysis @@@@@@@@@@@@@@@@@@@@@@@
  vector<int> xpl=m_par->q_planes("x"),upl=m_par->q_planes("u"),vpl=m_par->q_planes("v");
  vector<bool> planes_hit_tr(8,false),planes_hit_bg(8,false);
  int event=-1;//*** ASK ABOUT THIS!!!
  for(unsigned int j=0; j<m_par->setup.size(); j++){
    hdst_key key=track[j].key;
    if(Hits_Data_Set_Time.find(key)!=Hits_Data_Set_Time.end()) event = Hits_Data_Set_Time[key].event;
    if(Hits_Data_Set_Time[key].truth_nbg) planes_hit_tr[j]=true;
    else planes_hit_bg[j]=true;
  }

  int n_xpl_tr=0,n_xpl_bg=0,n_uvpl_tr=0,n_uvpl_bg=0;
  for(unsigned int xp=0; xp<xpl.size(); xp++){
    int plane=xpl[xp];
    n_xpl_tr+=planes_hit_tr[plane];
    n_xpl_bg+=planes_hit_bg[plane];
  }
  for(unsigned int up=0; up<upl.size(); up++){
    int plane=upl[up];
    n_uvpl_tr+=planes_hit_tr[plane];
    n_uvpl_bg+=planes_hit_bg[plane];
  }
  for(unsigned int vp=0; vp<vpl.size(); vp++){
    int plane=vpl[vp];
    n_uvpl_tr+=planes_hit_tr[plane];
    n_uvpl_bg+=planes_hit_bg[plane];
  }
  //FINISH ME!!!!!!!
  Event_Fit[event].X_hits_in_fit=n_xpl_tr;
  Event_Fit[event].UV_hits_in_fit=n_uvpl_tr;
  Event_Fit[event].bg_X_fit=n_xpl_bg;
  Event_Fit[event].bg_UV_fit=n_uvpl_bg;
  Event_Fit[event].fit_dtheta=Delta_Theta.getValue();

  for(unsigned int plane=0; plane<m_par->setup.size(); plane++){
    if(ROI.theta==-999 || Delta_Theta==-999) continue; //&& Delta_Theta_division~=-999
    hdst_key key=track[plane].key;
    int true_hit_event = Hits_Data_Set_Time[key].event;
    float32fixed<4> track_fit_theta = ROI.theta, track_fit_phi = ROI.phi;
    if(Generate_FPGA_Sample) Hits_Data_Set_Time[key].fit_fill(ROI.theta,ROI.phi,Delta_Theta,M_x_global,M_u_global,M_v_global,M_x_local,ROI.m_x,ROI.m_y,ROI.roi);
    else Hits_Data_Set_Time[key].fit_fill(ROI.theta,ROI.phi,Delta_Theta);
    fit = true;
    int roi=(Generate_FPGA_Sample?ROI.roi:0);
    Event_Fit[true_hit_event].fit_theta=track_fit_theta;
    Event_Fit[true_hit_event].fit_phi=track_fit_phi;
    Event_Fit[true_hit_event].fit_dtheta=Delta_Theta_division;
    Event_Fit[true_hit_event].fit_roi=roi;
    
    if(Hits_Data_Set_Time[key].truth_nbg) Event_Fit[true_hit_event].truth_planes_hit+=pow(10,plane);
    else Event_Fit[true_hit_event].bg_planes_hit+=pow(10,plane);
  }
  if(fit) n_fit++;
}

evFit_entry MMT_Fitter::fit_event(int event, vector<Hit>& track, vector<hdst_entry>& hdsts, int& nfit,vector<pair<double,double> >&mfits) const{
  if(debug)cout<<"Begin fit event!"<<endl;
  bool did_fit=false;
  int check=Filter_UV(track);
  vector<int> xpl=m_par->q_planes("x"),upl=m_par->q_planes("u"),vpl=m_par->q_planes("v");
  //---- Calc global slopes and local X slope -----
  float32fixed<2> M_x_global = Get_Global_Slope(track,"x"),M_u_global = (check>=10?float32fixed<2>(-999.):Get_Global_Slope(track,"u")),M_v_global = (check%10==1?float32fixed<2>(-999.):Get_Global_Slope(track,"v"));
  //----  Calc delta theta ----------
  //----- Calc ROI ----------
  ROI ROI = Get_ROI(M_x_global,M_u_global,M_v_global,track,mfits);
//   if(ROI.roi<1) cout <<"WHaaaaa!? roi of "<<ROI.roi<<endl;
  if(ROI.theta==-999){
    for(unsigned int i=0;i<track.size();i++)  track[i].print(); 
    cerr<<"SOMETHING IS OFF!  fit_event\n";
    exit(-999);
  }
  
  float32fixed<2> M_x_local = Get_Local_Slope(track,ROI.theta.getFloat(),ROI.phi.getFloat()),Delta_Theta_division = Get_Delta_Theta_division(M_x_local,M_x_global,1.), Delta_Theta = Get_Delta_Theta(M_x_local,M_x_global), dtheta_idl=Get_Delta_Theta_division(ideal_local_slope(track),M_x_global);
  if(debug&&abs(dtheta_idl-Delta_Theta_division)>2.e-3)m_par->fill0=true;
  if(debug)cout<<"Mxg="<<M_x_global.getValue()<<",Mug="<<M_u_global.getValue()<<",Mvg="<<M_v_global.getValue()<<",Mxl="<<M_x_local.getValue()<<",dth="<<Delta_Theta.getValue()<<endl;
  //@@@@@@@@ Begin Info Storage for Later Analysis @@@@@@@@@@@@@@@@@@@@@@@
  vector<bool> planes_hit_tr(8,false),planes_hit_bg(8,false);
  for(unsigned int ihit=0; ihit<track.size(); ihit++){
    int hdst_pos=find_hdst(hdsts,track[ihit].key);
    if(hdst_pos==-1) continue;
    if(hdsts[hdst_pos].truth_nbg) planes_hit_tr[track[ihit].info.plane]=true;
    else planes_hit_bg[track[ihit].info.plane]=true;
  }
  int n_xpl_tr=0,n_xpl_bg=0,n_uvpl_tr=0,n_uvpl_bg=0;
  for(int xp=0; xp<(int)xpl.size(); xp++){
    int plane=xpl[xp];
    n_xpl_tr+=planes_hit_tr[plane];n_xpl_bg+=planes_hit_bg[plane];
  }
  if(check<10){
    for(unsigned int up=0; up<upl.size(); up++){
      int plane=upl[up];
      n_uvpl_tr+=planes_hit_tr[plane];
      n_uvpl_bg+=planes_hit_bg[plane];
    }
  }
  if(check%10==0){
    for(unsigned int vp=0; vp<vpl.size(); vp++){
      int plane=vpl[vp];
      n_uvpl_tr+=planes_hit_tr[plane];
      n_uvpl_bg+=planes_hit_bg[plane];
    }
  }
//   return evFit_entry(event,ROI.theta,ROI.phi,Delta_Theta_division,ROI.roi,n_xpl_tr,n_uvpl_tr,n_xpl_bg,n_uvpl_bg,Delta_Theta);

  //FINISH ME!!!!!!!
  float32fixed<4> candtheta=ROI.theta,candphi=ROI.phi;
  /* I think this bit appears redundant but could end up being stupid; the fitter shouldn't care about CT stuff (beyond having min num hits to do fit, which the -999 or w/e is responsible for taking care of)
  bool xfail=(n_xpl_tr+n_xpl_bg<m_par->CT_x),uvfail= (n_uvpl_tr+n_uvpl_bg<m_par->CT_uv);
  if(debug)cout<<n_xpl_tr+n_xpl_bg<<" x hits"<<endl;
  if(xfail)candtheta=-999.;
  if(uvfail)candphi=-999.;
  */
  bool fitkill=(ROI.theta==-999 || Delta_Theta==-999||Delta_Theta==-4);// ||xfail||uvfail);
  if(debug)cout<<"HIT CODE: "<<track_to_index(track)<<endl;
  evFit_entry aemon(event,candtheta,candphi,Delta_Theta_division,ROI.roi,n_xpl_tr,n_uvpl_tr,n_xpl_bg,n_uvpl_bg,Delta_Theta,track_to_index(track));
  if(fitkill) return aemon;
//   if(debug)cout<<" aemon fit_theta="<<aemon.fit_theta.getFloat()<<", fit_phi="<<aemon.fit_phi.getFloat()<<"...";
  int nplanes=m_par->setup.size();
  for(int plane=0; plane<nplanes; plane++){
    if(track[plane].info.slope==-999) continue; //&& Delta_Theta_division~=-999
    int hdst_pos=find_hdst(hdsts,track[plane].key);
    if(hdst_pos==-1) continue;
    did_fit=true;
    hdsts[hdst_pos].fit_fill(ROI.theta,ROI.phi,Delta_Theta,M_x_global,M_u_global,M_v_global,M_x_local,ROI.m_x,ROI.m_y,ROI.roi);
    aemon.fit_hit_keys.push_back(track[plane].key);
//     if(debug)cout<<"hdst fit_theta="<<hdsts[hdst_pos].fit_theta<<"...";
    if(hdsts[hdst_pos].truth_nbg) aemon.truth_planes_hit+=pow(10,nplanes-plane-1);
    else aemon.bg_planes_hit+=pow(10,nplanes-plane-1);
  }
  if(did_fit) nfit++;
//   if(debug)cout<<"aemon has "<<aemon.fit_hit_keys.size()<<" keys"<<endl;
//   if(event==2997||event==2800) aemon.print();
  return aemon;
}

int MMT_Fitter::find_hdst(const vector<hdst_entry>& hdsts, const hdst_key& key) const{
  for(unsigned int i=0;i<hdsts.size();i++){
    if(hdsts[i].BC_time==key.BC_time&&hdsts[i].time==key.time&&hdsts[i].gtime==key.gtime&&hdsts[i].VMM_chip==key.VMM_chip) return i;
  }
  return -1;
}


int MMT_Fitter::Filter_UV(vector<Hit>& track) const{
  return 0;
  float32fixed<2> h=m_par->h, tolerance = h;//*2;  //Can be optimized...
  vector<int> u_planes=m_par->q_planes("u"),v_planes=m_par->q_planes("v");
  vector<Hit> u_hits=q_hits("u",track),v_hits=q_hits("v",track);
  bool pass_u=!u_hits.empty(),pass_v=!v_hits.empty();
  //if the difference in slope between the first and last u/v planes is too great don't pass, set track hits to zero
  if(pass_u){
//     if(squack)cout<<"diff="<<fabs(u_hits.front().slope-u_hits.back().slope)<<" should be |"<<u_hits.front().slope<<"-"<<u_hits.back().slope<<"|....tolerance is"<<tolerance<<endl;
    if((float32fixed<2>(u_hits.front().info.slope-u_hits.back().info.slope)).fabs()>tolerance) pass_u=false;
  }
  if(pass_v){
    if((float32fixed<2>(v_hits.front().info.slope-v_hits.back().info.slope)).fabs()>tolerance) pass_v=false;
  }
  int return_val=0;//of form (bool ubad, bool vbad), so 10 would be bad u, good v
  if(!pass_u){
    for(unsigned int iup=0;iup<u_planes.size();iup++) track[u_planes[iup]].info.slope=-999;
    return_val+=10;
  }
  if(!pass_v){
    for(unsigned int ivp=0;ivp<v_planes.size();ivp++) track[v_planes[ivp]].info.slope=-999;
    return_val+=1;
  }
  return return_val;
}

float32fixed<2> MMT_Fitter::Get_Global_Slope(const vector<Hit>& track, const string& type) const{
  vector<Hit> qhits=q_hits(type,track);
  float32fixed<2> sum=0.;
  if(qhits.size()==0)return -999;
  float32fixed<2> nhitdiv=1./qhits.size();
//   if(debug)cout<<"Calculating global slope of type "<<type<<"---adding "<< qhits.size() <<" slopes: ";
  for(int ihit=0;ihit<(int)qhits.size();ihit++){
    sum+=(qhits[ihit].info.slope*nhitdiv);
//     if(debug)cout<<qhits[ihit].info.slope.getValue()<<" ";
//     if(debug)cout<<"("<<qhits[ihit].info.slope.getValue()<<" x "<<nhitdiv.getValue()<<" = "<<float32fixed<2>(qhits[ihit].info.slope*nhitdiv).getValue()<<"; sum is now "<<sum.getValue()<<" ) ";
  }
//   if(debug) cout<<" for a slope of "<<sum.getValue()<<endl;
  return sum;
}

//CHANGE!
float32fixed<2> MMT_Fitter::Get_Local_Slope(const vector<Hit>& Track,double theta,double phi)const{
  vector<int> x_planes=m_par->q_planes("x"),ybin_hits(x_planes.size(),-1);
  int nxp=x_planes.size();
  for(int ipl=0; ipl<nxp; ipl++) ybin_hits[ipl]=((Track[x_planes[ipl]].info.slope==-999||Track[x_planes[ipl]].info.slope==-4)?-1:m_par->ybin(Track[x_planes[ipl]].info.y));
  bool hit=false;
  float32fixed<yzdex> yzsum=0;
  float mxlf=0;
  int xdex,ybin,which;m_par->key_to_indices(ybin_hits,xdex,ybin,which);
  if(xdex<0||ybin<0||which<0) return -999;
  float32fixed<zbardex> zbar=m_par->Ak_local_slim[xdex][ybin][which];
  float32fixed<bkdex>bk=m_par->Bk_local_slim[xdex][ybin][which];
  if(debug)cout<<"zbar is "<<zbar.getValue()<<", and bk is "<<bk.getValue()<<endl;
  int ebin=m_par->eta_bin(theta),pbin=m_par->phi_bin(phi);
  for(int ipl=0; ipl<nxp; ipl++){
    float32fixed<yzdex> z=Track[x_planes[ipl]].info.z,y=Track[x_planes[ipl]].info.y;
    if(ebin!=-1){
//       cout<<"Before....y="<<y.getFloat()<<"..."<<endl;
//        double zflt=m_par->z_nominal[x_planes[ipl]].getFloat(),x=zflt*tan(theta)*sin(phi),yup=zflt*tan(theta)*cos(phi)-m_par->ybases[x_planes[ipl]][0].getFloat();
//        double alpha=m_par->correct.rotate.Z(),stuff=((cos(alpha)-1.)*yup+x*sin(alpha))/store_const();
//       cout<<"("<<m_par->ymod[ebin][pbin][ipl].getFloat()<<","<<stuff<<")"<<endl;
      z=z+m_par->zmod[ebin][pbin][ipl];
      y=y+m_par->ymod[ebin][pbin][ipl];
//       cout<<"After....y="<<y.getFloat()<<",dy="<<y-Track[x_planes[ipl]].info.y<<"..."<<endl;
    }
//     cout<<theta<<"("<<ebin<<"--"<<2*atan(exp(-1.*m_par->m_etabins[ebin]))<<") "<<phi<<"("<<pbin<<"--"<<m_par->m_phibins[pbin]<<") "<<m_par->z_nominal[x_planes[ipl]].getFloat()<<" "<<m_par->ybases[x_planes[ipl]][0].getFloat()<<" "<<(ebin==-1?0.:m_par->ymod[ebin][pbin][ipl].getFloat())*store_const()<<endl;
    if(Track[x_planes[ipl]].info.slope==-999||Track[x_planes[ipl]].info.slope==-4) continue;
    hit=true;
    yzsum+=y*(z*zbar-1.);
    mxlf += bk.getFloat()*y.getFloat()*(z.getFloat()*zbar.getFloat()-1.);
//     if(debug)cout<<"yzsum="<<yzsum.getValue()<<", mxlf="<<mxlf<<endl;
  }
  float32fixed<2> mxl=float32fixed<2>(bk.getValue()*yzsum.getValue());
  if(!hit) {return float32fixed<2>(999);}
  if(log10(abs(mxl.getValue()-mxlf))>-3.&&xdex!=5&&xdex!=10)cout<<setprecision(20)<<"*****AVENGE ME! *******fixed: "<<mxl.getValue()<<", float:"<<mxlf<<",    (ERROR: "<<mxlf-mxl.getValue()<<")"<<endl;
  return mxl;
}

int MMT_Fitter::track_to_index(const vector<Hit>&track)const{
  vector<bool>hits(m_par->setup.size(),false);
  for(int ihit=0;ihit<(int)track.size();ihit++)hits[track[ihit].info.plane]=(hits[track[ihit].info.plane]?true:track[ihit].info.slope>-2.);
  return m_par->bool_to_index(hits);
}

double MMT_Fitter::ideal_local_slope(const vector<Hit>& Track)const{
  vector<vector<double> > z_hit;
  for(int i = 0; i<m_par->z_large.size(); i++){
    vector<double> temp;
    for(int j = 0; j<m_par->z_large[i].size(); j++)
      temp.push_back(m_par->z_large[i][j].getFloat());
    z_hit.push_back(temp);
  }
  vector<int> x_planes=m_par->q_planes("x");
  int nxp=x_planes.size();
  bool hit=false;
  double sum_xy=0,sum_y=0;
  double ak_idl=ideal_ak(Track),bk_idl=ak_idl*ideal_zbar(Track);
  for(int ipl=0; ipl<nxp; ipl++){
    double y=Track[x_planes[ipl]].info.y.getFloat(),z=Track[x_planes[ipl]].info.z.getFloat();
    if(y==-999) continue;
    hit=true;
    sum_xy += ak_idl*z*y;
    sum_y  += bk_idl*y;
//     if(debug)cout<<"...z="<<z<<",y="<<y<<",sum_y="<<sum_y<<",sum_xy="<<sum_xy;
  }
  if(!hit) return -10.;
  double ls_idl=sum_xy-sum_y;
//   if(debug)cout<<endl<<"ak_idl="<<ak_idl<<",bk_idl="<<bk_idl<<",sum_xy"<<sum_xy<<",sum_y="<<sum_y<<",ls_idl="<<ls_idl<<endl;
  return ls_idl;
}

double MMT_Fitter::ideal_z(const Hit& hit)const{
  int plane=hit.info.plane; 
  double tilt=(plane<4?m_par->correct.rotate.X():0),dz=(plane<4?m_par->correct.translate.Z():0),
    nominal=m_par->z_nominal[plane].getFloat(),y=hit.info.y.getFloat()-m_par->ybases[plane].front().getFloat(),z=nominal+dz+y*tan(tilt);
  return z;
}

double MMT_Fitter::ideal_ak(const vector<Hit>& Track)const{
  vector<int> x_planes=m_par->q_planes("x");//this tells us which planes are x planes
  int n_xplanes=x_planes.size(),hits=0;
  double sum_x=0,sum_xx=0;
  for(int ip=0; ip<n_xplanes; ip++){
    if(Track[x_planes[ip]].info.slope==-999) continue;
    double addme=ideal_z(Track[x_planes[ip]]);//z_hit[x_planes[plane]];
    hits++;//there's a hit
    sum_x  += addme;
    sum_xx += addme*addme;
    if(squack) cout<<"z["<<ip<<"]="<<addme<<", sum_x="<<sum_x<<", sum_xx="<<sum_xx<<", hits="<<hits<<endl;
  }
  double diff = hits*sum_xx-sum_x*sum_x;
  return hits/diff;
}

double MMT_Fitter::ideal_zbar(const vector<Hit>& Track)const{
  vector<int> x_planes=m_par->q_planes("x");//this tells us which planes are x planes
  double ztot=0;int nhit=0;
  for(unsigned int ip=0;ip<x_planes.size();ip++){
    if(Track[x_planes[ip]].info.slope==-999)continue;
    nhit++;ztot+=ideal_z(Track[x_planes[ip]]);
  }
  return ztot/nhit;
}

float32fixed<2> MMT_Fitter::Get_Delta_Theta(float32fixed<2> M_local, float32fixed<2> M_global) const{
  int region=-1;
  if(div_hack)return Get_Delta_Theta_division(M_local,M_global);
  float32fixed<2> LG = M_local * M_global;
  for(int j=0;j<number_LG_regions;j++){   //number_LG_regions
    if(LG <= DT_Factors_val(j,0)){
      region = j;
      break;
    }
  }
//   if(debug) cout<<"m_l="<<M_local.getFloat() <<", m_g="<<M_global.getFloat() <<", for dtheta="<< DT_Factors_val(region,1)*(M_local - M_global) <<endl;
  if(region==-1) return -999;
  return DT_Factors_val(region,1)*(M_local - M_global);
}

float32fixed<2> MMT_Fitter::DT_Factors_val(int i, int j) const{
  if(m_par->val_tbl){
    return m_par->DT_Factors[i][j];
  }
  if(j<0||j>1){
    cerr<<"DT_Factors only has two entries on the second index (for LG and mult_factor); you inputed an index of " << j << endl;
    exit(1);
  }
  if(i<0||i>=number_LG_regions){
    cerr<<"There are " << number_LG_regions << " in DT_Factors(_val); you inputed an index of " << i << endl;
    exit(1);
  }
  double a=1.;//not sure what this is for, so hard to choose fixed_point algebra
  if(j==0) return mult_factor_lgr(i,a,number_LG_regions,LG_min,LG_max);
  return LG_lgr(i,a,number_LG_regions,LG_min,LG_max);
}

float32fixed<2> MMT_Fitter::LG_lgr(int ilgr, double a, int number_LG_regions, float32fixed<2> _min, float32fixed<2> _max) const{
  a+=0;
  return _min+float32fixed<2>(ilgr/number_LG_regions)*(_max-_min);
}

float32fixed<2> MMT_Fitter::mult_factor_lgr(int ilgr, double a, int number_LG_regions, float32fixed<2> _min, float32fixed<2> _max) const{
  return float32fixed<2>(1./(a+LG_lgr(ilgr,a,number_LG_regions,_min,_max)/a));
}

float32fixed<2> MMT_Fitter::Get_Delta_Theta_division(float32fixed<2> M_local, float32fixed<2> M_global, float32fixed<4> a) const{
//delta_theta = theta_local_slope - theta_global_fit
// a=1;  // Really sin(phi), but I use small angles about phi=pi/2
//   if(debug) cout<<"M_l="<<M_local.getFloat()<<", M_g="<<M_global.getFloat()<<", for dtheta="<<(M_local.getValue() - M_global.getValue())/(a + (M_local.getValue()*M_global.getValue())/a.getFloat())<<endl;

  //we could use 2 bits for the numerator and 3 for the denominator, but then
  //fixed_point doesn't know how to do the algebra. Assume we know how to do
  //this division (we don't, efficiently, thus the method Get_Delta_Theta
  return float32fixed<2>((M_local - M_global)/(a + (M_local*M_global)/a.getFloat()));
}

vector<Hit> MMT_Fitter::q_hits(const string& type,const vector<Hit>& track) const{
  string setup(m_par->setup);
  if(setup.length()!=track.size()){
    cerr<<"Setup has length: "<<setup.length()<<", but there are "<<track.size()<<" hits in the track"<<endl;
    exit(2);
  }
  vector<int> qpl(m_par->q_planes(type));
  vector<Hit> q_hits;
  for(unsigned int ihit=0; ihit<qpl.size(); ihit++){
    if(track[qpl[ihit]].info.slope!=-999) q_hits.push_back(track[qpl[ihit]]);
  }
  return q_hits;
}

double MMT_Fitter::degtorad(double degree_value) const{
  return atan(1.)/45.*degree_value;
}
  
double MMT_Fitter::radtodeg(double radian_value) const{
  return 45./atan(1.)*radian_value;
}

//change this to take u and/or v out of the roi calculation
ROI MMT_Fitter::Get_ROI(float32fixed<2> M_x,float32fixed<2> M_u,float32fixed<2> M_v,const vector<Hit>&track,vector<pair<double,double> >&mfits) const{
  //M_* are all global slopes
  if(debug) cout<<"\nGet_ROI("<<M_x.getValue()<<","<<M_u.getValue()<<","<<M_v.getValue()<<"); ";

  //--- calc constants ------
  float32fixed<2> b=degtorad(m_par->stereo_degree.getFloat());
  float32fixed<7> A=1./tan(b.getFloat()),B=1./tan(b.getFloat());

  //---  slope conversion equations ---- 
  float32fixed<2> m_y = M_x;
  float32fixed<2> m_xu = A*M_u.getFloat() - B*m_y.getFloat(),m_xv = B*m_y.getFloat() - A*M_v.getFloat();

  //--- which slopes are truly present ----  
  //Note that bad slopes are not necessarily 0 as I often use -999 to denote something missing
  //we have -999 for M_u or M_v to denote that it didn't pass filtering
  int nu=1,nv=1;
  if(M_u<0||M_u==float32fixed<2>(-999)){
    m_xu = 0;nu=0;
  }
  if(M_v<0||M_v==float32fixed<2>(-999)){
    m_xv=0;nv=0;
  }
  if(nu==0&&nv==0) return ROI(-999,-999,-999,-999,-999);
  //--- average of 2 mx slope values ----if both u and v were bad, give it a -999 value to know not to use m_x
  //*** check to see if U and V are necessary for fit
  float32fixed<2> m_x = (nu+nv==0?0:(m_xv+m_xu)/(nu+nv));
  if(m_par->correct.translate.X()!=0&&m_par->correct.type==2){
    m_x+=phi_correct_factor(track)*m_par->correct.translate.X()/m_par->z_nominal[3].getFloat();
  }
  if(debug) cout<<"(b,A,B,my,mxu,mxv,mx)=("<<b.getFloat()<<","<<A.getFloat()<<","<<B.getFloat()<<","<<m_y.getFloat()<<","<<m_xu.getFloat()<<","<<m_xv.getFloat()<<","<<m_x.getValue()<<")\n";
  mfits.push_back(pair<double,double>(m_x.getFloat(),m_y.getFloat()));

  //Get m_x and m_y in parameterized values
  int a_x = round((m_x.getValue()-m_par->m_x_min.getValue())/m_par->h_mx.getValue()), a_y = round((m_y.getValue()-m_par->m_y_min.getValue())/m_par->h_my.getValue());
  // Generally, this offers a reality check or cut.  The only reason a slope
  // should be "out of bounds" is because it represents a weird UV combination
  // -- ie. highly background influenced
  if(a_y>m_par->n_y || a_y<0){
    cout <<"y slope (theta) out of bounds in Get_ROI....(a_x,a_y,m_par->n_x,m_par->n_y)=("<<a_x<<","<<a_y<<","<<m_par->n_x<<","<<m_par->n_y<<")"<<endl;
    return ROI(-999,-999,-999,-999,-999);
  }
  if(a_x>m_par->n_x || a_x<0){
    cout <<"x slope (phi) out of bounds in Get_ROI....(a_x,a_y,m_par->n_x,m_par->n_y)=("<<a_x<<","<<a_y<<","<<m_par->n_x<<","<<m_par->n_y<<")"<<endl;
    return ROI(-999,-999,-999,-999,-999);
  }
//   xent.push_back(a_x);yent.push_back(a_y);
  if(debug)cout<<"fv_angles...(a_x,a_y)=("<<a_x<<","<<a_y<<")"<<endl;
  double phicor=0.;
  if(m_par->correct.rotate.Z()!=0&&m_par->correct.type==2){
//     phicor=-1.*phi_correct_factor(track)*m_par->correct.rotate.Z();
    phicor=-0.2*m_par->correct.rotate.Z();
  }
  float32fixed<4> fv_theta=Slope_Components_ROI_theta(a_y,a_x), fv_phi=(m_x.getValue()==0?-999:Slope_Components_ROI_phi(a_y,a_x).getFloat()+phicor);
  if(debug)cout<<"fv_theta="<<fv_theta.getValue()<<", fv_phi="<<fv_phi.getValue()<<endl;

  //--- More hardware realistic approach but need fine tuning ----
  int roi = Rough_ROI_temp(fv_theta,fv_phi);

  //--- current "roi" which is not an actual roi but an approx phi and theta
  return ROI(fv_theta.getValue(),fv_phi.getValue(),m_x.getValue(),m_y.getValue(),roi);
}

double MMT_Fitter::phi_correct_factor(const vector<Hit>&track)const{
  if((m_par->correct.rotate.Z()==0&&m_par->correct.translate.X()==0)||m_par->correct.type!=2)return 0.;
  int nxmis=0,nx=0,numis=0,nu=0,nvmis=0,nv=0;
  double xpart=0.5,upart=0.5,vpart=0.5;
  string set=m_par->setup;
  for(int ihit=0;ihit<(int)track.size();ihit++){
    int n_pln=track[ihit].info.plane;
    bool ismis=n_pln<4;
    char pln=set[n_pln];
    if(pln=='x'||pln=='X'){nx++;nxmis+=ismis;}
    else if(pln=='u'||pln=='U'){nu++;numis+=ismis;}
    else if(pln=='v'||pln=='V'){nv++;nvmis+=ismis;}
  }
  if(nu==0&&nv==0)return 0.;
  if(nu==0)upart=0.;
  else if(nv==0)vpart=0.;
  else xpart=0.;
  return xpart*1.*nxmis/nx+upart*1.*numis/nu+vpart*1.*nvmis/nv;
}

float32fixed<4> MMT_Fitter::Slope_Components_ROI_val(int jy, int ix, int thetaphi) const{
  if(m_par->val_tbl){
    return m_par->Slope_to_ROI[jy][ix][thetaphi];
  }
  if(thetaphi<0||thetaphi>1){
    cerr<<"Slope_Components_ROI only has two entries on the third index (for theta and phi); you inputed an index of " << thetaphi << endl;
    exit(2);
  }
  if(thetaphi==0) return Slope_Components_ROI_theta(jy,ix);
  return Slope_Components_ROI_phi(jy,ix);
}

float32fixed<4> MMT_Fitter::Slope_Components_ROI_theta(int jy, int ix) const{
  //get some parameter information
  if(jy<0||jy>=m_par->n_y){
    cerr << "You picked a y slope road index of " << jy << " in Slope_Components_ROI_theta; there are only " << m_par->n_y << " of these.\n";
    if(jy>=m_par->n_y)jy=m_par->n_y-1;
    else jy=0;
//     exit(2);
  }
  if(ix<0||ix>=m_par->n_x){
    cerr << "You picked an x slope road index of " << ix << " in Slope_Components_ROI_theta; there are only " << m_par->n_x << " of these.\n";
    if(ix>=m_par->n_x)ix=m_par->n_x-1;
    else ix=0;
//     exit(2);
  }
  int xdex=ix,ydex=jy+1;
  if(xdex==0)xdex++;
  float32fixed<2> m_x=m_par->m_x_min+m_par->h_mx*xdex, m_y=m_par->m_y_min+m_par->h_my*ydex;
  float32fixed<4> theta=atan(sqrt(m_x*m_x+m_y*m_y));
//   cout<<"in slope componets roi theta, theta must be in ["<<m_par->minimum_large_theta<<","<<m_par->maximum_large_theta<<"]"<<endl;
  if(theta<m_par->minimum_large_theta || theta>m_par->maximum_large_theta){
//     cout << "Our theta of "<<theta<<" is not in ["<<m_par->minimum_large_theta<<","<<m_par->maximum_large_theta<<"]"<<endl;
    theta=0;
  }
  return theta;
}

float32fixed<4> MMT_Fitter::Slope_Components_ROI_phi(int jy, int ix) const{
  if(jy<0||jy>=m_par->n_y){
    cerr << "You picked a y slope road index of " << jy << " in Slope_Components_ROI_phi; there are only " << m_par->n_y << " of these.\n";
    if(jy>=m_par->n_y)jy=m_par->n_y-1;
    else jy=0;
//    exit(2);
  }
  if(ix<0||ix>=m_par->n_x){
    cerr << "You picked an x slope road index of " << ix << " in Slope_Components_ROI_phi; there are only " << m_par->n_x << " of these.\n";
    //right now we're assuming these are cases just on the edges and so put the values to the okay limits
    if(ix>=m_par->n_x)ix=m_par->n_x-1;
    else ix=0;
//     exit(2);
  }
  int xdex=ix,ydex=jy+1;
  float32fixed<2> m_x=m_par->m_x_min+m_par->h_mx*xdex, m_y=m_par->m_y_min+m_par->h_my*ydex;
  if(debug)cout<<"m_par->m_x_min+m_par->h_mx*xdex="<<m_par->m_x_min.getValue()<<"+"<<m_par->h_mx.getValue()<<"*"<<xdex<<"="<<m_x.getValue()<<", ";
  if(debug)cout<<"m_par->m_y_min+m_par->h_my*ydex="<<m_par->m_y_min.getValue()<<"+"<<m_par->h_my.getValue()<<"*"<<ydex<<"="<<m_y.getValue()<<", ";
  float32fixed<4> phi(atan2(m_x.getValue(),m_y.getValue()));//the definition is flipped from what you'd normally think
  if(debug)cout<<"for a phi of "<<phi.getValue()<<endl;
  if(phi<m_par->minimum_large_phi || phi>m_par->maximum_large_phi){
    if(debug) cout<<"Chucking phi of " << phi.getValue()<<" which registers as not in ["<<m_par->minimum_large_phi.getValue()<<","<<m_par->maximum_large_phi.getValue()<<"]"<<endl;
    phi=999;
  }
  return phi;
}

int MMT_Fitter::Rough_ROI_temp(float32fixed<4> theta, float32fixed<4> phi) const{
  //temporary function to identify areas of the wedge.
  float32fixed<4> minimum_large_theta=m_par->minimum_large_theta, maximum_large_theta=m_par->maximum_large_theta,minimum_large_phi=m_par->minimum_large_phi, maximum_large_phi=m_par->maximum_large_phi;
  int n_theta_rois=32, n_phi_rois=16;//*** ASK BLC WHAT THESE VALUES OUGHT TO BE!

  float32fixed<4> h_theta = (maximum_large_theta - minimum_large_theta)/n_theta_rois;
  float32fixed<4> h_phi = (maximum_large_phi - minimum_large_phi)/n_phi_rois;
  //how is this done in the FPGA? No such division, for sure
  double roi_t = ceil((theta - minimum_large_theta)/h_theta.getFloat());
  double roi_p = ceil((phi - minimum_large_phi)/h_phi.getFloat());

  if(theta<minimum_large_theta || theta>maximum_large_theta) roi_t = 0;
  if(phi<minimum_large_phi || phi>maximum_large_phi) roi_p = 0;
  int ret_val=roi_t * 1000 + roi_p;
  /*
  cout<<"ret_val:"<<ret_val<<endl;
  if(ret_val<1){
    cout << "ROI theta of "<<theta<<" is not in ["<<m_par->minimum_large_theta<<","<<m_par->maximum_large_theta<<"]"<<endl;
    cout << "ROI phi of "<<phi<<" is not in ["<<m_par->minimum_large_phi<<","<<m_par->maximum_large_phi<<"]"<<endl;
  }
  */
  return ret_val;
}

