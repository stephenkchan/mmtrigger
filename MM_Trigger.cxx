#include "MM_Trigger.h"
const bool sdebug=false;
const bool print=false;

MM_Trigger::MM_Trigger(string fname, par_par var, const string& tname)
  :MMT_Loader(fname,var,tname),m_find(MMT_Finder(m_par)),redo(var.redo){//m_find(MMT_Finder((new MMT_Parameters(var)))){
  if(sdebug) cout << "MM_T::building trigger" << endl;
  plotcounter=0;
  n_etabins=m_par->n_etabins;
  n_phibins=m_par->n_phibins;
  correct_bcid=2;bad_bcid=0;
  m_etabins=m_par->m_etabins;
  m_phibins=m_par->m_phibins;

  float maxtheta=m_par->maximum_large_theta.getFloat(),mintheta=m_par->minimum_large_theta.getFloat(),maxphi=m_par->maximum_large_phi.getFloat(),minphi=m_par->minimum_large_phi.getFloat();
  etalo=-log(tan(0.5*maxtheta));etahi=-log(tan(0.5*mintheta));
  int nbins=300;
  TH1::SetDefaultSumw2();

  //begin0 puts you at the place where the half setup starts
  //begin1 puts you at the place where the energy starts (setup is included in par_par)
  string delim="_",hark="Det";
  if(fname.find(hark+delim)==string::npos){
    delim="";hark="ValAlg";
  }
  int begin0=fname.find(hark+delim)+4,begin1=fname.find(delim,begin0+1)+1,sublen=fname.length()-begin1-5;//-begin0-5; the 5 is ".root"
  nom=fname.substr(begin1,sublen);
  double reslo=-0.005*pi(),reshi=0.005*pi(),preslo=-0.05*pi(),preshi=0.05*pi(),dthlo=-0.03*pi(),dthhi=0.03*pi();
//   m_dtheta=new TH1D("dtheta","dtheta",nbins,dthlo,dthhi);
  m_theta_tru_entpos=new TH1D("theta_tru_entpos","theta_tru_entpos",nbins,reslo,reshi);
  m_fit_the=new TH1D("fit_theta","fit_theta",nbins,mintheta,maxtheta);
  m_fit_phi=new TH1D("fit_phi","fit_phi",nbins,minphi,maxphi);
  m_fit_dth=new TH1D("fit_dtheta","fit_dtheta",nbins,reslo,reshi);
  m_hcode=new TH2D("hcode","hcode",44,-0.5,44.5,99,-49.5,49.5);//x is reco plane CT hit code; y is CT_reco hit code

  int bcidbins=m_par->setup.size()+1;//any of nplanes and the 0 planes case...
  int nsimbins=m_par->nsimmax_1d();
  m_crep_the=vector<vector<vector<TH1D*> > >(n_etabins,vector<vector<TH1D*> >(n_phibins,vector<TH1D*>(nsimbins,NULL)));
  m_crep_phi=vector<vector<vector<TH1D*> > >(n_etabins,vector<vector<TH1D*> >(n_phibins,vector<TH1D*>(nsimbins,NULL)));
  m_crep_dth=vector<vector<vector<TH1D*> > >(n_etabins,vector<vector<TH1D*> >(n_phibins,vector<TH1D*>(nsimbins,NULL)));
  
  for(int i=0;i<n_etabins;i++){
    TString estr(m_par->eta_str(i));
    m_all_charge.push_back(new TH1D("all_charge"+estr,"all_charge"+estr,nbins,0.,5.));
    m_fit_charge.push_back(new TH1D("fit_charge"+estr,"fit_charge"+estr,nbins,0.,5.));
    m_res_the.push_back(new TH1D("res_theta"+estr,"res_theta"+estr,nbins,reslo,reshi));
    m_res_phi.push_back(new TH1D("res_phi"+estr,"res_phi"+estr,nbins,preslo,preshi));
    m_res_dth.push_back(new TH1D("res_dtheta"+estr,"res_dtheta"+estr,nbins,dthlo,dthhi));
    m_res_mxl.push_back(new TH1D("res_mXl"+estr,"res_mXl"+estr,nbins,reslo,reshi));
     m_res_mx.push_back(new TH1D("res_mx"+estr,"res_mx"+estr,nbins,reslo,reshi));
     m_res_my.push_back(new TH1D("res_my"+estr,"res_my"+estr,nbins,0.2*reslo,0.2*reshi));
   m_res_mxuv.push_back(new TH1D("res_mxuv"+estr,"res_mxuv"+estr,nbins,0.7*reslo,0.7*reshi));
   m_bcid_eff.push_back(new TH2D("bcid_eff"+estr,"bcid_eff"+estr,bcidbins,-0.5,bcidbins*1.-0.5,bcidbins,-0.5,bcidbins*1.-0.5));
    for(int j=0;j<3;j++){
      TString oof;oof+=j;
      m_res_mu.push_back(new TH1D("res_mu"+oof+estr,"res_mu"+oof+estr,nbins,0.7*reslo,0.7*reshi));
      m_res_mv.push_back(new TH1D("res_mv"+oof+estr,"res_mv"+oof+estr,nbins,0.7*reslo,0.7*reshi));
    }
    for(int j=0;j<4;j++){
      TString oof;oof+=j;
      m_res_mxn.push_back(new TH1D("res_mx"+oof+estr,"res_mx"+oof+estr,nbins,0.7*reslo,0.7*reshi));
    }
    for(int j=0;j<n_phibins;j++){
      TString pstr(m_par->phi_str(j));
      for(int k=0;k<nsimbins;k++){
	TString xstr(m_par->index_to_hit_str(k));
	m_crep_the[i][j][k]=new TH1D("res_theta"+estr+pstr+xstr,"res_theta"+estr+pstr+xstr,nbins,reslo,reshi);
	m_crep_phi[i][j][k]=new TH1D("res_phi"+estr+pstr+xstr,"res_phi"+estr+pstr+xstr,nbins,preslo,preshi);
	m_crep_dth[i][j][k]=new TH1D("res_dtheta"+estr+pstr+xstr,"res_dtheta"+estr+pstr+xstr,nbins,dthlo,dthhi);
      }
    }
  }
  m_store_no=new TH1D("store_no","store_no",10,0.5,10.5);
  nev=get_nevents();nana=0;
  nev_ent_wedge=new TH1D("nev_ent_wedge","nev_ent_wedge",n_etabins,0.5,n_etabins+0.5);
  nev_ct=new TH1D("nev_ct","nev_ct",n_etabins,0.5,n_etabins+0.5);
  nev_fit_clean=new TH1D("nev_fit_clean","nev_fit_clean",n_etabins,0.5,n_etabins+0.5);
  nev_fit_bg=new TH1D("nev_fit_bg","nev_fit_bg",n_etabins,0.5,n_etabins+0.5);
  nev_bcid_pass=new TH1D("nev_bcid_pass","nev_bcid_pass",n_etabins,0.5,n_etabins+0.5);
//   m_find=MMT_Finder(m_par);

  int nxp=m_par->q_planes("x").size(),nuvp=m_par->q_planes("u").size()+m_par->q_planes("v").size();
  nev_tru_xuv=new TH3D("nev_tru_xuv","nev_tru_xuv",nxp+1,-0.5,nxp+0.5,nuvp+1,-0.5,nuvp+0.5,n_etabins,0.5,n_etabins+0.5);
  nev_fit_xuv=new TH3D("nev_fit_xuv","nev_fit_xuv",nxp+1,-0.5,nxp+0.5,nuvp+1,-0.5,nuvp+0.5,n_etabins,0.5,n_etabins+0.5);
  nev_rec_xuv=new TH3D("nev_rec_xuv","nev_rec_xuv",nxp+1,-0.5,nxp+0.5,nuvp+1,-0.5,nuvp+0.5,n_etabins,0.5,n_etabins+0.5);
  nev_rft_xuv=new TH3D("nev_rft_xuv","nev_rft_xuv",nxp+1,-0.5,nxp+0.5,nuvp+1,-0.5,nuvp+0.5,n_etabins,0.5,n_etabins+0.5);
  nev_rpl_xuv=new TH3D("nev_rpl_xuv","nev_rpl_xuv",nxp+1,-0.5,nxp+0.5,nuvp+1,-0.5,nuvp+0.5,n_etabins,0.5,n_etabins+0.5);
}

void MM_Trigger::analysis(int nevents){
  int evstart=0;//1e4+1;
  if(nevents+evstart<0||nevents+evstart>nev) nevents=nev;
  string dir=(file_name.find_last_of("/")==string::npos||file_name.find_last_of("/")==-1?".":file_name.substr(0,file_name.find_last_of("/"))),hist_printer_name=dir+"/printouts/mmt_hit_print"+m_par->param_par().print_pars()+nom+".txt";
  ofstream hit_printer(hist_printer_name.c_str());
  print_header(hit_printer);
  for(int iev=evstart; iev<nevents+evstart; iev++){
    if(iev%100==0) progress_bar(iev-evstart,nevents);
    if(process_event(iev,hit_printer))nana++;//i.e. if the event was in the correct wedge type
    if(debug) cout<<"After "<<iev<<endl;
  }
  cout<<endl;
}

void MM_Trigger::print_header(ofstream& print)const{
  print<<"File: "<<file_name<<endl;
  print<<"Alg config: "<<m_par->param_par().print_pars()<<endl;
  print<<"-----------Y BASES-----------"<<endl;
  for(int sta=0;sta<m_par->n_stations_eta;sta++){
    for(int pl=0;pl<(int)m_par->setup.size();pl++)print<<setw(10)<<m_par->ybases[pl][sta].getFloat();
    print<<endl;
  }
  print<<"-----------Z BASES-----------"<<endl;
  for(int bin=0;bin<m_par->ybins;bin++){
    for(int pl=0;pl<(int)m_par->setup.size();pl++)print<<setw(10)<<m_par->z_large[bin][pl].getFloat();
    print<<endl;
  }
  print<<"-----------------------------"<<endl;
}
void MM_Trigger::print_ev_header(ofstream& print,int ev)const{
  ostringstream first;first<<"%Ev"<<ev<<" t";
  print<<setw(12)<<first.str()<<setw(12)<<"VMM"<<setw(12)<<"plane"<<setw(12)<<"station"<<setw(12)<<"strip"<<setw(12)<<"slope"<<endl;
}
void MM_Trigger::print_hit(ofstream& print,const hdst_entry&hit)const{
  hdst_key key=hit.entry_key();hdst_info info=hit.entry_info(m_par);
  print<<setw(12)<<key.time<<setw(12)<<key.VMM_chip<<setw(12)<<info.plane<<setw(12)<<hit.station_eta<<setw(12)<<hit.strip<<setw(12)<<info.slope.getFloat()<<endl;
}
void MM_Trigger::print_event_summary(ofstream& print,int ev,double tru_th,double tru_ph,double tru_dth,double mx,double my,double mxl,double theta,double phi,double dtheta)const{
  print<<"--------EVENT "<<ev<<" SUMMARY: (tru_th,tru_ph,tru_dth,mx,my,mxl,theta,phi,dtheta)=("<<setw(12)<<tru_th<<","<<setw(12)<<tru_ph<<","<<setw(12)<<tru_dth<<","<<setw(12)<<mx<<","<<setw(12)<<my<<","<<setw(12)<<mxl<<","<<setw(12)<<theta<<","<<setw(12)<<phi<<","<<setw(12)<<dtheta<<")--------"<<endl;
}

bool MM_Trigger::process_event(int event,ofstream&hit_printer){
  m_par->fill0=false;
  clear_buffers(event);
  if(debug) cout << "START MMT::process_event #"<<event<<endl;
//   else nev_ct++;
  bool pass_cuts=load_event(event);
  evInf_entry boom(Event_Info.find(event)->second);
  if(print) boom.print();
  double tent=boom.theta_ent,tpos=boom.theta_pos,ppos=boom.phi_pos,pent=boom.phi_ent,dt=boom.dtheta;
  int ebin=m_par->eta_bin(tpos);
  bool bigrth=false;
  print_ev_header(hit_printer,event);
  vector<hdst_entry> hdsts(event_hdsts(event));
  for(int i=0;i<(int)hdsts.size();i++)print_hit(hit_printer,hdsts[i]);
  if(hdsts.empty()) {print_event_summary(hit_printer,event,tpos,ppos,dt);return false;}
  double mu=0,mv=0,mxr=0;int nu=0,nv=0;
  vector<int> x_pl=m_par->q_planes("x"),u_pl=m_par->q_planes("u"),v_pl=m_par->q_planes("v");vector<bool>pl_reco_hit(m_par->setup.size(),false);
  for(unsigned int ihit=0;ihit<hdsts.size();ihit++){
    pl_reco_hit[hdsts[ihit].plane]=true;
    m_all_charge[ebin]->Fill(hdsts[ihit].charge);
    if(hdsts[ihit].entry_key()==hdst_key()){
      cout<<"HIT!? "<<ihit<<endl;
      continue;
    }
    TVector3 truth=hdsts[ihit].truth;int pl=hdsts[ihit].plane;
    if(m_par->setup.substr(pl,1)=="u"){
      double slope=hdsts[ihit].entry_info(m_par).slope.getValue(),mtru=(1+tan(degtorad(m_par->stereo_degree.getFloat()))*tan(ppos))*cos(ppos)*tan(tpos);
      mu+=mtru;nu++;
      if(debug)cout<<" tru="<<mtru<<", ureco="<<slope<<endl;
      if(pl==u_pl.front()) m_res_mu[3*ebin]->Fill(slope-mtru);
      if(pl==u_pl.back())m_res_mu[3*ebin+1]->Fill(slope-mtru);
    }
    if(m_par->setup.substr(pl,1)=="v"){
      double slope=hdsts[ihit].entry_info(m_par).slope.getValue(),mtru=(1-tan(degtorad(m_par->stereo_degree.getFloat()))*tan(ppos))*cos(ppos)*tan(tpos);
      mv+=mtru;nv++;
      if(debug)cout<<" tru="<<mtru<<", vreco="<<slope<<endl;
      if(pl==v_pl.front()) m_res_mv[3*ebin]->Fill(slope-mtru);
      if(pl==v_pl.back())m_res_mv[3*ebin+1]->Fill(slope-mtru);
    }
    else{
      double slope=hdsts[ihit].entry_info(m_par).slope.getValue();int bin=4*ebin,plane=hdsts[ihit].plane;
      if(plane==1)bin++;
      else if(plane==6)bin+=2;
      else if(plane==7)bin+=3;
      m_res_mxn[bin]->Fill(slope-tan(tpos)*cos(ppos));
    }
  }
  int xrec=0,uvrec=0;
  for(int ipl=0;ipl<(int)x_pl.size();ipl++) xrec+=pl_reco_hit[x_pl[ipl]];
  for(int ipl=0;ipl<(int)u_pl.size();ipl++)uvrec+=pl_reco_hit[u_pl[ipl]];
  for(int ipl=0;ipl<(int)v_pl.size();ipl++)uvrec+=pl_reco_hit[v_pl[ipl]];
  nev_rpl_xuv->Fill(xrec,uvrec,ebin+1);
  if(nu==2&&nv==2&&hdsts.size()==8){
    double mutru=0.5*mu,mvtru=0.5*mv,mu_mv_tan=0.5*(mutru-mvtru)/tan(degtorad(m_par->stereo_degree.getFloat()));
    m_res_mxuv[ebin]->Fill(mu_mv_tan-tan(tpos)*sin(ppos));
  }
  int x_tru=Event_Info.at(event).N_X_hits,uv_tru=Event_Info.at(event).N_UV_hits;
  nev_tru_xuv->Fill(x_tru,uv_tru,ebin+1);
  if(!pass_cuts){
    if(debug) cout<<"Event "<<event<<" does not pass cuts\n";
    return !(boom.bad_wedge);
  }
  if(m_par->minimum_large_phi<pent&&m_par->maximum_large_phi>pent && m_par->minimum_large_theta<tent&&m_par->maximum_large_theta>tent){
    nev_ent_wedge->Fill(ebin+1);
  }
  if(debug)cout<<"Event "<<event<<"start!...tpos="<<tpos<<"...ppos="<<ppos<<"...tent="<<tent<<"...pent="<<pent<<endl;
//   if(xhit<m_par->CT_x||uvhit<m_par->CT_uv)cout<<"CT ("<<m_par->CT_x<<","<<m_par->CT_uv<<") filled for ebin+1 "<<ebin+1<<" m_tru_xuv was filled at "<<xhit<<","<<uvhit<<","<<ebin+1<<endl;
  m_theta_tru_entpos->Fill(tent-tpos);
  vector<evFit_entry> road_fits;
  vector<pair<double,double> >mfits;double mxl;
  map<pair<int,int>,finder_entry> road_tracks=m_find.finder_event(event,hdsts,road_fits,mfits,mxl);
  int max_hcode=0,min_hcode=0;
  for(int road=0;road<(int)road_fits.size();road++){
    if(road_fits[road].hcode>max_hcode)max_hcode=road_fits[road].hcode;//means it passed reco CT; cf MMT_Finder::Coincidence_Gate
    if(road_fits[road].hcode<min_hcode)min_hcode=road_fits[road].hcode;//means it passed reco CT; cf MMT_Finder::Coincidence_Gate
  }
  if(max_hcode>0)nev_ct->Fill(ebin+1);
  int hcodeyfill=(max_hcode>abs(min_hcode)?max_hcode:min_hcode);
  m_hcode->Fill(10*xrec+uvrec,hcodeyfill);
  if(xrec==4&&uvrec==4&&hcodeyfill<-20&&debug)print_rd_summary(road_tracks,event);
  nev_rec_xuv->Fill(max_hcode/10,max_hcode%10,ebin+1);
  for(int ihit=0;ihit<(int)hdsts.size();ihit++){
    if(hdsts[ihit].fit_theta.getValue()>0)m_fit_charge[ebin]->Fill(hdsts[ihit].charge);
  }
  m_res_mxl[ebin]->Fill(mxl-tan(tent)*cos(pent));
  for(int i=0;i<(int)mfits.size();i++){
    if(mfits[i].first-tan(tpos)*sin(ppos)<-0.004&&debug)cout<<"FLAG! (truth mx="<<tan(tpos)*sin(ppos)<<")"<<endl<<endl;
    double mutru=(1+tan(degtorad(m_par->stereo_degree.getFloat()))*tan(ppos))*cos(ppos)*tan(tpos),mvtru=(1-1*tan(degtorad(m_par->stereo_degree.getFloat()))*tan(ppos))*cos(ppos)*tan(tpos),mxtru_alt=0.5*(mutru-mvtru)/tan(degtorad(m_par->stereo_degree.getFloat()));
    m_res_mx[ebin]->Fill(mfits[i].first-mxtru_alt);//tan(tpos)*sin(ppos));//sqrt(nume/(1+1./denflip))*(mfits[i].first>0?1:-1));
    m_res_my[ebin]->Fill(mfits[i].second-tan(tpos)*cos(ppos));//sqrt(nume/(1+denflip)));
  }
  double mx=-999,my=-999,thet=-999,ph=-999,dthet=-999;
  if(!mfits.empty()){mx=mfits.front().first;my=mfits.front().second;}
  bool did_clean_fit=false,did_bg_fit=false;
  int xhit=0,uvhit=0;
  for(unsigned int i=0; i<road_fits.size(); i++){
    if(road_fits[i].fit_roi<=0)continue;
    if(sdebug||print) road_fits[i].print();
    double the(road_fits[i].fit_theta.getFloat()),phi(road_fits[i].fit_phi.getFloat()),dth(road_fits[i].fit_dtheta.getFloat());
    if(the<-4||phi<-4||dth<-4)continue;
    thet=road_fits[i].fit_theta.getFloat();ph=road_fits[i].fit_phi.getFloat();dthet=road_fits[i].fit_dtheta.getFloat();
    if(debug)cout<<"fit "<<i<<": theta="<<road_fits[i].fit_theta.getFloat()<<", phi="<<road_fits[i].fit_phi.getFloat()<<", dtheta="<<road_fits[i].fit_dtheta.getFloat()<<endl;
    int cebin=m_par->eta_bin(the),cpbin=m_par->phi_bin(phi),ht=road_fits[i].hcode;
    m_fit_the->Fill(the); m_fit_phi->Fill(phi); m_fit_dth->Fill(dth); 
    m_res_the[ebin]->Fill(the-tpos-(m_par->correct.type==3?m_par->crep_table[cebin][cpbin][ht][0]:0));
    m_res_dth[ebin]->Fill(dth-dt-(m_par->correct.type==3?m_par->crep_table[cebin][cpbin][ht][2]:0));
    if(debug) cout<<endl<<"m_res_the["<<ebin<<"] fill "<<the-tpos<<endl;
    m_res_phi[ebin]->Fill(phi-ppos-(m_par->correct.type==3?m_par->crep_table[cebin][cpbin][ht][1]:0));
    if(debug) cout<<endl<<"m_res_phi["<<ebin<<"] fill "<<phi-ppos<<endl;
    if(debug)cout<<"for (tpos,ppos,dt)=("<<tpos<<","<<ppos<<","<<dt<<")...BINS[E"<<cebin<<"][P"<<cpbin<<"][x"<<ht<<"]:"<<" the="<<the-tpos<<",phi="<<phi-ppos<<",dth="<<dth-dt<<endl;
    m_crep_the[cebin][cpbin][ht]->Fill(the-tpos);m_crep_phi[cebin][cpbin][ht]->Fill(phi-ppos);m_crep_dth[cebin][cpbin][ht]->Fill(dth-dt);
    vector<hdst_key> keys(road_fits[i].fit_hit_keys);double uquant=0,vquant=0;int uct=0,vct=0;int gut=0,schlecht=0,womp=0;
    map<int,int>bcid_map;
    for(int ihit=0;ihit<(int)keys.size();ihit++){
      int pos=m_find.find_hdst(hdsts,keys[ihit]);
      int plane=hdsts[pos].plane,bcid=hdsts[pos].BC_time-10*(hdsts[pos].event+1);
      bcid_map[bcid]++;
      if(bcid==correct_bcid)gut++;
      else{schlecht++;womp=bcid;}
      if(m_par->setup.substr(plane,1)=="u"){uquant+=hdsts[ihit].entry_info(m_par).slope.getValue();uct++;}
      if(m_par->setup.substr(plane,1)=="v"){vquant+=hdsts[ihit].entry_info(m_par).slope.getValue();vct++;}
    }
    m_bcid_eff[ebin]->Fill(gut,schlecht+gut);
    int maj_bunch=-999,maj_ct=0;
    for(map<int,int>::const_iterator it=bcid_map.begin();it!=bcid_map.end();++it){
      if(it->second>maj_ct){maj_bunch=it->first;maj_ct=it->second;}
    }
    if(maj_bunch==correct_bcid)nev_bcid_pass->Fill(ebin+1);
    else{
      cerr<<"Particularly bad track; majority BCID: "<<maj_bunch<<endl;
      bad_bcid++;
    }
    if(uct==2&&vct==2){
      m_res_mu[3*ebin+2]->Fill(0.5*uquant-(1+tan(degtorad(m_par->stereo_degree.getValue()))*tan(ppos))*tan(tpos)*cos(ppos));
      m_res_mv[3*ebin+2]->Fill(0.5*vquant-(1-tan(degtorad(m_par->stereo_degree.getValue()))*tan(ppos))*tan(tpos)*cos(ppos));
    }
    if(debug) cout<<endl<<"m_res_dth["<<ebin<<"] fill "<<dth<<" - "<<dt<<"="<<dth-dt<<endl;
    if(debug&&m_par->fill0)cout<<"FLAG! fit_dth="<<dth<<", tru="<<dt<<endl;
    if(road_fits[i].bg_X_fit==0&&road_fits[i].bg_UV_fit==0) did_clean_fit=true;
    else did_bg_fit=true;
    //       if(did_bg_fit) road_fits[i].print();
    xhit=road_fits[i].bg_X_fit+road_fits[i].X_hits_in_fit;
    uvhit=road_fits[i].bg_UV_fit+road_fits[i].UV_hits_in_fit;
  }
  print_event_summary(hit_printer,event,tpos,ppos,dt,mx,my,mxl,thet,ph,dthet);
  if(print){
    for(unsigned int i=0; i<hdsts.size(); i++) hdsts[i].print();
    cout<<endl;
  }
  if(did_clean_fit||did_bg_fit) {nev_fit_xuv->Fill(x_tru,uv_tru,ebin+1);nev_rft_xuv->Fill(xhit,uvhit,ebin+1);}
  if(did_bg_fit) nev_fit_bg->Fill(ebin+1);//nev_fit_bg++;
  else if(did_clean_fit) nev_fit_clean->Fill(ebin+1);//nev_fit_clean++;
  if(debug)cout<<"MMT::p_e::END EVENT #"<<event<<endl;
  return true;
}

void MM_Trigger::print_rd_summary(const map<pair<int,int>,finder_entry>&road_tracks,int ev){
  if(plotcounter>10)return;
  set_style();
  TString hname="hit_matrix"+nom+"_ev";hname+=ev;
  TH2D*hit_matrix=new TH2D(hname,hname,m_par->setup.size(),-0.5,m_par->setup.size()-0.5,m_find.get_roads(),-0.5,m_find.get_roads()-0.5);
  bool shouldbefine=false;
  for(int road=0;road<m_find.get_roads();road++){
    int nrdhit=0;
    for(int plane=0;plane<(int)m_par->setup.size();plane++){
      if(road_tracks.find(pair<int,int>(road,plane))==road_tracks.end())continue;
      hit_matrix->Fill(plane,road);nrdhit++;
    }
    if(nrdhit==8)shouldbefine=true;
  }
  TString pnomine=(shouldbefine?"FINE":"s");pnomine+=hname;
  hname+="c";TCanvas*c=new TCanvas(hname,hname,1000,1000);TLatex l=plot_latex();TString title="Event ";title+=ev;title+=" Slope Roads";
  hit_matrix->Draw("colz");
  l.DrawLatex(0.20,0.95,title);
  TString pdfname="/n/atlasfs/atlasdata/atlasdata1/stchan/mmtrigger/runtime_diag/"+pnomine+".pdf";
  c->Print(pdfname);
  delete hit_matrix;delete c;
  plotcounter++;
}

string MM_Trigger::filename(const string& filedir,const string& opttag)const{
  ostringstream oname; oname<<filedir<<(filedir.substr(filedir.length()-1)=="/"?"":"/")<<"mmt_raw"<<m_par->param_par().print_pars()<<"_"<<nom<<(opttag==""?"":opttag.substr(0,1).compare("_")==0?"":"_")<<opttag<<".root";
  return oname.str();
}

bool MM_Trigger::file_exists(const string&name)const{
  if(debug)cout<<"checking to see if "<<name<<" exists."<<endl;
  string path="./",file=name;
  if(file.find("/")!=string::npos){
    path=name.substr(0,name.find_last_of("/")+1);file=name.substr(name.find_last_of("/")+1);
  }
  DIR*red=opendir(path.c_str());struct dirent*andy;
  if(red==NULL)return false;
  while((andy=readdir(red))!=NULL){
    if(debug)cout<<"compare "<<andy->d_name<<" "<<file<<endl;
    if(string(andy->d_name)==file)return true;
  }
  return false;
}

void MM_Trigger::save_analysis(const string& filedir,int nevents,const string& opttag,bool pcrep){
  if(!redo && file_exists(filename(filedir,opttag))){
    cout<<"File "<<filename(filedir,opttag)<<" already exists."<<endl;
    return;
  }
  analysis(nevents);
  save_file(filedir,opttag,pcrep);
}

void MM_Trigger::save_file(const string& filedir,const string& opttag,bool pcrep){
  bool diag=m_par->diag&&debug;string oname=filename(filedir,opttag);
  if(m_par->correct.type==0)pcrep=true; 
  TFile *f=TFile::Open(oname.c_str(),"RECREATE");
  if(pcrep){
    ostringstream crep_nom;crep_nom<<filedir<<(filedir.substr(filedir.length()-1)=="/"?"":"/")<<"pcrep"<<m_par->param_par().print_pars()<<"_"<<nom<<(opttag==""?"":opttag.substr(0,1).compare("_")==0?"":"_")<<opttag<<".txt";
    ofstream crep(crep_nom.str().c_str());
    for(int i=0;i<n_etabins;i++){
      string estr=m_par->eta_str(i);
      for(int j=0;j<n_phibins;j++){
	string pstr=m_par->phi_str(i);
	for(int k=0;k<m_par->nsimmax_1d();k++)crep<<estr<<pstr<<m_par->index_to_hit_str(k)<<" "<<m_crep_the[i][j][k]->GetMean()<<" "<<m_crep_phi[i][j][k]->GetMean()<<" "<<m_crep_dth[i][j][k]->GetMean()<<endl;
      }
    }
    cout<<"Wrote sim correct results to "<<crep_nom.str()<<endl;
  }
  if(diag){
    int stations=m_par->n_stations_eta;
    for(int i=0;i<stations;i++){
      m_diff_xuv[i]->Write();
      if(debug)cout<<"STATION_ETA "<<i<<":";
      ostringstream name;name<<"global_strippos_zeros_es"<<i;TString a(name.str());
      TH1D *dummy=new TH1D(a,a,(int)m_par->setup.size(),0.5,(int)m_par->setup.size()+0.5);
      for(int j=0;j<(int)m_par->setup.size();j++){
	TString nom(gposx_strippos_l[i][j]->GetName());nom+="line";
	TF1 *line=new TF1(nom,"[0]*x+[1]",0,5000);
	TFitResultPtr yoink=gposx_strippos_l[i][j]->Fit(nom,"RQS");
	dummy->SetBinContent(j+1,line->GetParameter(1));
	if(debug)cout<<"\t("<<j<<") "<< line->GetParameter(1);
      }
      if(debug)cout<<endl;
      dummy->Write();
    }
    for(int i=0;i<stations;i++){
      if(debug)cout<<"STATION_ETA "<<i<<":";
      ostringstream name;name<<"gloloc_strippos_zeros_es"<<i;TString a(name.str());
      TH1D *dummy=new TH1D(a,a,(int)m_par->setup.size(),0.5,(int)m_par->setup.size()+0.5);
      for(int j=0;j<(int)m_par->setup.size();j++){
	TString nom(glposx_strippos_l[i][j]->GetName());nom+="line";
	TF1 *line=new TF1(nom,"[0]*x+[1]",0,5000);
	TFitResultPtr yoink=glposx_strippos_l[i][j]->Fit(nom,"RQS");
	dummy->SetBinContent(j+1,line->GetParameter(1));
	if(debug)cout<<"\t("<<j<<") "<< line->GetParameter(1);
      }
      if(debug)cout<<endl;
      dummy->Write();
    }
    for(int i=0;i<stations;i++){
      for(int j=0;j<(int)m_par->setup.size();j++) ybases_solved[i][j]->Write();
    }
    for(int i=0;i<stations;i++){
      for(int j=0;j<(int)m_par->setup.size();j++) ybases_pos[i][j]->Write();
    }
    for(int i=0;i<stations;i++){
      for(int j=0;j<(int)m_par->setup.size();j++) ybases_ent[i][j]->Write();
    }
    for(int i=0;i<stations;i++){
      for(int j=0;j<(int)m_par->setup.size();j++) strip_widths[i][j]->Write();
    }
  }
  m_hcode->Write();
  nev_ent_wedge->Write();nev_ct->Write();nev_fit_clean->Write();nev_fit_bg->Write();
  nev_tru_xuv->Write();nev_fit_xuv->Write();nev_rec_xuv->Write();nev_rft_xuv->Write();nev_rpl_xuv->Write();nev_bcid_pass->Write();
  m_fit_the->Write(); m_fit_phi->Write(); m_fit_dth->Write();
  for(int i=0;i<n_etabins;i++){
    m_res_the[i]->Write(); m_res_phi[i]->Write(); m_res_dth[i]->Write();m_res_mxl[i]->Write();m_res_mx[i]->Write();m_res_my[i]->Write();m_res_mxuv[i]->Write();
    m_bcid_eff[i]->Write();m_all_charge[i]->Write();m_fit_charge[i]->Write();
  }
  for(int i=0;i<(int)m_res_mu.size();i++)m_res_mu[i]->Write();
  for(int i=0;i<(int)m_res_mv.size();i++)m_res_mv[i]->Write();
  for(int i=0;i<(int)m_res_mxn.size();i++)m_res_mxn[i]->Write();
  m_theta_tru_entpos->Write();// m_dtheta->Write();
  m_store_no->SetBinContent(1,nana);
  m_store_no->SetBinContent(2,nev_ct->Integral());
  m_store_no->SetBinContent(3,nev_fit_clean->Integral());
  m_store_no->SetBinContent(4,nev_fit_bg->Integral());
  m_store_no->SetBinContent(5,gen_bg);
  m_store_no->SetBinContent(10,nev);
  m_store_no->Write();
  f->Close();
  cout<<nev_fit_clean->Integral()<<" fits ("<<bad_bcid<<" bad BCID) of "<<nev_ct->Integral()<<" passing events---saved output to "<<oname<<endl;
}

void MM_Trigger::progress_bar(int current, int total){
  double progress=100.*current/total;
  int ticks=progress/2,stars=50-ticks;
  cerr << "\r Progress: |";
  for(int i=0; i<ticks; i++) cerr <<"-";
  cerr << ">";
  for(int i=0; i<stars-1; i++) cerr <<"*";
  cerr << "| " << setprecision(4) << current << " of " << total<< " (" << progress << "%)      " << std::flush;
}

/*
int MM_Trigger::find_earliest_hdst(const map<int,hdst_entry>& hdst) const{
  int minbctime=hdst.begin()->second.BC_time, mintime=hdst.begin()->second.time, min_dex=0;
  for(map<int,hdst_entry>::iterator it=hdst.begin(); it!=hdst.end(); ++it){
    if(it->second.BC_time<minbctime){
      minbctime=it->second.BC_time, mintime=it->second.time, min_dex=it-map.begin();
    }
    else if(it->second.BC_time==minbctime && it->second.time<mintime){
      minbctime=it->second.BC_time, mintime=it->second.time, min_dex=it-map.begin();
    }
  }
  return min_dex;
}

void MM_Trigger::Run_File(const string& background, int nentries){
  //some variables
  bool Generate_FPGA_Sample=true;
  double H_values, UV_values, CT_values, CT_x_values, CT_uv_values,LG_max=0, LG_min=1000;
  int N_coincidence_threshold_met, N_Fit;

  double energy = 200;  //in GeV
  double stereodeg = 15; // 1.5deg = 15
  double charge = 1000, z_shift = 0;
  int entries=(nentries<0||nentries<tr->GetEntries()?tr->GetEntries:nentries);
  for(int ien=0; ien<entries; ien++){
    proess_event(ien);
  }
  if(background=="shift"||background=="bgoff"){
    Event_Info = Event_Info_Without_Background;
    Hits_Data_Set_Time = Hits_Data_Set_Time_Without_Background;
  }
  else if(background=="incoherent"){
    energy = 0;
    Event_Info = Event_Info_Background;
    Hits_Data_Set_Time = Master_Background_Hits;
    Hits_Data_Set_Time(:,10)=Hits_Data_Set_Time(:,1).*2+(Hits_Data_Set_Time(:,10)-Hits_Data_Set_Time(:,1)*10);
    Hits_Data_Set_Time(:,14)=Hits_Data_Set_Time(:,14)-100*Hits_Data_Set_Time(:,1)+Hits_Data_Set_Time(:,1)*2;
    Mimic_VMM_Chip_Deadtime(Hits_Data_Set_Time);
  }
  else if(background=="bgon"){
    Event_Info = Event_Info_With_Background;
    Hits_Data_Set_Time = Hits_Data_Set_Time_With_Background;
    //         load('Background_30000events_2BCunif_dr10', 'Master_Background_Hits')
    //         Event_Info = Event_Info_Without_Background;
    //         Hits_Data_Set_Time = Mimic_VMM_Chip_Deadtime(sortrows(...
    //             [Hits_Data_Set_Time_Without_Background;Master_Background_Hits],[10,14]));
    //         Hits_Data_Set_Time = Mimic_VMM_Chip_Deadtime(sortrows(Master_Background_Hits,[10,14]));
  }

  [n_events_possible,aa] = size(Event_Info);
  Event_Fit = [Event_Info(:,1),zeros(n_events_possible,15)];
  double h = 0.0009;  //2.5*10^(-4);  //0.0009;  //8.5*10^(-4);  //   // 4* 10^(-4);   // //4*10^(-4);
  double uv_error = 0.0035;  //  0.004;  // //0.0035;   //3.5*10^(-3);     //3.4*10^(-3);

  for(unsigned int iz=0; iz<m_par->z_large.size(); iz++) m_par->z_large[iz]+=z_shift;
  //tic;
  N=100000;
  int event_mark = Finder_Control(Hits_Data_Set_Time,BC_window,N);
  //Front_Filter_3(Hits_Data_Set_Time,BC_window); event_mark = 30000;  //,N);  //faster version by ~4x
  //toc;

  double H_values=h,UV_values=uv_error,CT_values=CT,CT_x_values=CT_x,CT_uv_values=CT_uv,Charge_values=charge,Z_shift_values = z_shift;
  //Event_Fit = Event_Fit(1:event_mark,:); Event_Info = Event_Info(1:event_mark,:);//
}

*/
