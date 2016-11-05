#include "MMT_Loader.h"
const bool ldebug=false&&debug,striphack=true;

MMT_Loader::MMT_Loader(string fname, par_par var, const string& tname)
  :m_par(new MMT_Parameters(var)),file_name(fname){
  n_bar=0;n_ent=0;
  if(ldebug) cout<<"MMT_L---parameters loaded\n";
  if(fname.find(".root")==string::npos) fname+=".root";
  if(ldebug) cout<<"MMT_L---want to load "<<fname<<endl;
  string dir=fname.substr(0,fname.find_last_of("/")),file=fname.substr(fname.find_last_of("/")+1);
  DIR *fdir=opendir(dir.c_str()); struct dirent *ent=NULL; bool found=false;
  if(fdir==NULL){
    cerr<<"Can't open: "<<dir<<endl; exit(0);
  }
  while((ent=readdir(fdir))){
    if(file.compare(ent->d_name)==0){
      found=true; break;
    }
  }
  if(!found){
    cerr<<"Couldn't find: "<<fname<<endl;
    exit(1);
  }
  f=TFile::Open(fname.c_str());
  if(f==0||!f->IsOpen()){
    cerr<<"Could not open file "<<fname<<endl;
    exit(0);
  }
  if(ldebug) cout<<"MMT_L---want to tree "<<(tname==""?"NSWHitsTree":tname.c_str())<<endl;
  f->GetObject((tname==""?"NSWHitsTree":tname.c_str()),tr);
  Init(tr);
  if(ldebug) cout<<"MMT_L---tree loaded"<<endl;
  VMM_deadtime = 100;  //(ns)
  num_VMM_per_plane = 1000;
  int estations=m_par->n_stations_eta;
  reset_VMMs();
  m_rand=new TRandom3();
  double r_min=m_par->H.getFloat(),r_max=(m_par->Hnom+m_par->h1)/cos(degtorad(m_par->wedge_opening_angle.getFloat())*0.5);//,dr = 10;  int r_zones = ceil((r_max - r_min)/dr);
  hit_rate=new TF1("get_r","x*pow(0.1*x,-2.125)",r_min,r_max);//changed...

  hit_rate_GHz=0.0986657;//we just did the integral in Mathematica
  tlo=-50.;thi=170.;
  double interval=thi-tlo;
  hit_ct_rand=new TF1("get_hits","TMath::PoissonI(x,[0])",0,4*interval*hit_rate_GHz);//out to four times expectation cause why not
  hit_ct_rand->SetParameter(0,interval*hit_rate_GHz);
  gen_bg=m_par->genbg;
  wedgey=m_par->H+m_par->h1*0.5; wedgex=0;
  uvxxmod=(m_par->setup.compare("xxuvuvxx")==0);
  gposx_strippos_l=vector<vector<TH2D*> >(estations,vector<TH2D*>());
  glposx_strippos_l=vector<vector<TH2D*> >(estations,vector<TH2D*>());
  ybases_solved=vector<vector<TH1D*> >(estations,vector<TH1D*>());
  ybases_pos=vector<vector<TH1D*> >(estations,vector<TH1D*>());
  ybases_ent=vector<vector<TH1D*> >(estations,vector<TH1D*>());
  strip_widths=vector<vector<TH1D*> >(estations,vector<TH1D*>());
  for(int i=1;i<=estations;i++){
    TString xuv("diff_xuv_es");xuv+=i;
    m_diff_xuv.push_back(new TH1D(xuv,xuv,4000,-3999.5,4000.5));
    for(int j=0;j<(int)m_par->setup.size();j++){
      ostringstream ech;ech<<"_es"<<i<<"_pl"<<j; string yar="ybases"+ech.str(),gar="swidths"+ech.str(),har="gposx_strippos_l"+ech.str(),asd="glposx_strippos_l"+ech.str(),ent="ybases_ent"+ech.str(),pos="ybases_pos"+ech.str();
      gposx_strippos_l[i-1].push_back(new TH2D(har.c_str(),har.c_str(),2500,0,3200,2500,0.,r_max*1.1));
      glposx_strippos_l[i-1].push_back(new TH2D(asd.c_str(),har.c_str(),2500,0,3200,2500,0.,r_max*1.1));
      ybases_solved[i-1].push_back(new TH1D(yar.c_str(),yar.c_str(),10000,-1.*r_max,r_max*1.1));
      ybases_pos[i-1].push_back(new TH1D(pos.c_str(),yar.c_str(),10000,-1.*r_max,r_max*1.1));
      ybases_ent[i-1].push_back(new TH1D(ent.c_str(),yar.c_str(),10000,-1.*r_max,r_max*1.1));
      strip_widths[i-1].push_back(new TH1D(gar.c_str(),gar.c_str(),10000,m_par->strip_width*0.5,m_par->strip_width*5.));
    }
  }
  if(ldebug) cout<<"MMT_L::built object" <<endl;
//   double yhi=m_par->H+m_par->h1,striphi=yhi/m_par->strip_width;
}

bool MMT_Loader::load_event(int event){
 if(ldebug)cout<<"MMT_L::load_event #"<<event<<endl;
  int nentries=tr->GetEntries();
  if(event<0||event>=nentries){
    cerr<<"event is "<<event<<" but should be in [0,"<<nentries<<"]...process_event abort\n";
    exit(2);
  }
  GetEntry(event);
  reset_VMMs();
  float phiEntry = 0, phiPosition = 0, etaEntry = 0, etaPosition = 0, chargeThreshold = m_par->chargeThreshold, theTheta = 0;//, avTheta = 0;
  int pdgEntry = 0;
  TLorentzVector thePart, theInfo;
  TVector3 vertex;int pdg=0;
  if(ldebug)cout<<"MMT_L::load_event---about to look at the vertex, and there are "<<TruthParticle_Barcode->size()<<" truth particle barcodes and "<<MuEntry_Particle_n<<" mu entry particles"<<endl;
  for(unsigned int j = 0; j<TruthParticle_Barcode->size(); j++){
    if(TruthParticle_Barcode->at(j) == 10001 && abs(TruthParticle_Pdg->at(j))==13){
      thePart.SetPtEtaPhiE(TruthParticle_Pt->at(j),TruthParticle_Eta->at(j),TruthParticle_Phi->at(j),TruthParticle_E->at(j));
      if(MuEntry_Particle_n>0&&j<MuEntry_Particle_n){
	pdg=TruthParticle_Pdg->at(j);
	phiEntry = MuEntry_Particle_Phi->at(j);
	etaEntry = MuEntry_Particle_Eta->at(j);
	phiPosition = MuEntry_Position_Phi->at(j);
	etaPosition = MuEntry_Position_Eta->at(j);
      }
      vertex=TVector3(TruthVertex_X->at(j),TruthVertex_Y->at(j),TruthVertex_Z->at(j));
    }
  }
//   athena_header head(thePart,TruthParticle_n,etaPosition,etaEntry,phiPosition,phiEntry, MuEntry_Particle_n,vertex);

  //mirrors the loop in Reshape_Generated_Events(C)
  if(ldebug)cout<<"MMT_L::load_event---init the fit and info entries"<<endl;
  evFit_entry fit; fit.athena_event=event;//-1;
  double theta_pos=atan(exp(-etaPosition))*2,theta_ent=atan(exp(-etaEntry))*2,phi_pos=phi_shift(phiPosition),phi_ent=phi_shift(phiEntry);
  evInf_entry primer(event,pdg,thePart.E(),thePart.Pt(),thePart.Eta(),etaPosition,etaEntry,phi_shift(thePart.Phi()),phi_pos,phi_shift(phiEntry),thePart.Theta(),theta_pos,theta_ent,theta_ent-theta_pos,TruthParticle_n,MuEntry_Particle_n,vertex);
//   primer.athena_event=event;//-1;
  vector<athena_entry> entries;

  //continuing on in Import_Athena_Events_expandedinfo_April
  if(ldebug)cout<<"MMT_L::load_event---loop through entries Digits_MM_truth*:"<<Digits_MM_truth_barcode->size()<<endl;
  for(unsigned int i = 0; i<Digits_MM_truth_barcode->size(); i++){
    string sname(Digits_MM_stationName->at(i));
    if(sname.compare("MML")!=0){
      if(debug)cout<<" The bad wedge name is "<<sname<<endl;
      primer.bad_wedge=true;break;//if not from the correct wedge type, ignore
    }
    int mult=Digits_MM_multiplet->at(i),gap=Digits_MM_gas_gap->at(i);
    //match to truth particle
    TLorentzVector truthPart;
    for(unsigned int j = 0; j<TruthParticle_Barcode->size(); j++){
      if(TruthParticle_Barcode->at(j) == Digits_MM_truth_barcode->at(i)) truthPart.SetPtEtaPhiE(TruthParticle_Pt->at(j),TruthParticle_Eta->at(j),TruthParticle_Phi->at(j),TruthParticle_E->at(j));
    }
    theTheta = truthPart.Theta();
    if(Digits_MM_stripGposY->at(i).size() == 0) continue;
    if(Digits_MM_stripGposY->at(i)[0]<-100 || Digits_MM_stripGposY->at(i)[0]>100 ) continue;
    int indForPos = -1, indforpos_tru=-1;
    float earliestTime = 100000, etime_tru=100000;
    for(unsigned int j = 0; j<Digits_MM_stripPosition->at(i).size(); j++){
      if(Digits_MM_time->at(i)[j]<earliestTime){
	if(Digits_MM_charge->at(i)[j]<chargeThreshold){
	  earliestTime=Digits_MM_time->at(i)[j];indForPos=j;
	}
      }
      if(Digits_MM_time->at(i)[j]<etime_tru){
	etime_tru=Digits_MM_time->at(i)[j];indforpos_tru=j;
      }
    }
    if(indForPos>=0){
      entries.push_back(athena_entry(mult,gap, Hits_MM_globalTime->at(i),Digits_MM_time->at(i)[indForPos],
				     TVector3(Digits_MM_truth_localPosX->at(i),Digits_MM_truth_localPosY->at(i),-999),//Digits_MM_truth_localPosZ->at(i)),
				     TVector3(Digits_MM_stripLposX->at(i)[indForPos],Digits_MM_stripLposY->at(i)[indForPos],-999),//Digits_MM_stripLposZ->at(i)[indForPos]),
				     TVector3(Digits_MM_stripGposX->at(i)[indForPos],Digits_MM_stripGposY->at(i)[indForPos],Digits_MM_stripGposZ->at(i)[indForPos]),
				     Digits_MM_charge->at(i)[indForPos],Digits_MM_stripPosition->at(i)[indForPos],Digits_MM_stationEta->at(i)));
    }
    if(indForPos!=indforpos_tru && indforpos_tru>=0){
      entries.push_back(athena_entry(mult,gap, Hits_MM_globalTime->at(i),Digits_MM_time->at(i)[indforpos_tru],
				     TVector3(Digits_MM_truth_localPosX->at(i),Digits_MM_truth_localPosY->at(i),-999),//Digits_MM_truth_localPosZ->at(i)),
				     TVector3(Digits_MM_stripLposX->at(i)[indforpos_tru],Digits_MM_stripLposY->at(i)[indforpos_tru],-999),//Digits_MM_stripLposZ->at(i)[indforpos_tru]),
				     TVector3(Digits_MM_stripGposX->at(i)[indforpos_tru],Digits_MM_stripGposY->at(i)[indforpos_tru],Digits_MM_stripGposZ->at(i)[indforpos_tru]),
				     Digits_MM_charge->at(i)[indforpos_tru],Digits_MM_stripPosition->at(i)[indforpos_tru],Digits_MM_stationEta->at(i)));
    }
  }
  if(ldebug) cout<<"Begin cuts...\n";
  int min_hits = 1,max_hits = 10000,nent=entries.size();
  //Number of hits cut
  if(!primer.bad_wedge)primer.pass_cut=true;//default is false
  if(nent<min_hits||nent>max_hits) primer.pass_cut=false;
  if(!primer.pass_cut&&debug)cout<<"event FAIL at max hit mark"<<endl;
  double theta_min = m_par->minimum_large_theta.getFloat(),theta_max =m_par->maximum_large_theta.getFloat(),phi_min = m_par->minimum_large_phi.getFloat(),phi_max = m_par->maximum_large_phi.getFloat();
  double tru_phi=phi_shift(thePart.Phi()),tru_theta=thePart.Theta();
  //Theta cut 
  if(tru_theta<theta_min||tru_theta>theta_max) primer.pass_cut=false;//*** do a theta cut?
  //Phi cut
  if(ldebug)cout<<"tru_phi is "<<tru_phi<<" which should be in ["<<phi_min<<","<<phi_max<<"]"<<endl;
  if(tru_phi<phi_min||tru_phi>phi_max) primer.pass_cut=false;
  if(!primer.pass_cut&&debug)cout<<"event FAIL at phi cut"<<endl;

  map<hdst_key,hdst_entry> targaryen;
  int fstation=0;
  for(unsigned int ient=0; ient<entries.size(); ient++){
    athena_entry examine=entries[ient];
    //DLM_NEW plane assignments
    //stated [3,2,1,0;7,6,5,4]
    int plane=(examine.multiplet-1)*4+examine.gas_gap-1;
    if(m_par->dlm_new)plane=examine.multiplet*4-examine.gas_gap;
//     if(debug) cout<<"SUBSTR CALL MMT_L--0...plane: "<<plane<<", multiplet: "<<examine.multiplet<<endl;
    int BC_id=ceil(examine.time/25.);
    TVector3 athena_tru(examine.strip_gpos.X(),examine.strip_gpos.Y()-examine.truth_lpos.Y(),examine.strip_gpos.Z());
    if(m_par->dlm_new){
      athena_tru.SetX(examine.strip_gpos.X()-examine.strip_lpos.X());
//       cerr<<"IT'S THE NEW TIME!"<<endl;
    }
    TVector3 athena_rec(examine.strip_gpos);
    //now store some variables BLC initializes; we might trim this down for efficiency later
    //the following line should be rather easy to one-to-one replace

    TVector3 truth(athena_tru.Y(),-athena_tru.X(),athena_tru.Z()), recon(athena_rec.Y(),-athena_rec.X(),athena_rec.Z());
    if(zplanes.find(truth.Z())==zplanes.end()) zplanes[truth.Z()]=1;
    else zplanes[truth.Z()]++;
    /*
    if(ldebug){ 
      cout<<"theta from tru="<<atan2(truth.Y(),truth.Z())<<", theta_pos="<<theta_pos<<"; phi from tru="<<atan2(truth.X(),truth.Y())<<", phi_pos="<<phi_pos<<endl;
      cout<<"truth prior: "; truth.Print();cout<<"recon prior: "; recon.Print();
    }
    */
    if(uvxxmod){
      xxuv_to_uvxx(truth,plane);xxuv_to_uvxx(recon,plane);
    }
//     if(ldebug){ cout<<"recon after:"<<m_par->setup.substr(plane,1)<<" "; recon.Print();cout<<endl;}

    //Cut on being in the intended wedge, i.e. no small wedge hits
    if(m_par->z_nominal.front()>truth.Z()||m_par->z_nominal.back()<truth.Z()) primer.pass_cut=false;
    if(!primer.pass_cut&&debug)cout<<"event FAIL at z cut with Z="<<truth.Z()<<" not in ["<<m_par->z_nominal.front().getFloat()<<","<<m_par->z_nominal.back().getFloat()<<"] "<<(m_par->z_nominal.front()>truth.Z())<<" "<<(m_par->z_nominal.back()<truth.Z())<<endl;

    double charge = examine.charge;//,  module_y_center = fabs(examine.strip_gpos.X());
    
    int strip_pos=examine.strip_pos, station=examine.eta_station,strip=strip_number(station,plane,strip_pos);  //theta_strip_id,module_y_center,plane); // true_x,true_y,plane);
    fstation=station;
    string schar=m_par->setup.substr(plane,1);
//     if(debug)cout<<"STATION_ETA="<<station<<", STRIP#="<<strip_pos<<", Y="<<recon.Y()<<" ("<<truth.Y()<<")"<<endl;
    if(m_par->H<recon.Y()&&m_par->diag){
      double width=m_par->strip_width.getFloat(),base=m_par->ybases[plane][station-1].getFloat(),yhere=recon.Y(),xhere=truth.X(),msl=0;//recon.Y();//(schar.compare("x")==0?recon.Y():truth.Y());
      if(schar=="u"||schar=="v")width/=cos(degtorad(m_par->stereo_degree.getFloat()));
      if(schar=="u"){
	msl=-tan(degtorad(m_par->stereo_degree.getFloat()));
      }
      if(schar=="v"){
	msl=tan(degtorad(m_par->stereo_degree.getFloat()));
      }
      ybases_solved[station-1][plane]->Fill(yhere-strip_pos*width-msl*xhere);
      double ypos=m_par->z_nominal[plane].getFloat()*tan(theta_pos),yent=m_par->z_nominal[plane]*tan(theta_ent);
      ybases_pos[station-1][plane]->Fill(ypos-strip_pos*width);ybases_ent[station-1][plane]->Fill(yent-strip_pos*width);
      strip_widths[station-1][plane]->Fill((yhere-base)/strip_pos);
      gposx_strippos_l[station-1][plane]->Fill(strip_pos,recon.Y());
      glposx_strippos_l[station-1][plane]->Fill(strip_pos,truth.Y());
    }
    int VMM_chip = Get_VMM_chip(strip);
    //we're doing everything by the variable known as "athena_event" to reflect C++ vs MATLAB indexing
    int btime=(event+1)*10+(BC_id-1);
    primer.NUV_bg_preVMM = 0;   //examine.gtime;
    double special_time = examine.time;// + (event+1)*200;//we just reset all the VMMs at an event's start right now, so the event time length bookeeping isn't necessary
    hdst_entry shaka(event,examine.gtime,charge,VMM_chip,plane,strip_pos,station,tru_theta,tru_phi,true,btime,special_time,truth,recon);//leave the rest of the info as 0's
    targaryen[shaka.entry_key()]=shaka;
//     if(debug){ cout<<"Filling targaryen slot: "; shaka.entry_key().print();}
//     keys.push_back(shaka.entry_key());
  }
  vector<hdst_key> keys,scheduled_for_deletion;
  int xhit=0,uvhit=0,strip_X_tot=0,strip_UV_tot=0;
  bool frontx=false,backx=false;
  vector<bool>plane_hit(m_par->setup.size(),false);
  for(map<hdst_key,hdst_entry>::iterator it=targaryen.begin();it!=targaryen.end();++it){
    plane_hit[it->second.plane]=true;
    if(it->second.charge<chargeThreshold)scheduled_for_deletion.push_back(it->first);
    else keys.push_back(it->first);
  }
  //remove hits that don't pass charge threshold but that we keep until this point for CT_tru purposes
  for(int ikill=0;ikill<(int)scheduled_for_deletion.size();ikill++)targaryen.erase(scheduled_for_deletion[ikill]);
  for(int ipl=0;ipl<(int)plane_hit.size();ipl++){
    if(plane_hit[ipl]){
      if(m_par->setup.substr(ipl,1)=="x"){
	xhit++;
	if(ipl<4)frontx=true;
	else backx=true;
      }
      else if(m_par->setup.substr(ipl,1)=="u"||m_par->setup.substr(ipl,1)=="v") uvhit++;
    }
  }
  primer.N_X_hits=xhit;
  primer.N_UV_hits=uvhit;
  //X and UV hits minumum cut
//   if(xhit<m_par->CT_x) primer.pass_cut=false;//return;
  //NEW guarantee at least one hit on each multiplet--guarantees sufficiently long lever arm for local slope calculation
//   if(!frontx||!backx)primer.pass_cut=false;
//   if(uvhit<m_par->CT_uv) primer.pass_cut=false;//return;
//   if(!primer.pass_cut&&debug)cout<<"event FAIL at CT_tru cut"<<endl;


  //---- Sort by time and apply VMM chip deadtime:Hits_Data_Set_Time_Raw = sortrows(Hits_Data_Set_Time_Raw,[10,14]);
  if(gen_bg){
//     vector<hdst_entry> napoleon(blc_Generate_Incoherent_Background(keys));
    vector<hdst_entry> napoleon(Generate_Incoherent_Background(keys));
    for(unsigned int i=0;i<napoleon.size();i++){
      hdst_key key=napoleon[i].entry_key();
      if(targaryen.find(key)==targaryen.end())targaryen[key]=napoleon[i];
    }
    if(debug)cout<<"Background generation complete."<<endl;
  }
  primer.N_hits_preVMM=targaryen.size();
  primer.N_hits_postVMM=0;
  map<hdst_key,hdst_entry>nicht_tot;
  for(map<hdst_key,hdst_entry>::iterator it=targaryen.begin();it!=targaryen.end();++it){
    if(!Mimic_VMM_Chip_Deadtime(it->second)){
      continue;
      if(debug){cout<<"VMM DEAD ERASE:"<<endl;it->first.print();}      
    }
    if(it->first==hdst_key())continue;
    nicht_tot[it->first]=it->second;
  }

  //*** FIGURE OUT WHAT TO DO WITH THE TIES--IS MIMIC VMM BUSTED? DO WE PLACE THE HIT REQUIREMENTS HERE? (PROBABLY)
  //the ties are a bookeeping curiosity--just overwrite, and the hit requirements do go here
  if(debug)cout<<"nicht_tot has "<<nicht_tot.size()<<" entries ("<<primer.N_hits_preVMM<<" pre deadtime).\n";
  for(map<hdst_key,hdst_entry>::iterator it=nicht_tot.begin(); it!=nicht_tot.end(); ++it){
    /*
    if(Hits_Data_Set_Time.find(it->first)!=Hits_Data_Set_Time.end()){
      cout<<"WE HAVE A TIE!"<<endl;
      it->second.print();
      Hits_Data_Set_Time[it->first].print();
      continue;
    }
    */
    if(it->first==hdst_key()){
      cout<<"ZEROED OUT KEY!"<<endl;
      it->second.print();
      continue;
    }
    int plane=it->second.plane;
    plane_hit[plane]=true;
    primer.N_hits_postVMM++;
//     if(Hits_Data_Set_Time.find(aegon)!=Hits_Data_Set_Time.end()) continue;
    Hits_Data_Set_Time[it->first]=it->second;
//     if(debug)it->second.print();
    if(m_par->setup.substr(plane,1).compare("x")==0){//if(debug)cout<<"ADD X STRIP VALUE "<<it->second.strip<<endl;
	strip_X_tot+=it->second.strip;
	
    }
    else{//if(debug)cout<<"ADD UV STRIP VALUE "<<it->second.strip<<endl;
      strip_UV_tot+=it->second.strip;
    }
//     if(debug)cout<<"(hit"<<rae<<",pl"<<rhaegar[rae].plane<<")...";
  }
  if(xhit==4&&uvhit==4){
    if(debug)cout<<"and so it's uv("<<strip_UV_tot<<") minus x("<<strip_X_tot<<")*0.25="<<0.25*(strip_UV_tot-strip_X_tot)<<endl;
    m_diff_xuv[fstation-1]->Fill(0.25*(strip_UV_tot-strip_X_tot));
  }
  //*** place any cuts on n_x, n_uv, n_postvmm here...
  Event_Info[event]=primer;
//   primer.print();
  if(debug){
    cout<<"HAHAHA\n";
    for(map<int,int>::iterator it=zplanes.begin();it!=zplanes.end();++it){
      ostringstream lar;lar<<it->first<<"("<<it->second<<")";
      cout<<setw(12)<<lar.str();
    }
    cout<<endl<<"And now the stored ones..."<<endl;
    for(int iz=0;iz<(int)m_par->z_nominal.size();iz++){
      cout<<setw(12)<<m_par->z_nominal[iz].getFloat();
    }
    cout<<endl;
  }
  int skip=(m_par->dlm_new?5000:500);
  bool do_print=debug&&false&&event%skip==0;
  if(nicht_tot.size()==m_par->setup.size()&&do_print){//&&primer.pass_cut){
    map<int,hdst_entry>fit_hits;
    map<hdst_key,hdst_entry>::iterator it=nicht_tot.begin();
    for(;it!=nicht_tot.end();++it) fit_hits[it->second.plane]=it->second;
    TString tag("20150330-es");tag+=(nicht_tot.begin()->second.station_eta);tag+="-";tag+=event;
    print_event(fit_hits,tag);
  }
  if(debug)cout<<"Event "<<event<<" did "<<(primer.pass_cut?"":"(NOT) ")<<"pass cuts."<<endl;
  return primer.pass_cut;
}

void MMT_Loader::print_event(const map<int,hdst_entry>&fit_hits,const TString& evtag)const{
  string printdir="/n/atlascode/backedup/stchan/mmtrigger/plots/evprint/evprint"+m_par->param_par().print_pars()+"_"+string(evtag.Data());
  vector<double> x_y,x_z,uv_y,uv_z,t_y,t_z;
  map<int,hdst_entry>::const_iterator it=fit_hits.begin();
  for(;it!=fit_hits.end();++it){
    hdst_info rec=it->second.entry_info(m_par);
    double y=it->second.recon.Y(),z=it->second.recon.Z(),zrec=rec.z.getFloat()*store_const()/*,yrec=rec.y*/,road=rec.slope.getFloat();
    t_y.push_back(y/z);t_z.push_back(z);
    if(m_par->setup[it->first]=='x'){
      x_y.push_back(road);x_z.push_back(zrec);
    }
    else{
      uv_y.push_back(road);uv_z.push_back(zrec);
    }
  }
  int nx=x_y.size(),nuv=uv_y.size();
  vector<double> x_ery(nx,m_par->x_error.getFloat()),uv_ery(nuv,m_par->uv_error.getFloat());//,x_erz(nx,5),uv_erz(nuv,5);
  set_style();
  TGraphErrors *xevnt=new TGraphErrors(x_z.size(),&(x_z.front()),&(x_y.front()),0,&(x_ery.front())),*uvevnt=new TGraphErrors(uv_z.size(),&(uv_z.front()),&(uv_y.front()),0,&(uv_ery.front()));
  TGraph *tru=new TGraph(t_y.size(),&(t_z.front()),&(t_y.front()));
       xevnt->SetLineColor(kRed);     xevnt->SetMarkerColor(kRed);        xevnt->SetMarkerStyle(kOpenCross);
    tru->SetLineColor(kSpring-5);  tru->SetMarkerColor(kSpring-5);         tru->SetMarkerStyle(kFullCircle);
  uvevnt->SetLineColor(kAzure+1);uvevnt->SetMarkerColor(kAzure+1);uvevnt->SetMarkerStyle(kFullTriangleDown);
  TString nomu(printdir),nomx(printdir),har(printdir);nomu+="_uv";nomx+="_x";
  xevnt->SetNameTitle(nomx,nomx);uvevnt->SetNameTitle(nomu,nomu);har+="_c";
  TCanvas *square=new TCanvas(har,har,500,500);
  TLatex l=plot_latex();har+="m";
  TMultiGraph *drawme=new TMultiGraph(har,har);
  drawme->Add(xevnt);drawme->Add(uvevnt);drawme->Add(tru);
  drawme->Draw("AP");
  if(debug)cout<<"Added graphs to multigraph....";
  drawme->GetXaxis()->SetLabelSize(0.03);
  drawme->GetYaxis()->SetLabelSize(0.03);
  drawme->GetXaxis()->SetTitle("z [mm]");
  drawme->GetYaxis()->SetTitle("slope");
  if(debug)cout<<"configured labels and titles...."<<endl;
  l.DrawLatex(0.20,0.90,evtag);
  square->Print((printdir+".pdf").c_str());
}

TStyle* MMT_Loader::set_style() const{
  TStyle *style = new TStyle("Plain","");//AtlasStyle();//
  style->SetNumberContours(100);
  style->SetPadRightMargin(0.10);
  style->SetPadLeftMargin(0.16);
  int n_pts=11;
//makes standard color palette slightly pastel (purple to red) with: crimson at red, lighter purple and 
  Double_t R[11]={1.00,0.60,0.20,0.20,0.40,0.20,0.20,0.60,1.00,1.00,0.80};
  Double_t G[11]={0.60,0.20,0.20,0.60,1.00,1.00,1.00,1.00,1.00,0.60,0.00};
  Double_t B[11]={1.00,1.00,1.00,1.00,1.00,0.60,0.20,0.20,0.20,0.20,0.00},length[11];
  for(int i=0;i<n_pts;i++)length[i]=0.1*i;
  TColor::CreateGradientColorTable(n_pts,length,R,G,B,100);
//   style->SetPalette(1,0);
  style->cd();
  return style;
}

TLatex MMT_Loader::plot_latex() const{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextSize(0.04);
//   l.SetTextFont(132);
  //label text color: black
  l.SetTextColor(kBlack);
  return l;
}

vector<hdst_entry> MMT_Loader::event_hdsts(int find_event) const{
  vector<hdst_entry> bolero;
  int fnd_entries=0;
  for(map<hdst_key,hdst_entry>::const_iterator entry=Hits_Data_Set_Time.begin(); entry!=Hits_Data_Set_Time.end(); ++entry){
    if(entry->second.event==find_event){
      if(entry->first==hdst_key()||entry->second.station_eta<1)continue;
      bolero.push_back(entry->second);
      fnd_entries++;
    }
    else if(fnd_entries>0) break;//we use the fact that maps store things according to the strict weak ordering of the key's comparison operator
  }
  return bolero;
}

void MMT_Loader::xxuv_to_uvxx(TVector3& hit,int plane)const{
  if(plane<4)return;
  else if(plane==4)hit_rot_stereo_bck(hit);//x to u
  else if(plane==5)hit_rot_stereo_fwd(hit);//x to v
  else if(plane==6)hit_rot_stereo_fwd(hit);//u to x
  else if(plane==7)hit_rot_stereo_bck(hit);//v to x
}

void MMT_Loader::hit_rot_stereo_fwd(TVector3& hit)const{
  double degree=degtorad(m_par->stereo_degree.getFloat());
  if(striphack) hit.SetY(hit.Y()*cos(degree));
  else{
    double xnew=hit.X()*cos(degree)+hit.Y()*sin(degree),ynew=-hit.X()*sin(degree)+hit.Y()*cos(degree);
    hit.SetX(xnew);hit.SetY(ynew);
  }
}

void MMT_Loader::hit_rot_stereo_bck(TVector3& hit)const{
  double degree=-degtorad(m_par->stereo_degree.getFloat());
  if(striphack) hit.SetY(hit.Y()*cos(degree));
  else{
    double xnew=hit.X()*cos(degree)+hit.Y()*sin(degree),ynew=-hit.X()*sin(degree)+hit.Y()*cos(degree);
    hit.SetX(xnew);hit.SetY(ynew);
  }
}

vector<hdst_key> MMT_Loader::event_hdst_keys(int find_event) const{
  vector<hdst_key> ravel;
  int fnd_entries=0;
  for(map<hdst_key,hdst_entry>::const_iterator entry=Hits_Data_Set_Time.begin(); entry!=Hits_Data_Set_Time.end(); ++entry){
    if(entry->second.event==find_event){
      ravel.push_back(entry->first);
      fnd_entries++;
    }
    else if(fnd_entries>0) break;//we use the fact that maps store things according to the strict weak ordering of the key's comparison operator
  }
  return ravel;
}

void MMT_Loader::clear_buffers(int event){
  int cl_event=event-2;
  if(cl_event<0) return;
//   cout<<"Clearing buffers in event "<<cl_event<<"...";
//   Event_Fit.erase(Event_Fit.find(cl_event));
  Event_Info.erase(Event_Info.find(cl_event));
//   cout<<"event_info done...";
  vector<hdst_key> kill_keys(event_hdst_keys(cl_event));
  for(unsigned int ik=0;ik<kill_keys.size();ik++) Hits_Data_Set_Time.erase(Hits_Data_Set_Time.find(kill_keys[ik]));
//   cout<<"hdst done..."<<endl;
}

void MMT_Loader::reset_VMMs(){
  int estations=m_par->n_stations_eta,npl(m_par->setup.size());
  VMM_chip_status=vector<vector<vector<bool> > >(num_VMM_per_plane,vector<vector<bool> >(npl,vector<bool>(estations,true)));
  VMM__chip_last_hit_time=vector<vector<vector<double> > >(num_VMM_per_plane,vector<vector<double> >(npl,vector<double>(estations,0)));
}

bool MMT_Loader::Mimic_VMM_Chip_Deadtime(hdst_entry& candy){//** ASK BLC IF THERE'S A WAY TO DO THIS SINGLE ENTRY
  if(candy.strip<0){
    candy.VMM_chip=0;//clear out non-hits  (i.e. empty events -- probably from cuts in generation)
    return false;
  }
  else{
    int VMM_chip=candy.VMM_chip-1,plane=candy.plane,station=candy.station_eta-1;double time=candy.time;  //true time in sub ns from simulation
    if(!VMM_chip_status[VMM_chip][plane][station]){ //is the chip active?
      if(debug)cout<<"DEAD CHIP!"<<endl;
      if(VMM__chip_last_hit_time[VMM_chip][plane][station] + VMM_deadtime <= time){
	//if not, should the chip be active?
	//if the dead time is past, reset the chip
	VMM__chip_last_hit_time[VMM_chip][plane][station]=0;
	VMM_chip_status[VMM_chip][plane][station]=true;
      }
      else{
// 	candy.VMM_chip=0;
	return false;
      }
    }
    if(VMM_chip_status[VMM_chip][plane][station]){ //is the chip active? 
      //if the chip is active, let the hit pass, shut off the chip, and save the time
      VMM_chip_status[VMM_chip][plane][station]=false;
      VMM__chip_last_hit_time[VMM_chip][plane][station]=time;
    }
  }
  return true;
}

double MMT_Loader::phi_shift(double athena_phi) const{
  return athena_phi+(athena_phi>=0?-1:1)*pi();
  return (-1.*athena_phi+(athena_phi>=0?1.5*pi():-0.5*pi()));
}
int MMT_Loader::Get_VMM_chip(int strip) const{  //Not Finished... Rough
  int strips_per_VMM = 64;
  return ceil(1.*strip/strips_per_VMM);
}

int MMT_Loader::strip_number(int station,int plane,int spos)const{
  assert(station>0&&station<=m_par->n_stations_eta);
  assert(plane>=0&&plane<(int)m_par->setup.size());
  bool do_auto=false;
  //if true do strip # (ceil(Y/strip_width); what's currently fed into the algorithm)  calculation based on evenly spaced eta assumption of stations
  double H=m_par->H.getFloat()/*,h=m_par->h1,z=m_par->z_nominal[plane],z0=m_par->z_nominal.front()*/,ybase=m_par->ybases[plane][station-1].getFloat();
  if(do_auto){
    //-log(tan(0.5(atan(y/z))))=eta
    //this is the even y spacing
    if(m_par->dlm_new) ybase=H+1100.*(station-1);
    else ybase=H+950.*(station-1);
    /*//this is the even eta spacing version
    double etalo=-log(tan(0.5*atan((h+H)/z))),etahi=-log(tan(0.5*atan(H/z))),inc=(etahi-etalo)/m_par->n_stations_eta;
    double this_eta=etalo+inc*(station-1);
    ybase=z*tan(2*atan(exp(-1.*this_eta)));
    */
  }
  double width=m_par->strip_width.getFloat(); string plane_char=m_par->setup.substr(plane,1);
//   if(plane_char.compare("u")==0||plane_char.compare("v")==0) width/=cos(degtorad(m_par->stereo_degree));
  int base_strip=ceil(ybase/width)+spos;
  return base_strip;
}

int MMT_Loader::Get_Strip_ID(double X,double Y,int plane,int&station) const{  //athena_strip_id,module_y_center,plane)  
  if(Y==-9999) return -1;
  if(ldebug)   cout<<"Strip (width="<<m_par->strip_width.getFloat()<<") for (X,Y,pl)=("<<X<<","<<Y<<","<<plane<<") is ";
  string setup(m_par->setup);
  double strip_width=m_par->strip_width.getFloat(), degree=degtorad(m_par->stereo_degree.getFloat());//,vertical_strip_width_UV = strip_width/cos(degree);
  int setl=setup.length();
  if(plane>=setl||plane<0){
    cerr<<"Pick a plane in [0,"<<setup.length()<<"] not "<<plane<<endl; exit(1);
  }
//   if(debug) cout<<"SUBSTR CALL MMT_L--2\n";
  string xuv=setup.substr(plane,1);
  if(xuv=="u"||xuv=="v"||xuv=="U"||xuv=="V")strip_width/=cos(degtorad(m_par->stereo_degree.getFloat()));
  station=-999;double yup=0;
  for(int ista=m_par->ybases[plane].size()-1;ista>=0;ista--){
    if(Y>=m_par->ybases[plane][ista].getFloat()){
      station=ista+1;//station indexing starts at 1 not 0
      yup=Y-m_par->ybases[plane][ista].getFloat();
      break;
    }
  }
  return yup*1./strip_width; 
}

void MMT_Loader::Init(TTree *t){
  if(!t){ cerr<<"Not a valid tree pointer in MMT_Loader::Init.\n";
    exit(1);
  }
  //set all the variables to zero before reading out of a new tree
  runNumber=0;eventNumber=0;TruthVertex_n=0;TruthParticle_n=0;
  TruthVertex_X=0;TruthVertex_Y=0;TruthVertex_Z=0;
  TruthParticle_Pt=0;TruthParticle_Eta=0;TruthParticle_Phi=0;TruthParticle_E=0;
  TruthParticle_Pdg=0;TruthParticle_Barcode=0;
  MuEntry_Particle_n=0;
  MuEntry_Particle_Pt=0;MuEntry_Particle_Eta=0;MuEntry_Particle_Phi=0;
  MuEntry_Position_Eta=0;MuEntry_Position_Phi=0;
  Hits_MM_n=0;
  Hits_MM_globalTime=0;
  Digits_MM=0;
  Digits_MM_multiplet=0;Digits_MM_gas_gap=0;
  Digits_MM_time=0; Digits_MM_charge=0;
  Digits_MM_stripPosition=0;
  Digits_MM_stationName=0;Digits_MM_stationEta=0;Digits_MM_stationPhi=0;
  Digits_MM_stripLposX=0;Digits_MM_stripLposY=0;//Digits_MM_stripLposZ=0;
  Digits_MM_stripGposX=0;Digits_MM_stripGposY=0;Digits_MM_stripGposZ=0;
  Digits_MM_truth_barcode=0;
  Digits_MM_truth_localPosX=0;Digits_MM_truth_localPosY=0;//Digits_MM_truth_localPosZ=0;
  Digits_MM_truth_XZ_angle=0;
  Digits_MM_stripForTrigger=0;
  Digits_MM_stripTimeForTrigger=0;
  t->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
  t->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
  t->SetBranchAddress("TruthVertex_n", &TruthVertex_n, &b_TruthVertex_n);
  t->SetBranchAddress("TruthVertex_X", &TruthVertex_X, &b_TruthVertex_X);
  t->SetBranchAddress("TruthVertex_Y", &TruthVertex_Y, &b_TruthVertex_Y);
  t->SetBranchAddress("TruthVertex_Z", &TruthVertex_Z, &b_TruthVertex_Z);
  t->SetBranchAddress("TruthParticle_n", &TruthParticle_n, &b_TruthParticle_n);
  t->SetBranchAddress("TruthParticle_Pt", &TruthParticle_Pt, &b_TruthParticle_Pt);
  t->SetBranchAddress("TruthParticle_Eta", &TruthParticle_Eta, &b_TruthParticle_Eta);
  t->SetBranchAddress("TruthParticle_Phi", &TruthParticle_Phi, &b_TruthParticle_Phi);
  t->SetBranchAddress("TruthParticle_E", &TruthParticle_E, &b_TruthParticle_E);
  t->SetBranchAddress("TruthParticle_Pdg", &TruthParticle_Pdg, &b_TruthParticle_Pdg);
  t->SetBranchAddress("TruthParticle_Barcode", &TruthParticle_Barcode, &b_TruthParticle_Barcode);
  t->SetBranchAddress("MuEntry_Particle_n", &MuEntry_Particle_n, &b_MuEntry_Particle_n);
  t->SetBranchAddress("MuEntry_Particle_Pt", &MuEntry_Particle_Pt, &b_MuEntry_Particle_Pt);
  t->SetBranchAddress("MuEntry_Particle_Eta", &MuEntry_Particle_Eta, &b_MuEntry_Particle_Eta);
  t->SetBranchAddress("MuEntry_Particle_Phi", &MuEntry_Particle_Phi, &b_MuEntry_Particle_Phi);
  t->SetBranchAddress("MuEntry_Position_Eta", &MuEntry_Position_Eta, &b_MuEntry_Position_Eta);
  t->SetBranchAddress("MuEntry_Position_Phi", &MuEntry_Position_Phi, &b_MuEntry_Position_Phi);
  t->SetBranchAddress("Hits_MM_n", &Hits_MM_n, &b_Hits_MM_n);
  t->SetBranchAddress("Hits_MM_globalTime", &Hits_MM_globalTime, &b_Hits_MM_globalTime);
  t->SetBranchAddress("Hits_MM_hitGlobalPositionX", &Hits_MM_hitGlobalPositionX, &b_Hits_MM_hitGlobalPositionX);
  t->SetBranchAddress("Hits_MM_hitGlobalPositionY", &Hits_MM_hitGlobalPositionY, &b_Hits_MM_hitGlobalPositionY);
  t->SetBranchAddress("Hits_MM_hitGlobalPositionZ", &Hits_MM_hitGlobalPositionZ, &b_Hits_MM_hitGlobalPositionZ);
  t->SetBranchAddress("Digits_MM", &Digits_MM, &b_Digits_MM);
  t->SetBranchAddress("Digits_MM_multiplet", &Digits_MM_multiplet, &b_Digits_MM_multiplet);
  t->SetBranchAddress("Digits_MM_gas_gap", &Digits_MM_gas_gap, &b_Digits_MM_gas_gap);
  t->SetBranchAddress("Digits_MM_time", &Digits_MM_time, &b_Digits_MM_time);
  t->SetBranchAddress("Digits_MM_charge", &Digits_MM_charge, &b_Digits_MM_charge);
  t->SetBranchAddress("Digits_MM_stripPosition", &Digits_MM_stripPosition, &b_Digits_MM_stripPosition);
  t->SetBranchAddress("Digits_MM_stationEta", &Digits_MM_stationEta, &b_Digits_MM_stationEta);
  t->SetBranchAddress("Digits_MM_stationPhi", &Digits_MM_stationPhi, &b_Digits_MM_stationPhi);
  t->SetBranchAddress("Digits_MM_stationName", &Digits_MM_stationName, &b_Digits_MM_stationName);
  t->SetBranchAddress("Digits_MM_stripLposX", &Digits_MM_stripLposX, &b_Digits_MM_stripLposX);
  t->SetBranchAddress("Digits_MM_stripLposY", &Digits_MM_stripLposY, &b_Digits_MM_stripLposY);
//   t->SetBranchAddress("Digits_MM_stripLposZ", &Digits_MM_stripLposZ, &b_Digits_MM_stripLposZ);
  t->SetBranchAddress("Digits_MM_stripGposX", &Digits_MM_stripGposX, &b_Digits_MM_stripGposX);
  t->SetBranchAddress("Digits_MM_stripGposY", &Digits_MM_stripGposY, &b_Digits_MM_stripGposY);
  t->SetBranchAddress("Digits_MM_stripGposZ", &Digits_MM_stripGposZ, &b_Digits_MM_stripGposZ);
  t->SetBranchAddress("Digits_MM_truth_barcode", &Digits_MM_truth_barcode, &b_Digits_MM_truth_barcode);
  t->SetBranchAddress("Digits_MM_truth_localPosX", &Digits_MM_truth_localPosX, &b_Digits_MM_truth_localPosX);
  t->SetBranchAddress("Digits_MM_truth_localPosY", &Digits_MM_truth_localPosY, &b_Digits_MM_truth_localPosY);
//   t->SetBranchAddress("Digits_MM_truth_localPosZ", &Digits_MM_truth_localPosZ, &b_Digits_MM_truth_localPosZ);
  t->SetBranchAddress("Digits_MM_truth_XZ_angle", &Digits_MM_truth_XZ_angle, &b_Digits_MM_truth_XZ_angle);
  t->SetBranchAddress("Digits_MM_stripForTrigger", &Digits_MM_stripForTrigger, &b_Digits_MM_stripForTrigger);
  t->SetBranchAddress("Digits_MM_stripTimeForTrigger", &Digits_MM_stripTimeForTrigger, &b_Digits_MM_stripTimeForTrigger);
}

void MMT_Loader::GetEntry(int i){
  if(ldebug) cout<<"MMT_L::GetEntry start!\n";
  b_runNumber->GetEntry(i); b_eventNumber->GetEntry(i);

  if(ldebug) cout<<"MMT_L::GetEntry before truth vertex!\n";
  b_TruthVertex_n->GetEntry(i);b_TruthVertex_X->GetEntry(i);b_TruthVertex_Y->GetEntry(i);b_TruthVertex_Z->GetEntry(i);

  if(ldebug) cout<<"MMT_L::GetEntry before truth particle!\n";
  b_TruthParticle_n->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before truth particle!\n";
  b_TruthParticle_Pt->GetEntry(i);b_TruthParticle_Eta->GetEntry(i);b_TruthParticle_Phi->GetEntry(i);b_TruthParticle_E->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before truth particle!\n";
  b_TruthParticle_Pdg->GetEntry(i);b_TruthParticle_Barcode->GetEntry(i);

  if(ldebug) cout<<"MMT_L::GetEntry before mu particle!\n";
  b_MuEntry_Particle_n->GetEntry(i);
  b_MuEntry_Particle_Pt->GetEntry(i);b_MuEntry_Particle_Eta->GetEntry(i);b_MuEntry_Particle_Phi->GetEntry(i);
  b_MuEntry_Position_Eta->GetEntry(i);b_MuEntry_Position_Phi->GetEntry(i);

  if(ldebug) cout<<"MMT_L::GetEntry before hits mm!\n";
  b_Hits_MM_n->GetEntry(i);b_Hits_MM_globalTime->GetEntry(i);

  if(ldebug) cout<<"MMT_L::GetEntry before digits mm0\n";
  b_Digits_MM->GetEntry(i);b_Digits_MM_multiplet->GetEntry(i);b_Digits_MM_gas_gap->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before digits mm time\n";
  b_Digits_MM_time->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before digits mm charge\n";
  b_Digits_MM_charge->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before digits mm strip pos\n";
  b_Digits_MM_stripPosition->GetEntry(i);
  b_Digits_MM_stationEta->GetEntry(i); b_Digits_MM_stationPhi->GetEntry(i); b_Digits_MM_stationName->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before digits mm lgpos\n";
  b_Digits_MM_stripLposX->GetEntry(i);b_Digits_MM_stripLposY->GetEntry(i);//b_Digits_MM_stripLposZ->GetEntry(i);
  b_Digits_MM_stripGposX->GetEntry(i);b_Digits_MM_stripGposY->GetEntry(i);b_Digits_MM_stripGposZ->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before digits mm truth\n";
  b_Digits_MM_truth_barcode->GetEntry(i);
  b_Digits_MM_truth_localPosX->GetEntry(i);b_Digits_MM_truth_localPosY->GetEntry(i);//b_Digits_MM_truth_localPosZ->GetEntry(i);
  b_Digits_MM_truth_XZ_angle->GetEntry(i);
  if(ldebug) cout<<"MMT_L::GetEntry before digits mm strip\n";
  b_Digits_MM_stripForTrigger->GetEntry(i);b_Digits_MM_stripTimeForTrigger->GetEntry(i);
}

vector<hdst_entry> MMT_Loader::Generate_Incoherent_Background(const vector<hdst_key>& keys){
  vector<hdst_entry> aegon;
  double phi_min=m_par->minimum_large_phi.getFloat(),dphi=m_par->maximum_large_phi.getFloat()-phi_min;
  if(keys.empty())return aegon;
  int event=keys.front().event;
  //so what do we need? a time and a position
  //the time we just uniformly plop in a given window for the "time" and not "gtime" variable
  //position is plane, theta, phi: pick a random plane, pick a random phi, pick an r (which is phi dependent for the limits
  double interval=thi-tlo,tbase=0.;//(event+1)*100;//time limits;the base stuff is in case legacy time keeping convention becomes important
  double gt_i=keys.front().gtime,charge=Hits_Data_Set_Time[keys.front()].charge;
  
  //calculate occupancy numbers
  int nhits=round(hit_ct_rand->GetRandom());
  for(int i=0;i<nhits;i++){
    int plane=floor(m_par->setup.size()*m_rand->Rndm());
    double stump=m_par->ybases[plane].front().getFloat(),top=m_par->h1.getFloat()+stump;
    double atime=tbase+interval*m_rand->Rndm();
    double phi=phi_min+dphi*m_rand->Rndm(),r=hit_rate->GetRandom(stump/cos(phi),top/cos(phi));
    double y=r*cos(phi),x=r*sin(phi),theta=atan(r/m_par->z_nominal[plane].getFloat());
    TVector3 truth(x,y,m_par->z_nominal[plane].getFloat()),recon(0,0,0);
    int station=-1,strip=Get_Strip_ID(x,y,plane,station),VMM_chip = Get_VMM_chip(strip),BC_id = ceil(atime/25),BC_time = (event+1)*10 + (BC_id-1);
    if(debug)cout<<"bg hit: (time,strip,plane,station,vmm)=("<<atime<<","<<strip<<","<<plane<<","<<station<<","<<VMM_chip<<")"<<endl;
    if(station<0||strip<0||VMM_chip<1)continue;
    aegon.push_back(hdst_entry(event,gt_i+atime,charge,VMM_chip,plane,strip,station,theta,phi,false,BC_time,atime,truth,recon));
//     if(debug)aegon.back().print();
  }
  return aegon;
}
/*
void MMT_Loader::Build_Occupancy_Table(){
  double r_min=m_par->H.getFloat(),r_max=m_par->H+m_par->h1,dr = 10;  int r_zones = ceil((r_max - r_min)/dr);
  double phi_min=m_par->minimum_large_phi.getFloat(),dphi=m_par->maximum_large_phi.getFloat()-phi_min;
  hit_rate=new TF1("get_r","pow(0.1*x,-2.125)",r_min,r_max);
  int BC_window=5;

  //hit_rate_per_cm2 I(r)~r^(-2.125), w/14.1 kHz/cm2\equiv I0 @ 100 cm, so I(r)=14.1 kHz/cm2 (r/1000 mm)^-2.125
  //use the dimensional variable u=r/(100 cm), dA=(dphi/2pi)*2pi r dr=(100 cm)^2*dphi u du
  //so we want tBC*int{dA I(r)}=25ns*dphi*I0*(1e4 cm2)int{du u^-1.125}=(2.5e8s)8*dphi*(1.41e8 s-1)*(u_min^-1/8-u_max^-1/8)
  hit_rate_per_BC=8*2.5*1.41*(m_par->wedge_opening_angle*pi()/180.)*(pow(0.001*r_min,-0.125)-pow(0.001*r_max,-0.125));
  Occupancy_Table.clear();
  Occupancy_Table=vector<vector<double> >(r_zones,vector<double>(3,0));
  for(int i=0;i<r_zones;i++){
    double r = r_min + i*dr + dr/2,r_lo=r-dr/2,r_hi=r+dr/2;
    Occupancy_Table[i][0]=r;
    double hit_rate_per_cm2 = 2.5*1e8*pow(r*0.1,-2.125);
    Occupancy_Table[i][2]=hit_rate_per_cm2;
    
    double Area = 0.5*0.01*(r_hi*r_hi-r_lo*r_lo)*dphi;
    double occupancy = hit_rate_per_cm2*Area/(40.e6)*BC_window;
    Occupancy_Table[i][1]=occupancy;
  }
}

vector<hdst_entry> MMT_Loader::blc_Generate_Incoherent_Background(const vector<hdst_key>& keys){
  vector<hdst_entry> aegon;
  double r_min=m_par->H.getFloat(),r_max=m_par->H.getFloat()+m_par->h1.getFloat(),dr = 10;  int r_zones = ceil((r_max - r_min)/dr);
  if((int)(Occupancy_Table.size())!=r_zones) Build_Occupancy_Table();
  double phi_min=m_par->minimum_large_phi.getFloat(),dphi=m_par->maximum_large_phi-phi_min;
  vector<double> z_nominal;
  for(int i = 0; i<m_par->z_nominal.size(); i++) z_nominal.push_back(m_par->z_nominal[i].getFloat());
  if(keys.empty())return aegon;
  int event=keys.front().event;
  int preBC=0,postBC=1;
//   int BC_0=keys.front().BC_time,
  int BC_window=postBC+preBC+1;
  double gt_i=keys.front().gtime-preBC*25.,charge=Hits_Data_Set_Time[keys.front()].charge;
  
  for(unsigned int plane=0;plane<m_par->setup.size();plane++){
    for(int i=0;i<r_zones;i++){
      double q=m_rand->Rndm();
      if(q<=Occupancy_Table[i][1]){
	double atime=gt_i+m_rand->Rndm()*BC_window*25.;//25 ns per bunch crossing
	double r=-0.5*Occupancy_Table[i][0]+m_rand->Rndm(),theta=atan(r/z_nominal[plane]),phi=phi_min+dphi*m_rand->Rndm();;
	double y=r*cos(phi-0.5*pi()),x=r*sin(0.5*pi()-phi);
	TVector3 truth(x,y,z_nominal[plane]),recon(0,0,0);
	int station=-1,strip=Get_Strip_ID(x,y,plane,station),VMM_chip = Get_VMM_chip(),BC_id = ceil(atime/25),BC_time = event*10 + (BC_id-1);
	aegon.push_back(hdst_entry(event,gt_i+atime,charge,VMM_chip,plane,strip,station,theta,phi,false,BC_time,atime,truth,recon));
	i--;
      }
    }
  }
  return aegon;
}
*/
