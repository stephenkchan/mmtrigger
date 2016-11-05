#include "MMT_Finder.h"

MMT_Finder::MMT_Finder(MMT_Parameters *par):MMT_Fitter(par){
  if(debug) cout<<"MMT_Find::building finder"<< endl;
  roads=ceil(1.*(m_par->slope_max-m_par->slope_min)/m_par->h.getFloat());//initialization, can use floats
  int nplanes=m_par->setup.size();
  if(debug) cout<<"MMT_Find::finder entries " << roads << " " << m_par->slope_max.getFloat() << " " << m_par->slope_min.getFloat() << " " << m_par->h.getFloat() << endl;
  Gate_Flags=vector<vector<double> >(roads,(vector<double>(2,0)));// sloperoad,   
  //plane, [10*xhits+uvhits,hit yes/no]//[hit yes/no, time_stamp]
  Finder=vector<vector<finder_entry> >(roads,(vector<finder_entry>(nplanes,finder_entry())));  //[strip,slope,hit_index];
  if(debug) cout<<"MMT_Find::built finder"<< endl;
}

int MMT_Finder::Finder_Control(map<hdst_key,hdst_entry>& data_hdst,int BC_window, map<int,evFit_entry>& Event_Fit, map<hdst_key,hdst_entry>& Hits_Data_Set_Time){  //LOOK AT EVENT #11//***ASK BLC ABOUT EVENT MARK!
  if(debug)cout<<"Running Finder...\n";
  /*
//-----used in other parts of simulation for record keeping------
// global N_coincidence_threshold_met N_Fit
//---------------
*/
  int N_coincidence_threshold_met = 0, event_mark = 30000, N_fitters = 1;  //how many fitters are possibly available
  max_age = BC_window;
  roads = ceil((m_par->slope_max - m_par->slope_min)/m_par->h.getFloat());

//   int n_total_hits=data_hdst.size(),w=(n_total_hits==0?0:data_hdst.front().size());
  int BC_clock_min = data_hdst.begin()->second.BC_time;
  int BC_clock_max = data_hdst.rbegin()->second.BC_time+4;

  vector<Hit> hits_to_write;
  for(int BC_clock=BC_clock_min; BC_clock<BC_clock_max+10; BC_clock++){   // FPGA_clock = 1:10^10    //
    int fitters_occupied = 0;  //how fits are requested this BC
    clock = BC_clock;    // clock = FPGA_clock;  
//     Receive_Hits(hits_to_write,mark,new_hits,data_hdst);//***
    //-----If there is a hit, write to the finder----
//     if(new_hits) Record_Hits_in_Finder(hits_to_write);
    Record_Hits_in_Finder(Receive_Hits(data_hdst));
    Check_Coincidence_Gates();
    
    //----Check for track candidates-----
    for(int road=0; road<roads; road++){  //add a mark for logic later
      if(fitters_occupied == N_fitters) break;
      if(Gate_Flags[road][0]>0 && Gate_Flags[road][1]>0){
	Fit_Gate(road,Event_Fit,Hits_Data_Set_Time); //sends candidate for local logic, fitting and clears road and neighbors
	N_coincidence_threshold_met = N_coincidence_threshold_met + 1;
	fitters_occupied = fitters_occupied + 1;
      }
    }
  
    for(int road=0; road<roads; road++){  //add a mark for logic later
      if(Gate_Flags[road][1]>0) Clear_Finder(road,"hit_expiring");
    }
  }
  //event_mark=A[mark][0];
  return event_mark;
}

void MMT_Finder::sort_hdsts(vector<hdst_entry>& hdsts) const{
  map<hdst_key,hdst_entry> sorter;
  for(unsigned int ient=0; ient<hdsts.size(); ient++){
    hdst_key key(hdsts[ient].BC_time,hdsts[ient].time,hdsts[ient].gtime,hdsts[ient].VMM_chip,hdsts[ient].event);
    sorter[key]=hdsts[ient];
  }
  hdsts.clear();
  for(map<hdst_key,hdst_entry>::iterator it=sorter.begin(); it!=sorter.end(); ++it) hdsts.push_back(it->second);
}

map<pair<int,int>,finder_entry> MMT_Finder::finder_event(int event,vector<hdst_entry>& hdsts,vector<evFit_entry>& fit,vector<pair<double,double> >&mfits,double&mxl) const{
//   if(debug) cout<<"Begin finder_event"<<event<<"..."<<hdsts.size()<<" hdst entries..."<<endl;
  fit=vector<evFit_entry>(roads,evFit_entry());
  const int nfit_max=1;
  map<pair<int,int>,finder_entry> evFinder;
  if(hdsts.empty()) return evFinder;
  //receive a hit (translate the hdst entry) and record it in the finder for the event
  //Receive_Hits and Record_Hits_in_Finder functionality
  for(int ihds=0; ihds<(int)hdsts.size(); ihds++)record_hit_evFinder(evFinder,hdsts[ihds].entry_hit(m_par));

  //Check_Coincidence_Gates stuff(); since hdsts isn't empty if this code is seen, we don't need to check "new_hits"
  int nplanes=m_par->setup.size();
  //the following container performs the function of Gate_Flags[rd][0]
  //for Gate_Flags[rd][1]: search to see if there is an entry for any evFinder(rd,(all planes))
  int fits_occupied=0;
  for(int road=0; road<roads; road++){
    vector<bool> plane_is_hit; vector<Hit> track; pair<int,int>key(road,0);
    for(int plane=0; plane<nplanes; plane++){
      key.second=plane;Hit hit;
      if(evFinder.count(key)){
	plane_is_hit.push_back(true);hit=evFinder.at(key).hit;
      }
      else plane_is_hit.push_back(false);
      track.push_back(hit);
    }
    int road_num=Coincidence_Gate(plane_is_hit);//it's a hit code finder_rd_cgate(evFinder,road,track);
//     if(debug&&!track.empty()&&road_num<-10) track.front().print_track(track);
    //check for track candidates
    if(debug)cout<<"road: "<<road<<" code; "<<road_num<<", fits: "<<fits_occupied<<endl;
    if(road_num>0){
      if(debug)cout<<"currently have done "<<fits_occupied<<" fits."<<endl;
      if(fits_occupied>=nfit_max) break;
      evFit_entry candidate=fit_event(event,track,hdsts,fits_occupied,mfits);
      if(abs(candidate.fit_dtheta.getValue())<m_par->dtheta_cut){
	fit[road]=candidate;
	mxl=Get_Local_Slope(track).getValue();
	fits_occupied++;
	if(debug) cout<<"roi of fit in road "<<road<<" is "<<fit[road].fit_roi<<", with theta="<<fit[road].fit_theta.getFloat()<<", phi="<<fit[road].fit_phi.getFloat()<<endl;
      }
    }
    fit[road].hcode=road_num;
  }
//   cout<<"xent size: "<<xent.size()<<endl;
  return evFinder;
}

void MMT_Finder::record_hit_evFinder(map<pair<int,int>,finder_entry>& evFinder, const Hit& hit) const{
  float32fixed<3> tol;
  float32fixed<3> h=m_par->h.getFloat();//transforming where the point goes
  float32fixed<3> slope=hit.info.slope.getFloat();//transforming where the point goes
  hdst_key hkey=hit.key;
  int plane=hit.info.plane;//,strip=hit.strip,station=hit.station_eta;
//   if(debug) cout<<"SUBSTR CALL MMT_Finder--m_par is "<<m_par<<"; looking at plane "<<plane<<" of "<<m_par->setup.length()<<endl;
  string plane_type=m_par->setup.substr(plane,1);
  if(plane_type=="x") tol=m_par->x_error.getFloat();
  else if(plane_type=="u"||plane_type=="v") tol=m_par->uv_error.getFloat();
  else return;//if it's an unsupported plane option, don't fill
  //---slope road boundaries based on hit_slope +/- tolerance---; if min or max is out of bounds, put it at the limit
  float32fixed<3> s_min = slope - tol, s_max = slope + tol;
  
  //how is this division done in the hardware?
  int road_max = round((s_max - m_par->slope_min)/h.getFloat());
  if(road_max>=roads){
//     cout<<"we wanted to fill from theta "<<atan(s_min)<<" to "<<atan(s_max)<<", but road_max theta is "<<atan(h*roads+m_par->slope_min)<<endl;
    road_max=roads-1;
  }
  //how is this division done in the hardware?
  int road_min = round((s_min - m_par->slope_min)/h.getFloat());
  if(road_min<0)road_min=0;
  
  //----fill buffer----
//   if(debug)cout<<"RECORD HIT IN PLANE "<<plane<<" IN ROADS ["<<road_min<<","<<road_max<<"]"<<endl;
  for(int road = road_min; road<=road_max; road++){
    pair<int,int> key(road,plane);
    if(evFinder.find(key)==evFinder.end()) evFinder[key]=finder_entry(true,clock,hit);
    else if(hkey<evFinder.find(key)->second.hit.key) evFinder[key]=finder_entry(true,clock,hit);
  }
}

void MMT_Finder::Check_Coincidence_Gates(){
  bool hit=false;
  for(int road=0; road<roads; road++){
    for(unsigned int j=0; j<Finder[road].size(); j++) hit|=Finder[road][j].is_hit;//*** ASK BLC!
  }
  if(!hit) return;
  int npl=m_par->setup.size();
  for(int road=0; road<roads; road++){
    for(int plane=0;plane<npl;plane++){
      double age = clock - Finder[road][plane].clock;
      if(age == max_age) Gate_Flags[road][1] = 1;
    }
    Gate_Flags[road][0] = Check_Road_for_Coincidence(road);//this value is
  }
}

int MMT_Finder::Check_Road_for_Coincidence(int road) const{
  int npl=m_par->setup.size();
  vector<bool> a(npl,0);
  for(int i=0; i<npl; i++) a[i]=Finder[road][i].is_hit;
  return Coincidence_Gate(a);
}

int MMT_Finder::finder_rd_cgate(const map<pair<int,int>,finder_entry>& finder,int road,vector<Hit>&track)const{
  track.clear();
  bool front=false,back=false;int nx=0,nuv=0;
  for(int pl=0;pl<m_par->setup.size();pl++){
    Hit hit;
    if(finder.count(pair<int,int>(road,pl))){
      char typ=m_par->setup.at(pl);
      if(typ=='x'||typ=='X'){
	nx++;
	if(pl<4)front=true;
	else back=true;
      }
      else if(typ=='u'||typ=='v'||typ=='U'||typ=='V')nuv++;
      hit=finder.at(pair<int,int>(road,pl)).hit;
    }
    track.push_back(hit);
  }
  bool xpass=nx>=m_par->CT_x,uvpass=nuv>=m_par->CT_uv,fbpass=front&&back;
  int value = 10*nx+nuv;
  if(!xpass||!uvpass){
    value*=-1;
    if(debug)cout<<"CG hit count fail with value: "<<value<<endl;
  }
  else if(!fbpass){
    if(value>0)value*=-1;
    value-=5;
    if(debug)cout<<"CG quadruplet fail with value: "<<value<<endl;
  }
  return value;
}

int MMT_Finder::Coincidence_Gate(const vector<bool>& plane_hits) const{
  //8 for eight detector planes
  if(plane_hits.size()!=8){
    cerr << "There are " << plane_hits.size() << " entries in plane_hits in (Build_)Coincidence_Gate(), not 8...aborting\n";
    exit(0);
  }
  //Might want to establish a heirarchy of gates
  //Eg, 4X+4UV > 3+3.   Also,
  int X_count=0,U_count=0,V_count=0,value=0;bool front=false,back=false;
  //search the string
  vector<int> u_planes=m_par->q_planes("u"), x_planes=m_par->q_planes("x"), v_planes=m_par->q_planes("v");
  for(unsigned int ip=0;ip<x_planes.size();ip++){
    if(plane_hits[x_planes[ip]]){
      X_count++;
      if(x_planes[ip]<4)front=true;
      else back=true;
    }
  }
  for(unsigned int ip=0;ip<u_planes.size();ip++) U_count+=plane_hits[u_planes[ip]];
  for(unsigned int ip=0;ip<v_planes.size();ip++) V_count+=plane_hits[v_planes[ip]];

  int UV_count = U_count + V_count;
  bool xpass=X_count>=m_par->CT_x,uvpass=UV_count>=m_par->CT_uv,fbpass=front&&back;
  value = 10*X_count+UV_count;
  if(!xpass||!uvpass){
    value*=-1;
    if(debug&&value<-10)cout<<"CG hit count fail with value: "<<value<<endl;
  }
  else if(!fbpass&&X_count+UV_count>0){
    if(value>0)value*=-1;
    value-=5;
    if(debug)cout<<"CG quadruplet fail with value: "<<value<<endl;
  }
  return value;
}

vector<Hit> MMT_Finder::Receive_Hits(map<hdst_key,hdst_entry>& data_hdst) const{
  /*
    old arguments were vector<Hit>& hits_to_write, int& mark, bool& new_hits,const map<hdst_key,hdst_entry>& data_hdst
  int n=0;
  new_hits = false;
  //   Hits_to_Write = [0,0,0,0];//what?
  // [n_total_hits,w] = size(Data_Set);
  //Change to Data_Set(mark,14) everywhere for fpga clock

  if(clock >= data_hdst[mark].BC_time){
    new_hits = true;
    for(int i=mark; i<hits_to_write.size(); i++){
      if(data_hdst[i].BC_time <= clock){
	n++;
	hits_to_write[n]=Translate_Hit(i,data_hdst[i]);
      }
      else mark=i;
    }
  }
  */
  vector<Hit> vulcan;
  for(map<hdst_key,hdst_entry>::iterator it=data_hdst.begin(); it!=data_hdst.end(); ++it) vulcan.push_back(Translate_Hit(it->second));
  return vulcan;
}

Hit MMT_Finder::Translate_Hit(const hdst_entry& data) const{//*** ASK BLC WHERE PLANE COMES FROM
  /*
  if(debug){
    int plane=data.plane;//,strip=data.strip,station=data.station_eta;
    if((plane==2||plane==4)){ cout<<"@@@@@@ U! "; data.truth.Print();}
    else if((plane==3||plane==5)){ cout<<"@@@@@@ V! "; data.truth.Print();}
    else { cout<<"$$$$$$$ X! "; data.truth.Print();}
  }
  */
  return data.entry_hit(m_par);
  return Hit(data.entry_key(),data.entry_info(m_par));
}

void MMT_Finder::Record_Hits_in_Finder(const vector<Hit>& hits_to_write){
  string setup=m_par->setup;
  for(unsigned int i=0; i<hits_to_write.size(); i++){
    float32fixed<2> tol;
//     if(debug) cout<<"SUBSTR CALL MMT_Finder--1\n";
    string plane_type=setup.substr(hits_to_write[i].info.plane,1);
    if(plane_type=="x") tol=m_par->x_error;
    else if(plane_type=="u"||plane_type=="v") tol=m_par->uv_error;
    else return;

    hdst_key key=hits_to_write[i].key;hdst_info info=hits_to_write[i].info;
    float32fixed<2> slope=hits_to_write[i].info.slope;
    int plane=hits_to_write[i].info.plane;
    
    //---slope road boundaries based on hit_slope +/- tolerance---
    float32fixed<2> s_min = slope - tol, s_max = slope + tol;
    
    //how is this done in the FPGA?
    int road_max = round((s_max - m_par->slope_min.getFloat())/m_par->h.getFloat());
    if(road_max>roads) road_max=roads;
    
    //how is this done in the FPGA?
    int road_min = round((s_min - m_par->slope_min.getFloat())/m_par->h.getFloat());
    if(road_min<0)road_min=0;
    
    //----fill buffer----
    for(int road = road_min; road<=road_max; road++){
      if(Finder[road][plane]==finder_entry()) Finder[road][plane]=finder_entry(true,clock,hits_to_write[i]);
      else{
	if(key<Finder[road][plane].hit.key) Finder[road][plane]=finder_entry(true,clock,hits_to_write[i]);
      }
    }
  }
}


//FIT TRACK CANDIDATES
//----------------------
void MMT_Finder::Fit_Gate(int road, map<int,evFit_entry>& Event_Fit, map<hdst_key,hdst_entry>& Hits_Data_Set_Time){
//----------------
//Add logic here to find the "best candidate" locally
//----------------
//No longer use "road" imported here and find a max road value
//search in main functionstill necessary to ensure there is a result
//*** ASK BLC WHY THERE WAS A LOOP OVER ROADS HERE
  int r=0,maxval=0;
  for(unsigned int i=0; i<Gate_Flags.size(); i++){
    if(Gate_Flags[i][0]>maxval){
      maxval=Gate_Flags[i][0]; r=i;
    }
  }
  if(Gate_Flags[r][1]<=0) Gate_Flags[r][0]=0;
  road = r;//*** DOES THIS VARIABLE REALLY NEED TO BE PASSED ON??? IT'S IN A FOR LOOP ABOVE
  vector<Hit> Track = Read_Finder_to_Track(road);
  set_roads(road);//*** really?
  Get_Fit(Track,Event_Fit,Hits_Data_Set_Time);
  Clear_Finder(road,"fit_requested");
}

vector<Hit> MMT_Finder::Read_Finder_to_Track(int road) const{
  vector<Hit> track;
  for(int plane=0; plane<8; plane++) track.push_back(Finder[road][plane].hit);
  return track;
}

void MMT_Finder::Clear_Finder(int road, string type){
  int neighbors = 1,road_min=road-neighbors,road_max=road+neighbors;
  if(type=="fit_requested"){
    if(road-neighbors<1){
      road_min = 1;
      road_max = road + neighbors;
    }
    else if(road+neighbors<roads){
      road_max = roads;
      road_min = road - neighbors;
    }
    for(int i=road_min; i<=road_max; i++){
      for(int plane=0; plane<8; plane++) Finder[i][plane]=finder_entry();
      for(int gf=0;gf<2;gf++) Gate_Flags[i][gf]=0;
    }
  }
  else if(type=="hit_expiring"){
    for(int plane=0; plane<8; plane++) Finder[road][plane]=finder_entry();
  }
  else{
    cout << "Clear_Finder invoked for an unsupported option: " << type << "...aborting.\n";
    exit(7);
  }
}
