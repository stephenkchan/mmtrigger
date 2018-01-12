function Get_Fit(Track)

global Hits_Data_Set_Time Event_Fit 
global Generate_FPGA_Sample LG_max LG_min %sometimes used when debugging
global N_Fit

fit=0;

%@@@@@@@@@@ Begin Track Fitting @@@@@@@@@@@@

%----- UV filter --------------
Track = Filter_UV(Track);

%---- Calc global slopes and local X slope -----
M_x_global = Get_Global_Slope(Track,'X');
M_u_global = Get_Global_Slope(Track,'U');
M_v_global = Get_Global_Slope(Track,'V');
M_x_local = Get_Local_Slope(Track);

%---------Debugging------------

LG=M_x_local*M_x_global;

if LG>LG_max
    LG_max=LG;
end
if LG<LG_min
    LG_min=LG;
end

%------------------------------

%----  Calc delta theta ----------
Delta_Theta_division = Get_Delta_Theta_division(M_x_local,M_x_global,1);
Delta_Theta = Get_Delta_Theta(M_x_local,M_x_global); %Delta without division operation                    %This should have been changed, but runs on tuning h had only diviosn for both slots

%----- Calc ROI ----------
ROI = Get_ROI(M_x_global,M_u_global,M_v_global);

%----- Abandon fit if ROI comes back as out of bounds ------
if ROI(1)==-999
    Track;
    'SOMETHING IS OFF!'
    return;
end





%@@@@@@@@ Begin Info Storage for Later Analysis @@@@@@@@@@@@@@@@@@@@@@@
planes_hit_t=zeros(1,8);
planes_hit_bg=zeros(1,8);
for j=1:8
    index = Track(3,j);
    if index~=0
        event = Hits_Data_Set_Time(index,1);
        if Hits_Data_Set_Time(index,9)==1
            planes_hit_t(j)=1;
        elseif Hits_Data_Set_Time(index,9)==0
            planes_hit_bg(j)=1;
        end
    end
    
end
N_X_planes_t = planes_hit_t(1)+planes_hit_t(2)+planes_hit_t(5)+planes_hit_t(6);
N_UV_planes_t = planes_hit_t(3)+planes_hit_t(4)+planes_hit_t(7)+planes_hit_t(8);
N_X_planes_bg = planes_hit_bg(1)+planes_hit_bg(2)+planes_hit_bg(5)+planes_hit_bg(6);
N_UV_planes_bg = planes_hit_bg(3)+planes_hit_bg(4)+planes_hit_bg(7)+planes_hit_bg(8);
Event_Fit(event,6:9) = [N_X_planes_t,N_UV_planes_t,N_X_planes_bg,N_UV_planes_bg];

Event_Fit(event,10) = Delta_Theta;

for j=1:8
    if Track(1,j)~=0 && ROI(1)~=-999 && Delta_Theta~=-999 %&& Delta_Theta_division~=-999
        true_hit_index = Track(3,j);
        true_hit_event = Hits_Data_Set_Time(true_hit_index,1);
        track_fit_theta = ROI(1);
        track_fit_phi = ROI(2);
        Hits_Data_Set_Time(true_hit_index,11:13) = [track_fit_theta,track_fit_phi,Delta_Theta];
        %         delta_theta_analysis = Delta_Theta_Analysis_new(Track); Hits_Data_Set_Time(true_hit_index,20:22) = delta_theta_analysis;
        
        fit = 1;
        
        if Generate_FPGA_Sample==1
            mx=ROI(3);
            my=ROI(4);
            roi=ROI(5);
            Hits_Data_Set_Time(true_hit_index,24:30) = [M_x_global,M_u_global,M_v_global,M_x_local,mx,my,roi];
        else
            roi=0;
        end
        
        Event_Fit(true_hit_event,2:5) = [track_fit_theta,track_fit_phi,Delta_Theta_division,roi];
        
        
        if Hits_Data_Set_Time(true_hit_index,9)==1;
            if j==1
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 1;
            elseif j==2
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 10;
            elseif j==3
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 100;
            elseif j==4
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 1000;
            elseif j==5
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 10^4;
            elseif j==6
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 10^5;
            elseif j==7
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 10^6;
            elseif j==8
                Event_Fit(true_hit_event,11) = Event_Fit(true_hit_event,11) + 10^7;
            end
        elseif Hits_Data_Set_Time(true_hit_index,9)==0;
            if j==1
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 1;
            elseif j==2
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 10;
            elseif j==3
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 100;
            elseif j==4
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 1000;
            elseif j==5
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 10^4;
            elseif j==6
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 10^5;
            elseif j==7
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 10^6;
            elseif j==8
                Event_Fit(true_hit_event,12) = Event_Fit(true_hit_event,12) + 10^7;
            end
        end
        
    end
end

if fit==1
N_Fit = N_Fit+1;
end

end



function output = Filter_UV(Track)
global h

tolerance = h*2;  %Can be optimized...

pass_u=1;
pass_v=1;
if abs(Track(2,3)-Track(2,7))>tolerance;
    pass_u=0;
elseif abs(Track(2,4)-Track(2,8))>tolerance;
    pass_v=0;
end
if pass_u==0 && pass_v==1
    Track(:,[3,7])=0;
elseif pass_u==1 && pass_v==0
    Track(:,[4,8])=0;
end

output = Track;

end

function M = Get_Global_Slope(Track,type)
global setup

n_x=0;
n_u=0;
n_v=0;
for plane=1:8
    switch setup(plane)
        case 'x'
            stereo=0;
        case 'u'
            stereo=1;
        case 'v'
            stereo=-1;
    end
    
    if stereo==0
        n_x = n_x+1;
        x_planes(1,n_x)=plane;
    elseif stereo==1
        n_u = n_u+1;
        u_planes(1,n_u)=plane;
    elseif stereo==-1
        n_v = n_v+1;
        v_planes(1,n_v)=plane;
    end
end


switch type
    case 'X'
        sum=0;
        N=0;
        for j=x_planes    %[1,2,7,8]
            if Track(1,j)~=0
                %sum = sum + Track(j)/z_hit(j);
                sum = sum + Track(2,j);
                N=N+1;
            end
        end
        if N>0
            M = sum/N;
        else
            M=-999;
        end
        
    case 'U'
        sum=0;
        N=0;
        for j=u_planes   %[3,6]
            if Track(1,j)~=0
                %sum = sum + Track(j)/z_hit(j);
                sum = sum + Track(2,j);
                N=N+1;
            end
        end
        if N>0
            M = sum/N;
        else
            M=-999;
        end
        
    case 'V'
        sum=0;
        N=0;
        for j=v_planes    %[4,5]
            if Track(1,j)~=0
                %sum = sum + Track(j)/z_hit(j);
                sum = sum + Track(2,j);
                N=N+1;
            end
        end
        if N>0
            M = sum/N;
        else
            M=-999;
        end
        
end

end

function slope = Get_Local_Slope(Track)

global A_local B_local
global z_large

z_hit = z_large;

planes = [1,2,5,6];

X=[0,0,0,0];
N=0;
if Track(1,planes(1))~=0
    N=N+1;
    X(1) = 1;
end
if Track(1,planes(2))~=0
    N=N+1;
    X(2) = 1;
end
if Track(1,planes(3))~=0
    N=N+1;
    X(3) = 1;
end
if Track(1,planes(4))~=0
    N=N+1;
    X(4) = 1;
end


if X==[1,1,1,1] %[1,2,7,8]
    k=1;
elseif X==[1,1,1,0] % [1,2,7]
    k=2;
elseif X==[1,1,0,1] %[1,2,8]
    k=3;
elseif X==[1,0,1,1]  %[1,7,8]
    k=4;
elseif X==[0,1,1,1]  %[2,7,8]
    k=5;
elseif X==[1,1,0,0]  %[1,2]
    k=6;
elseif X==[1,0,1,0]  %[1,7]
    k=7;
elseif X==[1,0,0,1]  %[1,8]
    k=8;
elseif X==[0,1,1,0]  %[2,7]
    k=9;
elseif X==[0,1,0,1]  %[2,8]
    k=10;
elseif X==[0,0,1,1]  %[7,8]
    k=11;
else
    'No X hits'
    slope=-10;
    return;
end

sum_xy=0;
sum_y=0;
for i=1:4
    if X(i)==1
        plane=planes(i);
        sum_xy = sum_xy + z_hit(plane)*Track(1,plane);
        sum_y = sum_y + Track(1,plane);
    end
end

slope = A_local(k)*sum_xy - B_local(k)*sum_y;

end

function Delta_Theta = Get_Delta_Theta(M_local,M_global)

global DT_Factors

%delta_theta = theta_local_slope - theta_global_fit

Mult_Factors = DT_Factors;

% if M_global<0 || M_local<0   %if M-local doesn't point back into the detector, abandon
%     M_global
%     M_local
%     Delta_Theta=-999;
%     return;
% end
region=0;
[max_factor,a] = size(DT_Factors);
LG = M_local * M_global;
for j=1:max_factor   %number_LG_regions
    if LG <= Mult_Factors(j,1)
        region = j;
        break;
    end
end

if region==0
    Delta_Theta = -999;
    return;
end
    

Delta_Theta = (M_local - M_global) * Mult_Factors(region,2);

end

function Delta_Theta = Get_Delta_Theta_division(M_local,M_global,a)
%delta_theta = theta_local_slope - theta_global_fit
% a=1;  % Really sin(phi), but I use small angles about phi=pi/2
Delta_Theta = (M_local - M_global) / (a + M_local*M_global/a);
end

function ROI = Get_ROI(M_x,M_u,M_v)
%M_? are all global slopes

global Slope_to_ROI h_mx h_my m_y_min stereo_degree

%--- calc constants ------
b=degtorad(stereo_degree);
A=1/sin(b);
B=cot(b);

%---  slope conversion equations ---- 
m_y = M_x;
m_xu = A*M_u - B*m_y;
m_xv = -A*M_v + B*m_y;

%--- which slopes are truly present ----  
%Note that bad slopes are not necessarily 0 as I often use -999 to denote something missing
nu=1;
nv=1;
if M_u<0
    m_xu = 0;
    nu=0;
end
if M_v<0
    m_xv=0;
    nv=0;
end

%this should not happen...
if nu+nv==0
    'slope out of bounds -- Error slopes'
    ROI=[-999,-999];
    return;
end

%--- average of 2 mx slope values -----
m_x = (m_xv+m_xu)/(nu+nv);

%degbugging stuff
% if abs(m_xu-m_xv)>0.05 && (nu+nv)==2
%     M_x
%     M_u
%     M_v
%     m_xu
%     m_xv
%     m_x
%     pause
% end


%generated table is only for +x coordinates...
if m_x<0
    m_x = -m_x;
    T=1;
else
    T=0;
end

%Get m_x and m_y in parameterized values   %NEED TO THINK ABOUT THESE!!!!!!
a_x = round(m_x/h_mx) + 1;
a_y = round((m_y-m_y_min)/h_my)+1;

[a1,a2,a3]=size(Slope_to_ROI);

% Generally, this offers a reality check or cut.  The only reason a slope
% should be "out of bounds" is because it represents a weird UV combination
% -- ie. highly background influenced
if a_y>a1 || a_x>a2 || a_y<1 || a_x<1
    'slope out of bounds'
    ROI=[-999,-999];
    return;
end

%Reference slope-to-roi table
if T==1
    Fit_Vector = [Slope_to_ROI(a_y,a_x,1),pi-Slope_to_ROI(a_y,a_x,2)];
else
    Fit_Vector = [Slope_to_ROI(a_y,a_x,1),Slope_to_ROI(a_y,a_x,2)];
end

%--- More hardware realistic approach but need fine tuning ----
roi = Rough_ROI_temp(Fit_Vector);

%--- current "roi" which is not an actual roi but an approx phi and theta
ROI = [Fit_Vector,m_x,m_y,roi];

end

function roi = Rough_ROI_temp(Fit_Vector)

%temporary function to identify areas of the wedge.

global minimum_large_theta maximum_large_theta minimum_large_phi maximum_large_phi
global n_theta_rois n_phi_rois

theta = Fit_Vector(1);
phi = Fit_Vector(2);

h_theta = (maximum_large_theta - minimum_large_theta)/n_theta_rois;
h_phi = (maximum_large_phi - minimum_large_phi)/n_phi_rois;

roi_t = ceil((theta - minimum_large_theta)/h_theta);
roi_p = ceil((phi - minimum_large_phi)/h_phi);

if theta<minimum_large_theta || theta>maximum_large_theta || phi<minimum_large_phi || phi>maximum_large_phi
    roi_t = 0;
    roi_p = 0;
end

roi = roi_t * 1000 + roi_p;

end