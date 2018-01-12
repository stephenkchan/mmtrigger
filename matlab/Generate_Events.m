
function Generated_Hits = Generate_Events(N)

clear Hits_Data_Set_Time Event_Info Event_Info
% switch athena_output
%     case 'yes'
%         Cartesian_Hits = Athena_Generated_Events(
% 
% 

global Lwedge Lplane_hit Theta Phi
global add_backgnd_tracks add_backgnd_salt_pep

N_tracks = 0;
true_count=0;
Table_of_Hits = zeros(1,9);
Cart_Hits_Record = zeros(1,5);

for event=1:2*N
    if true_count == N  %N_tracks == N
        break;
    end
    event;
    
    switch add_backgnd_tracks
        case 'yes'
            q=rand(1);
            if q>0.3   %0.2
                event_type='true';
                true=1;
            else
                event_type='backgnd';
                true=0;
            end
        case 'no'
            event_type='true';
            true=1;
    end
    
            
    True_Vector = Get_random_event(event_type);   %[0.4,pi/2+0.1];%0.15];   %
    Theta = True_Vector(1);
    Phi = True_Vector(2);
    
    
    switch add_backgnd_salt_pep    %THIS NEEDS TO BE GENERALIZED TO OTHER WEDGES
        case 'yes'
            q=rand(1);
            if q<1
            bsp1 = [N_tracks+1,1,1,1,round(random('unif',1,8)),round(random('unif',1000,7000)),Theta,Phi,0];
            bsp2 = [N_tracks+1,1,1,1,round(random('unif',1,8)),round(random('unif',1000,7000)),Theta,Phi,0];
            BSP = [bsp1;bsp2];
            BSP_true = 1;
            else
                BSP_true = 0;
            end
        case 'no'
            BSP_true = 0;
    end
    
    
    %Is there a potential hit?
    Identify_Wedge(Theta,Phi,event_type)
    if Lwedge~=-999;  %i.e. there is a hit
        
        %Generate the hit addresses
        wedge = Lwedge;
        wedge_type = 1;%'L';
        Hit_Addresses_and_Cart_Hits = Get_Simple_Event_Hit_Addresses(wedge,wedge_type,Lplane_hit);
        Hit_Addresses = Hit_Addresses_and_Cart_Hits(:,1:5);
        Cart_Hits = Hit_Addresses_and_Cart_Hits(:,6:10);
        if Hit_Addresses(1,1)~=-999
            N_tracks = N_tracks+1;
            [a,b] = size(Hit_Addresses);
            Theta_rec = Theta*ones(a,1);
            Phi_rec = Phi*ones(a,1);
            True = true*ones(a,1);
            A = N_tracks * ones(a,1);
            B = [A,Hit_Addresses];
            C = [B,Theta_rec,Phi_rec];
            D = [C,True];
            Table_of_Hits = [Table_of_Hits;D];
            Cart_Hits_Record = [Cart_Hits_Record;Cart_Hits];
            if BSP_true==1
            Table_of_Hits = [Table_of_Hits; BSP];
            Cart_Hits_Record = [Cart_Hits_Record;zeros(2,5)];
            end
            if true==1
            true_count = true_count + 1;
            end
        end
    end
end
Table_of_Hits(1,1) = true_count;
Table_of_Hits(1,2) = N_tracks;
Generated_Hits = [Table_of_Hits,Cart_Hits_Record];
end




%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%        Good and Background Tracks Generator
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function output = Get_random_event(event_type)

global H h1 z_large

theta_min = atan(H/z_large(1));
theta_max = atan((H+h1)/z_large(8));

 switch event_type
    case 'true'
        
        z=z_large(8);
        
        a_min = H;
        a_max = H+h1; 

       A = log(a_max/a_min);
        
    for i=1:1000000
        
        box_height = A/a_min;
        box_width = a_max - a_min;
        
        shift_height = A/a_max;
        
        a=random('unif',a_min,a_max);
        
        p = A/a;
        
        b=random('unif',0,box_height-shift_height);
        
        if b <= p-shift_height
            theta = atan(a/z);
            break;
        end
        
    end
    
theta = random('unif',theta_min+0.02,theta_max-0.02); %random('unif',0.14,0.55);  %0.126,0.55);
    

    case 'backgnd'   %NEED TO WORK ON THIS!!!!!
        
        z = z_large(8); % + rand(1)*4000;
        
        a_min = H;
        a_max = H+h1; 

       A = log(a_max/a_min);
        
    for i=1:1000000
        
        box_height = A/a_min;
        box_width = a_max - a_min;
        
        shift_height = A/a_max;
        
        a=random('unif',a_min,a_max);
        
        p = A/a;
        
        b=random('unif',0,box_height-shift_height);
        
        if b <= p-shift_height
            theta = atan(a/z);
            break;
        end
        
    end
    theta = random('unif',0.126,0.9);
 end
 
 
%theta = random('unif',0.2,0.5);  %0.126,0.55);
%    phi = random('unif',pi/2-0.2,pi/2+0.2);  %1.27,1.87);
    

phi = random('unif',1.27,1.87);
%       theta = 0.4;
%       phi = pi/2-.2;
%     
    output = [theta,phi];
    
end
    


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%        Is it  hit?  Which planes and wedge?
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%Is the given event on a wedge?
%Which wedge?
%Which planes are potentially hit?
%refering to chamber I, II, III,..., VIII
function output = Identify_Wedge(theta,phi,event_type)

global cartesian_hit_rot Lwedge Lplane_hit

global H w1 w2 w3 h1 h2 wedge_opening_angle z_large

%-------
lambda = (pi - degtorad(wedge_opening_angle))/2;  %tilt of side of wedges  (see paper 7)

%-------


%------- Which large wedge zone? ----------------------------

delta=0.1;  %restricting more tightly on the actual wedge

if phi>3*pi/8+delta && phi<5*pi/8-delta
    Lwedge = 1;
elseif phi>5*pi/8+delta && phi<7*pi/8-delta
    Lwedge = 2;
elseif phi>7*pi/8+delta && phi<9*pi/8-delta
    Lwedge = 3;
elseif phi>9*pi/8+delta && phi<11*pi/8-delta
    Lwedge = 4;
elseif phi>11*pi/8+delta && phi<13*pi/8-delta
    Lwedge = 5;
elseif phi>13*pi/8+delta && phi<15*pi/8-delta
    Lwedge = 6;
elseif phi>15*pi/8+delta && phi<pi/8-delta
    Lwedge = 7;
elseif 1*phi>pi/8+delta && phi<3*pi/8+delta
    Lwedge = 8;
else
    Lwedge = -999;
end

switch event_type
    case 'true'
        if Lwedge == 0
           return;
        else
           z_hit = z_large;
        end
        true_hits = Get_Cartesian_Track(theta,phi,z_hit);
    case 'backgnd'
        z_hit = z_large - 4000;  %3500;%+rand(1)*2000;
        true_hits = Get_Cartesian_Track(theta,phi,z_hit);
        true_hits(3,:) = z_large;
end
        

%true_hits = Get_Cartesian_Track(theta,phi,z_hit);
x_hit = true_hits(1,:);
y_hit = true_hits(2,:);



%----- Rotate coordinate frame from wedge zone ? to wedge zone #1 -----------------------

rot_phi = -(Lwedge-1)*2*pi/8;
vec = zeros(2,1);
for j=1:8
    vec = [x_hit(j);y_hit(j)];
    vec_rot = [cos(rot_phi),-sin(rot_phi);sin(rot_phi),cos(rot_phi)]*vec;
    x_hit(j) = vec_rot(1,1);
    y_hit(j) = vec_rot(2,1);
end

cartesian_hit_rot = [x_hit;y_hit;z_hit];

%---- Is it actually on a large wedge wedge? yes=1, no=0 ------------------

Lplanes_hit = zeros(1,8);
hit=0;
for j=1:8
    if y_hit(j)>H && y_hit(j)<(H+h2) && y_hit(j)>-tan(lambda)*x_hit(j)+H-w3/2*tan(lambda) && y_hit(j)>tan(lambda)*x_hit(j)+H-w3/2*tan(lambda)
        Lplane_hit(j) = 1;
        hit=hit+1;
    elseif y_hit(j)>H+h2 && y_hit(j)<H+h1 && x_hit(j)>-w2/2 && x_hit(j)<w2/2
        Lplane_hit(j) = 1;
        hit=hit+1;
    else
        Lplane_hit(j) = 0;
    end
    
%     %Add gaps...
%     if Lplane_hit(j) == 1
%         if j==1 || j==3 || j==5 || j==7
%             if y_hit(j)>=gap_odd_bot_1 && y_hit(j)<=gap_odd_bot_1+module_gap
%                 hit=hit-1;
%                 Lplane_hit(j) = 0;
%             elseif y_hit(j)>=gap_odd_bot_2 && y_hit(j)<=gap_odd_bot_2+module_gap
%                 hit=hit-1;
%                 Lplane_hit(j) = 0;
%             elseif y_hit(j)>=gap_odd_bot_3 && y_hit(j)<=gap_odd_bot_3+module_gap
%                 hit=hit-1;
%                 Lplane_hit(j) = 0;
%             end
%         elseif j==2 || j==4 || j==6 || j==8
%             if y_hit(j)>=gap_even_bot_1 && y_hit(j)<=gap_even_bot_1+module_gap
%                 hit=hit-1;
%                 Lplane_hit(j) = 0;
%             elseif y_hit(j)>=gap_even_bot_2 && y_hit(j)<=gap_even_bot_2+module_gap
%                 hit=hit-1;
%                 Lplane_hit(j) = 0;
%             elseif y_hit(j)>=gap_even_bot_3 && y_hit(j)<=gap_even_bot_3+module_gap
%                 hit=hit-1;
%                 Lplane_hit(j) = 0;
%             end
%         end
%     end
end

if hit==0
    Lwedge = 0;
end

% Lwedge
% Lplane_hit
% cartesian_hit_rot


end

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%        Cartesian Hit Coordinates (used by Identify_Wedge)
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function cartesian = Get_Cartesian_Track(theta,phi,z_hit)
r1=z_hit(1)/cos(theta);
r2=z_hit(2)/cos(theta);
r3=z_hit(3)/cos(theta);
r4=z_hit(4)/cos(theta);
r5=z_hit(5)/cos(theta);
r6=z_hit(6)/cos(theta);
r7=z_hit(7)/cos(theta);
r8=z_hit(8)/cos(theta);

x1=r1*sin(theta)*cos(phi);
x2=r2*sin(theta)*cos(phi);
x3=r3*sin(theta)*cos(phi);
x4=r4*sin(theta)*cos(phi);
x5=r5*sin(theta)*cos(phi);
x6=r6*sin(theta)*cos(phi);
x7=r7*sin(theta)*cos(phi);
x8=r8*sin(theta)*cos(phi);

y1=r1*sin(theta)*sin(phi);
y2=r2*sin(theta)*sin(phi);
y3=r3*sin(theta)*sin(phi);
y4=r4*sin(theta)*sin(phi);
y5=r5*sin(theta)*sin(phi);
y6=r6*sin(theta)*sin(phi);
y7=r7*sin(theta)*sin(phi);
y8=r8*sin(theta)*sin(phi);

x_hit = [x1,x2,x3,x4,x5,x6,x7,x8];
y_hit = [y1,y2,y3,y4,y5,y6,y7,y8];

cartesian = [x_hit;y_hit;z_hit];

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%        Hit Address (will call on Ionization below if requested)
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function output = Get_Simple_Event_Hit_Addresses(wedge,wedge_type,plane_hit)

clear track_addresses cart_hit

global cartesian_hit_rot
global stereo_degree setup strip_width add_ionization

degree=degtorad(stereo_degree);

track_addresses = -999*[1,1,1,1,1,1,1,1,1,1];

hits=0;
old=0;
for j=1:8
    
    if plane_hit~=0
        old=hits;
        
        switch add_ionization
            case 'yes'
                z_hit = cartesian_hit_rot(3,j);
                plane = j;   %CANNOT BE USED FOR BACKGROUND TRACKS YET
                Drifted_XY = Add_Ionization_Statistics(z_hit,plane);
                X = Drifted_XY(1);
                Y = Drifted_XY(2);
                BC = Drifted_XY(3);
                cart_hit = [cartesian_hit_rot(1,j),cartesian_hit_rot(2,j),cartesian_hit_rot(3,j),Drifted_XY(1),Drifted_XY(2)];
            case 'no'
                X = cartesian_hit_rot(1,j);
                Y = cartesian_hit_rot(2,j);
                BC = 1;
                cart_hit = [cartesian_hit_rot(1,j),cartesian_hit_rot(2,j),cartesian_hit_rot(3,j),cartesian_hit_rot(1,j),cartesian_hit_rot(2,j)];
        end
        
        switch setup(j)
            case 'x'
                strip_hit = ceil(Y*1/strip_width);
                hits=hits+1;
                
            case 'u'
                y_hit = X*sin(degree)+Y*cos(degree);
                strip_hit = ceil(y_hit*1/strip_width);
                hits=hits+1;
                
            case 'v'
                y_hit = -X*sin(degree)+Y*cos(degree);
                strip_hit = ceil(y_hit*1/strip_width);
                hits=hits+1;
                
        end
    end
    
    plane = j;
    
    if hits>old
        track_addresses(hits,:) = [BC,wedge_type,wedge,plane,strip_hit,cart_hit];
    end
    
end

output = track_addresses;

end


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%               ADD IONIZATION STATISTICS
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
function output = Add_Ionization_Statistics(z_hit,plane)

global Theta Phi ion_type s


exp_drift_gap = 5;
d=exp_drift_gap;   %exp_drift_gap;   %s;

switch ion_type
    case 'George'  %20deg
        P1=0.834;
        P2=0.981;
        P3=0.99;
        G=rand(1);
        if G<=P1
            BC=1;
        elseif G<=P2
            BC=2;
        elseif G<=P3
            BC=3;
        else
            BC=4;
        end
        
    case 'simple'
        mfp = (-1/4/cos(degtorad(20))*(1/log(1-0.834)+2/log(1-0.981)+3/log(1-0.99))*exp_drift_gap)/3;
       % mfp = 3.5311;  %(avg mfp calc from george above)
        L=d/4/cos(Theta);
        P1=1-exp(-L/mfp);
        P2=1-exp(-2*L/mfp);
        P3=1-exp(-3*L/mfp);
        G=rand(1);
        if G<=P1
            BC=1;
        elseif G<=P2
            BC=2;
        elseif G<=P3
            BC=3;
        else
            BC=4;
        end
        
    case 'geo'
        P1=0.25;
        P2=0.5;
        P3=0.75;
        G=rand(1);
        if G<=P1
            BC=1;
        elseif G<=P2
            BC=2;
        elseif G<=P3
            BC=3;
        else
            BC=4;
        end        
        
end


switch plane   %SHOULD MAKE THIS DYNAMIC
    case {1,3,5,7}
        z = z_hit + (BC-1)*d/4;
        r_b = (z + d/4)/cos(Theta);
        r_a = z/cos(Theta);
        
        %r = r_a+random('poiss',1/mfp);  %1/2.6);   %2.6mm is mfp for 90:10 Ar:CO2
        r = random('unif',r_a,r_b);
    case {2,4,6,8}
        z= z_hit - (BC-1)*d/4;
        r_a = z/cos(Theta);
        r_b = (z - d/4)/cos(Theta);
        
        %r = r_a-random('poiss',1/mfp);  %1/2.6);   %2.6mm is mfp for 90:10 Ar:CO2
        r = random('unif',r_b,r_a);
end

x = r*sin(Theta)*cos(Phi);
y = r*sin(Theta)*sin(Phi);

output = [x,y,BC];

end