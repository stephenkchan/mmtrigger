
function event_mark = Finder_Control(A,BC_window,N)  %LOOK AT EVENT #11
global clock max_age
global Gate_Flags

global slope_max slope_min h
global roads

global Data_Set

global Event_Info

'Running Finder...'
Data_Set = A;%(1:1000,:);

% h=0.1;
%-----used in other parts of simulation for record keeping------
global N_coincidence_threshold_met N_Fit
N_Fit = 0;
N_coincidence_threshold_met = 0;
%---------------

event_mark = 30000;

N_fitters = 1;  %how fitters are possibly available
max_age = BC_window;
roads = ceil((slope_max - slope_min)/h);

Build_Finder(roads);
Build_Coincidence_Gate();

[n_total_hits,w] = size(Data_Set);
mark=2;   %The first line of the data set is not a hit, but rather a storage line for data set information, which is no longer utilized
BC_clock_min = Data_Set(2,10);
BC_clock_max = Data_Set(n_total_hits,10)+4;

a=0;
for BC_clock = BC_clock_min:BC_clock_max+10   % FPGA_clock = 1:10^10    %
    
    fitters_occupied = 0;  %how fits are requested this BC
    
    clock = BC_clock;    % clock = FPGA_clock;  %
    
%      %only run requested number of events
%     events = sum(Event_Info(1:Data_Set(mark,1),15));  %15th entry generally means it passed cuts... not most reliable
%     if events>=N || mark==n_total_hits
%         event_mark = Data_Set(mark,1)
%         break;
%     else
%         %-----Read data----------
        [Hits_to_Write,mark,new_hits] = Receive_Hits(mark);
%     end
    
    %-----If there is a hit, write to the finder----
    'record';
    if new_hits == 1;
        Record_Hits_in_Finder(Hits_to_Write);
    end
    
    Check_Coincidence_Gates()
    
    %----Check for track candidates-----
    for road = 1:roads  %add a mark for logic later
        if fitters_occupied == N_fitters;
            break;
        end
        if Gate_Flags(road,1)>0 && Gate_Flags(road,2)>0
                Fit_Gate(road); %sends candidate for local logic, fitting and clears road and neighbors
                N_coincidence_threshold_met = N_coincidence_threshold_met + 1;
                fitters_occupied = fitters_occupied + 1;
        end
        
    end
    
    for road = 1:roads
            if Gate_Flags(road,2)>0
                Clear_Finder(road,'hit_expiring')
            end
    end
    
end


end


%INITIALIZE BUFFER AND COINCIDENCE GATE REFERENCE TABLE
%--------------------
function Build_Coincidence_Gate()
'Building coincidence gate reference...'
%Might want to establish a heirarchy of gates
%Eg, 4X+4UV > 3+3.   Also,

global setup CT_x CT_uv
global Coincidence_Gate

Coincidence_Gate = zeros(2,2,2,2,2,2,2,2);

for p1=0:1
    for p2=0:1
        for p3=0:1
            for p4=0:1
                for p5=0:1
                    for p6=0:1
                        for p7=0:1
                            for p8=0:1
                                
                                switch setup
                                    case 'xxuvxxuv'
                                        X_count = p1+p2+p5+p6;
                                        U_count = p3+p7;
                                        V_count = p4+p8;
                                        UV_count = U_count + V_count;
                                    case 'xxuvvuxx'
                                        X_count = p1+p2+p5+p6;
                                        U_count = p3+p7;
                                        V_count = p4+p8;
                                        UV_count = U_count + V_count;
                                end
                                        
                                if X_count>=CT_x && UV_count>=CT_uv
                                    value = CT_x*10 + CT_uv;
                                else
                                    value = 0;
                                end
                                
                                Coincidence_Gate(p1+1,p2+1,p3+1,p4+1,p5+1,p6+1,p7+1,p8+1) = value;
                                
                            end
                        end
                    end
                end
            end
        end
    end
end

end
function Build_Finder(roads)
'Building Finders...'
global Finder Finder_Detail Gate_Flags

Finder = zeros(roads,8,2);   % sloperoad, plane, [hit yes/no, time_stamp]
Finder_Detail = zeros(roads,8,3);  %[strip,slope,hit_index];
Gate_Flags = zeros(roads,2);  % [coincidence threshold, hit expiring in road]

end
%-------------------


%PROCESS INCOMING HIT ADDRESSES
%--------------------
function [Hits_to_Write,mark,new_hits] = Receive_Hits(mark)
global Data_Set clock

n=0;
new_hits = 0;
Hits_to_Write = [0,0,0,0];

[n_total_hits,w] = size(Data_Set);

%Change to Data_Set(mark,14) everywhere for fpga clock

if clock >= Data_Set(mark,10)
    new_hits = 1;
    for i=mark:n_total_hits
        if Data_Set(i,10) <= clock
            n=n+1;
            hit_index = i;
            [plane,strip,slope] = Translate_Hit(Data_Set(i,:));
            Hits_to_Write(n,:) = [hit_index,plane,strip,slope];
        else
            mark=i;
            break;
        end
    end
end

end
function [plane,strip,slope] = Translate_Hit(Hit)

global strip_width z_large

z_hit = z_large;

strip = Hit(6);
plane = Hit(5);
z = z_hit(plane);

y = strip*strip_width;
slope = y/z;

end
function Record_Hits_in_Finder(Hits_to_Write)
global clock
global Finder Finder_Detail
global slope_min
global roads
global h x_error uv_error setup

for i=1:size(Hits_to_Write,1)
    
    hit_index = Hits_to_Write(i,1);
    plane = Hits_to_Write(i,2);
    strip = Hits_to_Write(i,3);
    slope = Hits_to_Write(i,4);
    
    
    %-----tolerance for smearing hits across slope roads------
    switch setup(plane)
        case 'x'
            tol = x_error;
        case {'u','v'}
            tol = uv_error;
        otherwise
            'ERROR -- not a plane in writing hit!'
    end
    
    %---slope road boundaries based on hit_slope +/- tolerance---
    s_min = slope - tol;
    s_max = slope + tol;
    
    road_max = round((s_max - slope_min)/h);
    if road_max>roads
        road_max=roads;
    end
    road_min = round((s_min - slope_min)/h);
    if road_min<1
        road_min=1;
    end
    
    %----fill buffer----
    for road = road_min:road_max
        Finder(road,plane,:) = [1,clock];
        Finder_Detail(road,plane,:) = [strip,slope,hit_index];
    end
    
end

end
%----------------------


%CREATE FLAGS FOR READOUT ON NEXT CLOCK TICK
%----------------------
function Check_Coincidence_Gates()

global roads
global Finder Gate_Flags
global max_age clock

Gate_Flags = zeros(roads,2);

Counter = sum(Finder(:,:,1),2);

for road = 1:roads
    if Counter(road)>0
        for plane = 1:8
            age = clock - Finder(road,plane,2);
            if age == max_age
                Gate_Flags(road,2) = 1;
%             elseif age > max_age
%                 Gate_Flags(road,2) = 0;
            end
        end
        
        Gate_Flags(road,1) = Check_Road_for_Coincidence(road);
        
    end
end

end
function gate_coincidence = Check_Road_for_Coincidence(road)

global Finder Coincidence_Gate

a = Finder(road,:,1) + 1;
gate_coincidence = Coincidence_Gate(a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8));

end
%----------------------



%FIT TRACK CANDIDATES
%----------------------
function Fit_Gate(road)

global roads Gate_Flags
%----------------
%Add logic here to find the "best candidate" locally
%----------------


%No longer use "road" imported here and find a max road value
%search in main functionstill necessary to ensure there is a result
for i=1:roads
    [max_value,r] = max(Gate_Flags(:,1));
    if Gate_Flags(r,2)>0
        break;
    else
        Gate_Flags(r,1)=0;
    end
end
road = r;

Track = Read_Finder(road);
Get_Fit(Track);
Clear_Finder(road,'fit_requested');

end
function Track = Read_Finder(road)

global Finder_Detail

Track(1,:) = Finder_Detail(road,:,1);
Track(2,:) = Finder_Detail(road,:,2);
Track(3,:) = Finder_Detail(road,:,3);


% for j=1:8
%             if Track_Indexes_4(i,j,1)~=0
%                 Track(:,j) = Track_Indexes_4(i,j,:);
%             elseif Track_Indexes_3(i,j,1)~=0
%                 Track(:,j) = Track_Indexes_3(i,j,:);
%             elseif Track_Indexes_2(i,j,1)~=0
%                 Track(:,j) = Track_Indexes_2(i,j,:);
%             elseif Track_Indexes_1(i,j,1)~=0
%                 Track(:,j) = Track_Indexes_1(i,j,:);
%             end
%         end

end
function Clear_Finder(road,type)

global Finder Finder_Detail Gate_Flags roads
neighbors = 1;

switch type
    case 'fit_requested'
        if road-neighbors<1
            road_min = 1;
            road_max = road + neighbors;
        elseif road+neighbors<roads
            road_max = roads;
            road_min = road - neighbors;
        else
            road_max = road + neighbors;
            road_min = road - neighbors;
        end
        
        Finder(road_min : road_max ,:,:) = 0;
        Finder_Detail(road_min : road_max ,:,:) = 0;
        Gate_Flags(road_min : road_max ,:) = 0;
        
    case 'hit_expiring'
        for plane = 1:8
            Finder(road,plane,:) = 0;
            Finder_Detail(road,plane,:) = 0;
        end
end



%
% Finder = zeros(roads,8,2);
% Finder_Detail = zeros(roads,8,3);
% Gate_Flags = zeros(roads,2);

end
%----------------------























