% Use for sample sets generated post approx. August 8

function output = Import_Athena_Events_expandedinfo(filename,startline,endline,savefile)   %

clear Hits_Data_Set_Time Event_Info Event_Info_With_Background Event_Info_Without_Background Athena_Generated_Events Hits_Data_Set_Time_With_Background Hits_Data_Set_Time_Without_Background
global z_large
global Event_Info Event_Fit
global Hits_Data_Set_Time_With_Background Hits_Data_Set_Time_Without_Background


'Importing Athena Generated Events...'
filename

Set_Parameters();
debug=0;

%Hits_Data_Set_Time layout
%first row is info -- num true events, num back_events
%[event, BC_id(1:4), L/S, L#, plane, strip#, true_theta,
%   true_phi, true_or_backgnd, time(BC), fit_theta, fit_phi, delta_theta,...]

%true starting line (if 1) because file output makes first 2 lines invalid
if startline==1
    startline=startline+2;
end

%made a change to defaults in import function

Athena_Generated_Events = Reshape_Generated_Events(importfile(filename,startline,endline));


j=0;
[m,n]=size(Athena_Generated_Events);
for i=1:m
   multilayer = Athena_Generated_Events(i,2);
    plane_in_multilayer = Athena_Generated_Events(i,3);
    if multilayer==2
        plane = 4 + plane_in_multilayer;
    else
        plane = plane_in_multilayer;
    end
    

    
    event = Athena_Generated_Events(i,1)+1;
    %athena_strip_id = Athena_Generated_Events(i,4);
    %athena_truth_xz_angle = Athena_Generated_Events(i,4);
    wild_variable = Athena_Generated_Events(i,4);
    athena_time = Athena_Generated_Events(i,5);
    athena_true_local_x = Athena_Generated_Events(i,6);
    athena_true_local_y = Athena_Generated_Events(i,7);
    athena_recon_local_x = Athena_Generated_Events(i,8);
    athena_recon_local_y = Athena_Generated_Events(i,9);
    athena_recon_global_x = Athena_Generated_Events(i,10);
    athena_recon_global_y = Athena_Generated_Events(i,11);
    athena_recon_global_z = Athena_Generated_Events(i,12);
    athena_recon_charge = Athena_Generated_Events(i,13);
    
    athena_true_x = athena_recon_global_x;  % + athena_true_local_x;  %athena_recon_global_x;  %
    athena_true_y = athena_recon_global_y - athena_true_local_y;  %-athena_recon_global_y;  %
    athena_true_z = athena_recon_global_z;
    
    athena_recon_x = athena_recon_global_x; %athena_recon_global_x - athena_recon_local_x;  %athena_recon_global_x;  %
    athena_recon_y = athena_recon_global_y; %athena_recon_global_y - athena_recon_local_y;  %-athena_recon_global_y;  %
    athena_recon_z = athena_recon_global_z;
    
    
    athena_true_eta = Event_Info(event,4);
    athena_true_theta = atan(exp(-athena_true_eta))*2;
    athena_true_phi = Event_Info(event,6);
    athena_true_muposition_phi = Event_Info(event,21);
    athena_true_muentry_phi = Event_Info(event,22);
    athena_true_muposition_eta = Event_Info(event,17);
    athena_true_muentry_eta = Event_Info(event,18);
    athena_true_muposition_theta = atan(exp(-athena_true_muposition_eta))*2;
    athena_true_muentry_theta = atan(exp(-athena_true_muentry_eta))*2;
    if athena_true_muposition_phi>=0
        true_muposition_phi = pi/2-(athena_true_muposition_phi-pi);
    else
        true_muposition_phi = pi/2-(pi-abs(athena_true_muposition_phi));
    end
    if athena_true_muentry_phi>=0
        true_muentry_phi = pi/2-(athena_true_muentry_phi-pi);
    else
        true_muentry_phi = pi/2-(pi-abs(athena_true_muentry_phi));
    end
    Event_Info(event,19) = athena_true_muposition_theta;
    Event_Info(event,20) = athena_true_muentry_theta;
    Event_Info(event,21) = true_muposition_phi;
    Event_Info(event,22) = true_muentry_phi;
    
    athena_pt = Event_Info(event,3);
    athena_E = Event_Info(event,2);
    
    true_x = athena_true_y;
    true_y = -athena_true_x;   
    true_z = athena_true_z;
    recon_x = athena_recon_y;
    recon_y = -athena_recon_x;
    recon_z = athena_recon_z;
    
    charge = athena_recon_charge;
    
    module_y_center = abs(athena_recon_global_x);
    
    BC_id = ceil(athena_time/25);
    strip = Get_Strip_ID(recon_x,recon_y,plane);  %athena_strip_id,module_y_center,plane); % true_x,true_y,plane);  % 
    VMM_chip = Get_VMM_chip(strip);
    true_theta = athena_true_theta;
    
    if athena_true_phi>=0
        true_phi = pi/2-(athena_true_phi-pi);
    else
        true_phi = pi/2-(pi-abs(athena_true_phi));
    end
    BC_time = event*10 + (BC_id-1);
    
    
    Event_Info(event,14) = 0;   %wild_variable;  %athena_truth_xz_angle;
    
    Event_Info(event,5) = true_theta;    
    
    
    if debug==1
    if plane==3 && athena_true_theta<0.3%plane==1 %athena_true_phi<0  && plane==1;  %plane == 3 || plane==7 %
        j=j+1;
        events_to_plot_1(j,:)=[athena_true_local_x,athena_true_local_y];%,true_theta];  %[athena_true_y,true_phi+athena_true_phi];  %[athena_true_y,athena_true_phi]; % %[athena_true_x,athena_true_y]; %[recon_x,recon_y];  %[athena_phi,athena_true_y];
        events_to_plot_2(j,:)=[athena_recon_local_x,athena_recon_local_y];
        events_to_plot_3(j,:)=[athena_recon_global_x,athena_recon_global_y];
        events_to_plot_4(j,:)=[athena_true_x,athena_true_y];
        events_to_plot_5(j,:)=[athena_recon_x,athena_recon_y];
        events_to_plot_6(j,:)=[recon_x,recon_y];
        events_to_plot_7(j,:)=[true_x,true_y];
        events_to_plot_8(j,:)=[true_x,true_phi]; 
        events_to_plot_9(j,:)=[athena_true_y,abs(athena_true_phi)];
        events_to_plot_10(j,:)=[athena_true_y,abs(athena_true_phi)-true_phi];
    end
    end
    
    
    Hits_Data_Set_Time(i+1,:) = [event,wild_variable,charge,VMM_chip,plane,strip,true_theta,true_phi,1,BC_time,0,0,0];
    Time(i+1,:) = athena_time + event*100;
    Cartesian_Hit_Data(i+1,:) = [true_x,true_y,true_z,recon_x,recon_y];  
end
Hits_Data_Set_Time_Raw = [Hits_Data_Set_Time,Time,Cartesian_Hit_Data];

%----  Set Event_Info(:,15)=1 by default --------
Event_Info(:,15)=1;

%-------- Apply cuts to data set ------------
CT_x = 0;
CT_uv = 0;
theta_min = 0.1722;
theta_max = 0.5254;
phi_min = 1.3028;
phi_max = 1.8388;
min_hits = 1;
max_hits = 10000;
%Cut on being in the intended wedge, i.e. no small wedge hits
Hits_Data_Set_Time_Raw = Apply_Cut(Hits_Data_Set_Time_Raw,17,z_large(1),z_large(8));
%Theta cut %ALTER THIS FOR BACKGROUND!!!!!
%Hits_Data_Set_Time_Raw = Apply_Cut(Hits_Data_Set_Time_Raw,7,theta_min,theta_max);
%Phi cut
Hits_Data_Set_Time_Raw = Apply_Cut(Hits_Data_Set_Time_Raw,8,phi_min,phi_max);
%Number of hits cut
Hits_Data_Set_Time_Raw = Apply_Cut(Hits_Data_Set_Time_Raw,-1,min_hits,max_hits);
%X and UV hits minumum cut
Hits_Data_Set_Time_Raw = Apply_Cut(Hits_Data_Set_Time_Raw,0,CT_x,CT_uv);

%----------  Add zeros for matrix columns to be filled later ----------
[m,n]=size(Hits_Data_Set_Time_Raw);
Hits_Data_Set_Time_Raw(1,:) = zeros(1,n);
Hits_Data_Set_Time_Raw(:,11:13)=zeros(m,3);
Hits_Data_Set_Time_Raw(:,20:23)=zeros(m,4);
Hits_Data_Set_Time_Raw(:,24:30)=zeros(m,7);

%------ Assemble event info preVMM ----------
Get_Event_Info_preVMM(Hits_Data_Set_Time_Raw);

%---- Sort by time and apply VMM chip deadtime ----------------
Hits_Data_Set_Time_Raw = sortrows(Hits_Data_Set_Time_Raw,[10,14]);
Hits_Data_Set_Time_Raw = Mimic_VMM_Chip_Deadtime(Hits_Data_Set_Time_Raw);

%---------- Sort by event and assemble event info postVMM -------------
Get_Event_Info_postVMM(sortrows(Hits_Data_Set_Time_Raw,1));

%---- File without background --------
Hits_Data_Set_Time_Without_Background = Hits_Data_Set_Time_Raw;
Event_Info_Without_Background = Event_Info;

% %---------- Add background and apply VMM chip deadtime ------------------
% Hits_Data_Set_Time_Raw = Add_sp_backgnd(Hits_Data_Set_Time_Raw,8);
% 
% %--------- Time-order hits and apply VMM vhip deadtime --------------
% Hits_Data_Set_Time = sortrows(Hits_Data_Set_Time_Raw,[10,14]);
% Hits_Data_Set_Time = Mimic_VMM_Chip_Deadtime(Hits_Data_Set_Time);
% 
% %---------- File with background ---------------
% Hits_Data_Set_Time_With_Background = Hits_Data_Set_Time;
% Event_Info_With_Background = Event_Info;


%--------- Remove zero-hit events from Event_Info ------------
[a,~]=size(Event_Info);
for i=1:a
    if Event_Info(i,11)<0 || Event_Info(i,12)<0 || Event_Info(i,9)<1
        Event_Info(i,15)=0;
    end
end

%--------- Dummy file with all events to later store fitting results -----
Event_Fit_Blank = Event_Fit;

%-------- Save files ----------
% save(savefile,'Hits_Data_Set_Time_With_Background','Hits_Data_Set_Time_Without_Background','Event_Info_With_Background','Event_Info_Without_Background','Event_Fit_Blank');
save(savefile,'Hits_Data_Set_Time_Without_Background','Event_Info_Without_Background','Event_Fit_Blank');


output = Hits_Data_Set_Time_Without_Background; %Hits_Data_Set_Time;


if debug==1
figure(1);
scatter(events_to_plot_1(:,1),events_to_plot_1(:,2))
figure(2);
scatter(events_to_plot_2(:,1),events_to_plot_2(:,2))
figure(3);
scatter(events_to_plot_3(:,1),events_to_plot_3(:,2))
figure(4);
scatter(events_to_plot_4(:,1),events_to_plot_4(:,2))
figure(5);
scatter(events_to_plot_5(:,1),events_to_plot_5(:,2))
figure(6);
scatter(events_to_plot_6(:,1),events_to_plot_6(:,2))
figure(7);
scatter(events_to_plot_7(:,1),events_to_plot_7(:,2))
figure(8);
scatter(events_to_plot_8(:,1),events_to_plot_8(:,2))   
figure(9);
scatter(events_to_plot_9(:,1),events_to_plot_9(:,2))
figure(10);
scatter(events_to_plot_10(:,1),events_to_plot_10(:,2))
end


end


function VMM_chip=Get_VMM_chip(strip)  %Not Finished... Rough
strips_per_VMM = 64;
VMM_chip = ceil(strip/strips_per_VMM);
end


function strip = Get_Strip_ID(X,Y,plane)  %athena_strip_id,module_y_center,plane)  
global setup strip_width stereo_degree

degree=degtorad(stereo_degree);
vertical_strip_width_UV = strip_width/cos(degree);

switch setup(plane)
    case 'x'
        strip_hit = ceil(Y*1/strip_width); 
    case 'u'
        y_hit = X*sin(degree)+Y*cos(degree);
        strip_hit = ceil(y_hit*1/strip_width);
        %strip_hit = ceil(Y/vertical_strip_width_UV);
    case 'v'
        y_hit = -X*sin(degree)+Y*cos(degree);
        strip_hit = ceil(y_hit*1/strip_width);
        %strip_hit = ceil(Y/vertical_strip_width_UV);
end


% bottom_y=zeros(1,4);
% bottom_y(1)=982;
% bottom_y(2)=bottom_y(1)+930+20;
% bottom_y(3)=bottom_y(2)+930+20;
% bottom_y(4)=bottom_y(3)+930+20;
% 
% bottom_width=zeros(1,4);
% bottom_width(1)=582.3;
% bottom_width(2)=1145.11;
% bottom_width(3)=1707.91;
% bottom_width(4)=2270.72;
% 
% 

% if module_y_center == 1447
%     module=1;
% elseif module_y_center == 2397
%     module=2;
% elseif module_y_center == 3347
%     module=3;
% elseif module_y_center == 4239.5
%     module=4;
% else
%     'Error: unknown module center'
%     strip=0;
%     return;
% end
% 
% bottom_strip = ceil([bottom_y(module)/strip_width,(bottom_y(module)-bottom_width(module)*tan(degree))/vertical_strip_width_UV,(bottom_y(module)-bottom_width(module)*tan(degree))/vertical_strip_width_UV])
% %bottom_strip = ceil([bottom_y(module)/strip_width,bottom_y(module)/vertical_strip_width_UV,bottom_y(module)/vertical_strip_width_UV])
% 
% switch setup(plane)
%     case 'x'
%         plane_type=1;
%     case 'u'
%         plane_type=2;
%     case 'v'
%         plane_type=3;
% end
% 
% strip_hit = bottom_strip(plane_type) + athena_strip_id - 1;

strip=strip_hit;
end




function output = Reshape_Generated_Events(C)
'reshaping imported matrix...'
clear D

global Event_Info Event_Fit


Event_Info = zeros(1,35);

Event_Fit = zeros(1,10);

[m,n]=size(C);

a=0;
b=0;
c=0;
for i=1:m-1
    sprintf('reshaping %i/%i',i,m)
    if C(i,1)==-1
        athena_event = C(i,4);
        a=a+1;
        c=c+1;
        
        for j=i+2:m
            if C(j,1)==-2
                break;
            else
                b=b+1;
                D(b,:) = [athena_event,C(j,:)];
            end
        end

        Event_Info(a,1) = athena_event;
        Event_Info(a,2) = C(i+1,4);
        Event_Info(a,3) = C(i+1,1);
        Event_Info(a,4) = C(i+1,2);
        Event_Info(a,6) = C(i+1,3);
        
        Event_Info(a,16) = C(i+1,5);
        Event_Info(a,17) = C(i+1,6);
        Event_Info(a,18) = C(i+1,7);
        Event_Info(a,21) = C(i+1,8);
        Event_Info(a,22) = C(i+1,9);
        Event_Info(a,23) = C(i+1,10);
        
        Event_Info(a,31) = C(i+1,11);
        Event_Info(a,32) = C(i+1,12);
        Event_Info(a,33) = C(i+1,13);

        Event_Fit(a,1) = athena_event;
                
    end
end

output = D;
end


% function Get_Event_Info_preVMM(A)
% 
% global Event_Info
% 
% [m,mm]=size(A);
% [n,nn]=size(Event_Info);
% 
% N_hits_pre_VMM = 0;
% % planes_hit = zeros(1,8);    
% event = min(A(2:m,1));
% for i=2:m
%     
%     sprintf('accumulating event pre-fit pre-VMM info %i/%i',i,m)
%     
%     if A(i,1)==event
%         
%         N_hits_pre_VMM = N_hits_pre_VMM + 1;
%         
%     end
%         
%     if A(i,1)~=event || i==m
%         
%         Event_Info(event,9) = N_hits_pre_VMM;
%         
%         N_hits_pre_VMM = 0;
%         
%         if i==m
%             break;
%         else
%             i=i-1;
%             event = event + 1;
%         end
%     end
% 
% end
% end


% function Get_Event_Info_postVMM(A)
% 
% global Event_Info
% 
% [m,mm]=size(A);
% [n,nn]=size(Event_Info);
% 
% N_hits_post_VMM = 0;
% event = min(A(2:m,1));
% planes_hit = zeros(1,8);    
% 
% for i=2:m
%     
%     sprintf('accumulating event pre-fit post-VMM info %i/%i',i,m)
%     
%     if A(i,1)==event
%         
%         N_hits_post_VMM = N_hits_post_VMM + 1;
%         planes_hit(A(i,5)) = 1;
%         
%     end
%         
%     if A(i,1)~=event || i==m
%         N_X_planes = planes_hit(1)+planes_hit(2)+planes_hit(5)+planes_hit(6);
%         N_UV_planes = planes_hit(3)+planes_hit(4)+planes_hit(7)+planes_hit(8);
%         
%         Event_Info(event,10) = N_hits_post_VMM;
%         Event_Info(event,11:12) = [N_X_planes,N_UV_planes];
% 
%         N_hits_post_VMM = 0;
%         planes_hit = zeros(1,8);
%         
%         if i==m
%             break;
%         else
%             i=i-1;
%             event = event + 1;
%         end
%     end
% 
% end
% end

function Get_Event_Info_postVMM(A)

global Event_Info

[m,mm]=size(A);
[n,nn]=size(Event_Info);

% planes_hit = zeros(1,8);
event_min = min(A(2:m,1));
event_max = max(A(2:m,1));
for event=event_min:event_max
    
    N_hits_post_VMM = 0;
    planes_hit = zeros(1,8);
    
    for i=2:m
        if A(i,1)==event
            planes_hit(A(i,5)) = 1;
            N_hits_post_VMM = N_hits_post_VMM + 1;
            
        end
        
    end
    
    N_X_planes = planes_hit(1)+planes_hit(2)+planes_hit(5)+planes_hit(6);
    N_UV_planes = planes_hit(3)+planes_hit(4)+planes_hit(7)+planes_hit(8);
    Event_Info(event,10) = N_hits_post_VMM;
    Event_Info(event,11:12) = [N_X_planes,N_UV_planes];
    
end

end


function Get_Event_Info_preVMM(A)

global Event_Info

[m,mm]=size(A);
[n,nn]=size(Event_Info);

% planes_hit = zeros(1,8);
event_min = min(A(2:m,1));
event_max = max(A(2:m,1));
for event=event_min:event_max
    
    N_hits_pre_VMM = 0;
    for i=2:m
        if A(i,1)==event
            
            N_hits_pre_VMM = N_hits_pre_VMM + 1;
%             true_phi = A(i,8);
            
        end
        
    end
    
    Event_Info(event,9) = N_hits_pre_VMM;
%     Event_Info(event,6) = true_phi;
    
end

end





function Generated_Events = importfile(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   DIGISXXUV200GEV15DEG1000ENEW = IMPORTFILE(FILENAME) Reads data from
%   text file FILENAME for the default selection.
%
%   DIGISXXUV200GEV15DEG1000ENEW = IMPORTFILE(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   digisxxuv200GeV15deg1000enew =
%   importfile('digis_xxuv_200GeV_15deg_1000e_new.txt', 2, 269395);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2014/04/30 03:08:12

%% Initialize variables.
delimiter = ' ';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = [dataArray{:,1:end-1}];
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end


%% Replace specified string with -1.0
R = cellfun(@(x) ischar(x) && strcmp(x,'<<'),raw); % Find non-numeric cells
raw(R) = {-1.0}; % Replace non-numeric cells

%% Replace specified string with -2.0
R = cellfun(@(x) ischar(x) && strcmp(x,'>>'),raw); % Find non-numeric cells
raw(R) = {-2.0}; % Replace non-numeric cells

%% Replace blank cells with 0.0
R = cellfun(@(x) isempty(x) || (ischar(x) && all(x==' ')),raw);
raw(R) = {0.0}; % Replace blank cells

%% Replace non-numeric cells with 0.0
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {0.0}; % Replace non-numeric cells

%% Create output variable
Generated_Events = cell2mat(raw);
end


% function Generated_Events = importfile(filename, startRow, endRow)
% %IMPORTFILE Import numeric data from a text file as a matrix.
% %   DIGISXXUV200GEV = IMPORTFILE(FILENAME) Reads data from text file
% %   FILENAME for the default selection.
% %
% %   DIGISXXUV200GEV = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
% %   from rows STARTROW through ENDROW of text file FILENAME.
% %
% % Example:
% %   digisxxuv200GeV = importfile('digis_xxuv_200GeV.txt', 1, 2656);
% %
% %    See also TEXTSCAN.
% 
% % Auto-generated by MATLAB on 2013/09/25 01:30:09
% 
% %% Initialize variables.
% delimiter = ' ';
% if nargin<=2
%     startRow = 1;
%     endRow = inf;
% end
% 
% %% Read columns of data as strings:
% % For more information, see the TEXTSCAN documentation.
% formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';
% 
% %% Open the text file.
% fileID = fopen(filename,'r');
% 
% %% Read columns of data according to format string.
% % This call is based on the structure of the file used to generate this
% % code. If an error occurs for a different file, try regenerating the code
% % from the Import Tool.
% dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
% for block=2:length(startRow)
%     frewind(fileID);
%     dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
%     for col=1:length(dataArray)
%         dataArray{col} = [dataArray{col};dataArrayBlock{col}];
%     end
% end
% 
% %% Close the text file.
% fclose(fileID);
% 
% %% Convert the contents of columns containing numeric strings to numbers.
% % Replace non-numeric strings with NaN.
% raw = [dataArray{:,1:end-1}];
% numericData = NaN(size(dataArray{1},1),size(dataArray,2));
% 
% for col=[1,2,3,4,5,6,7,8,9,10,11,12]
%     % Converts strings in the input cell array to numbers. Replaced non-numeric
%     % strings with NaN.
%     rawData = dataArray{col};
%     for row=1:size(rawData, 1);
%         sprintf('column %i/12, row %i/%i',col,row,endRow)
%         % Create a regular expression to detect and remove non-numeric prefixes and
%         % suffixes.
%         regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
%         try
%             result = regexp(rawData{row}, regexstr, 'names');
%             numbers = result.numbers;
%             
%             % Detected commas in non-thousand locations.
%             invalidThousandsSeparator = false;
%             if any(numbers==',');
%                 thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
%                 if isempty(regexp(thousandsRegExp, ',', 'once'));
%                     numbers = NaN;
%                     invalidThousandsSeparator = true;
%                 end
%             end
%             % Convert numeric strings to numbers.
%             if ~invalidThousandsSeparator;
%                 numbers = textscan(strrep(numbers, ',', ''), '%f');
%                 numericData(row, col) = numbers{1};
%                 raw{row, col} = numbers{1};
%             end
%         catch me
%         end
%     end
% end
% 
% 
% %% Replace specified string with -1.0
% R = cellfun(@(x) ischar(x) && strcmp(x,'<<'),raw); % Find non-numeric cells
% raw(R) = {-1.0}; % Replace non-numeric cells
% 
% %% Replace specified string with -2.0
% R = cellfun(@(x) ischar(x) && strcmp(x,'>>'),raw); % Find non-numeric cells
% raw(R) = {-2.0}; % Replace non-numeric cells
% 
% %% Replace blank cells with 0.0
% R = cellfun(@(x) isempty(x) || (ischar(x) && all(x==' ')),raw);
% raw(R) = {0.0}; % Replace blank cells
% 
% %% Replace non-numeric cells with 0.0
% R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
% raw(R) = {0.0}; % Replace non-numeric cells
% 
% %% Create output variable
% Generated_Events = cell2mat(raw);
% end
