% function Tuning_Simple()

clear all
close all

global Event_Info Hits_Data_Set_Time
global Generate_FPGA_Sample
global h uv_error x_error CT_x CT_uv CT
global Event_Info_Without_Background Event_Info_With_Background
global Hits_Data_Set_Time_Without_Background Hits_Data_Set_Time_With_Background
global Event_Fit Event_Fit_Tuning
global H_values UV_values CT_values CT_x_values CT_uv_values
global LG_max LG_min
global z_large

LG_max = 0;
LG_min = 1000;

global N_coincidence_threshold_met N_Fit

global computer

% clear Event_Fit_Tuning

computer = 'clark';

Set_Parameters();
Table_Generators();
Generate_FPGA_Sample=1;

%---- Set files, folders, and parameters --------

background =  'bgon';  %

openfolder = sprintf('/Users/%s/Desktop/eventfiles/allevents/',computer);
% openfolder = sprintf('/Users/%s/Desktop/eventfiles/multiple/',computer);
% openfolder_addon = '/Users/clark/Desktop/eventfiles/tuning/';
% savefolder = '/Users/clark/Desktop/eventfiles/tuning/';
savefolder = sprintf('/Users/%s/Desktop/eventfiles/runs/',computer);

energy = 200;  %in GeV
stereodeg = 15; % 1.5deg = 15
charge = 1000;
z_shift = 0;

openfile = sprintf('%s%iGeV_%ideg_%ie_new.mat',openfolder,energy,stereodeg,charge);

load(openfile)

switch background
    case 'shift'
        openfile = sprintf('/Users/%s/Desktop/eventfiles/shiftz.mat',computer);
        load(openfile)
        Event_Info = Event_Info_Without_Background;
        Hits_Data_Set_Time = Hits_Data_Set_Time_Without_Background;
    case 'incoherent'
        openfile = sprintf('/Users/%s/Desktop/eventfiles/Incoherent_Background.mat',computer);
        energy = 0;
        load(openfile);
        Event_Info = Event_Info_Background;
        Hits_Data_Set_Time = Master_Background_Hits;
        Hits_Data_Set_Time(:,10)=Hits_Data_Set_Time(:,1).*2+(Hits_Data_Set_Time(:,10)-Hits_Data_Set_Time(:,1)*10);
        Hits_Data_Set_Time(:,14)=Hits_Data_Set_Time(:,14)-100*Hits_Data_Set_Time(:,1)+Hits_Data_Set_Time(:,1)*2;
        Mimic_VMM_Chip_Deadtime(Hits_Data_Set_Time);
    case 'bgoff'
        load(openfile);
        Event_Info = Event_Info_Without_Background;
        Hits_Data_Set_Time = Hits_Data_Set_Time_Without_Background;
    case 'bgon'
        load(openfile)
        Event_Info = Event_Info_With_Background;
        Hits_Data_Set_Time = Hits_Data_Set_Time_With_Background;
        %         load('Background_30000events_2BCunif_dr10', 'Master_Background_Hits')
        %         Event_Info = Event_Info_Without_Background;
        %         Hits_Data_Set_Time = Mimic_VMM_Chip_Deadtime(sortrows(...
        %             [Hits_Data_Set_Time_Without_Background;Master_Background_Hits],[10,14]));
        % %         Hits_Data_Set_Time = Mimic_VMM_Chip_Deadtime(sortrows(Master_Background_Hits,[10,14]));
end



[n_events_possible,aa] = size(Event_Info);



Event_Fit = [Event_Info(:,1),zeros(n_events_possible,15)];

h = 0.0009;  %2.5*10^(-4);  %0.0009;  %8.5*10^(-4);  %   % 4* 10^(-4);   % %4*10^(-4);
uv_error = 0.0035;  %  0.004;  % %0.0035;   %3.5*10^(-3);     %3.4*10^(-3);

CT_x = 3; %3; %
CT_uv = 3;  %3; %
CT = CT_x + CT_uv;
x_error = h*0.5;
BC_window = 2;

z_large = z_large + z_shift;

[m,n]=size(Hits_Data_Set_Time);
Hits_Data_Set_Time(:,11:13)=zeros(m,3);
Hits_Data_Set_Time(:,20:23)=zeros(m,4);
Hits_Data_Set_Time(:,24:30)=zeros(m,7);

tic;
N=100000;
event_mark = Finder_Control(Hits_Data_Set_Time,BC_window,N);
%Front_Filter_3(Hits_Data_Set_Time,BC_window); event_mark = 30000;  %,N);  %faster version by ~4x
toc;

H_values = h;
UV_values = uv_error;
CT_values = CT;
CT_x_values = CT_x;
CT_uv_values = CT_uv;
Charge_values = charge;
Z_shift_values = z_shift;



N_coincidence_threshold_met
N_Fit

% LG_max
% LG_min


switch background
    case 'shift'
        savefile = sprintf('%sshiftz.mat',savefolder);
    case 'incoherent'
        savefile = sprintf('%sincoherent.mat',savefolder);
    otherwise
        savefile = sprintf('%s%iGeV_%ideg_%ie_%s_%i%i_%i.mat',savefolder,energy,stereodeg,charge,background,CT_x,CT_uv,BC_window)
end
Event_Fit = Event_Fit(1:event_mark,:);
Event_Info = Event_Info(1:event_mark,:);
event_mark
save(savefile,'Hits_Data_Set_Time','H_values','UV_values','Event_Fit','Event_Info','CT_values','CT_x_values','CT_uv_values','Charge_values','Z_shift_values')
