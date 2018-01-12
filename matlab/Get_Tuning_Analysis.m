function Get_Tuning_Analysis()

clear all
% clear Tuning_Analysis

global Event_Info Event_Fit_Tuning  % Event_Fit Event_Fit_Blank
global H_values UV_values CT_values CT_x_values CT_uv_values
global z_large strip_width
global computer intrinsic
global Data_Trimmed
global event_mark
Set_Parameters();

computer = 'clark';

intrinsic = 1;

tunetype = 'run';
background = 'bgon';  %'bgon';
data=44;  %XhitsUVhits

compare = 'compare_on';

close all

if data==33
max_hits= 100;
min_hits = 6;
min_X_hits = 3;
min_UV_hits = 3;
elseif data==44
max_hits= 8;
min_hits = 8;
min_X_hits = 4;
min_UV_hits = 4;
elseif data==32
max_hits= 100;
min_hits = 5;
min_X_hits = 3;
min_UV_hits = 2;
elseif data==22
max_hits= 100;
min_hits = 4;
min_X_hits = 2;
min_UV_hits = 2;
end


energy = 200;  %in GeV
stereodeg = 15; % 1.5deg = 15
charge = 1000; %1000;

theta_min = 0.1722;
theta_max = 0.5254;
theta_bin_size = (theta_max-theta_min)/3;
phi_min = 1.3028;
phi_max = 1.8388;




switch tunetype
    case 'run'
        CT_x = 3;
        CT_uv = 3;
        charge = 1000;
        BC_window=2;
        z_shift=50;
        
        openfolder = sprintf('/Users/%s/Desktop/eventfiles/allevents/',computer);
%         openfolder = sprintf('/Users/%s/Desktop/eventfiles/multiple/',computer);
        openfolder2 = sprintf('/Users/%s/Desktop/eventfiles/runs/',computer);
        savefolder = sprintf('/Users/%s/Desktop/eventfiles/runs/analysis/',computer);
%         
       openfile = sprintf('%s%iGeV_%ideg_%ie_new.mat',openfolder,energy,stereodeg,charge);
       openfile2 = sprintf('%s%iGeV_%ideg_%ie_%s_%i%i_%i.mat',openfolder2,energy,stereodeg,charge,background,CT_x,CT_uv,BC_window)
%         openfile = sprintf('/Users/clark/Desktop/eventfiles/shiftz.mat');
%         openfile2 = sprintf('%sshiftz.mat',openfolder2);
%         
        load(openfile,'Event_Info_Without_Background','Event_Info_With_Background')
        %         openfile2 = sprintf('%s%iGeV_%ideg_%ie_%s_%i%i.mat',openfolder2,energy,stereodeg,charge,background,CT_x,CT_uv);
        load(openfile2,'H_values','UV_values','Event_Fit','CT_values','CT_x_values','CT_uv_values');
        savefile = sprintf('%s%iGeV_%ideg_%ie_%s.mat',savefolder,energy,stereodeg,charge,background);
        
        Event_Fit_Tuning = Event_Fit;
        Tune = [theta_min+theta_bin_size/2,theta_min+3*theta_bin_size/2,theta_max-theta_bin_size/2];
                
        [event_mark,u] = size(Event_Fit)
        
        %         switch compare
        %             case 'compare_on'
        %            openfile3 = sprintf('%s%iGeV_%ideg_%ie_%s_%i%i_%i.mat',openfolder2,energy,stereodeg,charge,'bgoff',CT_x,CT_uv,BC_window);
        %         load(openfile3,'Event_Fit');
        %
        %
        
        
    otherwise
        openfolder = sprintf('/Users/%s/Desktop/eventfiles/',computer);
        openfolder2 = sprintf('/Users/%s/Desktop/eventfiles/tuning/',computer);
        savefolder = sprintf('/Users/%s/Desktop/eventfiles/tuning/analysis/',computer);
        
        openfile = sprintf('%s%iGeV_%ideg_%ie.mat',openfolder,energy,stereodeg,charge);
        load(openfile,'Event_Info_Without_Background','Event_Info_With_Background')
        
        openfile2 = sprintf('%s%iGeV_%ideg_%ie_%s_tune%s.mat',openfolder2,energy,stereodeg,charge,background,tunetype);
        load(openfile2,'Event_Fit_Tuning','H_values','UV_values','Charge_values','Z_shift_values');
        
        savefile = sprintf('%s%iGeV_%ideg_%ie_%s_tune%s.mat',savefolder,energy,stereodeg,charge,background,tunetype);
        
        switch tunetype
            case 'h'
                Tune = H_values;
            case 'uv'
                Tune = UV_values;
            case 'charge'
                Tune = [0:500:10000];  %Charge_values;
            case 'z_shift'
                Tune = Z_shift_values;
        end
end




% switch background
%     case 'bgoff'
%         Event_Info = Event_Info_Without_Background;
%
%     case 'bgon'
%         Event_Info = Event_Info_With_Background;
% end

% Event_Info = Fix_Event_Info_Phi(Event_Info,6);
% Event_Info = Fix_Event_Info_Phi(Event_Info,21);
% Event_Info = Fix_Event_Info_Phi(Event_Info,22);
% Event_Info = add_to_event_info(Event_Info);
% [a,b]=size(Event_Info);
% if a<30000
%     c=30000-a;
%     Event_Info = [Event_Info;zeros(c,b)];
% end

[a,b,N]=size(Event_Fit_Tuning);
if a<30000
    c=30000-a;
    
    Event_Fit_Tuning(a+1:event_mark,:,:) = 0;
end

events=event_mark;


% for i=1:N
%
%     Event_Fit = Event_Fit_Tuning(:,:,i);
%     Event_Analysis = [Event_Info(:,1),zeros(event_mark,15)];
%
%     theta_error = Event_Fit(:,2) - Event_Info(:,5);
%     phi_error = Event_Fit(:,3) - Event_Info(:,6);
%     dtheta_error = Event_Fit(:,10) - Event_Info(:,7);
%
%     Event_Analysis(:,2:4) = [theta_error,phi_error,dtheta_error];
%
%     Event_Analysis_Tuning(:,:,i) = Event_Analysis;
%
% end


Tuning_Analysis = zeros(N,20);


for i=1:N
    
    switch tunetype
        case 'charge'
            charge = Tune(i);
    end
    
    Event_Info = Event_Info_Alterations(charge,background,openfile,event_mark);
    [rr,ss]=size(Event_Info)
    if ss<33
        Event_Info = [Event_Info,zeros(event_mark,9)];
    end
    
    
    Event_Fit = Event_Fit_Tuning(:,:,i);
    Event_Analysis = [Event_Info(:,1),zeros(event_mark,15)];
    if intrinsic==1
        theta_error = Event_Fit(:,2) - Event_Info(:,19);  %Event_Fit(:,2) - Event_Info(:,5); %   %
        phi_error = Event_Fit(:,3) - Event_Info(:,21);   %Event_Fit(:,3) - Event_Info(:,6);
    elseif intrinsic==0
        theta_error = Event_Fit(:,2) - Event_Info(:,5);
        phi_error = Event_Fit(:,3) - Event_Info(:,6);
    end
    dtheta_error = Event_Fit(:,10) - Event_Info(:,7);
    X_fitted = Event_Fit(:,6);
    UV_fitted = Event_Fit(:,7);
    background_X_fitted = Event_Fit(:,8);
    background_UV_fitted = Event_Fit(:,9);
    Event_Analysis(:,2:4) = [theta_error,phi_error,dtheta_error];
    Event_Analysis(:,10:13) = [X_fitted,UV_fitted,background_X_fitted,background_UV_fitted];
    
    
    qual_events=0;
    nxt0=0;
    nxt1=0;
    nxt2=0;
    nxt3=0;
    nxt4=0;
    nuvt0=0;
    nuvt1=0;
    nuvt2=0;
    nuvt3=0;
    nuvt4=0;
    nxb0=0;
    nxb1=0;
    nxb2=0;
    nxb3=0;
    nxb4=0;
    nuvb0=0;
    nuvb1=0;
    nuvb2=0;
    nuvb3=0;
    nuvb4=0;
    
    n_events_with_bg = 0;
    bg_track_hits=0;
    true_track_hits=0;
    
    n_fit=0;
    n_not_fit=0;
    
    n_fit_good=0;
    
    n_fit1=0;
    n_not_fit1=0;
    n_fit2=0;
    n_not_fit2=0;
    n_fit3=0;
    n_not_fit3=0;
    
    Theta_Error=[0,0];
    Theta_Error1=[0,0];
    Theta_Error2=[0,0];
    Theta_Error3=[0,0];
    Phi_Error=[0,0];
    Phi_Error1=[0,0];
    Phi_Error2=[0,0];
    Phi_Error3=[0,0];
    Dtheta_Error=[0,0];
    Dtheta_Error1=[0,0];
    Dtheta_Error2=[0,0];
    Dtheta_Error3=[0,0];
    
    Phi=[0,0];
    Theta=[0,0];
    Dtheta=[0,0];
    Dtheta_division=[0,0];
    
    drift_sum=0;
    slope_spread_X_sum=0;
    slope_spread_UV_sum=0;
    
    q=0;
    for event=1:events
        if Event_Info(event,15)==1 %passed importing cuts which are not restrictive on X and UV
            if Event_Info(event,16)==1  % number of truth particles
                if Event_Info(event,23)==1  %number mu enter NSW
                    if Event_Info(event,9)>=min_hits && Event_Info(event,9)<=max_hits %number of hits pre-VMM and w/o bg
                        if Event_Info(event,11)>=min_X_hits  %minimum number of X hits in event
                            if Event_Info(event,12)>=min_UV_hits %minimum number of UV hits in event
                                if 1==1 %Event_Analysis(event,3)>-0.8;  %abs(Event_Info(event,6)-pi/2)>0.1   % Event_Fit(event,6)+Event_Fit(event,7)==8 %
                                    qual_events = qual_events + 1;
                                    Event_Analysis(event,5) = 1;
                                    
%                                     
                                    if Event_Fit(event,2)==0  %abs(Event_Analysis(event,3))>1 %WHYYYYYYDWQF
                                        q=q+1;
                                        event;
%                                         Event_Analysis(event,3);
%                                         Event_Info(event,9);
                                        D(q,:)= [Event_Info(event,6),Event_Info(event,5),Event_Info(event,25),Event_Info(event,26)];
%                                         %[Event_Info(event,6),Event_Fit(event,3),Event_Info(event,5),Event_Fit(event,2)]; %Event_Info(event,5)];
                                        %                                     Event_Analysis(event,12)*100+Event_Analysis(event,13)
                                    end
                                    
                                    theta = Event_Info(event,19); %Event_Info(event,5);
                                    theta_error = Event_Analysis(event,2);
                                    phi_error = Event_Analysis(event,3);
                                    dtheta_error = Event_Analysis(event,4);
                                    if Event_Fit(event,2)>0 && Event_Fit(event,3)>0 %&& phi_error<0.02 && theta_error<0.02 && dtheta_error<0.01;
                                        n_fit = n_fit + 1;
                                        Event_Analysis(event,6) = 1;
                                        wild_1 = Event_Info(event,9);
                                        wild_2 = Event_Info(event,25);
                                        wild_3 = Event_Info(event,24);
                                        Theta_Error(n_fit,:) = [event,theta_error];
                                        Phi_Error(n_fit,:) = [event,phi_error];
                                        Dtheta_Error(n_fit,:) = [event,dtheta_error];
                                        Wild(n_fit,:) = [wild_1,wild_2,wild_3];
                                        
                                        
                                        fit_phi = Event_Fit(event,3);
                                        fit_theta = Event_Fit(event,2);
                                        fit_dtheta = Event_Fit(event,10);
                                        fit_dtheta_division = Event_Fit(event,4);
                                        Phi(n_fit,:) = [event,fit_phi];
                                        Theta(n_fit,:) = [event,fit_theta];
                                        Dtheta(n_fit,:) = [event,fit_dtheta];
                                        Dtheta_division(n_fit,:) = [event,fit_dtheta_division];
                                        
                                        if theta<(theta_min+theta_bin_size)
                                            n_fit1 = n_fit1 + 1;
                                            Theta_Error1(n_fit1,:) = [event,theta_error];
                                            Phi_Error1(n_fit1,:) = [event,phi_error];
                                            Dtheta_Error1(n_fit1,:) = [event,dtheta_error];
                                        elseif theta<(theta_min+2*theta_bin_size)
                                            n_fit2 = n_fit2 + 1;
                                            Theta_Error2(n_fit2,:) = [event,theta_error];
                                            Phi_Error2(n_fit2,:) = [event,phi_error];
                                            Dtheta_Error2(n_fit2,:) = [event,dtheta_error];
                                        elseif theta<(theta_min+3*theta_bin_size)
                                            n_fit3 = n_fit3 + 1;
                                            Theta_Error3(n_fit3,:) = [event,theta_error];
                                            Phi_Error3(n_fit3,:) = [event,phi_error];
                                            Dtheta_Error3(n_fit3,:) = [event,dtheta_error];
                                        else
                                            'WTF'
                                            pause;
                                        end
                                        if theta_error<0.02 && phi_error<0.02
                                            n_fit_good = n_fit_good + 1;
                                            Event_Analysis(event,7) = 1;
                                        end
                                    else
                                        n_not_fit = n_not_fit + 1;
                                        if theta<(theta_min+theta_bin_size)
                                            n_not_fit1 = n_not_fit1 + 1;
                                        elseif theta<(theta_min+2*theta_bin_size)
                                            n_not_fit2 = n_not_fit2 + 1;
                                        else
                                            n_not_fit3 = n_not_fit3 + 1;
                                        end
                                    % end
                                    
                                    
                                    %True hits in fit
                                    if Event_Fit(event,6)==0
                                        nxt0 = nxt0 + 1;
                                    elseif Event_Fit(event,6)==1
                                        nxt1 = nxt1 + 1;
                                    elseif Event_Fit(event,6)==2
                                        nxt2 = nxt2 + 1;
                                    elseif Event_Fit(event,6)==3
                                        nxt3 = nxt3 + 1;
                                    elseif Event_Fit(event,6)==4
                                        nxt4 = nxt4 + 1;
                                    end
                                    if Event_Fit(event,7)==0
                                        nuvt0 = nuvt0 + 1;
                                    elseif Event_Fit(event,7)==1
                                        nuvt1 = nuvt1 + 1;
                                    elseif Event_Fit(event,7)==2
                                        nuvt2 = nuvt2 + 1;
                                    elseif Event_Fit(event,7)==3
                                        nuvt3 = nuvt3 + 1;
                                    elseif Event_Fit(event,7)==4
                                        nuvt4 = nuvt4 + 1;
                                    end
                                    
                                    %Background hits in fit
                                    if Event_Fit(event,8)==0
                                        nxb0 = nxb0 + 1;
                                    elseif Event_Fit(event,8)==1
                                        nxb1 = nxb1 + 1;
                                    elseif Event_Fit(event,8)==2
                                        nxb2 = nxb2 + 1;
                                    elseif Event_Fit(event,8)==3
                                        nxb3 = nxb3 + 1;
                                    elseif Event_Fit(event,8)==4
                                        nxb4 = nxb4 + 1;
                                    end
                                    if Event_Fit(event,9)==0
                                        nuvb0 = nuvb0 + 1;
                                    elseif Event_Fit(event,9)==1
                                        nuvb1 = nuvb1 + 1;
                                    elseif Event_Fit(event,9)==2
                                        nuvb2 = nuvb2 + 1;
                                    elseif Event_Fit(event,9)==3
                                        nuvb3 = nuvb3 + 1;
                                    elseif Event_Fit(event,9)==4
                                        nuvb4 = nuvb4 + 1;
                                    end
                                    
                                    if Event_Fit(event,8)~=0 || Event_Fit(event,9)~=0
                                        n_events_with_bg = n_events_with_bg + 1;
                                    end
                                    
                                    end
                                    
                                    true_track_hits = true_track_hits + Event_Fit(event,6) + Event_Fit(event,7);
                                    bg_track_hits = bg_track_hits + Event_Fit(event,8) + Event_Fit(event,9);
                                    
                                    
                                    drift_sum = drift_sum + Event_Info(event,24);
                                    slope_spread_X_sum = slope_spread_X_sum + Event_Info(event,25);
                                    slope_spread_UV_sum = slope_spread_UV_sum + Event_Info(event,26);
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    qual_events1 = n_fit1 + n_not_fit1;
    qual_events2 = n_fit2 + n_not_fit2;
    qual_events3 = n_fit3 + n_not_fit3;
    
    Tuning_Analysis(i,1:14) = [qual_events,n_fit,Get_Statistics(Theta_Error(:,2)),Get_Statistics(Phi_Error(:,2)),Get_Statistics(Dtheta_Error(:,2))];
    Tuning_Analysis(i,15:28) = [qual_events1,n_fit1,Get_Statistics(Theta_Error1(:,2)),Get_Statistics(Phi_Error1(:,2)),Get_Statistics(Dtheta_Error1(:,2))];
    Tuning_Analysis(i,29:42) = [qual_events2,n_fit2,Get_Statistics(Theta_Error2(:,2)),Get_Statistics(Phi_Error2(:,2)),Get_Statistics(Dtheta_Error2(:,2))];
    Tuning_Analysis(i,43:56) = [qual_events3,n_fit3,Get_Statistics(Theta_Error3(:,2)),Get_Statistics(Phi_Error3(:,2)),Get_Statistics(Dtheta_Error3(:,2))];
    
    Tuning_Analysis(i,57:66) = [nxt0,nxt1,nxt2,nxt3,nxt4,nuvt0,nuvt1,nuvt2,nuvt3,nuvt4];
    Tuning_Analysis(i,67:76) = [nxb0,nxb1,nxb2,nxb3,nxb4,nuvb0,nuvb1,nuvb2,nuvb3,nuvb4];
    
    Tuning_Analysis(i,77) = n_fit_good;
    Tuning_Analysis(i,78:80) = [n_events_with_bg,true_track_hits,bg_track_hits];
    
    Tuning_Analysis(i,81) = drift_sum / qual_events;
    Tuning_Analysis(i,82) = slope_spread_X_sum / qual_events;
    Tuning_Analysis(i,83) = slope_spread_UV_sum / qual_events;
    
end


n_fit

fit = n_fit/qual_events
good = n_fit_good/n_fit     %qual_events

qual_events

[n_fit1/qual_events1,n_fit2/qual_events2,n_fit3/qual_events3]

sum(Event_Analysis(:,6))/sum(Event_Analysis(:,5))

%--Debugging-----------
for i=1:1
debug=0;
if debug==1
    % figure(1); scatter(Wild(:,2),Phi_Error(:,2))
    % figure(2); scatter(Wild(:,1),Wild(:,2))
    %figure(3); scatter(Wild(:,3),Wild(:,2))
    figure(34); hold on; scatter(D(:,1),D(:,2));hold off;
    figure(35); hold on; hist(D(:,3));hold off;
    figure(36); hold on; hist(D(:,4));hold off;
%     figure(5); scatter(D(:,1),D(:,3))
    % figure(6); scatter(D(:,3),D(:,2))
end

debug=0;
if debug==1
    % scatter(Event_Info(:,6),Event_Info(:,5))
    
    figure(1)
    scatter(Event_Info(:,5),Event_Analysis(:,3));
    figure(2)
    scatter(Event_Info(:,5),Event_Analysis(:,2));
    
    figure(3)
    scatter(Event_Info(:,6),Event_Analysis(:,3));
    
    figure(4)
    scatter(Event_Info(:,6),Event_Analysis(:,2));
    
    figure(5)
    m=0;
    [a,b]=size(Event_Info);
    for q=1:a
        if Event_Info(q,15)==1
            m=m+1;
            D(m,:)=[Event_Info(q,6),Event_Info(q,5)];
        end
    end
    % scatter(Event_Info(:,6),Event_Info(:,5))
    scatter(D(:,1),D(:,2))
end
end
%-------------------------

clear Dtheta
n_fit=0;
debug2=1;
if debug2==1
    a=0;
    b=0;
    for i=1:event_mark
        if Event_Analysis(i,6)==1
            n_fit=n_fit+1;
            Dtheta(n_fit,:)=[i,Event_Fit(i,10),Event_Info(i,33)];
        end
    end
end
figure(1223)
scatter(Dtheta(:,3),Dtheta(:,2))
figure(1224)
hist(Dtheta(:,3))
figure(1225)
hist(Dtheta(:,2))



% close all
switch tunetype
    case 'run'
                 Get_Run_Plots(Tune,Tuning_Analysis,Theta_Error,Phi_Error,Dtheta_Error,n_fit)
        %         Get_Efficiency_Plot(Event_Analysis,5,1000,800)
        %         Get_Efficiency_Plot(Event_Analysis,3,100,900)
        
        
         Get_Run_Hist(Tune,Tuning_Analysis,Theta_Error,Phi_Error,Dtheta_Error,n_fit)
                Get_Data_Set_Properties_Without_Background(Event_Info)
%                  Background_Hit_Number_Influence(Event_Analysis,12)
                Event_Analysis(:,14)=Event_Analysis(:,10)-Event_Analysis(:,12);
                 Background_Hit_Number_Influence(Event_Analysis,14)
        
      %  Get_Efficiency_Plot(Event_Analysis,33,10,9000)
        %
        %         [mt,st] = normfit(Theta(:,2));
        %         [mp,sp] = normfit(Phi(:,2));
        %         [md,sd] = normfit(Dtheta(:,2));
        %         [mdd,sdd] = normfit(Dtheta_division(:,2));
        %         dtheta_rms = rms(Dtheta(:,2))
        %         dtheta_division_rms = rms(Dtheta_division(:,2))
        %         figure(65)
        %         histfit(Dtheta(:,2))
        %         figure(66)
        %         histfit(Dtheta_division(:,2))
        %         Dtheta_options(Dtheta, Dtheta_division)
        
    case 'h'
        Get_H_Plots(Tune,Tuning_Analysis)
    case 'uv'
        Get_UV_Plots(Tune,Tuning_Analysis)
    case 'charge'
        Get_Charge_Plots(Tune,Tuning_Analysis)
end


j=0;
for i=1:event_mark
    if Event_Info(i,9)~=0
        j=j+1;
    end
end

% j=0;
% for i=1:event_mark
%     if Event_Analysis(i,5)==1
%         j=j+1;
%         histo(j) = Event_Info(i,3)/1000;
%     end
% end
% figure(786)
% hist(histo,100)



j=0;
for i=1:event_mark
    if Event_Analysis(i,5)==1 && Event_Info(i,25)~=0  %&& Event_Info(i,25)<0.1 && Event_Analysis(i,6)==1
        j=j+1;
        histo(j) = Event_Info(i,25);
    end
end
figure(786)
hist(histo,0:0.000001:0.05)
j;
dis=mean(histo);








%---------------------


switch tunetype
    case 'run'
        save(savefile,'Tuning_Analysis');  %'Did not save!'
    otherwise
        save(savefile,'Tuning_Analysis');
end

end


function Dtheta_options(Dtheta, Dtheta_division)

f=3;
count=0;
count_d=0;

[n,a]=size(Dtheta);
[n_d,a]=size(Dtheta_division);

r = rms(Dtheta);
r_d = rms(Dtheta_division);
m = mean(Dtheta);
m_d = mean(Dtheta_division);



for i=1:n
    if Dtheta(i,2)> 0.0075 || Dtheta(i,2)< -0.0075
        count = count+1;
    end
    
    if Dtheta_division(i,2)>0.0075 || Dtheta_division(i,2)<-0.0075
        count_d = count_d+1;
    end
end
a= count/n
b=count_d/n




end


function Get_Efficiency_Plot(Event_Analysis,bin_parameter_index,n_bins,plot_index)

global Event_Info


[n_events,n]=size(Event_Analysis)


%/////////   Define histogram bins   ///////////////
bin_parameter_min = 100000;
bin_parameter_max = -9999;
for i=1:size(Event_Info)
    if Event_Analysis(i,5)==1
        if Event_Info(i,bin_parameter_index)<bin_parameter_min
            bin_parameter_min = Event_Info(i,bin_parameter_index);
            a=i;
        end
        if Event_Info(i,bin_parameter_index)>bin_parameter_max
            bin_parameter_max = Event_Info(i,bin_parameter_index);
            b=i;
        end
    end
end

bin_width = (bin_parameter_max-bin_parameter_min)/n_bins;
Hist_Table=zeros(n_bins,6);

for i=1:n_bins
    Hist_Table(i,1)=(bin_parameter_min+bin_width/2)+(i-1)*bin_width;
end

%//////////  Bin the events //////////////
for i=1:n_events
    if Event_Analysis(i,5)==1
        parameter=Event_Info(i,bin_parameter_index);
        for j=1:n_bins
            if parameter<(Hist_Table(j,1) + bin_width/2)
                Hist_Table(j,2)=Hist_Table(j,2)+1;
                if Event_Analysis(i,6)==1
                    Hist_Table(j,3)=Hist_Table(j,3)+1;
                    Hist_Table(j,5)=Hist_Table(j,5)+Event_Analysis(i,3);
                end
                if Event_Analysis(i,7)==1
                    Hist_Table(j,4)=Hist_Table(j,4)+1;
                end
                break;
            end
        end
    end
end


figure(plot_index);
hold on;
scatter(Hist_Table(:,1),Hist_Table(:,3)./Hist_Table(:,2));  %fit?
scatter(Hist_Table(:,1),Hist_Table(:,4)./Hist_Table(:,2));  %good?
scatter(Hist_Table(:,1),Hist_Table(:,5)./Hist_Table(:,3));   %avg error
hold off;


end

function Get_Charge_Plots(Tune,Tuning_Analysis)
global z_large strip_width

figure(1);
hold on;
scatter(Tune,Tuning_Analysis(:,13))
% scatter(Tune,Tuning_Analysis(:,81)/100)
scatter(Tune,Tuning_Analysis(:,82))
scatter(Tune,Tuning_Analysis(:,83))
hold off;

figure(2);
hold on;
scatter(Tune,Tuning_Analysis(:,77)./Tuning_Analysis(:,1))
% scatter(Tune,Tuning_Analysis(:,81)/100)
scatter(Tune,Tuning_Analysis(:,82))
scatter(Tune,Tuning_Analysis(:,83))
hold off;


figure(3);
hold on;
scatter(Tune,Tuning_Analysis(:,2))
hold off;



end

%
% function Get_UV_Plots(Tune,Tuning_Analysis)  %proportions
% global z_large strip_width
%
% Plot_Tuning_uv(Tune,Tuning_Analysis(:,66)./Tuning_Analysis(:,1),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65))./Tuning_Analysis(:,1),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65)+Tuning_Analysis(:,64))./Tuning_Analysis(:,1));
% hold on;
% scatter(Tune,Tuning_Analysis(:,13)*100,'x')
% hold off;
%
% Plot_Tuning_uv(Tune*z_large(1)/strip_width,Tuning_Analysis(:,66)./Tuning_Analysis(:,1),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65))./Tuning_Analysis(:,1),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65)+Tuning_Analysis(:,64))./Tuning_Analysis(:,1));
% hold on;
% scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,13)*100,'x')
% hold off;
%
% figure(11);
% hold on;
% scatter(Tune,Tuning_Analysis(:,77)./Tuning_Analysis(:,1))
% scatter(Tune,Tuning_Analysis(:,13)*100)
% hold off;
%
%
% Plot_Tuning_uv(Tune,Tuning_Analysis(:,16)./Tuning_Analysis(:,15),Tuning_Analysis(:,30)./Tuning_Analysis(:,29),Tuning_Analysis(:,44)./Tuning_Analysis(:,43));
% hold on;
% scatter(Tune,Tuning_Analysis(:,12)*100,'x')
% hold off;
%
%
% end
%
% function Get_H_Plots(Tune,Tuning_Analysis) %proportions
% global z_large strip_width
%
% Plot_Tuning_h(Tune,Tuning_Analysis(:,61)./Tuning_Analysis(:,1),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60))./Tuning_Analysis(:,1),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60)+Tuning_Analysis(:,59))./Tuning_Analysis(:,1));
% hold on;
% scatter(Tune,Tuning_Analysis(:,13)*100,'x')
% hold off;
%
% Plot_Tuning_h(Tune*z_large(1)/strip_width,Tuning_Analysis(:,61)./Tuning_Analysis(:,1),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60))./Tuning_Analysis(:,1),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60)+Tuning_Analysis(:,59))./Tuning_Analysis(:,1));
% hold on;
% scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,13)*100,'x')
% hold off;
%
% figure(11);
% hold on;
% scatter(Tune,Tuning_Analysis(:,77)./Tuning_Analysis(:,1))
% scatter(Tune,Tuning_Analysis(:,12)*100)
% hold off;
%
% Plot_Tuning_h(Tune,Tuning_Analysis(:,16)./Tuning_Analysis(:,15),Tuning_Analysis(:,30)./Tuning_Analysis(:,29),Tuning_Analysis(:,44)./Tuning_Analysis(:,43));
% hold on;
% scatter(Tune,Tuning_Analysis(:,12)*100,'x')
% hold off;
%
% end



function Get_UV_Plots(Tune,Tuning_Analysis)  %no proportions
global z_large strip_width

Plot_Tuning_uv(Tune,Tuning_Analysis(:,66),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65)),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65)+Tuning_Analysis(:,64)));
hold on;
scatter(Tune,Tuning_Analysis(:,13)*100,'x')
hold off;
 
% Plot_Tuning_uv(Tune,Tuning_Analysis(:,66),(Tuning_Analysis(:,65)),(Tuning_Analysis(:,64)));
% hold on;
% scatter(Tune,Tuning_Analysis(:,13)*100,'x')
% hold off;

Plot_Tuning_uv(Tune*z_large(1)/strip_width,Tuning_Analysis(:,66),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65)),(Tuning_Analysis(:,66)+Tuning_Analysis(:,65)+Tuning_Analysis(:,64)));
hold on;
scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,13)*100,'x')
hold off;

figure(11);
hold on;
scatter(Tune,Tuning_Analysis(:,77))
scatter(Tune,Tuning_Analysis(:,13)*100)
hold off;

figure(12);
hold on;
scatter(Tune,Tuning_Analysis(:,2)./Tuning_Analysis(:,1))
scatter(Tune,Tuning_Analysis(:,16)./Tuning_Analysis(:,15))
scatter(Tune,Tuning_Analysis(:,30)./Tuning_Analysis(:,29))
scatter(Tune,Tuning_Analysis(:,44)./Tuning_Analysis(:,43))
hold off;

figure(111)
hold on;
scatter(Tune,Tuning_Analysis(:,2))
scatter(Tune,Tuning_Analysis(:,78))
hold off;

Plot_Tuning_uv(Tune,Tuning_Analysis(:,16),Tuning_Analysis(:,30),Tuning_Analysis(:,44));
hold on;
scatter(Tune,Tuning_Analysis(:,12)*100,'x')
hold off;


end

function Get_H_Plots(Tune,Tuning_Analysis) %no proportions
global z_large strip_width

Plot_Tuning_h(Tune,Tuning_Analysis(:,61),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60)),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60)+Tuning_Analysis(:,59)));
hold on;
scatter(Tune,Tuning_Analysis(:,13)*100,'x')
hold off;
 
% Plot_Tuning_h(Tune,Tuning_Analysis(:,61),(Tuning_Analysis(:,61)),(Tuning_Analysis(:,59)));
% hold on;
% scatter(Tune,Tuning_Analysis(:,13)*100,'x')
% hold off;

% Plot_Tuning_h(Tune*z_large(1)/strip_width,Tuning_Analysis(:,61),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60)),(Tuning_Analysis(:,61)+Tuning_Analysis(:,60)+Tuning_Analysis(:,59)));
% hold on;
% scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,13)*100,'x')
% hold off;
% 
% figure(11);
% hold on;
% scatter(Tune,Tuning_Analysis(:,77)./Tuning_Analysis(:,1))
% scatter(Tune,Tuning_Analysis(:,12)*100)
% hold off;
% 
% 
figure(12);
hold on;
scatter(Tune,Tuning_Analysis(:,2)./Tuning_Analysis(:,1))
scatter(Tune,Tuning_Analysis(:,16)./Tuning_Analysis(:,15))
scatter(Tune,Tuning_Analysis(:,30)./Tuning_Analysis(:,29))
scatter(Tune,Tuning_Analysis(:,44)./Tuning_Analysis(:,43))
hold off;
% 
% 
% figure(112)
% hold on;
% scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,2))
% scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,78))
% hold off;
% figure(111)
% hold on;
% scatter(Tune,Tuning_Analysis(:,2))
% scatter(Tune,Tuning_Analysis(:,78))
% hold off;
% 
% 
% 
% Plot_Tuning_h(Tune,Tuning_Analysis(:,16),Tuning_Analysis(:,30),Tuning_Analysis(:,44));
% hold on;
% scatter(Tune,Tuning_Analysis(:,12)*100,'x')
% hold off;

figure(111)
hold on;
scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,2))
scatter(Tune*z_large(1)/strip_width,Tuning_Analysis(:,78))
hold off;

figure(121)
hold on
scatter(Tune,Tuning_Analysis(:,2)./Tuning_Analysis(:,1))
scatter(Tune,Tuning_Analysis(:,12)*100,'x')
hold off

end



function Get_Run_Hist(Tune,Tuning_Analysis,Theta_Error,Phi_Error,Dtheta_Error,n_fit)


%------------------------
Stats = Get_Statistics(Theta_Error(:,2));
mean=Stats(1);
sig=Stats(2);
rms=Stats(3);
tails=Stats(4);
figure(2001)
hold on;
hist(Theta_Error(:,2),200,-0.04,0.0)
% scatter([-0.04+0.0001:0.0002:-0.00-0.0001],h)
title('Theta Error','FontSize',18);
xlabel('Theta Error (mrad)','FontSize',16);
ylabel('Number of Events','FontSize',16);
text_mean = sprintf('mean = %irad',mean);
text_sig = sprintf('sig = %irad',sig);
text_rms = sprintf('rms = %irad',rms);
text_tails = sprintf('tails = %i%',tails);
text_events = sprintf('events fitted = %i',n_fit);
annotation('textbox','String',{text_mean,text_sig,text_rms,text_tails,text_events},'FitBoxToText','on','FontSize',16);

%------------------------
Stats = Get_Statistics(Phi_Error(:,2));
mean=Stats(1);
sig=Stats(2);
rms=Stats(3);
tails=Stats(4);
figure(2002)
hist(Phi_Error(:,2),1000)
title('Phi Error','FontSize',18);
xlabel('Phi Error (mrad)','FontSize',16);
ylabel('Number of Events','FontSize',16);
text_mean = sprintf('mean = %irad',mean);
text_sig = sprintf('sig = %irad',sig);
text_rms = sprintf('rms = %irad',rms);
text_tails = sprintf('tails = %i%',tails);
text_events = sprintf('events fitted = %i',n_fit);
annotation('textbox','String',{text_mean,text_sig,text_rms,text_tails,text_events},'FitBoxToText','on','FontSize',16);


%------------------------
Stats = Get_Statistics(Dtheta_Error(:,2));
mean=Stats(1);
sig=Stats(2);
rms=Stats(3);
tails=Stats(4);
figure(2003)
hist(Dtheta_Error(:,2),1000)
title('Delta Theta Error','FontSize',18);
xlabel('Delta Theta Error (mrad)','FontSize',16);
ylabel('Number of Events','FontSize',16);
text_mean = sprintf('mean = %irad',mean);
text_sig = sprintf('sig = %irad',sig);
text_rms = sprintf('rms = %irad',rms);
text_tails = sprintf('tails = %i%',tails);
text_events = sprintf('events fitted = %i',n_fit);
annotation('textbox','String',{text_mean,text_sig,text_rms,text_tails,text_events},'FitBoxToText','on','FontSize',16);




% figure(8745)
% hold on
% T = Get_Bin_Data(Dtheta_Error,33,30)
% scatter(T(:,1),T(:,2)./500./1000)
% scatter(T(:,1),T(:,3))
% scatter(T(:,1),T(:,4))
% hold off
% % global Event_Info
% % C = Event_Info(Dtheta_Error(:,1),33);
% % scatter(C,Dtheta_Error(:,2))

end

function Get_Run_Plots(Tune,Tuning_Analysis,Theta_Error,Phi_Error,Dtheta_Error,n_fit)

global Data_Trimmed

% figure(1005);
% Get_Statistics(Dtheta_Error(:,2));
% % Dtheta_Error_Trimmed = Move_Trimmed_Data(0,'import');
% histfit(Data_Trimmed)



%----------
bin_parameter_index = 5;
n_bins = 100;



%     figure(1000);
%     hold 'on';
%     errorbar(Tune,[Tuning_Analysis(1,19),Tuning_Analysis(1,33),Tuning_Analysis(1,47)],...
%         [Tuning_Analysis(1,19),Tuning_Analysis(1,33),Tuning_Analysis(1,47)]./...
%         sqrt([Tuning_Analysis(1,16),Tuning_Analysis(1,30),Tuning_Analysis(1,44)]),'b')
%     errorbar(Tune,[Tuning_Analysis(1,23),Tuning_Analysis(1,37),Tuning_Analysis(1,51)],...
%         [Tuning_Analysis(1,27),Tuning_Analysis(1,41),Tuning_Analysis(1,55)]./...
%         sqrt([Tuning_Analysis(1,16),Tuning_Analysis(1,30),Tuning_Analysis(1,44)]),'g')
%     errorbar(Tune,[Tuning_Analysis(1,27),Tuning_Analysis(1,41),Tuning_Analysis(1,55)],...
%         [Tuning_Analysis(1,27),Tuning_Analysis(1,41),Tuning_Analysis(1,55)]./...
%         sqrt([Tuning_Analysis(1,16),Tuning_Analysis(1,30),Tuning_Analysis(1,44)]),'r')
%     hold 'off';
%
%
%
%     figure(1001);
%     hold 'on';
%     errorbar(Tune,[Tuning_Analysis(1,17),Tuning_Analysis(1,31),Tuning_Analysis(1,45)],...
%         [Tuning_Analysis(1,18),Tuning_Analysis(1,32),Tuning_Analysis(1,46)],'b')
%     errorbar(Tune,[Tuning_Analysis(1,21),Tuning_Analysis(1,35),Tuning_Analysis(1,49)],...
%         [Tuning_Analysis(1,22),Tuning_Analysis(1,36),Tuning_Analysis(1,50)],'g')
%     errorbar(Tune,[Tuning_Analysis(1,25),Tuning_Analysis(1,39),Tuning_Analysis(1,53)],...
%         [Tuning_Analysis(1,26),Tuning_Analysis(1,40),Tuning_Analysis(1,54)],'r')
%     hold 'off';
%
%


base_index = 0;
bin_parameter_index = 5;
n_bins = 10;

% %--------------
% figure(base_index+3);
% hold on;
% title('Tails Removed -- Mean (with STD) v. Theta');
% xlabel('Theta (radians)');
% ylabel('Mean (radians)');
%
% Hist_Table = Get_Bin_Data(Theta_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,3),Hist_Table(:,4),'b');
% tails_theta = Hist_Table(:,6);
%
% Hist_Table = Get_Bin_Data(Phi_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,3),Hist_Table(:,4),'g');
% tails_phi = Hist_Table(:,6);
%
% Hist_Table = Get_Bin_Data(Dtheta_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,3),Hist_Table(:,4),'r');
% tails_dtheta = Hist_Table(:,6);
%
% legend('Theta','Phi','DTheta');
%
% annotation('textbox','String',{'proportion in theta tails',tails_theta},'FitBoxToText','on');
% annotation('textbox','String',{'proportion in phi tails',tails_phi},'FitBoxToText','on');
% annotation('textbox','String',{'proportion in dtheta tails',tails_dtheta},'FitBoxToText','on');
%
% hold off;
% %-------------------

%-------------
figure(base_index+12);
hold on;
title('Tails Removed -- Error STD v. Theta');  % -- no averaging up front');
xlabel('Theta (radians)');
ylabel('Error Standard Deviation (mrad)');

Hist_Table = Get_Bin_Data(Theta_Error,bin_parameter_index,n_bins);
errorbar(Hist_Table(:,1),Hist_Table(:,4).*1000,Hist_Table(:,4).*1000./sqrt(Hist_Table(:,2)),'b');
tails_theta = Hist_Table(:,6);

Hist_Table = Get_Bin_Data(Phi_Error,bin_parameter_index,n_bins);
errorbar(Hist_Table(:,1),Hist_Table(:,4).*1000,Hist_Table(:,4).*1000./sqrt(Hist_Table(:,2)),'g');
tails_phi = Hist_Table(:,6);

Hist_Table = Get_Bin_Data(Dtheta_Error,bin_parameter_index,n_bins);
errorbar(Hist_Table(:,1),Hist_Table(:,4).*1000,Hist_Table(:,4).*1000./sqrt(Hist_Table(:,2)),'r');
tails_dtheta = Hist_Table(:,6);

legend('Theta','Phi','DTheta');

annotation('textbox','String',{'proportion in theta tails',tails_theta},'FitBoxToText','on');
annotation('textbox','String',{'proportion in phi tails',tails_phi},'FitBoxToText','on');
annotation('textbox','String',{'proportion in dtheta tails',tails_dtheta},'FitBoxToText','on');

% annotation('textbox','String',{'proportion in tails',tails,sprintf('core %i rms fitted',core)},'FitBoxToText','on');
hold off;
%---------------


% %-------------
% figure(base_index+13);
% hold on;
% title('Tails Removed -- Error STD v. Momentum');  % -- no averaging up front');
% xlabel('Theta (radians)');
% ylabel('Error RMS (radians)');
%
% bin_parameter_index = 3;
%
% Hist_Table = Get_Bin_Data(Theta_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,4),Hist_Table(:,4)./sqrt(Hist_Table(:,2)),'b');
% tails_theta = Hist_Table(:,6);
%
% Hist_Table = Get_Bin_Data(Phi_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,4),Hist_Table(:,4)./sqrt(Hist_Table(:,2)),'g');
% tails_phi = Hist_Table(:,6);
%
% Hist_Table = Get_Bin_Data(Dtheta_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,4),Hist_Table(:,4)./sqrt(Hist_Table(:,2)),'r');
% tails_dtheta = Hist_Table(:,6);
%
% legend('Theta','Phi','DTheta');
%
% annotation('textbox','String',{'proportion in theta tails',tails_theta},'FitBoxToText','on');
% annotation('textbox','String',{'proportion in phi tails',tails_phi},'FitBoxToText','on');
% annotation('textbox','String',{'proportion in dtheta tails',tails_dtheta},'FitBoxToText','on');
%
% % annotation('textbox','String',{'proportion in tails',tails,sprintf('core %i rms fitted',core)},'FitBoxToText','on');
% hold off;
% %---------------

% %-------------
% figure(base_index+2);
% hold on;
% title('Tails Removed -- Error RMS v. Theta');  % -- no averaging up front');
% xlabel('Theta (radians)');
% ylabel('Error RMS (radians)');
%
% bin_parameter_index = 5;
%
% Hist_Table = Get_Bin_Data(Theta_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,5),Hist_Table(:,5)./sqrt(Hist_Table(:,2)),'b');
% tails_theta = Hist_Table(:,6);
%
% Hist_Table = Get_Bin_Data(Phi_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,5),Hist_Table(:,5)./sqrt(Hist_Table(:,2)),'g');
% tails_phi = Hist_Table(:,6);
%
% Hist_Table = Get_Bin_Data(Dtheta_Error,bin_parameter_index,n_bins);
% errorbar(Hist_Table(:,1),Hist_Table(:,5),Hist_Table(:,5)./sqrt(Hist_Table(:,2)),'r');
% tails_dtheta = Hist_Table(:,6);
%
% legend('Theta','Phi','DTheta');
%
% annotation('textbox','String',{'proportion in theta tails',tails_theta},'FitBoxToText','on');
% annotation('textbox','String',{'proportion in phi tails',tails_phi},'FitBoxToText','on');
% annotation('textbox','String',{'proportion in dtheta tails',tails_dtheta},'FitBoxToText','on');
%
% % annotation('textbox','String',{'proportion in tails',tails,sprintf('core %i rms fitted',core)},'FitBoxToText','on');
% hold off;
% %---------------


% if n_bins==1
% figure(base_index+3);
% hold 'on'
% errorbar((theta_max-theta_min)/2,Tuning_Analysis(1,3),Tuning_Analysis(1,4))
% errorbar((theta_max-theta_min)/2,Tuning_Analysis(1,7),Tuning_Analysis(1,8))
% errorbar((theta_max-theta_min)/2,Tuning_Analysis(1,11),Tuning_Analysis(1,12))
% hold 'off'
% figure(base_index+2);
% hold 'on'
% errorbar((theta_max-theta_min)/2,Tuning_Analysis(1,5),Tuning_Analysis(1,5)/sqrt(Tuning_Analysis(1,2)))
% errorbar((theta_max-theta_min)/2,Tuning_Analysis(1,9),Tuning_Analysis(1,9)/sqrt(Tuning_Analysis(1,2)))
% errorbar((theta_max-theta_min)/2,Tuning_Analysis(1,13),Tuning_Analysis(1,13)/sqrt(Tuning_Analysis(1,2)))
% hold 'off'
% end

end

function Get_Data_Set_Properties_Without_Background(Event_Info)
global event_mark
%Hits per event
n=200;
N_Hits_per_Event = zeros(n,3);
for i=1:event_mark
    n_hits = Event_Info(i,10);
    if n_hits>0
        N_Hits_per_Event(n_hits,2) = N_Hits_per_Event(n_hits,2) + 1;
    end
end
for i=1:n
    N_Hits_per_Event(i,1)=i;
    N_Hits_per_Event(i,3)=sum(N_Hits_per_Event(1:i,2));
end
figure(123)
scatter(N_Hits_per_Event(:,1),N_Hits_per_Event(:,3)/sum(N_Hits_per_Event(:,2)))

figure(124)
bar(N_Hits_per_Event(:,1),N_Hits_per_Event(:,2))



%IP spread
n=0;
for event = 1:event_mark
    if Event_Info(event,15)==1
        n=n+1;
        IPxyz(n,:) = Event_Info(event,31:33);
    end
end
figure(342)
hold on;
histfit(IPxyz(:,1)), histfit(IPxyz(:,2))
hold off;
figure(343)
histfit(IPxyz(:,3))

[mu_x,sig_x]=normfit(IPxyz(:,1))
[mu_y,sig_y]=normfit(IPxyz(:,2))
[mu_z,sig_z]=normfit(IPxyz(:,3))

    

% figure(125)
%insert drift time plot

end

function output = Get_Bin_Data(B,bin_parameter_index,n_bins) %,core,base_index)

global Event_Info

%B must have one entry per event

[n_events,mm]=size(B);
B = [B,zeros(n_events,1)];  %grow B by 1 column for indexing histogram

%/////////   Define histogram bins   ///////////////
bin_parameter_min = 100000;
bin_parameter_max = -9999;
for i=1:size(Event_Info)
    if Event_Info(i,15)==1
        if Event_Info(i,bin_parameter_index)<bin_parameter_min
            bin_parameter_min = Event_Info(i,bin_parameter_index);
            a=i;
        end
        if Event_Info(i,bin_parameter_index)>bin_parameter_max
            bin_parameter_max = Event_Info(i,bin_parameter_index);
            b=i;
        end
    end
end


bin_width = (bin_parameter_max-bin_parameter_min)/n_bins;
Hist_Table=zeros(n_bins,6);

for i=1:n_bins
    Hist_Table(i,1)=(bin_parameter_min+bin_width/2)+(i-1)*bin_width;
end

%//////////  Bin the events //////////////
for i=1:n_events
    parameter=Event_Info(B(i,1),bin_parameter_index);
    for j=1:n_bins
        if parameter<(Hist_Table(j,1) + bin_width/2)
            Hist_Table(j,2)=Hist_Table(j,2)+1;
            B(i,mm+1) = j;
            break;
        end
    end
end

%/////////  Calc statistics  //////////////
for j=1:n_bins
    clear Points_for_Fit
    a=0;
    for i=1:n_events
        if B(i,mm+1) == j;
            a=a+1;
            Points_for_Fit(a) = B(i,2);
        end
    end
    if a==0
        Points_for_Fit(1) = 0;
    end
    
    Hist_Table(j,3:6) = Get_Statistics(Points_for_Fit);
    
end

output = Hist_Table;

return;

%//////   CREATE FIGURES    ////////////

figure(base_index);
hold on;
title('Error RMS v. Theta');  % -- no averaging up front');
xlabel('Theta (radians)');
ylabel('Error RMS (radians)');
errorbar(Hist_Table(:,1),Hist_Table(:,3),Hist_Table(:,3)./sqrt(Hist_Table(:,2)),'b');
errorbar(Hist_Table(:,1),Hist_Table(:,4),Hist_Table(:,4)./sqrt(Hist_Table(:,2)),'g');
errorbar(Hist_Table(:,1),Hist_Table(:,5),Hist_Table(:,5)./sqrt(Hist_Table(:,2)),'r');
legend('Theta','Phi','DTheta');
hold off;

figure(base_index+1);
hold on;
title('Mean (with STD) v. Theta');  % -- no averaging up front');
xlabel('Theta (radians)');
ylabel('Mean (radians)');
errorbar(Hist_Table(:,1),Hist_Table(:,6),Hist_Table(:,7),'b');
errorbar(Hist_Table(:,1),Hist_Table(:,8),Hist_Table(:,9),'g');
errorbar(Hist_Table(:,1),Hist_Table(:,10),Hist_Table(:,11),'r');
legend('Theta','Phi','DTheta');
hold off;


tails=(Hist_Table(:,2)-Hist_Table(:,12))./Hist_Table(:,2);

figure(base_index+2);
hold on;
title('Tails Removed -- Error RMS v. Theta');  % -- no averaging up front');
xlabel('Theta (radians)');
ylabel('Error RMS (radians)');
errorbar(Hist_Table(:,1),Hist_Table(:,13),Hist_Table(:,13)./sqrt(Hist_Table(:,12)),'b');
errorbar(Hist_Table(:,1),Hist_Table(:,14),Hist_Table(:,14)./sqrt(Hist_Table(:,12)),'g');
errorbar(Hist_Table(:,1),Hist_Table(:,15),Hist_Table(:,15)./sqrt(Hist_Table(:,12)),'r');
legend('Theta','Phi','DTheta');
annotation('textbox','String',{'proportion in tails',tails,sprintf('core %i rms fitted',core)},'FitBoxToText','on');
hold off;

figure(base_index+3);
hold on;
title('Tails Removed -- Mean (with STD) v. Theta');  % -- no averaging up front');
xlabel('Theta (radians)');
ylabel('Mean (radians)');
errorbar(Hist_Table(:,1),Hist_Table(:,16),Hist_Table(:,17),'b');
errorbar(Hist_Table(:,1),Hist_Table(:,18),Hist_Table(:,19),'g');
errorbar(Hist_Table(:,1),Hist_Table(:,20),Hist_Table(:,21),'r');
legend('Theta','Phi','DTheta');
annotation('textbox','String',{tails},'FitBoxToText','on');
hold off;




% 
% %//////   SAVE FIGURES    ////////////
% fileout=sprintf('/Users/clark/Desktop/Figures/%s_%i',base_index);
% saveas(figure(base_index),fileout);
% 
% fileout=sprintf('/Users/clark/Desktop/Figures/%s_%i',base_index+1);
% saveas(figure(base_index+1),fileout);
% 
% fileout=sprintf('/Users/clark/Desktop/Figures/%s_%i',base_index+2);
% saveas(figure(base_index+2),fileout);
% 
% fileout=sprintf('/Users/clark/Desktop/Figures/%s_%i',base_index+3);
% saveas(figure(base_index+3),fileout);

end



function output = Get_Statistics(Data0,returndata)

global Data_Trimmed

Data_Trimmed = Data0;

core = 3;
type = 'sig';
repeat = 1;

Data = Data0;

if length(Data)<100
    'too few events!';
    output=[0,0,0,0];
    return;
end

[mean,sig] = normfit(Data);
% rms = sqrt(sum(Data-mean)/length(Data));
rms = sqrt(mean^2+sig^2);
% rms = sqrt(sum(Data.^2)/length(Data));
percent_tails = 0;


if core>0
    for i=0:repeat
        if i>1
            Data = Data1;
            clear Data1=0;
        end
        switch type
            case 'sig'
                lower_limit = mean - core*sig;
                upper_limit = mean + core*sig;
            case 'rms'
                lower_limit = -core*rms;
                upper_limit = core*rms;
        end
        k=0;
        for j=1:length(Data)
            if Data(j)>lower_limit && Data(j)<upper_limit
                k=k+1;
                Data1(k) = Data(j);
            end
        end
        [mean,sig] = normfit(Data1);
        %         rms = sqrt(sum(Data-mean)/length(Data));
        rms = sqrt(mean^2+sig^2);
        %         rms = sqrt(sum(Data1.^2)/length(Data1));
    end
    
    percent_tails = (length(Data0)-length(Data1))/length(Data0);
    
    %     Move_Trimmed_Data(Data1,'export');
    
    Data_Trimmed = Data1;
    
end

output = [mean,sig,rms,percent_tails];
end

% function output = Move_Trimmed_Data(Data1,direction)
%     switch direction
%         case 'export'
%             output = 1;
%         case 'import'
%             output = Data1;
%     end
% end

function output = add_dtheta_event_info(A)  %dtheta truth is calculated here!!!
global intrinsic
if intrinsic==1
    A(:,7) = A(:,20) - A(:,19);
else
    A(:,7) = A(:,20) - A(:,5);
end

output = A;

end

function output=Do_altered_cuts(A,level)

CT_x = 0;
CT_uv = 0;
theta_min = 0.1722;
theta_max = 0.5254;
phi_min = 1.3028;
phi_max = 1.8388;
min_hits = 1;
max_hits = 10000;


[n,m]=size(A);
if level==1
    for i=1:n
        if A(i,9)==0 %|| A(i,5)<theta_min || A(i,5)>theta_max %|| A(i,6)==0
            A(i,15)=0;
        end
        %     if A(i,11)<=0 || A(i,12)<=0 || A(i,9)<=1
        %         A(i,15)=0;
        %     end
    end
end

if level==2
    for i=1:n
        if A(i,6)<phi_min || A(i,6)>phi_max
            A(i,15)=0;
        end
    end
end

if level==3
    for i=1:n
        if A(i,21)<phi_min || A(i,21)>phi_max
            A(i,15)=0;
        end
    end
end

if level==4
    for i=1:n
%         if A(i,5)<theta_min || A(i,5)>theta_max
%             A(i,15)=0;
%         end
        if A(i,19)<theta_min || A(i,19)>theta_max
            A(i,15)=0;
        end
    end
end

output = A;

end


% function output = Fix_Event_Info_Phi(A,index)
% [n,m]=size(A);
% for i=1:n
%     athena_phi = A(i,index);
%     if A(i,15)==1
%         if athena_phi>0
%             true_phi = pi/2-(athena_phi-pi);
%         else %if athena_phi<0
%             true_phi = pi/2-(pi-abs(athena_phi));
%         end
%         A(i,index)=true_phi;
%     end
% end
% output = A;
% end


function output = Fix_Event_Info_Phi(A,index)
[n,m]=size(A);
for i=1:n
    athena_phi = A(i,index);
    if A(i,15)==1
        if athena_phi>0
            true_phi = pi/2+(pi-athena_phi);
        elseif athena_phi<0
            true_phi = pi/2-(pi+athena_phi);
        end
        A(i,index)=true_phi;
    end
end
output = A;
end

function output = Event_Info_Alterations(charge,background,openfile,event_mark)
global strip_width z_large

energy=200;
stereodeg=15;

% openfolder = sprintf('/Users/%s/Desktop/eventfiles/',computer);
% openfile = sprintf('%s%iGeV_%ideg_%ie.mat',openfolder,energy,stereodeg,charge);
load(openfile);  %,'Event_Info_Without_Background','Event_Info_With_Background','Hits_Data_Set_Time_Without_Background','Hits_Data_Set_Time_With_Background')

switch background
    case 'bgoff'
        Event_Info = Event_Info_Without_Background;
    case 'bgon'
        Event_Info = Event_Info_With_Background;
end

Event_Info = Event_Info(1:event_mark,:);
event_mark

Event_Info = Do_altered_cuts(Event_Info,1);
Event_Info = Fix_Event_Info_Phi(Event_Info,6);
% Event_Info = Do_altered_cuts(Event_Info,2);
Event_Info = Fix_Event_Info_Phi(Event_Info,21);
Event_Info = Fix_Event_Info_Phi(Event_Info,22);
% Event_Info = Do_altered_cuts(Event_Info,3);
Event_Info = add_dtheta_event_info(Event_Info);

[a,b]=size(Event_Info)
if a<event_mark
    c=event_mark-a;
    Event_Info = [Event_Info;zeros(c,b)];
end


%-----------------

% for q=1:a
%     if Event_Info(q,21)>2
%         Event_Info(q,21)=Event_Info(q,21)-pi/2;
%     end
% end

% m=0;
% D=[0,0];
% [a,b]=size(Event_Info);
% for q=1:a
%     if Event_Info(q,15)==1
%         m=m+1;
%         D(m,:)=[Event_Info(q,21)-Event_Info(q,22),Event_Info(q,5)];
%     end
% end
% % scatter(Event_Info(:,6),Event_Info(:,5))
% scatter(D(:,1),D(:,2))


%-----------------------


% switch plane_type
%     case 'x'
%         vertical_strip_width = strip_width;
%     case 'uv'
%         vertical_strip_width = vertical_strip_width_UV;
%     otherwise
%         'NOT A CHOICE'
%         pause;
% end

A=Hits_Data_Set_Time_Without_Background;
[m,n]=size(A);
mark=2;
events_max = event_mark;  %[events_max,dummy] = size(Event_Info_Without_Background);
B=zeros(event_mark,1);
C=zeros(event_mark,2);
for event=1:events_max
    N_hits=0;
    drift_time_sum=0;
    j=0;
    k=0;
    Hit_Slopes_X = 0;
    Hit_Slopes_UV = 0;
    for i=mark:m
        if A(i,1) == event
            N_hits = N_hits + 1;
            drift_time_sum = drift_time_sum + (A(i,14)-100*event);   %conversion of drift time from spacing
            if A(i,5)==1 || A(i,5)==2 || A(i,5)==5 || A(i,5)==6
                j=j+1;
                Hit_Slopes_X(j) = A(i,6)*strip_width/z_large(A(i,5));
            elseif A(i,5)==3 || A(i,5)==4 || A(i,5)==7 || A(i,5)==8
                k=k+1;
                Hit_Slopes_UV(k) = A(i,6)*strip_width/z_large(A(i,5));
            end
            
%             Event_Info(event,6)=A(i,8);
            
        elseif A(i,1) > event
            mark = i;
            break;
        end
    end
    if N_hits>0
        drift_time_avg = drift_time_sum / N_hits;
        B(event) = drift_time_avg;
        if j>1
            %C(event,1) = mean(Hit_Slopes_X);
            C(event,1) = max(Hit_Slopes_X)-min(Hit_Slopes_X);
        end
        if k>1
            %C(event,2) = mean(Hit_Slopes_UV);
            C(event,2) = max(Hit_Slopes_UV)-min(Hit_Slopes_UV);
        end
    end
end

[bbb,bbbb]=size(B)
Event_Info(:,24) = B;
Event_Info(:,25:26) = C;

Event_Info = Do_altered_cuts(Event_Info,2);
Event_Info = Do_altered_cuts(Event_Info,3);
Event_Info = Do_altered_cuts(Event_Info,4);


output = Event_Info;

end


function Background_Hit_Number_Influence(Event_Analysis,hit_type_index)
global event_mark

figure(456)
hold on;
for i=1:3
    error_index=1+i;
    n_events = event_mark;
    hits0=0;
    hits1=0;
    hits2=0;
    hits3=0;
    hits4=0;
    Hits0=0;
    Hits1=0;
    Hits2=0;
    Hits3=0;
    Hits4=0;
    for event=1:n_events
        if Event_Analysis(event,6)==1
            if Event_Analysis(event,hit_type_index)==0
                hits0=hits0+1;
                Hits0(hits0)=Event_Analysis(event,error_index);
            elseif Event_Analysis(event,hit_type_index)==1
                hits1=hits1+1;
                Hits1(hits1)=Event_Analysis(event,error_index);
            elseif Event_Analysis(event,hit_type_index)==2
                hits2=hits2+1;
                Hits2(hits2)=Event_Analysis(event,error_index);
            elseif Event_Analysis(event,hit_type_index)==3
                hits3=hits3+1;
                Hits3(hits3)=Event_Analysis(event,error_index);
            elseif Event_Analysis(event,hit_type_index)==4
                hits4=hits4+1;
                Hits4(hits4)=Event_Analysis(event,error_index);
            end
        end
    end
    % [mean,sig,rms,percent_tails];
    Hits0_Stats=[0,0,0,0];
    Hits1_Stats=[0,0,0,0];
    Hits2_Stats=[0,0,0,0];
    Hits3_Stats=[0,0,0,0];
    Hits4_Stats=[0,0,0,0];
    if hits0>0
        Hits0_Stats=Get_Statistics(Hits0);
    end
    if hits1>0
        Hits1_Stats=Get_Statistics(Hits1);
    end
    if hits2>0
        Hits2_Stats=Get_Statistics(Hits2);
    end
    if hits3>0
        Hits3_Stats=Get_Statistics(Hits3);
    end
    if hits4>0
        Hits4_Stats=Get_Statistics(Hits4);
    end
    Hit_Stats = [Hits0_Stats(2),Hits1_Stats(2),Hits2_Stats(2),Hits3_Stats(2),Hits4_Stats(2)];
    
    scatter([0,1,2,3,4],Hit_Stats)
    % errorbar([0,1,2,3,4],Hit_Stats,Hit_Stats./sqrt([hits0,hits1,hits2,hits3,hits4]))
end
hold off;
figure(457)
hold on;
bar([0,1,2,3,4],[hits0,hits1,hits2,hits3,hits4]./sum([hits0,hits1,hits2,hits3,hits4]));
hold off;
hits0
end

function Get_Background_Plot(openfile)
global strip_width

load(openfile);
background = 'bgon';
charge = 1000;

A = Hits_Data_Set_Time_With_Background;

Event_Info=Event_Info_Alterations(charge,background,openfile);

[m,n] = size(A);
j=0;
% n_events=0;
% event=-1;
for i=2:m
    if A(i,9)==0
        vmm = A(i,4);
        vmm_center = ((vmm-1)*64+32)*strip_width;
        plane = A(i,5);
        strip = A(i,6);
        y = strip*strip_width;
        j=j+1;
        Y(j) = vmm_center/10;  %vmm;  %y;
    end
    %     if A(i,1)~=event
    %         event = A(i,1);
    %         n_events = n_events+1;
    %     end
end
j
N = sum(Event_Info(:,15))
y_min=min(Y);
y_max=max(Y);
n_bins = ceil((y_max-y_min)/(64*strip_width));
hist(Y,300); %n_bins)
hold on;
title('Background model output','FontSize',18);
xlabel('y [mm] (distance from beam along wedge symmetry axis)','FontSize',16);
ylabel('number of background hits in data sample','FontSize',16);
text = sprintf('%i background hits randomly generated for %i events',j,N);
annotation('textbox','String',text,'FitBoxToText','on');

% Event=zeros(1,event_mark);
% Event_bg=zeros(1,event_mark);
% for i=2:m
%     if A(i,9)==1
%         Event(A(i,1))=1;
%     end
%     if A(i,9)==0
%         Event_bg(A(i,1))=Event_bg(A(i,1))+1;
%     end
% end
% NN = sum(Event)
% NNN = sum(Event_bg)
% A(1:100,1:14)/1000
% Event_diff = Event - Event_bg;
% true_only_events=0;
% bg_only_events=0;
% for i=1:event_mark
%     if Event_diff(i)<0
%         bg_only_events = bg_only_events + 1;
%     elseif Event_diff(i)>0
%         true_only_events = true_only_events + 1;
%     end
% end
% true_only_events
% bg_only_events




end  %NO LONGER NECESSARY WITH NEW BG MODEL