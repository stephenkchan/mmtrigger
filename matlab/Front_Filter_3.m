function Front_Filter_3(A,BC_window)   %FIFO and coincidence checking
%A is data set of hits and BC_window is depth of FIFO

global BC_clock
global slope_min slope_max h regions
global Slope_Roads_1 Slope_Roads_2 Slope_Roads_3 Slope_Roads_4 Track_Indexes_1 Track_Indexes_2 Track_Indexes_3 Track_Indexes_4

% global Event_Fit
% %quick fix because originally saved this variable as too small
% [qq,ww]=size(Event_Fit);
% if ww==10
%     Event_Fit = [Event_Fit,zeros(qq,6)];
% end

global N_coincidence_threshold_met
N_coincidence_threshold_met = 0;

global BC_Timer BC_Timer_Activate  %used to control waiting window for hits
BC_Timer=0;
BC_Timer_Activate=0;

'Fitting Events...'

regions=ceil((slope_max - slope_min)/h);

%FIFO (effective size is not necessarily 4 deep depending on BC_window)
%slope_Roads_# form the coincidence FIFO
Slope_Roads_1 = zeros(regions,8);
Slope_Roads_2 = zeros(regions,8);
Slope_Roads_3 = zeros(regions,8);
Slope_Roads_4 = zeros(regions,8);
%follow the coincidence FIFO, but serve to keep track of hit identification
Track_Indexes_1 = zeros(regions,8,3);
Track_Indexes_2 = zeros(regions,8,3);
Track_Indexes_3 = zeros(regions,8,3);
Track_Indexes_4 = zeros(regions,8,3);

%How many total hits to interrogaate?
[m,ll]=size(A);

mark=2;   %The first line of the data set is not a hit, but rather a storage line for data set information, which is no longer utilized
BC_clock_min = A(2,10);
BC_clock_max = A(m,10)+4;

for BC_clock = BC_clock_min:BC_clock_max
    
    %BC_Timer and BC_Timer_Activate serve to efffectively shrink the FIFO depth, or hit collection window, from 4 to
    %3,2,or 1 BC's
    if BC_Timer_Activate==1  %activated when first slope is stored and deactivated when slopes are read
        BC_Timer = BC_Timer+1;
    end
    
    %move FIFO along...
    Slope_Roads_4 = Slope_Roads_3;
    Slope_Roads_3 = Slope_Roads_2;
    Slope_Roads_2 = Slope_Roads_1;
    Slope_Roads_1 = zeros(regions,8);
    Track_Indexes_4(:,:,:) = Track_Indexes_3(:,:,:);
    Track_Indexes_3(:,:,:) = Track_Indexes_2(:,:,:);
    Track_Indexes_2(:,:,:) = Track_Indexes_1(:,:,:);
    Track_Indexes_1(:,:,:) = zeros(regions,8,3);
    
    %pulls hits at the current BC clock time
    if BC_clock == A(mark,10)
        for i=mark:m
            if A(i,10)<=BC_clock
                Fill_Slope_Roads_1(A(i,:),i);
            else
                mark=i;
                break;
            end
        end
    end
    
    %Checks FIFO for tracks and acts accordingly
    if BC_Timer==BC_window
        Track=Check_Counters(); %check fifo for coincidence
        if sum(Track(1,:))~=0 %i.e. there are hits meeting coincidence thresholds)
%             Get_Fit(Track); %fit track to determine roi and dtheta
        end
        BC_Timer=0;
        BC_Timer_Activate=0;
    end
    
    
end


end



function Fill_Slope_Roads_1(B,hit_index) %records hits in FIFO

global Slope_Roads_1 Track_Indexes_1
global slope_min slope_max h regions x_error uv_error
global z_large setup strip_width vertical_strip_width_UV
q=regions;
z_hit = z_large;   %NEED TO MAKE DYNAMIC

%activate BC_Timer if not done so already, and advance by one
global BC_Timer BC_Timer_Activate
if BC_Timer_Activate==0
    BC_Timer = BC_Timer+1;
    BC_Timer_Activate=1;
end


[n,l]=size(B);
for i=1:n %originally intended to fill many at once -- currently n==1 always, i.e. one hit at a time
    
    strip = B(i,6);
    plane = B(i,5);
    z = z_hit(plane);
    
    slope = strip*strip_width/z;   %could be more effective to treat UV and X strip_widths as different -- NEED TO CHECK
    
    
    %what type of plane?
    switch setup(plane)
        case 'x'
            stereo=0;
            %             slope = strip*strip_width/z;
        case 'u'
            stereo=1;
            %             slope = strip*vertical_strip_width_UV/z;
        case 'v'
            stereo=-1;
            %             slope = strip*vertical_strip_width_UV/z;
        otherwise
            stereo=999;
    end
    
    %determine slope tolerance depending on plane type
    if stereo==1 || stereo==-1   %plane==3 || plane==4 || plane==7 || plane==8
        err = uv_error;
    elseif stereo==0
        err = x_error;
%         if slope>0.33
%             err = x_error*1.5;
%         end
    else
        'ERROR'
        return;
    end
    
    %slope road boundaries based on hit_slope +/- tolerance
    slope_p = slope + err;
    slope_m = slope - err;
    s_p = round((slope_p - slope_min)/h);
    if s_p>q
        s_p=q;
    end
    s_m = round((slope_m - slope_min)/h);
    if s_m<1
        s_m=1;
    end
    
    %fill multiple slope roads based on tolerance
    for k=s_m:s_p
        Slope_Roads_1(k,plane)=1;
        Track_Indexes_1(k,plane,:) = [strip,slope,hit_index];
    end
    
end

end


function output = Check_Counters()  %checks FIFO for coincidence and clears
global  setup Slope_Roads_1 Slope_Roads_2 Slope_Roads_3 Slope_Roads_4 Track_Indexes_1 Track_Indexes_2 Track_Indexes_3 Track_Indexes_4 regions
global CT CT_x CT_u CT_v CT_uv
global N_coincidence_threshold_met
q=regions;
coder.extrinsic('clear')
clear Track count a b c d
Track=zeros(3,8);

%counting coincidence numbers
count=sum(Slope_Roads_1,2)+sum(Slope_Roads_2,2)+sum(Slope_Roads_3,2)+sum(Slope_Roads_4,2);
count_x = sum(Slope_Roads_1(:,[1,2,5,6]),2)+sum(Slope_Roads_2(:,[1,2,5,6]),2)+sum(Slope_Roads_3(:,[1,2,5,6]),2)+sum(Slope_Roads_4(:,[1,2,5,6]),2);
count_uv = sum(Slope_Roads_1(:,[3,4,7,8]),2)+sum(Slope_Roads_2(:,[3,4,7,8]),2)+sum(Slope_Roads_3(:,[3,4,7,8]),2)+sum(Slope_Roads_4(:,[3,4,7,8]),2);

for w=1:q
    if count(w,1)>=CT && count_x(w,1)>=CT_x && count_uv(w,1)>=CT_uv  %i.e. threshold met
        
        %check to see if there is a newarby row with a higher coincidence count
        neighbers_to_check=100;
        i=w; %primary index I work with going forward
        if q-w>=neighbers_to_check
            r=neighbers_to_check;
        else
            r=q-w;
        end
        max_total=count(w,1);
        max_x=count_x(w,1);
        max_uv=count_uv(w,1);
        for k=w:w+r
            if count(k,1)>=max_total && count_x(k,1)>=max_x && count_uv(k,1)>=max_uv
                i=k;
                max_total=count(k,1);
                max_x=count_x(k,1);
                max_uv=count_uv(k,1);
            end
        end
        
        for j=1:8
            if Track_Indexes_4(i,j,1)~=0
                Track(:,j) = Track_Indexes_4(i,j,:);
            elseif Track_Indexes_3(i,j,1)~=0
                Track(:,j) = Track_Indexes_3(i,j,:);
            elseif Track_Indexes_2(i,j,1)~=0
                Track(:,j) = Track_Indexes_2(i,j,:);
            elseif Track_Indexes_1(i,j,1)~=0
                Track(:,j) = Track_Indexes_1(i,j,:);
            end
        end
        
        %             %$$$$$$$$$$$$  DEBUGGING TOOL   $$$$$$$$$$$
        %             for j=1:8
        %             if Track(3,4)~=0 && A(Track(3,4),1)==181
        %             Slope_Roads_4(i-2:i+2,:)
        %             Slope_Roads_3(i-2:i+2,:)
        %             Slope_Roads_2(i-2:i+2,:)
        %             Slope_Roads_1(i-2:i+2,:)
        %             pause;
        %             break;
        %             end
        %             end
        %             %$$$$$$$$$$$$     $$$$$$$$$$$
        
        %------ confiming coincidence thresholds met -----  REDUNDANT ARTIFACT!!!
        g=0;
        for j=1:8
            if Track(1,j)~=0
                g=g+1;
            end
        end
        if g<CT
            Track=zeros(3,8);
            break;
        end
        gx=0;
        gu=0;
        gv=0;
        for j=1:8
            if Track(1,j)~=0
                switch setup(j)
                    case 'x'
                        gx=gx+1;
                    case 'u'
                        gu=gu+1;
                    case 'v'
                        gv=gv+1;
                end
            end
        end
        if (gu+gv) < CT_uv
            Track=zeros(3,8);
            break;
        end
        if gx < CT_x || gu < CT_u || gv < CT_v
            Track=zeros(3,8);
            break;
        end
        %---------------------
        
        
        
        %served to erase only one part of the FIFO rather than the whole
        %FIFO once a track was discovered
        %             err = ceil(x_error/h);
        %             a=i-err;
        %             b=i+err;
        %             if b>q
        %                 b=q;
        %             end
        %             if a<1
        %                 a=1;
        %             end
        
        %             for k=1:q  %k=a:b
        %                 Slope_Roads_1(k,:)=zeros(1,8);
        %                 Slope_Roads_2(k,:)=zeros(1,8);
        %                 Slope_Roads_3(k,:)=zeros(1,8);
        %                 Slope_Roads_4(k,:)=zeros(1,8);
        %                 Track_Indexes_1(k,:,:)=zeros(1,8,3);
        %                 Track_Indexes_2(k,:,:)=zeros(1,8,3);
        %                 Track_Indexes_3(k,:,:)=zeros(1,8,3);
        %                 Track_Indexes_4(k,:,:)=zeros(1,8,3);
        %             end
        
        
        
        %Currently, clearing the entire FIFO.
        Slope_Roads_1=zeros(regions,8);
        Slope_Roads_2=zeros(regions,8);
        Slope_Roads_3=zeros(regions,8);
        Slope_Roads_4=zeros(regions,8);
        Track_Indexes_1=zeros(regions,8,3);
        Track_Indexes_2=zeros(regions,8,3);
        Track_Indexes_3=zeros(regions,8,3);
        Track_Indexes_4=zeros(regions,8,3);
        
        %keeping a tally for analysis
        N_coincidence_threshold_met = N_coincidence_threshold_met+1;

        output = Track;  %Output the track for fitting
        return;

    end
end

output = Track;  %this is a null track because coincidence thresholds were not satisfied

end


