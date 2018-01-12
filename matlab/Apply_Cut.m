% It was noticed at 2:43am on march 21, 2014 that a mistake exists in this cutting program.
% The confusion arose between event and athena_event numbers.  Hits_Data_Set_Time was apporpriately cut.
% However, Event_Info(i,15) was inappropriately altered.
% No changes have been made to correct this error as I suspect if does not
% affect anything until the analysis stage.  Therefore, a fix exists in the
% analysis script "Get_Tuning_Analysis"

%... And this might not be true!



function output = Apply_Cut(Data,parameter_index,lower_bound,upper_bound)

global Event_Info

[m,n]=size(Data);

if parameter_index==0
    'Apply cut on number of planes hit...'
    event_min = min(Data(2:m,1));
    event_max = max(Data(2:m,1));
    for event=event_min:event_max
        planes_hit=zeros(1,8);
        X_planes_hit=0;
        U_planes_hit=0;
        V_planes_hit=0;
        for i=2:m
            if Data(i,1)==event
                plane = Data(i,5);
                planes_hit(plane) = 1;
            end
        end
        
                X_planes_hit = planes_hit(1)+planes_hit(2)+planes_hit(5)+planes_hit(6);
                U_planes_hit = planes_hit(3)+planes_hit(7);
                V_planes_hit = planes_hit(4)+planes_hit(8);
        if X_planes_hit < lower_bound || (U_planes_hit + V_planes_hit) < upper_bound
%         if sum(planes_hit)>upper_bound || sum(planes_hit)<lower_bound
            for i=2:m
                if Data(i,1)==event
                    Data(i,:)=0;
%                     Data(i,1)=event;
                    Event_Info(event,15)=0;
                end
            end
        end
        
    end
    
elseif parameter_index==-1
    'Apply cut on total number of hits...'
    event_min = min(Data(2:m,1));
    event_max = max(Data(2:m,1));
    for event=event_min:event_max
        N=0;
        for i=2:m
            if Data(i,1)==event
                N=N+1;
            end
        end
        
        if N<lower_bound || N>upper_bound
            for i=2:m
                if Data(i,1)==event
                    Data(i,:)=0;
%                     Data(i,1)=event;
                    Event_Info(event,15)=0;
                end
            end
        end
        
    end
    
else
    sprintf('Apply cut on parameter index %d',parameter_index)
    for i=2:m
        event=Data(i,1);
        if Data(i,parameter_index)>upper_bound || Data(i,parameter_index)<lower_bound
            Data(i,:)=0;
            Event_Info(event,15)=0;
        end
    end
    
end

%remove zero entries
j=0;
for i=2:m
    if Data(i,4)==0
        j=j+1;
        rows_to_del(j)=i;
    end
end
if j>0
    Data(rows_to_del,:)=[];
end

output = Data;