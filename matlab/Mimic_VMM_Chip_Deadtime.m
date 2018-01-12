function Hits_Data_Set_Time = Mimic_VMM_Chip_Deadtime(Hits_Data_Set_Time_Full)

clear rows_to_del Hits_Data_Set_Time

'Mimicing VMM chip 1st-hit deadtime...'

sortrows(Hits_Data_Set_Time_Full,[10,14]);  %10);

VMM_deadtime = 100;  %(ns)

num_VMM_per_plane = 1000;

VMM_chip_status = ones(num_VMM_per_plane,8);
VMM__chip_last_hit_time = zeros(num_VMM_per_plane,8);

[m,n]=size(Hits_Data_Set_Time_Full);

A=Hits_Data_Set_Time_Full;
for i=2:m
    if A(i,6)==0  %clear out non-hits  (i.e. empty events -- probably from cuts in genertion)
        A(i,4)=0;
    else
        VMM_chip = A(i,4);
        plane = A(i,5);
        time = A(i,14);  %true time in sub ns from simulation
        
        if VMM_chip_status(VMM_chip,plane)==0  %is the chip active?
            if VMM__chip_last_hit_time(VMM_chip,plane) + VMM_deadtime <= time  %if not, should the chip be active?
                VMM__chip_last_hit_time(VMM_chip,plane)=0;
                VMM_chip_status(VMM_chip,plane)=1;
            end
        end
        
        if VMM_chip_status(VMM_chip,plane)==1  %is the chip active?
            VMM_chip_status(VMM_chip,plane)=0;
            VMM__chip_last_hit_time(VMM_chip,plane)=time;
        else
            A(i,4)=0;
        end
        
    end
end

%Get rid of hits that would not have registered
j=0;
for i=2:m
    if A(i,4)==0
        j=j+1;
        rows_to_del(j)=i;
    end
end
if j>0
A(rows_to_del,:)=[];
end
Hits_Data_Set_Time = A;
    