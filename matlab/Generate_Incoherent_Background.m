function Generate_Incoherent_Background(index)

global strip_width z_large
global H h1

Set_Parameters();

'Add Salt-and-Pepper Background...'

figure(1);
hold on;

theta_min = 0.1722;
theta_max = 0.5254;
phi_min = 1.3028;
phi_max = 1.8388;

n_events = 1000;

BC_window=2;

%calculate occupancy numbers
r_min=H;
r_max=H+h1;
dr = 10;
r_zones = ceil((r_max - r_min)/dr);
dphi = phi_max - phi_min;
Occupancy_Table = zeros(r_zones,3);
for i=1:r_zones
    r = r_min + (i-1)*dr + dr/2;
    hit_rate_per_cm2 = 2.5*10^8*(r*0.1)^(-2.125);
    
    Area = 0.5*(((r+dr/2)*0.1)^2-((r-dr/2)*0.1)^2)*dphi;
    occupancy = hit_rate_per_cm2*Area/(40*10^6)*BC_window;

%     occupancy = -50*dphi*((1/(r*0.1)^(-0.125)-1/((r+dr)*0.1)^(-0.125)));



     
%     Area = (((r+dr/2)*0.1)^2-((r-dr/2)*0.1)^2)/cot(0.268);
%     psi = pi/2 - phi_min;
%     Area = 2*tan(psi)*dr*r + tan(psi)*dr^2;
%     occupancy = hit_rate_per_cm2*Area/(40*10^6);
% 
%     occupancy = 3.35*(log((r+dr)*0.1)-log(r*0.1));   %3.35*(log((r+dr/2)*0.1)-log((r-dr/2)*0.1));
   
    
    
    Occupancy_Table(i,:) = [r,occupancy,hit_rate_per_cm2];
end
Occupancy_Table(:,2)
Occupancy_Table(:,1)
scatter(Occupancy_Table(:,1),Occupancy_Table(:,2))

b=0;
for event=((index-1)*1000+1):((index-1)*1000+n_events)
    event
    for plane = 1:8
        for i=1:r_zones
            q = rand(1);
            if q <= Occupancy_Table(i,2)
                r = random('unif',Occupancy_Table(i,1)-dr/2,Occupancy_Table(i,1)+dr/2);
                phi = random('unif',phi_min,phi_max);
                y = r*cos(phi-pi/2);
                x = r*sin(pi/2-phi);
                theta = atan(r/z_large(plane));
                
                athena_time = random('unif',0,BC_window*25);
                recon_y = 0;
                recon_x = 0;
                true_x = x;
                true_y = y;
                true_z = z_large(plane);
                
                strip = Get_Strip_ID(x,y,plane);
                VMM_chip = Get_VMM_chip(strip);
                BC_id = ceil(athena_time/25);
                BC_time = event*10 + (BC_id-1);
                truth = 0;
                
                b=b+1;
                Hits(1,:) = [event,r,0,VMM_chip,plane,strip,theta,phi,truth,BC_time,0,0,0];
                Time(1,:) = athena_time + event*100;
                Cart_Hit(1,:) = [true_x,true_y,true_z,recon_x,recon_y];
                
                A(b,:) = [Hits,Time,Cart_Hit];
                
                i=i-1;  %no reason you can't have multiple hits
            end
        end
    end
end

A=[A,zeros(b,11)];


b
Background_Hits=sortrows(A,[10,14]);
filename = sprintf('Background_%ievents_%iBCunif_dr%i_%i',n_events,BC_window,dr,index);
save(filename)

hist(A(:,2),100)

scatter(Occupancy_Table(:,1),Occupancy_Table(:,3))
hold on;
title('Incoherent Background Model','FontSize',18);
xlabel('Radius from beam (mm)','FontSize',16);
ylabel('Hit Rate (Hz/cm^2)','FontSize',16);
% text = ;
% annotation('textbox','String',text,'FitBoxToText','on');

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
    case 'v'
        y_hit = -X*sin(degree)+Y*cos(degree);
        strip_hit = ceil(y_hit*1/strip_width);
end

strip=strip_hit;
end
