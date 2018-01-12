function Table_Generators()

'Generating Tables...'

%Reference for local slope least squares fit
Local_Slope_A_B();

%cartesian slopes to phi and theta reference table
Slope_Components_ROI();

%Delta theta reference table to avoid division operation
Delta_theta_optimization_LG();

end



function Local_Slope_A_B()
clear A_local B_local
global A_local B_local
global z_large setup strip_width

z_hit = z_large;

switch setup
    case 'xxuvvuxx'
        planes = [1,2,7,8];
    case 'xuxvxuxv'
        planes = [1,3,5,7];
    case 'xxuvxxuv'
        planes = [1,2,5,6];
    case 'xxvuxxvu'
        planes = [1,2,5,6];
end


A_local = zeros(1,11);
B_local = zeros(1,11);

for k=1:11
    
    if k==1
        N=4;
        x=[1,2,3,4];
    elseif k==2
        N=3;
        x=[1,2,3];
    elseif k==3
        N=3;
        x=[1,2,4];
    elseif k==4
        N=3;
        x=[1,3,4];
    elseif k==5
        N=3;
        x=[2,3,4];
    elseif k==6
        N=2;
        x=[1,2];
    elseif k==7
        N=2;
        x=[1,3];
    elseif k==8
        N=2;
        x=[1,4];
    elseif k==9
        N=2;
        x=[2,3];
    elseif k==10
        N=2;
        x=[2,4];
    elseif k==11
        N=2;
        x=[3,4];
    end
    
    sum_x=0;
    sum_xx=0;
    for i=1:N
        j=x(i);
        sum_x = sum_x + z_hit(planes(j));
        sum_xx = sum_xx + z_hit(planes(j))^2;
    end
    
    D = N * sum_xx - sum_x^2;
    
    A_local(k) = N / D * strip_width;
    B_local(k) = sum_x / D * strip_width;
    
end

format short e
A_local
B_local
z_large

end


function Slope_Components_ROI()
clear Slope_to_ROI M
global H h1 w2 z_large
global Slope_to_ROI h_mx h_my m_y_min
global minimum_large_theta maximum_large_theta minimum_large_phi maximum_large_phi

z_hit=z_large;

%these are mm/mm i.e. not units of strip #
m_y_max = (H+h1)/z_hit(1); %2;
m_y_min = H/z_hit(8);  %
m_x_max = (w2/2)/z_hit(1); %2; %
m_x_min = (-w2/2)/z_hit(1); %-2; %



n_x = ceil(m_x_max/h_mx)+1;
n_y = ceil((m_y_max - m_y_min)/h_my);

for i=1:n_x
    i;
   for j=1:n_y
       
       m_x = (i-1)*h_mx;
       m_y = m_y_min+j*h_my;
       
       theta = atan(sqrt(m_x^2+m_y^2));
       phi = atan(m_y/m_x);
       
       
       if theta<minimum_large_theta || theta>maximum_large_theta || phi<minimum_large_phi || phi>maximum_large_phi
            theta = 0;
            phi = 0;
       end
       
       
       M(j,i,1) = theta;
       M(j,i,2) = phi;
       
   end
end

Slope_to_ROI = M;

end


function Delta_theta_optimization_LG()

global DT_Factors
global computer
global minimum_large_theta maximum_large_theta minimum_large_phi maximum_large_phi

a=1;   %sin(pi/2+degtorad(28/4));  %sin(phi);

% Delta_Theta = (M_local - M_global) / (a + M_local*M_global/a);

theta_max = maximum_large_theta;   %0.56;  %Check these!!!
theta_min = minimum_large_theta;  %0.125;
M_global_max = tan(theta_max);
M_global_min = tan(theta_min);

LG_min = 0;  %ie.e no neg slopes
LG_max = 0.5;  %
number_LG_regions = 256;
LG_region_width = (LG_max - LG_min)/number_LG_regions;

for i=1:number_LG_regions
    
    LG = LG_min + (LG_region_width * (i-1));
    mult_factor = 1 / (a + LG/a);
    
    Mult_Factors(i,:) = [LG,mult_factor];
    
end

DT_Factors = Mult_Factors;


%Everything past here is debugging nonsense


D_limit = 20 * 10^(-2);
DM_cut_200GeV = 1000; %0.07;
N=1000;
D_theta_res = zeros(1,N);
D_theta_truth = zeros(1,N);
D_theta_calc = zeros(1,N);
DM = zeros(1,N);
for k=1:N
    
    M_global = random('unif',M_global_min,M_global_max);
    D_theta = random('unif',D_limit*-1,1*D_limit);  %random('unif',20*10^(-3),D_limit);  %random('unif',D_limit,100*10^(-3));
    phi = random('unif',minimum_large_phi,maximum_large_phi);
    M_local = (sin(phi)*D_theta + M_global)/(1-D_theta*M_global/sin(phi));
    
    if (M_local - M_global) <= DM_cut_200GeV
        D_theta_truth(k) = (M_local - M_global) / (sin(phi) + M_local*M_global/sin(phi));
        if M_global >= M_global_min && M_global <= M_global_max
            for j=1:number_LG_regions
                if LG >= Mult_Factors(j,1)
                    region = j;
                    break;
                end
            end
        end
        D_theta_calc(k) = (M_local - M_global) * Mult_Factors(region,2);
        DM(k) = M_local - M_global;
    else
        k=k-1;
    end
    
end

D_theta_res = D_theta_calc - D_theta_truth;

% %hist(D_theta,100);
% % histfit(D_theta_res,1000);  %round(sqrt(N)))
% figure(3);
% hist(D_theta_res,1000);
% figure(4);
% hist(D_theta_calc,1000);
% % figure(3);
% % hist(DM,1000);

% savefile = sprintf('/Users/%s/Desktop/Delta_Theta_Mult_Factors_LG.txt',computer);
% dlmwrite(savefile,DT_Factors);

end
