

function Set_Parameters()

global w1 w2 w3 h1 h2 h3 H L wedge_opening_angle strip_width stereo_degree z_large setup
global h_mx h_my m_y_min

global slope_min slope_max h x_error uv_error CT CT_x CT_u CT_v CT_uv

global minimum_large_theta maximum_large_theta minimum_large_phi maximum_large_phi

global n_theta_rois n_phi_rois

global mid_plane_large_X mid_plane_large vertical_strip_width_UV mid_plane_large_UV


'Setting Parameters...'



%y = vertical distance from beamline
%z = into wedge along beamline
%x = horizontal distance from beam looking down

%//////////  Define the large wedge /////////////////
w1=2220;   %top
w2=2540.3;  %determined by 33deg angle   %shelves part
w3=582.3;  %bottom
h1=3665; %how tall wedge is at w1, ie total height
h3=360;
h2=3665-h3; %height at w2

wedge_opening_angle = 33;  %degree

z_large = [7478,7489,7521,7532,7604,7615,7647,7658];

mid_plane_large = mean(z_large);
mid_plane_large_X = (z_large(1)+z_large(2)+z_large(5)+z_large(6))/4;
mid_plane_large_UV = (z_large(3)+z_large(4)+z_large(7)+z_large(8))/4;

H=982; %bottom of wedge to the beamline
L=7500; %distance from IP to front side of wedge

strip_width = 0.445;  % 0.5;
stereo_degree = 1.5;  %in degrees!

degree=degtorad(stereo_degree);
vertical_strip_width_UV = strip_width/cos(degree);

%setup = ['x','x','u','v','v','u','x','x'];
setup = ['x','x','u','v','x','x','u','v'];
%setup = ['x','x','v','u','x','x','v','u'];
% setup = ['x','u','x','v','x','u','x','v'];

%////// TABLE GENERATORS ///////////////  %size of cartesian steps
h_mx = 0.0001; %0.0001; % 0.005;  %0.0001;
h_my =  0.0001; %0.0001; %0.005;  %0.0001;
%///////////////////////////////////////






%////////  for cut applications  /////////////////
tol = 0;  %0.02;
minimum_large_theta = atan((cos(pi/2-degtorad(wedge_opening_angle/2))*w3+H)/z_large(1))+tol;
maximum_large_theta = atan((H+h1)/z_large(8))-tol;
minimum_large_phi = pi/2-degtorad(wedge_opening_angle/2)+tol;
maximum_large_phi = pi/2+degtorad(wedge_opening_angle/2)-tol;
%/////////////////////////////////////////////////



%/////   Rough ROI  ////////////
n_theta_rois = 32;
n_phi_rois = 16;
%////////////////////////////





%///// Front Filter /////////////
theta_max = maximum_large_theta;
theta_min = minimum_large_theta;
slope_max = tan(theta_max);
slope_min = tan(theta_min);
%--- the following slope road and thresholds are often re-defined immediately in any master script running the file
h = 2.5*10^(-4);  %10^(-3);
CT_x = 3;
CT_u = 0;
CT_v = 0;
CT_uv = 2;
CT = CT_x + CT_uv;
x_error = h*0.5;
uv_error = 3.5*10^(-3);
BC_window = 2;
%////////////////////////////////



%/////   Stuff sometimes worth seeing ////////////////
theta_max = maximum_large_theta;   %0.56;  %Check these!!!
theta_min = minimum_large_theta;  %0.125;
M_global_max = tan(theta_max);
M_global_min = tan(theta_min);
atan((H)/z_large(8));
atan((H+h1)/z_large(1));
%/////////////////////////////////////////////////////

