% Test of ellipsoid_fit.m function
clc
clear 
close all

% Noisy ellipsoid params %
% Semi principal axes
ax= 20;  
bx = 30; 
cx = 50;

% Centre  
xc = 10; 
yc = 100; 
zc = 10;

% Orientation (3-2-1 Euler)
yaw = pi/4; 
pitch = pi/3; 
roll = pi/7;

% Noise factor
noise_factor = 5;

% Generate test ellipoid %
% x	=	a*cosu*sinv 	
% y	=	b*sinu*sinv	
% z	=	c*cosv
% u in [0,2pi) and v in [0,pi]
[u, v] = meshgrid(0:0.3:pi*2, 0:0.3:pi);
x = ax*cos(u).*cos(v);
y = bx*cos(u).*sin(v);
z = cx*sin(u);

% Meshgrid to vector
x=x(:); y=y(:); z=z(:); 
xyz = [x y z];
Q = dcm_321_euler(yaw, pitch, roll);

% Rotate using DCM (321)
for i_iters = 1: length(x)
     new = Q*xyz(i_iters,:)';
     xyz(i_iters,:) = new'; 
end

% Move centre after rotation
x = xc + xyz(:,1);
y = yc + xyz(:,2);
z = zc + xyz(:,3);

% Add noise to generated points
x = x + noise_factor*rand(size(x));
y = y + noise_factor*rand(size(y));
z = z + noise_factor*rand(size(z));

% Ellipsoid fit algoritm
v = ellipsoid_fit(x,y,z);

% Plot result
plot3(x,y,z,'.','MarkerSize',15); hold on; 
plot_ellipsoid(v);
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
title('Ellipsoid fit');
grid on;
