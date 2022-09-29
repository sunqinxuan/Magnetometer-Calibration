% Test of ellipsoid_fit.m function
clc
clear 
close all

%%%%% Ellipsoid parameters for test %%%%%
% Semi principal axes
ax= 20; 
bx = 30; 
cx = 50;  

% Centre
xc = 10; 
yc = 100; 
zc = 10;  
centre = [xc,yc,zc];

% Rotation angles (3-2-1)
yaw = pi/4; 
pitch = pi/3; 
roll = pi/7;  

% Signal to noise ratio
SNR = 1;

%%%%% Dont change beyond this point %%%%%
% Generate test points
% x	=	a*cosu*sinv 	
% y	=	b*sinu*sinv	
% z	=	c*cosv
% u in [0,2pi) and v in [0,pi]
[u, v] = meshgrid(0:0.3:pi*2,0:0.3:pi);
x = ax*cos(u).*cos(v);
y = bx*cos(u).*sin(v);
z = cx*sin(u);

% Meshgrid to vector
x=x(:); y=y(:); z=z(:); 
xyz = [x y z];

% Rotate using DCM (321)
C = dcm321Euler(yaw,pitch,roll);
for i_iters = 1: length(x)
     new = C*xyz(i_iters,:)';
     xyz(i_iters,:) = new'; 
end

% Move centre after rotation
x = xc + xyz(:,1);
y = yc + xyz(:,2);
z = zc + xyz(:,3);

% Add noise to generated points
x = x + SNR*rand(size(x));
y = y + SNR*rand(size(y));
z = z + SNR*rand(size(z));

% Ellipsoid fit algoritm
v = ellipsoid_fit(x,y,z);

% Plot and camera settings
% set(gca,'NextPlot','add', 'Visible','off'); view(59,13); hold on;
plot_ellipsoid(v);
plot3(x,y,z,'.','MarkerSize',15); hold on; % Input points
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Ellipsoid fitting');
grid on;
