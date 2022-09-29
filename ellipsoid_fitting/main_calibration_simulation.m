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
 
% Unpack ellipsoid coefficients
a = v(1); b = v(2); c = v(3);
f = v(4); g = v(5); h = v(6); 
p = v(7); q = v(8); r = v(9); 
d = v(10); 

% Coordinate frame transformation i.e diagonalize M 
Q =[a, h, g; h, b, f; g, f, c]; % Original ellipsoid matrix 
u = [p, q, r]';
k = d;

[evec, eval]=eig(Q); % Compute eigenvectors matrix
rotation = evec'; % DCM = eigenvectors matrix
offset = - Q \ u; % Eqn(21)
eval = -eval
M_ = evec'*Q*evec; % Diagonalize M

% Coefficients of the ellipsoid in new frame
% Note the ellipsoid is not rotating in this new frame so f, g and h = 0 
pqr_ = [p,q,r]*C;   
a_ = M_(1,1);
b_ = M_(2,2);
c_ = M_(3,3);
p_ = pqr_(1);
q_ = pqr_(2);
r_ = pqr_(3);
d_ = d;

% Semi principal axes (Still no rotation)
ax_ = sqrt(p_^2/a_^2 + q_^2/(a_*b_) + r_^2/(a_*c_) - d_/a_)
bx_ = sqrt(p_^2/(a_*b_) + q_^2/b_^2 + r_^2/(b_*c_) - d_/b_)
cx_ = sqrt(p_^2/(a_*c_) + q_^2/(b_*c_) + r_^2/c_^2 - d_/c_)
gain = [1/ax_, 0, 0; 0, 1/bx_, 0; 0,0,1/cx_];

for i_iters = 1:length(x)
    % Sensor data
    h_hat = [x(i_iters); y(i_iters); z(i_iters)]; 
    
    % Calibration, Eqn(11)
    h = gain*rotation*(h_hat - offset);
    
    % Calibrated values
    x_cal(i_iters) = h(1);
    y_cal(i_iters) = h(2);
    z_cal(i_iters) = h(3);
end

% Plot and camera settings
% set(gca,'NextPlot','add', 'Visible','off'); view(59,13); hold on;
figure; grid on; hold on;
plot_ellipsoid(v);
plot3(x,y,z,'.','MarkerSize',15); hold on; % Input points
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
title('Ellipsoid fitting');

figure; 
plot3(x_cal,y_cal,z_cal,'.','MarkerSize',15); hold on;
title({'After magnetometer calibration','(Normalized to unit circle)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal; grid on;


