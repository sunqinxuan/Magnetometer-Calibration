% Test of ellipsoid_fit.m function
% 2020/09/30

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

% Simulate sensor measurements %
% x_m	=	a*cosu*sinv 	
% y_m	=	b*sinu*sinv	
% z_m	=	c*cosv
% u in [0,2pi) and v in [0,pi]
[u, v] = meshgrid(0:0.3:pi*2,0:0.3:pi);
x_m = ax*cos(u).*cos(v);
y_m = bx*cos(u).*sin(v);
z_m = cx*sin(u);

% Meshgrid to vector
x_m=x_m(:); y_m=y_m(:); z_m=z_m(:); 
xyz = [x_m y_m z_m];
Q = dcm_321_euler(yaw, pitch, roll);

% Rotate using DCM (321)
for i_iters = 1: length(x_m)
     new = Q*xyz(i_iters,:)';
     xyz(i_iters,:) = new'; 
end

% Move centre after rotation
x_m = xc + xyz(:,1);
y_m = yc + xyz(:,2);
z_m = zc + xyz(:,3);

% Add noise to generated points
x_m = x_m + noise_factor*rand(size(x_m));
y_m = y_m + noise_factor*rand(size(y_m));
z_m = z_m + noise_factor*rand(size(z_m));

% Ellipsoid fit algoritm
v = ellipsoid_fit(x_m,y_m,z_m);
 
% Unpack ellipsoid coefficients
a = v(1); b = v(2); c = v(3);
f = v(4); g = v(5); h = v(6); 
p = v(7); q = v(8); r = v(9); 
d = v(10); 

% Coordinate frame transformation i.e diagonalize M 
M =[a, h, g; h, b, f; g, f, c]; % Original ellipsoid matrix 
u = [p, q, r]';
k = d;

[evec, eval]=eig(M); % Compute eigenvectors matrix
rotation = evec'; % DCM = eigenvectors matrix
eval = -eval;
M_ = evec'*M*evec; % Diagonalize M

% Coefficients of the ellipsoid in new frame
% Note the ellipsoid is not rotating in this new frame so f, g and h = 0 
pqr_ = [p,q,r]*evec;   
a_ = M_(1,1);
b_ = M_(2,2);
c_ = M_(3,3);
p_ = pqr_(1);
q_ = pqr_(2);
r_ = pqr_(3);
d_ = d;

% Semi principal axes (Still no rotation)
ax_ = sqrt(p_^2/a_^2 + q_^2/(a_*b_) + r_^2/(a_*c_) - d_/a_);
bx_ = sqrt(p_^2/(a_*b_) + q_^2/b_^2 + r_^2/(b_*c_) - d_/b_);
cx_ = sqrt(p_^2/(a_*c_) + q_^2/(b_*c_) + r_^2/c_^2 - d_/c_);

offset = - M \ u; % Eqn(21)
gain = [1/ax_, 0, 0; 0, 1/bx_, 0; 0,0,1/cx_];
matrix = gain*rotation;

% Calibration %
% Memory to calibrated readings 
x_hat = zeros(length(x_m),1); 
y_hat = zeros(length(x_m),1); 
z_hat = zeros(length(x_m),1);
for i_iters = 1:length(x_m)
    % Sensor data
    h_hat = [x_m(i_iters); y_m(i_iters); z_m(i_iters)]; 
    
    % Calibration, Eqn(11)
    h = matrix*(h_hat - offset);
    
    % Calibrated values
    x_hat(i_iters) = h(1);
    y_hat(i_iters) = h(2);
    z_hat(i_iters) = h(3);
end

% Visualization %
% Sensor readings and ellipoid fit
scatter3(x_m, y_m, z_m, 'fill', 'MarkerFaceColor', 'red'); hold on; 
plot_ellipsoid(v); 
title({'Before magnetometer calibration', '(Ellipsoid fit)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

% After calibrations
figure;
scatter3(x_hat, y_hat, z_hat, 'fill', 'MarkerFaceColor', 'blue'); hold on;
plot_sphere([0,0,0]', 1);
title({'After magnetometer calibration', '(Normalized to unit sphere)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;
