% Batch 3-axis compass calibration using least squares ellipsoid fitting
%
% References:
%   [1] Renaudin - Complete Triaxis Magnetometer Calibration in the 
%                   Magnetic Domain (2010)
%
% This implementation relies on the theory explained in
%   https://teslabs.com/articles/magnetometer-calibration/ 
%
% 2020/06/03
clc
clear
close all

% Import raw magnetometer readings
file = 'raw_data.txt'; 

% SWITCH = 0 Calibrate to unit circle
% SWITCH = 1 Do not change radius
SWITCH = 0;

% Uncalibrated readings
raw = importdata(file);
x_hat = raw(:,1); 
y_hat = raw(:,2); 
z_hat = raw(:,3);

% Calibrated data vectors 
X_cal = zeros(length(x_hat),1); 
Y_cal = zeros(length(x_hat),1); 
Z_cal = zeros(length(x_hat),1);

% Ellipsoid fit
% ax^2 + by^2 + cz^2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0
% v = [a, b, c, f, g, h, p, q, r, d]' (in the paper k = -d)
% Q = [a h g; h b f; g f c]
% u = [p, q, r]'
v = ellipsoid_fit(x_hat,y_hat,z_hat);

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
pqr_ = [p,q,r]*evec;   
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

Ainv  = gain*rotation;
for i_iters = 1:length(x_hat)
    % Sensor data
    h_hat = [x_hat(i_iters); y_hat(i_iters); z_hat(i_iters)]; 
    
    % Calibration, Eqn(11)
    h = gain*rotation*(h_hat - offset);
    
    % Calibrated values
    X_cal(i_iters) = h(1);
    Y_cal(i_iters) = h(2);
    Z_cal(i_iters) = h(3);
end

% Plot uncalibrated data
##subplot(1,2,1);
plot_ellipsoid(v); 
##hold on;

scatter3(x_hat,y_hat,z_hat,'fill','MarkerFaceColor','red');
title({'Before magnetometer calibration','(Ellipsoid fitted)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

##% Plot calibrated data
##subplot(1,2,2);
##plot_sphere([0,0,0],radius);
##hold on;

figure;
scatter3(X_cal,Y_cal,Z_cal,'fill','MarkerFaceColor','blue');
if SWITCH == 0
    title({'After magnetometer calibration','(Normalized to unit circle)'});
else
    title({'After magnetometer calibration'});
end
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

% Print calibration matrices
fprintf('3D magnetometer calibration based on ellipsoid fitting');
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('\nThe calibration equation to be implemented:') 
fprintf('\n\t\t\t\th_hat = M*(h_m - b) \nWhere,')
fprintf('\nh_m   = Measured sensor data vector');
fprintf('\nh_hat = Calibrated sensor data vector');
fprintf('\n\nM =\n'); disp(Ainv);
fprintf('\nb =\n'); disp(offset);
