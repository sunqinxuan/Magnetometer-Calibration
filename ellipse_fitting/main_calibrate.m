% 2022/09/21

clc
clear
close all

% TODO
% 1. Consistancy of phi for different quadrant
% 2. Pack matrix and scale to single transformation matrix

% Noisy ellipse params
a = 6;      % Semi-major axis
b = 12;     % Semi-minor axis
xc = -20;   % Center X-coordinate  
yc = 15;    % Center Y-coordinate
phi = pi/4; % Angle wrt. X-axis
noise = 4;  % Noise factor

% Generate noisy ellipse
t = 0:0.1:2*pi;
x_raw = xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi) + noise*rand(size(t));
y_raw = yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi) + noise*rand(size(t));

% Fit ellipse
v = ellipse_fit(x_raw', y_raw');
[a_, b_, xc_, yc_, phi_] = ellipse_conic_to_parametric(v);

% Calibration matric and bias
[matrix, scale, bias] = ellipse_to_circle(v);
x_cal = x_raw - bias(1);
y_cal = y_raw - bias(2);
temp = [x_cal; y_cal]'*matrix;
x_cal = temp(:,1).*scale;
y_cal = temp(:,2);

% Ellipse and circle for visualization
t_ = 0:0.1:2*pi;
x_circ = b_*cos(t);
y_circ = b_*sin(t);
x_ellipse = xc_ + a_*cos(t_)*cos(phi_) - b_*sin(t_)*sin(phi_);
y_ellipse = yc_ + a_*cos(t_)*sin(phi_) + b_*sin(t_)*cos(phi_);

% Visualization
figure; hold on; axis equal; grid on;
plot(x_raw, y_raw, 'o'); 
plot(x_ellipse, y_ellipse);
plot(x_cal, y_cal, 'o')
plot(x_circ, y_circ);
xlabel("X"); ylabel("Y");
title("2D magnetometer calibration simulation");
legend("Sensor input", "Ellipse fit", "Calibration", "Output circle");
