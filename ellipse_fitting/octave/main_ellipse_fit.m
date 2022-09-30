% 2022/09/11

clc
clear
close all

% Noisy ellipse params
a = 6;      % Semi-major axis
b = 12;     % Semi-minor axis
xc = -20;   % Center X-coordinate  
yc = 15;    % Center Y-coordinate
phi = pi/4; % Angle wrt. X-axis
noise = 4;  % Noise factor

% Generate noisy ellipse
t = 0:0.1:2*pi;
x = xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi) + noise*rand(size(t));
y = yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi) + noise*rand(size(t));

% Fit ellipse
v = ellipse_fit(x', y');
[a_, b_, xc_, yc_, phi_] = ellipse_conic_to_parametric(v);

% Output ellipse
t_ = 0:0.01:2*pi;
x_ = xc_ + a_*cos(t_)*cos(phi_) - b_*sin(t_)*sin(phi_);
y_ = yc_ + a_*cos(t_)*sin(phi_) + b_*sin(t_)*cos(phi_);

% Visualization
figure; hold on; axis equal; grid on;
plot(x, y, 'o'); 
plot(x_, y_);
legend("Input", "Ellipse fit");
xlabel("X"); ylabel("Y");
