% 2022/10/05

clc
clear
close all

% Import raw magnetometer readings
file = 'sensor_data.csv'; 

% Import sensor readings
raw = importdata(file);
x_m = raw(:,1); 
y_m = raw(:,2);

% Fit ellipse
v = ellipse_fit(x_m, y_m);
[a_, b_, xc_, yc_, phi_] = ellipse_conic_to_parametric(v);

% Calibration matric and bias
[matrix, scale, bias] = ellipse_to_circle(v);
x_hat = x_m - bias(1);
y_hat = y_m - bias(2);
temp = [x_hat'; y_hat']'*matrix;
x_hat = temp(:,1).*scale(1);
y_hat = temp(:,2).*scale(2);
csvwrite("calibrated_data.csv", [x_hat'; y_hat']');

% Ellipse and circle for visualization
t_ = 0:0.1:2*pi;
x_circ = cos(t_);
y_circ = sin(t_);
x_ellipse = xc_ + a_*cos(t_)*cos(phi_) - b_*sin(t_)*sin(phi_);
y_ellipse = yc_ + a_*cos(t_)*sin(phi_) + b_*sin(t_)*cos(phi_);

% Visualization %
% Before calibration
figure;
plot(x_m, y_m, 'o', 'MarkerFaceColor','red'); 
hold on; axis equal; grid on;
plot(x_ellipse, y_ellipse, 'color', [0 0.4470 0.7410]);
xlabel("X"); ylabel("Y");
title({"Before calibation", "(Ellipse fit)"});

% After calibration
figure;
plot(x_circ, y_circ, 'color', 'red');
hold on; axis equal; grid on;
plot(x_hat, y_hat, 'o', 'MarkerFaceColor', [0 0.4470 0.7410]);
xlabel("X"); ylabel("Y");
title({"After calibation", "(Normalized to unit circle)"});
