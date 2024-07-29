clc
clear
close all

addpath('.\data')
addpath('..\m_IGRF')

data_original_filename = 'Flt1003_train.h5';
time = datenum([2020 6 29]);
lines={1003.02,1003.04,1003.08};
load('model_Flt1002.mat');

% data_original_filename = 'Flt1007_train.h5';
% time = datenum([2020 7 7]);
% lines={1007.02,1007.06};
% load('model_Flt1006.mat');

[x_m,y_m,z_m,mag_earth_intensity]=loadMITData(data_original_filename, lines, time);

%%
% R=[0,1,0;1,0,0;0,0,1];
% matrix=matrix*R;
residual_h_m=zeros(size(x_m));
residual_h_hat=zeros(size(x_m));

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

%     residual_h_m(i_iters)=abs(norm(h_hat)-mag_earth(i_iters));
%     residual_h_hat(i_iters)=abs(norm(h)-mag_earth(i_iters));

    residual_h_m(i_iters)=abs(norm(h_hat)-mag_earth_intensity);
    residual_h_hat(i_iters)=abs(norm(h)-mag_earth_intensity);
end

residual_h_m_mean=mean(residual_h_m);
residual_h_hat_mean=mean(residual_h_hat);

%%
figure;
% Visualization %
% Sensor readings and ellipoid fit
scatter3(x_m, y_m, z_m, 'fill', 'MarkerFaceColor', 'red'); hold on; 
plot_ellipsoid(v,'r'); 

% After calibrations
% figure;
scatter3(x_hat, y_hat, z_hat, 'fill', 'MarkerFaceColor', 'blue'); hold on;
plot_sphere([0,0,0]', mag_earth_intensity);
% title({'After magnetometer calibration', '(Normalized to unit sphere)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;
legend('before calibration: measured data','before calibration: fitted ellipsoid','after calibration: calibrated data','after calibration: sphere');

% Print calibration params
fprintf('3D magnetometer calibration based on ellipsoid fitting');
fprintf('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~');
fprintf('\nThe calibration equation to be implemented:') 
fprintf('\n\t\t\t\th_hat = M*(h_m - b) \nWhere,')
fprintf('\nh_m   = Measured sensor data vector');
fprintf('\nh_hat = Calibrated sensor data vector');
fprintf('\n\nM =\n'); disp(matrix);
fprintf('\nb =\n'); disp(offset);

fprintf('\nresidual_h_m_mean ='); disp(residual_h_m_mean);
fprintf('\nresidual_h_hat_mean ='); disp(residual_h_hat_mean);
