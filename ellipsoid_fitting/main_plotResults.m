clc
clear
close all

%% 1st flight

% data_UDAU=readData_UDAU('.\data\uav\1\mag13\UDAU_data.txt');
% save_mat_name='model_uav_1_mag13.mat';
% x_m=data_UDAU(1688:9391,2);
% y_m=data_UDAU(1688:9391,3);
% z_m=data_UDAU(1688:9391,4);

%% 2nd flight

% data_UDAU=readData_UDAU('.\data\uav\2\mag13\UDAU_data.txt');
% save_mat_name='model_uav_2_mag13.mat';
% x_m=data_UDAU(1:7754,2);
% y_m=data_UDAU(1:7754,3);
% z_m=data_UDAU(1:7754,4);

%% 3rd flight

data_csv=readData_csv('.\data\uav\3\mag13\Mag13_2024-07-10-09-48-03.csv');
save_mat_name='model_uav_3_mag13.mat';
x_m=table2array(data_csv(1107:10071,2));
y_m=table2array(data_csv(1107:10071,3));
z_m=table2array(data_csv(1107:10071,4));

%%
mag_m=[];
for i=1:size(x_m,1)
    m=[x_m(i),y_m(i),z_m(i)];
%     sum=sum+norm(m);
    mag_m(i)=norm(m);
end
% mag_earth_intensity=sum/size(x_m,1); 
mag_earth_intensity=mean(mag_m);


%%

% data_original_filename = 'UDAU_data.txt';
% time = datenum([2020 6 20]); 
% lines={1002.02,1002.20};
% 
% % data_original_filename = 'Flt1006_train.h5';
% % time = datenum([2020 7 6]); 
% % lines={1006.04,1006.04,1006.08};
% 
% [x_m,y_m,z_m,mag_earth_intensity]=loadMITData(data_original_filename, lines, time);

%%

% Ellipsoid fit
% ax^2 + by^2 + cz^2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0
% v = [a, b, c, f, g, h, p, q, r, d]' (in the paper k = -d)
% M = [a h g; h b f; g f c]
% u = [p, q, r]'
v = ellipsoid_fit(x_m, y_m, z_m);

[matrix,offset]=calculateCalibCoeffs(v,mag_earth_intensity);

% cell_str=strsplit(data_original_filename,'_');
% save_mat_name=['model_',cell_str{1,1},'.mat'];
save(save_mat_name,'matrix','offset','v');

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
% title({'Before magnetometer calibration', '(Ellipsoid fit)'});
% xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
% axis equal;

% After calibrations
% figure;
scatter3(x_hat, y_hat, z_hat, 'fill', 'MarkerFaceColor', 'blue'); hold on;
v2 = ellipsoid_fit(x_hat, y_hat, z_hat);
plot_ellipsoid(v2,'b'); 
plot_sphere([0,0,0]', mag_earth_intensity);
% title({'After magnetometer calibration', '(Normalized to unit sphere)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;
legend('before calibration: measured data','before calibration: fitted ellipsoid','after calibration: calibrated data','after calibration: fitted ellipsoid','after calibration: sphere');

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
