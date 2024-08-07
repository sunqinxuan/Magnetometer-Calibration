clc
clear
close all

%% 1st flight

% save_mat_name='model_uav_1_sensor3dmgx5.mat';

% data_mag13=readData_UDAU('.\data\uav\1\mag13\UDAU_data.txt');
% x_mag13=data_mag13(1688:9391,2);
% y_mag13=data_mag13(1688:9391,3);
% z_mag13=data_mag13(1688:9391,4);
% 
% data_3dmgx5=readData_csv('.\data\uav\1\sensor3dmgx5\Sensor3dmgx5_2024-07-04-10-02-21.csv');
% x_m=rmmissing(table2array(data_3dmgx5(315:10186,8)))*1e5;
% y_m=rmmissing(table2array(data_3dmgx5(315:10186,9)))*1e5;
% z_m=rmmissing(table2array(data_3dmgx5(315:10186,10)))*1e5;

%% 2nd flight

% save_mat_name='model_uav_2_sensor3dmgx5.mat';
% 
% data_mag13=readData_UDAU('.\data\uav\2\mag13\UDAU_data.txt');
% x_mag13=data_mag13(1:7754,2);
% y_mag13=data_mag13(1:7754,3);
% z_mag13=data_mag13(1:7754,4);
% 
% data_3dmgx5=readData_csv('.\data\uav\2\sensor3dmgx5\Sensor3dmgx5_2024-07-04-10-17-18.csv');
% x_m=rmmissing(table2array(data_3dmgx5(1:7500,8)))*1e5;
% y_m=rmmissing(table2array(data_3dmgx5(1:7500,9)))*1e5;
% z_m=rmmissing(table2array(data_3dmgx5(1:7500,10)))*1e5;

%% 3rd flight

save_mat_name='model_uav_3_sensor3dmgx5.mat';

data_mag13=readData_csv('.\data\uav\3\mag13\Mag13_2024-07-10-09-48-03.csv');
x_mag13=table2array(data_mag13(1107:10071,2));
y_mag13=table2array(data_mag13(1107:10071,3));
z_mag13=table2array(data_mag13(1107:10071,4));

data_3dmgx5=readData_csv('.\data\uav\3\sensor3dmgx5\Sensor3dmgx5_2024-07-10-09-48-03.csv');
x_m=rmmissing(table2array(data_3dmgx5(670:10700,8)))*1e5;
y_m=rmmissing(table2array(data_3dmgx5(670:10700,9)))*1e5;
z_m=rmmissing(table2array(data_3dmgx5(670:10700,10)))*1e5;

%%
mag_mag13=[];
for i=1:size(x_mag13,1)
    m=[x_mag13(i),y_mag13(i),z_mag13(i)];
%     sum=sum+norm(m);
    mag_mag13(i)=norm(m);
end
% mag_earth_intensity=sum/size(x_m,1); 
mag_earth_intensity=mean(mag_mag13);

mag_3dmgx5=[];
for i=1:size(x_m,1)
    m=[x_m(i),y_m(i),z_m(i)];
%     sum=sum+norm(m);
    mag_3dmgx5(i)=norm(m);
end
% mag_earth_intensity=sum/size(x_m,1); 
% mag_earth_intensity=mean(mag_3dmgx5);


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
