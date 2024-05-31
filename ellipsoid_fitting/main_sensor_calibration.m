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

addpath('.\data')
addpath('..\m_IGRF')

% data_original_filename = 'Flt1002_train.h5';
% time = datenum([2020 6 20]); 
% lines={1002.02,1002.20};

data_original_filename = 'Flt1006_train.h5';
time = datenum([2020 7 6]); 
lines={1006.04,1006.04,1006.08};

tt=[];
x_m=[];
y_m=[];
z_m=[];
mag_earth=[];
for i=1:size(lines,2)
    [tt_tmp,x_m_tmp,y_m_tmp,z_m_tmp,mag_earth_tmp]=readH5File(data_original_filename, lines{i}, time);
    tt=[tt;tt_tmp];
    x_m=[x_m;x_m_tmp];
    y_m=[y_m;y_m_tmp];
    z_m=[z_m;z_m_tmp];
    mag_earth=[mag_earth;mag_earth_tmp];
end

% line_number = 1002.20; 
% [tt_2,x_m_2,y_m_2,z_m_2,mag_earth_2]=readH5File(data_original_filename, line_number, time);

% tt=[tt_1;tt_2];
% x_m=[x_m_1;x_m_2];
% y_m=[y_m_1;y_m_2];
% z_m=[z_m_1;z_m_2];
% mag_earth=[mag_earth_1;mag_earth_2];

% figure;
% plot(tt,mag_earth,'k');hold on;
% plot(tt,mag_1_uc,'r');hold on;
% plot(tt,mag_1_c,'g');hold on;
% plot(tt,mag_1_dc,'b');hold on;

%%

mag_earth_intensity=mean(mag_earth);

% Import raw magnetometer readings
% file = 'sensor_data.txt'; 

% Import sensor readings
% raw = importdata(file);
% x_m = raw(:,1); 
% y_m = raw(:,2); 
% z_m = raw(:,3);
% x_m=flux_b_x;
% y_m=flux_b_y;
% z_m=flux_b_z;

% Ellipsoid fit
% ax^2 + by^2 + cz^2 + 2fyz + 2gxz + 2hxy + 2px + 2qy + 2rz + d = 0
% v = [a, b, c, f, g, h, p, q, r, d]' (in the paper k = -d)
% M = [a h g; h b f; g f c]
% u = [p, q, r]'
v = ellipsoid_fit(x_m, y_m, z_m);

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
matrix = gain*rotation*mag_earth_intensity;

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

cell_str=strsplit(data_original_filename,'_');
save_mat_name=['model_',cell_str{1,1},'.mat'];
save(save_mat_name,'matrix','offset','v');

residual_h_m_mean=mean(residual_h_m);
residual_h_hat_mean=mean(residual_h_hat);

% Visualization %
% Sensor readings and ellipoid fit
scatter3(x_m, y_m, z_m, 'fill', 'MarkerFaceColor', 'red'); hold on; 
plot_ellipsoid(v); 
% title({'Before magnetometer calibration', '(Ellipsoid fit)'});
% xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
% axis equal;

% After calibrations
% figure;
scatter3(x_hat, y_hat, z_hat, 'fill', 'MarkerFaceColor', 'blue'); hold on;
plot_sphere([0,0,0]', mag_earth_intensity);
% title({'After magnetometer calibration', '(Normalized to unit sphere)'});
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal;

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
