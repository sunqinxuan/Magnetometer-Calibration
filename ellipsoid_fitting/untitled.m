close all;
clc;clear;

addpath('.\data')

% anomaly_map_file='Canada_MAG_RES_200m.hdf5';
% anomaly_map = h5read(anomaly_map_file,'/map');
% anomaly_map_xx=h5read(anomaly_map_file,'/xx');
% anomaly_map_yy=h5read(anomaly_map_file,'/yy');
% 
% imshow(anomaly_map,[]);
% 
% map=load("map.txt");
% figure;
% imshow(map',[]);

load('model_Flt1002.mat');

data_original_filename = 'Flt1002_train.h5';
time = datenum([2020 6 20]); 
lines={1002.02,1002.20};
[x_m,y_m,z_m,mag_earth_intensity]=loadMITData(data_original_filename, lines, time);


% Unpack ellipsoid coefficients
a = v(1); b = v(2); c = v(3);
f = v(4); g = v(5); h = v(6); 
p = v(7); q = v(8); r = v(9); 
d = v(10); 

% Coordinate frame transformation i.e diagonalize M 
M =[a, h, g; h, b, f; g, f, c]; % Original ellipsoid matrix 
u = [p, q, r]';
k = d;

fprintf('\n\nM =\n'); disp(M);
fprintf('\nu = '); disp(u);
fprintf('\nk = '); disp(k);

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

fprintf('\n\nmatrix =\n'); disp(matrix);
fprintf('\noffset = '); disp(offset);