function [x_m,y_m,z_m,mag_earth_intensity]=loadMITData(data_original_filename, lines, time)

tt=[];
x_m=[];
y_m=[];
z_m=[];
mag_earth=[];
map_idx_x=[];
map_idx_y=[];
for i=1:size(lines,2)
    [tt_tmp,x_m_tmp,y_m_tmp,z_m_tmp,mag_earth_tmp,map_idx_x_tmp,map_idx_y_tmp]=readH5File(data_original_filename, lines{i}, time);
    tt=[tt;tt_tmp];
    x_m=[x_m;x_m_tmp];
    y_m=[y_m;y_m_tmp];
    z_m=[z_m;z_m_tmp];
    mag_earth=[mag_earth;mag_earth_tmp];
    map_idx_x=[map_idx_x;map_idx_x_tmp];
    map_idx_y=[map_idx_y;map_idx_y_tmp];
end

anomaly_map= h5read('Canada_MAG_RES_200m.hdf5','/map');
max_value=max(max(anomaly_map));
img_traj=anomaly_map;
for i=1:size(map_idx_x,1)
    img_traj(map_idx_y(i),map_idx_x(i))=max_value;
    img_traj(map_idx_y(i)-1,map_idx_x(i))=max_value;
    img_traj(map_idx_y(i),map_idx_x(i)-1)=max_value;
    img_traj(map_idx_y(i)+1,map_idx_x(i))=max_value;
    img_traj(map_idx_y(i),map_idx_x(i)+1)=max_value;
end
% figure;
% imshow(img_traj,[]);

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

mag_earth_intensity=mean(mag_earth);
