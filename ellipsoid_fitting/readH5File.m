function [tt,x_m,y_m,z_m,mag_earth,map_idx_x,map_idx_y]=readH5File(file_name, line_number, time)

% data_info = h5info(data_original_filename);
data_line = h5read(file_name,'/line');
i1 = find(data_line==line_number, 1 );
i2 = find(data_line==line_number, 1, 'last' );

tt=readDataField(file_name,'/tt',i1,i2);

x_m=readDataField(file_name,'/flux_c_x',i1,i2);
y_m=readDataField(file_name,'/flux_c_y',i1,i2);
z_m=readDataField(file_name,'/flux_c_z',i1,i2);

% flux_b_t=readDataField(data_original_filename,'/flux_b_t',i1,i2);
% mag_1_uc=readDataField(data_original_filename,'/mag_1_uc',i1,i2);
% mag_1_c=readDataField(data_original_filename,'/mag_1_c',i1,i2);
% mag_1_dc=readDataField(data_original_filename,'/mag_1_dc',i1,i2);
% mag_1_igrf=readDataField(data_original_filename,'/mag_1_igrf',i1,i2);

% figure;
% plot(tt,mag_1_uc,'r');hold on;
% plot(tt,mag_1_c,'g');hold on;
% plot(tt,mag_1_dc,'b');hold on;
% % plot(tt,mag_1_igrf,'k');hold on;

% w_x=readH5File(data_original_filename,'/roll_rate',i1,i2);
% w_y=readH5File(data_original_filename,'/pitch_rate',i1,i2);
% w_z=readH5File(data_original_filename,'/yaw_rate',i1,i2);

utm_x=readDataField(file_name,'/utm_x',i1,i2);
utm_y=readDataField(file_name,'/utm_y',i1,i2);

[mag_anomaly,map_idx_x,map_idx_y]=read_anomaly_map('Canada_MAG_RES_200m.hdf5',utm_x,utm_y);

mag_diurnal=readDataField(file_name,'/diurnal',i1,i2);

baro=readDataField(file_name,'/baro',i1,i2);
lat=readDataField(file_name,'/lat',i1,i2);
lon=readDataField(file_name,'/lon',i1,i2);

mag_earth=zeros(size(lat));
gh = loadigrfcoefs(time);
for i=1:size(lat,1)
    latitude=lat(i);
    longitude=lon(i);
    altitude=baro(i)*1e-3;
    [Bx, By, Bz] = igrf(gh, latitude, longitude, altitude, 'geod');
    mag_earth(i)=norm([Bx,By,Bz])+mag_anomaly(i)+mag_diurnal(i);
end


%%
function y=readDataField(file_name, field_name, i1,i2)

y = h5read(file_name,field_name);
y = y(i1:i2,:);