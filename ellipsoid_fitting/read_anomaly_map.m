function [mag_anomaly,indices_x,indices_y]=read_anomaly_map(anomaly_map_file,utm_x,utm_y)


anomaly_map = h5read(anomaly_map_file,'/map');
anomaly_map_xx=h5read(anomaly_map_file,'/xx');
anomaly_map_yy=h5read(anomaly_map_file,'/yy');

% img_traj=anomaly_map;
% max_value=max(max(anomaly_map));

x=utm_x(1);
y=utm_y(1);
idx_x=0;
idx_y=0;
mag_anomaly=zeros(size(utm_x));
indices_x=zeros(size(utm_x));
indices_y=zeros(size(utm_x));
for i=1:size(anomaly_map_xx,1)
    if abs(x-anomaly_map_xx(i))<100
        for j=1:size(anomaly_map_yy,1)
            if abs(y-anomaly_map_yy(j))<100
%                 img_traj(j,i)=max_value;
                indices_x(1)=i;
                indices_y(1)=j;
                mag_anomaly(1)=anomaly_map(j,i);
                idx_x=i;
                idx_y=j;
%                 fprintf('x=');disp(x);
%                 fprintf('y=');disp(y);
%                 fprintf('anomaly=');disp(mag_anomaly(1));
                break;
            end
        end
    end
end

if idx_x==0 || idx_y==0
    fprintf('cannot find the utm location!');
else
    for i=2:size(utm_x,1)
        x=utm_x(i);
        y=utm_y(i);
        for ii=idx_x-20:idx_x+20
            for jj=idx_y-20:idx_y+20
                if abs(x-anomaly_map_xx(ii))<100 && abs(y-anomaly_map_yy(jj))<100
%                     img_traj(jj,ii)=max_value;
                    indices_x(i)=ii;
                    indices_y(i)=jj;
                    mag_anomaly(i)=anomaly_map(jj,ii);
                    idx_x=ii;
                    idx_y=jj;
%                     disp(i);
%                     fprintf('x=');disp(x);
%                     fprintf('y=');disp(y);
%                     fprintf('anomaly=');disp(mag_anomaly(i));
                    break;
                end
            end
        end
    end
end

