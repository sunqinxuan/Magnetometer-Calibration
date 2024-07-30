function plotResults(x_m,y_m,z_m,x_hat,y_hat,z_hat,mag_earth_intensity)

figure;
% Visualization %
% Sensor readings and ellipoid fit
scatter3(x_m, y_m, z_m, 'fill', 'MarkerFaceColor', 'red'); hold on; 
v = ellipsoid_fit(x_m, y_m, z_m);
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
legend('before calibration: measured data','before calibration: fitted ellipsoid','after calibration: calibrated data','after calibration: fitted ellipsoid','sphere');

