function [x_hat,y_hat,z_hat]=applyModel(x_m,y_m,z_m,mag_earth_intensity,matrix,offset)

residual_h_m=zeros(size(x_m));
residual_h_hat=zeros(size(x_m));

x_hat = zeros(length(x_m),1); 
y_hat = zeros(length(x_m),1); 
z_hat = zeros(length(x_m),1);

for i_iters = 1:length(x_m)
    % Sensor data
    h_m = [x_m(i_iters); y_m(i_iters); z_m(i_iters)]; 
    % Calibration, Eqn(11)
    h_hat = matrix*(h_m - offset);
    % Calibrated values
    x_hat(i_iters) = h_hat(1);
    y_hat(i_iters) = h_hat(2);
    z_hat(i_iters) = h_hat(3);
    % residuals
    residual_h_m(i_iters)=abs(norm(h_m)-mag_earth_intensity);
    residual_h_hat(i_iters)=abs(norm(h_hat)-mag_earth_intensity);
end

residual_h_m_mean=mean(residual_h_m);
residual_h_hat_mean=mean(residual_h_hat);

fprintf('\nresidual_h_m_mean ='); disp(residual_h_m_mean);
fprintf('\nresidual_h_hat_mean ='); disp(residual_h_hat_mean);