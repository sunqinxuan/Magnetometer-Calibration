function [v] = ellipse_fit(x, y)
% Output:
%     v = [a, b, c, d, e, f]'
%       ax^2 + bxy + cy^2 + dx + ey + f = 0  
%
% 2022/09/1
  
% Design matrix 
D1 = [x.*x, x.*y, y.*y];    % Quadratic part
D2 = [x, y, ones(size(x))]; % Linear part

% Scatter matrix0 
S1 = D1' * D1; % Quatratic part
S2 = D1' * D2; % Combined part
S3 = D2' * D2; % Linear part
  
% Reduced scatter matrix, M
T = - inv(S3) * S2';
M = S1 + S2 * T;
M = [M(3, :)./2; - M(2, :); M(1, :)./2];

% Eolve eigensystem
[evec, eval] = eig(M);

% Eigenvector for min. pos. eigenvalue
cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :).*evec(2, :); % evaluate a’Ca
a1 = evec(:, find(cond > 0));

% Coefficients of ellipse
v = [a1; T * a1];
end 
