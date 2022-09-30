function [M, gain, bias] = ellipse_to_circle(v)
% Generate a transformation matrix and a bias vector that changes arbitrary 
% ellipse into a circle centered at the origin.
%
% Inputs:
%   v = Coefficients of ellipse in conic form as vector [a, b, c, d, e, f].
%       ax^2 + bxy + cy^2 + dx + ey + f = 0
%
% 2022/09/21
 
a = v(1); b = v(2); c = v(3); d = v(4); e = v(5); f = v(6);

% Matrix form of ellipse
[a_, b_, xc, yc, phi] = ellipse_conic_to_parametric(v);

c = cos(phi);
s = sin(phi);
M = [c, s; -s, c]';
gain = [1/a_, 1/b_];
bias = [xc, yc]';
end
