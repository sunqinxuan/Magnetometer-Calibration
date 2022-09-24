function [a_, b_, xc, yc, phi] = ellipse_conic_to_parametric(v)
% Convert the coefficients of ellipse in the conic form to the parameters
% of the ellipse.
%
% Inputs:
%   v = Coefficients of ellipse in conic form as vector [a, b, c, d, e, f].
%       ax^2 + bxy + cy^2 + dx + ey + f = 0
% 
% Outputs:
%   a_  = Semi-major axis 
%   b_  = Semi-minor axis 
%   xc  = X-coordinate of center 
%   yc  = Y-coordinate of center 
%   phi = Rotation of the ellipse wrt. X-axis (wrt. semi-major axis) 
%
% Reference:
%  Loan - Using the Ellipse to Fit and Enclose Data Points
%
% 2022/09/21

a = v(1); b = v(2); c = v(3); d = v(4); e = v(5); f = v(6);
M0 = [f, d/2, e/2; d/2, a, b/2; e/2, b/2, c];
M = M0(2:3, 2:3);
lambda = eig(M);

if(~(abs(lambda(1)-a) <= abs(lambda(1)-c)))
  temp = lambda(2);
  lambda(2) = lambda(1);
  lambda(1) = temp;
end

a_ = sqrt(-det(M0)/(det(M)*lambda(1)));
b_ = sqrt(-det(M0)/(det(M)*lambda(2)));
temp = 4*a*c - b^2;
xc = (b*e - 2*c*d)/temp;
yc = (b*d - 2*a*e)/temp;
phi = acot((a-c)/b)/2;

if(b_ > a_)
  temp = b_;
  b_ = a_;
  a_ = temp;
  phi = pi/2 + phi;
end
end
