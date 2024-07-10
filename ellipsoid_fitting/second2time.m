function [h,m,s]=second2time(t)

h=floor(t/3600);
m=mod(t,3600)/60;
s=(m-floor(m))*60;