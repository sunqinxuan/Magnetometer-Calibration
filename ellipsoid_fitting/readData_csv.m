function y=readData_csv(name)

opts = detectImportOptions(name);
opts.DataLines = [2,Inf];
y = readtable(name, opts);
% N=size(y,1);
% y=y(1:N/2,:);