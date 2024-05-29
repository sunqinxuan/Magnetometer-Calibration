function y=readH5File(file_name, field_name, i1,i2)

y = h5read(file_name,field_name);
y = y(i1:i2,:);