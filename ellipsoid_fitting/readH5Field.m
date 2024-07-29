function data=readH5Field(file_name, line_number, field_name)

data_line = h5read(file_name,'/line');
i1 = find(data_line==line_number, 1 );
i2 = find(data_line==line_number, 1, 'last' );

data = h5read(file_name,field_name);
data = data(i1:i2,:);
