function y=readData_UDAU(filename)

a=load(filename);
N=size(a,1);

fid=fopen(filename,'r');
data=[];
y=[];
for i=1:N
    % diurnal data
    data{i,1}=fscanf(fid,'%f',1);
    data{i,2}=fscanf(fid,'%f',1);
    data{i,3}=fscanf(fid,'%f',1);
    tmp=fscanf(fid,'%f',1);
    tmp=fscanf(fid,'%f',1);
    tmp=fscanf(fid,'%f',1);
    tmp=fscanf(fid,'%s',1);
    % timestamp
    data{i,4}=fscanf(fid,'%s',1);
    str_time=split(data{i,4},':');
    hour=str2double(str_time{1});
    if hour<10
        hour=hour+12;
    end
    minute=str2double(str_time{2});
    second=str2double(str_time{3});
    time=hour*3600+minute*60+second;
    y(i,:)=[time,data{i,1},data{i,2},data{i,3}];
end
fclose(fid);
