clc; clear; close all; format longg;

cd 'D:\GPBL SIT\Project\Final Project\Parameter Japan';
d = dir('D:\GPBL SIT\Project\Final Project\Parameter Japan');
file_name = {d.name};
allRt = [];
for i = 3:41
    
    file = xlsread(file_name{i});
    allRt(i-2,:) = smooth(file(:,8)',5);
    
end

cd 'D:\GPBL SIT\Project\Final Project';
% serial_date = datestr(now, 'dd-mm-yy_HH-MM-SS');
filename_save = strjoin({'allRt','.xlsx'});
filename_save = filename_save(~isspace(filename_save));
xlswrite(filename_save,allRt,'sheet1');

data = xlsread('allRt.xlsx');

J = isoutlier(data, 'mean');
total = sum(J(:)~=0);
i=1;

while total>=1
    for j=1:294
        replaced = J(:,j);
        replacer = ~replaced;
        data(replaced, j) = mean(data(replacer,j));
    end
    J = isoutlier(data, 'mean');
    total = sum(J(:)~=0);
    i=i+1;
    disp(i);
end

writematrix(data, 'D:\GPBL SIT\Project\Final Project\RtClean.xlsx');