clear all
close all
%--------------------------------------------------------------------------

% Folder name  - provide folder name for which you want to convert the data
% Folder='28.4.21-DCS only test';

Folder='20211116-12';
%--------------------------------------------------------------------------

Directory=strcat(pwd, '\',Folder,'\');

files_temp=dir(Directory);
filesFlags = ~[files_temp.isdir];
subFiles = files_temp(filesFlags);

for f=1:size(subFiles,1)
    temp=subFiles(f).name;
    ix=strfind(temp,'_');
    if ix>0
        filename=temp(1:ix(1)-1);
    end
end

subFiles2=0;
for f=1:size(subFiles,1)
    temp=subFiles(f).name;
    ix=strfind(temp,'csv');
    if ix>0
        subFiles2=subFiles2+1;
    end
end

% Save as a mat -----------------------------------------------------------

file=strcat(Directory, filename,'_tau.csv');
Data_tau=load (file);

for f=1:1:subFiles2-1
    file=strcat(Directory, filename,num2str(f),'.csv');
    temp=load (file);
    Data(f,:,:)=temp(1:4,:);
end

clearvars -except Directory Data_tau Data
filename=strcat(Directory,'Data.mat'); 
save(filename,'Data','Data_tau');

%--------------------------------------------------------------------------

disp('Data export done')

