clear all;
close all;
%--------------------------------------------------------------------------

% Folder name  - provide folder name for which you want to convert the data
Folder='20220309 - 8';

%--------------------------------------------------------------------------

Directory = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\TNQT Pulsatility study\06062022\Pressure test_20220606\";

fprintf("Start....\n")

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

% Save as a mat -----------------------------------------------------------

file=strcat(Directory, filename,'_tau.csv');
Data_tau=load (file);

for f=1:1:size(subFiles,1)-1
    file=strcat(Directory, filename,num2str(f),'.csv');
    temp=load(file);
    Data(f,:,:)=temp(1:4,:);
    clear temp
end

clearvars -except Directory Data_tau Data
filename=strcat(Directory,'Data.mat'); 
save(filename,'Data','Data_tau');
fprintf("Complete!!\n")
%--------------------------------------------------------------------------


