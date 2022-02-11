clear all
close all
%--------------------------------------------------------------------------



%--------------------------------------------------------------------------

Directory="D:\Jignesh\MSc Western Uni\Research MSc\Codes\Western-MSc\Data\DCS\Farah\Farah Data 2_10_2022\DCS\20220210 - 4\";

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
    Data(f,:,:)=load (file);
end

clearvars -except Directory Data_tau Data
filename=strcat(Directory,'Data.mat'); 
save(filename,'Data','Data_tau');

%--------------------------------------------------------------------------


