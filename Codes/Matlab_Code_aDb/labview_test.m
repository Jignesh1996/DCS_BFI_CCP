Directory = "C:\Users\Jignesh\OneDrive - The University of Western Ontario\Research\Data\Pressure_modulation\Leena\Leena_Probe_pressure1\"



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


file=strcat(Directory, filename,num2str(i),'.csv');
temp=load(file);

adb = standalone_dcs(temp(1:4,:),Data_tau) ;

%     clear temp
