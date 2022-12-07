clear all;
close all;

Directory = "D:\Jignesh\MSc Western Uni\Research\Data\20220706\"; %Directory of the source file, 1 level up.

fprintf("Start....\n")


topLevelFolder = Directory; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'

cur_dir = topLevelFolder;
% Get a list of all files and folders in this folder.
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..



for i=1:size(subFolderNames,2)
    work_dir = strcat(Directory,subFolderNames(i));
    fold_name = subFolderNames(i);
    files_temp=dir(work_dir);
    filesFlags = ~[files_temp.isdir];
    subFiles = files_temp(filesFlags);
    
    for f=1:size(subFiles,1)
        temp=subFiles(f).name;
        ix=strfind(temp,'Data.mat');
        if ix>0
            filename=strcat(work_dir,"\Data.mat");
            break;
        end
    end
    
    new_dir = "D:\"; % Directory of destination, where the file should be copied
    dest = strcat(new_dir,fold_name,"_Data.mat");
    source = filename;
    copyfile(source,dest)
end
fprintf("Complete!!\n")
%--------------------------------------------------------------------------

