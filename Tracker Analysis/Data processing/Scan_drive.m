% first run conda.bat, cd ~\Code\Python\Drive_API, 
% activate DRIVE_API and python quickstart.py.
%
% then compy the 'Drive_files.txt to this working directory
clear
clc

f = fopen('Drive_files.txt','r');



Files = {};
while size(Files,1) == 0 ||  ischar(Files{end})
    Files{end+1} = fgetl(f);
end
fclose(f);


Filepath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
if ~exist(fullfile(Filepath,'to_copy'),'file')
    mkdir(fullfile(Filepath,'to_copy'))
end
oldfiles = dir(fullfile(Filepath,'to_copy'));
for i = 1:size(oldfiles,1)
    if ~oldfiles(i).isdir
        delete(fullfile(oldfiles(i).folder, oldfiles(i).name))
    end
end

copy_files = 1;


% Check if all dat files are uploaded
dat_files = dir(fullfile(Filepath, '*.dat'));
Not_found = {};
for i = 1:size(dat_files,1)
   if isempty(find(strcmp(dat_files(i).name, Files)))
       Not_found{end+1} = dat_files(i).name;       
   end    
end

if isempty(Not_found)
    fprintf('All .dat files are uploaded!\n')
else
    fprintf('Not all .dat files are uploaded:\n')
    for i = 1:size(Not_found, 2)
        fprintf('\t%s\n',Not_found{i})
        if copy_files
            copyfile(fullfile(Filepath,Not_found{i}), fullfile(Filepath,'to_copy',Not_found{i}))
        end
    end
end


% Check if all metafiles are uploaded
meta_files = dir(fullfile(Filepath, '*.mat'));
Not_found = {};
for i =1 :size(meta_files, 1)
    if length(meta_files(i).name) == 14 
        if isempty(find(strcmp(meta_files(i).name, Files)))
            Not_found{end+1} = meta_files(i).name;            
        end
    end
end

if isempty(Not_found)
    fprintf('All metafiles files are uploaded!\n')
else
    fprintf('Not all metafiles files are uploaded:\n')
    for i = 1:size(Not_found, 2)
        fprintf('\t%s\n',Not_found{i})
        
        if copy_files
            copyfile(fullfile(Filepath,Not_found{i}), fullfile(Filepath,'to_copy',Not_found{i}))
        end
    end
end

% Check if all manual Annotations are uploaded
ann_files = dir(fullfile(Filepath, '*Annotations.mat'));
Not_found = {};
for i = 1:size(ann_files, 1)
    if isempty(find(strcmp(ann_files(i).name, Files)))
        Not_found{end+1} = ann_files(i).name;
    end
end

if isempty(Not_found)
    fprintf('All MANUAL annotation files are uploaded!\n')
else
    fprintf('Not all MANUAL files are uploaded:\n')
    for i = 1:size(Not_found, 2)
        fprintf('\t%s\n',Not_found{i})
        if copy_files
            copyfile(fullfile(Filepath,Not_found{i}), fullfile(Filepath,'to_copy',Not_found{i}))
        end
    end
end


% Check if all Tracker Annotations are uploaded
ann_files = dir(fullfile(Filepath, '*Annotations_Tracker.mat'));
Not_found = {};
for i = 1:size(ann_files, 1)
    if isempty(find(strcmp(ann_files(i).name, Files)))
        Not_found{end+1} = ann_files(i).name;
    end
end

if isempty(Not_found)
    fprintf('All TRACKER annotation files are uploaded!\n')
else
    fprintf('Not all TRACKER files are uploaded:\n')
    for i = 1:size(Not_found, 2)
        fprintf('\t%s\n',Not_found{i})
        if copy_files
            copyfile(fullfile(Filepath,Not_found{i}), fullfile(Filepath,'to_copy',Not_found{i}))
        end
    end
end


% Check if all compiled files are uploaded
comp_files = dir(fullfile(Filepath, '*compiled.mat'));
Not_found = {};
for i = 1:size(comp_files, 1)
    if isempty(find(strcmp(comp_files(i).name, Files)))
        Not_found{end+1} = comp_files(i).name;
    end
end

if isempty(Not_found)
    fprintf('All COMPILED annotation files are uploaded!\n')
else
    fprintf('Not all COMPILED files are uploaded:\n')
    for i = 1:size(Not_found, 2)
        fprintf('\t%s\n',Not_found{i})
        if copy_files
            copyfile(fullfile(Filepath,Not_found{i}), fullfile(Filepath,'to_copy',Not_found{i}))
        end
    end
end



    

