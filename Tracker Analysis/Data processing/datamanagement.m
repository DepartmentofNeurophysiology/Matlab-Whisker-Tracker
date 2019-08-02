% DATA MANAGEMENT
% Scavange input folder (db.path_huma_clacked) for videos with manual
% annotations, copy available videos to a single directory
% (db.path_human_clicked_store) and save videos in common format
% (M(mouse)_R(session)_(trial). Save excel file with overview of all videos
%
% For all available videos, find 'general.mat' files and extract gapwidth
%%
clc
clear

% Manage data for human vs tracker comparison
db.path_human_clicked_store = 'F:\Studie\Stage Neurobiologie\Videos\Human clicked';
db.path_human_clicked_paper = 'F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
db.base_name = 'F:\Studie\Stage Neurobiologie\Videos\VideoDatabase';

db.name_manual_annotations = '_Annotations.mat';
db.name_video = '.dat';
db.name_metadata = '.mat';
db.name_tracker = '_Annotations_Tracker.mat';
db.name_compiled = '_compiled.mat';

db.unique = 'M%2d_R%02d_%02d';

% Index available videos
files = dir( fullfile(db.path_human_clicked_store,'**','*.dat'));

if exist(fullfile( db.base_name, 'overview.xlsx'), 'file')
    delete( fullfile( db.base_name, 'overview.xlsx' ));
end

copy_data_to_path = 1;



%% COPY DATA
fprintf('Checking for available videos, convert them to standard dataformat (%s)\n',db.unique)
fprintf('%d files found in data store:\n', size(files,1));

for i = 1:size(files,1)
    token = regexp(files(i).folder, '\', 'split');
    session = token{end};
    mouse = token{end-1};
    
    tid = regexp(files(i).name,'Data_(?<tid>\d+)','names');
    sid = regexp(session,'R(?<sid>\d+)','names');
    mid = regexp(mouse,'Mouse (?<mid>\d+)','names');
    
    fname = sprintf([db.unique '%s'],str2double(mid.mid), str2double(sid.sid),...
        str2double(tid.tid), db.name_video);
    
    fprintf('\t%s -',fname);
    
    vid_file_in = fullfile( files(i).folder, files(i).name);
    man_file_in = fullfile( files(i).folder, ...
        [files(i).name(1:end-4) db.name_manual_annotations]);
    meta_file_in= fullfile( files(i).folder, ...
        [files(i).name(1:end-4) db.name_metadata]);
    
          
    vid_file_out = sprintf([db.unique '%s'], str2double( mid.mid), ...
        str2double( sid.sid), str2double(tid.tid), db.name_video);
    man_file_out = sprintf([db.unique '%s'], str2double( mid.mid), ...
        str2double( sid.sid), str2double(tid.tid), db.name_manual_annotations);
    meta_file_out = sprintf([db.unique '%s'], str2double( mid.mid), ...
        str2double( sid.sid), str2double(tid.tid), db.name_metadata);
  
    
    
    if ~exist( fullfile( db.path_human_clicked_paper, vid_file_out ), 'file')
        fprintf(' does not exist');
        if copy_data_to_path
            fprintf(' - copying');
            copyfile(vid_file_in, fullfile( db.path_human_clicked_paper, vid_file_out));
            fprintf(' - done\n');
        else
            fprintf(' - not copying\n');
        end      
        
    else
        fprintf(' allready exists\n');
    end
    
     if ~exist( fullfile( db.path_human_clicked_paper, man_file_out ), 'file')
         if copy_data_to_path
             copyfile(man_file_in, fullfile( db.path_human_clicked_paper, man_file_out));
         end
     end
     
     if ~exist( fullfile( db.path_human_clicked_paper , meta_file_out ), 'file')
         if copy_data_to_path
             copyfile(meta_file_in, fullfile( db.path_human_clicked_paper, meta_file_out));
         end
     end
     
    
    outname{i,1} = fname;  
end




%% check for general .mat files
% copy gapwidth into an array and store in output folder
fprintf('\n\n checking for ''general.mat'' files:\n');

addpath(genpath('F:\Studie\Stage Neurobiologie\Controller')) % path to 'Controller path'


gapwidth = cell( size(files,1),2);

files = dir( fullfile( db.path_human_clicked_paper, '*.dat'));

for i = 1:size(files,1)
    
    fprintf('\t%s',files(i).name);
    gapwidth{i,1} = files(i).name;
    
    tokens = regexp( files(i).name, 'M(?<mouse>\d+)_R(?<session>\d+)_(?<trial>\d+).dat','names');
    gen_path = fullfile(db.path_human_clicked_store, ...
        sprintf('Mouse %s',tokens.mouse),sprintf('R%d',str2double(tokens.session)));
    
    if exist( fullfile( gen_path, 'General.mat'), 'file' )
        gapdata = load( fullfile( gen_path, 'General.mat'));
        trialnr = str2double( tokens.trial);
        gw = gapdata.CGSave.Paradigm.Trials(trialnr).Distance;
        gapwidth{i,2} = gw;   
        
        fprintf(': gapwidth - %2.2f\n', gw);
    else
        gapwidth{i,2} = NaN;
        fprintf(': gapwidth - NaN\n');
    end
    
    
end

fprintf('saving file...')
save( fullfile( db.path_human_clicked_paper, 'gapwidths'),'gapwidth')
fprintf('done!\n');



%% 


% Manage data for Single vs Multi whisker
db.path_single_vs_multi = 'F:\Studie\Stage Neurobiologie\Videos\Whisker Deprivation mouse34';
db.path_single_vs_multi_paper = 'F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi';
db.name_video = '.dat';
db.name_metadata = '.mat';

files = dir( fullfile(db.path_single_vs_multi,'*','*.dat'));

mid.mid = '34';
fprintf('Checking for available videos, convert them to standard dataformat (%s)\n',db.unique)
fprintf('%d files found in data store:\n', size(files,1));


if ~exist(db.path_single_vs_multi_paper, 'file')
    mkdir(db.path_single_vs_multi_paper)
end

for i  = 1:size(files, 1)
    token = regexp(files(i).folder, '\', 'split');
    session = token{end};
     tid = regexp(files(i).name,'Data_(?<tid>\d+)','names');
    sid = regexp(session,'R(?<sid>\d+)','names');
    fname = sprintf([db.unique '%s'],str2double(mid.mid), str2double(sid.sid),...
        str2double(tid.tid), db.name_video);
   fprintf('\t%s -',fname);
    
    
   vid_file_in = fullfile( files(i).folder, files(i).name);
    meta_file_in= fullfile( files(i).folder, ...
        [files(i).name(1:end-4) db.name_metadata]);
    
     vid_file_out = sprintf([db.unique '%s'], str2double( mid.mid), ...
        str2double( sid.sid), str2double(tid.tid), db.name_video);
    meta_file_out = sprintf([db.unique '%s'], str2double( mid.mid), ...
        str2double( sid.sid), str2double(tid.tid), db.name_metadata);
    
    
      if ~exist( fullfile( db.path_single_vs_multi_paper, vid_file_out ), 'file')
        fprintf(' does not exist');
        if copy_data_to_path
            fprintf(' - copying');
            copyfile(vid_file_in, fullfile( db.path_single_vs_multi_paper, vid_file_out));
            fprintf(' - done\n');
        else
            fprintf(' - not copying\n');
        end      
        
    else
        fprintf(' allready exists\n');
    end
       
     if ~exist( fullfile( db.path_single_vs_multi_paper , meta_file_out ), 'file')
         if copy_data_to_path
             copyfile(meta_file_in, fullfile( db.path_single_vs_multi_paper, meta_file_out));
         end
     end
     
     outname{end+1,1} = fname;
     
end

fprintf('saving excel file...');
xlswrite( fullfile( db.base_name, 'overview.xlsx' ), outname);
fprintf(' done\n');

