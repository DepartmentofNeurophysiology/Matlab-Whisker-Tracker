% This path contains functions used to analyse and visualize
% MWT tracker results. Example usage:

% Before analyzing the results (as stored in *_Annotations_Tracker.mat
% files), the data is first preprocessed using compiledata:

compiledata('file', 'full_filename', 'data', {'Tracker'})

% or do a whole directory at once:

compiledata('path', 'full_path_name', 'data', {'Tracker'})

% The function saves a *_compiled.mat file in the directory with 
% processed data

% the results can be visualized and saved using PRINT_VIDEO:

PRINT_VIDEO('dPath', 'path_to_videos', 'Name' , 'filename_compiled.mat',...
    'dNose',1,'dTclean',1,'dTtouch',1,'dExp',1)


% The remaining functions are used to generate figures from the data
% Figures/Figure_tracker_processing
% Tracker Performanc/Figure_Parameter_Validation
% Single vs Multi/Figure_behaviour_1
%
% and support those scripts

