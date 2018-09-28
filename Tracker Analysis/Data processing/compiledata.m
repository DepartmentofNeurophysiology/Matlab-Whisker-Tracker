function compiledata(varargin)
%COMPILEDATA
%
% Adds all available tracking data into an 'Annotations' struct and saves
% output into the videopath.
%
%
% OPTIONAL INPUT NAME-VALUE PAIRS:
%
%   'path' - [string] datapath to scan for videofiles
%   'overwrite' - [0 or 1] overwrite previously saved Annotations 
%                  (default: 0)   
%   'file' - [string] a single filename for compilation, leave empty to
%   scan a directory (default: '')
%   'fileindex' - [array] index for files to compile (ie. [1,3,5] will
%           compile files [1, 3, 5] as found in the scanning path,
%           (default: all)
%   'data' - [cell] multiple entries to add for compilation
%   ('Tracker','Manual','Janelia'), (ie {'Tracker','Manual'} will add only
%   MWT and Manual data, (default: 'Tracker;)
%
%
% OUTPUT:
% 
% Annotations struct:
%   .Output (raw tracker data)
%   .Settings (settings used in tracking)
%   .Tracker (processed tracking data, with new fields:)
%       .Parameters - parametrized (raw) tracking data
%       .gapinfo    - dimensions of gap
%       .dist_nose_target - distance nose to target
%       .exploring  - (not used) flag marking frames with mouse exploring
%       .Traces_clean- processed traces (noise filtered, spline fit)
%       .Parameters_clean - parametrized (clean) tracking data
%       .Touch      - Detected touch
%
%
% EXAMPLE USAGE:
%
%   single-file, tracker only compilation with overwrite:
%       compiledata('file','C:\Videos\example.dat','data',{'Tracker'}, 'overwrite', 1)
%
%   multi-file, tracker/janelia L
%       compiledata('path','C:\Videos\','data',{'Tracker',Janelia'})
%
%
%
%
%% Parse input

p = inputParser;
p.CaseSensitive = 0;
defaultPath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
defaultOverwrite = 0;
defaultFile = '';
defaultTrackingData = {'Tracker'};
defaultFileIndex = 'all';


addParameter(p,'Path', defaultPath);
addParameter(p,'Overwrite',defaultOverwrite);
addParameter(p,'File',defaultFile);
addParameter(p,'Data',defaultTrackingData);
addParameter(p,'FileIndex', defaultFileIndex);


parse(p,varargin{:});


if ~isempty(p.Results.File)
    idx = find(p.Results.File == '\',1,'last');
    Files(1).folder = p.Results.File(1:idx-1);
    Files(1).name = p.Results.File(idx+1:end);
    
    
    if isdir(p.Results.File)
        fprintf('Specify filename, a directory was provided: \n%s\n', p.Results.File)
        return
    elseif ~exist( fullfile( Files(1).folder, Files(1).name), 'file')
        fprintf('File ''%s'' does not exis\nt', fullfile( Files(1).folder, Files(1).name))
        return
    end
    
    FilesToAdd = 1;
    
    
else
    Files = dir(fullfile(p.Results.Path,'*.mat'));
    
    if isstr(p.Results.FileIndex)
        switch(p.Results.FileIndex)
            case 'all'
                FilesToAdd = 1:size(Files,1);
            otherwise
                fprintf('Fileindex string not recognized\n')
                return
        end
    else
        FilesToAdd = p.Results.FileIndex;
    end
    
   
end






%% Process data

for file_index = 1:length(FilesToAdd)
    PathName = Files( FilesToAdd( file_index )).folder;
    BaseName = Files( FilesToAdd( file_index )).name(1:end-4);
    
    
    fprintf('(%2d/%2d) Saving tracking data ',file_index,length(FilesToAdd))
    Annotations = [];
    
    % Meta data
    meta_file = fullfile( PathName, [BaseName '.mat']);
    
    % Tracker data
    tracker_file = fullfile( PathName, [BaseName '_Annotations_Tracker.mat']);
    
    % Manual data
    manual_file = fullfile( PathName, [BaseName '_Annotations.mat']);
       
    % Output file
    save_file = fullfile( PathName, [BaseName '_compiled.mat']);
    disp_name = [BaseName '_compiled.mat'];
    
    fprintf('into %s', disp_name);
    
    if exist(save_file, 'file') && ~p.Results.Overwrite
        fprintf(' - Allready compiled\n');
        continue
    end
    
    if exist(meta_file, 'file') & 0
        MetaData = load(meta_file);
        MetaData.Data.TimeMS = (MetaData.Data.Time(:,1) - MetaData.Data.Time(1,1))*1000;        
    else
        fprintf(' - ERROR: metadata not found\n');
        %continue
    end
    
    fprintf('\n\t Structs added')    
    if any(strcmp('Tracker', p.Results.Data)) && exist(tracker_file, 'file')        
        tracker_structs = who('-file', tracker_file);
        if any(strcmp('Output',tracker_structs)) & any(strcmp('Settings',tracker_structs))
            
            load(tracker_file,'Output','Settings') 
            Settings.Video(1) = PathName(1);
            
            % Copy tracker output to new Tracker struct
            %Tracker.MetaData = MetaData.Data;
            Tracker.Objects = Output.Objects;
            Tracker.Direction = Output.Direction;
            Tracker.Nose = Output.Nose;
            Tracker.Headvec = Output.AngleVector; % Angle of head w.r.t. frame
            Tracker.Origins = Output.Origins; % Tracking seeds
            Tracker.Traces = Output.Traces;
            Tracker.Parameters = getParams(Tracker, 'raw');
            
            % Find specifications of gap
            Tracker.gapinfo = detectGap(Tracker.Objects);
            [Tracker.dist_nose_target, Tracker.exploring] = getDistTarget(Tracker);
            
            % Parameterize Traces
            % Tracker.Parameters = getParams(Tracker,'raw');
            
            % Clean Traces
            Tracker.Traces_clean = CleanTraces(Tracker, 1);
            Tracker.Parameters_clean = getParams(Tracker, 'clean');
            
            % Detect Touch
            %[Tracker.Touch, ~] = DetectTouch(Tracker, Settings);
            
            Annotations.Output = Output;
            Annotations.Settings = Settings;
            Annotations.Tracker = Tracker;
            fprintf(' - Tracker');
            
        else
            fprintf(' - Tracker file incomplete')
        end
        
    elseif any(strcmp('Tracker', p.Results.Data))
        fprintf(' - Tracker not found');
    end
    
    
    
    if any(strcmp('Manual', p.Results.Data)) && exist(manual_file ,'file')
        manual_vars = who('-file', manual_file);
        
        if any(strcmp('CurvesByFrame', manual_vars))
            load(manual_file, 'CurvesByFrame')

            % Store unprocessed manual notations
            Manual.RawNotations = CurvesByFrame;

            % Convert manual notations format to tracker format
            converted_data = ConvertAnnotations(Manual.RawNotations);
            Manual.Objects = Output.Objects; % Manual tracking data does not contain object detection
            Manual.Nose = Output.Nose; % Manual tracking data does not contain nose
            Manual.Headvec = Output.AngleVector; % Manual tracking data does not contain headangle
            Manual.Traces = converted_data.Traces;
            Manual.Parameters = getParams(Manual, 'raw');
            Manual.Labels = converted_data.Labels.Full;
            Manual.Label_names = converted_data.Labels.Names;
            Manual.Touch = converted_data.Touch;
            
            Annotations.Manual = Manual;
            fprintf(' - Manual');
        else
            fprintf(' - Manual file incomplete')            
        end
    elseif any(strcmp('Manual', p.Results.Data))
        fprintf(' - Manual not found');       
    end
    
    
    
    save(save_file , 'Annotations')
    fprintf(' - data saved\n')
    
end






































