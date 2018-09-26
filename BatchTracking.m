% This script does:
%   - Show a GUI to select videos for tracking
%   - Export tracking results to 'outpath' (see makeSettings.m) under a new
%   name ('Video_##')
%   - Saves tracker performance in 'tracker_log.txt'
%   - Maintains a table with tracked videos in 'outpath'

%% Initialize
clear
close all
clc

addpath(genpath(pwd))

% Select video files
[PathName,Files,Extension] = BatchProcessing;
for i = 1:size(Files,1)
    ImportNames{i,1} = fullfile(PathName,[Files{i} Extension]);
end

% Check if parameter setup has been ran before
if ~exist('Settings\Settings.mat','file')
    disp('Run ParameterSetup first')
end

makeSettings;

%%
% Open table with tracking results in outpath
if exist(fullfile(Settings.outpath,'Tracked_Videos.mat'), 'file')
    load(fullfile(Settings.outpath,'Tracked_Videos.mat'))
    nfiles = size(Files,1);
    nvids = size(Tracked_Videos,1);
    for i = 1:size(Files,1)
        ExportNames{i,1} = [ImportNames{i,1}(1:end-4) '_Annotations_Tracker.mat'];
        
        
        
        Tracked_Videos = [Tracked_Videos; {ImportNames{i},ExportNames{i},...
            NaN, NaN, NaN, NaN}];
        FilesIdx(i) = size(Tracked_Videos,1);
    end
else
    for i = 1:size(Files,1)
        ExportNames{i,1} = [ImportNames{i,1}(1:end-4) '_Annotations_Tracker.mat'];
        TotalTime(i,1) = NaN;
        FrameCount(i,1) = NaN;
        FPS(i,1) = NaN;
        TraceCount(i,1) = NaN;
        FilesIdx(i) = i;
    end
    Tracked_Videos = table(ImportNames,ExportNames,TotalTime,FrameCount,...
        FPS,TraceCount);
end

% Add new session to log file
time = clock;




%% Track videos





for i = 1:size(Files,1)
    
    % clear
    % Use try/catch so that single errors do not stop tracking in other
    % videos
    
    
    % Print video details
    time_start = clock;
    full_name = Tracked_Videos.ImportNames{FilesIdx(i)};
    slash_idx = find(full_name == '\',1,'last');
    
    
    % Generate settings for file tracking
    Settings.PathName = full_name(1:slash_idx-1);
    Settings.FileName = full_name(slash_idx+1:end);
    Settings.Video = fullfile(Settings.PathName,Settings.FileName);
    Settings.batch_mode = 1;
    Settings.outpath = Settings.PathName;
    
    % Load video details
    if Settings.use_external_specfile
        try
            m_file = fullfile(Settings.PathName,Settings.FileName);
            m_file(end-2) = 'm';
            load(m_file)
            Settings.Video_width = Data.Resolution(1);
            Settings.Video_heigth = Data.Resolution(2);
            Settings.Nframes = Data.NFrames;
        catch
            disp('Make sure to turn of --use external specfile-- in settings or update the section loading the video specifications.')
        end
        
        
    else
        try
            Settings.Video = fullfile(Settings.PathName, Settings.FileName);
            Video_object = VideoReader(Settings.Video);
            Settings.Video_width = Video_object.Height;
            Settings.Video_heigth = Video_object.Width;
            Settings.Nframes = floor(Video_object.Duration * Video_object.FrameRate);
            Settings.Video_object = Video_object;
        catch
            disp('The Video is not readable with ''VideoReader''')
        end
        
    end
    
    
    
    % Track video
    [Output.Objects,~] = ObjectDetection(Settings);
    Output.gapinfo = detectGap(Output.Objects);
    
    Output = TrackNose(Settings, Output);
    
    
    
    % Initialize empty variables
    n_frames = Settings.Nframes;
    n_tracked = 0;
    Traces = cell(n_frames,1);
    Origins = cell(n_frames,1);
    h = waitbar(0,'Tracking Video -');
    time_buffer_size = 20; % Size of buffer to calculate processing speed
    timestamps = zeros(1,time_buffer_size);
    
    frame_idx = CostumFrameSelection(Settings, Output);
   
    
    for framenr = 1:n_frames
        
        if frame_idx(framenr)
            % Track frame
            Settings.Current_frame = framenr;
            Output = TrackFrame(Settings, Output);
            Traces{framenr} = Output.Traces;
            Origins{framenr} = Output.Origins;
            
        else
            Traces{framenr} = {};
            Origins{framenr} = {};
            continue
        end
        
        
        % Update GUI variables
        n_tracked = n_tracked+1;
        time = clock;
        timestamps = circshift(timestamps,-1);
        timestamps(end) = time(4)*3600 + time(5)*60 + time(6);
        elapsed_time = timestamps(end) - timestamps(1);
        n_frames_left = n_frames - n_tracked;
        track_speed = time_buffer_size/elapsed_time;
        time_left = n_frames_left/track_speed;
        bar_string = sprintf('Tracking video - %d/%d \n@%1.2fFPS   Time left: %4.0fs',framenr,n_frames,track_speed,time_left);
        h.Children.Title.String = bar_string;
        
        
        waitbar(n_tracked/numel(find(frame_idx)));
        
        
    end
    
    close(h)
    Output.Traces = Traces;
    Output.Origins = Origins;
    
    
    % Store tracking resuts
    save(Tracked_Videos.ExportNames{FilesIdx(i)},'Output','Settings')
    
    
    % Print tracking stats
    time_end = clock;
    n_empty = numel(find(isnan(Output.Nose(:,1))));
    n_empty = 0;
    n_tot = Settings.Nframes;
    n_tracked=  n_tot-n_empty;
    tot_time = time_end(4)*3600+time_end(5)*60+time_end(6) - ...
        time_start(4)*3600 - time_start(5)*60 - time_start(6) ;
    fps = n_tracked/tot_time;
    n_traces = 0;
    for q = 1:size(Output.Traces,1)
        n_traces = n_traces + size(Output.Traces{q},2);
    end
    Tracked_Videos.TotalTime(FilesIdx(i)) = tot_time;
    Tracked_Videos.FrameCount(FilesIdx(i)) = n_tot;
    Tracked_Videos.FPS(FilesIdx(i)) = fps;
    Tracked_Videos.TraceCount(FilesIdx(i)) = n_tracked;
    
    
    
    
    
    
end

save(fullfile(Settings.outpath,'Tracked_Videos'),'Tracked_Videos')

fprintf('Finished!\n')