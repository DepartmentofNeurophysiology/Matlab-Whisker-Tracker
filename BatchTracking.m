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

% Check if parameter setup has been ran before
if ~exist('Settings\Settings.mat','file')
    disp('Run ParameterSetup first')
end

makeSettings;


%% Track videos

if Settings.use_parfor
poolobj = gcp;
end


for i = 1:size(Files,1)
   
   
    time_start = clock;
         
    % Generate settings for file tracking  
    Settings.Video = fullfile(PathName, [Files{i} Extension]);
    Settings.batch_mode = 1;
    Settings = getMetaData(Settings);

    
    
    
    % Track video
    [Output.Objects,~, Output.Edges] = ObjectDetection(Settings);
    Output.gapinfo = detectGap(Output.Objects);
    
    Output = TrackNose(Settings, Output);
    
    timer = logspeed([], 20);
    
    % Variables for tracking
    frame_idx = CostumFrameSelection(Settings, Output);
    frame_idx = find(frame_idx);
    Traces = cell(Settings.Nframes,1);

   
   %% 
    if ~Settings.use_parfor    
        tic
        h = waitbar(0,'Tracking Video -');
        count = 0;
        for ii = frame_idx


            % Track frame
            Settings.Current_frame = ii;
            Output.Traces = TrackFrame(Settings, Output);
            Traces{ii} = Output.Traces;

            count = count+1;


            % Update GUI variables   
            timer = logspeed(timer, []);

            time_left = (length(frame_idx) - count) /timer.speed;

            if ~isnan(timer.speed)
                bar_string = sprintf('Tracking video - %d/%d \n%1.2fFPS   Time left: %4.0fs',...
                    count,length(frame_idx),timer.speed,time_left);
            else
                bar_string = sprintf('Tracking video - %d/%d \n   FPS:   Time left:    s', count, length(frame_idx));
            end
            h.Children.Title.String = bar_string;        

            waitbar(count/ length(frame_idx));
            

        end
        close(h)
        Output.ProcessingTime = toc;


    elseif Settings.use_parfor       
        
       tic
       ppm = ParforProgMon('Tracking Whiskers ', length(frame_idx));
        TempTraces = cell(1, length(frame_idx));
        parfor ii = 1:length(frame_idx)
            loopsettings = Settings;
            loopsettings.Current_frame = frame_idx(ii);
            TempTraces{ii} = TrackFrame(loopsettings, Output);   
            ppm.increment();          
        end
        ppm.delete();
        for ii = 1:length(frame_idx)
            Traces{frame_idx(ii)} = TempTraces{ii};            
        end
        Output.ProcessingTime = toc;
        
    end




    Output.Traces = Traces;    
    
    Settings.ExportName = [Settings.Video(1:end-4) '_Annotations_Tracker.mat'];
    % Store tracking resuts
    save( Settings.ExportName,'Output','Settings')   
 
   
   
end


fprintf('Finished!\n')