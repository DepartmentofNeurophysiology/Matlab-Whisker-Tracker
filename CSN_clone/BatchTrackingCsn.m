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
    ParameterSetup;
end

load('Settings\Settings.mat')

%%
setFlag = 0;
i = 1;
while setFlag == 0
    Settings.Video = fullfile(PathName, [Files{i} Extension]);
    Settings.Current_frame = 100;
    frame = LoadFrame(Settings);
    ROI = setROI(frame);
    
    figure(1)
    clf
    
    
    for j = 1:25
        Settings.Video = fullfile(PathName, [Files{j} Extension]);
        Settings.Current_frame = 100;
        frame = LoadFrame(Settings);
        subplot(5,5,j)
        imshow(frame)
        hold on
        scatter(ROI(:,1), ROI(:,2), 'g','filled')
    end
    
    res = inputdlg('satisfied');
    if res{1} == 'y'
        setFlag = 1;
    end
    %}
end
close all



%% Track videos


False_videos = {};

for i = 1:size(Files,1)
    disp([i, size(Files,1)])
    
    time_start = clock;
    
    % Generate settings for file tracking
    Settings.Video = fullfile(PathName, [Files{i} Extension]);
    Settings.batch_mode = 1;
    Settings = getMetaData(Settings);
    Settings.track_nose = 1;
    Settings.DefaultDirection = 'Up';

    
    % Track video
    Results = getBackground(Settings);
    Settings.Current_frame=1;
    test_frame = LoadFrame(Settings);
    
    
   
    Output.ObjectsNative = Results.Objects;
    Output.Edges = Results.Edges;
    Output.gapinfo = Results.gapinfo;
    Output.ROI = ROI;
    
    x = Output.ROI(:,1);
    y = Output.ROI(:,2);
    Objects = Results.Objects;
    ObjectsTouch = Results.Objects;
    
    ymean1 = round(mean(y(1:2)));
    ymean2 = round(mean(y(3:4)));
    nrows = ymean2-ymean1;
    quadrantsize = 0.25*nrows;
    ymax = ymean1 + round(quadrantsize);
    Objects(1:ymean1,:) = 1;
    ObjectsTouch(1:ymean1,:) = 0;
   
    Objects(ymax:end,:) = 1;
    ObjectsTouch(ymax:end,:) = 0;
    
    xmean1 = round(mean([x(1), x(4)]));
    xmean2 = round(mean(x(2:3)));
    

    Objects(:,1:xmean1) = 1;
    Objects(:,xmean2:end) = 1;
    ObjectsTouch(:,1:xmean1) = 0;
    ObjectsTouch(:,xmean2:end) = 0;
    Output.Objects = Objects;
    Output.ObjectsTouch = ObjectsTouch;
    
    if Settings.track_nose
       Output = TrackNose(Settings, Output);
    end
  %% 
    timer = logspeed([], 20);
    
    % Variables for tracking
    frame_idx = CostumFrameSelection(Settings, Output);
    frame_idx = find(frame_idx);
    Traces = cell(Settings.Nframes,1);
    
    if Settings.use_parfor
        poolobj = gcp;
    end
    
    
        %%
        if ~Settings.use_parfor
            tic
            h = waitbar(0,'Tracking Video -');
            count = 0;
            for ii = frame_idx
                
                
                % Track frame
                Settings.Current_frame = ii;
               
                [Traces{ii}, ~] = TrackFrame(Settings, Output);
                
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
                [TempTraces{ii}, ~] = TrackFrame(loopsettings, Output);
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
        
        tname = split(Settings.Video,'\');
        dpath = tname{1};
        for j = 2:length(tname)-1
            dpath = [dpath '\' tname{j}];
        end
        
       compiledata('file',Settings.Video,'data',{'Tracker'},'overwrite',1)
       makeFileForPython( [Settings.Video(1:end-4) '_compiled.mat'])
      %   PRINT_VIDEO('dPath',dpath,'FileName',Settings.FileName(1:end-4),...
      %       'dTclean',1,'dNose',1,'dTtouch',1,'ROI',1,'FrameSelect','full','dExp',1)
        
    
        
        
   %     False_videos{end+1} = Settings.Video;
    
end


fprintf('Finished!\n')