function [res] = timeMWT(varargin)
% Time processing speed of MWT by running the listed videos
p = inputParser;

addParameter(p,'files',[]);

parse(p, varargin{:});

if isempty(p.Results.files)
    disp('This function requires the name-value pair: ''files'' - ''filelist''')
    return
end

files = p.Results.files;
for i = 1:size(files,2)
    timer(i).name =  files{i};
    timer(i).start = clock;
    
    makeSettings;
    slash_idx = find(files{i} == '\', 1,'last');
    Settings.PathName = files{i}(1:slash_idx-1);
    Settings.FileName = files{i}(slash_idx+1:end);
    Settings.Video = fullfile(Settings.PathName, Settings.FileName);
    Settings.batch_mode = 1;
    m_file = Settings.Video;
    m_file(end-2) = 'm';
    load(m_file);
    Settings.Video_width = Data.Resolution(1);
    Settings.Video_heigth = Data.Resolution(2);
    Settings.Nframes = Data.NFrames;  
    
    timer(i).MWTsetup = clock;
    
    [Output.Objects,~] = ObjectDetection(Settings);       
    Output.gapinfo = detectGap(Output.Objects);  
    
    timer(i).MWTobjectdetection = clock;
    
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
    
    timer(i).MWTtracetracking = clock;
end

res = timer;
    
    
    
    
    
    
    
    