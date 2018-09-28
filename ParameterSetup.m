%% Initialize
% clear workspace, add path, make settings for tracking
clear
clc
addpath(genpath(pwd))
makeSettings;


%% Load Video data and corresponding metadata

[Settings.FileName, Settings.PathName] = uigetfile(fullfile(Settings.default_video_path,'*.*'),'Select video file');
Settings.Video = fullfile(Settings.PathName,Settings.FileName);


% Metadata required for tracking are the Video width and heigth and the
% number of frames in the video. Those are native properties of video
% objects, but not in costum data files like our .dat format. In our case
% we load a second file containing metadata.

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
        return
    end
elseif Settings.video_extension == '.mat'
       m = matfile(Settings.Video);
       data = m.movf(1,1);
       Settings.Video_width = size(data.cdata,1);
       Settings.Video_heigth = size(data.cdata,2);
       Settings.Nframes = size(m.movf, 2);
     
       
       
    
else
    try
        Video_object = VideoReader(Settings.Video);
        Settings.Video_width = Video_object.Height;
        Settings.Video_heigth = Video_object.Width;
        Settings.Nframes = floor(Video_object.Duration * Video_object.FrameRate);
        Settings.Video_object = Video_object;
    catch
        disp('The Video is not readable with ''VideoReader''')
        return
    end    
end



%% Data Tracking and display

% Track Objects and Nose
[Output.Objects,Settings.object_threshold] = ObjectDetection(Settings);
Output = TrackNose(Settings, Output);

% Open a GUI to display results of tracking on a single frame
Settings.definite_settings = 0;
while Settings.definite_settings == 0
    Settings = FrameSelection(Settings);
    Settings = TrackingParameters(Settings, Output);
end


save(fullfile(pwd,'Settings','Settings.mat'),'Settings')


