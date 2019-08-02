clear
close all
clc

addpath(genpath(pwd))

path = 'H:\Data\Julien\Highspeed video\16_Oct';

vid = 'Acq_A_020';



load('Settings\Settings.mat')

Settings.Video = fullfile(path, vid, 'ImgA.avi');
Settings.batch_mode = 0;
Settings = getMetaData(Settings);
 Settings.DefaultDirection = 'Up';
 
Results = getBackground(Settings);
Settings.Current_frame=1;
frame = LoadFrame(Settings);
ROI = setROI(frame);


    % Generate settings for file tracking
 
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
    Objects(1:ymean1,:) = 1;
    ObjectsTouch(1:ymean1,:) = 0;
    ymean2 = round(mean(y(3:4)));
    Objects(ymean2:end,:) = 1;
    ObjectsTouch(ymean2:end,:) = 0;
    
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

for i = 1:Settings.Nframes
    Settings.Current_frame = i;
    frame= LoadFrame(Settings);
    imagesc(frame)
    colormap gray
    hold on
    scatter(Output.Nose(i,2), Output.Nose(i,1), 'r', 'filled')
    %scatter(Annotations.Tracker.Nose(i,2), Annotations.Tracker.Nose(i, 1), 'b', 'filled')
    scatter(data.Output.Nose(i, 2), data.Output.Nose(i, 1), 'g', 'filled')
    hold off
end