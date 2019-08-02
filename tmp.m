viddir = 'G:\HighspeedVideo\16_Oct';
trials = dir(viddir);

for i = 60
    if trials(i).isdir
        if ~isfile(fullfile(trials(i).folder, trials(i).name, 'touches.mat'))
            continue
        end
        load(fullfile(trials(i).folder, trials(i).name, 'touches.mat'))
        load(fullfile(trials(i).folder, trials(i).name, 'ImgA_Annotations_Tracker.mat'))
        firsttouch = find(touch,1,'first');
        lasttouch = find(touch, 1, 'last');
        
        if ~isempty(firsttouch)
            Settings.Current_frame = firsttouch;
            Settings.Video = fullfile(viddir, trials(i).name, 'ImgA.avi');
            f = LoadFrame(Settings);
          
            fr = figure('Units','normalized','Position',[0 0 1 1]);
            ax = axes(fr,'Units','points','Position',[0,0m]);
            imagesc(ax,f)
            
             Settings.Current_frame = firsttouch-1;
             subplot(2,2,2)
             imshow(f)
            
            Settings.Current_frame = lasttouch;
            f = LoadFrame(Settings);            
            subplot(2,2,3)
            imshow(f)
            
            Settings.Current_frame = lasttouch+1;
            f = LoadFrame(Settings);
            subplot(2,2,4)
            imshow(f)
            
        end
            
        
    end
end

%%

i = 10;
load(fullfile(trials(i).folder, trials(i).name, 'ImgA_Annotations_Tracker.mat'))
Settings.Video = fullfile(viddir, trials(i).name, 'ImgA.avi');

for j=  1:10
    Settings.Video = fullfile(viddir, trials(i).name, 'ImgA.avi');
    Settings.Current_frame = j;
    f = LoadFrame(Sett);
    figure
    imshow(f)
end



