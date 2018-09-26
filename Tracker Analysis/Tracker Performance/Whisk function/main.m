clc
clear


datadir = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
files = dir(fullfile(datadir,'*.dat'));



% create a video file for each file
for i = 1:2
    try
    disp(i)
   metafile = fullfile( files(i).folder, [files(i).name(1:end-4) '.mat']);
   meta = load(metafile);
   Settings = [];
   Settings.Video = fullfile( files(i).folder, files(i).name );
   Settings.Video_width = meta.Data.Dims(1);
   Settings.Video_heigth = meta.Data.Dims(2);   
   nframes = meta.Data.NFrames;
   
    % Setup figure
    fig = figure(1);
    set(gcf,'position',[100 100 round(1*Settings.Video_heigth) ...
        round(1*Settings.Video_width)]);
    set(gcf,'Units','pixels')
    set(gca,'Units','normalized')
    set(gca,'Position',[0 0 1 1])
    ax = gca;
    colormap(ax,'gray')
   
   outfile = fullfile( files(i).folder, [files(i).name(1:end-4)]);
   
   if exist([outfile '.avi'],'file')
       continue
   end
   
   
   
   vidout = VideoWriter(outfile,'Motion JPEG AVI');
   open(vidout);
   h = waitbar(0, '');
   for j = 1:nframes
       Settings.Current_frame = j;
       frame = LoadFrame(Settings);
       imagesc(ax, frame);
       
       
        frameout = getframe;
        
        writeVideo(vidout, frameout.cdata);
       
       waitbar(j/nframes)
   end
   close(h)
    
    end
end




files = dir( fullfile( datadir,'*.avi'));
fprintf('Starting to track:')
for i = 1:size(files,1)
    try
        
        
    fprintf('\n\t- video %s ...',files(i).name)
    fname = [files(i).name(1:end-4) '.whiskers'];
    
    if exist(fullfile(files(i).folder,fname), 'file')
        continue
    end
    
    if length(files(i).name) > 14 && strcmp(files(i).name(end-12:end-4), 'Annotated')
        continue
    end
    
       if length(files(i).name) > 14 && strcmp(files(i).name(end-14:end-4), 'and_tracker')
        continue
       end
    
          if length(files(i).name) > 14 && strcmp(files(i).name(end-18:end-4), 'and_tracker_raw')
        continue
    end
    
   sysstr = sprintf('cd %s && trace %s %s',files(i).folder,files(i).name,fname);
   [status{i},cmout{i}] = system( sysstr );
    fprintf('done!')
    end
end


%%


for i = 1:2%:size(files,1)
    
      
    if length(files(i).name) > 14 && strcmp(files(i).name(end-12:end-4), 'Annotated')
        continue
    end
   
   vidfile = fullfile( files(i).folder, [files(i).name(1:end-4) '.dat'] );
   clear('Settings')
   makeSettings;
    metafile = fullfile( files(i).folder, [files(i).name(1:end-4) '.mat']);
   meta = load(metafile);
   Settings.Video = vidfile;
   Settings.Video_width = meta.Data.Dims(1);
   Settings.Video_heigth = meta.Data.Dims(2);  
   Settings.batch_mode = 1;
   Settings.Nframes = meta.Data.NFrames; 
   
   [Objects, ~] = ObjectDetection(Settings);

   gapinfo = detectGap(Objects);
   o.Objects = Objects;
   o.gapinfo = gapinfo;
   o = TrackNose(Settings, o);
   
   
    wfile = fullfile( files(i).folder, [ files(i).name(1:end-4) '.whiskers']);
    mfile = fullfile( files(i).folder, [files(i).name(1:end-4) '.mat'] );
    
    MetaData = load(mfile);
    Janelia = ConvertJanelia( wfile );
    Janelia.MetaData = MetaData.Data;
    Janelia.Objects = Objects;
    Janelia.gapinfo = o.gapinfo;
    Janelia.Direction = o.Direction;
    Janelia.Nose = o.Nose;
    Janelia.Headvec = o.AngleVector;    
    Janelia.Traces_clean = CleanJanelia(Janelia);
    Janelia.Parameters = getParams(Janelia, 'clean');
    Janelia.Labels = getLabels(Janelia);
    [Janelia.Touch, ~] = DetectTouch(Janelia, Settings);
    
    
    manfile = fullfile(files(i).folder, [files(i).name(1:end-4) '_Annotations.mat']);
    if exist(manfile, 'file')
        load(manfile)
        % Store unprocessed manual notations
        Manual.RawNotations = CurvesByFrame;
        
        % Convert manual notations format to tracker format
        converted_data = ConvertAnnotations(Manual.RawNotations);
        Manual.Objects = Objects; % Manual tracking data does not contain object detection
        Manual.Nose = o.Nose; % Manual tracking data does not contain nose
        Manual.Headvec = o.AngleVector; % Manual tracking data does not contain headangle
        Manual.Traces = converted_data.Traces;
        Manual.Parameters = getParams(Manual, 'raw');
        %Manual.Labels = converted_data.Labels.Full;
        %Manual.Label_names = converted_data.Labels.Names;
        Manual.Touch = converted_data.Touch;
        
    else
        Manual = [];
    end
    save( fullfile(files(i).folder, [files(i).name(1:end-4) '_Annotations_Janelia']), 'Janelia','Manual')
    
end









