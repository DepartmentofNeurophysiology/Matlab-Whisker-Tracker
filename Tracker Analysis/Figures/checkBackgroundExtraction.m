clear
path = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';

vids = dir(fullfile(path,'*.dat'));

for i = 1:size(vids,1)
    clear Settings
    makeSettings;
    Settings.Video = fullfile( vids(i).folder, vids(i).name);
    Settings.PathName = vids(i).folder;
    Settings.FileName = vids(i).name;
    mfile = Settings.Video;
    mfile(end-2) = 'm';
    load(mfile)
    Settings.Video_width = Data.Resolution(1);
    Settings.Video_heigth = Data.Resolution(2);
    Settings.Nframes = Data.NFrames;
    Settings.batch_mode = 1;
    
    [Objects{i}, ~] = ObjectDetection(Settings);
    ginf{i} = detectGap(Objects{i});
    
end


%%
figure(1)
mkdir(fullfile(path, 'background'))
for i = 1:size(Objects,2)
    clf
    clear Settings
    makeSettings;
    Settings.Video = fullfile( vids(i).folder, vids(i).name);
    Settings.PathName = vids(i).folder;
    Settings.FileName = vids(i).name;
    mfile = Settings.Video;
    mfile(end-2) = 'm';
    load(mfile)
    Settings.Video_width = Data.Resolution(1);
    Settings.Video_heigth = Data.Resolution(2);
    Settings.Nframes = Data.NFrames;
    Settings.Current_frame = round(Settings.Nframes*0.8);
    frame = LoadFrame(Settings);
    
    
    imagesc(frame)
    colormap gray
    
    
    [ia,ib] = find(edge(Objects{i}));
    hold on
    scatter(ib,ia,1,'r','filled')
    
    if ginf{i}.main_ax == 1
        line([ginf{i}.edge_1 ginf{i}.edge_1], [0 Settings.Video_width],'color','g','LineWidth',2)
        line([ginf{i}.edge_2 ginf{i}.edge_2], [0 Settings.Video_width],'color','g','LineWidth',2)

        
    elseif ginf{i}.main_ax == 2
        line([0 Settings.Video_heigth], [ginf{i}.edge_1 ginf{i}.edge_1],'color','g','LineWidth',2)
        line([0 Settings.Video_heigth], [ginf{i}.edge_2 ginf{i}.edge_2],'color','g','LineWidth',2)

    end
    
    saveas(gcf,fullfile(path,'background',[vids(i).name(1:end-4) '.png']))
end