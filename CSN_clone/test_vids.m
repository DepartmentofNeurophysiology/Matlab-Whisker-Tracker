%% Select 50 random videos

%{
vidfol = 'G:\HighspeedVideo';
sess = {'21-10','22-10','23-10'};

filesout = {};
for i = 1:3
    files = dir(fullfile(vidfol, sess{i},'*\*.avi'));
    for j = 1:length(files)
        filesout{end+1} = files(j);
    end    
end
%

nvids = 50;
vididx = sort(randperm(length(filesout), nvids));

vidfiles = {};
vidnames = {};
for i = 1:length(vididx)
    vidfiles{end+1} = fullfile(filesout{vididx(i)}.folder, ...
        filesout{vididx(i)}.name);
    names = split(filesout{vididx(i)}.folder,'\');
    vidnames{end+1} = [names{end-1} '_' names{end} '.avi'];        
end

%
outfolder = 'C:\TempVids';
mkdir(outfolder);

for i = 1:length(vidfiles) 
    copyfile(vidfiles{i}, fullfile(outfolder,vidnames{i}))
end

%}



%%


clear
clc

addpath(genpath('C:\code\MWT'))

load('C:\code\MWT\Settings\Settings.mat')
vidfolder = 'C:\TempVids';

Videos = dir(fullfile(vidfolder, '*.avi'));


%% Test a single video
vididx = 1;

Settings.Video = fullfile( Videos(vididx).folder, Videos(vididx).name);
Settings.batch_mode = 1;
Settings = getMetaData(Settings);
Results = getBackground(Settings);

fidx = 1;
Settings.Current_frame = 1;
frame = LoadFrame(Settings);

%%
figure(1)
subplot(1,3,1)
imshow(Results.Objects)

subplot(1,3,2)
imshow(Results.Edges)

subplot(1,3,3)
imshow(frame)

%%

Objects = setROI(Results.Objects,frame);


%%
Results.Objects = Objects;
Settings.DefaultDirection = 'Up';
O = TrackNose(Settings,Results);



%%
figure(1)
for i = 1:Settings.Nframes
    Settings.Current_frame = i;
    f = LoadFrame(Settings);
    clf
    imshow(f)
    hold on
    scatter(O.Nose(i,2), O.Nose(i,1), 'r', 'filled')
    quiver(O.Nose(i,2), O.Nose(i,1), O.AngleVector(i,2)*100, O.AngleVector(i,1)*100, 'r')
    pause(0.01)
end


%%
Settings.track_nose = 0;
Settings.Dilation = 30;
Data.Objects = Objects;
Data.Edges = Results.Edges;
Data.Nose = [NaN, NaN];

tic
[T, M] = TrackFrame(Settings, Data);
toc

figure(1)
clf
imshow(frame)
hold on
for i = 1:size(T, 2)
    plot(T{i}(:,2), T{i}(:,1), 'r')
end










