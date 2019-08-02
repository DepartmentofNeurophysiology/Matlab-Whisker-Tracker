function p_touch()
%%
[FileName, PathName] = uigetfile('*.avi');
compiledata('file',fullfile(PathName, FileName),...
    'data',{'Tracker'},'overwrite',1)
%
manfile = fullfile(PathName, [FileName(1:end-4) '_touches.mat']);
load(manfile)

mwtfile = fullfile(PathName, [FileName(1:end-4) '_compiled.mat']);
load(mwtfile)

%%
mwtTouch = Annotations.Tracker.Touch;
manFlag = TouchFlag;

mwtres = zeros(size(mwtTouch));


for distidx = 1:10
mwtFlag = zeros(size(mwtTouch,2),1);

for i = 1:length(mwtFlag)
    if any(mwtTouch{distidx, i})
        mwtres(distidx, i) = 1;
    end
end
end
mwtres(end+1,:) = manFlag;

figure(1)
clf
imagesc(mwtres)

keyboard
%%

figure(2)
clf

fnr = 170;
Settings = Annotations.Settings;
Settings.Current_frame = fnr;
frame = LoadFrame(Settings);

imshow(frame)
hold on

Traces = Annotations.Tracker.Traces{fnr};
for t = 1:size(Traces,2)
    plot(Traces{t}(:,2), Traces{t}(:,1), 'r')
end


%%

figure(2)
clf

fnr = 170;
Settings = Annotations.Settings;
Settings.Current_frame = fnr;
frame = LoadFrame(Settings);

imshow(Annotations.Output.ObjectsTouch)
hold on

Traces = Annotations.Tracker.Traces_clean{fnr};
for t = 1:size(Traces,2)
    plot(Traces{t}(:,2), Traces{t}(:,1), 'r')
end


%%
h = 500;
f = figure(1);
clf(f)
f.Position = [100 100 h h*(Settings.Video_heigth/Settings.Video_width)];
ax = axes(f);

vidout = VideoWriter(fullfile(PathName, [FileName(1:end-4) 'res.avi']));
vidout.FrameRate = 30;
open(vidout)
for i = 1:Settings.Nframes
    Settings.Current_frame = i;
    frame = LoadFrame(Settings);
    
    imagesc(ax, frame);
    colormap(ax,'gray')
    axis(ax,'off')
    if manFlag(i) == 1
        hold on
        rectangle('Position',[10 50 60 50],'FaceColor',[0 1 0])
        hold off
    end
    if mwtres(5,i) == 1
        hold on
        rectangle('Position',[340 50 60 50],'FaceColor',[0 1 0])
        hold off
    end
    text(10,30,'MAN','color','r','fontsize',20)
    text(340,30,'MWT','color','r','fontsize',20)

        
    
    frameout = getframe(f);
    writeVideo(vidout, frameout.cdata);
end
    

save(fullfile(PathName, [FileName(1:end-4) '_touches']), 'TouchFlag')
close(f)   
close(vidout)