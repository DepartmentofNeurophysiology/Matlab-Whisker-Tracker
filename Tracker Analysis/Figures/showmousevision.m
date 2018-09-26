clear
load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_compiled.mat')

Settings = Annotations.Settings;
Settings.PathName(1) = 'D';
Settings.Video(1) = 'D';
%%

figure(1)
colormap('gray')
frame_disp = zeros(301,301);
for i = 1:Settings.Nframes
clf
Settings.Current_frame = i;
frame = LoadFrame(Settings);
imagesc(frame)
hold on
for j = 1:size(Annotations.Tracker.Traces_clean{i},2)
    t = Annotations.Tracker.Traces_clean{i}{j};
    plot(t(:,2), t(:,1), 'r')


end
hold off

fr = getframe(gcf);
framesave{i} = fr.cdata;
end

%%

for i = 1500:2500
    if isempty(framesave{i})
        continue
    end
    
    

vec = Annotations.Tracker.Headvec(i,:);
angle = atan2d(vec(1),vec(2))-90;


cut.x1 = 150;
cut.x2 = 300;
cut.y1 = 200;
cut.y2 = 200;

nose = round(Annotations.Tracker.Nose(i,:));

x1 = nose(1)-cut.x1;
if x1 < 1
    x1 = 1;
end

x2 = nose(1) + cut.x2;
if x2 > size(framesave{i},1) 
    x2 = size(framesave{i},1);
end

y1 = nose(2) - cut.y1;
if y1 < 1
    y1 = 1;
end

y2 = nose(2) + cut.y2;
if y2 > size(framesave{i},2)
    y2 = size(framesave{i},2);
end

frame_cut = framesave{i}(x1:x2,y1:y2);
frame_cut = imrotate(frame_cut, angle,'crop');

frame_disp = zeros(cut.x1+cut.x2+1,cut.y1+cut.y2+1);
dx = size(frame_disp,1) - size(frame_cut,1);
dy = size(frame_disp,2) - size(frame_cut,2);
if dx < 1 && dy < 1
    frame_disp = frame_cut;
elseif dx < 1 & dy > 0    
frame_disp(: ,...
    round(dy/2):round(dy/2)-1+size(frame_cut,2)) = frame_cut;
elseif dx > 0 & dy < 1
    frame_disp(round(dx/2):round(dx/2)-1+size(frame_cut,1),...
    :) = frame_cut;
elseif dx > 0 & dy > 0
    frame_disp(round(dx/2):round(dx/2)-1+size(frame_cut,1),...
    round(dy/2):round(dy/2)-1+size(frame_cut,2)) = frame_cut;
end

imagesc(frame_disp)

text(10,10,num2str(i),'color','r','BackgroundColor','k')
drawnow

pause(0.05)
end
