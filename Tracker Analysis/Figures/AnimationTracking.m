%% Animation tracker explanation
clear
clc
close all

cc = cbrewer('div','RdBu',12);
colors.originprofile = cc(10,:);
colors.origins = cc(3,:);
colors.square = cc(3,:);
colors.red = cc(3,:);


cc = cbrewer('div','RdYlGn',12);
colors.green = cc(10,:);



datapth = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
video = 'M46_R01_02.dat';
load(fullfile(datapth, [video(1:end-4) '_compiled.mat']))
Settings = Annotations.Settings;

Settings.Current_frame = 300;


frame = LoadFrame(Settings);
f_width = Settings.Video_heigth;
f_heigth = Settings.Video_width;

ewidth = 200;

f = figure('Units','pixels','position',[150 150 f_width+ewidth f_heigth]);
ax = axes(gcf,'Units','normalized','Position',[0 0 (f_width/(f_width+ewidth)) 1]);


gifname = 'temp.gif';



%%
% Frame 1 - Raw frame
h = imshow(frame);
colormap gray
makeText(ax,f_width,f_heigth,'pre')


%%
frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'Loopcount', inf, 'DelayTime',0.5)




for alpha = 0:0.1:0.95   
    set(h,'AlphaData',ones(size(frame))*(1-alpha))
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)
    
end


%%

Objects = gray2ind(Annotations.Output.Objects);



%%


% Frame 2 - Extracted objects
cla
imshow(Objects);
[a,b] = find(Objects);
hold on
scatter(b,a,3,'r','MarkerFaceColor', colors.red, 'MarkerEdgeColor', colors.red)
makeText(ax,f_width,f_heigth,'obj')

frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.5)


%}

%%



Silhouette = zeros(size(Objects));
Silhouette( find(frame <= Settings.Silhouettethreshold) ) = 1; %#ok<*FNDSB>
Silhouette(find(Objects)) = 0;
dispSil = Silhouette;
dispSil(dispSil == 0) = 2;
dispSil(dispSil == 1) = 0.1;


%%

% Frame 3 - Silhouette extraction
cla
imshow(frame);
colormap gray
hold on
s = scatter(b,a,3,'r','MarkerFaceColor', colors.red, 'MarkerEdgeColor', colors.red,...
        'MarkerFaceAlpha',0, 'MarkerEdgeAlpha',0);
makeText(ax,f_width,f_heigth,'sil')
    
for alpha = 0.05:0.05:0.4  
    set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
    drawnow
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)
end

cla
imshow(frame);
colormap gray
hold on
scatter(b,a,3,'r','MarkerFaceColor', colors.red, 'MarkerEdgeColor', colors.red,...
    'MarkerFaceAlpha',1, 'MarkerEdgeAlpha',1)
makeText(ax,f_width,f_heigth,'sil')


drawnow
frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)
disp(alpha)

[d,e] = find(Silhouette);
[fr,g] = find(Silhouette==0);


cla
imshow(frame);
colormap gray
hold on
scatter(b,a,3,'r','MarkerFaceColor', colors.red, 'MarkerEdgeColor', colors.red,...
    'MarkerFaceAlpha',1, 'MarkerEdgeAlpha',1)
s1 =  scatter(e,d,20,'w','filled','MarkerFaceAlpha',0,'MarkerEdgeAlpha', 0);
s2 =   scatter(g,fr,20,'k','filled','MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
makeText(ax,f_width,f_heigth,'sil')

for i = 0.01:0.01:0.15
    
    
    set(s1, 'MarkerFaceAlpha', i ,'MarkerEdgeAlpha',i)
    set(s2, 'MarkerFaceAlpha', i, 'MarkerEdgeAlpha', i)
    
    drawnow
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)
end


cla
imshow(Silhouette)
makeText(ax,f_width,f_heigth,'sil')

frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.2)

%}

%%


Shapes.normal = Silhouette;
Shapes.normal = imerode(Silhouette,strel('diamond',2));
%Data.Silhouette( find(Data.Objects) ) = 0;
SilhouetteSmall = imerode(Silhouette,strel('diamond',5));
SE = strel('disk',Settings.Dilationsize);
Silhouette = imdilate(SilhouetteSmall,SE);
TraceOriginIDX = find(edge(Silhouette));
TraceOriginIDX = unique(TraceOriginIDX);

dispSil = Silhouette;
dispSil(dispSil == 0) = 2;
dispSil(dispSil == 1) = 0.1;


%%


% Frame 4 - Dilation
[a,b] = find(~Shapes.normal & Silhouette);

for alpha = 0.01:0.01:0.2
    cla
    imshow(Shapes.normal)
    hold on
    scatter(b,a,'w','filled','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha', alpha)
    makeText(ax,f_width,f_heigth,'edg')

    drawnow
    
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)
    
end

[e1,e2] = find(edge(Silhouette));

for alpha = 0.05:0.05:0.6
    cla
    imshow(Silhouette)
    hold on
    scatter(e2, e1, 3, 'MarkerFaceColor', colors.originprofile, 'MarkerEdgeColor', colors.originprofile,...
        'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
        makeText(ax,f_width,f_heigth,'edg')

    drawnow
    
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)
end


cla
imshow(Silhouette)
hold on
scatter(e2, e1, 10, 'MarkerFaceColor', colors.originprofile, 'MarkerEdgeColor', colors.originprofile,...
    'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1)
    makeText(ax,f_width,f_heigth,'edg')

drawnow

frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.05)

%%

Shapes.large = Silhouette;
Shapes.small = SilhouetteSmall;
TraceOriginIDX = find(edge(Silhouette));
TraceOriginIDX = unique(TraceOriginIDX);
[T1,T2] = ind2sub(size(frame),TraceOriginIDX);
unmarked = 2:length(TraceOriginIDX);
marked = 1;
temp(1,:) = [T1(1),T2(2)];
for j = 1:length(T1)
    dist = sqrt(sum( ([T1(unmarked),T2(unmarked)] - [temp(end,1),temp(end,2)]).^2,2));
    id = find(dist == min(dist));
    temp = [temp;T1(unmarked(id)),T2(unmarked(id))];
    marked = [marked,unmarked(id)];
    unmarked = [unmarked(1:id-1),unmarked(id+1:end)];
end
TraceOriginIDX = sub2ind(size(frame),temp(:,1),temp(:,2));



OriginProfile =  frame(TraceOriginIDX);

% Lowpass filter before peakdetection
LP = medfilt1(OriginProfile,500);
OriginProfile = abs(OriginProfile-LP); %*
Settings.Origin_threshold = 0.02;
[~,OriginsIDX] = findpeaks(OriginProfile,'Threshold',...
    Settings.Origin_threshold);




%% Frame 6 - Raw origins





for alpha = 0:0.05:0.4
    
    cla
    imshow(frame)
    hold on
    scatter(temp(:,2), temp(:,1), 3, 'MarkerFaceColor',colors.originprofile, 'MarkerEdgeColor',colors.originprofile)
    
    scatter(temp(OriginsIDX,2),temp(OriginsIDX,1), 10 ,'MarkerFaceColor', colors.origins,...
        'MarkerEdgeColor', colors.origins,'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha', alpha)
        makeText(ax,f_width,f_heigth,'fndsd')

    
    
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end




%%

Origins = [];
[Origins(:,1),Origins(:,2)] = ind2sub(size(frame),...
    TraceOriginIDX(OriginsIDX));
Shapes.Objects = Objects;
Shapes.Frame = frame;
Shapes.Tracked = zeros(size(Objects));
for idx = 1:size(Origins,1)
    Settings.object_tick = 2;
    [flag,~] = ValidatePoint(Settings, Shapes, Origins(idx,:));
    if ~flag
        Origins(idx,1:2) = NaN;
    end
end


% Track traces
% Rank origins based on their intensity values
id = find(~isnan(Origins(:,1)));
Origins = Origins(id,:);
id = sub2ind(size(frame),Origins(:,1),Origins(:,2));
Origins(:,3) = frame(id);
Origins = sortrows(Origins,3);

flag = 1;
tick = 1;
if size(Origins,1) > 1
    while flag
        if Origins(tick+1,:) == Origins(tick,:)
            if tick < size(Origins,1)-1
                Origins = [Origins(1:tick,:);Origins(tick+2:end,:)];
            elseif tick < size(Origins,1)
                Origins = [Origins(1:tick,:)];
            end
            
        else
            tick = tick+1;
        end
        
        if tick == size(Origins,1)
            flag = 0;
        end
    end
end




%%

for alpha = 0.4:-0.05:0
    
    cla
    imshow(frame)
    hold on
    scatter(temp(:,2), temp(:,1), 3, 'MarkerFaceColor',colors.originprofile,...
        'MarkerEdgeColor',colors.originprofile,'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha', alpha)
    
    scatter(temp(OriginsIDX,2),temp(OriginsIDX,1), 10 ,'MarkerFaceColor', colors.origins,...
        'MarkerEdgeColor', colors.origins,'MarkerFaceAlpha',alpha, 'MarkerEdgeAlpha', alpha)
    scatter(Origins(:,2), Origins(:,1), 10 ,'MarkerFaceColor', colors.green,...
        'MarkerEdgeColor', colors.green,'MarkerFaceAlpha',0.4-alpha, 'MarkerEdgeAlpha', 0.4-alpha+0.1)
        makeText(ax,f_width,f_heigth,'fndsd')

    
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end


cla
imshow(frame)
hold on

scatter(Origins(:,2), Origins(:,1), 10 ,'MarkerFaceColor', colors.green,...
    'MarkerEdgeColor', colors.green,'MarkerFaceAlpha',1, 'MarkerEdgeAlpha', 1)


frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)

%%
Oid = 15;
Origin = Origins(Oid,1:2);
[Trace, Shapes] = TrackTrace(Settings, Shapes, Origin);
dy = max(Trace(:,1)) - min(Trace(:,1));
dx = max(Trace(:,2)) - min(Trace(:,2));
testTrace = Trace;
if dy > dx
    Y1 = min(Trace(:,1))-20;
    Y2 = max(Trace(:,1))+20;
    dy = Y2-Y1;
    ratio = f_width/f_heigth;
    
    dx = dy*ratio;
    X1 = (max(Trace(:,2))+min(Trace(:,1)))/2 - dx/2;
    X2 = (max(Trace(:,2))+min(Trace(:,1)))/2 + dx/2;
    
    
else
end


%%

% Frame 8
xlim([0 f_width])
ylim([0 f_heigth])
for alpha = 0.05:0.05:0.4
    
    cla
    imshow(frame)
    hold on
    
    scatter(Origins(:,2), Origins(:,1), 20 ,'MarkerFaceColor', colors.green,...
        'MarkerEdgeColor', colors.green,'MarkerFaceAlpha',1, 'MarkerEdgeAlpha', 1)
    
    r = rectangle('Position',[X1 Y1 X2-X1 Y2-Y1],'FaceColor',[colors.square alpha],'EdgeColor',[colors.square, alpha*2],...
        'LineWidth',3);
    
     makeText(ax,f_width,f_heigth,'trcks')

    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end

nsteps = 10;
dy1 = Y1/nsteps;
dy2 = (f_heigth-Y2)/nsteps;

dx1 = X1/nsteps;
dx2 = (f_width-X2)/nsteps;

dalpha = alpha/nsteps;
for i = 1:nsteps
    ylim([0+dy1*i f_heigth-dy2*i])
    xlim([0+dx1*i f_width-dx2*i])
    set(r, 'FaceColor', [colors.square alpha-(dalpha*(i))])
    pause(0.1)
    
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end



cla
imshow(frame)
hold on
     makeText(ax,f_width,f_heigth,'trcks')

scatter(Origins(:,2), Origins(:,1), 40 ,'MarkerFaceColor', colors.green,...
    'MarkerEdgeColor', 'k' ,'MarkerFaceAlpha',1, 'MarkerEdgeAlpha', 1)
   
r = rectangle('Position',[X1 Y1 X2-X1 Y2-Y1],'EdgeColor',[colors.square, alpha*2],...
        'LineWidth',3);
frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)

%%



xlim([X1 X2])
ylim([Y1 Y2])


swidth = 50;
x1 = Origin(2)-swidth;
x2 = Origin(2)+swidth;
y1 = Origin(1)- swidth/ratio;
y2 = Origin(1)+ swidth/ratio;


%%


for alpha = 0.4:-0.05:0
    cla 
    imshow(frame)
    scatter(Origins(:,2), Origins(:,1), 40, 'MarkerFaceColor', colors.green, ...
        'MarkerEdgeColor', 'k','MarkerFaceAlpha', alpha,'MarkerEdgeAlpha', alpha)
    rectangle('Position',[x1 y1 x2-x1 y2-y1], 'FaceColor',[colors.square 0.4-alpha],'EdgeColor', colors.square,...
    'LineWidth',3)
    scatter(Origin(2), Origin(1), 50,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
        'MarkerFaceAlpha', 1-(alpha*0.2),'MarkerEdgeAlpha', 1-(alpha*2))
         makeText(ax,f_width,f_heigth,'trcks')

    drawnow
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)

end

    





%%

ref = [0 1];
frame = Shapes.Frame;
SL = Shapes.large;


Shapes.Tracked = zeros(size(Shapes.Tracked));

%Track towards snout
first_stepsize = 3; % initial stepsize is always small


Trace = Origin;
Cx = ceil(Trace(1,1) + first_stepsize*sind(1:10:360)); % x-coordinates
Cy = ceil(Trace(1,2) + first_stepsize*cosd(1:10:360)); % y-coordinates
Cidx = sub2ind(size(frame), Cx, Cy);                   % linear indices
Cidx = unique(Cidx);
% only allow indices pointing toward shape
keep_idx = find(SL(Cidx));


Cidx = Cidx(keep_idx);
[Cx, Cy] = ind2sub(size(frame), Cidx);




xlim([X1 X2])
ylim([Y1 Y2])


% Frame 11 - First ROI
cla
imshow(frame)
hold on
s = scatter(Origin(2), Origin(1), 50,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
        'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
r = rectangle('Position',[x1 y1 x2-x1 y2-y1], 'FaceColor',[colors.square 0.3],'EdgeColor', colors.square,...
    'LineWidth',3);
         makeText(ax,f_width,f_heigth,'trcks')



nsteps = 10;
dy1 = (y1-Y1)/nsteps;
dy2 = (Y2-y2)/nsteps;

dx1 = (x1-X1)/nsteps;
dx2 = (X2-x2)/nsteps;

dalpha = 0.3/nsteps;
for i = 1:nsteps
    ylim([Y1+dy1*i Y2-dy2*i])
    xlim([X1+dx1*i X2-dx2*i])
    
    set(r, 'FaceColor',[colors.square 0.3-dalpha*i])
          

    
   drawnow 
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end


cla
imshow(frame)
hold on
s = scatter(Origin(2), Origin(1), 50,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1);
r = rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor', colors.square,...
    'LineWidth',3);
         makeText(ax,f_width,f_heigth,'trcks')

frameout = getframe(f);
im = frame2im(frameout);
[imind, cm] = rgb2ind(im, 256);
imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)

s = scatter(Cy, Cx, 40, 'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0, 'MarkerEdgeAlpha',0);


for alpha = 0.05:0.1:1
    
    set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
   drawnow 
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end








% Extract intensity in ROI
Proi = frame(Cidx);

% Find pixel with minimum value

if ~isempty(Proi)
    pidx = Cidx( find( Proi == min(Proi),1,'first'));
else
    pidx = [];
end
newpt = zeros(1,2);
if ~isempty(pidx)
    [newpt(1), newpt(2)] = ind2sub(size(frame), pidx);
else
    newpt = [1,1]; % This values will be rejected
end

[ptflag, Settings] = ValidatePoint( Settings, Shapes, newpt);
if ptflag
    Trace(end+1,:) = newpt;
end

first_stepsize = 5;



s = scatter(Trace(end,2), Trace(end,1),50, 'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
for alpha = 0:0.2:1
      set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
      
      
   drawnow 
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end








while ptflag
    
    %% Frame 13
    cla
    imshow(frame)
    hold on
    scatter(Trace(:,2), Trace(:,1), 50,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k')
    scatter(Origin(2), Origin(1), 70,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k')
    
         makeText(ax,f_width,f_heigth,'trcks')

    
    %%
    
    vt = Trace(end,:) - Trace(end-1,:);
    Angle = atan2d(vt(1), vt(2));
    HalfCircle = -30:30;
    Theta = HalfCircle + Angle;
    
    Cx = ceil(Trace(end,1) + first_stepsize*sind(Theta));
    Cy = ceil(Trace(end,2) + first_stepsize*cosd(Theta));
    
  s = scatter(Cy, Cx, 40, 'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0, 'MarkerEdgeAlpha',0);


    for alpha = 0.05:0.1:1

        set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
       drawnow 
        frameout = getframe(f);
        im = frame2im(frameout);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
    end
    


    Cidx = unique(sub2ind(size(frame),Cx,Cy));
    % Extract intensity in ROI
    Proi = frame(Cidx);
    
    % Find pixel with minimum value
    
    if ~isempty(Proi)
        pidx = Cidx( find( Proi == min(Proi),1,'first'));
    else
        pidx = [];
    end
    
    if ~isempty(pidx)
        [newpt(1), newpt(2)] = ind2sub(size(frame), pidx);
    else
        newpt = [1,1]; % This values will be rejected
    end
    
    [ptflag, Settings] = ValidatePoint( Settings, Shapes, newpt);
    if ptflag
        Trace(end+1,:) = newpt;
    end
    
    
    
    %% Frame 14
    s = scatter(Trace(end,2), Trace(end,1),50, 'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
        'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
    for alpha = 0:0.2:1
        set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
        
        
        drawnow
        frameout = getframe(f);
        im = frame2im(frameout);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
    end
    %%
    % Frames Toward Snout
    
    
    
end


%%
Trace = flip(Trace, 1);


if size(Trace,1) > 2
    xdata= Trace(:,1);
    ydata = Trace(:,2);
    px = polyfit(1:length(xdata),xdata',1);
    py = polyfit(1:length(ydata),ydata',1);
    
    
    xdataf = polyval(px,[1, length(xdata)]);
    ydataf = polyval(py,[1, length(ydata)]);
    
    vt = [xdataf(2)-xdataf(1), ydataf(2)-ydataf(1)];
    Angle = atan2d(vt(1), vt(2));
    HalfCircle = Settings.circle_start:Settings.circle_end;
    Theta = HalfCircle + Angle;
    Cx = ceil(Trace(end,1) + first_stepsize*sind(Theta));
    Cy = ceil(Trace(end,2) + first_stepsize*cosd(Theta));
    
    
else
    Cx = ceil(Trace(1,1) + first_stepsize*sind(1:10:360)); % x-coordinates
    Cy = ceil(Trace(1,2) + first_stepsize*cosd(1:10:360)); % y-coordinates
    Cidx = sub2ind(size(frame), Cx, Cy);                   % linear indices
    Cidx = unique(Cidx);
    % only allow indices pointing toward shape
    keep_idx = find(~SL(Cidx));
    
    Cidx = Cidx(keep_idx);
    [Cx, Cy] = ind2sub(size(frame), Cidx);
end



Cidx = unique(sub2ind(size(frame),Cx,Cy));

Proi = frame(Cidx);

% Find pixel with minimum value
if ~isempty(Proi)
    pidx = Cidx( find( Proi == min(Proi),1,'first'));
else
    pidx = [];
end

if ~isempty(pidx)
    [newpt(1), newpt(2)] = ind2sub(size(frame), pidx);
else
    newpt = [1,1]; % This values will be rejected
end

[ptflag, Settings] = ValidatePoint( Settings, Shapes, newpt);
if ptflag
    Trace(end+1,:) = newpt;
end

HalfCircle = Settings.circle_start:Settings.circle_end;


noiseflag = 0;


%% 15


% Frame - First point toward tip
cla
imshow(frame)
hold on
scatter(Trace(:,2), Trace(:,1), 50,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k')
scatter(Origin(2), Origin(1), 70,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k')


s = scatter(Cy, Cx, 40, 'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0, 'MarkerEdgeAlpha',0);

         makeText(ax,f_width,f_heigth,'trckt')

for alpha = 0.05:0.1:1
    
    set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
    drawnow
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end



s = scatter(Trace(end,2), Trace(end,1),50, 'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
for alpha = 0:0.2:1
    set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
    
    
    drawnow
    frameout = getframe(f);
    im = frame2im(frameout);
    [imind, cm] = rgb2ind(im, 256);
    imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.1)
end
%%

s0 = 50;
send = max(testTrace(:,1)) - min(testTrace(:,1));
ntotrack = size(testTrace,1) - size(Trace,1);
ds = (send-s0)/ntotrack;
tick = 1;

%
while ptflag
    vt = Trace(end,:) - Trace(end-1,:);
    
    Angle = atan2d(vt(1), vt(2));
    
    if size(Trace,1) > 9
        HalfCircle = Settings.circle_start-10:Settings.circle_end + 10; % loose up roi near tip
    end
    
    Theta = HalfCircle + Angle;
    Cx = ceil(Trace(end,1) + Settings.stepsize*sind(Theta));
    Cy = ceil(Trace(end,2) + Settings.stepsize*cosd(Theta));
    Cidx = unique(sub2ind(size(frame),Cx,Cy));
    
    % Extract intensity in ROI
    Proi = frame(Cidx);
    
    % Find pixel with minimum value
    if ~isempty(Proi)
        lowestpoint = find( Proi == min(Proi), 1, 'first');
        lowestpoints = find( Proi - Proi(lowestpoint) < 0.01);
        
        pidx = Cidx( lowestpoints );
    else
        pidx = [];
    end
    
    if length(pidx) > 1 & size(Trace,1) > 3% if multiple pixels are competing
        vtrace = Trace(end-3,:) - Trace(end,:);
        atrace = atan2d(vtrace(2), vtrace(1));
        
        apt = []; % find point with smallest angle w.r.t. current trace
        for ii = 1:length(pidx)
            [newpt(1), newpt(2)] = ind2sub(size(frame),pidx(ii));
            vpt = Trace(end-3,:) - newpt;
            apt(ii) = abs( atrace- atan2d(vpt(2),vpt(1)) );
        end
        
        lowestval = find( apt == min(apt), 1, 'first');
        pidxout = pidx( lowestval );
        [newpt(1), newpt(2)] = ind2sub( size(frame), pidxout );
        
    elseif ~isempty(pidx)
        [newpt(1), newpt(2)] = ind2sub(size(frame), pidx(1));
    else
        newpt = [1,1]; % This values will be rejected
    end
    
    Settings.current_length = size(Trace,1);
    [ptflag, Settings] = ValidatePoint( Settings, Shapes, newpt);
    if ptflag
        Trace(end+1,:) = newpt;
        
        if noiseflag > 0
            noiseflag =0;
        end
        
        
        
        
    elseif Settings.stop.noise
        if noiseflag < 0
            Trace(end+1,:) = newpt;
            noiseflag = noiseflag + 1;
            ptflag = 1;
        end
        
    else % Use extrapolation to scan for whiskerpoints further away
        % To pass objects/noise around the tip
        
        if size(Trace,1) > 5 & Settings.stop.object
            xdata = Trace(:,1);
            ydata = Trace(:,2);
            
            px = polyfit(1:length(xdata), xdata',2);
            py = polyfit(1:length(ydata), ydata',2);
            fitax = 1:length(xdata)+Settings.extrapolationsize;
            loop_trace = [];
            loop_trace(:,1) = round(polyval(px, fitax));
            loop_trace(:,2) = round(polyval(py, fitax));
            
            check = [];
            for i = 1:size(loop_trace,1)
                [ptflag, Settings] = ValidatePoint(Settings, Shapes, loop_trace(i,:));
                if ptflag
                    check(i) = 1;
                else
                    check(i) = 0;
                end
            end
            
            newidx = length(xdata) + find(check(length(xdata) + 1:end), 1, 'first');
            if ~isempty(newidx)
                Trace(end+1,1:2) = loop_trace(newidx,:);
                ptflag = 1;
            else
                ptflag = 0;
                
            end
        else
            ptflag =0;
        end
        
    end
    
    swidth = s0+(ds*tick);
    tick = tick+1;
    x1 = mean(Trace(:,2))-swidth;
    x2 = mean(Trace(:,2))+swidth;
    y1 = mean(Trace(:,1))- swidth/ratio;
    y2 = mean(Trace(:,1))+ swidth/ratio;
    
    
    %%  Frames towards tip 16
    cla
    imshow(frame)
        xlim([x1 x2])
    ylim([y1 y2])
    
    hold on
    scatter(Trace(:,2), Trace(:,1), 50,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k')
    scatter(Origin(2), Origin(1), 70,'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k')
             makeText(ax,f_width,f_heigth,'trckt')

    s = scatter(Cy, Cx, 40, 'MarkerFaceColor', colors.red, 'MarkerEdgeColor', 'k',...
    'MarkerFaceAlpha',0, 'MarkerEdgeAlpha',0);


    for alpha = 0:0.25:1

        set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
        drawnow
        frameout = getframe(f);
        im = frame2im(frameout);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.005)
    end
    


    s = scatter(Trace(end,2), Trace(end,1),50, 'MarkerFaceColor', colors.green, 'MarkerEdgeColor', 'k',...
        'MarkerFaceAlpha',0,'MarkerEdgeAlpha',0);
    for alpha = 0:0.25:1
        set(s, 'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)


        drawnow
        frameout = getframe(f);
        im = frame2im(frameout);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, gifname, 'gif', 'WriteMode', 'append', 'DelayTime',0.005)
    end
    

    
end


%%
Output = TrackFrame(Settings, Annotations.Output);
TracesFrame = Output.Traces;

figure(1)
axes();
cla
imshow(frame)
hold on
for i = 1:size(TracesFrame,2)
    t=  TracesFrame{i};
    plot(t(:,2), t(:,1), 'b')
    text(t(end,2), t(end,1), num2str(i))
    
end

xlim([0 f_width])
ylim([0 f_heigth])

x = 1;



%%





function makeText(ax, f_width, f_heigth, name)
tmargx = 10;
tmargy = 70;
tmargyin = 15;

t.pre = text(ax,f_width+tmargx,f_heigth-tmargy,'Preprocessed frame','Units','pixels','color','k');
t.obj = text(ax,f_width+tmargx,f_heigth-tmargy-tmargyin,'Extract objects','Units','pixels','color','k');
t.sil = text(ax,f_width+tmargx,f_heigth-tmargy-tmargyin*2,'Extract Silhouette', 'Units','pixels','color','k');
t.edg = text(ax,f_width+tmargx,f_heigth-tmargy-tmargyin*3,'Define seed ROI', 'Units','pixels','color','k');
t.fndsd = text(ax,f_width+tmargx, f_heigth-tmargy-tmargyin*4,'Track trace seeds','Units','pixels','color','k');
t.int = text(ax,f_width+tmargx, f_heigth-tmargy-tmargyin*5,'Track individual traces:','Units','pixels','color','k');
t.trcks = text(ax,f_width+tmargx, f_heigth-tmargy-tmargyin*6,'   -Converge to snout','Units','pixels','color','k');
t.trckt = text(ax,f_width+tmargx, f_heigth-tmargy-tmargyin*7,'   -Converge to tip','Units','pixels','color','k');
set(t.(name), 'FontWeight','bold','Color',[1 0 0])

if strcmp(name,'trcks') || strcmp(name,'trckt')
    set(t.int, 'FontWeight','bold','Color',[1 0 0])
end
end
