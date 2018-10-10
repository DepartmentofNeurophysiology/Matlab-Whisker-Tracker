

clear 
clc
close all
warning('off')





video = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02.dat';
comp = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_compiled.mat';

load(comp)

Settings = Annotations.Settings;
Settings.Current_frame = 300;
Frame = LoadFrame(Settings);
Objects = Annotations.Tracker.Objects;
Edges = Annotations.Tracker.Edges;
Nose = Annotations.Tracker.Nose(300,:);

ax_heigth = 150;
ax_width = ax_heigth*(Settings.Video_heigth/Settings.Video_width);

marg_x = 50;
x_space = 10;
marg_y = 50;
y_space = 10;

fig_width = 3*ax_width+2*marg_x+2*x_space;
fig_heigth = 3*ax_heigth+2*marg_y+2*y_space;

f = figure('Units','points','Position',[10 100 fig_width fig_heigth]);

x = marg_x;
y = fig_heigth - marg_y - ax_heigth;  
ax1 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]); %#ok<*SAGROW>
x = x + x_space + ax_width;
ax2 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax22 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax23 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax24 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
x = x + x_space+ ax_width;
ax3 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax30 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax31 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax32 = axes(f,'Units','points','Position',[x y ax_width 40]);
ax33 = axes(f,'Units','points','Position',[x y 20 ax_heigth]);
set(ax33,'view',[90,-90])


x = marg_x;
y = y - ax_heigth - y_space;  
ax4 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax40 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax41 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
x = x + x_space + ax_width;
ax5 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax51 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);

x = x + x_space + ax_width;
ax6 = axes(f, 'Units', 'points', 'Position', [x y ax_width ax_heigth]);

y = y - ax_heigth - y_space;
x = marg_x;
ax7 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax71 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);

x = x + x_space + ax_width;
ax8 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);






hold(ax33,'on')









cmap = cbrewer('seq','Reds',20);
c.Objects =[0 0 0; cmap(17,:)];

cmap = cbrewer('div','RdBu',20);

c.Nose = cmap(17,:);
c.Roi = c.Nose;
c.Roi2 = cmap(19,:);
c.Red = cmap(5,:);
cmap = cbrewer('div','PuOr',20);
c.Snout = cmap(5,:);
c.Snout2 = cmap(2,:);
c.Interest = cmap(15,:);
cmap = cbrewer('seq','Greens',20);
c.Trace = cmap(10,:);
c.Trace2 = cmap(15,:);
    

    
% Ax1 - Raw Frame

imagesc(ax1, Frame)
colormap('gray')
caxis(ax1,[0 1])
set(ax1, 'Visible','off')

% Preprocessing

imagesc(ax2, Frame)
colormap(ax2, 'gray')
caxis(ax2,[0 1])
set(ax2, 'Visible', 'off')

img = imagesc(ax22, Objects);
colormap(ax22, c.Objects)
set(ax22, 'Visible','off')
img.AlphaData = 0.5*Objects;

Frame_binary = Frame;
Frame_binary(Frame_binary > Settings.Silhouettethreshold ) = 0;
Frame_binary(find( Frame_binary )) = 1; %#ok<*FNDSB>
Frame_binary(find( Objects )) = 0;

frame_normal = Frame_binary;
frame = imdilate( Frame_binary, strel('diamond',5));
frame = imerode(frame, strel('diamond',45));
mask = imdilate(frame, strel('diamond',40));
frame = frame_normal.*mask;
fe = edge(Frame_binary);
Frame_edge = imdilate(fe, strel('diamond',4));


img = imagesc(ax23, frame);
colormap(ax23, c.Snout)
set(ax23, 'Visible', 'off')
img.AlphaData = 0.5.*frame;

img = imagesc(ax24, Edges);
colormap(ax24, c.Objects)
set(ax24, 'Visible', 'off');
img.AlphaData = Edges;


% Nose Tracking
imagesc(ax3, Frame)
colormap(ax3, 'gray')
caxis(ax3, [0 5])
set(ax3, 'Visible' ,'off')

zoomheigth = 200;
zoomwidth = zoomheigth*(Settings.Video_heigth/Settings.Video_width);
X0 = 150;
Y0 = 250;
xlim(ax3, [X0 X0+zoomwidth])
ylim(ax3, [Y0 Y0+zoomheigth])

img = imagesc(ax31, Frame_binary);
colormap(ax31, c.Snout)
set(ax31, 'Visible','off')
img.AlphaData = 0.5*Frame_binary;
xlim(ax31,[X0 X0+zoomwidth])
ylim(ax31, [Y0 Y0+zoomheigth])

img = imagesc(ax30, Frame_edge);
colormap(ax30, c.Snout2)
set(ax30, 'Visible','off')
img.AlphaData = 0.5*Frame_edge;
xlim(ax30,[X0 X0+zoomwidth])
ylim(ax30, [Y0 Y0+zoomheigth])

area(ax32,1:640,sum(Frame_binary,1),'FaceColor',c.Snout,'FaceAlpha',0.5)
line(ax3, [Nose(2) Nose(2)], [0 Settings.Video_width], 'Color',c.Snout,'LineStyle',':')
set(ax32, 'Visible','off')
xlim(ax32, [X0 X0+zoomwidth])
ylim(ax32, [Y0 Y0+zoomheigth])

area(ax33,1:512,flip(sum(Frame_binary,2)),'FaceColor',c.Snout,'FaceAlpha',0.5)
set(ax33, 'Visible','off')
%ylim(ax33, [0 512])
line(ax3, [0 Settings.Video_heigth], [Nose(1) Nose(1)],...
    'Color', c.Snout, 'LineStyle',':','LineWidth',1)
xlim(ax33, [512-(Y0+zoomheigth) 512-Y0])
hold(ax31, 'on')
scatter(ax31, Nose(2), Nose(1), 'MarkerFaceColor', c.Snout, 'MarkerEdgeColor','k')

theta = 1:1:360;
Cx = round(Nose(2) + 40*sind(theta));
Cx = [Cx, Cx+1];
Cy = round(Nose(1) + 40*cosd(theta));
Cy = [Cy, Cy+1]; 
scatter(ax31, Cx, Cy, 6,'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
scatter(ax31, Cx, Cy, 3,'MarkerFaceColor', c.Nose, 'MarkerEdgeColor', c.Nose)
IDX = find(Cx > 1 & Cx < size(frame,2) & Cy > 1 & Cy< size(frame,1));
C = sub2ind( size(frame), Cy(IDX), Cx(IDX) );
PTS = [];
PTS(:,1) = find( fe( C ));
[PTS(1:length(PTS),2),PTS(1:length(PTS),1)] = ind2sub(size(frame),C(PTS));


% - clean entries to find points on edge only
for id = 1:size(PTS,1)
    dist = [];
    if isnan(PTS(id,1))
        continue
    end
    for ii = 1:size(PTS,1)
        if id == ii || isnan(PTS(ii,1))
            continue
        end
        dist(ii) = finddist(PTS(id,:)',PTS(ii,:));
        if dist(ii) < 25
            PTS(ii,:) = [NaN,NaN];
        end
        
    end
end
PTS = PTS( ~isnan(PTS(:,1)),:);

line(ax31, PTS(1:2,1), PTS(1:2,2), 'Color', 'k','LineWidth',2.5)
line(ax31, PTS(1:2,1), PTS(1:2,2), 'Color', c.Objects(2,:),'LineWidth',2)
scatter(ax31, PTS(:,1), PTS(:,2), 30,'MarkerFaceColor', c.Objects(2,:), 'MarkerEdgeColor', 'k')




% ROI for tracking initiation
imagesc(ax4, Frame)
colormap(ax4, 'gray')
caxis(ax4, [0 1])
set(ax4, 'Visible' ,'off')

xlim(ax4, [X0 X0+zoomwidth])
ylim(ax4, [Y0 Y0+zoomheigth])


SLN = zeros(size(Objects)); % Variable for silhouette (normal shape)
SLN( Frame <= Settings.Silhouettethreshold ) = 1;
SLN( Objects == 1 ) = 0;


SLS = imerode(SLN, strel('diamond', 5)); % Small silhouette
SLL = imdilate(SLS, strel('diamond', Settings.Dilationsize));

K = 5;
L = 1;
CHl = conv2(Frame, ones(K,K)./K^2, 'same'); % Low pass filter frame with large kernel
CHs = conv2(Frame, ones(L,L)./L^2, 'same');  % Low pass filter frame with small kernel
CH = zeros(size(Objects));
CH(0.5*CHs./CHl < 0.485) = 1; % Matrix with potential valid pixels


FrameHeight = size(Frame, 1);
FrameWidth = size(Frame, 2);


Frame(1,1) = 999; % Default 'dead' pixel, to refer to if tracking should
% be interrupted.
CH(Objects == 1) = 0;
CH(Edges == 1) = 0.5;
CH(SLL == 1) = 0;
CH(SLN == 1) = 0;


RoiInd = find(edge(SLL));
f_edge = imdilate(edge(SLL), strel('diamond',2));

img = imagesc(ax40, SLL);
colormap(ax40, c.Snout);
set(ax40, 'Visible','off')
img.AlphaData = 0.5*SLL;
xlim(ax40,[X0 X0+zoomwidth])
ylim(ax40, [Y0 Y0+zoomheigth])

img = imagesc(ax41, f_edge);
colormap(ax41, c.Roi);
set(ax41, 'Visible','off')
img.AlphaData = f_edge;
xlim(ax41,[X0 X0+zoomwidth])
ylim(ax41, [Y0 Y0+zoomheigth])


% Segmentations
cla(ax5)
imagesc(ax5, Frame)
colormap(ax5, 'gray')
caxis(ax5, [0 1])
set(ax5, 'Visible' ,'off')

zoomheigth = 100;
zoomwidth = zoomheigth*(Settings.Video_heigth/Settings.Video_width);
X0 = 250;
Y0 = 280;
xlim(ax5, [X0 X0+zoomwidth])
ylim(ax5, [Y0 Y0+zoomheigth])



set(ax51, 'Visible', 'off')

xlim(ax51, [X0 X0+zoomwidth])
ylim(ax51, [Y0 Y0+zoomheigth])


[RoiSub(:,1), RoiSub(:,2)] = ind2sub(size(Frame), RoiInd);
RS2 = RoiSub;
RS2(:,1) = RS2(:,1)+1;
RoiSub = [RoiSub; RS2];
todo = 1:size(RoiSub,1);
todo(1) = NaN;
did(1:length(RoiInd)) = NaN;
did(1) = 1;
TempVar = zeros(length(RoiInd), 2);
TempVar(1,:) = RoiSub(1,:);

for i = 2:size(RoiSub,1)
    todo = todo(~isnan(todo));
    dist = sqrt(sum( (RoiSub(todo,:) - TempVar(end,:)).^2, 2));
    [~, idx] = min(dist);
    TempVar(i,:) = RoiSub(todo(idx),:);
    did(i) = todo(idx);
    todo(idx) = NaN;
end
TempVar(TempVar(:,1) < 1, :) = NaN;
TempVar(TempVar(:,1) > FrameHeight, :) = NaN;
TempVar(TempVar(:,2) < 1, :) = NaN;
TempVar(TempVar(:,2) > FrameWidth,:) = NaN;
TempVar = TempVar(~isnan(TempVar(:,1)),:);

RoiSub = TempVar;

%scatter(RoiSub(:,2), RoiSub(:,1), 'm','filled')

RoiInd = sub2ind(size(Frame), TempVar(:,1), TempVar(:,2));

SeedProfile = Frame(RoiInd); % Intensity profile over ROI

LP = medfilt1(SeedProfile, 50);
SeedProfile = abs(SeedProfile-LP); % Subtract local background

if isempty(SeedProfile); TracesOut = {}; return; end

[~, SeedInd] = findpeaks(SeedProfile, 'MinPeakDistance',3);
Seeds = RoiSub(SeedInd,:);
d = sqrt( sum(  (Seeds-Nose).^2,2 ));
Seeds = Seeds( d <=150 , :);



img = imagesc(ax51, SLL);
colormap(ax51, c.Snout);
set(ax51, 'Visible','off')
img.AlphaData = 0.5*SLL;
xlim(ax51,[X0 X0+zoomwidth])
ylim(ax51, [Y0 Y0+zoomheigth])

hold(ax51, 'on')
[fa, fb] = find(edge(SLL));
%scatter(ax51, fb, fa, 20, 'MarkerFaceColor', c.Roi, 'MarkerEdgeColor', c.Roi)
scatter(ax51, Seeds(:,2) ,Seeds(:,1), 16, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
scatter(ax51, Seeds(:,2) ,Seeds(:,1), 11, 'MarkerFaceColor', c.Trace, 'MarkerEdgeColor',c.Trace)

hold(ax41, 'on')
scatter(ax41, Seeds(:,2) ,Seeds(:,1), 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
scatter(ax41, Seeds(:,2) ,Seeds(:,1), 5, 'MarkerFaceColor', c.Interest , 'MarkerEdgeColor',c.Interest )

bufferSize = 100; % max nr of points on a trace


Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';
H = -15:15;
HH = -30:30;
trace_finished = zeros(1, size(Seeds,1));
trace_edges = zeros(1, size(Seeds,1));

i = 2;

Theta = [];
if i == 2
    ThetaRoi = repmat(1:20:360, [size(Seeds,1), 1]); % Roi bounded by angle range
    %N = repmat(Centre,[size(Traces,2) 1]);
    %vt = Traces(:,:,i-1) - N';
    %ThetaRoi = atan2d(vt(2,:), vt(1,:))' + [-60:20:60];
else
    vt = Traces(:,:,i-1)- Traces(:,:,i-2);
    ThetaRoi = atan2d(vt(2,:), vt(1,:))' + H;
    
end

Theta(1,:,:) = [3*cosd(ThetaRoi)]; % ROI per previous point
Theta(2,:,:) = [3*sind(ThetaRoi)];

LastPoints = repmat(Traces(:,:,i-1), [1 , 1, size(Theta, 3)]);

RoiSub = round(LastPoints + Theta);

% Remove indices outside of frame range
RoiSub(RoiSub < 1) = 1;
RoiSub(1,RoiSub(1,:,:) > FrameHeight) = 1;
RoiSub(2,RoiSub(2,:,:) > FrameWidth) = 1;
RoiInd = squeeze(sub2ind(size(Frame), RoiSub(1,:,:), RoiSub(2,:,:)));


[a,b] = find(SLL(RoiInd));
%RoiInd(SLL(RoiInd) == 1) = 1;
for j = 1:size(RoiInd, 1)
    
    
    idx = b(a == j);
    if ~any(diff(idx)>1) & min(idx) > 2 & max(idx) < size(RoiInd,2)-1
        idx = [[-2 -1]'+idx(1); idx; [1 2]'+idx(end)];
    elseif ~any(diff(idx) > 1) & min(idx) == 2 & max(idx) < size(RoiInd, 2) - 1
        idx = [-1+idx(1); idx; [1 2]'+idx(end); size(RoiInd,2)];
    elseif ~any(diff(idx) > 1) & min(idx) == 1 & max(idx) < size(RoiInd, 2) - 1
        idx = [idx; [1 2]'+idx(end);[-1 0]'+size(RoiInd,2)];
        
    elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(RoiInd,2)
        idx = [ [1 2]'; [-2 -1]'+idx(1); idx];
    elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(RoiInd,2)-1
        idx = [ [1]'; [-2 -1]'+idx(1); idx; size(RoiInd,2)];
        
    elseif numel(find(diff(idx) > 1)) == 1 & min(idx) == 1 & max(idx) == size(RoiInd,2)
        id = find(diff(idx) > 1);
        idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
        
    elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) < size(RoiInd,2) - 1
        idx = [ [-2 -1]'+min(idx); idx ; [1 2]'+idx(end)];
        
    elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) == size(RoiInd,2) - 1
        idx = [1; [-2 -1]'+min(idx); idx ; [1 ]'+idx(end)];
        
        
    elseif numel(find(diff(idx) > 1)) > 1 & min(idx) == 1 & max(idx) == size(RoiInd,2)
        [~, id] = max(diff(idx));
        idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
        
        
    else
        idx =idx;
        
        
    end
    
    RoiInd(j,idx) = 1;
    
end

for j = 1:size(RoiInd, 1)
    [e,f] = ind2sub(size(Frame), RoiInd(j,:));
    
    scatter(ax51,f,e , 5, 'MarkerFaceColor', c.Roi, 'MarkerEdgeColor', c.Roi)
end

Iroi = Frame(RoiInd);
if size(Seeds,1) == 1
    [~, PointSub] = min(Iroi);
    PointInd = PointSub;
else
    [~, PointSub] = min(Iroi, [], 2);
    PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
end


[Traces(1, :, i), Traces(2,:,i)] = ind2sub(size(Frame), RoiInd(PointInd));
             
for j = 1:size(Traces, 2)
    scatter(ax51, Traces(2,j,i), Traces(1,j,i), 16, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    scatter(ax51, Traces(2,j,i), Traces(1,j,i), 11, 'MarkerFaceColor', c.Interest, 'MarkerEdgeColor', c.Interest)
       
end
%%

zoomheigth = 70;
zoomwidth = zoomheigth*(Settings.Video_heigth/Settings.Video_width);
X0 = 230;
Y0 = 420;
imagesc(ax7, Frame)
colormap(ax7, 'gray')
caxis(ax7, [0 1])
set(ax7, 'Visible', 'off')
xlim(ax7,[X0 X0+zoomwidth])
ylim(ax7, [Y0 Y0+zoomheigth])
hold(ax7, 'on')

img = imagesc(ax71, Edges);
colormap(ax71, c.Objects);
img.AlphaData = Edges;
set(ax71, 'Visible', 'off')
xlim(ax71,[X0 X0+zoomwidth])
ylim(ax71, [Y0 Y0+zoomheigth])
%%



% Segmentations

imagesc(ax6, CH)
colormap(ax6, 'gray')
caxis(ax6, [0 1])
set(ax6, 'Visible' ,'off')

zoomheigth = 30;
zoomwidth = zoomheigth*(Settings.Video_heigth/Settings.Video_width);
X0 = 300;
Y0 = 290;
xlim(ax6,[X0 X0+zoomwidth])
ylim(ax6, [Y0 Y0+zoomheigth])

hold(ax6, 'on')
alf = 1;
didinterpol = 0;
for i = 2:bufferSize
    
    if ~any(any(Traces(:,:, i-1) ~= 1)) % If no points were assigned in previous loop
        break
    end
    
    Theta = [];
    if i == 2
        ThetaRoi = repmat(1:20:360, [size(Seeds,1), 1]); % Roi bounded by angle range
        %N = repmat(Centre,[size(Traces,2) 1]);
        %vt = Traces(:,:,i-1) - N';
        %ThetaRoi = atan2d(vt(2,:), vt(1,:))' + [-60:20:60];
    else
        vt = Traces(:,:,i-1)- Traces(:,:,i-2);
        
        for qq = 1:size(vt, 2)
            if i > 2 
                p1 = Traces(:,qq,i-1);
                p2 = p1 + vt(:,qq); 
           
                line(ax6, [p1(2) p2(2)], [p1(1) p2(1)], 'Color', [c.Roi2 alf],'LineWidth',3)
            end
        end
        ThetaRoi = atan2d(vt(2,:), vt(1,:))' + H;
        
    end
    
    
    
    Theta(1,:,:) = [3*cosd(ThetaRoi)]; % ROI per previous point
    Theta(2,:,:) = [3*sind(ThetaRoi)];
    
    LastPoints = repmat(Traces(:,:,i-1), [1 , 1, size(Theta, 3)]);
    
    RoiSub = round(LastPoints + Theta);
    
    % Remove indices outside of frame range
    RoiSub(RoiSub < 1) = 1;
    RoiSub(1,RoiSub(1,:,:) > FrameHeight) = 1;
    RoiSub(2,RoiSub(2,:,:) > FrameWidth) = 1;
    RoiInd = squeeze(sub2ind(size(Frame), RoiSub(1,:,:), RoiSub(2,:,:)));
    
    n_cut = 2;
    if i == 2
        [a,b] = find(SLL(RoiInd));
        %RoiInd(SLL(RoiInd) == 1) = 1;
        for j = 1:size(RoiInd, 1)
            
            
            idx = b(a == j);
            if ~any(diff(idx)>1) & min(idx) > 2 & max(idx) < size(RoiInd,2)-1
                idx = [[-2 -1]'+idx(1); idx; [1 2]'+idx(end)];
            elseif ~any(diff(idx) > 1) & min(idx) == 2 & max(idx) < size(RoiInd, 2) - 1
                idx = [-1+idx(1); idx; [1 2]'+idx(end); size(RoiInd,2)];
            elseif ~any(diff(idx) > 1) & min(idx) == 1 & max(idx) < size(RoiInd, 2) - 1
                idx = [idx; [1 2]'+idx(end);[-1 0]'+size(RoiInd,2)];
                
            elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(RoiInd,2)
                idx = [ [1 2]'; [-2 -1]'+idx(1); idx];
            elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(RoiInd,2)-1
                idx = [ [1]'; [-2 -1]'+idx(1); idx; size(RoiInd,2)];
                
            elseif numel(find(diff(idx) > 1)) == 1 & min(idx) == 1 & max(idx) == size(RoiInd,2)
                id = find(diff(idx) > 1);
                idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
                
            elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) < size(RoiInd,2) - 1
                idx = [ [-2 -1]'+min(idx); idx ; [1 2]'+idx(end)];
                
            elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) == size(RoiInd,2) - 1
                idx = [1; [-2 -1]'+min(idx); idx ; [1 ]'+idx(end)];
                
                
            elseif numel(find(diff(idx) > 1)) > 1 & min(idx) == 1 & max(idx) == size(RoiInd,2)
                [~, id] = max(diff(idx));
                idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
                
                
            else
                idx =idx;
                
                
            end
            
            RoiInd(j,idx) = 1;
            
            
            
            
            
        end
        
    end
    
    if i > 2
        for qq = 1:size(RoiInd, 1)          
       
            
            [a,b] = ind2sub(size(Frame), RoiInd(qq,:));
            scatter(ax6, b, a, 32, 'MarkerFaceColor', c.Roi, 'MarkerEdgeColor', c.Roi,...
                'MarkerFaceAlpha',alf,'MarkerEdgeAlpha',alf)
        end
    end
    
    
    Iroi = Frame(RoiInd);
    if size(Seeds,1) == 1
        [~, PointSub] = min(Iroi);
        PointInd = PointSub;
    else
        [~, PointSub] = min(Iroi, [], 2);
        PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
    end
    
    
    [Traces(1, :, i), Traces(2,:,i)] = ind2sub(size(Frame), RoiInd(PointInd));
    
    trace_finished(CH(RoiInd(PointInd)) == 0) = trace_finished(CH(RoiInd(PointInd)) == 0) + 1;
    trace_edges(CH(RoiInd(PointInd)) == 0.5) = trace_edges(CH(RoiInd(PointInd)) == 0.5) + 1;

    
    if i < 3
        trace_finished(CH(RoiInd(PointInd)) == 0.5) =  trace_finished(CH(RoiInd(PointInd)) == 0.5) + 1;
    elseif any(CH(RoiInd(PointInd)) == 0.5)
        nan_traces = find(CH(RoiInd(PointInd)) == 0.5);
        
       
        
        for j = 1:length(nan_traces)
            
          
            t = squeeze(Traces(:, nan_traces(j), :));
            if didinterpol == 0
            scatter(ax7, t(2,:), t(1,:), 32, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')            
            scatter(ax7, t(2,:), t(1,:), 16, 'MarkerFaceColor', c.Trace, 'MarkerEdgeColor', c.Trace)
            end
            i1 = find(t(1,:) == 1, 1, 'first');
            i2 = find(t(2,:) == 1, 1, 'first');
            if ~isempty([i1 i2])
                t = t(:,1:max([i1 i2])-1);
            end
            
            if i < 10
                rawax = 1:i;
                fitax = 1:0.5:i+10;
            else
                rawax = 1:10;
                fitax = 1:0.5:i+10;
            end
            
            px = polyfit(rawax, t(1,end-length(rawax)+1:end), 1);
            py = polyfit(rawax, t(2,end-length(rawax)+1:end), 1);
            
            tfit = [];
            tfit(1,:) = round(polyval(px, fitax));
            tfit(2,:) = round(polyval(py, fitax));
            
            throwidx = [find(tfit(1,:) < 1), find(tfit(1,:) > FrameHeight)];
            throwidx = [throwidx, find(tfit(2,:) < 1), find(tfit(2,:) > FrameWidth)];
            
            if ~isempty(throwidx)
                tfit(:, throwidx) = NaN;
                tfit = tfit(:, ~isnan(tfit(1,:)));
            end
            
            tfitInd = sub2ind(size(Frame), tfit(1,:), tfit(2,:));
            ptInd = find(CH(tfitInd) == 0.5,1,'last');
            
            if didinterpol == 0
                hold(ax71, 'on')
                plot(ax71, tfit(2,:), tfit(1,:), 'color',c.Roi2,'LineWidth',2)
                scatter(ax71, tfit(2, ptInd), tfit(1, ptInd), 32,'MarkerFaceColor','k','MarkerEdgeColor','k')
                scatter(ax71, tfit(2, ptInd), tfit(1, ptInd), 16,'MarkerFaceColor',c.Interest,'MarkerEdgeColor',c.Interest)
            end
            
            if ptInd < size(tfit,2)-5
                LastPt = tfit(:, ptInd+1);
                %ThetaLoop = squeeze(Theta(:, nan_traces(j), :));
                
                vt = tfit(:, end) - tfit(:, end-5);
                T1 = atan2d(vt(2,1), vt(1,1)) + HH;
                T2 = [];
                T2(1,:) = 6*cosd(T1)';
                T2(2,:) = 6*sind(T1)';
                
                
                
                RoiSub = round(LastPt + T2);
                RoiSub(RoiSub < 1) = 1;
                RoiSub(1,RoiSub(1,:) > FrameHeight) = 1;
                RoiSub(2,RoiSub(2,:) > FrameWidth) = 1;
                RoiInd = squeeze(sub2ind(size(Frame), RoiSub(1,:), RoiSub(2,:)));
                
                scatter(ax71, RoiSub(2,:), RoiSub(1,:), 'MarkerFaceColor', c.Roi, 'MarkerEdgeColor', c.Roi)
                Iroi = Frame(RoiInd);
                
                [~, PointSub] = min(Iroi, [], 2);
                PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
                [a,b] = ind2sub(size(Frame), RoiInd(PointInd));
                
                if didinterpol == 0
                    scatter(ax71, b, a, 32, 'MarkerFaceColor','k','MarkerEdgeColor','k')
                    scatter(ax71, b, a, 16, 'MarkerFaceColor',c.Interest,'MarkerEdgeColor',c.Interest)
                end
                
                didinterpol = 1;
                if CH(a,b) == 1
                    Traces(1, nan_traces(j), i) = a;
                    Traces(2, nan_traces(j), i) = b;
                    trace_edges(nan_traces(j)) = 0;
                    
                else
                    trace_finished( nan_traces(j) ) = trace_finished( nan_traces(j) )  + 1;
                end
                
            end
        end
        
    end
    
    
    Traces(1:2,trace_finished>2,i) = 1;
    
    idx = find(trace_edges > 10 & Traces(1,:,i)>1);
    if ~isempty(idx)
        for j = 1:length(idx)         
            n_delete = trace_edges(idx(j));
            Traces(:, idx(j), i-n_delete:i) = 1;
        end
    end
    
  
  

end



%plot(ax6, squeeze(Traces(2,105,1:9)), squeeze(Traces(1,105,1:9)), 'color',c.Red,'LineWidth',2)

for j = 1:size(Traces,2)
  
    scatter(ax6, squeeze(Traces(2,j,:)), squeeze(Traces(1,j,:)), 32,'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor','k','MarkerFaceAlpha',alf,'MarkerEdgeAlpha',alf)
    scatter(ax6, squeeze(Traces(2,j,:)), squeeze(Traces(1,j,:)), 25,'MarkerFaceColor', c.Trace, ...
        'MarkerEdgeColor', c.Trace,'MarkerFaceAlpha',alf,'MarkerEdgeAlpha',alf)

end



Traces(Traces == 1) = NaN;



TOut = {};
TracesOut = {};
Or = [];
for i = 1:size(Traces,2)
    t = squeeze(Traces(:,i,:));
    idx = find(~isnan(t(1,:)),1,'last');
    t = t(:, 1:idx);
    if size(t,2) > 10
        TOut{end+1} = t';
        Or(end+1,:) = t(:,1)';
    end  
    
    
end


if ~isempty(Or)
    dx = abs(Or(:,1) - Or(:,1)');
    dy = abs(Or(:,2) - Or(:,2)');
    d = sqrt( (dx+dy).^2 );



    flag = ones(1, size(TOut, 2));
    for i = 1:length(flag)
        d(i,i) = NaN;
        idx = find(d(i,:) < 10);
        if ~isempty(idx)
            t1 = TOut{i};

            for j = 1:length(idx)
                t2 = TOut{idx(j)};
               dx = abs(t1(:,1) - t2(:,1)');
                dy = abs(t1(:,2) - t2(:,2)');
                dnew = sqrt( (dx+dy).^2);
                val = min(dnew,[],2);
                dist = mean(min(dnew,[],2));

                if dist < 3
                    flag(i) = 0;
                    d(:,i) = NaN;
                end    
            end
        end
    end


    for i = 1:length(flag)
        if flag(i)
            TracesOut{end+1} = TOut{i};
        end
    end
end






%%
cla(ax8)
imagesc(ax8, Frame)
colormap(ax8, 'gray')
caxis([0 1])
hold(ax8, 'on')
set(ax8, 'Visible', 'off')
for i = 1:size(TracesOut, 2)
    t = TracesOut{i};
    plot(ax8, t(:,2), t(:,1),'color',c.Trace2)  
end

zoomheigth = 400;
zoomwidth = zoomheigth*(Settings.Video_heigth/Settings.Video_width);
X0 = 0;
Y0 = 100;
xlim(ax8,[X0 X0+zoomwidth])
ylim(ax8, [Y0 Y0+zoomheigth])

% Objects interaction









