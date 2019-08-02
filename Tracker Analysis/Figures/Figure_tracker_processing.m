clear
clc
close all
warning('off')

video = 'F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02.dat';
comp = 'F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_compiled.mat';

load(comp)
Settings = Annotations.Settings;
Settings.Video(1) = 'F';
cmap = makeColor();

nframes = Settings.Nframes;
idx = round(linspace(1,nframes,100));
Frames = zeros(512, 640, length(idx));
for i = 1:length(idx)
    Settings.Current_frame = idx(i)-1;
    Frames(:,:,i) = LoadFrame(Settings);
end
sumFrames = sum(Frames,3);
meanFrames = mean(Frames,3);


Settings.Current_frame = 300;



%%
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
% ax24 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
x = x + x_space+ ax_width;
ax3 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax30 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
% ax31 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
% ax32 = axes(f,'Units','points','Position',[x y ax_width 40]);
% ax33 = axes(f,'Units','points','Position',[x y 20 ax_heigth]);


x = marg_x;
y = y - ax_heigth - y_space;
ax4 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax40 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
% ax41 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
x = x + x_space + ax_width;

ax5 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
%ax51 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);

x = x + x_space + ax_width;
ax6 = axes(f, 'Units', 'points', 'Position', [x y ax_width ax_heigth]);
ax61 = axes(f, 'Units', 'points', 'Position', [x y ax_width ax_heigth]);

y = y - ax_heigth - y_space;
x = marg_x;
ax7 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax71 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);

x = x + x_space + ax_width;
ax8 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);
ax81 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);


x = x + x_space + ax_width;
ax9 = axes(f, 'Units', 'points', 'Position', [x y ax_width ax_heigth]);



%% Axes 1 - raw frame
Frame = LoadFrame(Settings);
Frame = im2double(abs(Frame));

if Settings.doGaussian
    ng = Settings.Gaussian_kernel_size;
    kfil = zeros(1, ng);
    kfil(1:ceil(ng/2)) = 1:ceil(ng/2);
    kfil(ceil(ng/2)+1:end) = max(kfil)-1:-1:1;
    Kgaus = repmat(kfil,[ng, 1]) + repmat(kfil,[ng, 1])';
    Kgaus = Kgaus./ sum(sum(Kgaus));
    Frame = conv2(Frame, Kgaus, 'same');
end



Frame = adapthisteq(Frame, 'NBins',256,'NumTiles',[5 5]);
minval = min(min(Frame));
Frame = Frame - minval;
Frame = Frame./max(max(Frame));

imagesc(ax1, Frame)
set(ax1, 'Visible','off')
colormap(ax1, 'gray')

rectangle(ax1, 'Position',[0 0 70 70], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax1, 15, 35, 'A','FontSize',16)


%% Axes 2 - Background extraction
Background = sumFrames./100;

Background = Background - min(min(Background));
Background = Background ./ max(max(Background));


Background(Background > Settings.Background_threshold) = 0;
Objects = zeros(size(Background));

Objects(Background > 0) = 1;
KL = conv2(meanFrames, ...
    ones(Settings.Edges_kernel_large,Settings.Edges_kernel_large)./...
    Settings.Edges_kernel_large,'same');
KS = conv2(meanFrames, ...
    ones(Settings.Edges_kernel_small, Settings.Edges_kernel_small)./...
    Settings.Edges_kernel_small^2, 'same');

H = KS./KL;
H = H-min(min(H));
H = H./max(max(H));

tax = 0:0.01:1;
for i = 1:length(tax)
    count(i) = numel(find(H < tax(i))); %#ok<AGROW>
end
[~, idx] = max(diff(count));
edge_threshold = tax(idx)-0.05;

Edges = zeros(size(Objects));

%Edges(H < Settings.Edges_threshold) = 1;
Edges(H < edge_threshold) = 1;


Res = Costumbackground(Objects, Edges);
Objects = Res.Objects;
Edges = Res.Edges;



imagesc(ax2, sumFrames./100)
colormap(ax2, 'gray')
set(ax2, 'Visible', 'off')

img = imagesc(ax22, Objects);
img.AlphaData = 0.5.*Objects;
colormap(ax22, cmap.Objects);
set(ax22, 'Visible', 'off')

img = imagesc(ax23, Edges);
img.AlphaData = 1.*Edges;
colormap(ax23, cmap.Objects)
set(ax23, 'Visible', 'off')

rectangle(ax23, 'Position',[0 0 70 70], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax23, 15, 35, 'B','FontSize',16)

%% Axes 3 - Shape
SLN = zeros(size(Objects));
SLN( Frame <= Settings.Shape_threshold) = 1;
SLN(Objects == 1) = 0;

SLS = imerode(SLN, strel('diamond', 5)); % Small silhouette
SLL = imdilate(SLS, strel('diamond', Settings.Dilation));


imagesc(ax3, Frame)
colormap(ax3, 'gray')
set(ax3, 'Visible', 'off')

img = imagesc(ax30, SLL);
colormap(ax30, cmap.Snout)
set(ax30, 'Visible', 'off')
img.AlphaData = 0.5.*SLL;

rectangle(ax30, 'Position',[0 0 70 70], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax30, 15, 35, 'C','FontSize',16)

%% Axes 4 - Nose
cla(ax40)
FM = imdilate(SLN, strel('diamond', 5));
FM = imerode(FM, strel('diamond', 45));
FM = imdilate(FM, strel('diamond', 55));

FrameNose = SLN.*FM;
f_sum = sum(FrameNose ,2)./ size(FrameNose, 2);
X = find( f_sum > 0, 1, 'last');
Y = round( mean( find( FrameNose(X,:)), 'omitnan' ));


FrameNose2 = edge(FrameNose);
theta = 0:1:360;
% Get headangle:
% - find non-zero entries in circle around nose
Cx = round(Y + 40*sind(theta));
Cx = [Cx, Cx+1]; %#ok<AGROW>
Cy = round(X + 40*cosd(theta));
Cy = [Cy, Cy+1]; %#ok<AGROW>
IDX = find(Cx > 1 & Cx < size(FrameNose2,2) & Cy > 1 & Cy< size(FrameNose2,1));

C = sub2ind( size(FrameNose2), Cy(IDX), Cx(IDX) );
PTS = [];
PTS(:,1) = find( FrameNose2( C ));
[PTS(1:length(PTS),2),PTS(1:length(PTS),1)] = ind2sub(size(FrameNose2),C(PTS));


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

PTS = PTS(1:2,:); % Lazy
Vp = PTS(2,:) - PTS(1,:);

% Headangle described as vector, normal to Vp
AngleVector(i,1:2) = [-Vp(1),Vp(2)];

% Normalise
AngleVector(i,1:2) = AngleVector(i,1:2)...
    ./sqrt(sum(AngleVector(i,1:2).^2));


if AngleVector(i,1) > 0
    AngleVector(i,:) = -AngleVector(i,:);
end




imagesc(ax4, Frame)
colormap(ax4, 'gray')
set(ax4, 'Visible', 'off')
hold(ax4, 'on')

img = imagesc(ax40, FrameNose);
colormap(ax40, cmap.Snout)
img.AlphaData = 0.5.*FrameNose;
set(ax40, 'Visible', 'off')
hold(ax40, 'on')
scatter(ax40,Y,X, 'MarkerFaceColor', cmap.Nose, 'MarkerEdgeColor', 'k')
scatter(ax40, Cx,Cy, 1, 'MarkerFaceColor', cmap.Roi,'MarkerEdgeColor',cmap.Roi)
line(ax40, PTS(2,:), PTS(1,:), 'color', cmap.Interest)

scatter(ax40, PTS(:,1), PTS(:,2), 20, 'MarkerFaceColor', cmap.Interest, 'MarkerEdgeColor','k')

zoom_heigth = 340;
zoom_width = zoom_heigth*(Settings.Video_heigth/Settings.Video_width);
zoom_base_x = 60;
zoom_base_y = 120;
xlim(ax4, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax4, [zoom_base_y zoom_base_y+zoom_heigth])
xlim(ax40, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax40, [zoom_base_y zoom_base_y+zoom_heigth])

rectangle(ax40, 'Position',[0+zoom_base_x 0+zoom_base_y 50 50], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax40, 10+zoom_base_x, 25+zoom_base_y, 'D','FontSize',16)

%% Tracker threshold

K =  Settings.Trace_kernel_large;
L =  Settings.Trace_kernel_small;
CHl = conv2(Frame, ones(K,K)./K^2, 'same'); % Low pass filter frame with large kernel
CHs = conv2(Frame, ones(L,L)./L^2, 'same');  % Low pass filter frame with small kernel
CH = zeros(size(Objects));
CH(0.5*CHs./CHl < Settings.Trace_threshold) = 1; % Matrix with potential valid pixels
CH(Objects == 1) = 0;
CH(Edges == 1) = 0.5;
CH(SLL == 1) = 0;
CH(SLN == 1) = 0;

imagesc(ax5, CH)
colormap(ax5, 'gray')
set(ax5, 'Visible', 'off')


rectangle(ax5, 'Position',[0 0 70 70], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax5, 15, 35, 'E','FontSize',16)





%% SEED roi





% Find tracing seeds
RoiInd = find(edge(SLL));

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
FrameHeight = size(Frame, 1);
FrameWidth = size(Frame, 2);


TempVar(TempVar(:,1) < 1, :) = NaN;
TempVar(TempVar(:,1) > FrameHeight, :) = NaN;
TempVar(TempVar(:,2) < 1, :) = NaN;
TempVar(TempVar(:,2) > FrameWidth,:) = NaN;
TempVar = TempVar(~isnan(TempVar(:,1)),:);

RoiSub = TempVar;

%scatter(RoiSub(:,2), RoiSub(:,1), 'm','filled')

RoiInd = sub2ind(size(Frame), TempVar(:,1), TempVar(:,2));

SeedProfile = double(Frame(RoiInd)); % Intensity profile over ROI

LP = medfilt1(SeedProfile, 50);
SeedProfile = abs(SeedProfile-LP); % Subtract local background



Seeds = RoiSub(find(CH(RoiInd)==1),:);



imagesc(ax6, Frame)
colormap(ax6, 'gray')
set(ax6, 'Visible' ,'off')

img = imagesc(ax61, SLL);
colormap(ax61, cmap.Snout)
set(ax61, 'Visible', 'off')
img.AlphaData = 0.3*SLL;
hold(ax61, 'on')

scatter(ax61, RoiSub(:,2), RoiSub(:,1), 5, 'MarkerFaceColor', cmap.Roi, 'MarkerEdgeColor', cmap.Roi)
scatter(ax61, Seeds(:,2), Seeds(:,1), 10, 'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', cmap.Trace)



rectangle(ax61, 'Position',[0 0 70 70], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax61, 15, 35, 'F','FontSize',16)



%% Tracker initiation

cla(ax7)
cla(ax71)

imagesc(ax7, Frame)
colormap(ax7, 'gray')
set(ax7, 'Visible', 'off')

hold(ax7, 'on')

img = imagesc(ax71, SLL);
colormap(ax71, cmap.Snout);
set(ax71, 'Visible', 'off')
img.AlphaData = 0.3*SLL;
hold(ax71, 'on')

zoom_heigth = 15;
zoom_width = zoom_heigth*(Settings.Video_heigth/Settings.Video_width);
zoom_base_x = 307;
zoom_base_y = 288;
xlim(ax7, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax7, [zoom_base_y zoom_base_y+zoom_heigth])
xlim(ax71, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax71, [zoom_base_y zoom_base_y+zoom_heigth])


TraceToPlot = 126;

% Track Traces
bufferSize = 100; % max nr of points on a trace


Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';
H = -15:15;
HH = -30:30;
trace_finished = zeros(1, size(Seeds,1));
trace_edges = zeros(1, size(Seeds,1));

scatter(ax71, Seeds(TraceToPlot,2), Seeds(TraceToPlot,1), 15, ...
    'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', cmap.Trace)


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
        ThetaRoi = atan2d(vt(2,:), vt(1,:))' + H;
        
        
        if i == 3
            
            
            qq = vt(:, TraceToPlot);
            
            line(ax71, [Traces(2, TraceToPlot, i-1), Traces(2, TraceToPlot,i-1)+qq(2)],...
                [Traces(1, TraceToPlot, i-1) Traces(1, TraceToPlot,i-1)+qq(1)],...
                'color', cmap.Interest,'LineWidth',3)
            scatter(ax71, Traces(2, TraceToPlot, i-1), Traces(1, TraceToPlot, i-1), 20, ...
                'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', cmap.Trace)
        elseif i ==  4
   
            line(ax71, [Traces(2, TraceToPlot, i-1), Traces(2, TraceToPlot,i-1)+qq(2)],...
                [Traces(1, TraceToPlot, i-1) Traces(1, TraceToPlot,i-1)+qq(1)],...
                'color', cmap.Interest,'LineWidth',3,'LineStyle',':')
            scatter(ax71, Traces(2, TraceToPlot, i-1), Traces(1, TraceToPlot, i-1), 20, ...
                'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', cmap.Trace)
        end
    end
    
    
    
    Theta(1,:,:) = [3*cosd(ThetaRoi)]; % ROI per previous point
    Theta(2,:,:) = [3*sind(ThetaRoi)];
    
    
    
    LastPoints = repmat(Traces(:,:,i-1), [1 , 1, size(Theta, 3)]);
    
    RoiSub = round(LastPoints + Theta);
    
    if i < 5
        scatter(ax7, RoiSub(2, TraceToPlot, :), RoiSub(1, TraceToPlot, :), 15,...
            'MarkerFaceColor', 'k' , 'MarkerEdgeColor', 'k','Clipping','on')
    end
    
    
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
    Roisub2 = [];
    [Roisub2(:,1), Roisub2(:,2)] = ind2sub(size(Frame), RoiInd(TraceToPlot,:));
    
    if i < 5
        scatter(ax7, Roisub2(:,2), Roisub2(:,1), 35,...
            'MarkerFaceColor', cmap.Roi , 'MarkerEdgeColor', cmap.Roi2)
    end
    
    Iroi = Frame(RoiInd);
    Iroi(CH(RoiInd) == 0) = NaN;
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
                Iroi = Frame(RoiInd);
                
                [~, PointSub] = min(Iroi, [], 2);
                PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
                [a,b] = ind2sub(size(Frame), RoiInd(PointInd));
                
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
    

    
    if i < 4
        scatter(ax71, Traces(2, TraceToPlot, i), Traces(1, TraceToPlot, i), 20,...
            'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', cmap.Trace)
    end
    
    
end
Traces(Traces == 1) = NaN;

rectangle(ax71, 'Position',[0+zoom_base_x 0+zoom_base_y 2.2 2.2], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax71, zoom_base_x+0.4, zoom_base_y+1, 'G','FontSize',16)


%% Axes 8 - Edge Interaction




TraceToPlot = 109;


% Track Traces
bufferSize = 100; % max nr of points on a trace


Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';
H = -15:15;
HH = -30:30;
trace_finished = zeros(1, size(Seeds,1));
trace_edges = zeros(1, size(Seeds,1));

scatter(ax71, Seeds(TraceToPlot,2), Seeds(TraceToPlot,1), 15, ...
    'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', cmap.Trace)


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
    Roisub2 = [];
    [Roisub2(:,1), Roisub2(:,2)] = ind2sub(size(Frame), RoiInd(TraceToPlot,:));
    

    
    Iroi = Frame(RoiInd);
    Iroi(CH(RoiInd) == 0) = NaN;
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
            i1 = find(t(1,:) == 1, 1, 'first');
            i2 = find(t(2,:) == 1, 1, 'first');
            if ~isempty([i1 i2])
                t = t(:,1:max([i1 i2])-1);
            end
            
            if nan_traces(j) == TraceToPlot
                printdata.t = t;
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
            
            if nan_traces(j) == TraceToPlot
               printdata.tfit = tfit; 
               
            end
            
            throwidx = [find(tfit(1,:) < 1), find(tfit(1,:) > FrameHeight)];
            throwidx = [throwidx, find(tfit(2,:) < 1), find(tfit(2,:) > FrameWidth)];
            
            if ~isempty(throwidx)
                tfit(:, throwidx) = NaN;
                tfit = tfit(:, ~isnan(tfit(1,:)));
            end
            
            tfitInd = sub2ind(size(Frame), tfit(1,:), tfit(2,:));
            ptInd = find(CH(tfitInd) == 0.5,1,'last');
            
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
                Iroi = Frame(RoiInd);
                
             
                
                [~, PointSub] = min(Iroi, [], 2);
                PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
                [a,b] = ind2sub(size(Frame), RoiInd(PointInd));
                
                   if nan_traces(j) == TraceToPlot
                    printdata.ptInd = ptInd;
                    printdata.Roisub = RoiSub;
                    printdata.pt = [a,b];
                end
                
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
Traces(Traces == 1) = NaN;


%%


cla(ax8)
cla(ax81)
imagesc(ax8, Frame)
colormap(ax8, 'gray')
set(ax8, 'Visible', 'off')
hold(ax8, 'on')

img = imagesc(ax81, Edges);
colormap(ax81, cmap.Objects)
set(ax81, 'Visible', 'off')
img.AlphaData = 0.6.*Edges;
hold(ax81, 'on')

plot(ax81, printdata.tfit(2,[1 20]), printdata.tfit(1,[1 20]), 'color', cmap.Interest,...
    'LineStyle', '-','LineWidth',2)
plot(ax81, printdata.tfit(2,[1 50]), printdata.tfit(1,[1 50]), 'color', cmap.Interest,...
    'LineStyle', '--','LineWidth',2)

t = squeeze(Traces(:, TraceToPlot, :));
s= scatter(ax81, t(2,:), t(1,:), 20, 'MarkerFaceColor', [cmap.Trace], ...
    'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9);

scatter(ax81, printdata.t(2,:), printdata.t(1,:), 40, ...
    'MarkerFaceColor', cmap.Trace, 'MarkerEdgeColor', 'k')
scatter(ax81, printdata.t(2,end), printdata.t(1,end), 40, ...
    'MarkerFaceColor', cmap.Roi, 'MarkerEdgeColor', 'k')


scatter(ax81, printdata.Roisub(2,:), printdata.Roisub(1,:), 20, ...
    'MarkerFaceColor', cmap.Roi, 'MarkerEdgeColor', cmap.Roi2)

scatter(ax81, printdata.pt(2), printdata.pt(1), 40, 'MarkerFaceColor', cmap.Trace, ...
    'MarkerEdgeColor', 'k')


zoom_heigth = 40;
zoom_width = zoom_heigth*(Settings.Video_heigth/Settings.Video_width);
zoom_base_x = 235;
zoom_base_y = 450;
xlim(ax8, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax8, [zoom_base_y zoom_base_y+zoom_heigth])
xlim(ax81, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax81, [zoom_base_y zoom_base_y+zoom_heigth])


rectangle(ax81, 'Position',[0+zoom_base_x 0+zoom_base_y 5.8 5.8], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax81, zoom_base_x+1.2, zoom_base_y+2.2, 'H','FontSize',16)

%% full results



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

cla(ax9)
imagesc(ax9, Frame)
set(ax9, 'Visible', 'off')
colormap(ax9, 'gray')
hold(ax9, 'on')
T  = Annotations.Tracker.Traces_clean{Settings.Current_frame};
for i = 1:size(T, 2)
    t = T{i};
    plot(ax9, t(:,2), t(:,1), 'color', cmap.Trace2)
end


zoom_heigth = 400;
zoom_width = zoom_heigth*(Settings.Video_heigth/Settings.Video_width);
zoom_base_x = 30;
zoom_base_y = 100;
xlim(ax9, [zoom_base_x zoom_base_x+zoom_width])
ylim(ax9, [zoom_base_y zoom_base_y+zoom_heigth])

pos = [140 180 40 40];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor', cmap.Objects(2,:))
pos = [360 225 40 40];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor', cmap.Objects(2,:))
pos = [257 385 40 40];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor', cmap.Objects(2,:))
pos = [142 314 40 40];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor', cmap.Objects(2,:))




rectangle(ax9, 'Position',[0+zoom_base_x 0+zoom_base_y 58 58], 'FaceColor', cmap.Gray,'EdgeColor','none')
text(ax9, zoom_base_x+22, zoom_base_y+22, 'I','FontSize',16)







