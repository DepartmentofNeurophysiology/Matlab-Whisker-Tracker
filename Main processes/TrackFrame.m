function Traces = TrackFrame(Settings, Output)
%% Traces = TrackFrame(Settings, Objects)
% Track raw traces in a frame using user defined traces


Objects = Output.Objects;
Edges = Output.Edges;

%% Proprocessing
% Three frames are required for tracking:
% - Dilated binarized mouse shape, for seed ROI
% - Lowpass filtered from to detect average pixel value
% - Lowpass filtered from to estimate pixel value

Frame = LoadFrame(Settings);

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

figure(1)
clf
imshow(Frame)
hold on



% Find tracing seeds

RoiInd = find(edge(SLL));

if isempty(RoiInd); return; end
clear RoiSub
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

if isempty(SeedProfile); return; end

[~, SeedInd] = findpeaks(SeedProfile, 'Threshold', 0.01);
Seeds = RoiSub(SeedInd,:);

dist =sqrt(sum( (Seeds - Output.Nose(Settings.Current_frame,:)).^2,2));
Cpoints = Seeds(dist<200,:);
Centre =( max(Cpoints,[],1)+min(Cpoints,[],1)) / 2;

% Track Traces
bufferSize = 100; % max nr of points on a trace


Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';
H = -15:15;

trace_finished = zeros(1, size(Seeds,1));


clc
for i = 2:bufferSize
    
    if ~any(any(Traces(:,:, i-1) ~= 1)) % If no points were assigned in previous loop
        break
    end
    
    Theta = [];
    if i == 2
        %ThetaRoi = repmat(1:20:360, [size(Seeds,1), 1]); % Roi bounded by angle range
        N = repmat(Centre,[size(Traces,2) 1]);
        vt = Traces(:,:,i-1) - N';
        ThetaRoi = atan2d(vt(2,:), vt(1,:))' + [-60:20:60];
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
    
    Iroi = Frame(RoiInd);    

    [~, PointSub] = min(Iroi, [], 2);
    PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
    
    [Traces(1, :, i), Traces(2,:,i)] = ind2sub(size(Frame), RoiInd(PointInd));
    
    trace_finished(CH(RoiInd(PointInd)) == 0) = trace_finished(CH(RoiInd(PointInd)) == 0) + 1;

    
    if i < 10
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
           
           rawax = 1:6;
           fitax = 1:0.5:5+10;
           
           px = polyfit(rawax, t(1,end-5:end), 1);
           py = polyfit(rawax, t(2,end-5:end), 1);
           
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
           
           if ptInd < length(fitax)-5
               LastPt = tfit(:, ptInd+1);
               ThetaLoop = squeeze(Theta(:, nan_traces(j), :));
               RoiSub = round(LastPt + ThetaLoop);
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
                   
               else
                   trace_finished( nan_traces(j) ) = trace_finished( nan_traces(j) )  + 1;
               end
           
           end
        end
        
    end
        
    
    Traces(1:2,trace_finished>2,i) = 1;
end
Traces(Traces == 1) = NaN;





scatter(Centre(2), Centre(1), 'b', 'filled')
scatter(Output.Nose(Settings.Current_frame,2), Output.Nose(Settings.Current_frame,1), 'g','filled')

%scatter(Seeds(:,2), Seeds(:,1), 'g','filled')
for i = 1:size(Traces,2)
    t = squeeze(Traces(:,i,:));
    idx = find(~isnan(t(1,:)),1,'last');
    t = t(:, 1:idx);
    if size(t,2) >10
    plot(t(2,:), t(1,:), 'b')
    end
   % text( t(2,end), t(1,end), num2str(i))
end









