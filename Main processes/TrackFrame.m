function TracesOut = TrackFrame(Settings, Output)
%% Traces = TrackFrame(Settings, Objects)
% Track raw traces in a frame using user defined traces
TracesOut = {};
Objects = Output.Objects;
Edges = Output.Edges;

%% Proprocessing
% Three frames are required for tracking:
% - Dilated binarized mouse shape, for seed ROI
% - Lowpass filtered from to detect average pixel value
% - Lowpass filtered from to estimate pixel value

Frame = LoadFrame(Settings);
Frame = im2double(Frame);
Frame = imadjust(Frame, [],[], Settings.Gamma);


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


SLN = zeros(size(Objects)); % Variable for silhouette (normal shape)
SLN( Frame <= Settings.Shape_threshold ) = 1;
SLN( Objects == 1 ) = 0;


SLS = imerode(SLN, strel('diamond', 5)); % Small silhouette
SLL = imdilate(SLS, strel('diamond', Settings.Dilation));

K = Settings.Trace_kernel_large;
L = Settings.Trace_kernel_small;
CHl = conv2(Frame, ones(K,K)./K^2, 'same'); % Low pass filter frame with large kernel
CHs = conv2(Frame, ones(L,L)./L^2, 'same');  % Low pass filter frame with small kernel
CH = zeros(size(Objects));
CH(0.5*CHs./CHl < Settings.Trace_threshold) = 1; % Matrix with potential valid pixels


FrameHeight = size(Frame, 1);
FrameWidth = size(Frame, 2);

Frame(1,1) = 999; % Default 'dead' pixel, to refer to if tracking should
% be interrupted.
CH(Objects == 1) = 0;
CH(Edges == 1) = 0.5;
%CH(SLL == 1) = 0;
CH(SLN == 1) = 0;


% Find tracing seeds




RoiInd = find(edge(SLL));

if isempty(RoiInd); TraceOut = {}; return; end %#ok<NASGU>
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

if isempty(SeedProfile); TracesOut = {}; return; end

if Settings.FullRoi == 0
    [~, SeedInd] =  findpeaks(SeedProfile, 'MinPeakDistance',3,'MinPeakProminence',Settings.Seed_threshold);
    Seeds = RoiSub(SeedInd,:);
elseif Settings.FullRoi == 1
    Seeds = RoiSub(find(CH(RoiInd) == 1),:);
end


%
if Settings.track_nose
    N = Output.Nose(Settings.Current_frame,:);
    d = sqrt( sum(  (Seeds-N).^2,2 ));
    Seeds = Seeds( d <=200 , :);
end





%%
%dist =sqrt(sum( (Seeds - Output.Nose(Settings.Current_frame,:)).^2,2));
%Cpoints = Seeds(dist<200,:);
%Centre = median(Cpoints,1);

% Track Traces
bufferSize = 100; % max nr of points on a trace


Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';
H = -15:15;
HH = -30:30;
trace_finished = zeros(1, size(Seeds,1));
trace_edges = zeros(1, size(Seeds,1));



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

%{
figure(1)
clf
imshow(Frame)
hold on

%{
[a,b] = ind2sub(size(Frame), X);
for i = 1:200
scatter(b(i,:),a(i,:),'r','filled')
end
%}

%{
[a,b] = ind2sub(size(Frame), RoiInd);
for i = 800
scatter(b(i,:), a(i,:), 'r','filled')
end
%}
%
scatter(Seeds(:,2), Seeds(:,1), 'g','filled')

scatter(Centre(2), Centre(1), 'g','filled')
for i = 1:size(Traces,2)
    t = squeeze(Traces(:,i,:));
    idx = find(~isnan(t(1,:)),1,'last');
    t = t(:, 1:idx);
    if size(t,2) > 10
        plot(t(2,:), t(1,:), 'b')
end

  
end
%}
