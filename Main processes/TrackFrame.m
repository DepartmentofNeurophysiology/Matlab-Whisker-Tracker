function [TracesOut, Misc] = TrackFrame(Settings, Data)
%% Traces = TrackFrame(Settings, Objects)
% Track raw traces in a frame using user defined traces
% Input:
% - Settings (struct with fields:)
%   - doGaussian (binary, apply gaussian filter)
%   - Gamme (scalar, 0-~, gamma rescaling
%   - Gaussian_kernel_size, (scalar, size of gaussian kernel)
%   - Shape_threshold (scalar, [0-1], threshold for shape extraction
%   - Dilation (integer, number of steps for dilation for seed detection
%   - Trace_kernel_large (integer, size of box kernel)
%   - Trace_kernel_small (integer, size of delta kernel)
%   - Trace_threshold (scalar, threshold of whisker detection)
%   - track_nose (binary, filter using nose distance)
%   - Current_frame (integer, current frame)
%
% - Data (struct with fields:)
%   - Objects (binary matrix marking environmental objects)
%   - Edges   (binary matrix marking edges)
%   - Nose    (n x 2) integer matrix, containing nose position per frame
%
%
% Output:
% - TracesOut (1 x n) cell, containing [2 x m] arrays for n traces of
%   length m
% - Misc (struct with optional fields for external use)


%% Preprocessing
TracesOut = {};
Objects = Data.Objects;
Edges = Data.Edges;

Frame = LoadFrame(Settings);
Frame = im2double(Frame);
Frame = imadjust(Frame, [],[], Settings.Gamma);


if Settings.doGaussian
    % Make gaussian kernel
    kfil = zeros(1, Settings.Gaussian_kernel_size);
    kfil(1:ceil(Settings.Gaussian_kernel_size/2)) = 1:ceil(Settings.Gaussian_kernel_size/2);
    kfil(ceil(Settings.Gaussian_kernel_size/2)+1:end) = max(kfil)-1:-1:1;
    Kgaus = repmat(kfil,[Settings.Gaussian_kernel_size, 1]) + repmat(kfil,[Settings.Gaussian_kernel_size, 1])';
    Kgaus = Kgaus./ sum(sum(Kgaus));
    
    % Convolute frame with gaussian
    Frame = conv2(Frame, Kgaus, 'same');
end


% Enhance contrast, normalize data
Frame = adapthisteq(Frame, 'NBins',256,'NumTiles',[5 5]);
minval = min(min(Frame));
Frame = Frame - minval;
Frame = Frame./max(max(Frame));


%% Detect whisker like edges


% Get mouse silhouettes
SLN = zeros(size(Objects)); % Silhouette Normal
SLN( Frame <= Settings.Shape_threshold ) = 1;
SLN( Objects == 1 ) = 0;

SLS = imerode(SLN, strel('diamond', 5)); % Silhouette Small
SLL = imdilate(SLS, strel('diamond', Settings.Dilation)); % Silhouette Large

% Filter edge like pixels
Box_kernel = conv2(Frame, ones(Settings.Trace_kernel_large,Settings.Trace_kernel_large)...
    ./Settings.Trace_kernel_large^2, 'same'); %  % Box kernel
Delta_kernel = conv2(Frame, ones(Settings.Trace_kernel_small,Settings.Trace_kernel_small)...
    ./Settings.Trace_kernel_small^2, 'same');  % Delta kernel

Control = zeros(size(Objects));
Control(0.5*Delta_kernel./Box_kernel < Settings.Trace_threshold) = 1; % Modified gaussian of laplacian


FrameHeight = size(Frame, 1);
FrameWidth = size(Frame, 2);

% Setup Frame and control matrix
Frame(1,1) = 999; % Default 'dead' pixel, to refer to if tracking should be interrupted.
Control(Objects == 1) = 0;
Control(Edges == 1) = 0.5;
Control(SLN == 1) = 0;





%% Detect Seeds


Roi = find(edge(SLL));
if isempty(Roi); TraceOut = {}; return; end %#ok<NASGU>
[Roi_sub(:,1), Roi_sub(:,2)] = ind2sub(size(Frame), Roi);
RS2 = Roi_sub;

% Add a pixel to the right to deal with rounding errors
RS2(:,1) = RS2(:,1)+1;
Roi_sub = [Roi_sub; RS2];
Roi_sub(Roi_sub(:,1) > FrameHeight,:) = 1;
Roi_sub(Roi_sub(:,2) > FrameWidth,:) = 1;
Misc.SeedRoi = Roi_sub;
%
Roi = sub2ind(size(Frame), Roi_sub(:,1), Roi_sub(:,2));
Seeds = Roi_sub(find(Control(Roi) == 1),:); %#ok<FNDSB>


if Settings.track_nose
    N = Data.Nose(Settings.Current_frame,:);
    d = sqrt( sum(  (Seeds-N).^2,2 ));
    Seeds = Seeds( d <=200 , :);
end




%% Whisker tracking

% Track Traces
bufferSize = 100; % max nr of points on a trace
Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';


Finished = zeros(1, size(Seeds,1)); % Keep track of finished seeds
trace_on_edges = zeros(1, size(Seeds,1)); % Keep track of traces on edges




for i = 2:bufferSize
    
    
    
    if ~any(any(Traces(:,:, i-1) ~= 1)) % If no points were assigned in previous loop
        break
    end
    
    %
    %
    % FIND ROI
    %
    %
    % Generate arcs per seed/last tracked point
    Theta = [];
    if i == 2
        ThetaRoi = repmat(1:20:360, [size(Seeds,1), 1]); % Roi bounded by angle range
    else
        vt = Traces(:,:,i-1)- Traces(:,:,i-2);
        ThetaRoi = atan2d(vt(2,:), vt(1,:))' + [-15:15]; %#ok<*NBRAK>
    end
    Theta(1,:,:) = [3*cosd(ThetaRoi)]; % ROI per previous point
    Theta(2,:,:) = [3*sind(ThetaRoi)];
    
    % Define new ROI
    LastPoints = repmat(Traces(:,:,i-1), [1 , 1, size(Theta, 3)]);
    Roi_sub = round(LastPoints + Theta);
    
    % Remove indices outside of frame range
    Roi_sub(Roi_sub < 1) = 1;
    Roi_sub(1,Roi_sub(1,:,:) > FrameHeight) = 1;
    Roi_sub(2,Roi_sub(2,:,:) > FrameWidth) = 1;
    Roi = squeeze(sub2ind(size(Frame), Roi_sub(1,:,:), Roi_sub(2,:,:)));    
    if i == 2
        Roi = cut_roi(Roi, SLL);
    end
    
    %
    %
    % FIND LOCAL MINIMA
    %
    %
    Iroi = Frame(Roi);
    Iroi(Control(Roi) == 0) = NaN;
    if size(Seeds,1) == 1
        [~, minima] = min(Iroi);
        Local_min = minima;
    else
        [~, minima] = min(Iroi, [], 2);
        Local_min = sub2ind(size(Roi), 1:length(minima), minima');
    end
    
    
    %
    %
    % ASSIGN NEW POINTS
    %
    %
    [Traces(1, :, i), Traces(2,:,i)] = ind2sub(size(Frame), Roi(Local_min));
    
    % If a new found point does not pass control, +1
    Finished(Control(Roi(Local_min)) == 0) = Finished(Control(Roi(Local_min)) == 0) + 1;
    trace_on_edges(Control(Roi(Local_min)) == 0.5) = trace_on_edges(Control(Roi(Local_min)) == 0.5) + 1;
    
   
    %
    %
    % LINEAR EXTRAPOLATION
    %
    %
    if i < 3
        Finished(Control(Roi(Local_min)) == 0.5) =  Finished(Control(Roi(Local_min)) == 0.5) + 1;
    elseif any(Control(Roi(Local_min)) == 0.5)
        
        
        
        
        traces_on_edge = find(Control(Roi(Local_min)) == 0.5);
        
        
        for j = 1:length(traces_on_edge)
            
            
            % Get trace
            trace = squeeze(Traces(:, traces_on_edge(j), :));
            i1 = find(trace(1,:) == 1, 1, 'first');
            i2 = find(trace(2,:) == 1, 1, 'first');
            if ~isempty([i1 i2])
                trace = trace(:,1:max([i1 i2])-1);
            end
            
            % Fit linear polynomial
            if i < 10
                rawax = 1:i;
                fitax = 1:0.5:i+10;
            else
                rawax = 1:10;
                fitax = 1:0.5:i+10;
            end
            
            px = polyfit(rawax, trace(1,end-length(rawax)+1:end), 1);
            py = polyfit(rawax, trace(2,end-length(rawax)+1:end), 1);
            
            tfit = [];
            tfit(1,:) = round(polyval(px, fitax));
            tfit(2,:) = round(polyval(py, fitax));
            
            % Keep values within frame range
            throwidx = [find(tfit(1,:) < 1), find(tfit(1,:) > FrameHeight)];
            throwidx = [throwidx, find(tfit(2,:) < 1), find(tfit(2,:) > FrameWidth)];            
            if ~isempty(throwidx)
                tfit(:, throwidx) = NaN;
                tfit = tfit(:, ~isnan(tfit(1,:)));
            end
            
            % Detect last point of fit on edge
            tfitInd = sub2ind(size(Frame), tfit(1,:), tfit(2,:));
            ptInd = find(Control(tfitInd) == 0.5,1,'last');
            
            if ptInd < size(tfit,2)-5
                
                
                % Use last point as seed
                LastPt = tfit(:, ptInd+1);
                vt = tfit(:, end) - tfit(:, end-5);
                
                % Create arc for ROI
                T1 = atan2d(vt(2,1), vt(1,1)) + [ -30:30;];
                T2 = [];
                T2(1,:) = 6*cosd(T1)';
                T2(2,:) = 6*sind(T1)';                
                Roi_sub = round(LastPt + T2);
                Roi_sub(Roi_sub < 1) = 1;
                Roi_sub(1,Roi_sub(1,:) > FrameHeight) = 1;
                Roi_sub(2,Roi_sub(2,:) > FrameWidth) = 1;
                Roi = squeeze(sub2ind(size(Frame), Roi_sub(1,:), Roi_sub(2,:)));
                
                % Find local minima
                Iroi = Frame(Roi);                
                [~, minima] = min(Iroi, [], 2);
                Local_min = sub2ind(size(Roi), 1:length(minima), minima');
                [a,b] = ind2sub(size(Frame), Roi(Local_min));
                
                % Validate new point
                if Control(a,b) == 1
                    Traces(1, traces_on_edge(j), i) = a;
                    Traces(2, traces_on_edge(j), i) = b;
                    trace_on_edges(traces_on_edge(j)) = 0;                    
                else
                    Finished( traces_on_edge(j) ) = Finished( traces_on_edge(j) )  + 1;
                end
                
            end
        end
        
    end
    
    % Stop tracking traces that are finished
    Traces(1:2,Finished>2,i) = 1;
    
    % Delete points on objects/edges
    idx = find(trace_on_edges > 10 & Traces(1,:,i)>1);
    if ~isempty(idx)
        for j = 1:length(idx)
            n_delete = trace_on_edges(idx(j));
            Traces(:, idx(j), i-n_delete:i) = 1;
        end
    end
    
    
end


Traces(Traces == 1) = NaN;


% Reshape results
TOut = {};
TracesOut = {};
Or = [];
for i = 1:size(Traces,2)
    trace = squeeze(Traces(:,i,:));
    idx = find(~isnan(trace(1,:)),1,'last');
    trace = trace(:, 1:idx);
    if size(trace,2) > 10
        TOut{end+1} = trace';
        Or(end+1,:) = trace(:,1)';
    end
end


% Remove double tracked traces
Misc.CH = Control;
Control(Edges == 1) = 0;
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
            idx = sub2ind(size(Frame), TOut{i}(:,1), TOut{i}(:,2));
            idx(isnan(idx)) = 1;
            
            if sum(Control(idx))/length(idx) > 0.8
                TracesOut{end+1} = TOut{i};
            end
        end
    end
end



Misc.Seeds = Seeds;
Misc.RoiSub = Roi_sub;
Misc.SLL = SLL;
Misc.SLN = SLN;

end



function Roi = cut_roi(Roi,SLL)
[a,b] = find(SLL(Roi));
for j = 1:size(Roi, 1)
    idx = b(a == j);
    if ~any(diff(idx)>1) & min(idx) > 2 & max(idx) < size(Roi,2)-1 %#ok<*AND2>
        idx = [[-2 -1]'+idx(1); idx; [1 2]'+idx(end)]; %#ok<*AGROW>
    elseif ~any(diff(idx) > 1) & min(idx) == 2 & max(idx) < size(Roi, 2) - 1
        idx = [-1+idx(1); idx; [1 2]'+idx(end); size(Roi,2)];
    elseif ~any(diff(idx) > 1) & min(idx) == 1 & max(idx) < size(Roi, 2) - 1
        idx = [idx; [1 2]'+idx(end);[-1 0]'+size(Roi,2)];        
    elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(Roi,2)
        idx = [ [1 2]'; [-2 -1]'+idx(1); idx];
    elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(Roi,2)-1
        idx = [ [1]'; [-2 -1]'+idx(1); idx; size(Roi,2)];        
    elseif numel(find(diff(idx) > 1)) == 1 & min(idx) == 1 & max(idx) == size(Roi,2)
        id = find(diff(idx) > 1);
        idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];        
    elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) < size(Roi,2) - 1
        idx = [ [-2 -1]'+min(idx); idx ; [1 2]'+idx(end)];        
    elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) == size(Roi,2) - 1
        idx = [1; [-2 -1]'+min(idx); idx ; [1 ]'+idx(end)];     
    elseif numel(find(diff(idx) > 1)) > 1 & min(idx) == 1 & max(idx) == size(Roi,2)
        [~, id] = max(diff(idx));
        idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
    else
        idx =idx; %#ok<ASGSL>
    end
    Roi(j,idx) = 1;
end
end
