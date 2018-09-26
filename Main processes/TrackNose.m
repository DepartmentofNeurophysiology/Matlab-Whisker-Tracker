function Output = TrackNose(Settings,Output)
%%
% Tracks nose position and headangle:
%   - Estimate movement axis (X/Y) and direction (Up/Down)
%   - Track nose and headangle
%%

% Extract variables
nframes = Settings.Nframes;
stepsize = Settings.nose_interval;
Base(1:numel(1:stepsize:nframes),1:2) = NaN;
Frames_filtered = cell(numel(1:stepsize:nframes),1);


% Extract centre of mass (Base) for frames
store_idx = 1;
h = waitbar(0,'Tracking Nose');
for framenr = 1:stepsize:nframes
    
    % Load frame
    Settings.Current_frame = framenr;
    frame_in = LoadFrame(Settings);
    frame = frame_in;
    
    % binarize frame
    frame( frame > Settings.Silhouettethreshold ) = 0;
    frame( find( frame )) = 1;  %#ok<*FNDSB>
    frame( find( Output.Objects )) = 0;
    
    % Apply mask to filter objects/cable
    frame_normal = frame; 
    frame = imdilate( frame, strel('diamond',5));
    frame = imerode(frame, strel('diamond',45));
    mask = imdilate(frame, strel('diamond',40));    
    frame = frame_normal.*mask; 
    
    
    Frames_filtered{store_idx} = frame;  
    
    % Only continue if mouse is present in frame
    if numel(find(frame)) < 3000       
        Frames_filtered{store_idx}(:) = 0;        
        waitbar(framenr/nframes)
         store_idx = store_idx+1;
        continue
    end
    
    % Store CoM
    sx = sum(frame,1);
    sy = sum(frame,2);
    Base(store_idx,:) = [mean(find(sy)) mean(find(sx))];
    store_idx = store_idx+1;
    waitbar(framenr/nframes)
end
close(h);
ax_default = 'Y'; % in our data mice only move along Y, verify correct movement

% Logic to extract movement axis
dB = diff(Base,1);
dB(abs(dB) > 20) = NaN;
dB = abs(sum(dB,1,'omitnan'));
if dB(1) > dB(2)
    ax_dir = 'Y';
    base_idx = 1;
    f_idx = 2;
  
elseif dB(2) >= dB(1)
    ax_dir = 'X';
    base_idx = 2;
    f_idx=1;
   
end

if ~strcmp(ax_default,  ax_dir)
    disp('Measured direction axis does not match default axis...')
    disp('continue with default axis')
end

% Logic to extract movement direction
xax = 1:stepsize:nframes;
idx = find(~isnan(Base(:, base_idx)));
base_slope = polyfit(xax(idx), Base(idx, base_idx)',1);
if base_slope(1) < 0
    Direction = 'Up';
elseif base_slope(1) > 0
    Direction = 'Down';
end



%% Nose and headangle tracking

Nose_raw = [];
theta = 1:1:360;

Nose_raw(1:size(Frames_filtered,1), 1:2) = NaN;
AngleVector(1:size(Frames_filtered,1) ,1:2) = NaN;

for i = 1:size(Frames_filtered,1)
    if isnan(Base(i,1))  
        continue
    end
    
    % Find nose as point closest to target platform
    frame = Frames_filtered{i};
    f_sum = sum(frame ,f_idx)./ size(frame, f_idx);
    switch(Direction)
        case 'Up'
            X = find( f_sum > 0, 1, 'first');
            
        case 'Down'
            X = find( f_sum > 0, 1, 'last');
    end
    Y = round( mean( find( frame(X,:)), 'omitnan' ));
    Nose_raw(i,:) = [X, Y];    
    frame = edge(frame);
    
    
    % Get headangle:
    % - find non-zero entries in circle around nose
    Cx = round(Nose_raw( i, 2) + 40*sind(theta));
    Cx = [Cx, Cx+1]; %#ok<AGROW>
    Cy = round(Nose_raw( i, 1) + 40*cosd(theta));
    Cy = [Cy, Cy+1]; %#ok<AGROW>
    IDX = find(Cx > 1 & Cx < size(frame,2) & Cy > 1 & Cy< size(frame,1));
   
    C = sub2ind( size(frame), Cy(IDX), Cx(IDX) ); 
    PTS = [];
    PTS(:,1) = find( frame( C ));
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
    
    
    
    
    if ~isempty(PTS) && size(PTS,1) == 2
        PTS = PTS(1:2,:); % Lazy
        Vp = PTS(2,:) - PTS(1,:);
        
        % Headangle described as vector, normal to Vp
        AngleVector(i,1:2) = [-Vp(1),Vp(2)];
        
        % Normalise
        AngleVector(i,1:2) = AngleVector(i,1:2)...
            ./sqrt(sum(AngleVector(i,1:2).^2));
        
        
        % Determine if facing right direction (towards tail)
        switch(Direction)
            case 'Up'
                if AngleVector(i,1) < 0
                    AngleVector(i,:) = -AngleVector(i,:);
                end
            case 'Down'
                if AngleVector(i,1) > 0
                    AngleVector(i,:) = -AngleVector(i,:);
                end
            case 'Left'
                if AngleVector(i,2) > 0
                    AngleVector(i,:) = -AngleVector(i,:);
                end
            case 'Right'
                if AngleVector(i,2) < 0
                    AngleVector(i,:) = -AngleVector(i,:);
                end
        end
        
        
    else
        AngleVector(i,1:2) = NaN;
    end    
   
end


%% Fit data
track_idx = 1:stepsize:nframes;
keep = find(~isnan( Nose_raw(:,1)));
track_idx = track_idx(keep);
fitax = track_idx(1):track_idx(end);
Nose(1:nframes,1:2) = NaN;
Nose(fitax,1) = spline(track_idx , Nose_raw( keep, 1), fitax);
Nose(fitax,2) = spline(track_idx , Nose_raw( keep, 2), fitax);


track_idx = 1:stepsize:nframes;
keep = find(~isnan(AngleVector(:,1)));
track_idx = track_idx(keep);
fitax = track_idx(1): track_idx(end);
AngleVec(1:nframes,1:2) = NaN;
AngleVec(fitax,1) = spline(track_idx, AngleVector( keep, 1), fitax);
AngleVec(fitax,2) = spline(track_idx, AngleVector( keep, 2), fitax);


Output.Direction = Direction;
Output.Nose = Nose;
Output.AngleVector = AngleVec;






