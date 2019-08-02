function Output = TrackNose(Settings,Output)
%%
% Tracks nose position and headangle:
%   - Estimate movement axis (X/Y) and direction (Up/Down)
%   - Track nose and headangle


%%


if isfield(Settings,'DefaultDirection') && ~isempty(Settings.DefaultDirection)
    Default_direction = Settings.DefaultDirection;
else
    Default_direction = [];
    if exist(fullfile(Settings.PathName, 'Selected_frames.mat'), 'file')
        dataIn = load(fullfile(Settings.PathName, 'Selected_frames.mat'));
        id = [];
        for i = 1:size(dataIn.Output, 2)
            if strcmp(Settings.FileName, dataIn.Output(i).Video)
                id = i;
            end
        end

        if ~isempty(id)
            if isfield(dataIn.Output(id), 'Direction')
                if ~isempty(dataIn.Output(id).Direction)
                    if dataIn.Output(id).Direction == 1
                        Default_direction = 'Up';
                    elseif dataIn.Output(id).Direction == 2
                        Default_direction = 'Down';
                    end
                end
            end

        end
    end
end



% Extract variables
nframes = Settings.Nframes;
stepsize = 5;%Settings.nose_interval;
frameidx = 1:stepsize:nframes;
Base(1:length(frameidx),1:2) = NaN;

% Extract centre of mass (Base) for frames
O = Output.Objects;
O = imdilate(O,strel('diamond',3));


h = waitbar(0,'Tracking Nose');
for i = 1:length(frameidx)
    
    % Load Frame
    Settings.Current_frame = frameidx(i);
    Frame = LoadFrame(Settings);
    Frame = im2double( abs(Frame));
    
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
    
    SLN = zeros(size(Frame));
    SLN(Frame <= Settings.Shape_threshold) = 1;    
    SLN(O==1) = 0;
    
%     FM = imdilate(SLN, strel('diamond', 5));
%     FM = imerode(FM, strel('diamond', 45));
%     FM = imdilate(Frame, strel('diamond', 40));
    
    Frame = SLN;

    
    % Only continue if mouse is present in Frame
    if numel(find(Frame)) < 5  
        waitbar(i/(2*length(frameidx)))
        continue
    end
    
    % Store CoM
    sx = sum(Frame,1);
    sy = sum(Frame,2);
    Base(i,:) = [mean(find(sy)) mean(find(sx))];   
    waitbar(i/(2*length(frameidx)))
end


%%
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
    ax_dir = 'Y';
    base_idx = 1;
    f_idx = 2;
end

% Logic to extract movement direction
xax = 1:stepsize:nframes;
idx = find(~isnan(Base(:, base_idx)));
base_slope = polyfit(xax(idx), Base(idx, base_idx)',1);
if base_slope(1) < 0
    Direction = 'Up';
elseif base_slope(1) > 0
    Direction = 'Down';
elseif base_slope(1) == 0
    Direction = 'Up';
end

if ~isempty(Default_direction)
    fprintf('Old direction: %s , New direction: %s\n', Direction, Default_direction);
    Direction = Default_direction;
end



%% Nose and headangle tracking

Nose_raw = [];
theta = 1:1:360;

Nose_raw(1:length(frameidx), 1:2) = NaN;
AngleVector(1:length(frameidx) ,1:2) = NaN;

for i = 1:length(frameidx)
    if isnan(Base(i,1))  
        continue
    end
    
    % Load Frame
    Settings.Current_frame = frameidx(i);
    Frame = LoadFrame(Settings);
    Frame = im2double( abs(Frame));
    
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
    
    SLN = zeros(size(Frame));
    SLN(Frame <= Settings.Shape_threshold) = 1;
    SLN(O == 1) = 0;
    
    FM = imdilate(SLN, strel('diamond', 1));
    %FM = imerode(FM, strel('diamond', 45));
    %FM = imdilate(FM, strel('diamond', 55));
    
    Frame = SLN.*FM;  
    
    f_sum = sum(Frame ,f_idx)./ size(Frame, f_idx);
    switch(Direction)
        case 'Up'
            X = find( f_sum > 0, 1, 'first');
            
        case 'Down'
            X = find( f_sum > 0, 1, 'last');
    end
    
    
    Y = round( mean( find( Frame(X,:)), 'omitnan' ));
    Nose_raw(i,:) = [X, Y];    
    Frame = edge(Frame);
        
    % Get headangle:
    % - find non-zero entries in circle around nose
    Cx = round(Nose_raw( i, 2) + 40*sind(theta));
    Cx = [Cx, Cx+1]; %#ok<AGROW>
    Cy = round(Nose_raw( i, 1) + 40*cosd(theta));
    Cy = [Cy, Cy+1]; %#ok<AGROW>
    IDX = find(Cx > 1 & Cx < size(Frame,2) & Cy > 1 & Cy< size(Frame,1));
   
    C = sub2ind( size(Frame), Cy(IDX), Cx(IDX) ); 
    PTS = [];
    PTS(:,1) = find( Frame( C ));
    [PTS(1:length(PTS),2),PTS(1:length(PTS),1)] = ind2sub(size(Frame),C(PTS));    
  
   
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
    
        waitbar((length(frameidx)+i)/(2*length(frameidx)))

    
   
end

close(h);
%% Fit data
n = 20;
track_idx = 1:stepsize:nframes;
keep = find(~isnan( Nose_raw(:,1)));
if isempty(keep)
    Nose(1:nframes,1:2) = NaN;
    AngleVec(1:nframes,1:2) = NaN;
    Output.Direction = Direction;
    Output.Nose = Nose;
    Output.AngleVector = AngleVec;
    return
end
       
    
track_idx = track_idx(keep);
fitax = track_idx(1):track_idx(end);
Nose(1:nframes,1:2) = NaN;
sfitX = spline(track_idx , Nose_raw( keep, 1), fitax);
sfitY = spline(track_idx , Nose_raw( keep, 2), fitax);
Nose(fitax,1) = filter( ones(1,n)/n, 1, sfitX);
Nose(fitax,2) = filter( ones(1,n)/n, 1, sfitY);

AngleVec(1:nframes,1:2) = NaN;
Araw = AngleVector;

meanX = mean(Araw(:,1),'omitnan');
stdX = std(Araw(:,1), 'omitnan');

meanY = mean(Araw(:,2), 'omitnan');
stdY = std(Araw(:,2), 'omitnan');

Arawf = Araw;
Arawf( abs(Araw(:,1)-meanX) > 3*stdX   ,:) = NaN;
Arawf( abs(Araw(:,2)-meanY) > 3*stdY   ,:) = NaN;

track_idx = 1:stepsize:nframes;
dTheta = diff(Arawf,1);
dT =diff(track_idx)';
slope(:,1) = dTheta(:,1)./dT;
slope(:,2) = dTheta(:,2)./dT;
slope(end+1,:) = 0;

keep = find(~isnan(Nose_raw(:,1)) & ~isnan(Arawf(:,1)) & abs(slope(:,1)) < 1 & abs(slope(:,2))<1);
track_idx = 1:stepsize:nframes;
track_idx = track_idx(keep);
if numel(find(~isnan(AngleVector))) > 4
    
    sfitX = spline(track_idx, Arawf(keep,1), fitax);
    sfitY = spline(track_idx, Arawf(keep,2), fitax);

    AngleVec(fitax,1) = filter( ones(1,n)/n, 1, sfitX );
    AngleVec(fitax,2) = filter( ones(1,n)/n, 1, sfitY );
end

%%



Output.Direction = Direction;
Output.Nose = Nose;
Output.AngleVector = AngleVec;









%%

