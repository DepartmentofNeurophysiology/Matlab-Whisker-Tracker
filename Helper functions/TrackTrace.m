function [Trace, Shapes] = TrackTrace(Settings, Shapes, Origin)
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
%[Cx, Cy] = ind2sub(size(frame), Cidx);

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





%%

% Repeat untill snout is found
while ptflag
    vt = Trace(end,:) - Trace(end-1,:);
    Angle = atan2d(vt(1), vt(2));
    HalfCircle = -30:30;
    Theta = HalfCircle + Angle;
    
    Cx = ceil(Trace(end,1) + first_stepsize*sind(Theta));
    Cy = ceil(Trace(end,2) + first_stepsize*cosd(Theta));
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
    
   
end


%{
figure(1)
clf
imagesc(frame)
colormap('gray')
hold on


scatter(Origin(2), Origin(1), 'b','filled')
scatter(Cy, Cx, 'y','filled')
for i = 1:size(Trace,1)
    scatter(Trace(i,2), Trace(i,1), 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
    text(Trace(i,2), Trace(i,1), num2str(i))
end
scatter(newpt(2), newpt(1), 'g','filled')
%}


%
%
% Track towards the tip
Trace = flip(Trace,  1);

% Start with a linear fit on previous found point to estimate direction

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


%%


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
%%
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
    
    
end

if noiseflag > 0
    remove_idx = noiseflag;
    Trace = Trace(1:end-remove_idx,:);
end

%}
%{

figure(1)
clf
imagesc(frame)
colormap('gray')
hold on


scatter(Origin(2), Origin(1), 'b','filled')
scatter(Cy, Cx, 'y','filled')
for i = 1:size(Trace,1)
    scatter(Trace(i,2), Trace(i,1), 'MarkerFaceColor','r', 'MarkerEdgeColor','r')
    
end
scatter(newpt(2), newpt(1), 'g','filled')


%}



