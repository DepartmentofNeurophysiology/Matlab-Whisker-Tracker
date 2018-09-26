function Output = TrackFrame(Settings, Output)
%Data = TrackFrame(Data,Settings) returns a field .Traces in Data
% containing all tracked traces of the current frame using the intensity
% convergence algorithm:
% -	Expand mouse silhouette
% -	Detect origins using peakdetection
% - Repeatedly apply 'FindPoint' untill all origins have been traced
%% Initialize

Tracked = zeros(size(Output.Objects));
Traces = [];
Origins = [];
Objects = Output.Objects;
frame = LoadFrame(Settings);



% Find Trace origins
Silhouette = zeros(size(Objects));
Silhouette( find(frame <= Settings.Silhouettethreshold) ) = 1; %#ok<*FNDSB>
Silhouette(find(Objects)) = 0;
Shapes.normal = Silhouette;
Shapes.normal = imerode(Silhouette,strel('diamond',2));
%Data.Silhouette( find(Data.Objects) ) = 0;
SilhouetteSmall = imerode(Silhouette,strel('diamond',5));
SE = strel('disk',Settings.Dilationsize);
Silhouette = imdilate(SilhouetteSmall,SE);
TraceOriginIDX = find(edge(Silhouette));
TraceOriginIDX = unique(TraceOriginIDX);

if isempty(TraceOriginIDX)
    return
end

% Arange originIDX acoording to shape
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
%* due to artefactis in background extraction, the origin profile looks
%funky sometimes. By subtraciting a lowpass filter and taking the absolute
%value this is corrected for.

if ~isempty(OriginProfile)
    [~,OriginsIDX] = findpeaks(OriginProfile,'Threshold',...
        Settings.Origin_threshold);
else
    OriginsIDX = [];
end

Shapes.small = SilhouetteSmall;
Shapes.large = Silhouette;
Shapes.Objects = Objects;
Shapes.Frame = frame;
Shapes.Tracked = Tracked;




%%
if ~isempty(OriginsIDX)
    
    
    [Origins(:,1),Origins(:,2)] = ind2sub(size(frame),...
        TraceOriginIDX(OriginsIDX));
    
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
    
    
    
    %Data.Origins = flip(Data.Origins,1);
    
    
    
    traceidx = 1;

    Settings.object_tick = 0;   
    for idx = 1:size(Origins,1)
      
        [Trace, Shapes] = TrackTrace(Settings, Shapes, Origins(idx,1:2));
        
        
        
        if ValidateTrace(Settings,Trace)
            Traces{traceidx} = Trace;
            traceidx = traceidx+1;
            
            
           
            
        else % cleanup tracking record
            for i = 1:size(Trace,1)
                Shapes.Tracked(Trace(i,1)-Settings.stepsize+2:Trace(i,1)+Settings.stepsize-2,...
                    Trace(i,2)-Settings.stepsize+2:Trace(i,2) + Settings.stepsize-2) = 1;
            end
            
        end
        
    end
    
    
end

Output.Origins = Origins;
Output.Traces = Traces;

%{

figure(1)
clf
imagesc(frame)
colormap('gray')
hold on



for i = 1:size(Traces,2)
    plot(Traces{i}(:,2), Traces{i}(:,1),'r')
end
scatter(Origins(:,2), Origins(:,1), 'b','filled')
for i = 1:size(Origins,1)
    text(Origins(i,2), Origins(i,1)+5,num2str(i))
end
%}

