function [Touch, TouchFiltered] = DetectTouch(Traces, Edges, edgeIDX)
% Detect Traces endings within the target platform edge in
% 'Edges', ridge specified as 'edgeIDX'

%%
edge_width = 25; % Roi width of detection
max_dist = 5:5:50; % distance tip-edge (max) to detect touch 
ndist = length(max_dist);

% Detect touch only around target platform

mask = zeros(size(Edges));
nframes = max(size(Traces));

Touch = cell(ndist, nframes);
TouchPoints = cell(ndist,nframes);


if ~isnan(edgeIDX(1))
    for i = 1:length(edgeIDX)

        if edgeIDX(i) <= edge_width | edgeIDX(i)+edge_width > size(mask,1)
            continue
        end
    mask(edgeIDX(i)-edge_width:edgeIDX(i)+edge_width,:) = 1;
    Edges = Edges.*mask;
    end
end

[Opts(:,1), Opts(:,2)] = find(Edges);


h = waitbar(0, 'Detecting Touch');


for i = 1:nframes   
    
    for distidx = 1:ndist
        looptouch = zeros(1, size(Traces{i}, 2));
        looppoints = [];
        for j = 1:size(Traces{i}, 2)
            t = Traces{i}{j};
            if isempty(t)
                continue
            end

            dist = Opts - t(end,:);       
            dist = sqrt( sum( dist.^2, 2));
            [~, tpt] = find(dist <= max_dist(distidx));

            if ~isempty(tpt)
                looptouch(j) = 1;
                
                looppoints(end+1,:) = Traces{i}{j}(end,:);
            end
        end


        Touch{distidx,i} = looptouch;
        TouchPoints{distidx,i} = looppoints;
    end
    
    waitbar(i/nframes)
end

close(h);
  



%%
MatchNext = cell(size(Touch));
MatchPrev = cell(size(Touch));
TouchFiltered = cell(size(Touch));
for i = 2:size(MatchNext, 2)
    
    looptouch = zeros(size(Touch{i}));
    
    
    if ~isempty(TouchPoints{i})
        MatchNext{i} = zeros(1,size(TouchPoints{i},1));
        MatchPrev{i} = zeros(1,size(TouchPoints{i},1));
        
        if ~isempty(TouchPoints{i+1})
            for j = 1:size(TouchPoints{i}, 1)
                dx = TouchPoints{i}(j,1) - TouchPoints{i+1}(:,1);
                dy = TouchPoints{i}(j,2) - TouchPoints{i+1}(:,2);
                d = sqrt( dx.^2 + dy.^2);
                [val, id] = min(d);
                if val < 20
                    MatchNext{i}(j) = id;
                end
            end
        end
        
        if ~isempty(TouchPoints{i-1})
            for j = 1:size(TouchPoints{i},1)
                dx = TouchPoints{i}(j,1) - TouchPoints{i-1}(:,1);
                dy = TouchPoints{i}(j,2) - TouchPoints{i-1}(:,2);
                d = sqrt( dx.^2 + dy.^2);
                [val, id] = min(d);
                if val < 20
                    MatchPrev{i}(j) = id;
                end
            end
        end
        
        touchidx = find(Touch{i});
        idx = find(MatchNext{i} | MatchPrev{i});
        for j  = 1:length(idx)
            looptouch(touchidx(idx(j))) = 1;
        end
        
        TouchFiltered{i} = looptouch;
        
        
    end
end


