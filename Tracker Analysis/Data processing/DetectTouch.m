function Touch = detectTouch(Traces, Edges, edgeIDX)
% Detect Traces endings within the target platform edge in
% 'Edges', ridge specified as 'edgeIDX'

%%
edge_width = 5; % Roi width of detection
max_dist = 10; % distance tip-edge (max) to detect touch 


% Detect touch only around target platform

mask = zeros(size(Edges));
nframes = max(size(Traces));
mask(edgeIDX-edge_width:edgeIDX+edge_width,:) = 1;
Edges = Edges.*mask;

[Opts(:,1), Opts(:,2)] = find(Edges);

Touch = cell(1, nframes);
h = waitbar(0, 'Detecting Touch');

for i = 1:nframes    
    looptouch = zeros(1, size(Traces{i}, 2));
    for j = 1:size(Traces{i}, 2)
        t = Traces{i}{j};
        if isempty(t)
            continue
        end
        
        dist = Opts - t(end,:);       
        dist = sqrt( sum( dist.^2, 2));
        [~, tpt] = find(dist <= max_dist);
        
        if ~isempty(tpt)
            looptouch(j) = 1;
        end
    end
    
    Touch{i} = looptouch;
end

close(h);
    



