function Labels = getLabels(Janelia)

%% Assign side labels
nframes = Janelia.MetaData.NFrames;
Labels.side = cell( nframes, 1);

for i = 1:nframes
    
    % use corrected whisker angle to assign side 
    ntraces = size( Janelia.Traces_clean{i}, 2);
    if ntraces == 0
        continue
    end
    
    Labels.side{i} = zeros(1, ntraces);
    Labels.side{i}( Janelia.Parameters{i}(:, 10) <= 0 ) = 1; 
    Labels.side{i}( Janelia.Parameters{i}(:, 10) > 0 ) = 2;
    
    if any(Labels.side{i} == 0)
        idx = find(Labels.side{i} == 0);
        if any(~isnan(Janelia.Parameters{i}(idx,6)))
            keyboard
        end
    end
    
end

%% Assign 'column' label

% assign # whiskers on rostral and distal side
n_rost = 3;
n_dist = 5;

Labels.dist = cell( nframes, 1);

for i = 1:nframes
   
    ntraces = size( Janelia.Traces_clean{i}, 2);
    if ntraces == 0
        continue
    end
    
    Labels.dist{i} = zeros(1, ntraces);
    
    
    
    % assign left side
    left_idx = find( Labels.side{i} == 1);
    
    px = Janelia.Parameters{i}(left_idx, 9); % whisekr tip x
    py = Janelia.Parameters{i}(left_idx, 10); % whisker tip y
    dist_nose = sqrt( sum( [px,py].^2,2));
    [~, sort_idx] = sort(dist_nose);    
    
    %{
    if length(left_idx) > n_dist+n_rost
        Labels.dist{i}(left_idx( sort_idx(end-n_dist+1:end))) = 3;
        Labels.dist{i}(left_idx( sort_idx(1:n_rost))) = 1;
        Labels.dist{i}(left_idx( sort_idx(n_rost+1:end-n_dist))) = 2;
    elseif length(left_idx) > n_dist
        Labels.dist{i}(left_idx( sort_idx(end-n_dist+1:end))) = 3;
        Labels.dist{i}(left_idx( sort_idx(1:end-n_dist))) = 2;
    else
        Labels.dist{i}(left_idx) = 3;
    end
    %}
    n_frame = length(left_idx);
    n_rost = floor(n_frame/2);   
    Labels.dist{i}(left_idx( sort_idx(1:n_rost) ) ) = 2;
    Labels.dist{i}(left_idx( sort_idx(n_rost+1:end) )) = 3;
    
    % assign right side
    right_idx = find( Labels.side{i} == 2);
    
   px = Janelia.Parameters{i}(right_idx, 9); % whisekr tip x
    py = Janelia.Parameters{i}(right_idx, 10); % whisker tip y
    dist_nose = sqrt( sum( [px,py].^2,2));
    [~, sort_idx] = sort(dist_nose);  
    
    
    %{
    if length(right_idx) > n_dist+n_rost
        Labels.dist{i}(right_idx( sort_idx(end-n_dist+1:end))) = 3;
        Labels.dist{i}(right_idx( sort_idx(1:n_rost))) = 1;
        Labels.dist{i}(right_idx( sort_idx(n_rost+1:end-n_dist))) = 2;
    elseif length(right_idx) > n_dist
        Labels.dist{i}(right_idx( sort_idx(end-n_dist+1:end))) = 3;
        Labels.dist{i}(right_idx( sort_idx(1:end-n_dist))) = 2;
    else
        Labels.dist{i}(right_idx) = 3;
    end
    %}
    
    n_frame = length(right_idx);
    n_rost = floor(n_frame/2);
    Labels.dist{i}(right_idx( sort_idx(1:n_rost))) = 2;
    Labels.dist{i}(right_idx( sort_idx(1:n_rost))) = 3;
    
    
    
    
end

