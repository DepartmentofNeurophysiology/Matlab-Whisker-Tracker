function [Traces_clean] = CleanTraces(Tracker,measuredNose)
%%
if nargin == 1
    measuredNose = 1;
end
Settings = makeAnalyseSettings;

nframes = size( Tracker.Traces,1);
tracked_frames = ones(1,nframes);

Nose = Tracker.Nose;
frame_heigth = size(Tracker.Objects, 1);
frame_width = size(Tracker.Objects,2);

if measuredNose
    % Filter points with nose at border
    for i = 1:nframes
        if isnan(Nose(i,1)) || isnan(Nose(i,2))
            tracked_frames(i) = 0;
            continue
        end
        
        if Nose(i,1) < Settings.min_nose_dist || ...
                Nose(i,2) > frame_heigth - Settings.min_nose_dist
            tracked_frames(i) = 0;
            continue
        end
        
        if Nose(i,2) < Settings.min_nose_dist || ...
                Nose(i,2) > frame_width - Settings.min_nose_dist
            tracked_frames(i) = 0;
        end
        
    end
    
    
    
    % Filter frames too large acceleration on nose
    da = diff(Nose,1);
    da = sqrt( sum( da.^2,2));
    da = medfilt1(da,3);
    da(end+1) = 0;
    tracked_frames(da > Settings.max_accel) = 0;
    
    % Filter frames with too few neighbouring tracked frames
    frame_density = conv2(tracked_frames, ones(1,10)./10,'same');
    tracked_frames(frame_density < Settings.frame_density) = 0;
end
%
% Fill small frame gaps


for i = 1:length(tracked_frames)
    if tracked_frames(i) == 0
        idx = find( tracked_frames(i+1:end), 1, 'first');
        if idx <= Settings.max_gap_fill
            tracked_frames(i:i+idx) = 1;
        end
    end
end




Traces = Tracker.Traces;
Parameters = Tracker.Parameters;
% Create boolean cell with tag to keep/reject per trace
Valid_trace = cell(1,nframes);
for i = 1:nframes
    if tracked_frames(i)
        Valid_trace{i} = ones(1, size(Traces{i},2)); % keep all traces by default
    else
        Valid_trace{i} = zeros(1, size(Traces{i},2)); % reject all traces by default
    end
end


[pts(:,1), pts(:,2)] = find(Tracker.Objects);
h = waitbar(0, 'filtering traces');
% Filter based on minimum length and distance from nose
for i = 1:nframes
    
    n_traces = length(Valid_trace{i});
    for j = 1:n_traces
        
        
        
        l = Parameters{i}(j,7);
        d = sqrt(sum( (Parameters{i}(j,9:10)).^2));
        
        
        if l <= Settings.min_trace_length
            Valid_trace{i}(j) = 0;
            
        end
        
        
        
        
        if measuredNose
            if d > Settings.max_dist_trace_nose
                Valid_trace{i}(j) = 0;
            end
        end
        
    end
    
    
    % Check for duplicate traces
    
    duplicates = zeros(n_traces, n_traces);
    for idx_1 = 1:n_traces
        for idx_2 = idx_1+1:n_traces
            
            trace_1 = Traces{i}{idx_1};
            trace_2 = Traces{i}{idx_2};
            
            % set shortest trace on 2nd dimension
            if size(trace_1,1) < size(trace_2,1)
                trace_dist = sqrt( (trace_1(:,1)'- trace_2(:,1)).^2 + ...
                    (trace_1(:,2)'-trace_2(:,2)).^2);
                
            elseif size(trace_1,1) >= size(trace_2,1)
                trace_dist = sqrt( (trace_2(:,1)'- trace_1(:,1)).^2 + ...
                    (trace_2(:,2)'-trace_1(:,2)).^2);
                
            end
            
            trace_dist_mean = sum( min( trace_dist, [], 1)) / size(trace_dist,1);
            
            if trace_dist_mean < 1
                duplicates(idx_1, idx_2) = 1;
               % duplicates(idx_2, idx_1) = 1;
            end
            
        end
    end
    
    for j = 1:length(Valid_trace{i})
        if any( duplicates( j, :))
            matched_traces = [j, find(duplicates(j,:))];
            
            %select longest of matched traces
            length_traces = [];
            for k = 1:length(matched_traces)
                length_traces(k) = size(Traces{i}{matched_traces(k)},1);
            end
            [~, maxidx] = max(length_traces);
            
            for k = 1:length(length_traces)
                if k ~= maxidx
                    Valid_trace{i}(matched_traces(k)) = 0;
                end
            end
        end
    end
    
    waitbar(i/nframes)
end
close(h)

% Create a new trace set with fits of the measured traces
Traces_clean = cell(nframes,1);
Parameters_clean = cell(1, nframes);

h = waitbar(0,'Fitting polynomials...');


for i = 1:nframes
    
  
    
    
    nclean = numel(find(Valid_trace{i}));
    idx = find(Valid_trace{i});
    Tsave = cell(1, nclean);
    Psave = zeros(nclean, size(Parameters{i},2));
    
    for j = 1:length(idx)
        trace = Traces{i}{idx(j)};
        params = Parameters{i}(idx(j),:);
        
        % old method fitting with splines
       %fitax = linspace(1, size(trace,1), Settings.n_spline_samples);
        %tsplx = spline( 1:size(trace,1), trace(:,1), fitax);
        %tsply = spline( 1:size(trace,1), trace(:,2), fitax);                
        %Tsave{j} = [tsplx; tsply]';   
        
        t1 = [];
        t2 = [];
        t3 = [];
        
        % fit 1/2/3 degree polynomial to find which has least error
        fitax = [1:length(trace)]';
        px1 = polyfit(fitax, trace(:,1), 1);
        py1 = polyfit(fitax, trace(:,2), 1);
        t1(:,1) = polyval(px1, fitax);
        t1(:,2) = polyval(py1, fitax);
        
        px2 = polyfit(fitax, trace(:,1), 2);
        py2 = polyfit(fitax, trace(:,2), 2);
        t2(:,1) = polyval(px2, fitax);
        t2(:,2) = polyval(py2, fitax);
        
        px3 = polyfit(fitax, trace(:,1), 3);
        py3 = polyfit(fitax, trace(:,2), 3);
        t3(:,1) = polyval(px3, fitax);
        t3(:,2) = polyval(py3, fitax);
        
        r(1) = getRSlocal(trace(:,1), t1(:,1), 1) + getRSlocal(trace(:,2), t1(:,2), 1);
        r(2) = getRSlocal(trace(:,1), t2(:,1), 2) + getRSlocal(trace(:,2), t2(:,2), 2);
        r(3) = getRSlocal(trace(:,1), t3(:,1), 3) + getRSlocal(trace(:,2), t3(:,2), 3);
        
        [~, keep_idx] = max(r);
        
        switch(keep_idx)
            case 1
                Tsave{j} = t1;
            case 2
                Tsave{j} = t2;
            case 3 
                Tsave{j} = t3;
        end
        
    
        Psave(j,:) = params;
        
    end
    
    Traces_clean{i} = Tsave;
    Parameters_clean{i} = Psave;
    
    waitbar(i/nframes)
end

close(h)
end


%%


function rscor = getRSlocal(Y, Yfit, fitdeg)
SSres = sum( (Y-Yfit).^2);
SStot = sum( (Y-mean(Y)).^2);
Rs = 1 - SSres/SStot;
rscor = 1-(1-Rs)*(length(Y)-1)/(length(Y)-fitdeg-1);
end


