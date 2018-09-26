function Traces_clean  = CleanJanelia(Janelia)
%%
%clc
%clear
%load('E:\Studie\Stage Neurobiologie\Videos\Human clicked\Mouse 46\R1\Data_15_Annotations_Janelia.mat')
%vidfile = 'E:\Studie\Stage Neurobiologie\Videos\Human clicked\Mouse 46\R1\Data_15.dat';

%%
Settings = makeAnalyseSettings;
nframes = size(Janelia.Traces,1);
Nose = Janelia.Nose;
frame_heigth = size(Janelia.Objects, 1);
frame_width = size(Janelia.Objects, 2);

%%
tracked_frames = ones(1,nframes);
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
        continue
    end
    
    
    
end


da = diff(Nose, 1);
da = sqrt( sum( da.^2,2));
da = medfilt1(da, 3);
da(end+1) = 0;
tracked_frames(da > Settings.max_accel) = 0;

frame_density = conv2(tracked_frames, ones(1,10)./10,'same');
tracked_frames(frame_density < Settings.frame_density) = 0;

for i = 1:length(tracked_frames)
    if tracked_frames(i) == 0
        idx = find( tracked_frames(i+1:end), 1 ,'first');
        if idx <= Settings.max_gap_fill
            tracked_frames(i:i+idx) = 1;
        end
    end
end

Traces = Janelia.Traces;

%%
Parameters = Janelia.Parameters;

min_length = 30;
max_length = 300;

Valid_trace = cell(1,nframes);
for i = 1:nframes
    if tracked_frames(i)
        Valid_trace{i} = ones(1, size(Traces{i}, 2));
    else
        Valid_trace{i} = zeros(1, size(Traces{i}, 2));
    end
end
[pts(:,1), pts(:,2)] = find(Janelia.Objects);
object_idx = sub2ind(size(Janelia.Objects),pts(:,1), pts(:,2));

h = waitbar(0, 'filtering traces');
for i = 1:nframes
    
    n_traces = length(Valid_trace{i});
   for j = 1:n_traces
       
       t = round(Traces{i}{j});  
       
       if isempty(t) | any(any(t == 0))
           Valid_trace{i}(j) = 0;
           continue
       end
       
       trace_idx = sub2ind( size(Janelia.Objects), t(:,1), t(:,2));
       o_count = sum(Janelia.Objects(trace_idx));
       
       if o_count/length(trace_idx) > 0.1
           Valid_trace{i}(j) = 0;
       end
       
       
       dist_to_nose_1 = sqrt( sum( (t(1,:)-Nose(i,:)).^2));
       dist_to_nose_2 = sqrt( sum( (t(end,:) - Nose(i,:)).^2));
       
       if dist_to_nose_2 < dist_to_nose_1
           t = flip(t,1);
           Traces{i}{j} = t;
           dist_to_nose = dist_to_nose_2;
       else
           dist_to_nose = dist_to_nose_1;
       end
       
       
       dt = diff(t,1);
       l = sum(sqrt( sum( dt.^2, 2)));
       
       if l<= min_length | l>= max_length
           Valid_trace{i}(j) = 0;
       end
       
     
       if dist_to_nose > Settings.max_dist_trace_nose
           Valid_trace{i}(j) = 0;
       end
    
     
       
   end   
   
   
   % check for duplicate traces
    duplicates = zeros(n_traces, n_traces);
    for idx_1 = 1:n_traces
        for idx_2 = idx_1+1:n_traces
            
            trace_1 = Traces{i}{idx_1};
            trace_2 = Traces{i}{idx_2};
            
            if isempty(trace_1) | isempty(trace_2) |...
                    Valid_trace{i}(idx_1) == 0 | Valid_trace{i}(idx_2) == 0
                continue
            end
            
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







%%

Traces_clean = cell(nframes,1);
Parameters_clean = cell(1, nframes);

h = waitbar(0, 'Fitting spline...');
for i = 1:nframes
    nclean = numel(find( Valid_trace{i}));
    idx = find(Valid_trace{i});
    Tsave = cell(1, nclean);
    Psave = zeros(nclean, size(Parameters{i},2));
    
    for j = 1:length(idx)
        trace = Traces{i}{idx(j)};
        params = Parameters{i}(idx(j), :);
        
        fitax = linspace(1, size(trace,1), Settings.n_spline_samples);
        tsplx = spline( 1:size(trace,1), trace(:,1), fitax);
        tsply = spline( 1:size(trace,1), trace(:,2), fitax);
                
        Tsave{j} = [tsplx; tsply]';
        Psave(j,:) = params;
    end
    
    Traces_clean{i} = Tsave;
    Parameters_clean{i} = Psave;
    
    waitbar(i/nframes)
end
close(h)
%%

%{
Raw_Traces = Traces;
Settings.Video = [metafile(1:end-4) '.dat'];
Settings.Video_width = Janelia.MetaData.Dims(1);
Settings.Video_heigth = Janelia.MetaData.Dims(2);

% Setup figure
fig = figure(1);
set(gcf,'position',[100 100 round(1*Settings.Video_heigth) ...
    round(1*Settings.Video_width)]);
set(gcf,'Units','pixels')
set(gca,'Units','normalized')
set(gca,'Position',[0 0 1 1])
ax = gca;
colormap(ax,'gray')
for j = 1210:nframes
    
    if any(isnan(Nose(j,:)))
        continue
    end
    
    Settings.Current_frame = j;
    frame = LoadFrame(Settings);
    imagesc(ax, frame)
    
    hold(ax, 'on')
    
    
    
    
    for k = 1:size(Raw_Traces{j}, 2)
        t = round(Raw_Traces{j}{k});
        if ~isempty(t)
            plot(t(:,2), t(:,1), 'r')
            text(double(t(end,2)), double(t(end,1)), num2str(k))
        end
    end
    for k = 1:size(Traces_clean{j}, 2)
        t = Traces_clean{j}{k};
        plot(t(:,2), t(:,1), 'g')
        
    end
    
    
    
    hold(ax, 'off')
    drawnow
    
    
end
%}
