clear
close all


print_manual_param = 0;
print_tracker_param = 1;
print_manual_touch = 0;

print_fig = 1;
param_show = 6;

% Load files
datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
Files = dir(fullfile(datapath,'*_Annotations_Tracker.mat'));

if ~exist(fullfile(datapath,'THETA'),'file')
    mkdir( fullfile( datapath, 'THETA'))
end


for file_idx =  18%: size(Files,1)
    close all
    loadfile = fullfile(Files(file_idx).folder, Files(file_idx).name);
    load(loadfile)
    metafile = [loadfile(1:end-24) '.mat'];
    mData = load(metafile);
    Settings = makeAnalyseSettings;
    colors = makeColor();
    c1 = colors.manual_dark;
    c2 = colors.tracker_dark;
    
    
    % Get time axis
    time_ax = mData.Data.Time(:,1)*1000 + mData.Data.Time(:,2); % time in ms
    time_ax = time_ax - time_ax(1);
    
    % Find frames of interest
    nframes = size(Manual.Traces,1);
    tracked_frames = zeros(1, nframes);
    for i = 1:nframes
        if ~isempty(Manual.Traces{i})
            tracked_frames(i) = 1;
        end
    end
    x1 = find(tracked_frames,1,'first');
    x2 = find(tracked_frames,1,'last');
    
    
    % Setup parameter extraction
    manual_labels = Manual.Label_names;
    nlabels = size(manual_labels, 2);
    param_manual(1:nframes,1:nlabels) = NaN;
    param_tracker = [];
    
    
    
    % Extract parameter
    for i = x1:x2
        ntraces = size(Manual.Parameters{i},1);
        for j = 1:ntraces
            if isempty(Manual.Labels{i})
                continue
            end
            idx = find(strcmp(Manual.Labels{i}{j}, manual_labels));
            param_manual(i,idx) = Manual.Parameters{i}(j,param_show);
        end
        
        
        
        
        ntraces = size(Tracker.Parameters_clean{i},1);
        if ntraces> 0
            param_tracker(end+1:end+ntraces, 1:3) =[ ones(ntraces,1)*time_ax(i), Tracker.Parameters_clean{i}(:,param_show) ,ones(ntraces,1)*i];
        end
    end
    
    if isempty(param_tracker)
        disp('No parameters found:')
        disp(Files(file_idx).name)
        continue
    end
    
    % Visualize result
    f = figure;
    f.Units = 'normalized';
    f.Position = [0 0 1 1];
    hold on    
    if print_tracker_param
      
            scatter(param_tracker(:,1), param_tracker(:,2), 'MarkerFaceColor',c2,'MarkerEdgeColor',c2)
        
    end    
    if print_manual_param
        plot(time_ax(x1:x2),param_manual(x1:x2,:),'r','LineWidth',1.5)
    end    
    ylim([-180 180])    
    if print_manual_touch
        for i = x1:x2
            if ~isempty(Manual.Touch.pt{i})
                line([time_ax(i) time_ax(i)],[-600 600],'color','g')
                for j =1 :size(Manual.Touch.label{i},2)
                    disp(Manual.Touch.label{i}{j})
                end
            end            
        end
    end
    
    
    if print_fig
        saveas(gcf, fullfile(datapath,'THETA',[Files(file_idx).name(1:end-24) '_raw.png']))
    end
    
    
    
    Par = Tracker.Parameters_clean;
    
    boundary = 0;
    
    
    Par_left = [];
    Par_right = [];
    for i = x1:x2
        
        if isempty( Par{i} )
            Par_left( : , i-x1+1) = NaN;
            Par_right( :, i-x1+1) = NaN;
            continue
        end
        
        
        
        par_loop = Par{i}(:, param_show);
        
        idx_left = find(par_loop <= boundary);
        idx_right = find(par_loop > boundary);
        
        
        par_left_size(i-x1+1) = length(idx_left);
        Par_left( 1: length(idx_left), i-x1+1 ) = par_loop( idx_left );
        
        par_right_size(i-x1+1) = length( idx_right);
        Par_right( 1: length( idx_right), i-x1+1 ) = par_loop( idx_right );
        
        
    end
    
    for i = x1:x2
        
        if i-x1+1 <= length(par_left_size) &&  par_left_size(i-x1+1) > 0
            Par_left( par_left_size(i-x1+1)+1 :end, i-x1+1) = NaN;
        else
            Par_left(:, i-x1+1) = NaN;
        end
        
        if i-x1+1 <= length(par_right_size) && par_right_size(i-x1+1) > 0
            Par_right( par_right_size(i-x1+1)+1:end, i-x1+1) = NaN;
        else
            Par_right( :, i-x1+1) = NaN;
        end
    end
    
    
    filtersize = 5;
    
    mean_left = mean(Par_left, 1, 'omitnan');
    mean_right = mean(Par_right, 1, 'omitnan');
    
    mean_left_filt = medfilt1( mean_left, filtersize);
    mean_right_filt = medfilt1( mean_right, filtersize);
    
    plot( time_ax(x1:x2), mean_left_filt, 'b','LineWidth',2)
    plot( time_ax(x1:x2), mean_right_filt, 'b', 'LineWidth', 2)
    
    
    
    
    
    
    nstartpoints = numel(find(~isnan( Par_left(:,1))));
    
    
    
    TFILT_LEFT = {};
    PL = Par_left;
    value_trace = zeros(1, size(PL, 2));
    same_val = 3;
    same_val_clean = 3;
    
    while any(any(~isnan(PL)))
        value_trace =[];
        value_trace_filt = [];
        
        % Assign lowest values in array as first trace
        for i = 1: size(PL, 2)
            idx = find( PL(:,i) <= min( PL(:,i) ) + same_val );
            value_trace(i) = mean( PL(idx,i));
        end
        
        % filter signal
        xax = time_ax(x1:x2);
        f_width = 2;
        t_const = median( diff( time_ax ));
        f_width_time = f_width * t_const;
        t_ax = time_ax( x1:x2 );
        
        for i = 1:length( value_trace )
            idx = find( t_ax > t_ax(i) - f_width_time & t_ax < t_ax(i) + f_width_time);
            value_trace_filt(i) = mean( value_trace( idx ), 'omitnan' );
        end
        
        TFILT_LEFT{end+1} = value_trace_filt;
        
        % Clean parameter data
        for i = 1:length( value_trace_filt)
            idx_omit = find( PL(:,i) < value_trace_filt(i) - same_val_clean );
            idx_assign = find( PL(:,i) >= value_trace_filt(i) - same_val_clean & ...
                PL(:,i) < value_trace_filt(i) + same_val_clean);
            
            PL(idx_omit, i) = NaN;
            PL(idx_assign, i) = NaN;
            
        end
        
    end
    
    
    PR = Par_right;
    TFILT_RIGHT = {};
    while any(any(~isnan(PR)))
        value_trace = [];
        
        % Assign lowest values in array as first trace
        for i = 1: size( PR, 2)
            idx = find( PR(:,i) >= max( PR(:,i) ) - same_val );
            value_trace(i) = mean( PR(idx,i));
        end
        
        % filter signal
        xax = time_ax(x1:x2);
        f_width = 2;
        t_const = median( diff( time_ax ));
        f_width_time = f_width * t_const;
        t_ax = time_ax( x1:x2 );
        
        for i = 1:length( value_trace )
            idx = find( t_ax > t_ax(i) - f_width_time & t_ax < t_ax(i) + f_width_time);
            value_trace_filt(i) = mean( value_trace( idx ), 'omitnan' );
        end
        
        TFILT_RIGHT{end+1} = value_trace_filt;
        
        % Clean parameter data
        for i = 1:length( value_trace_filt)
            idx_omit = find( PR(:,i) > value_trace_filt(i) + same_val_clean );
            idx_assign = find( PR(:,i) >= value_trace_filt(i) - same_val_clean & ...
                PR(:,i) < value_trace_filt(i) + same_val_clean);
            
            PR(idx_omit, i) = NaN;
            PR(idx_assign, i) = NaN;
            
        end
        
    end
    
    
    
    
    % Visualize result
    f = figure;
    f.Units = 'normalized';
    f.Position = [0 0 1 1];
    hold on
    cc = hsv(size(TFILT_LEFT,2)); 
    for i = 1:size( PL, 2)
        scatter( ones(1, size(PL,1))*xax(i), PL(:, i), 'b', 'filled')
    end
    for i = 1:size(TFILT_LEFT, 2)
        plot(xax, TFILT_LEFT{i}, 'color', cc(i,:) )
    end    
    cc = hsv(size(TFILT_RIGHT,2));    
    for i = 1:size(TFILT_RIGHT, 2)
        plot(xax, TFILT_RIGHT{i}, 'color', cc(i,:) )
    end
    
    if print_fig
        saveas(gcf, fullfile(datapath,'THETA',[Files(file_idx).name(1:end-24) '_colored.png']))
    end
    
    
    f = figure;
    f.Units = 'normalized';
    f.Position = [0 0 1 1];
    hold on
    cc = hsv(size(TFILT_LEFT,2));    
    for i = 1:size(TFILT_LEFT, 2)
        t = TFILT_LEFT{i};
        r = numel(find( isnan(t)))/ numel(t);
        if r < 0.3
            plot(xax, TFILT_LEFT{i}, 'color', 'b','LineWidth',1.5)
        end
    end    
    cc = hsv(size(TFILT_RIGHT,2));    
    for i = 1:size(TFILT_RIGHT, 2)
        t = TFILT_RIGHT{i};
        r = numel(find( isnan(t)))/ numel(t);
        if r < 0.3
            plot(xax, TFILT_RIGHT{i}, 'color', 'b','LineWidth',1.5)
        end
    end    
    % scatter(param_tracker(:,1), param_tracker(:,2), 'MarkerFaceColor',c2,'MarkerEdgeColor',c2)    
    plot(time_ax(x1:x2),param_manual(x1:x2,:),'r','LineWidth',1.5)
    
    if print_fig
        saveas(gcf, fullfile(datapath,'THETA',[Files(file_idx).name(1:end-24) '_filtered.png']))
    end
    
    
    
    
end
















