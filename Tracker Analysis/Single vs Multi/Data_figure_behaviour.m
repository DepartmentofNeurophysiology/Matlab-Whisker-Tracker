filepath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi';
files = dir(fullfile(filepath, '*compiled.mat'));

multi_idx = 0:39;
single_idx = 40:48;


h = waitbar(0, 'Loading data');


% Table for data storage
table_id = 1;

Dist = [];
nTouch = [];
Video = {};
Type = [];
Width = [];

resolution = 0.128; % mm/pxl

%%
for file_idx = 1:size(files, 1)
    
    % Load data
    load(fullfile(files(file_idx).folder, files(file_idx).name))
    
    % Extract parameters of interest
    dist = Annotations.Tracker.dist_nose_target; % distance nose to target
    touch = Annotations.Tracker.Touch;           % detected touch
    nose = Annotations.Tracker.Nose;
    traces = Annotations.Tracker.Traces;
    nframes = size(dist, 1); % number of frames
    gapinfo = Annotations.Tracker.gapinfo;
    
  
   
    warning('off')
    
    % Extract data
    for i = 1:nframes
        if ~isnan(nose(i,:)) & ~isempty(traces{i})
            
            % distance to target
            Dist(table_id,1) = dist(i)*resolution;
            
            % number of touches
            nTouch(table_id,1) = numel(find(touch{i} == 1));
   
            
            % video name
            Video{table_id,1} = files(file_idx).name;
            
            % type
            if ismember(str2double(files(file_idx).name(6:7)), multi_idx)
                Type{table_id,1} = 'Multi';
            elseif ismember(str2double(files(file_idx).name(6:7)), single_idx)
                Type{table_id,1} = 'Single';
            end
            
            % Width
            Width(table_id,1) = abs(gapinfo.edge_1 - gapinfo.edge_2)*resolution;
            
            table_id = table_id + 1;
        end
    end
    
    % Store datafrequency
    freq(file_idx) = 1000/median(diff(Annotations.Tracker.MetaData.TimeMS));
    
    waitbar(file_idx/size(files, 1))
end

close(h)



%%
t = table(Dist,nTouch,Width,Video,Type,'VariableNames',{'Dist','nTouch','Width','Video','Type'});
dsize=  2;
dist_bins = -25:dsize:50;
t.Dist_bin = discretize(t.Dist, dist_bins);
for i = 1:size(t, 1)
    if ~isnan(t.Dist_bin(i))
    t.Dist_bin(i) = dist_bins(t.Dist_bin(i)) + dsize/2;
    end
end


t_sorted = varfun(@sum, t, 'InputVariables','nTouch','GroupingVariables',{'Type','Video','Width'});
t_sorted.Properties.VariableNames{4} = 'Duration';
t_sorted.Duration = t_sorted.Duration* (1/ median(freq));
keep_idx = find(...
    t_sorted.sum_nTouch > 0 & ...
    t_sorted.Width > 5 & ...
    t_sorted.Duration < 2);

t_filtered = t_sorted(keep_idx,:);


%% AX 1


Data.table = t;
Data.table_sorted = t_sorted;
Data.table_filtered = t_filtered;

Data.ax1.single_width = t_filtered.Width( strcmp(t_filtered.Type, 'Single'));
Data.ax1.single_touch = t_filtered.sum_nTouch( strcmp(t_filtered.Type, 'Single'));
Data.ax1.single_fit_p = polyfit(Data.ax1.single_width, Data.ax1.single_touch, 1);
Data.ax1.single_fit_ax = 20:60;
Data.ax1.single_fit = polyval(Data.ax1.single_fit_p, Data.ax1.single_width);
Data.ax1.single_fit_plot = polyval(Data.ax1.single_fit_p, Data.ax1.single_fit_ax);
Data.ax1.single_SSE = sum( (Data.ax1.single_touch - mean(Data.ax1.single_touch)).^2);
Data.ax1.single_SSTO = sum( (Data.ax1.single_touch - Data.ax1.single_fit).^2);
Data.ax1.single_RS = abs(1- Data.ax1.single_SSE/ Data.ax1.single_SSTO);

Data.ax1.multi_width = t_filtered.Width( strcmp(t_filtered.Type, 'Multi'));
Data.ax1.multi_touch = t_filtered.sum_nTouch( strcmp(t_filtered.Type, 'Multi'));
Data.ax1.multi_fit_p = polyfit(Data.ax1.multi_width, Data.ax1.multi_touch, 1);
Data.ax1.multi_fit_ax = 20:60;
Data.ax1.multi_fit = polyval(Data.ax1.multi_fit_p, Data.ax1.multi_width);
Data.ax1.multi_fit_plot = polyval(Data.ax1.multi_fit_p, Data.ax1.multi_fit_ax);
Data.ax1.multi_SSE = sum( (Data.ax1.multi_touch - mean(Data.ax1.multi_touch)).^2);
Data.ax1.multi_SSTO = sum( (Data.ax1.multi_touch - Data.ax1.multi_fit).^2);
Data.ax1.multi_RS = abs(1-Data.ax1.multi_SSE/ Data.ax1.multi_SSTO);


%% AX2  

Data.ax2.single_width = t_filtered.Width( strcmp(t_filtered.Type, 'Single'));
Data.ax2.single_duration = t_filtered.Duration( strcmp(t_filtered.Type, 'Single'));
Data.ax2.single_fit_p = polyfit(Data.ax2.single_width, Data.ax2.single_duration, 1);
Data.ax2.single_fit_ax = 20:60;
Data.ax2.single_fit = polyval(Data.ax2.single_fit_p, Data.ax2.single_width);
Data.ax2.single_fit_plot = polyval(Data.ax2.single_fit_p, Data.ax2.single_fit_ax);
Data.ax2.single_SSE = sum( (Data.ax2.single_duration - mean(Data.ax2.single_duration)).^2);
Data.ax2.single_SSTO = sum( (Data.ax2.single_duration - Data.ax2.single_fit).^2);
Data.ax2.single_RS = abs(1-Data.ax2.single_SSE/Data.ax2.single_SSTO);

Data.ax2.multi_width = t_filtered.Width(strcmp( t_filtered.Type, 'Multi'));
Data.ax2.multi_duration = t_filtered.Duration( strcmp(t_filtered.Type, 'Multi'));
Data.ax2.multi_p = polyfit(Data.ax2.multi_width, Data.ax2.multi_duration, 1);
Data.ax2.multi_fit_ax = 20:60;
Data.ax2.multi_fit = polyval(Data.ax2.multi_p, Data.ax2.multi_width);
Data.ax2.multi_fit_plot = polyval(Data.ax2.multi_p, Data.ax2.multi_fit_ax);
Data.ax2.multi_SSE = sum( (Data.ax2.multi_duration - mean(Data.ax2.multi_duration)).^2);
Data.ax2.multi_SSTO = sum( (Data.ax2.multi_duration - Data.ax2.multi_fit).^2);
Data.ax2.multi_RS = abs(1-Data.ax2.multi_SSE/Data.ax2.multi_SSTO);


%%
% bins for distance binarization

t_sorted = varfun(@sum, t, 'InputVariables', 'nTouch', 'GroupingVariables', {'Dist_bin','Type'});

t_sorted.Properties.VariableNames{3} = 'Duration';
t_sorted.Duration =  t_sorted.Duration * (1/median(freq));
t_filtered = t_sorted;

tt = varfun(@mean, t, 'InputVariables', 'Width', 'GroupingVariables',{'Video','Type'});
n_videos_single = numel(find(strcmp(tt.Type, 'Single')));
n_videos_multi = numel(find(strcmp(tt.Type, 'Multi')));



Data.ax3.table = t_sorted;
sidx = find(strcmp(t_filtered.Type, 'Single'));
Data.ax3.single_dist = t_filtered.Dist_bin( sidx);
Data.ax3.single_touch = t_filtered.sum_nTouch( sidx)./ n_videos_single;
Data.ax4.single_dist = t_filtered.Dist_bin( sidx);
Data.ax4.single_duration = t_filtered.Duration(sidx)./ n_videos_single;

midx = find(strcmp(t_filtered.Type, 'Multi'));
Data.ax3.multi_dist = t_filtered.Dist_bin(midx);
Data.ax3.multi_touch = t_filtered.sum_nTouch(midx)./ n_videos_multi;
Data.ax4.multi_dist = t_filtered.Dist_bin(midx);
Data.ax4.multi_duration = t_filtered.Duration(midx)./ n_videos_multi;




save(fullfile(filepath, 'DataFigBeh'), 'Data')





