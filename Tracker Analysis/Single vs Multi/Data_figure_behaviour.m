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

Data.ax1.multi_width = t_filtered.Width( strcmp(t_filtered.Type, 'Multi'));
Data.ax1.multi_touch = t_filtered.sum_nTouch( strcmp(t_filtered.Type, 'Multi'));


%%

% Data.fit.single.xax = Data.table_sorted_1.mean_Width(sidx);
% Data.fit.single.data = Data.table_sorted_1.GroupCount(sidx);
% Data.fit.single.p = polyfit( Data.fit.single.xax, Data.fit.single.data, 1);
% Data.fit.single.showax = 25:55;
% Data.fit.single.showfit = polyval(Data.fit.single.p, Data.fit.single.showax);
% Data.fit.single.Rfit = polyval(Data.fit.single.p, Data.fit.single.xax);
% Data.fit.single.SSE = sum( ( Data.fit.single.data - mean(Data.fit.single.data)).^2);
% Data.fit.single.SSTO = sum( (Data.fit.single.data - Data.fit.single.Rfit).^2);
% Data.fit.single.RS = abs(1 - Data.fit.single.SSE/ Data.fit.single.SSTO);
% 
% 
% Data.fit.multi.xax = Data.table_sorted_1.mean_Width(midx);
% Data.fit.multi.data = Data.table_sorted_1.GroupCount(midx);
% Data.fit.multi.p = polyfit( Data.fit.multi.xax, Data.fit.multi.data, 1);
% Data.fit.multi.showax = 25:55;
% Data.fit.multi.showfit = polyval(Data.fit.multi.p, Data.fit.multi.showax);
% Data.fit.multi.Rfit = polyval(Data.fit.multi.p, Data.fit.multi.xax);
% Data.fit.multi.SSE = sum( (Data.fit.multi.data - mean(Data.fit.multi.data)).^2);
% Data.fit.multi.SSTO  = sum( (Data.fit.multi.data - Data.fit.multi.Rfit).^2);
% Data.fit.multi.RS = abs(1-Data.fit.multi.SSE/ Data.fit.multi.SSTO);
% 

%% AX2  

Data.ax2.single_width = t_filtered.Width( strcmp(t_filtered.Type, 'Single'));
Data.ax2.single_duration = t_filtered.Duration( strcmp(t_filtered.Type, 'Single'));

Data.ax2.multi_width = t_filtered.Width(strcmp( t_filtered.Type, 'Multi'));
Data.ax2.multi_duration = t_filtered.Duration( strcmp(t_filtered.Type, 'Multi'));


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





