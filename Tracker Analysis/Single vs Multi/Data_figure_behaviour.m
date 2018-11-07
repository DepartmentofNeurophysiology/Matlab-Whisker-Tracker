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

t = table(Dist,nTouch,Width,Video,Type,'VariableNames',{'Dist','nTouch','Width','Video','Type'});

t = t(t.Width > 15,:);


a = varfun(@mean, t, 'InputVariables','Width','GroupingVariables',{'Type','Video'});
sidx=  find(strcmp(a.Type,'Single'));
midx = find(strcmp(a.Type,'Multi'));


Data.table = t;
Data.table_sorted_1 = a;


sidx = find(strcmp(Data.table_sorted_1.Type, 'Single') & Data.table_sorted_1.GroupCount < 1000);
midx = find(strcmp(Data.table_sorted_1.Type, 'Multi') & Data.table_sorted_1.GroupCount < 1000);

Data.show.single.xax = Data.table_sorted_1.mean_Width(sidx);
Data.show.single.data = Data.table_sorted_1.GroupCount(sidx);
Data.show.multi.xax = Data.table_sorted_1.mean_Width(midx);
Data.show.multi.data = Data.table_sorted_1.GroupCount(midx);



Data.fit.single.xax = Data.table_sorted_1.mean_Width(sidx);
Data.fit.single.data = Data.table_sorted_1.GroupCount(sidx);
Data.fit.single.p = polyfit( Data.fit.single.xax, Data.fit.single.data, 1);
Data.fit.single.showax = 25:55;
Data.fit.single.showfit = polyval(Data.fit.single.p, Data.fit.single.showax);
Data.fit.single.Rfit = polyval(Data.fit.single.p, Data.fit.single.xax);
Data.fit.single.SSE = sum( ( Data.fit.single.data - mean(Data.fit.single.data)).^2);
Data.fit.single.SSTO = sum( (Data.fit.single.data - Data.fit.single.Rfit).^2);
Data.fit.single.RS = abs(1 - Data.fit.single.SSE/ Data.fit.single.SSTO);


Data.fit.multi.xax = Data.table_sorted_1.mean_Width(midx);
Data.fit.multi.data = Data.table_sorted_1.GroupCount(midx);
Data.fit.multi.p = polyfit( Data.fit.multi.xax, Data.fit.multi.data, 1);
Data.fit.multi.showax = 25:55;
Data.fit.multi.showfit = polyval(Data.fit.multi.p, Data.fit.multi.showax);
Data.fit.multi.Rfit = polyval(Data.fit.multi.p, Data.fit.multi.xax);
Data.fit.multi.SSE = sum( (Data.fit.multi.data - mean(Data.fit.multi.data)).^2);
Data.fit.multi.SSTO  = sum( (Data.fit.multi.data - Data.fit.multi.Rfit).^2);
Data.fit.multi.RS = abs(1-Data.fit.multi.SSE/ Data.fit.multi.SSTO);




%%
% bins for distance binarization
dsize=  2;
dist_bins = -25:dsize:50;
t.Dist_bin = discretize(t.Dist, dist_bins);



res = [];
for i = 1:length(dist_bins)
    ids = find(t.Dist_bin == i & strcmp(t.Type,'Single'));
    idm = find(t.Dist_bin == i & strcmp(t.Type,'Multi'));
    res(i,1) = dist_bins(i) + dsize/2;
    res(i,2) = mean(t.nTouch(ids));
    res(i,4) = std(t.nTouch(ids));
    res(i,3) = mean(t.nTouch(idm));
    res(i,5) = std(t.nTouch(idm));
end
res(isnan(res)) = 0;
res(:,4) = res(:,4)./sum(res(:,2));
res(:,5) = res(:,5)./sum(res(:,3));
res(:,2) = res(:,2)./sum(res(:,2));
res(:,3) = res(:,3)./sum(res(:,3));

Data.dist_v_ntouch = res;


res = [];
for i = 1:length(dist_bins)
     ids = find(t.Dist_bin == i & strcmp(t.Type,'Single'));
    idm = find(t.Dist_bin == i & strcmp(t.Type,'Multi'));
    res(i,1) = dist_bins(i) + dsize/2;
    res(i,2) = numel(find(ids));
    res(i,3) = numel(find(idm));
end
res(:,2) = res(:,2)./sum(res(:,2));
res(:,3) = res(:,3)./sum(res(:,3));

Data.dist_v_count = res;

save(fullfile(filepath, 'DataFigBeh'), 'Data')





