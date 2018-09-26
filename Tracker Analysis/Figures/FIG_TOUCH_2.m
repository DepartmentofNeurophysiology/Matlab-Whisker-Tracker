clear


datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
widthfile = 'gapwidths.mat';

gw = load( fullfile( datapath, widthfile) );

Files = dir(fullfile(datapath,'*_Annotations_Tracker.mat'));

min_nr_frames = 100;



ncontacts_tracker = [];

ncontacts_manual = [];
resolution(1:size(Files,1)) = NaN;


h = waitbar(0, 'Loading data...');

for i = 1:size(Files,1)
    
    
    loadfile = fullfile(Files(i).folder, Files(i).name);
    load(loadfile)    
    resolution(i) = gw.gapwidth{i,2} / (Tracker.gapinfo.edge_2 - Tracker.gapinfo.edge_1);
    waitbar(i/size(Files,1))
end
res_mm = mean(resolution);
%%
for i = 1:size(Files,1)
     loadfile = fullfile(Files(i).folder, Files(i).name);
    load(loadfile)  
    for j = 1:size(Tracker.Touch,2)
        
        if isnan(Tracker.dist_nose_target(j))
            continue
        end
        
        if ~isempty(Tracker.Touch{j})
        n_touch_frame = sum(Tracker.Touch{j});
        else
        n_touch_frame = 0;
        end
        if n_touch_frame ~= 0
            ncontacts_tracker(end+1,1:2) = [abs(Tracker.dist_nose_target(j))*res_mm, n_touch_frame];
            
        end
    end
    
    for j = 1:size(Manual.Touch.pt,2)
        if isnan(Tracker.dist_nose_target(j))
            continue
        end
        n_touch_frame = 0;
        if ~isempty(Manual.Touch.pt{j})
            n_touch_frame = size( Manual.Touch.pt{j},1);
        end
        
        if n_touch_frame ~= 0
            ncontacts_manual(end+1,1:2) = [abs(Tracker.dist_nose_target(j))*res_mm, n_touch_frame];
        end
    end
    
    waitbar(i/size(Files,1))
    
end


%%
Tracker.distance = ncontacts_tracker(:,1);
Tracker.ncontacts = ncontacts_tracker(:,2);
Manual.distance = ncontacts_manual(:,1);
Manual.ncontacts = ncontacts_manual(:,2);

max_dist = max([Tracker.distance; Manual.distance]);
max_contacts = 8;% max([Tracker.ncontacts; Manual.ncontacts]);

binsize_dist = 2.5;
bin_edges = 0:binsize_dist: max_dist-mod(max_dist, binsize_dist) + binsize_dist;
Tracker.distance = discretize( Tracker.distance, bin_edges );
Manual.distance = discretize( Manual.distance, bin_edges );

max_dist_bin = max([Tracker.distance; Manual.distance]);

tracker_table = table(Tracker.distance, Tracker.ncontacts);
tracker_table_sorted = varfun(@GroupCount, tracker_table,'GroupingVariables',{'Var1','Var2'});

tracker_image = zeros(max_contacts, max_dist_bin);
for i = 1:size(tracker_table_sorted)
    if tracker_table_sorted.Var2(i) > max_contacts
        continue
    end
    
    tracker_image( tracker_table_sorted.Var2(i), ...
        tracker_table_sorted.Var1(i)) = ...
        tracker_table_sorted.GroupCount(i);
    
end
tracker_image = tracker_image ./ sum( sum( tracker_image));
tracker_image = flip( flip( tracker_image,1),2);

f = figure(1);
clf(f);
ax = axes(f);

imagesc(tracker_image);
colormap gray
caxis([0 0.25])

xtag = (0:max_dist_bin)*binsize_dist;
xpos = (0:max_dist_bin)+0.5;
ax.XTick = xpos;
ax.XTickLabel = flip(xtag);
xlabel('dist from target (mm)');

ytag = 1:max_contacts;
ypos = 1:max_contacts;
yidx = 2:2:length(ypos);


ax.YTick = ypos(yidx-1);
ax.YTickLabel = flip(ytag(yidx));

title('Tracker data')




%

manual_table = table(Manual.distance, Manual.ncontacts);
manual_table_sorted = varfun(@GroupCount, manual_table,'GroupingVariables',{'Var1','Var2'});

manual_image = zeros(max_contacts, max_dist_bin);
for i = 1:size(tracker_table_sorted)
    if manual_table_sorted.Var2(i) > max_contacts
        continue
    end
    
    manual_image(manual_table_sorted.Var2(i), ...
        manual_table_sorted.Var1(i)) = ...
        manual_table_sorted.GroupCount(i);
    
end
manual_image = manual_image ./ sum( sum( manual_image));
manual_image = flip( flip( manual_image,1),2);

f = figure(2);
clf(f);
ax = axes(f);

imagesc(manual_image);
colormap gray
caxis([0 0.25])

xtag = (0:max_dist_bin)*binsize_dist;
xpos = (0:max_dist_bin)+0.5;
ax.XTick = xpos;
ax.XTickLabel = flip(xtag);
xlabel('dist from target (mm)');

ytag = 1:max_contacts;
ypos = 1:max_contacts;
yidx = 2:2:length(ypos);


ax.YTick = ypos(yidx-1);
ax.YTickLabel = flip(ytag(yidx));


title('Manual data')






























