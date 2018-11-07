function findResolution()
% Find video resolution 
%%
path = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
load(fullfile(path, 'gapwidths.mat'))

h = waitbar(0, 'Loading data...');
for i = 1:size(gapwidth, 1)
    fname = gapwidth{i,1}(1:end-4);
    width = gapwidth{i,2};
    
    if exist(fullfile(path,[fname '_compiled.mat']),'file')
        load(fullfile(path, [fname '_compiled.mat']))
        npixgap = abs(Annotations.Tracker.gapinfo.edge_1 - Annotations.Tracker.gapinfo.edge_2);
        resolution(i) = width/npixgap;
    else
        resolution(i) = NaN;
    end
    waitbar(i/size(gapwidth, 1))
end
close(h)

fprintf('Resolution is %2.8f mm/pixel\n', mean(resolution,'omitnan'))
    
    