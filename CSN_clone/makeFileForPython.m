function makeFileForPython(infile)
%%
load(infile)

%% get touches
T = Annotations.Tracker.Touch;
touch = zeros(size(T,2),1);
sensidx = 4;
for idx = 1:size(T,2)
    touch(idx) = sum(T{sensidx,idx});
 
end

%% get noseposition
NoseX = Annotations.Tracker.Nose(:,1);
NoseY = Annotations.Tracker.Nose(:,2);
AngleX = Annotations.Tracker.Headvec(:,1);
AngleY = Annotations.Tracker.Headvec(:,2);
Angle = atand(AngleY./AngleX);

%%
res = split(infile, '\');
outpath = res{1};
for i = 2:size(res,1)-1
    outpath = fullfile(outpath, res{i});
end
disp(outpath)
save( fullfile(outpath, 'MWTforPython'), 'touch','NoseX','NoseY','AngleX','AngleY', 'Angle')
