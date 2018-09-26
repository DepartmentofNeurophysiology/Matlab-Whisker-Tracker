function dist = finddist(Point,Trace)
% Find distance between a point and a trace
npoints = size(Trace,1);

for i = 1:npoints
    dist(i) = sqrt( sum( (Point' - Trace(i,:)).^2));
end

dist = min(dist);