function dist = finddist(Point,Trace)
% dist = finddist(Point, Trace)
%
% Returns shortest distance between a point and a trace


npoints = size(Trace,1);
dist(1:npoints) = NaN;
for i = 1:npoints
    dist(i) = sqrt( sum( (Point' - Trace(i,:)).^2));
end

dist = min(dist);