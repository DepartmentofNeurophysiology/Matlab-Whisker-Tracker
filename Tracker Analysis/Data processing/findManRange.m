function range = findManRange(Manual)

t = Manual.Traces;
tracked = zeros(1, size(t,1));
for i = 1:size(t, 1)
    if ~isempty(t{i})
       tracked(i) = 1; 
    end
end
range(1) = find(tracked, 1 ,'first');
range(2) = find(tracked, 1, 'last');