function Output = ConvertJanelia(filename)
%%
%filename = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02.whiskers';
data = LoadWhiskers(filename);
nframes = data(end).time + 1;
Traces = cell(nframes, 1);

for i = 1:size(data,1)
   Traces{data(i).time+1}{data(i).id+1} = [data(i).y, data(i).x]; 
end

Output.Traces = Traces;
