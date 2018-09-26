function Output = ConvertAnnotations(CurvesByFrame)
%%


% Verify number of clicks to process
nclicks = 0;
for i = 1:size(CurvesByFrame, 1)
    nclicks = nclicks + size(CurvesByFrame{i},2);
end


% Find unique whiskers
uniquewhiskers = {};
for i = 1:size(CurvesByFrame , 1)
    
    if ~isempty(CurvesByFrame{i})
        for j = 1:size(CurvesByFrame{i},2)
            label = sprintf('%s%d',CurvesByFrame{i}{j}{5},CurvesByFrame{i}{j}{6});
            
            
            flag = 1;
            for k = 1:length(uniquewhiskers)
                if strcmp(uniquewhiskers{k},label)
                    flag = 0;
                end
            end
            
            if flag
                uniquewhiskers{end+1} = label;
            end
        end
    end
end
nwhiskers = size( uniquewhiskers,2);

Traces = cell(size(CurvesByFrame,1), 1);
Labels = cell(size(CurvesByFrame,1), 1);
Labelsfull = cell(1, size(CurvesByFrame,1));
Clusters = cell(1, size(CurvesByFrame,1));

% Extract clicking data
nprocessed = 0;
for i = 1:size( CurvesByFrame, 1)
    if ~isempty(CurvesByFrame{i})
        Traces{i} = cell(1, nwhiskers);
        tick = 1;
        for j = 1:size(CurvesByFrame{i},2)
            
            switch(CurvesByFrame{i}{j}{4})
                case 'track'
                    label = sprintf('%s%d',CurvesByFrame{i}{j}{5},CurvesByFrame{i}{j}{6});
                    tid = find(strcmp(uniquewhiskers, label));
                    Traces{i}{tid}(end+1,1:2) = [CurvesByFrame{i}{j}{2}, CurvesByFrame{i}{j}{1}];
                    Labels{i}(tid) = label(1);
                    Clusters{i}(tid) = label(2);
                    Labelsfull{i}{tid} = label;
                    
                    nprocessed = nprocessed + 1;
                    
                case 'touch'
                    Touch.label{i}{tick} = sprintf('%s%d',CurvesByFrame{i}{j}{5},CurvesByFrame{i}{j}{6});
                    Touch.pt{i}(tick,:)= [CurvesByFrame{i}{j}{2}, CurvesByFrame{i}{j}{1}];
                    label = sprintf('%s%d',CurvesByFrame{i}{j}{5},CurvesByFrame{i}{j}{6});
                    tick = tick+1;
                    
                    nprocessed = nprocessed + 1;
                    
                otherwise
                    keyboard
            end
            
            
        end
    end
end





if nprocessed ~= nclicks
    keyboard
end

% Fit a quadratic funtion trough the traces
Tracesfit = cell(size(CurvesByFrame, 1), 1);

for i = 1:size(Traces, 1)
    if ~isempty(Traces{i})
        for j = 1:size(Traces{i},2)
            if ~isempty(Traces{i}{j})
                t = Traces{i}{j};
                npts = size(t,1);
                rawax = 1:npts;
                fitax = linspace(1,npts,90);
                
                px = polyfit(rawax, t(:,1)', 2);
                py = polyfit(rawax, t(:,2)', 2);
                
                tfit(:,1) = polyval(px, fitax);
                tfit(:,2) = polyval(py, fitax);
                
                Tracesfit{i}{j} = tfit;
            else
                Tracesfit{i}{j} = [];
            end
        end
        
    end
end


% Find unique labels
unique_labels = {};
for i = 1:size(Labelsfull,2)
    
    for j = 1:size(Labelsfull{i},2)
        l = Labelsfull{i}{j};
        
        if isempty(unique_labels)
            unique_labels{1} = l;
        end
        
        if ~any( strcmp(l ,unique_labels) )
            unique_labels{end+1} = l;
            
        end
        
    end
    
    
end


Output.TracesSmall = Traces;
Output.Traces = Tracesfit;
Output.Labels.Side = Labels;
Output.Labels.Clusters = Clusters;
Output.Labels.Full = Labelsfull;
Output.Labels.Names = unique_labels;
Output.Touch = Touch;

