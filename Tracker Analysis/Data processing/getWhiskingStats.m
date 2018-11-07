function Data = getWhiskingStats(Annotations)
A = getAngles(Annotations);
names = {'r_min','r_max','l_min','l_max'};

Mr_min = [];
Mr_max = [];
Ml_min = [];
Ml_max = [];
for i = 1:length(names)
    
    t = []; %#ok<*NASGU>
    eval(sprintf('T%s = getStats(A.Tracker.%s_filtered);',names{i},names{i} ))
    match = [];
    duration = [];
    amplitude = [];
    if isfield(A, 'Manual') & eval(sprintf('~isempty(T%s)',names{i})) & eval(sprintf('~isnan(T%s.Type(1))',names{i}))
        eval(sprintf('M%s = getStats(A.Manual.%s_filtered);', names{i},names{i}))
        eval(sprintf('tR = matchCycle(T%s.Onset(T%s.Type == 1), M%s(M%s.Type == 1,:));',names{i},names{i},names{i},names{i}));
        eval(sprintf('tP = matchCycle(T%s.Onset(T%s.Type == 2), M%s(M%s.Type == 2,:));',names{i},names{i},names{i},names{i}));
        eval(sprintf('match(T%s.Type == 1,1) = tR.Match;',names{i}))
        eval(sprintf('duration(T%s.Type == 1,1) = tR.Duration;', names{i}))
        eval(sprintf('amplitude(T%s.Type == 1,1) = tR.Amplitude;', names{i}))

        eval(sprintf('match(T%s.Type == 2,1) = tP.Match;',names{i}))
        eval(sprintf('duration(T%s.Type == 2,1) = tP.Duration;',names{i}))
        eval(sprintf('amplitude(T%s.Type == 2,1) = tP.Amplitude;',names{i}))
        if ~isempty(match)
            eval(sprintf('T%s.Match = match;',names{i}))
            eval(sprintf('T%s.MDuration = duration;', names{i}))
            eval(sprintf('T%s.MAmplitude = amplitude;',names{i}))
        else
            eval(sprintf('T%s.Match = NaN;', names{i}))
            eval(sprintf('T%s.MDuration = NaN;',names{i}))
            eval(sprintf('T%s.MAmplitude = NaN;', names{i}))
        end
    else
        eval(sprintf('T%s.Match = NaN;', names{i}))
        eval(sprintf('T%s.MDuration = NaN;',names{i}))
        eval(sprintf('T%s.MAmplitude = NaN;', names{i}))
    end
    
end

Data.Tr_min = Tr_min;
Data.Tr_max = Tr_max;
Data.Tl_min = Tl_min;
Data.Tl_max = Tl_max;

if isfield(A, 'Manual')  
    Data.Mr_min = Mr_min;
    Data.Mr_max = Mr_max;
    Data.Ml_min = Ml_min;
    Data.Ml_max = Ml_max;
end

end


function Res = matchCycle(TrackerOnsets, M) %#ok<*DEFNU>
ManualOnsets = M.Onset;


nOnsets = length(TrackerOnsets);
Res.Match(1:nOnsets) = NaN;
for i = 1:nOnsets
    t = TrackerOnsets(i);
    dt = ManualOnsets-t;
    [val, id] = min(abs(dt));
    if val < 10
        Res.Match(i) = id;
    end
end

Res.Amplitude(1:nOnsets) = NaN;
Res.Duration(1:nOnsets) = NaN;
for i = 1:length(Res.Match)
    if ~isnan(Res.Match(i))
        Res.Amplitude(i) = M.Amplitude(Res.Match(i));
        Res.Duration(i) = M.Duration(Res.Match(i));
    end
    
end

end




function Stats = getStats(data)
if ~any(~isnan(data))
    Stats = table(NaN,NaN,NaN,NaN,'VariableNames', {'Type','Duration','Amplitude','Onset'});
    return
end

if mean(data, 'omitnan') > 0
    sgn = 1;
elseif mean(data, 'omitnan') < 0
    sgn = -1;
end

[~, peaks] = findpeaks(-sgn*data, 'MinPeakDistance', 3);
[~, troghs] = findpeaks(sgn*data, 'MinPeakDistance', 3);

R = zeros(1, length(data));
R(peaks) = 1;
R(troghs) = 2;

dt = (1/300)*1000;


Stats = [];
for i = 1:length(R)
    if R(i) > 0
        id1 = R(i);
        val = find(R(i+1:end) > 0,1,'first');
        id2 = R(i+val);
        
        if id1 ~= id2
            if id1 == 1
                type = 1; % Retraction
            elseif id1 == 2
                type = 2; % Protraction
            end
            
            duration = val*dt;
            amplitude = data(i+val) - data(i);
            onset = i;
            t = table(type,duration,amplitude,onset,'VariableNames',...
                {'Type','Duration','Amplitude','Onset'});
            
            if isempty(Stats)
                Stats = t;
            else
                Stats = [Stats; t]; %#ok<*AGROW>
            end
        end
    end
    
    
end

end
