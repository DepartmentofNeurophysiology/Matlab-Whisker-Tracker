function Angles = getAngles(Annotations)
%%

kernel_width = 2;
kernel_size = 5;
kernel = normpdf(-kernel_size:kernel_size,0,kernel_width);
filter = kernel./ sum(kernel);




if isfield(Annotations,'Tracker')
    P = Annotations.Tracker.Parameters_clean;
    N = Annotations.Tracker.Nose;
    H = Annotations.Tracker.Headvec;
    T = Annotations.Tracker.Traces_clean;
    
    nframes = max(size(T));
    
    Angles.Tracker.r_min(1:nframes) = NaN;
    Angles.Tracker.r_max(1:nframes) = NaN;
    Angles.Tracker.l_min(1:nframes) = NaN;
    Angles.Tracker.l_max(1:nframes) = NaN;
    
    for i = 1:nframes
        if isempty(P{i})
            continue
        end
        
        
        n = N(i,:);
        h = H(i,:);
        if h(1) == 0
            h(1) = 0.0000001;
        elseif h(2) == 0
            h(2) = 0.0000001;
        end
        y = @(x)  x.*h(1)/h(2) + (n(1)-n(2).*h(1)/h(2));
        
        if strcmp(Annotations.Tracker.Direction, 'Down')
            if h(1)/h(2) > 0
                sgn = 1;
            else
                sgn = -1;
            end
        elseif strcmp(Annotations.Tracker.Direction, 'Up')
            if h(1)/h(2) > 0
                sgn = -1;
            else
                sgn = 1;
            end
        end
        
        l = sgn.*(y(P{i}(:,2)) - P{i}(:,1));
        Angles.lout{i} = l;
        
        
        l_angles = P{i}(l>=0,6);
        l_weights = P{i}(l>=0,7);
        
        [~, l_sort] = sort(l_angles);
        n_left = length(l_angles);
        if mod(n_left,2) == 0
            l_max_idx = l_sort(1:n_left/2);
            l_min_idx = l_sort(n_left/2+1:end);
        elseif mod(n_left,2) == 1
            l_max_idx = l_sort(1:ceil(n_left/2));
            l_min_idx = l_sort(ceil(n_left/2)+1:end);
        end
        
        if isfield(Annotations, 'deprived') && Annotations.deprived == 1
            l_max_idx = [l_min_idx; l_max_idx];
            l_min_idx = [];
        end
        
        Angles.Tracker.l_min(i) = sum( l_angles(l_min_idx).*l_weights(l_min_idx))/ sum(l_weights(l_min_idx));
        Angles.Tracker.l_max(i) = sum( l_angles(l_max_idx).*l_weights(l_max_idx))/ sum(l_weights(l_max_idx));
        
        r_angles = P{i}(l<0,6);
        r_weights = P{i}(l<0,7);
        
        [~, r_sort] = sort(r_angles);
        n_right = length(r_angles);
        if mod(n_right, 2) == 0
            r_min_idx = r_sort(1:n_right/2);
            r_max_idx = r_sort(n_right/2+1:end);
        elseif mod(n_right, 2) == 1
            r_min_idx = r_sort(1:ceil(n_right/2));
            r_max_idx = r_sort(ceil(n_right/2)+1:end);
        end
        
        if isfield(Annotations, 'deprived') && Annotations.deprived == 1
            r_max_idx = [ r_min_idx; r_max_idx];
            r_min_idx = [];
        end
        
        Angles.Tracker.r_min(i) = sum( r_angles(r_min_idx).*r_weights(r_min_idx))/ sum(r_weights(r_min_idx));
        Angles.Tracker.r_max(i) = sum( r_angles(r_max_idx).*r_weights(r_max_idx))/ sum(r_weights(r_max_idx));
        
        
    end
    
    Angles.Tracker.r_min = fillGaps(Angles.Tracker.r_min);
    Angles.Tracker.r_max = fillGaps(Angles.Tracker.r_max);
    Angles.Tracker.l_min = fillGaps(Angles.Tracker.l_min);
    Angles.Tracker.l_max = fillGaps(Angles.Tracker.l_max);
    
    
    switch(Annotations.Output.Direction)
        case 'Down'
            Angles.Tracker.r_min_filtered = conv(Angles.Tracker.r_min, filter, 'same');
            %[Angles.Tracker.r_min_peaks, Angles.Tracker.r_min_troghs] = addPeaks(-Angles.Tracker.r_min_filtered);
            Angles.Tracker.r_max_filtered = conv(Angles.Tracker.r_max, filter, 'same');
            %[Angles.Tracker.r_max_peaks, Angles.Tracker.r_max_troghs] = addPeaks(-Angles.Tracker.r_max_filtered);
            Angles.Tracker.l_min_filtered = conv(Angles.Tracker.l_min, filter, 'same');
            %[Angles.Tracker.l_min_peaks, Angles.Tracker.l_min_troghs] = addPeaks(Angles.Tracker.l_min_filtered);
            Angles.Tracker.l_max_filtered = conv(Angles.Tracker.l_max, filter, 'same');
            %[Angles.Tracker.l_max_peaks, Angles.Tracker.l_max_troghs] = addPeaks(Angles.Tracker.l_max_filtered);
        case 'Up'
            Angles.Tracker.r_min_filtered = conv(Angles.Tracker.r_min, filter, 'same');
            %[Angles.Tracker.r_min_peaks, Angles.Tracker.r_min_troghs] = addPeaks(Angles.Tracker.r_min_filtered);
            Angles.Tracker.r_max_filtered = conv(Angles.Tracker.r_max, filter, 'same');
            %[Angles.Tracker.r_max_peaks, Angles.Tracker.r_max_troghs] = addPeaks(Angles.Tracker.r_max_filtered);
            Angles.Tracker.l_min_filtered = conv(Angles.Tracker.l_min, filter, 'same');
            %[Angles.Tracker.l_min_peaks, Angles.Tracker.l_min_troghs] = addPeaks(-Angles.Tracker.l_min_filtered);
            Angles.Tracker.l_max_filtered = conv(Angles.Tracker.l_max, filter, 'same');
            %[Angles.Tracker.l_max_peaks, Angles.Tracker.l_max_troghs] = addPeaks(-Angles.Tracker.l_max_filtered);
    end
    
    
    if isfield(Annotations, 'max_slope')
        tr = Annotations.max_slope;
        idx_r_min = find(abs(diff( Angles.Tracker.r_min_filtered)) >= tr);
        idx_r_max = find(abs(diff( Angles.Tracker.r_max_filtered)) >= tr);
        idx_l_min = find(abs(diff( Angles.Tracker.l_min_filtered)) >= tr);
        idx_l_max = find(abs(diff( Angles.Tracker.l_max_filtered)) >= tr);
        
        if ~isempty(idx_r_min)
            idx_r_min = [idx_r_min idx_r_min(end)+1];
        end
        if ~isempty(idx_r_max)
            idx_r_max = [idx_r_max idx_r_max(end)+1];
        end
        if ~isempty(idx_l_min)
            idx_l_min = [idx_l_min idx_l_min(end)+1];
        end
        if ~isempty(idx_l_max)
            idx_l_max = [idx_l_max idx_l_max(end)+1];
        end
        
        Angles.Tracker.r_min_filtered( idx_r_min) = NaN;
        Angles.Tracker.r_max_filtered( idx_r_max) = NaN;
        Angles.Tracker.l_min_filtered( idx_l_min) = NaN;
        Angles.Tracker.l_max_filtered( idx_l_max) = NaN;
    end
    
    
end



%%



if isfield(Annotations, 'Manual')
    L = Annotations.Manual.Labels;
    P = Annotations.Manual.Parameters;
    
    nframes = size(L,2);
    
    Manual.r_min(1:nframes) = NaN;
    Manual.r_max(1:nframes) = NaN;
    Manual.l_min(1:nframes) = NaN;
    Manual.l_max(1:nframes) = NaN;
    
    for i = 1:nframes
        if isempty(L{i})
            continue
        end
        lidx = [];
        l = L{i};
        for j = 1:length(l)
            if ~isempty(l{j}) && l{j}(1) == 'L'
                lidx(end+1) = j;
            end
        end
        
        
        if ~isempty(lidx)
            l_angles = P{i}(lidx, 6);
            Manual.l_min(i) = max(l_angles);
            Manual.l_max(i) = min(l_angles);
        end
        
        
        ridx = [];
        for j = 1:length(l)
            if ~isempty(l{j}) && l{j}(1) == 'R'
                ridx(end+1) = j;
            end
        end
        
        if ~isempty(ridx)
            r_angles = P{i}(ridx,  6);
            Manual.r_min(i) = min(r_angles);
            Manual.r_max(i) = max(r_angles);
        end
        
    end
    %Angles.Manual.r_min = fillGaps(Angles.Manual.r_min);
    %Angles.Manual.r_max = fillGaps(Angles.Manual.r_max);
    %Angles.Manual.l_min = fillGaps(Angles.Manual.l_min);
    %Angles.Manual.l_max = fillGaps(Angles.Manual.l_max);
    switch(Annotations.Output.Direction)
        case 'Up'
            Manual.r_min_filtered = conv(Manual.r_min, filter, 'same');
            [Manual.r_min_peaks, Manual.r_min_troghs] = addPeaks(-Manual.r_min_filtered);
            Manual.r_max_filtered = conv(Manual.r_max, filter, 'same');
            [Manual.r_max_peaks, Manual.r_max_troghs] = addPeaks(-Manual.r_max_filtered);
            Manual.l_min_filtered = conv(Manual.l_min, filter, 'same');
            [Manual.l_min_peaks, Manual.l_min_troghs] = addPeaks(Manual.l_min_filtered);
            Manual.l_max_filtered = conv(Manual.l_max, filter, 'same');
            [Manual.l_max_peaks, Manual.l_max_troghs] = addPeaks(Manual.l_max_filtered);
        case 'Down'
            Manual.r_min_filtered = conv(Manual.r_min, filter, 'same');
            [Manual.r_min_peaks, Manual.r_min_troghs] = addPeaks(Manual.r_min_filtered);
            Manual.r_max_filtered = conv(Manual.r_max, filter, 'same');
            [Manual.r_max_peaks, Manual.r_max_troghs] = addPeaks(Manual.r_max_filtered);
            Manual.l_min_filtered = conv(Manual.l_min, filter, 'same');
            [Manual.l_min_peaks, Manual.l_min_troghs] = addPeaks(-Manual.l_min_filtered);
            Manual.l_max_filtered = conv(Manual.l_max, filter, 'same');
            [Manual.l_max_peaks, Manual.l_max_troghs] = addPeaks(-Manual.l_max_filtered);
    end
    
    
    
    
    
    %%
    
    names = {'r_min','l_min','r_max','l_max'};
    
    for j = 1:4
        
        eval(sprintf('data = Angles.Tracker.%s_filtered;', names{j}))
        
        for i = 1:4
            eval(sprintf('sameidx = 1:length(Manual.%s_filtered);',names{i}))
            eval(sprintf('idx = find(~isnan(data(sameidx)) & ~isnan(Manual.%s_filtered(sameidx)));', names{i}))
            eval(sprintf('d(i) = sum(abs(data(sameidx(idx)) - Manual.%s_filtered(sameidx(idx))));', names{i}))
        end
        
        [~, id] = min(d);
        
        eval(sprintf('Angles.Manual.%s = Manual.%s;',names{j}, names{id}))
        eval(sprintf('Angles.Manual.%s_filtered = Manual.%s_filtered;',names{j}, names{id}))
        eval(sprintf('Angles.Manual.%s_peaks = Manual.%s_peaks;',names{j}, names{id}))
        eval(sprintf('Angles.Manual.%s_troghs = Manual.%s_troghs;',names{j}, names{id}))
        
        
    end
    
    
end

end


function [peaks,troghs] = addPeaks(data)
[~, peaks] = findpeaks(data, 'MinPeakDistance', 10);
[~, troghs] = findpeaks(-data, 'MinPeakDistance', 10);
end

function data = fillGaps(data)
%%

start = find(~isnan(data), 1, 'first');
last= find(~isnan(data), 1,'last');


keep_nan = zeros(1, length(data));

for i = start+1:last-1
    if isnan(data(i))
        d_next = data(i+1:end);
        idx = find(~isnan(d_next), 1, 'first');
        
        if idx > 15
            keep_nan(i:i+idx) = 1;
        end
        
        if i - start > 2
            startval = mean(data(i-2:i-1));
        else
            startval = data(i-1);
        end
        
        
        if ~isnan(data(i+idx+1))
            endval = mean(data(i+idx:i+idx+1));
        else
            endval = data(i+idx);
        end
        
        n_to_fill = idx;
        for j = 1:idx
            
            if keep_nan(i+j-1)
                continue
            end
            fillval = startval + j*(endval-startval)/n_to_fill;
            data(i+j-1) = fillval;
        end
        
        
        
    end
end

end

