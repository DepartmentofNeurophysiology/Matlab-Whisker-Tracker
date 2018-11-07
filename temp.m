filepath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi';
files = dir(fullfile(filepath, '*compiled.mat'));
multi_idx = 0:39;
single_idx = 40:48;

deprived = [];
nn = 0;
h = waitbar(0, 'Loading data...');
for i = 1:size(files, 1)
    load(fullfile(files(i).folder, files(i).name))
    
    Annotations.max_slope = 15;
       if ismember(str2double(files(i).name(6:7)), multi_idx)
                Annotations.deprived = 0;
            elseif ismember(str2double(files(i).name(6:7)), single_idx)
                Annotations.deprived = 1;
        end
        
    A = getAngles(Annotations);
    
    Stats.rmin = getpeaks(A.Tracker.r_min_filtered);
    Stats.rmax = getpeaks(A.Tracker.r_max_filtered);
    Stats.lmin = getpeaks(A.Tracker.l_min_filtered);
    Stats.lmax = getpeaks(A.Tracker.l_max_filtered);
    
    if i == 1
        Allcycles = [ Stats.rmin.table; Stats.rmax.table; Stats.lmin.table; Stats.lmax.table];
        if ismember(str2double(files(i).name(6:7)), multi_idx)
                deprived = zeros(size(Allcycles,1),1);
            elseif ismember(str2double(files(i).name(6:7)), single_idx)
                deprived = ones(size(Allcycles,1), 1);
        end
        
        
        
    else
        Allcycles = [Allcycles; Stats.rmin.table; Stats.rmax.table; Stats.lmin.table; Stats.lmax.table];
        di = size(Allcycles, 1) - size(deprived, 1);
         if ismember(str2double(files(i).name(6:7)), multi_idx)
                deprived = [deprived; zeros(di,1)];
            elseif ismember(str2double(files(i).name(6:7)), single_idx)
                deprived = [deprived; ones(di, 1)];
        end
    end
    waitbar(i/size(files,1))
    
%     
%     if deprived(end) == 1 & nn < 11
%         figure()
%         hold on
%         plot(A.Tracker.r_min_filtered, 'r')
%         plot(A.Tracker.r_max_filtered, 'r')
%         plot(A.Tracker.l_min_filtered, 'b')
%         plot(A.Tracker.l_max_filtered, 'b')
%         nn = nn + 1;
%     end
%     
    
end
close(h)
Allcycles.deprived = deprived;

%%











function stats = getpeaks(data)
warning('off')
t = table([],[],[],[],[],'VariableNames',{'t0','t1','a0','a1','type'});

if numel(find(~isnan(data))) < 10
    stats.table = t;
    stats.peaks = [];
    stats.troghs = [];
    return
end


if mean(data, 'omitnan') > 0
    sgn = 1;
elseif mean(data, 'omitnan') < 0
    sgn = -1;
end

[~, stats.peaks] = findpeaks(-sgn*data, 'MinPeakDistance', 10);
[~, stats.troghs] = findpeaks(sgn*data, 'MinPeakDistance', 10);

peakidx = zeros(size(data));
peakidx(stats.peaks) = 1;
peakidx(stats.troghs) = 2;

table_id = 1;

for i = 1:length(peakidx)
    
    if peakidx(i) ~= 0
        t0 = i;
        a0 = data(t0);
        type = peakidx(i);
        
        id = find(peakidx(i+1:end) > 0,1,'first');
        t1 = i+id;
        a1 = data(t1);
        
        if peakidx(t1) ~= type
            t.t0(table_id, 1) = t0;
            t.t1(table_id, 1) =t1;
            t.a0(table_id, 1) = a0;
            t.a1(table_id, 1) = a1;
            t.type(table_id, 1) = type;
            table_id = table_id + 1;
        end
    end
end

stats.table = t;




end


