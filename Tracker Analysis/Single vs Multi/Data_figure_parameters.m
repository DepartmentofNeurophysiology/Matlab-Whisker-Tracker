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
    
    dist = Annotations.Tracker.dist_nose_target;
    Stats.rmin = getpeaks(A.Tracker.r_min_filtered);
    Stats.rmin.table.dist = dist(Stats.rmin.table.t0);
    Stats.rmax = getpeaks(A.Tracker.r_max_filtered);
    Stats.rmax.table.dist = dist(Stats.rmax.table.t0);
    Stats.lmin = getpeaks(A.Tracker.l_min_filtered);
    Stats.lmin.table.dist = dist(Stats.lmin.table.t0);
    Stats.lmax = getpeaks(A.Tracker.l_max_filtered);
    Stats.lmax.table.dist = dist(Stats.lmax.table.t0);
    
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



Allcycles.dt = (Allcycles.t1 - Allcycles.t0).* (1000/333);
tsize = 10;
tbins = 0:tsize:100;
Allcycles.dtbin = discretize(Allcycles.dt, tbins);
for i = 1:size(Allcycles, 1)
    if ~isnan(Allcycles.dtbin(i))
        Allcycles.dtbin(i) = tbins(Allcycles.dtbin(i)) + tsize/2;
    end
end


Allcycles.da = abs(Allcycles.a1) - abs(Allcycles.a0);
Allcycles.speed =abs( Allcycles.da)./(Allcycles.dt/1000);


keep_idx = ones(size(Allcycles, 1),1);
keep_idx( Allcycles.dt >= 100 ) = 0;
keep_idx( Allcycles.type == 2 & Allcycles.da > -5) = 0;
keep_idx( Allcycles.type == 1 & Allcycles.da < 5) = 0;
keep_idx( Allcycles.dist < 50) = 0;
keep_idx( Allcycles.speed < 50 ) = 0;
keep_idx( Allcycles.speed > 1500 ) = 0;

a = Allcycles(find(keep_idx == 1),:);

Data.All = Allcycles;
Data.All_filtered = a;
%%

I.deprived_pro = find(a.type == 1 & a.deprived == 1);
I.deprived_ret = find(a.type == 2 & a.deprived == 1);
I.full_pro = find(a.type == 1 & a.deprived == 0);
I.full_ret = find(a.type == 2 & a.deprived == 0);



tbins = unique(a.dt);
nbins = length(tbins);

Tbin.deprived_pro = zeros(nbins, 1);
Tbin.deprived_ret = zeros(nbins, 1);
Tbin.full_pro = zeros(nbins, 1);
Tbin.full_ret = zeros(nbins, 1);


for i = 1:nbins
    Tbin.deprived_pro(i) = numel(find(a.dt(I.deprived_pro) == tbins(i)));% & a.dt(I.deprived_pro) <= tbins(i+1)));
    Tbin.deprived_ret(i) = numel(find(a.dt(I.deprived_ret) ==  tbins(i)));% & a.dt(I.deprived_ret) <= tbins(i+1)));
    Tbin.full_pro(i) =     numel(find(a.dt(I.full_pro) == tbins(i)));%& a.dt(I.full_pro) <= tbins(i+1)));
    Tbin.full_ret(i) =     numel(find(a.dt(I.full_ret) == tbins(i)));% & a.dt(I.full_ret) <= tbins(i+1)));
end

Tbin.deprived_pro = Tbin.deprived_pro./numel(I.deprived_pro);
Tbin.deprived_ret = Tbin.deprived_ret./numel(I.deprived_ret);
Tbin.full_pro = Tbin.full_pro./numel(I.full_pro);
Tbin.full_ret = Tbin.full_ret./numel(I.full_ret);

xax = tbins;




Data.Index = I;
Data.Tbin = Tbin;
Data.Tbin_ax = xax;


%%

abins = 0:5:90;
nbins = length(abins) - 1;


Abin.deprived_pro = zeros(nbins, 1);
Abin.deprived_ret = zeros(nbins, 1);
Abin.full_pro = zeros(nbins, 1);
Abin.full_ret = zeros(nbins, 1);

a.da = abs(a.da);
for i = 1:nbins
    Abin.deprived_pro(i) = numel(find(a.da(I.deprived_pro) > abins(i) & a.da(I.deprived_pro) <= abins(i+1)));
    Abin.deprived_ret(i) = numel(find(a.da(I.deprived_ret) > abins(i) & a.da(I.deprived_ret) <= abins(i+1)));
    Abin.full_pro(i) = numel(find(a.da(I.full_pro) > abins(i) & a.da(I.full_pro) <= abins(i+1)));
    Abin.full_ret(i) = numel(find(a.da(I.full_ret) > abins(i) & a.da(I.full_ret) <= abins(i+1)));
end


Abin.deprived_pro = Abin.deprived_pro./numel(I.deprived_pro);
Abin.deprived_ret = Abin.deprived_ret./numel(I.deprived_ret);
Abin.full_pro = Abin.full_pro./numel(I.full_pro);
Abin.full_ret = Abin.full_ret./numel(I.full_ret);

xax = abins(1:end-1) + 2.5;


Data.Abin = Abin;
Data.Abin_ax = xax;

%%

sbins = 0:50:1500;
nbins = length(sbins) - 1;

Sbin.deprived_pro = zeros(nbins, 1);
Sbin.deprived_ret = zeros(nbins, 1);
Sbin.full_pro = zeros(nbins, 1);
Sbin.full_ret = zeros(nbins, 1);

for i = 1:nbins
    Sbin.deprived_pro(i) = numel(find( a.speed(I.deprived_pro) > sbins(i) & a.speed(I.deprived_pro) <= sbins(i+1)));
    Sbin.deprived_ret(i) = numel(find( a.speed(I.deprived_ret) > sbins(i) & a.speed(I.deprived_ret) <= sbins(i+1)));
    Sbin.full_pro(i) =numel(find( a.speed(I.full_pro) > sbins(i) & a.speed(I.full_pro) <= sbins(i+1)));
    Sbin.full_ret(i) = numel(find( a.speed(I.full_ret) > sbins(i) & a.speed(I.full_ret) <= sbins(i+1)));
end


Sbin.deprived_pro = Sbin.deprived_pro./numel(I.deprived_pro);
Sbin.deprived_ret = Sbin.deprived_ret./numel(I.deprived_ret);
Sbin.full_pro = Sbin.full_pro./numel(I.full_pro);
Sbin.full_ret = Sbin.full_ret./numel(I.full_ret);

xax = sbins(1:end-1) + 25;

Data.Sbin = Sbin;
Data.Sbin_ax = xax;


save(fullfile(filepath, 'DataFigPar'), 'Data')
    




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


