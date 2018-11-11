
%% Data axes 1 - EXAMPLE FRAME

example_file = 1;
example_frame = 1326;
Datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
comp_files = dir(fullfile(Datapath, '*compiled.mat'));

%%
load(fullfile(Datapath, 'Frames_with_nose.mat'))
load(fullfile(comp_files(example_file).folder, comp_files(example_file).name))

Settings = Annotations.Settings;
Settings.Current_frame = example_frame;
Frame = LoadFrame(Settings);
Frame = Frame-min(min(Frame));
Frame = Frame./max(max(Frame));
Data.Frame = Frame;

P = Annotations.Tracker.Parameters_clean;
N = Annotations.Tracker.Nose;
H = Annotations.Tracker.Headvec;
T = Annotations.Tracker.Traces_clean;

nframes = size(T,1);

r_min(1:nframes) = NaN;
r_max(1:nframes) = NaN;
l_min(1:nframes) = NaN;
l_max(1:nframes) = NaN;

i = example_frame;
n = N(i,:);
h = H(i,:);
y = @(x)  x.*h(1)/h(2) + (n(1)-n(2).*h(1)/h(2));

if h(1)/h(2) > 0
    sgn = 1;
else
    sgn = -1;
end

l = sgn.*(y(P{i}(:,2)) - P{i}(:,1));
l_idx = find(l<=0);
r_idx = find(l>0);
l_angles = P{i}(l<=0,6);
l_weights = P{i}(l<=0,7);

[~, l_sort] = sort(l_angles);
n_left = length(l_angles);
if mod(n_left,2) == 0
    l_max_idx = l_sort(1:n_left/2);
    l_min_idx = l_sort(n_left/2+1:end);
elseif mod(n_left,2) == 1
    l_max_idx = l_sort(1:ceil(n_left/2));
    l_min_idx = l_sort(ceil(n_left/2)+1:end);
end


r_angles = P{i}(l>0,6);
r_weights = P{i}(l>0,7);

[~, r_sort] = sort(r_angles);
n_right = length(r_angles);
if mod(n_right, 2) == 0
    r_min_idx = r_sort(1:n_right/2);
    r_max_idx = r_sort(n_right/2+1:end);
elseif mod(n_right, 2) == 1
    r_min_idx = r_sort(1:ceil(n_right/2));
    r_max_idx = r_sort(ceil(n_right/2)+1:end);
end

T = T{example_frame};
Data.TCRflag = zeros(1, size(T,2));
for i = 1:size(T,2)
    if ismember(i, l_idx(l_min_idx)) | ismember(i, r_idx(r_min_idx))
        Data.TCRflag(i) = 1;
    end
end
Data.TTraces = T;

Data.MTraces = Annotations.Manual.Traces{example_frame};
L = Annotations.Manual.Labels;
P = Annotations.Manual.Parameters;
i = example_frame;

lidx = [];
l = L{i};
for j = 1:length(l)
    if ~isempty(l{j}) && l{j}(1) == 'L'
        lidx(end+1) = j;
    end
end

if ~isempty(lidx)
    l_angles = P{i}(lidx, 6);
    [~,idma] = max(l_angles);
    [~,idmi] = min(l_angles);
end
Data.MCRflag(lidx(idma)) = 1;
Data.MCRflag(lidx(idmi)) = 2;

ridx = [];
for j = 1:length(l)
    if ~isempty(l{j}) && l{j}(1) == 'R'
        ridx(end+1) = j;
    end
end

if ~isempty(ridx)
    r_angles = P{i}(ridx,  6);
    [~,idma] = min(r_angles);
    [~,idmi]  = max(r_angles);
end
Data.MCRflag(ridx(idma)) = 1;
Data.MCRflag(ridx(idmi)) = 2;


%% AX2 - Thetas
A = getAngles(Annotations);
IDX = find(~isnan(A.Tracker.r_min_filtered) & ~isnan(A.Tracker.r_max_filtered) & ...
    ~isnan(A.Tracker.l_min_filtered) & ~isnan(A.Tracker.l_max_filtered) & ...
    ~isnan(A.Manual.r_min_filtered) & ~isnan(A.Manual.r_max_filtered) & ...
    ~isnan(A.Manual.l_min_filtered) & ~isnan(A.Manual.r_max_filtered));
nframes = size(Annotations.Output.Traces,1);
frange = zeros(1, nframes);
pairs = Output(example_file).Pairs;
for j = 1:size(pairs, 2)
    frange(pairs{j}(1):pairs{j}(2)) = 1;
end

Data.IDX = IDX;
Data.Angles = A;

dt = 1000/300;
Data.xax = [1:size(A.Tracker.r_min,2)]*dt;
% Data.corr.l_min = compareThetas(Data.Angles.Tracker.l_min_filtered, Data.Angles.Tracker.l_min_peaks, ...
%     Data.Angles.Manual.l_min_filtered, Data.Angles.Manual.l_min_peaks, frange);
% Data.corr.l_max = compareThetas(Data.Angles.Tracker.l_max_filtered, Data.Angles.Tracker.l_max_peaks, ...
%     Data.Angles.Manual.l_max_filtered, Data.Angles.Manual.l_max_peaks, frange);
% Data.corr.r_min = compareThetas(Data.Angles.Tracker.r_min_filtered, Data.Angles.Tracker.r_min_peaks, ...
%     Data.Angles.Manual.r_min_filtered, Data.Angles.Manual.r_min_peaks,frange);
% Data.corr.r_max = compareThetas(Data.Angles.Tracker.r_max_filtered, Data.Angles.Tracker.r_max_peaks, ...
%     Data.Angles.Manual.r_max_filtered, Data.Angles.Manual.r_max_peaks,frange);

Data.corr.l_min = compareThetas(Data.Angles.Tracker.l_min_filtered, [], ...
    Data.Angles.Manual.l_min_filtered, [], frange);
Data.corr.l_max = compareThetas(Data.Angles.Tracker.l_max_filtered, [], ...
    Data.Angles.Manual.l_max_filtered,[], frange);
Data.corr.r_min = compareThetas(Data.Angles.Tracker.r_min_filtered, [], ...
    Data.Angles.Manual.r_min_filtered,[],frange);
Data.corr.r_max = compareThetas(Data.Angles.Tracker.r_max_filtered, [], ...
    Data.Angles.Manual.r_max_filtered,[],frange);
%% AX3 - Correlation
h = waitbar(0, 'doing things');
for i = 1:size(comp_files,1)
    load(fullfile(comp_files(i).folder, comp_files(i).name))
    nframes = size(Annotations.Output.Traces,1);
    frange = zeros(1, nframes);
    pairs = Output(i).Pairs;
    for j = 1:size(pairs, 2)
        frange(pairs{j}(1):pairs{j}(2)) = 1;
    end
    
    if numel(find(frange)) < 25
        Out(i).l_min = NaN;
        Out(i).l_max = NaN;
        Out(i).r_min = NaN;
        Out(i).r_max = NaN;
        continue
    end
    
    Angles = getAngles(Annotations);
    
    
%     
%     
%     Out(i).l_min = compareThetas(Angles.Tracker.l_min_filtered, Angles.Tracker.l_min_peaks, ...
%         Angles.Manual.l_min_filtered, Angles.Manual.l_min_peaks, frange);
%     Out(i).l_max = compareThetas(Angles.Tracker.l_max_filtered, Angles.Tracker.l_max_peaks, ...
%         Angles.Manual.l_max_filtered, Angles.Manual.l_max_peaks, frange);
%     Out(i).r_min = compareThetas(Angles.Tracker.r_min_filtered, Angles.Tracker.r_min_peaks, ...
%         Angles.Manual.r_min_filtered, Angles.Manual.r_min_peaks,frange);
%     Out(i).r_max = compareThetas(Angles.Tracker.r_max_filtered, Angles.Tracker.r_max_peaks, ...
%         Angles.Manual.r_max_filtered, Angles.Manual.r_max_peaks,frange);
%        
    Out(i).l_min = compareThetas(Angles.Tracker.l_min_filtered, [], ...
        Angles.Manual.l_min_filtered, [], frange);
    Out(i).l_max = compareThetas(Angles.Tracker.l_max_filtered, [], ...
        Angles.Manual.l_max_filtered, [], frange);
    Out(i).r_min = compareThetas(Angles.Tracker.r_min_filtered, [], ...
        Angles.Manual.r_min_filtered, [],frange);
    Out(i).r_max = compareThetas(Angles.Tracker.r_max_filtered,[], ...
        Angles.Manual.r_max_filtered, [],frange);
    
    waitbar(i/size(comp_files,1))
end

close(h)
Data.Correlation.Rlmin(1:size(Out, 2)) = NaN;
Data.Correlation.Rlmax(1:size(Out, 2)) = NaN;
Data.Correlation.Rrmin(1:size(Out, 2)) = NaN;
Data.Correlation.Rrmax(1:size(Out, 2)) = NaN;

for i =1:size(Out, 2)
    if isstruct(Out(i).l_min)
        Data.Correlation.Rlmin(i) = Out(i).l_min.R;
        Data.Correlation.Rlmax(i) = Out(i).l_max.R;
        Data.Correlation.Rrmin(i) = Out(i).r_min.R;
        Data.Correlation.Rrmax(i) = Out(i).r_max.R;
    end
end


Data.Correlation.Total = [...
    Data.Correlation.Rlmin,...
    Data.Correlation.Rlmax,...
    Data.Correlation.Rrmin,...
    Data.Correlation.Rrmax];


%% AX4 - theta inset
snaptime = 4410;

data = Data.Angles.Tracker.r_max_filtered;
if mean(data, 'omitnan') > 0
    sgn = 1;
elseif mean(data, 'omitnan') < 0
    sgn = -1;
end
[~, peaks] = findpeaks(-sgn*data, 'MinPeakDistance',3);
[~, troghs] = findpeaks(sgn*data, 'MinPeakDistance', 3);
I_ret = peaks*dt - snaptime;
id = find(I_ret > 0,1,'first');
Data.TPeak = peaks(id);


I_ret = troghs*dt - snaptime;
id = find(I_ret > 0,1,'first');
Data.Ttrogh= troghs(id);
Data.Ttrogh2= troghs(id+1);

%%
data = Data.Angles.Manual.r_max_filtered;
if mean(data, 'omitnan') > 0
    sgn = 1;
elseif mean(data, 'omitnan') < 0
    sgn = -1;
end
[pval, peaks] = findpeaks(-sgn*data, 'MinPeakDistance',3);
[p, troghs] = findpeaks(sgn*data, 'MinPeakDistance', 3);
I_ret = peaks*dt - snaptime;
id = find(I_ret > 0,1,'first');
Data.MPeak = peaks(id);

I_ret = troghs*dt - snaptime;
id = find(I_ret > 0,1,'first');
Data.Mtrogh= troghs(id);
Data.Mtrogh2 = troghs(id+1);



%% AX5-8 measured parameters

h = waitbar(0,'gathering data');

for i = 1:size(comp_files, 1)
    
     C = [Data.Correlation.Rlmin(i), Data.Correlation.Rlmax(i), Data.Correlation.Rrmin(i), Data.Correlation.Rrmax(i)];
%     if numel(find(C < 0.5)) > 1
%         disp(i)
%         continue
%     end
%     
    
    load(fullfile(comp_files(i).folder, comp_files(i).name))
    
    p = getWhiskingStats(Annotations);
    
    if i == 1
        Pars = p;
    else
%         if istable(p.Tr_min)
%             Pars.Tr_min = [Pars.Tr_min; p.Tr_min];
%         end
%         if istable(p.Tr_max)
%             Pars.Tr_max = [Pars.Tr_max; p.Tr_max];
%         end
%         if istable(p.Tl_min)
%             Pars.Tl_min = [Pars.Tl_min; p.Tl_min];
%         end
%         if istable(p.Tl_max)
%             Pars.Tl_max = [Pars.Tl_max; p.Tl_max];
%         end
%         
%         if ~isempty(p.Mr_min)
%             Pars.Mr_min = [Pars.Mr_min; p.Mr_min];
%         end
%         if ~isempty(p.Mr_max)
%             Pars.Mr_max = [Pars.Mr_max; p.Mr_max];
%         end
%         if~isempty(p.Ml_min)
%             Pars.Ml_min = [Pars.Ml_min; p.Ml_min];
%         end
%         if~isempty(p.Ml_max)
%             Pars.Ml_max = [Pars.Ml_max; p.Ml_max];
%         end
%         
%         
        
        
        if istable(p.Tr_min) & Data.Correlation.Rrmin(i) > 0.5
            Pars.Tr_min = [Pars.Tr_min; p.Tr_min];
        end
        if istable(p.Tr_max) & Data.Correlation.Rrmax(i) > 0.5
            Pars.Tr_max = [Pars.Tr_max; p.Tr_max];
        end
        if istable(p.Tl_min) & Data.Correlation.Rlmin(i) > 0.5
            Pars.Tl_min = [Pars.Tl_min; p.Tl_min];
        end
        if istable(p.Tl_max) & Data.Correlation.Rlmax(i) > 0.5
            Pars.Tl_max = [Pars.Tl_max; p.Tl_max];
        end
        
        if ~isempty(p.Mr_min)  & Data.Correlation.Rrmin(i) > 0.5
            Pars.Mr_min = [Pars.Mr_min; p.Mr_min];
        end
        if ~isempty(p.Mr_max)  & Data.Correlation.Rrmax(i) > 0.5
            Pars.Mr_max = [Pars.Mr_max; p.Mr_max];
        end
        if~isempty(p.Ml_min) & Data.Correlation.Rlmin(i) > 0.5
            Pars.Ml_min = [Pars.Ml_min; p.Ml_min];
        end
        if~isempty(p.Ml_max)  & Data.Correlation.Rlmax(i) > 0.5
            Pars.Ml_max = [Pars.Ml_max; p.Ml_max];
        end
    end
   
   waitbar(i/size(comp_files,1))
end

close(h)


names = {'r_min','r_max','l_min','l_max'};
minAmplitude = 0;
maxAmplitude = 60;
maxDuration = 60;
for i = 1:4
    eval(sprintf('I_ret = find(~isnan(Pars.T%s.Match)  & Pars.T%s.Type == 1 & Pars.T%s.Amplitude > minAmplitude & Pars.T%s.Amplitude < maxAmplitude & Pars.T%s.MAmplitude > minAmplitude & Pars.T%s.MAmplitude < maxAmplitude & Pars.T%s.Duration < maxDuration & Pars.T%s.MDuration < maxDuration);', names{i},names{i},names{i}, names{i}, names{i}, names{i},names{i},names{i}))
    eval(sprintf('I_pro = find(~isnan(Pars.T%s.Match)  & Pars.T%s.Type == 2 & Pars.T%s.Amplitude < -minAmplitude & Pars.T%s.Amplitude > -maxAmplitude & Pars.T%s.MAmplitude < -minAmplitude & Pars.T%s.MAmplitude > -maxAmplitude & Pars.T%s.Duration < maxDuration & Pars.T%s.MDuration < maxDuration);', names{i},names{i},names{i}, names{i}, names{i}, names{i},names{i},names{i}))

    eval(sprintf('l_ret = table(Pars.T%s.Amplitude(I_ret), Pars.T%s.MAmplitude(I_ret), Pars.T%s.Duration(I_ret), Pars.T%s.MDuration(I_ret));',names{i},names{i}, names{i}, names{i}))
    eval(sprintf('l_pro = table(Pars.T%s.Amplitude(I_pro), Pars.T%s.MAmplitude(I_pro), Pars.T%s.Duration(I_pro), Pars.T%s.MDuration(I_pro));',names{i},names{i}, names{i}, names{i}))

    if i == 1
        r_Pro = l_pro;
        r_Ret = l_ret;
    else
        r_Pro = [r_Pro; l_pro];
        r_Ret = [r_Ret; l_ret];
    end

end

Data.Protraction_table = r_Pro;
Data.Protraction_table.Properties.VariableNames = {'Tamplitude','Mamplitude','Tduration','Mduration'};
Data.Retraction_table = r_Ret;
Data.Retraction_table.Properties.VariableNames = {'Tamplitude','Mamplitude','Tduration','Mduration'};


Amplitude_bins = 0:5:maxAmplitude;

t = Data.Protraction_table;
g = table(discretize(abs(t.Tamplitude),Amplitude_bins), discretize(abs(t.Mamplitude), Amplitude_bins),...
    'VariableNames',{'Tamplitude','Mamplitude'});
a = varfun(@numel, g,'GroupingVariables'...
    ,{'Tamplitude','Mamplitude'});

I = zeros(length(Amplitude_bins)-1, length(Amplitude_bins)-1);
for i = 1:size(a, 1)
    I(a.Mamplitude(i), a.Tamplitude(i)) = a.GroupCount(i);
end
Data.histImProA = I./max(max(I));
b = polyfit(abs(t.Mamplitude), abs(t.Tamplitude),1);
yfit = polyval(b,abs(t.Mamplitude));
SSE = sum( (abs(t.Tamplitude) - mean(abs(t.Tamplitude))).^2 );
SSTO = sum( (abs(t.Tamplitude) - yfit).^2);
Data.histImProA_Rs = abs(1 - SSE/SSTO);







t = Data.Retraction_table;
g = table(discretize(abs(t.Tamplitude),Amplitude_bins), discretize(abs(t.Mamplitude),Amplitude_bins),...
    'VariableNames',{'Tamplitude','Mamplitude'});
a = varfun(@numel, g,'GroupingVariables'...
    ,{'Tamplitude','Mamplitude'});

I = zeros(length(Amplitude_bins)-1, length(Amplitude_bins)-1);
for i = 1:size(a, 1)
    I(a.Mamplitude(i), a.Tamplitude(i)) = a.GroupCount(i);
end
Data.histImRetA = I./max(max(I));
b = polyfit(abs(t.Mamplitude), abs(t.Tamplitude),1);
yfit = polyval(b,abs(t.Mamplitude));
SSE = sum( (abs(t.Tamplitude) - mean(abs(t.Tamplitude))).^2 );
SSTO = sum( (abs(t.Tamplitude) - yfit).^2);
Data.histImRetA_Rs = abs(1 - SSE/SSTO);
Data.Amplitude_bins = Amplitude_bins;



Duration_bins = 0:5:maxDuration;

t = Data.Protraction_table;
g = table(discretize(abs(t.Tduration),Duration_bins), discretize(abs(t.Mduration),Duration_bins),...
    'VariableNames',{'Tduration','Mduration'});
a = varfun(@numel, g,'GroupingVariables'...
    ,{'Tduration','Mduration'});

I = zeros(length(Duration_bins)-1, length(Duration_bins)-1);
for i = 1:size(a, 1)
    I(a.Mduration(i), a.Tduration(i)) = a.GroupCount(i);
end
b = polyfit(t.Mduration, t.Tduration, 1);
yfit = polyval(b, t.Mduration);
SSE = sum( (t.Tduration - mean(t.Tduration)).^2 );
SSTO = sum( (t.Tduration - yfit).^2);
Data.histImProD_Rs = abs(1- SSE/SSTO);
Data.histImProD =  I./max(max(I));

t = Data.Retraction_table;
g = table(discretize(abs(t.Tduration),Duration_bins), discretize(abs(t.Mduration), Duration_bins),...
    'VariableNames',{'Tduration','Mduration'});
a = varfun(@numel, g,'GroupingVariables'...
    ,{'Tduration','Mduration'});

I = zeros(length(Duration_bins)-1, length(Duration_bins)-1);
for i = 1:size(a, 1)
    I(a.Mduration(i), a.Tduration(i)) = a.GroupCount(i);
end
Data.histImRetD =  I./max(max(I));
b = polyfit(t.Mduration, t.Tduration, 1);
yfit = polyval(b, t.Mduration);
SSE = sum( (t.Tduration - mean(t.Tduration)).^2 );
SSTO = sum( (t.Tduration - yfit).^2);
Data.histImRetD_Rs = abs(1- SSE/SSTO);
Data.Duration_bins = Duration_bins;



Speed_bins = 0:100:2000;
TrackerSpeedPro = 1000*abs(Data.Protraction_table.Tamplitude)./ Data.Protraction_table.Tduration;
ManualSpeedPro = 1000*abs(Data.Protraction_table.Mamplitude)./Data.Protraction_table.Mduration;
g = table(discretize(TrackerSpeedPro,Speed_bins), discretize(ManualSpeedPro, Speed_bins),...
    'VariableNames',{'Tspeed','Mspeed'});
Speed_bins = 0:100:2000;
a = varfun(@numel, g,'GroupingVariables'...
    ,{'Tspeed','Mspeed'});

I = zeros(length(Speed_bins)-1, length(Speed_bins)-1);
for i = 1:size(a, 1)
    I(a.Mspeed(i), a.Tspeed(i)) = a.GroupCount(i);
end
b = polyfit(ManualSpeedPro ,TrackerSpeedPro, 1);
yfit = polyval(b, ManualSpeedPro);
SSE = sum( (TrackerSpeedPro - mean(TrackerSpeedPro)).^2);
SSTO = sum( (TrackerSpeedPro - yfit).^2);
Data.histImProS_Rs = abs(1-SSE/SSTO);
Data.histImProS = I./max(max(I));


TrackerSpeedRet = 1000*abs(Data.Retraction_table.Tamplitude)./ Data.Retraction_table.Tduration;
ManualSpeedRet = 1000*abs(Data.Retraction_table.Mamplitude)./Data.Retraction_table.Mduration;
g = table(discretize(TrackerSpeedRet,Speed_bins), discretize(ManualSpeedRet, Speed_bins),...
    'VariableNames',{'Tspeed','Mspeed'});
Speed_bins = 0:100:2000;
a = varfun(@numel, g,'GroupingVariables'...
    ,{'Tspeed','Mspeed'});

I = zeros(length(Speed_bins)-1, length(Speed_bins)-1);
for i = 1:size(a, 1)
    I(a.Mspeed(i), a.Tspeed(i)) = a.GroupCount(i);
end
Data.histImRetS =  I./max(max(I));
b = polyfit(ManualSpeedRet ,TrackerSpeedRet, 1);
yfit = polyval(b, ManualSpeedRet);
SSE = sum( (TrackerSpeedRet - mean(TrackerSpeedRet)).^2);
SSTO = sum( (TrackerSpeedRet - yfit).^2);
Data.histImRetS_Rs = abs(1-SSE/SSTO);
Data.Speed_bins = Speed_bins;

%% AX 11 - Touch Data

h = waitbar(0, 'doing touch');



for i =  1:size(comp_files,1)
     load(fullfile(comp_files(i).folder, comp_files(i).name))
     
     TrackerTouch = Annotations.Tracker.TouchFiltered;
     
     nframes = size(Annotations.Tracker.Touch, 2);
     TTouch = [];
     MTouch = [];

     for j = 1:nframes  
         
         if j <= size(Annotations.Manual.Touch.pt,2) && ...
             ~isempty(Annotations.Manual.RawNotations{j}) && ...
                 ~isempty(Annotations.Tracker.Traces_clean{j})
             
            TTouch(end+1) = numel(find(TrackerTouch{j}));
            MTouch(end+1) = size(Annotations.Manual.Touch.pt{j},1);
            
       
         end
         
         
      
     end
   
     
     if i == 1
         TOUCH = table(TTouch', MTouch','VariableNames',{'Tracker','Manual'});
     else
         looptouch = table(TTouch', MTouch','VariableNames',{'Tracker','Manual'});
         TOUCH  = [TOUCH;looptouch];
     end
    
    waitbar(i/size(comp_files,1))
    
end
close(h)



I_ret = find( TOUCH.Tracker ~= 0 | TOUCH.Manual ~= 0);
TOUCH = TOUCH(I_ret,:);
Data.TOUCH = TOUCH;


%maxval = max([TOUCH.Tracker;TOUCH.Manual]);
maxval = 8;
nbins = maxval+1;

TOUCHCOUNT = zeros(maxval+1,maxval+1);
for i = 0:maxval
    for j = 0:maxval
        TOUCHCOUNT(i+1,j+1) = numel(find(TOUCH.Tracker == i & TOUCH.Manual == j));
    end
end

Data.TOUCHCOUNT = TOUCHCOUNT./sum(sum(TOUCHCOUNT));
xtickax = 1:2:maxval+1;
Data.TOUCHtickax = xtickax;
Data.TOUCHticklabel = 0:2:maxval;


b = polyfit(Data.TOUCH.Manual, Data.TOUCH.Tracker, 1);
yfit = polyval(b, Data.TOUCH.Manual);
SSE = sum( (Data.TOUCH.Tracker - mean(Data.TOUCH.Tracker)).^2);
SSTO = sum( (Data.TOUCH.Tracker - yfit).^2);
Data.touch_Rs = abs(1-SSE/SSTO);

Data.touch_diff = Data.TOUCH.Tracker - Data.TOUCH.Manual;


save(fullfile(Datapath,'Data_Figure_Par_Eval'),'Data')




