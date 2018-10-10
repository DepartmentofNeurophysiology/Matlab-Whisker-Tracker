
%% Data axes 1 - EXAMPLE FRAME

example_file = 1;
example_frame = 1326;
Datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
comp_files = dir(fullfile(Datapath, '*compiled.mat'));
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
Data.corr.l_min = compareThetas(Data.Angles.Tracker.l_min_filtered, Data.Angles.Tracker.l_min_peaks, ...
    Data.Angles.Manual.l_min_filtered, Data.Angles.Manual.l_min_peaks, frange);
Data.corr.l_max = compareThetas(Data.Angles.Tracker.l_max_filtered, Data.Angles.Tracker.l_max_peaks, ...
    Data.Angles.Manual.l_max_filtered, Data.Angles.Manual.l_max_peaks, frange);
Data.corr.r_min = compareThetas(Data.Angles.Tracker.r_min_filtered, Data.Angles.Tracker.r_min_peaks, ...
    Data.Angles.Manual.r_min_filtered, Data.Angles.Manual.r_min_peaks,frange);
Data.corr.r_max = compareThetas(Data.Angles.Tracker.r_max_filtered, Data.Angles.Tracker.r_max_peaks, ...
    Data.Angles.Manual.r_max_filtered, Data.Angles.Manual.r_max_peaks,frange);


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
    
    
    
    
    Out(i).l_min = compareThetas(Angles.Tracker.l_min_filtered, Angles.Tracker.l_min_peaks, ...
        Angles.Manual.l_min_filtered, Angles.Manual.l_min_peaks, frange);
    Out(i).l_max = compareThetas(Angles.Tracker.l_max_filtered, Angles.Tracker.l_max_peaks, ...
        Angles.Manual.l_max_filtered, Angles.Manual.l_max_peaks, frange);
    Out(i).r_min = compareThetas(Angles.Tracker.r_min_filtered, Angles.Tracker.r_min_peaks, ...
        Angles.Manual.r_min_filtered, Angles.Manual.r_min_peaks,frange);
    Out(i).r_max = compareThetas(Angles.Tracker.r_max_filtered, Angles.Tracker.r_max_peaks, ...
        Angles.Manual.r_max_filtered, Angles.Manual.r_max_peaks,frange);
    
    
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

idx = Data.Angles.Tracker.r_max_peaks*dt - snaptime;
id = find(idx > 0,1,'first');
Data.TPeak = Data.Angles.Tracker.r_max_peaks(id);
idx = Data.Angles.Tracker.r_max_troghs*dt - snaptime;
id = find(idx > 0,1,'first');
Data.Ttrogh= Data.Angles.Tracker.r_max_troghs(id);
Data.Ttrogh2= Data.Angles.Tracker.r_max_troghs(id+1);


idx = Data.Angles.Manual.r_max_peaks*dt - snaptime;
id = find(idx > 0,1,'first');
Data.MPeak = Data.Angles.Manual.r_max_peaks(id);

idx = Data.Angles.Manual.r_max_troghs*dt - snaptime;
id = find(idx > 0,1,'first');
Data.Mtrogh= Data.Angles.Manual.r_max_troghs(id);
Data.Mtrogh2 = Data.Angles.Manual.r_max_troghs(id+1);



%% AX5-8 measured parameters


h = waitbar(0,'gathering data');

for i =   1:size(comp_files,1)
    %%
    % Load file
    load(fullfile(comp_files(i).folder, comp_files(i).name))
    if isempty(Annotations)
        disp(i)
        continue
    end
    A = getAngles(Annotations);
    
    if ~isfield(A, 'Manual')
        disp(i)
        continue
    end
    [Mrmi, Trmi] = getWhiskStats(A.Manual.r_min_filtered, A.Manual.r_min_peaks, A.Manual.r_min_troghs,...
        A.Tracker.r_min_filtered, A.Tracker.r_min_peaks, A.Tracker.r_min_troghs);
    [Mlmi, Tlmi] = getWhiskStats(-A.Manual.l_min_filtered, A.Manual.l_min_peaks, A.Manual.l_min_troghs,...
        -A.Tracker.l_min_filtered, A.Tracker.l_min_peaks, A.Tracker.l_min_troghs);
    [Mrma, Trma] = getWhiskStats(A.Manual.r_max_filtered, A.Manual.r_max_peaks, A.Manual.r_max_troghs,...
        A.Tracker.r_max_filtered, A.Tracker.r_max_peaks, A.Tracker.r_max_troghs);
    [Mlma, Tlma] = getWhiskStats(-A.Manual.l_max_filtered, A.Manual.l_max_peaks, A.Manual.l_max_troghs,...
        -A.Tracker.l_max_filtered, A.Tracker.l_max_peaks, A.Tracker.l_max_troghs);
    
    if i == 1
        LeftMinPro = table(...
            Mlmi.protraction_amplitude', Mlmi.protraction_duration', Mlmi.protraction_speed',...
            Tlmi.protraction_amplitude', Tlmi.protraction_duration', Tlmi.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        RightMinPro = table(...
            Mrmi.protraction_amplitude', Mrmi.protraction_duration', Mrmi.protraction_speed',...
            Trmi.protraction_amplitude', Trmi.protraction_duration', Trmi.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        LeftMaxPro = table(...
            Mlma.protraction_amplitude', Mlma.protraction_duration', Mlma.protraction_speed',...
            Tlma.protraction_amplitude', Tlma.protraction_duration', Tlma.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        RightMaxPro = table(...
            Mrma.protraction_amplitude', Mrma.protraction_duration', Mrma.protraction_speed',...
            Trma.protraction_amplitude', Trma.protraction_duration', Trma.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        
        LeftMinCon = table(...
            Mlmi.contraction_amplitude', Mlmi.contraction_duration', Mlmi.contraction_speed',...
            Tlmi.contraction_amplitude', Tlmi.contraction_duration', Tlmi.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        RightMinCon = table(...
            Mrmi.contraction_amplitude', Mrmi.contraction_duration', Mrmi.contraction_speed',...
            Trmi.contraction_amplitude', Trmi.contraction_duration', Trmi.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        LeftMaxCon = table(...
            Mlma.contraction_amplitude', Mlma.contraction_duration', Mlma.contraction_speed',...
            Tlma.contraction_amplitude', Tlma.contraction_duration', Tlma.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        RightMaxCon = table(...
            Mrma.contraction_amplitude', Mrma.contraction_duration', Mrma.contraction_speed',...
            Trma.contraction_amplitude', Trma.contraction_duration', Trma.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        
    else
        LeftMinProNew = table(...
            Mlmi.protraction_amplitude', Mlmi.protraction_duration', Mlmi.protraction_speed',...
            Tlmi.protraction_amplitude', Tlmi.protraction_duration', Tlmi.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        RightMinProNew = table(...
            Mrmi.protraction_amplitude', Mrmi.protraction_duration', Mrmi.protraction_speed',...
            Trmi.protraction_amplitude', Trmi.protraction_duration', Trmi.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        LeftMaxProNew = table(...
            Mlma.protraction_amplitude', Mlma.protraction_duration', Mlma.protraction_speed',...
            Tlma.protraction_amplitude', Tlma.protraction_duration', Tlma.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        RightMaxProNew = table(...
            Mrma.protraction_amplitude', Mrma.protraction_duration', Mrma.protraction_speed',...
            Trma.protraction_amplitude', Trma.protraction_duration', Trma.protraction_speed',...
            'VariableNames',{'Man_pro_a','Man_pro_d','Man_pro_s','Trck_pro_a','Trck_pro_d','Trck_pro_s'});
        
        LeftMinConNew = table(...
            Mlmi.contraction_amplitude', Mlmi.contraction_duration', Mlmi.contraction_speed',...
            Tlmi.contraction_amplitude', Tlmi.contraction_duration', Tlmi.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        RightMinConNew = table(...
            Mrmi.contraction_amplitude', Mrmi.contraction_duration', Mrmi.contraction_speed',...
            Trmi.contraction_amplitude', Trmi.contraction_duration', Trmi.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        LeftMaxConNew = table(...
            Mlma.contraction_amplitude', Mlma.contraction_duration', Mlma.contraction_speed',...
            Tlma.contraction_amplitude', Tlma.contraction_duration', Tlma.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        RightMaxConNew = table(...
            Mrma.contraction_amplitude', Mrma.contraction_duration', Mrma.contraction_speed',...
            Trma.contraction_amplitude', Trma.contraction_duration', Trma.contraction_speed',...
            'VariableNames',{'Man_Con_a','Man_Con_d','Man_Con_s','Trck_Con_a','Trck_Con_d','Trck_Con_s'});
        
        LeftMinPro = [LeftMinPro; LeftMinProNew];
        RightMinPro = [RightMinPro; RightMinProNew];
        LeftMaxPro = [LeftMaxPro; LeftMaxProNew];
        RightMaxPro = [RightMaxPro; RightMaxProNew];
        
        LeftMinCon = [LeftMinCon; LeftMinConNew];
        RightMinCon = [RightMinCon; RightMinConNew];
        LeftMaxCon = [LeftMaxCon; LeftMaxConNew];
        RightMaxCon = [RightMaxCon; RightMaxConNew];
        
        
    end
    
    
    waitbar(i/size(comp_files,1))
    
end
close(h)


%%

% Filter noise data
MAXAMP = 180;
MAXDUR = 200;
MAXSPEED = 2500;

names = {'LeftMinCon','LeftMaxCon','RightMinCon','RightMaxCon'};

for i = 1:length(names)
    eval(...
    sprintf(['idx = find('...
        '%s.Man_Con_a < MAXAMP & ' ...
        '%s.Man_Con_a > 0 &'...
        '%s.Man_Con_d < MAXDUR & ' ...
        '%s.Man_Con_d > 0 &'...
        '%s.Man_Con_s < MAXSPEED & '...
        '%s.Trck_Con_s > 0 & '...
        '%s.Trck_Con_a < MAXAMP & ' ...
        '%s.Trck_Con_a > 0 & '...
        '%s.Trck_Con_d < MAXDUR & ' ...
        '%s.Trck_Con_d > 0 &'...
        '%s.Trck_Con_s < MAXSPEED & '...
        '%s.Trck_Con_s > 0' ...
        ');'],names{i},names{i},names{i},...
        names{i},names{i},names{i},...
        names{i},names{i},names{i},...
        names{i},names{i},names{i})...
        )    
    eval(sprintf('%s = %s(idx,:);',names{i},names{i}))
end

names = {'LeftMinPro','LeftMaxPro','RightMinPro','RightMaxPro'};
for i = 1:length(names)
    eval(...
    sprintf(['idx = find('...
        '%s.Man_pro_a > -MAXAMP & ' ...
        '%s.Man_pro_a < 0 &'...
        '%s.Man_pro_d < MAXDUR & ' ...
        '%s.Man_pro_d > 0 &'...
        '%s.Man_pro_s > -MAXSPEED & '...
        '%s.Trck_pro_s < 0 & '...
        '%s.Trck_pro_a > -MAXAMP & ' ...
        '%s.Trck_pro_a < 0 & '...
        '%s.Trck_pro_d < MAXDUR & ' ...
        '%s.Trck_pro_d > 0 &'...
        '%s.Trck_pro_s > -MAXSPEED & '...
        '%s.Trck_pro_s < 0' ...
        ');'],names{i},names{i},names{i},...
        names{i},names{i},names{i},...
        names{i},names{i},names{i},...
        names{i},names{i},names{i})...
        )    
    eval(sprintf('%s = %s(idx,:);',names{i},names{i}))
end        

Data.LeftMinCon = LeftMinCon;
Data.LeftMaxCon = LeftMaxCon;
Data.RightMinCon = RightMinCon;
Data.RightMaxCon = RightMaxCon;
Data.LeftMinPro = LeftMinPro;
Data.LeftMaxPro = LeftMaxPro;
Data.RightMinPro = RightMinPro;
Data.RightMaxPro = RightMaxPro;



%%
maxval =80;
nbins = 0.3*maxval+1;
bins = linspace(0,maxval,nbins);

BIN_AMPLITUDE_PRO = zeros(length(bins)-1, length(bins)-1);
for i = 1:length(bins)-1
    for j = 1:length(bins)-1
        BIN_AMPLITUDE_PRO(i,j) = numel(find( abs(LeftMaxPro.Man_pro_a) >= bins(j) & abs(LeftMaxPro.Man_pro_a) < bins(j+1) & ...
            abs(LeftMaxPro.Trck_pro_a) >= bins(i) & abs(LeftMaxPro.Trck_pro_a) < bins(i+1))) + BIN_AMPLITUDE_PRO(i,j);
        BIN_AMPLITUDE_PRO(i,j) = numel(find( abs(LeftMinPro.Man_pro_a >= bins(j)) & abs(LeftMinPro.Man_pro_a) < bins(j+1) & ...
            abs(LeftMinPro.Trck_pro_a) >= bins(i) & abs(LeftMinPro.Trck_pro_a) < bins(i+1))) + BIN_AMPLITUDE_PRO(i,j);
        BIN_AMPLITUDE_PRO(i,j) = numel(find( abs(RightMaxPro.Man_pro_a) >= bins(j) & abs(RightMaxPro.Man_pro_a) < bins(j+1) & ...
            abs(RightMaxPro.Trck_pro_a) >= bins(i) & abs(RightMaxPro.Trck_pro_a) < bins(i+1))) + BIN_AMPLITUDE_PRO(i,j);
        BIN_AMPLITUDE_PRO(i,j) = numel(find( abs(RightMinPro.Man_pro_a) >= bins(j) & abs(RightMinPro.Man_pro_a) < bins(j+1) & ...
            abs(RightMinPro.Trck_pro_a) >= bins(i) & abs(RightMinPro.Trck_pro_a) < bins(i+1))) + BIN_AMPLITUDE_PRO(i,j);

     end
end
Data.BIN_AMPLITUDE_PRO = BIN_AMPLITUDE_PRO./max(max(BIN_AMPLITUDE_PRO));


BIN_AMPLITUDE_CON = zeros(length(bins)-1, length(bins)-1);
for i = 1:length(bins)-1
    for j = 1:length(bins)-1
        
        BIN_AMPLITUDE_CON(i,j) = numel(find( abs(LeftMaxCon.Man_Con_a) >= bins(j) & abs(LeftMaxCon.Man_Con_a) < bins(j+1) & ...
            abs(LeftMaxCon.Trck_Con_a) >= bins(i) & abs(LeftMaxCon.Trck_Con_a) < bins(i+1))) + BIN_AMPLITUDE_CON(i,j);
        BIN_AMPLITUDE_CON(i,j) = numel(find( abs(LeftMinCon.Man_Con_a >= bins(j)) & abs(LeftMinCon.Man_Con_a) < bins(j+1) & ...
            abs(LeftMinCon.Trck_Con_a) >= bins(i) & abs(LeftMinCon.Trck_Con_a) < bins(i+1))) + BIN_AMPLITUDE_CON(i,j);
        BIN_AMPLITUDE_CON(i,j) = numel(find( abs(RightMaxCon.Man_Con_a) >= bins(j) & abs(RightMaxCon.Man_Con_a) < bins(j+1) & ...
            abs(RightMaxCon.Trck_Con_a) >= bins(i) & abs(RightMaxCon.Trck_Con_a) < bins(i+1))) + BIN_AMPLITUDE_CON(i,j);
        BIN_AMPLITUDE_CON(i,j) = numel(find( abs(RightMinCon.Man_Con_a) >= bins(j) & abs(RightMinCon.Man_Con_a) < bins(j+1) & ...
            abs(RightMinCon.Trck_Con_a) >= bins(i) &abs( RightMinCon.Trck_Con_a) < bins(i+1))) + BIN_AMPLITUDE_CON(i,j);
     end
end
Data.BIN_AMPLITUDE_CON = BIN_AMPLITUDE_CON./max(max(BIN_AMPLITUDE_CON));



stepsize = nbins/4;
Data.BIN_AMPLITUDExtickax  = [0:stepsize:size(BIN_AMPLITUDE_CON,1)+1];
Data.BIN_AMPLITUDExtickax(1) = 1;
Data.BIN_AMPLITUDExtickax(end) = Data.BIN_AMPLITUDExtickax(end)-1;
stepsize = 20;
Data.BIN_AMPLITUDExticklabels = 0:stepsize:maxval;

%%
maxval = 80;
nbins = 0.3*maxval+1;
bins = linspace(0,maxval,nbins);

BIN_DURATION_PRO = zeros(length(bins)-1, length(bins)-1);
for i = 1:length(bins)-1
    for j = 1:length(bins)-1
        BIN_DURATION_PRO(i,j) = numel(find( LeftMaxPro.Man_pro_d >= bins(j) & LeftMaxPro.Man_pro_d < bins(j+1) & ...
            LeftMaxPro.Trck_pro_d >= bins(i) & LeftMaxPro.Trck_pro_d < bins(i+1))) + BIN_DURATION_PRO(i,j);
        BIN_DURATION_PRO(i,j) = numel(find( LeftMinPro.Man_pro_d >= bins(j) & LeftMinPro.Man_pro_d < bins(j+1) & ...
            LeftMinPro.Trck_pro_d >= bins(i) & LeftMinPro.Trck_pro_d < bins(i+1))) + BIN_DURATION_PRO(i,j);
        BIN_DURATION_PRO(i,j) = numel(find( RightMaxPro.Man_pro_d >= bins(j) & RightMaxPro.Man_pro_d < bins(j+1) & ...
            RightMaxPro.Trck_pro_d >= bins(i) & RightMaxPro.Trck_pro_d < bins(i+1))) + BIN_DURATION_PRO(i,j);
        BIN_DURATION_PRO(i,j) = numel(find( RightMinPro.Man_pro_d >= bins(j) & RightMinPro.Man_pro_d < bins(j+1) & ...
            RightMinPro.Trck_pro_d >= bins(i) & RightMinPro.Trck_pro_d < bins(i+1))) + BIN_DURATION_PRO(i,j);        

    end
end

Data.BIN_DURATION_PRO = BIN_DURATION_PRO./max(max(BIN_DURATION_PRO));



BIN_DURATION_CON = zeros(length(bins)-1, length(bins)-1);
for i = 1:length(bins)-1
    for j = 1:length(bins)-1
        
        BIN_DURATION_CON(i,j) = numel(find( LeftMaxCon.Man_Con_d >= bins(j) & LeftMaxCon.Man_Con_d < bins(j+1) & ...
            LeftMaxCon.Trck_Con_d >= bins(i) & LeftMaxCon.Trck_Con_d < bins(i+1))) + BIN_DURATION_CON(i,j);
        BIN_DURATION_CON(i,j) = numel(find( LeftMinCon.Man_Con_d >= bins(j) & LeftMinCon.Man_Con_d < bins(j+1) & ...
            LeftMinCon.Trck_Con_d >= bins(i) & LeftMinCon.Trck_Con_d < bins(i+1))) + BIN_DURATION_CON(i,j);
        BIN_DURATION_CON(i,j) = numel(find( RightMaxCon.Man_Con_d >= bins(j) & RightMaxCon.Man_Con_d < bins(j+1) & ...
            RightMaxCon.Trck_Con_d >= bins(i) & RightMaxCon.Trck_Con_d < bins(i+1))) + BIN_DURATION_CON(i,j);
        BIN_DURATION_CON(i,j) = numel(find( RightMinCon.Man_Con_d >= bins(j) & RightMinCon.Man_Con_d < bins(j+1) & ...
            RightMinCon.Trck_Con_d >= bins(i) & RightMinCon.Trck_Con_d < bins(i+1))) + BIN_DURATION_CON(i,j);
    end
end

Data.BIN_DURATION_CON = BIN_DURATION_CON./max(max(BIN_DURATION_CON));


stepsize = nbins/4;
xtickax  = 0:stepsize:size(BIN_DURATION_CON,1)+1;
xtickax(1) = 1;
xtickax(end) = xtickax(end)-1;
Data.BIN_DURATIONxtickax = xtickax;
stepsize = 20;
Data.BIN_DURATIONxticklabels = 0:stepsize:maxval;











%%
maxval = 2000;

bins = linspace(0,maxval,nbins);

BIN_SPEED_PRO = zeros(length(bins)-1, length(bins)-1);
for i = 1:length(bins)-1
    for j = 1:length(bins)-1
        
        BIN_SPEED_PRO(i,j) = numel(find( abs(LeftMaxPro.Man_pro_s >= bins(j)) &  abs(LeftMaxPro.Man_pro_s) < bins(j+1) & ...
             abs(LeftMaxPro.Trck_pro_s) >= bins(i) &  abs(LeftMaxPro.Trck_pro_s) < bins(i+1))) + BIN_SPEED_PRO(i,j);
        BIN_SPEED_PRO(i,j) = numel(find(  abs(LeftMinPro.Man_pro_s) >= bins(j) & abs(LeftMinPro.Man_pro_s) < bins(j+1) & ...
            abs(LeftMinPro.Trck_pro_s) >= bins(i) & abs(LeftMinPro.Trck_pro_s) < bins(i+1))) + BIN_SPEED_PRO(i,j);
        BIN_SPEED_PRO(i,j) = numel(find( abs(RightMaxPro.Man_pro_s) >= bins(j) & abs(RightMaxPro.Man_pro_s) < bins(j+1) & ...
            abs(RightMaxPro.Trck_pro_s) >= bins(i) & abs(RightMaxPro.Trck_pro_s) < bins(i+1))) + BIN_SPEED_PRO(i,j);
        BIN_SPEED_PRO(i,j) = numel(find( abs(RightMinPro.Man_pro_s) >= bins(j) & abs(RightMinPro.Man_pro_s) < bins(j+1) & ...
            abs(RightMinPro.Trck_pro_s) >= bins(i) & abs(RightMinPro.Trck_pro_s) < bins(i+1))) + BIN_SPEED_PRO(i,j);
    end
end
Data.BIN_SPEED_PRO = BIN_SPEED_PRO./max(max(BIN_SPEED_PRO));


BIN_SPEED_CON = zeros(length(bins)-1, length(bins)-1);
for i = 1:length(bins)-1
    for j = 1:length(bins)-1
        
       BIN_SPEED_CON(i,j) = numel(find( abs(LeftMaxCon.Man_Con_s >= bins(j)) &  abs(LeftMaxCon.Man_Con_s) < bins(j+1) & ...
             abs(LeftMaxCon.Trck_Con_s) >= bins(i) &  abs(LeftMaxCon.Trck_Con_s) < bins(i+1))) + BIN_SPEED_CON(i,j);
        BIN_SPEED_CON(i,j) = numel(find(  abs(LeftMinCon.Man_Con_s) >= bins(j) & abs( LeftMinCon.Man_Con_s) < bins(j+1) & ...
            abs( LeftMinCon.Trck_Con_s) >= bins(i) & abs( LeftMinCon.Trck_Con_s) < bins(i+1))) + BIN_SPEED_CON(i,j);
        BIN_SPEED_CON(i,j) = numel(find(  abs(RightMaxCon.Man_Con_s) >= bins(j) &  abs(RightMaxCon.Man_Con_s) < bins(j+1) & ...
             abs(RightMaxCon.Trck_Con_s) >= bins(i) &  abs(RightMaxCon.Trck_Con_s) < bins(i+1))) + BIN_SPEED_CON(i,j);
        BIN_SPEED_CON(i,j) = numel(find(  abs(RightMinCon.Man_Con_s) >= bins(j) &  abs(RightMinCon.Man_Con_s) < bins(j+1) & ...
             abs(RightMinCon.Trck_Con_s) >= bins(i) &  abs(RightMinCon.Trck_Con_s) < bins(i+1))) + BIN_SPEED_CON(i,j);
    end
end
Data.BIN_SPEED_CON = BIN_SPEED_CON./max(max(BIN_SPEED_CON));


stepsize = nbins/2;
xtickax  = 0:stepsize:size(BIN_SPEED_CON,1)+1;
xtickax(1) = 1;
xtickax(end) = xtickax(end)-1;
Data.BIN_SPEEDxtickax = xtickax;
stepsize = 1000;
Data.BIN_SPEEDxticklabels = 0:stepsize:maxval;


%% AX 11 - Touch Data

h = waitbar(0, 'doing touch');

for i =  1:size(comp_files,1)
     load(fullfile(comp_files(i).folder, comp_files(i).name))
     
     
     nframes = size(Annotations.Tracker.Touch, 2);
     TTouch = [];
     TTouch(1:nframes) = 0;
     MTouch = [];
     MTouch(1:nframes) = 0;
     MTouch_auto = [];
     MTouch_auto(1:nframes) = 0;
     for j = 1:nframes        
        TTouch(j) = numel(find(Annotations.Tracker.Touch{j})); 
        if j <= size(Annotations.Manual.Touch.pt,2)
            MTouch(j) = size(Annotations.Manual.Touch.pt{j},1);
        end
        MTouch_auto(j) = numel(find(Annotations.Manual.Touch_auto{j}));         
     end
     
     if i == 1
         TOUCH = table(TTouch', MTouch', MTouch_auto','VariableNames',{'Tracker','Manual','ManualAuto'});
     else
         looptouch = table(TTouch', MTouch', MTouch_auto','VariableNames',{'Tracker','Manual','ManualAuto'});
         TOUCH  = [TOUCH;looptouch];
     end
    
    waitbar(i/37)
    
end
close(h)



idx = find( TOUCH.Tracker ~= 0 | TOUCH.Manual ~= 0);
TOUCH = TOUCH(idx,:);
Data.TOUCH = TOUCH;


maxval = max([TOUCH.Tracker;TOUCH.Manual;TOUCH.ManualAuto]);
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



save(fullfile(Datapath,'Data_Figure_Par_Eval'),'Data')




