%Parametersetup
%
%This GUI generates a 'Settings.mat' file that is used
%as input for tracking. Various parameters can be adjusted
%to optimize tracking results.


if exist('Settings\Settings.mat','file')
    tempvar.f = figure('Units','Points','Position',[500 500 230 100],'NumberTitle','off','Name','Load Settings');
    tempvar.b1 = uicontrol(tempvar.f,'Units','points','Position',[10 20 100 50],...
        'Style','pushbutton','Callback',@(src,eventdata)reportLoad(src,eventdata,tempvar.f,'y'),...
        'String','Yes');
    tempvar.b2 = uicontrol(tempvar.f,'Units','points','Position',[120 20 100 50],...
        'Style','pushbutton','Callback',@(src,eventdata)reportLoad(src,eventdata,tempvar.f,'n'),...
        'String','No');
    tempvar.text = uicontrol(tempvar.f,'Units','points','Position',[0 80 230 10],...
        'Style','text','String','Load previously saved settings?');
    uiwait;
    resp = guidata(tempvar.f);
    close(tempvar.f);
    clear tempvar
    
end

if strcmp(resp,'y')
    load('Settings\Settings.mat')
    handles.show.FullRoi = Settings.FullRoi;
else
    Settings.Gaussian_kernel_size = 3;
    Settings.doGaussian = 1;
    Settings.Gamma = 1;
    Settings.Background_threshold = 0.5;
    Settings.Edges_kernel_large = 10;
    Settings.Edges_kernel_small = 6;
    Settings.Edges_threshold = 0.35;
    Settings.Shape_threshold = 0.21;
    Settings.Dilation = 5;
    Settings.Seed_threshold = 0.05;
    Settings.Trace_threshold = 0.49;
    Settings.Trace_kernel_large = 5;
    Settings.Trace_kernel_small = 1;
    Settings.video_extension = '.dat';
    Settings.FullRoi = 0;
    
    handles.show.FullRoi = 0;
end


handles.State.folder = uigetdir('Select video folder');
%handles.State.folder = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
handles.State.videos = scanfiles(handles.State.folder, Settings.video_extension, 10, '');
handles.State.current_video = 1;
if isempty(handles.State.videos)
    disp('No videos found...')
    return
end
Settings.Video = fullfile(handles.State.videos(handles.State.current_video).folder,...
    handles.State.videos(handles.State.current_video).name);
handles.State.Video = Settings.Video;
handles.State.current_frame = 1;


Settings = getMetaData(Settings);


handles.Settings = Settings;

handles.prop.GK = [1 20 1];
handles.prop.GM = [0 5 0.01];
handles.prop.BG = [0 1 0.01];
handles.prop.EKL = [1 100 1];
handles.prop.EKS = [1 100 1];
handles.prop.ET = [0 1 0.01];
handles.prop.SHT = [0 1 0.01];
handles.prop.DL = [0 40 1];
handles.prop.SET = [0 1 0.01];
handles.prop.TT = [0 1 0.01];
handles.prop.TKL = [0 20 1];
handles.prop.TKS = [0 20 1];

handles.show.Objects = 1;
handles.show.Edges = 1;
handles.show.Threshold = 1;
handles.cmap = makeColor();



fig_heigth = 600;
fig_width = 800;

f = figure('Units','points','Position',[600 0 fig_width fig_heigth],'NumberTitle','off','Name','ParameterSetup');



ratio = 512/640;

axwidth = 0.35*fig_heigth;
axheigth = axwidth*ratio;

x_offset = 0.28*fig_width;
y_offset = 0.03*fig_heigth;

x_space = 0.1*fig_width;
y_space = 0.1*fig_heigth;

% Axes 1

prop = handles.prop;

handles.foldertext = uicontrol(f,'Units','points','Position',[0.03*fig_width 0.6*fig_heigth 180 100],...
    'Style','text','String','','HorizontalAlignment','left');

videos = cell(1,size(handles.State.videos,1));
for i = 1:size(handles.State.videos,1)
    videos{i} = handles.State.videos(i).name;
end

handles.videomenu = uicontrol(f, 'Units','points','Position',[0.03*fig_width 0.6*fig_heigth 180 30],...
    'Style','popupmenu','String',videos,'Callback',@(src,eventdata)selectVideo(src,eventdata,f));
uicontrol(f,'Units','points','Position',[0.03*fig_width 0.6*fig_heigth+30 180 10],'Style','text','String','Video Select',...
    'Horizontalalignment','center','TooltipString','Select video for example tracking')

handles.frameslider = uicontrol(f, 'Units','points','Position',[0.03*fig_width 0.6*fig_heigth-20 180 15],...
    'Style','slider','Value',1,'Callback',@(src,eventdata)selectFrame(src,eventdata,f),...
    'Min',1,'Max',Settings.Nframes,'Sliderstep',ones(1,2)*(1/numel(1:Settings.Nframes)));
uicontrol(f,'Units','points','Position',[0.03*fig_width 0.6*fig_heigth 180 10],'Style','text','String','Frame Select',...
    'Horizontalalignment','center','TooltipString','Select frame for example tracking')

handles.frametext = uicontrol(f,'Units','points','Position',[0.03*fig_width 0.6*fig_heigth-40 180 10],...
    'Style','text','String','Current frame: 1','HorizontalAlignment','left');


uicontrol(f,'Units','points','Position',[0.03*fig_width 0.6*fig_heigth-100 100 30],...
    'Style','pushbutton','Callback',@(src,eventdata)exportSettings(src,eventdata,f),...
    'String','Save','TooltipString','Save current settings');

uicontrol(f,'Units','points','Position',[0.03*fig_width 0.6*fig_heigth-180 100 30],...
    'Style','pushbutton','Callback',@(src,eventdata)exitFig(src,eventdata,f),...
    'String','Quit','TooltipString','Save current settings');


x = x_offset;
y = fig_heigth - y_offset - axwidth;

handles.axes1 = axes('Units','points','Position',[x y axwidth axheigth]);

handles.slider1 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-25,axheigth,10],'Style','slider',...
    'Value',Settings.Gaussian_kernel_size,'Callback',@(src,eventdata)readValInt(src,eventdata,f,'Gaussian_kernel_size'),...
    'Min',prop.GK(1),'Max',prop.GK(2),'Sliderstep',ones(1,2)*(1/numel(prop.GK(1):prop.GK(3):prop.GK(2))));

uicontrol(f,'Units','points','Position',[x y-25 0.2*axwidth 10],'Style','text','String','BG. t.',...
    'Horizontalalignment','Left','TooltipString','Preprocessing: Gaussian kernel size')


x = x_offset + x_space + axheigth;
handles.axes2 = axes('Units','points','Position',[x y axwidth axheigth]);
handles.axes3 = axes('Units','points','Position',[x y axwidth axheigth]);
handles.slider2 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-25,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Background_threshold,'Callback',@(src,eventdata)readVal2D(src,eventdata,f,'Background_threshold'),...
    'Min',prop.BG(1),'Max',prop.BG(2),'Sliderstep',ones(1,2)*(1/numel(prop.BG(1):prop.BG(3):prop.BG(2))));
uicontrol(f,'Units','points','Position',[x y-25 0.2*axwidth 10],'Style','text','String','BG. t.',...
    'Horizontalalignment','Left','TooltipString','Background detection: Threshold')

handles.slider3 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-40,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Edges_kernel_small,'Callback',@(src,eventdata)readValInt(src,eventdata,f,'Edges_kernel_small'),...
    'Min',prop.EKS(1),'Max',prop.EKS(2),'Sliderstep',ones(1,2)*(1/numel(prop.EKS(1):prop.EKS(3):prop.EKS(2))));
uicontrol(f,'Units','points','Position',[x y-40 0.2*axwidth 10],'Style','text','String','EG. ks.',...
    'Horizontalalignment','Left','TooltipString','Edge detection: Laplacian centre size')

handles.slider4 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-55,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Edges_kernel_large,'Callback',@(src,eventdata)readValInt(src,eventdata,f,'Edges_kernel_large'),...
    'Min',prop.EKL(1),'Max',prop.EKL(2),'Sliderstep',ones(1,2)*(1/numel(prop.EKL(1):prop.EKL(3):prop.EKL(2))));
uicontrol(f,'Units','points','Position',[x y-55 0.2*axwidth 10],'Style','text','String','EG. kl.',...
    'Horizontalalignment','Left','TooltipString','Edge detection: Laplacian neighbourhood size')

handles.slider5 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-70,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Edges_threshold,'Callback',@(src,eventdata)readVal2D(src,eventdata,f,'Edges_threshold'),...
    'Min',prop.ET(1),'Max',prop.ET(2),'Sliderstep',ones(1,2)*(1/numel(prop.ET(1):prop.ET(3):prop.ET(2))));
uicontrol(f,'Units','points','Position',[x y-70 0.2*axwidth 10],'Style','text','String','EG. t.',...
    'Horizontalalignment','Left','TooltipString','Edge detection: Threshold')

handles.checkbox2 = uicontrol(f, 'Units','points','Position',[x+axwidth+20 y+axheigth-20 80 10],...
    'Style','checkbox','String','show Objects','Value',1,'Callback',@(src,eventdata)readBox(src,eventdata,f,'Objects'));
handles.checkbox3 = uicontrol(f, 'Units','points','Position',[x+axwidth+20 y+axheigth-40 80 10],...
    'Style','checkbox','String','show Edges','Value',1,'Callback',@(src,eventdata)readBox(src,eventdata,f,'Edges'));


x = x_offset;
y = fig_heigth - y_offset - 2*axwidth - y_space;
handles.axes4 = axes('Units','points','Position',[x y axwidth axheigth]);
handles.axes5 = axes('Units','points','Position',[x y axwidth axheigth]);
handles.slider6 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-25,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Shape_threshold,'Callback',@(src,eventdata)readVal2D(src,eventdata,f,'Shape_threshold'),...
    'Min',prop.TT(1),'Max',prop.TT(2),'Sliderstep',ones(1,2)*(1/numel(prop.TT(1):prop.TT(3):prop.TT(2))));
uicontrol(f,'Units','points','Position',[x y-25 0.2*axwidth 10],'Style','text','String','SE. t.',...
    'Horizontalalignment','Left','TooltipString','Shape extraction: Threshold')

handles.slider7 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-40,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Dilation,'Callback',@(src,eventdata)readValInt(src,eventdata,f,'Dilation'),...
    'Min',prop.DL(1),'Max',prop.DL(2),'Sliderstep',ones(1,2)*(1/numel(prop.DL(1):prop.DL(3):prop.DL(2))));
uicontrol(f,'Units','points','Position',[x y-40 0.2*axwidth 10],'Style','text','String','SE. dil.',...
    'Horizontalalignment','Left','TooltipString','Shape extraction: Dilation')

handles.slider8 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-55,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Seed_threshold,'Callback',@(src,eventdata)readVal2D(src,eventdata,f,'Seed_threshold'),...
    'Min',prop.SET(1),'Max',prop.SET(2),'Sliderstep',ones(1,2)*(1/numel(prop.SET(1):prop.SET(3):prop.SET(2))));
uicontrol(f,'Units','points','Position',[x y-55 0.2*axwidth 10],'Style','text','String','SD. t.',...
    'Horizontalalignment','Left','TooltipString','Seed detection: Threshold')
handles.checkbox4 = uicontrol(f, 'Units', 'points','Position',[x-60 y+axheigth-20 80 10],...
    'Style','checkbox','String','full ROI','Value',handles.Settings.FullRoi,'Callback',@(src,eventdata)readBox(src,eventdata,f,'FullRoi'));


x = x_offset + x_space + axheigth;
handles.axes6 = axes('Units','points','Position',[x y axwidth 0.8*axwidth]);
handles.slider9 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-25,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Trace_threshold,'Callback',@(src,eventdata)readVal2D(src,eventdata,f,'Trace_threshold'),...
    'Min',prop.SHT(1),'Max',prop.SHT(2),'Sliderstep',ones(1,2)*(1/numel(prop.SHT(1):prop.SHT(3):prop.SHT(2))));
uicontrol(f,'Units','points','Position',[x y-25 0.2*axwidth 10],'Style','text','String','TT. t.',...
    'Horizontalalignment','Left','TooltipString','Trace tracking: Threshold')

handles.slider10 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-40,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Trace_kernel_small,'Callback',@(src,eventdata)readValInt(src,eventdata,f,'Trace_kernel_small'),...
    'Min',prop.TKS(1),'Max',prop.TKS(2),'Sliderstep',ones(1,2)*(1/numel(prop.TKS(1):prop.TKS(3):prop.TKS(2))));
uicontrol(f,'Units','points','Position',[x y-40 0.2*axwidth 10],'Style','text','String','TT. ks.',...
    'Horizontalalignment','Left','TooltipString','Trace tracking: Laplacian neighbourhood size')

handles.slider11 = uicontrol(f,'Units','points','Position',[x+0.2*axwidth,y-55,0.8*axwidth,10],'Style','slider',...
    'Value',Settings.Trace_kernel_large,'Callback',@(src,eventdata)readValInt(src,eventdata,f,'Trace_kernel_large'),...
    'Min',prop.TKL(1),'Max',prop.TKL(2),'Sliderstep',ones(1,2)*(1/numel(prop.TKL(1):prop.TKL(3):prop.TKL(2))));
uicontrol(f,'Units','points','Position',[x y-55 0.2*axwidth 10],'Style','text','String','TT. kl.',...
    'Horizontalalignment','Left','TooltipString','Trace tracking: Laplacian centre size')

handles.checkbox1 = uicontrol(f, 'Units', 'points','Position',[x+axwidth+20 y+axheigth-20 80 10],...
    'Style','checkbox','String','show threshold','Value',1,'Callback',@(src,eventdata)readBox(src,eventdata,f,'Threshold'));



guidata(f, handles)

readFrames(f)
readFrame(f)
updatePath(handles.foldertext, handles.State.folder)
updateFigure(f)


function files = scanfiles(path,extension, nlayers, format)
files = [];

for i = 1:nlayers
    str = path;
    
    for j = 1:i-1
        str = [str '\*'];
    end
    str = [str '\*' extension];
    
    loopfiles = dir(str);
    if ~isempty( format )
        for j = 1:size(loopfiles,1)
            token = regexp(loopfiles(j).name,format,'names');
            if ~isempty(token)
                files = [files; loopfiles(j)];
                
            end
        end
    else
        for j = 1:size(loopfiles,1)
            ext = loopfiles(j).name(end-3:end);
            if strcmp(ext, extension)
                files = [files;loopfiles(j)];
            end
        end
        %files = [files; loopfiles];
    end
end
end

function readFrames(f)
handles = guidata(f);
Settings = handles.Settings;
%fid = fopen(handles.State.Video,'r');
%fdim = [Settings.Video_width, Settings.Video_heigth];
warning('off')
nframes = Settings.Nframes;
idx = round(linspace(1,nframes,100));
Frames = zeros(512, 640, length(idx));
for i = 1:length(idx)
    Settings.Current_frame = idx(i)-1;
    Frames(:,:,i) = LoadFrame(Settings);
end
handles.sumFrames = sum(Frames,3);
handles.meanFrames = mean(Frames,3);

guidata(f, handles);
end

function readFrame(f)
handles = guidata(f);
Settings = handles.Settings;
Settings.Current_frame = handles.State.current_frame;
handles.FrameIn = LoadFrame(Settings);
guidata(f,handles);

end

function selectVideo(h,~,f)
handles = guidata(f);
handles.State.current_video = get(h,'Value');
handles.State.Video = fullfile(handles.State.videos(handles.State.current_video).folder,...
    handles.State.videos(handles.State.current_video).name);
handles.Settings.Video = handles.State.Video;
handles.Settings = getMetaData(handles.Settings);
handles.State.current_frame = 1;
set(handles.frameslider,'Max',handles.Settings.Nframes);
set(handles.frameslider,'Sliderstep',ones(1,2)*(1/numel(1:handles.Settings.Nframes)))
set(handles.frametext,'String',sprintf('Current frame: %d',handles.State.current_frame))
set(handles.frameslider,'Value',handles.State.current_frame)
guidata(f, handles)
readFrame(f)
readFrames(f)
updateFigure(f)
end

function selectFrame(h,~,f)
handles=guidata(f);
handles.State.current_frame = round(get(h, 'Value'));
guidata(f, handles)
readFrame(f)
updateFigure(f)
set(handles.frametext,'String',sprintf('Current frame: %d',handles.State.current_frame))

end

function reportLoad(~,~,f,res)
guidata(f,res);
uiresume;
end


function readBox(h,~,f,field)
handles = guidata(f);

if get(h,'Value') == 1
    handles.show.(field) = 1;
else
    handles.show.(field) = 0;
end

if strcmp(field, 'FullRoi')
    handles.Settings.FullRoi = handles.show.FullRoi;
end

guidata(f,handles);
updateFigure(f);
end

function readValInt(h,~,f,field)
handles = guidata(f);
val = round(get(h, 'Value'));
if mod(val,2) == 0
    val = val+1;
end
handles.Settings.(field) = val;
guidata(f, handles)
updateFigure(f)
%fprintf('%s - %f\n',field, handles.Settings.(field))
end

function readVal2D(h,~,f,field)
handles = guidata(f);
handles.Settings.(field) = 100*get(h, 'Value')/100;
guidata(f, handles)
updateFigure(f)
%fprintf('%s - %f\n',field, handles.Settings.(field))

end

function updatePath(obj,folder)
str = 'Current Path:\n';
tokens = regexp(folder,'\','split');
for j = 1:size(tokens,2)
    n(j) = length(tokens{j});
end

j = 1;
while j <= length(n)
    
    count = n(j);
    if count > 20
        str = [str tokens{j} '\\\n'];
        j = j+1;
    else
        offset = j;
        j = j+1;
        while count <= 20 & j <= length(n)
            count = count+n(j);
            j = j + 1;
        end
        
        
        for k = offset:j-2
            str = [str tokens{k} '\\'];
        end
        str = [str tokens{j-1} '\\' '\n'];
        
    end
    
    
end

str = sprintf(str);
set(obj,'String',str);


end


function exportSettings(~,~,f)

handles = guidata(f);
Settings = handles.Settings;

Settings.format = 'M(?<MOUSE>\d+)_R(?<SESSION>\d+)_(?<TRIAL>\d+).dat';
Settings.use_parfor = 1;
Settings.track_nose = 1;
Settings.costum_background = 1;
Settings.outpath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
Settings.n_background_samples = 30; % number of sample frames to extract background
Settings.frame_select = 'use_file';
Settings.FrameFile = 'Frames_with_nose.mat';
Settings.circle_start = -25;
Settings.circle_end = 25;
Settings.stepsize = 5;
Settings.extrapolationsize = 7;
Settings.dist_from_edge = 5;
Settings.minimum_traclength = 8;

Settings = orderfields(Settings); %#ok<NASGU>
save('Settings\Settings.mat','Settings')
end

function exitFig(~,~,f)
tempvar.f = figure('Units','Points','Position',[500 500 230 100],'NumberTitle','off','Name','Exit');
tempvar.b1 = uicontrol(tempvar.f,'Units','points','Position',[10 20 100 50],...
    'Style','pushbutton','Callback',@(src,eventdata)reportLoad(src,eventdata,tempvar.f,'y'),...
    'String','Yes');
tempvar.b2 = uicontrol(tempvar.f,'Units','points','Position',[120 20 100 50],...
    'Style','pushbutton','Callback',@(src,eventdata)reportLoad(src,eventdata,tempvar.f,'n'),...
    'String','No');
tempvar.text = uicontrol(tempvar.f,'Units','points','Position',[0 80 230 10],...
    'Style','text','String','Save updated settings?');

uiwait
resp = guidata(tempvar.f);
close(tempvar.f);


if strcmp(resp,'y')
    exportSettings('','',f);
end


close(f)
return

end





function updateFigure(f)
handles = guidata(f);


cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes3)
cla(handles.axes4)
cla(handles.axes5)
cla(handles.axes6)


c = handles.cmap;
Settings = handles.Settings;
show = handles.show;
%% Panel 1
Frame = im2double( abs(handles.FrameIn));

Frame = imadjust(Frame, [],[], Settings.Gamma);

if Settings.doGaussian
    ng = Settings.Gaussian_kernel_size;
    kfil = zeros(1, ng);
    kfil(1:ceil(ng/2)) = 1:ceil(ng/2);
    kfil(ceil(ng/2)+1:end) = max(kfil)-1:-1:1;
    Kgaus = repmat(kfil,[ng, 1]) + repmat(kfil,[ng, 1])';
    Kgaus = Kgaus./ sum(sum(Kgaus));
    Frame = conv2(Frame, Kgaus, 'same');
end



Frame = adapthisteq(Frame, 'NBins',256,'NumTiles',[5 5]);
minval = min(min(Frame));
Frame = Frame - minval;
Frame = Frame./max(max(Frame));


imagesc(handles.axes1, Frame)
colormap(handles.axes1, 'gray')
axis(handles.axes1, 'off')



%% Panel 2
Background = handles.sumFrames./100;

Background = Background - min(min(Background));
Background = Background ./ max(max(Background));


Background(Background > Settings.Background_threshold) = 0;
Objects = zeros(size(Background));

Objects(Background > 0) = 1;
KL = conv2(handles.meanFrames, ...
    ones(Settings.Edges_kernel_large,Settings.Edges_kernel_large)./...
    Settings.Edges_kernel_large,'same');
KS = conv2(handles.meanFrames, ...
    ones(Settings.Edges_kernel_small, Settings.Edges_kernel_small)./...
    Settings.Edges_kernel_small^2, 'same');

H = KS./KL;
H = H-min(min(H));
H = H./max(max(H));

tax = 0:0.01:1;
for i = 1:length(tax)
    count(i) = numel(find(H < tax(i))); %#ok<AGROW>
end
[~, idx] = max(diff(count));
edge_threshold = tax(idx)-0.02;

Edges = zeros(size(Objects));

%Edges(H < Settings.Edges_threshold) = 1;
Edges(H < edge_threshold) = 1;


Res = Costumbackground(Objects, Edges);
Objects = Res.Objects;
Edges = Res.Edges;

imagesc(handles.axes2, Frame);
colormap(handles.axes2, 'gray')
axis(handles.axes2, 'off')
cla(handles.axes3);

if show.Objects
    img = imagesc(handles.axes3, Objects);
    colormap(handles.axes3, c.Objects);
    img.AlphaData = 0.5*Objects;
    set(handles.axes3,'Visible','off');
end

if show.Edges
    hold(handles.axes3, 'on');
    img = imagesc(handles.axes3, Edges);
    colormap(handles.axes3, c.Objects);
    img.AlphaData = 1*Edges;
    set(handles.axes3, 'Visible','off')
end

%% Panel 3
SLN = zeros(size(Objects));
SLN( handles.FrameIn <= Settings.Shape_threshold) = 1;
SLN(Objects == 1) = 0;

SLS = imerode(SLN, strel('diamond', 5)); % Small silhouette
SLL = imdilate(SLS, strel('diamond', Settings.Dilation));
%SLN = imdilate(SLN, strel('diamond', 1));

Frame(1,1) = 999; % Default 'dead' pixel, to refer to if tracking should
% be interrupted.

K =  Settings.Trace_kernel_large;
L =  Settings.Trace_kernel_small;
CHl = conv2(Frame, ones(K,K)./K^2, 'same'); % Low pass filter frame with large kernel
CHs = conv2(Frame, ones(L,L)./L^2, 'same');  % Low pass filter frame with small kernel
CH = zeros(size(Objects));
CH(0.5*CHs./CHl < Settings.Trace_threshold) = 1; % Matrix with potential valid pixels
CH(Objects == 1) = 0;
CH(Edges == 1) = 0.5;
%CH(SLL == 1) = 0;
CH(SLN == 1) = 0;


% Find tracing seeds
RoiInd = find(edge(SLL));
if isempty(RoiInd); TraceOut = {}; return; end %#ok<NASGU>

[RoiSub(:,1), RoiSub(:,2)] = ind2sub(size(Frame), RoiInd);
RS2 = RoiSub;
RS2(:,1) = RS2(:,1)+1;
RoiSub = [RoiSub; RS2];
todo = 1:size(RoiSub,1);
todo(1) = NaN;
did(1:length(RoiInd)) = NaN;
did(1) = 1;
TempVar = zeros(length(RoiInd), 2);
TempVar(1,:) = RoiSub(1,:);

for i = 2:size(RoiSub,1)
    todo = todo(~isnan(todo));
    dist = sqrt(sum( (RoiSub(todo,:) - TempVar(end,:)).^2, 2));
    [~, idx] = min(dist);
    TempVar(i,:) = RoiSub(todo(idx),:);
    did(i) = todo(idx);
    todo(idx) = NaN;
end
FrameHeight = size(Frame, 1);
FrameWidth = size(Frame, 2);


TempVar(TempVar(:,1) < 1, :) = NaN;
TempVar(TempVar(:,1) > FrameHeight, :) = NaN;
TempVar(TempVar(:,2) < 1, :) = NaN;
TempVar(TempVar(:,2) > FrameWidth,:) = NaN;
TempVar = TempVar(~isnan(TempVar(:,1)),:);

RoiSub = TempVar;

%scatter(RoiSub(:,2), RoiSub(:,1), 'm','filled')

RoiInd = sub2ind(size(Frame), TempVar(:,1), TempVar(:,2));

SeedProfile = double(Frame(RoiInd)); % Intensity profile over ROI

LP = medfilt1(SeedProfile, 50);
SeedProfile = abs(SeedProfile-LP); % Subtract local background



if isempty(SeedProfile); TracesOut = {}; return; end %#ok<NASGU>

handles.Settings.FullRoi = handles.show.FullRoi;
if handles.Settings.FullRoi == 0
    [~, SeedInd] = findpeaks(SeedProfile, 'MinPeakDistance',3,'MinPeakProminence',Settings.Seed_threshold);
    Seeds = RoiSub(SeedInd,:);
elseif handles.Settings.FullRoi == 1
    %[~, SeedInd] = findpeaks(SeedProfile, 'MinPeakDistance',1,'Threshold',0);
    
    Seeds = RoiSub(find(CH(RoiInd)==1),:);
end


imagesc(handles.axes4, handles.FrameIn)
colormap(handles.axes4, 'gray')
img = imagesc(handles.axes5, SLL);
colormap(handles.axes5, c.Snout)
img.AlphaData = 0.5*SLL;
set(handles.axes4, 'Visible', 'off')
set(handles.axes5, 'Visible','off')
hold(handles.axes5, 'on')
scatter(handles.axes5, RoiSub(:,2), RoiSub(:,1),5, ...
    'MarkerFaceColor',c.Roi,'MarkerEdgeColor',c.Roi)

scatter(handles.axes5, Seeds(:,2), Seeds(:,1), 5, ...
    'MarkerFaceColor',c.Trace,'MarkerEdgeColor',c.Trace)



%% Panel 4
% Track Traces
bufferSize = 100; % max nr of points on a trace

Traces = ones(2, size(Seeds,1), bufferSize);
Traces(:, :, 1) = Seeds(:, 1:2)';
H = -15:15;
HH = -30:30;
trace_finished = zeros(1, size(Seeds,1));
trace_edges = zeros(1, size(Seeds,1));




for i = 2:bufferSize
    
    if ~any(any(Traces(:,:, i-1) ~= 1)) % If no points were assigned in previous loop
        break
    end
    
    Theta = [];
    if i == 2
        ThetaRoi = repmat(1:20:360, [size(Seeds,1), 1]); % Roi bounded by angle range
        %N = repmat(Centre,[size(Traces,2) 1]);
        %vt = Traces(:,:,i-1) - N';
        %ThetaRoi = atan2d(vt(2,:), vt(1,:))' + [-60:20:60];
    else
        vt = Traces(:,:,i-1)- Traces(:,:,i-2);
        ThetaRoi = atan2d(vt(2,:), vt(1,:))' + H;
        
    end
    
    
    
    Theta(1,:,:) = [3*cosd(ThetaRoi)]; %#ok<*NBRAK> % ROI per previous point
    Theta(2,:,:) = [3*sind(ThetaRoi)];
    
    LastPoints = repmat(Traces(:,:,i-1), [1 , 1, size(Theta, 3)]);
    
    RoiSub = round(LastPoints + Theta);
    
    % Remove indices outside of frame range
    RoiSub(RoiSub < 1) = 1;
    RoiSub(1,RoiSub(1,:,:) > FrameHeight) = 1;
    RoiSub(2,RoiSub(2,:,:) > FrameWidth) = 1;
    RoiInd = squeeze(sub2ind(size(Frame), RoiSub(1,:,:), RoiSub(2,:,:)));
    
    if i == 2
        [a,b] = find(SLL(RoiInd));
        %RoiInd(SLL(RoiInd) == 1) = 1;
        for j = 1:size(RoiInd, 1)
            
            
            idx = b(a == j);
            if ~any(diff(idx)>1) & min(idx) > 2 & max(idx) < size(RoiInd,2)-1 %#ok<*AND2>
                idx = [[-2 -1]'+idx(1); idx; [1 2]'+idx(end)];
            elseif ~any(diff(idx) > 1) & min(idx) == 2 & max(idx) < size(RoiInd, 2) - 1
                idx = [-1+idx(1); idx; [1 2]'+idx(end); size(RoiInd,2)];
            elseif ~any(diff(idx) > 1) & min(idx) == 1 & max(idx) < size(RoiInd, 2) - 1
                idx = [idx; [1 2]'+idx(end);[-1 0]'+size(RoiInd,2)];
                
            elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(RoiInd,2)
                idx = [ [1 2]'; [-2 -1]'+idx(1); idx];
            elseif ~any(diff(idx) > 1) & min(idx) > 2  & max(idx) == size(RoiInd,2)-1
                idx = [ [1]'; [-2 -1]'+idx(1); idx; size(RoiInd,2)];
                
            elseif numel(find(diff(idx) > 1)) == 1 & min(idx) == 1 & max(idx) == size(RoiInd,2)
                id = find(diff(idx) > 1);
                idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
                
            elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) < size(RoiInd,2) - 1
                idx = [ [-2 -1]'+min(idx); idx ; [1 2]'+idx(end)];
                
            elseif numel(find(diff(idx) > 1)) == 1 & min(idx) >2 & max(idx) == size(RoiInd,2) - 1
                idx = [1; [-2 -1]'+min(idx); idx ; [1 ]'+idx(end)];
                
                
            elseif numel(find(diff(idx) > 1)) > 1 & min(idx) == 1 & max(idx) == size(RoiInd,2)
                [~, id] = max(diff(idx));
                idx = [idx(1:id); [1 2]'+idx(id); [-2 -1]'+idx(id+1);idx(id+1:end)];
                
                
            else
                idx =idx; %#ok<ASGSL>
                
                
            end
            
            RoiInd(j,idx) = 1;
            
            
            
            
            
        end
        
    end
    
    Iroi = Frame(RoiInd);
    Iroi(CH(RoiInd) == 0) = NaN;
    if size(Seeds,1) == 1
        [~, PointSub] = min(Iroi);
        PointInd = PointSub;
    else
        [~, PointSub] = min(Iroi, [], 2);
        PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
    end
    
    
    [Traces(1, :, i), Traces(2,:,i)] = ind2sub(size(Frame), RoiInd(PointInd));
    
    trace_finished(CH(RoiInd(PointInd)) == 0) = trace_finished(CH(RoiInd(PointInd)) == 0) + 1;
    trace_edges(CH(RoiInd(PointInd)) == 0.5) = trace_edges(CH(RoiInd(PointInd)) == 0.5) + 1;
    
    
    
    if i < 3
        trace_finished(CH(RoiInd(PointInd)) == 0.5) =  trace_finished(CH(RoiInd(PointInd)) == 0.5) + 1;
    elseif any(CH(RoiInd(PointInd)) == 0.5)
        nan_traces = find(CH(RoiInd(PointInd)) == 0.5);
        
        
        
        for j = 1:length(nan_traces)
            
            
            t = squeeze(Traces(:, nan_traces(j), :));
            i1 = find(t(1,:) == 1, 1, 'first');
            i2 = find(t(2,:) == 1, 1, 'first');
            
         
            if ~isempty([i1 i2])
                t = t(:,1:max([i1 i2])-1);
            end
            
            if i < 10
                rawax = 1:i;
                fitax = 1:0.5:i+10;
            else
                rawax = 1:10;
                fitax = 1:0.5:i+10;
            end
            
            
            if size(t, 2) >= length(rawax)
                px = polyfit(rawax, t(1,end-length(rawax)+1:end), 1);
                py = polyfit(rawax, t(2,end-length(rawax)+1:end), 1);
                
                tfit = [];
                tfit(1,:) = round(polyval(px, fitax));
                tfit(2,:) = round(polyval(py, fitax));
                
                throwidx = [find(tfit(1,:) < 1), find(tfit(1,:) > FrameHeight)];
                throwidx = [throwidx, find(tfit(2,:) < 1), find(tfit(2,:) > FrameWidth)];
                
                if ~isempty(throwidx)
                    tfit(:, throwidx) = NaN;
                    tfit = tfit(:, ~isnan(tfit(1,:)));
                end
                
                tfitInd = sub2ind(size(Frame), tfit(1,:), tfit(2,:));
                ptInd = find(CH(tfitInd) == 0.5,1,'last');
                
                if ptInd < size(tfit,2)-5
                    LastPt = tfit(:, ptInd+1);
                    %ThetaLoop = squeeze(Theta(:, nan_traces(j), :));
                    
                    vt = tfit(:, end) - tfit(:, end-5);
                    T1 = atan2d(vt(2,1), vt(1,1)) + HH;
                    T2 = [];
                    T2(1,:) = 6*cosd(T1)';
                    T2(2,:) = 6*sind(T1)';
                    
                    
                    
                    RoiSub = round(LastPt + T2);
                    RoiSub(RoiSub < 1) = 1;
                    RoiSub(1,RoiSub(1,:) > FrameHeight) = 1;
                    RoiSub(2,RoiSub(2,:) > FrameWidth) = 1;
                    RoiInd = squeeze(sub2ind(size(Frame), RoiSub(1,:), RoiSub(2,:)));
                    Iroi = Frame(RoiInd);
                    
                    [~, PointSub] = min(Iroi, [], 2);
                    PointInd = sub2ind(size(RoiInd), 1:length(PointSub), PointSub');
                    [a,b] = ind2sub(size(Frame), RoiInd(PointInd));
                    
                    if CH(a,b) == 1
                        Traces(1, nan_traces(j), i) = a;
                        Traces(2, nan_traces(j), i) = b;
                        trace_edges(nan_traces(j)) = 0;
                        
                    else
                        
                        trace_finished( nan_traces(j) ) = trace_finished( nan_traces(j) )  + 1;
                    end
                    
                end
            end
        end
        
    end
    
    
    
    Traces(1:2,trace_finished>2,i) = 1;
    
    
    
    idx = find(trace_edges > 10 & Traces(1,:,i)>1);
    if ~isempty(idx)
        for j = 1:length(idx)
            n_delete = trace_edges(idx(j));
            Traces(:, idx(j), i-n_delete:i) = 1;
        end
    end
    
    
    
end

Traces(Traces == 1) = NaN;



TOut = {};
TracesOut = {};
Or = [];
for i = 1:size(Traces,2)
    t = squeeze(Traces(:,i,:));
    idx = find(~isnan(t(1,:)),1,'last');
    t = t(:, 1:idx);
    if size(t,2) > 10
        TOut{end+1} = t';
        Or(end+1,:) = t(:,1)';
    end
    
    
end


if ~isempty(Or)
    dx = abs(Or(:,1) - Or(:,1)');
    dy = abs(Or(:,2) - Or(:,2)');
    d = sqrt( (dx+dy).^2 );
    
    
    
    flag = ones(1, size(TOut, 2));
    for i = 1:length(flag)
        d(i,i) = NaN;
        idx = find(d(i,:) < 10);
        if ~isempty(idx)
            t1 = TOut{i};
            
            for j = 1:length(idx)
                t2 = TOut{idx(j)};
                dx = abs(t1(:,1) - t2(:,1)');
                dy = abs(t1(:,2) - t2(:,2)');
                dnew = sqrt( (dx+dy).^2);
                %val = min(dnew,[],2);
                dist = mean(min(dnew,[],2));
                
                if dist < 3
                    flag(i) = 0;
                    d(:,i) = NaN;
                end
            end
        end
    end
    
    
    for i = 1:length(flag)
        if flag(i)
            TracesOut{end+1} = TOut{i}; %#ok<*AGROW>
        end
    end
end

if ~handles.show.Threshold
    imagesc(handles.axes6, handles.FrameIn)
else
    imagesc(handles.axes6, CH);
end
colormap(handles.axes6, 'gray')
hold(handles.axes6, 'on')
set(handles.axes6, 'Visible', 'off')
for i = 1:size(TracesOut, 2)
    t = TracesOut{i};
    plot(handles.axes6, t(:,2), t(:,1), 'color',c.Red,'LineWidth',3)
    plot(handles.axes6, t(:,2), t(:,1), 'color',c.Trace,'LineWidth',2)
end



end






