function varargout = verify_touch(varargin)
%VERIFY_TOUCH
%
% validate touch detection for Manual, MWT and Whisk data.
%
%   GUI:
%  figure(1) - display target platform with detected touches
%  figure(2) - legend for touch data:
%       - manual touches (red)
%       - MWT touches (blue)
%       - Whisk touches (green)
%
%   skips empty frames
%   for each frame fill in the 4 fields:
%   - # Touch - total nr of touches in frame
%   - # Correct Man - total nr of correct manual notated touches
%   - # Correct MWT - total nr of correct MWT notated touches
%   - # Correct Jan - total nr of correct Whisk notated touches
%
%   After filling the fields, press 'SUBMIT' (submit won't allow
%   empty fields). Data is saved after each frame, on startup, check
%   for previous saved inspection, ask to load and continue last
%   session.
%
%  Use QUIT to exit the GUI
%
%
%
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @verify_touch_OpeningFcn, ...
                   'gui_OutputFcn',  @verify_touch_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before verify_touch is made visible.
function verify_touch_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to verify_touch (see VARARGIN)

% Choose default command line output for verify_touch
handles.output = hObject;

% Select '_Annotations_Tracker' file
[handles.FileName, handles.PathName] = uigetfile('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\*_compiled.mat');
load_data = load(fullfile( handles.PathName, handles.FileName));
handles.Data = load_data.Annotations;
handles.SaveName = fullfile( handles.PathName, [handles.Data.Settings.FileName(1:end-4) '_Evaluation_touch.mat']);
view_heigth = 100;
max_heigth = size(handles.Data.Output.Objects, 1);

switch(handles.Data.Output.Direction)
    case 'Down'
        edge_y = handles.Data.Output.gapinfo.edge_2;        
        y1 = edge_y - view_heigth;
        if edge_y > max_heigth - view_heigth
            y2 = max_heigth;
        else
            y2 = edge_y + view_heigth;
        end
        
    case 'Up'
        edge_y = handles.Data.Output.gapinfo.edge_1;
        if edge_y < view_heigth
            y1 = 1;
        else
            y1 = edge_y-view_heigth;
        end
        y2 = edge_y + view_heigth;
end

dy = y2-y1;
handles.disp_inf.y1 = y1;
handles.disp_inf.y2 = y2;

handles.fig = figure(1);
set(gcf,'position',[100 100 round(1.5*handles.Data.Settings.Video_heigth) 1.5*dy]);
set(gcf,'Units','pixels')
set(gca,'Units','normalized')
set(gca,'Position',[0 0 1 1])
handles.ax = gca;
colormap(handles.ax,'gray')

handles.disp_inf.nframes = handles.Data.Tracker.MetaData.NFrames;


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes verify_touch wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = verify_touch_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

Settings = handles.Data.Settings;
nframes = Settings.Nframes;
disp_inf.nframes = nframes;

handles.quit = 0;

Manual = handles.Data.Manual;
Tracker = handles.Data.Tracker;

Evaluate.total_touch(1:nframes) = NaN;
Evaluate.n_manual_correct(1:nframes) = NaN;
Evaluate.n_tracker_correct(1:nframes) = NaN;


handles.uival.n_total = NaN;
handles.uival.n_correct_man = NaN;
handles.uival.n_correct_mwt = NaN;


current_frame = 1;
frame_flag = 1;

set(handles.n_total, 'String','-')
set(handles.n_correct_man, 'String','-')
set(handles.n_correct_mwt, 'String','-')


f2 = figure(2);
f2.Units = 'points';
f2.Position = [200 200 120 150];
ax2 = axes(f2);
ax2.Units = 'points';
ax2.Position = [0 0 120 150];
cla(ax2)
hold(ax2, 'on')
scatter(ax2, 20, 150, 60,'r','filled','Marker','^')
scatter(ax2, 20, 100, 60,'b','filled','Marker','o')
scatter(ax2, 20, 50, 60, 'g','filled','Marker','square')
ylim(ax2,[0 200])
xlim(ax2,[0 60])
text(ax2, 25, 150,'Manual')
text(ax2, 25, 100,'Tracker')


guidata(hObject, handles)

if exist(handles.SaveName, 'file')
    inpt = inputdlg('Load previous data (Y/N)?');
    if strcmp(inpt, 'Y')
        load(handles.SaveName);
        current_frame = find(~isnan(Evaluate.total_touch),1,'last') + 1;
        if isempty(current_frame)
            current_frame = 1;
        end
    end
end


while frame_flag == 1
    disp_inf.current_frame = current_frame;
    
    if size(Manual.Touch.pt,2) > current_frame 
        if isempty( Manual.Touch.pt{current_frame}) 
            if isempty(Tracker.Touch{current_frame}) || ...
                ~any(find(Tracker.Touch{current_frame}))
                current_frame = current_frame + 1;
                continue
            end
        end
    end
    
    
    
    update_disp_inf(handles, disp_inf);
    
    
    Settings.Current_frame = current_frame;
    frame = LF(Settings);
    snap = frame(handles.disp_inf.y1:handles.disp_inf.y2,:);
    imagesc(handles.ax, snap);
    hold(handles.ax, 'on');
        
    if current_frame <= size(Manual.Touch.pt,2)
        n_manual = size( Manual.Touch.pt{current_frame}, 1);
        for i = 1:n_manual
            pt = Manual.Touch.pt{current_frame}(i,:);
            pt(1) = pt(1)-handles.disp_inf.y1;
            scatter(handles.ax, pt(2), pt(1), 'r','filled','Marker','^')
        end
    end
    
    t_pts = find( Tracker.Touch{current_frame} );
    for i = 1:length(t_pts)
        pt = Tracker.Traces_clean{current_frame}{t_pts(i)}(end,:);
        pt(1) = pt(1)-handles.disp_inf.y1;
        scatter(handles.ax, pt(2), pt(1), 'b','filled','Marker','o')
    end
    
    uiwait(handles.figure1);
        
    
    handles = guidata(hObject); 
    if handles.quit == 1 
        close(f2)
        close(handles.fig)
        close(handles.figure1)
        return
    end
    
   

    Evaluate.total_touch( current_frame ) = handles.uival.n_total;
    Evaluate.n_manual_correct( current_frame) = handles.uival.n_correct_man;
    Evaluate.n_tracker_correct(current_frame) = handles.uival.n_correct_mwt;
    disp(handles.uival)
    handles.uival.n_total = NaN;
    handles.uival.n_correct_man = NaN;
    handles.uival.n_correct_mwt = NaN;
    
    
    guidata(hObject, handles)
    set(handles.n_total, 'String','-')
    set(handles.n_correct_man, 'String','-')
    set(handles.n_correct_mwt, 'String','-')
    
    
    save(handles.SaveName, 'Evaluate')
    
    current_frame = current_frame + 1;

end



save(handles.SaveName, 'Evaluate')









varargout{1} = handles.output;



function update_disp_inf(handles, disp_inf)
disp_infstr = sprintf('Frame %d/%d',disp_inf.current_frame, disp_inf.nframes);
set(handles.disp_text,'String',disp_infstr)




function n_total_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to n_total (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uival.n_total = str2double(get(hObject, 'String'));
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of n_total as text
%        str2double(get(hObject,'String')) returns contents of n_total as a double


% --- Executes during object creation, after setting all properties.
function n_total_CreateFcn(hObject, ~, ~)
% hObject    handle to n_total (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_correct_mwt_Callback(hObject, ~, handles)
% hObject    handle to n_correct_mwt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uival.n_correct_mwt = str2double(get(hObject, 'String'));
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of n_correct_mwt as text
%        str2double(get(hObject,'String')) returns contents of n_correct_mwt as a double


% --- Executes during object creation, after setting all properties.
function n_correct_mwt_CreateFcn(hObject, ~, ~)
% hObject    handle to n_correct_mwt (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in submit.
function submit_Callback(~, ~, handles)
% hObject    handle to submit (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA
clc

%disp_inf(handles.uival)

if isnan(handles.uival.n_total)
    disp('Enter # Touch field')
    return
end

if isnan(handles.uival.n_correct_man)
    disp('Enter # Correct Man field')
    return
end

if isnan(handles.uival.n_correct_mwt)
    disp('Enter # Correct MWT field')
    return
end


uiresume(handles.figure1)




% --- Executes on button press in quit.
function quit_Callback(hObject, ~, handles)
% hObject    handle to quit (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.quit = 1;
uiresume
guidata(hObject, handles)


function n_correct_man_Callback(hObject, ~, handles)
% hObject    handle to n_correct_man (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uival.n_correct_man = str2double(get(hObject, 'String'));
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of n_correct_man as text
%        str2double(get(hObject,'String')) returns contents of n_correct_man as a double


% --- Executes during object creation, after setting all properties.
function n_correct_man_CreateFcn(hObject, ~, ~)
% hObject    handle to n_correct_man (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









function frame = LF(Settings)
% Load a frame from videofile 'Video' (fullpath+extension) at frameposition
% 'framenr'.
Video = Settings.Video;
Width = Settings.Video_width;
Heigth = Settings.Video_heigth;
framenr = Settings.Current_frame;
%%


%{
if ~strcmp(filetype,'.dat')
    frame = read(Video,framenr);
    frame = im2double(frame(:,:,1));
    frame = (frame - min(min(frame))) ./ max(max(frame));
 
    h = fspecial('gaussian',10);
    frame = imfilter(frame,h);
    h = fspecial('laplacian');
    f = imfilter(frame,h);
    frame = frame - f;
else
      
    
    % In the case of .dat files, use hardcoded resolution  
    f = fopen(Video,'r');
    fWidth = 512;
    fHeight = 640;
    fdim = [fWidth, fHeight];
    fseek(f,(framenr-1)*fdim(1)*fdim(2),'bof');
    frame = fread(f,fdim,'*uint8');
    fclose(f);
    frame = im2double(frame);
    frame = (frame - min(min(frame))) ./ max(max(frame));
    h = fspecial('gaussian',10);
    frame = imfilter(frame,h);
    h = fspecial('laplacian');
    f = imfilter(frame,h);
    frame = frame - f;

end
%}


%%
dotposition = find(Video == '.',1,'last');
extension = Video(dotposition:end);

switch(extension)
    
    % Add costum read functions here
    case '.dat'
        f = fopen(Video,'r');
        fdim = [Width, Heigth];
        fseek(f,(framenr-1)*fdim(1)*fdim(2),'bof');
        frame = fread(f,fdim,'*uint8');
        fclose(f);
        frame = im2double(frame);
        
        
    otherwise
        
        if isfield(Settings,'Video_object') % Specified in 1st section in ParameterSetup
            frame = read(Settings.Video_object, framenr);
            frame = rgb2gray(frame(:,:,:));
            frame = im2double(frame);
            
            
        else
            fprintf('Extension not supported: %s#n',extension)
        
        end
end

frame = (frame - min(min(frame))) ./ max(max(frame));
h = fspecial('gaussian',10);
frame = imfilter(frame,h);
h = fspecial('laplacian');
f = imfilter(frame,h);
frame = frame-f;
