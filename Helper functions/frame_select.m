function varargout = frame_select(varargin)
% FRAME_SELECT MATLAB code for frame_select.fig
%      FRAME_SELECT, by itself, creates a new FRAME_SELECT or raises the existing
%      singleton*.
%
%      H = FRAME_SELECT returns the handle to a new FRAME_SELECT or the handle to
%      the existing singleton*.
%
%      FRAME_SELECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRAME_SELECT.M with the given input arguments.
%
%      FRAME_SELECT('Property','Value',...) creates a new FRAME_SELECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before frame_select_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to frame_select_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help frame_select

% Last Modified by GUIDE v2.5 26-Oct-2018 21:02:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @frame_select_OpeningFcn, ...
                   'gui_OutputFcn',  @frame_select_OutputFcn, ...
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


% --- Executes just before frame_select is made visible.
function frame_select_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to frame_select (see VARARGIN)

% Load previously saved settings
if ~exist('Settings\Settings.mat','file')
    makeSettings;
else
    
    % Overwrite default settings with settings from previous session
    load('Settings\Settings.mat')
end

PathName = uigetdir('Select video parent directory');
Settings.batch_dir = PathName;
handles.Settings = Settings;

handles.vid_files = scanfiles(PathName, Settings.video_extension, 10, Settings.format);
handles.PathName = PathName;
handles.default_dir= 1;
set(handles.checkbox2, 'Value', 1);

if exist(fullfile(PathName,'Selected_frames.mat'))
    DataIn = load(fullfile(PathName, 'Selected_frames.mat'));
    handles.Output = DataIn.Output;
end


% Choose default command line output for frame_select
handles.output = hObject;

for i = 1:size(handles.vid_files, 1)
    handles.disp_names.String{i} = handles.vid_files(i).name;
end



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes frame_select wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = frame_select_OutputFcn(hObject, ~, handles) 



if isfield(handles, 'Output')
    file_index = [];
    for i = 1:size(handles.Output,2)
        if isempty(handles.Output(i).Video)
            file_index = i;
            break
        end
    end
    
    if isempty(file_index)
        if size(handles.Output,2) < size(handles.vid_files,1)
            file_index = size(handles.Output,2) + 1;
        else
            file_index = 1;
        end
    end
else    
    file_index = 1;
end


handles.break_flag = 0;
while 1
    
    metafile = fullfile(handles.vid_files(file_index).folder, [handles.vid_files(file_index).name(1:end-4) '.mat']);
    handles.metadata = load(metafile);
    handles.Settings.Video = fullfile(handles.vid_files(file_index).folder, handles.vid_files(file_index).name);
    handles.Settings.Vid_file = handles.vid_files(file_index).name;
    handles.Settings.file_index = file_index;
    handles.Settings.Video_width = handles.metadata.Data.Resolution(1);
    handles.Settings.Video_heigth = handles.metadata.Data.Resolution(2);
    handles.Settings.NFrames = handles.metadata.Data.NFrames;
    handles.Current_direction = handles.default_dir;
    updateBoxes(handles);
    
    guidata(hObject, handles)
    Reset(hObject, handles, 'init');    
    
    showFrame(handles)
    uiwait(handles.figure1);

    handles = guidata(handles.figure1);

    if handles.break_flag
        break
    end
     file_index = handles.next_vid;
end
saveResults(handles.PathName, handles.Output)


varargout{1} = handles.output;
close(handles.figure1)

function updateBoxes(handles)
if handles.Current_direction == 1
    set(handles.checkbox2, 'Value', 1)
else
    set(handles.checkbox2, 'Value', 0)
end

if handles.Current_direction == 2
    set(handles.checkbox3, 'Value', 1)
else
    set(handles.checkbox3, 'Value', 0)
end

if handles.Current_direction == 3
    set(handles.checkbox4, 'Value', 1)
else
    set(handles.checkbox4, 'Value', 0)
end

if handles.Current_direction == 4
    set(handles.checkbox5, 'Value', 1)
else
    set(handles.checkbox5, 'Value', 0)
end



function Reset(hObject, handles, state)
handles.Select.start = [];
handles.Select.end = [];
handles.Select.pairs = {};
handles.Current_frame = 1;
set(handles.slide_select,'Value',1)
set(handles.slide_select,'Min',1)
set(handles.slide_select,'Max', handles.metadata.Data.NFrames)
set(handles.edit_current_frame,'String','1')
handles.pairs = makePairs(handles);
handles.Current_frame = 1;
set(handles.disp_names,'Value',handles.Settings.file_index)

if isfield(handles, 'Output') & size(handles.Output,2) > handles.Settings.file_index ...
        & ~isempty(handles.Output(handles.Settings.file_index).Video) & strcmp(state, 'init')
    
    pairs = handles.Output.Pairs;
    for i = 1:size(pairs,2)
        handles.Select.start(end+1) = pairs{i}(1);
        handles.Select.end(end+1) = pairs{i}(2);
        handles.pairs{i} = pairs{i};
    end
    
    handles.Current_direction = handles.Output.Direction;
    

end

UpdateBar(handles)
updateBoxes(handles)


guidata(hObject, handles)

function UpdateBar(handles)
nframes = handles.metadata.Data.NFrames;
cla(handles.ax_framenrs)
rectangle(handles.ax_framenrs,'Position',[1 0 nframes 1],'FaceColor',[1 1 1],'EdgeColor',[0 0 0])
axis(handles.ax_framenrs, 'off')

hold(handles.ax_framenrs, 'on')
clc
fprintf('PAIRS:\n')

for i = 1:size(handles.pairs,2)
   
    start = handles.pairs{i}(1);
    stop = handles.pairs{i}(2);  
     fprintf('%d\t%d\n',start, stop)
    %rectangle(handles.ax_framenrs,'Position',[start 0 stop 0.5 ],'FaceColor',[0 0 0 0.5])
    line(handles.ax_framenrs,[start start],[0 1],'color','g','LineWidth',2)
    line(handles.ax_framenrs,[stop stop],[0 1],'color','r','LineWidth',2)
end
hold(handles.ax_framenrs, 'off')
xlim(handles.ax_framenrs,[1 nframes])
ylim(handles.ax_framenrs,[0 1])


function Pairs = makePairs(handles)
Pairs = {};
Start = handles.Select.start;
End = handles.Select.end;

if ~isempty(Start) & ~isempty(End)
    for i = 1:length(Start)
        diff = End - Start(i);
        diff(diff<0) = NaN;
        [~, idx] = min(diff);     
        
        
        if any(~isnan(diff)) & ~isempty(idx) & ~isnan(idx)               
            Pairs{end+1} = [Start(i) End(idx)];
        else
            Pairs{end+1} = [Start(i) handles.metadata.Data.NFrames];
        end
    end
    
elseif ~isempty(Start) & isempty(End)
    Pairs{1} = [min(Start) handles.metadata.Data.NFrames];
    
elseif isempty(Start) & ~isempty(End)
    Pairs{1} = [1 max(End)];
    
elseif isempty(Start) & isempty(End)
    Pairs{1} = [1 handles.metadata.Data.NFrames];
end



% --- Executes during object creation, after setting all properties.
function ax_disp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ax_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate ax_disp


% --- Executes on button press in add_start.
function add_start_Callback(hObject, eventdata, handles)
handles.Select.start(end+1) = handles.Current_frame;
handles.pairs = makePairs(handles);
UpdateBar(handles)
guidata(hObject, handles)


% --- Executes on button press in add_end.
function add_end_Callback(hObject, eventdata, handles)
handles.Select.end(end+1) = handles.Current_frame;
handles.pairs = makePairs(handles);
UpdateBar(handles)
guidata(hObject, handles)



function edit_current_frame_Callback(hObject, ~, handles)
handles.Current_frame = str2double(get(hObject,'String'));
set(handles.slide_select,'Value',handles.Current_frame)
showFrame(handles)
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function edit_current_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_current_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in prev_video.
function prev_video_Callback(hObject, eventdata, handles)
% hObject    handle to prev_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
id = handles.Settings.file_index;
handles.Output(id).Video = handles.Settings.Vid_file;
handles.Output(id).Pairs = handles.pairs;
handles.Output(id).frames_to_track = zeros(1, handles.Settings.NFrames);
handles.Output(id).Direction  = handles.Current_direction;
for i = 1:size(handles.pairs,2)
    handles.Output(id).frames_to_track(...
        handles.pairs{i}(1):handles.pairs{i}(2)) = 1;
end
saveResults(handles.PathName, handles.Output)
handles.next_vid = id - 1;
guidata(handles.figure1, handles)
uiresume(handles.figure1)

% --- Executes on button press in next_video.
function next_video_Callback(hObject, ~, handles)
% hObject    handle to next_video (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
id = handles.Settings.file_index;
handles.Output(id).Video = handles.Settings.Vid_file;
handles.Output(id).Pairs = handles.pairs;
handles.Output(id).frames_to_track = zeros(1, handles.Settings.NFrames);
handles.Output(id).Direction = handles.Current_direction;
for i = 1:size(handles.pairs,2)
    handles.Output(id).frames_to_track(...
        handles.pairs{i}(1):handles.pairs{i}(2)) = 1;
end
saveResults(handles.PathName, handles.Output)
handles.next_vid = id + 1;
guidata(handles.figure1, handles)
uiresume(handles.figure1)





% --- Executes on button press in remove.
function remove_Callback(hObject, ~, handles)
Reset(hObject, handles, 'new');


function files = scanfiles(path,extension, nlayers, format)
files = [];

for i = 1:nlayers
    str = path;
    
    for j = 1:i-1
        str = [str '\*'];
    end
    str = [str '\*' extension];
    
    loopfiles = dir(str);
    disp(loopfiles)
    if ~isempty( format )
        for j = 1:size(loopfiles,1)
            token = regexp(loopfiles(j).name,format,'names');
            if ~isempty(token)
               files = [files; loopfiles(j)];
                
            end
        end
    else
        files = [files; loopfiles];
    end
end


% --- Executes on slider movement.
function slide_select_Callback(hObject, ~, handles) %#ok<*DEFNU>
handles.Current_frame = round(get(hObject, 'Value'));
set(handles.edit_current_frame,'String',num2str(handles.Current_frame))
showFrame(handles)
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slide_select_CreateFcn(hObject, ~, ~)
% hObject    handle to slide_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function showFrame(handles)
if isfield(handles, 'Current_frame')
handles.Settings.Current_frame = handles.Current_frame;
else
    handles.Settings.Current_frame = 1;
end
frame = LoadFrame(handles.Settings);
colormap(handles.ax_disp, 'gray')
imagesc(handles.ax_disp, frame)
axis(handles.ax_disp, 'off')


function saveResults(path, Output) %#ok<INUSD>
save(fullfile(path,'Selected_frames'), 'Output')

% --- Executes on button press in quit.
function quit_Callback(hObject, ~, handles)
handles.break_flag = 1;
guidata(hObject, handles)
uiresume(handles.figure1)


% --- Executes on selection change in disp_names.
function disp_names_Callback(hObject, ~, handles)
id = handles.Settings.file_index;
handles.Output(id).Video = handles.Settings.Vid_file;
handles.Output(id).Pairs = handles.pairs;
handles.Output(id).frames_to_track = zeros(1, handles.Settings.NFrames);
for i = 1:size(handles.pairs,2)
    handles.Output(id).frames_to_track(...
        handles.pairs{i}(1):handles.pairs{i}(2)) = 1;
end
handles.next_vid = get(hObject,'Value');
guidata(hObject, handles)
uiresume(handles.figure1)


% --- Executes during object creation, after setting all properties.
function disp_names_CreateFcn(hObject, ~, ~)
% hObject    handle to disp_names (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    handles.Current_direction = 1;
    updateBoxes(handles);
    guidata(hObject, handles)
end
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    handles.Current_direction = 2;
    updateBoxes(handles);
    guidata(hObject, handles);
end

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    handles.Current_direction = 3;
    updateBoxes(handles);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') == 1
    handles.Current_direction = 4;
    updateBoxes(handles);
    guidata(hObject, handles);
end
% Hint: get(hObject,'Value') returns toggle state of checkbox5
