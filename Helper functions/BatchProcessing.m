function varargout = BatchProcessing(varargin)
% BATCHPROCESSING MATLAB code for BatchProcessing.fig
%      BATCHPROCESSING, by itself, creates a new BATCHPROCESSING or raises the existing
%      singleton*.
%
%      H = BATCHPROCESSING returns the handle to a new BATCHPROCESSING or the handle to
%      the existing singleton*.
%
%      BATCHPROCESSING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATCHPROCESSING.M with the given input arguments.
%
%      BATCHPROCESSING('Property','Value',...) creates a new BATCHPROCESSING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BatchProcessing_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BatchProcessing_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BatchProcessing

% Last Modified by GUIDE v2.5 12-Jul-2018 12:12:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BatchProcessing_OpeningFcn, ...
                   'gui_OutputFcn',  @BatchProcessing_OutputFcn, ...
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


% --- Executes just before BatchProcessing is made visible.
function BatchProcessing_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BatchProcessing (see VARARGIN)



% Load previously saved settings
if ~exist('Settings\Settings.mat','file')
    makeSettings;
else
    
    % Overwrite default settings with settings from previous session
    load('Settings\Settings.mat')
    ot = Settings.object_threshold;
    dl = Settings.Dilationsize;
    ort = Settings.Origin_threshold;
    tt = Settings.trace_threshold;
    
    makeSettings;
    Settings.object_threshold = ot;
    Settings.Dilationsze = dl;
    Settings.Origin_threshold = ort;
    Settings.trace_threshold = tt;
end

PathName = uigetdir('Select video parent directory');
%PathName = 'E:\Studie\Stage Neurobiologie\Videos\Mouse 47';
Settings.batch_dir = PathName;
handles.Settings = Settings;

% Scan for files
vid_files = scanfiles(PathName, Settings.video_extension, 10, Settings.format);
track_files = scanfiles(PathName, '_Annotations_Tracker.mat', 10, '');


for i = 1:length(vid_files)
    disp_names{i} = fullfile(vid_files(i).folder, vid_files(i).name); %#ok<*AGROW>
    disp_names{i} = disp_names{i}(length(PathName)+1:end-4);
end

for i = 1:size(track_files,1)
    track_files(i).vid = [track_files(i).name];
end

handles.listbox1.String = {};
handles.listbox2.String = {};

for i = 1:length(vid_files)
    tracked = 0;    
    for j = 1:size(track_files,1)
        t= fullfile(track_files(j).folder, track_files(j).name);
        v = fullfile(vid_files(i).folder, [vid_files(i).name(1:end-4) '_Annotations_Tracker.mat']);
        
        if strcmp(t,v)
            tracked = 1;
            
        end
    end
    
    
    if tracked
        %idx = length(handles.listbox2.String)+1;
        disp_names{i} = [disp_names{i} ' !'];
        vid_files(i).tracked = 1;
        handles.listbox2.String{end+1} = disp_names{i};
        
    elseif ~tracked
        %idx = length(handles.listbox1.String)+1;
        vid_files(i).tracke = 0;
        handles.listbox1.String{end+1} = disp_names{i};
    end
end

handles.listbox1.Max = 100;
handles.listbox2.Max = 100;

% Choose default command line output for BatchProcessing
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);


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
        files = [files; loopfiles];
    end    
end


% UIWAIT makes BatchProcessing wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BatchProcessing_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.Settings.batch_dir;

string_out = handles.listbox1.String;
for i = 1:size(string_out,1)
    if strcmp(string_out{i}(end),'!')
        string_out{i} = string_out{i}(1:end-2);
    end
end


varargout{2} = string_out;
varargout{3} = handles.Settings.video_extension;
close(handles.figure1)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, ~, handles) %#ok<*DEFNU>
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.selected_to_track = get(hObject, 'Value');
guidata(hObject, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, ~, ~)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, ~, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.selected_not_track = get(hObject, 'Value');
guidata(hObject, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns listbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox2


% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, ~, ~)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, ~, handles)
if isfield(handles,'selected_to_track')
stringout = handles.listbox1.String(handles.selected_to_track);
handles.listbox2.String(end+1:end+size(stringout,1)) = stringout;

str = handles.listbox1.String;

new_str = [];
for i = 1:size(str,1)
    if ~ismember(i, handles.selected_to_track)
        new_str{end+1} = str{i};
    end
end





handles.listbox1.String = str2;
end
guidata(hObject, handles)



% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, ~, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'selected_not_track')
stringout = handles.listbox2.String(handles.selected_not_track);
handles.listbox1.String(end+1:end+size(stringout,1)) = stringout;
str = handles.listbox2.String;

new_str = [];
for i = 1:size(str,1)
    if ~ismember(i, handles.selected_not_track)
        new_str{end+1} = str{i};
    end
end


handles.listbox2.String = new_str;
end

guidata(hObject, handles)


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, ~, handles)
stringout = handles.listbox1.String;
handles.listbox2.String(end+1:end+size(stringout,1)) = stringout;
handles.listbox1.String = {};
guidata(hObject,handles)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(~, ~, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
