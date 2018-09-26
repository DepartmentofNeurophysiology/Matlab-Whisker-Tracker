function varargout = FrameSelection(varargin)
% FRAMESELECTION MATLAB code for FrameSelection.fig
%      FRAMESELECTION, by itself, creates a new FRAMESELECTION or raises the existing
%      singleton*.
%
%      H = FRAMESELECTION returns the handle to a new FRAMESELECTION or the handle to
%      the existing singleton*.
%
%      FRAMESELECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FRAMESELECTION.M with the given input arguments.
%
%      FRAMESELECTION('Property','Value',...) creates a new FRAMESELECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FrameSelection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FrameSelection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FrameSelection

% Last Modified by GUIDE v2.5 07-Jul-2018 21:59:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FrameSelection_OpeningFcn, ...
                   'gui_OutputFcn',  @FrameSelection_OutputFcn, ...
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


% --- Executes just before FrameSelection is made visible.
function FrameSelection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FrameSelection (see VARARGIN)
% Choose default command line output for FrameSelection
handles.output = hObject;
handles.mTrackSettings = varargin{1};
Settings = handles.mTrackSettings;
nframes = Settings.Nframes;
Settings.Current_frame = round(nframes/2);
frame = LoadFrame(Settings);


imagesc(handles.axes1, frame)
colormap('gray')
axis(handles.axes1, 'off');
set(handles.slider1,'Max',nframes);
set(handles.slider1,'Value',Settings.Current_frame);

handles.mTrackSettings = Settings;
guidata(hObject, handles);

% UIWAIT makes FrameSelection wait for user response (see UIRESUME)
 uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FrameSelection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.mTrackSettings;
close(handles.figure1);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Settings = handles.mTrackSettings;
Settings.Current_frame = round(get(hObject, 'Value'));
frame = LoadFrame(Settings);

imagesc(handles.axes1, frame)
colormap('gray')
axis(handles.axes1, 'off');

handles.mTrackSettings = Settings;
guidata(hObject, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
