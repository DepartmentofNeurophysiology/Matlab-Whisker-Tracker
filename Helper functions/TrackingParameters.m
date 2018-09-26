function varargout = TrackingParameters(varargin)
% TRACKINGPARAMETERS MATLAB code for TrackingParameters.fig
%      TRACKINGPARAMETERS, by itself, creates a new TRACKINGPARAMETERS or raises the existing
%      singleton*.
%
%      H = TRACKINGPARAMETERS returns the handle to a new TRACKINGPARAMETERS or the handle to
%      the existing singleton*.
%
%      TRACKINGPARAMETERS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACKINGPARAMETERS.M with the given input arguments.
%
%      TRACKINGPARAMETERS('Property','Value',...) creates a new TRACKINGPARAMETERS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TrackingParameters_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TrackingParameters_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TrackingParameters

% Last Modified by GUIDE v2.5 07-Jul-2018 23:51:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TrackingParameters_OpeningFcn, ...
                   'gui_OutputFcn',  @TrackingParameters_OutputFcn, ...
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


% --- Executes just before TrackingParameters is made visible.
function TrackingParameters_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to TrackingParameters (see VARARGIN)

% Choose default command line output for TrackingParameters
handles.output = hObject;
handles.mTrackSettings = varargin{1};
handles.Output = varargin{2};
% Update handles structure
guidata(hObject, handles);


set(handles.Dilationtext, 'String' , sprintf('%2.0f',handles.mTrackSettings.Dilationsize))
set(handles.edit2, 'String', sprintf('%1.2f',handles.mTrackSettings.Origin_threshold))
set(handles.stoptr, 'String', sprintf('%1.2f', handles.mTrackSettings.trace_threshold))



showParams(handles)

% UIWAIT makes TrackingParameters wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TrackingParameters_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.mTrackSettings;
close(handles.figure1);

function showParams(handles)
frame = LoadFrame(handles.mTrackSettings);
Output = TrackFrame(handles.mTrackSettings, handles.Output);


Origins = Output.Origins;
Traces = Output.Traces;



imagesc(handles.axes1, frame)
colormap('gray')
axis('off')

hold(handles.axes1, 'on');
scatter(Origins(:,2),  Origins(:,1), 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k')

for i = 1:size(Traces,2)
    plot(Traces{i}(:,2), Traces{i}(:,1),'r')
end



function Dilationtext_Callback(hObject, eventdata, handles)
% hObject    handle to Dilationtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mTrackSettings.Dilationsize = ...
    str2double(get(hObject,'String'));
showParams(handles)
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of Dilationtext as text
%        str2double(get(hObject,'String')) returns contents of Dilationtext as a double


% --- Executes during object creation, after setting all properties.
function Dilationtext_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Dilationtext (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EdgeOverlay.
function EdgeOverlay_Callback(hObject, eventdata, handles)
% hObject    handle to EdgeOverlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
flag = get(hObject, 'Value');
if flag
object_edge = edge(handles.Output.Objects);
[x,y] = find(object_edge);
hold(handles.axes1, 'on')
scatter(y,x,2,'MarkerFaceColor','y','MarkerEdgeColor','y')
else
showParams(handles)
end

% Hint: get(hObject,'Value') returns toggle state of EdgeOverlay


% --- Executes on button press in Confirm.
function Confirm_Callback(hObject, eventdata, handles)
% hObject    handle to Confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mTrackSettings.definite_settings = 1;
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in swtichframe.
function swtichframe_Callback(hObject, eventdata, handles)
% hObject    handle to swtichframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mTrackSettings.definite_settings = 0;
guidata(hObject, handles)
uiresume(handles.figure1);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mTrackSettings.Origin_threshold = str2double(get(hObject,'String'));
showParams(handles)
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function stoptr_Callback(hObject, eventdata, handles)
% hObject    handle to stoptr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mTrackSettings.trace_threshold = ...
    str2double(get(hObject, 'String'));
showParams(handles)
guidata(hObject, handles)
% Hints: get(hObject,'String') returns contents of stoptr as text
%        str2double(get(hObject,'String')) returns contents of stoptr as a double


% --- Executes during object creation, after setting all properties.
function stoptr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stoptr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
