function varargout = ObjectDetection(varargin)
% OBJECTDETECTION MATLAB code for ObjectDetection.fig
%      OBJECTDETECTION, by itself, creates a new OBJECTDETECTION or raises the existing
%      singleton*.
%
%      H = OBJECTDETECTION returns the handle to a new OBJECTDETECTION or the handle to
%      the existing singleton*.
%
%      OBJECTDETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OBJECTDETECTION.M with the given input arguments.
%
%      OBJECTDETECTION('Property','Value',...) creates a new OBJECTDETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ObjectDetection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ObjectDetection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ObjectDetection

% Last Modified by GUIDE v2.5 07-Jul-2018 19:08:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ObjectDetection_OpeningFcn, ...
                   'gui_OutputFcn',  @ObjectDetection_OutputFcn, ...
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


% --- Executes just before ObjectDetection is made visible.
function ObjectDetection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ObjectDetection (see VARARGIN)
% Choose default command line output for ObjectDetection
handles.output = hObject;
handles.mTrackSettings = varargin{1};
Settings = handles.mTrackSettings;

set(handles.checkbox1, 'Value', Settings.costum_background)
h = waitbar(0, 'Loading video samples');
frameidx = round(linspace(1,Settings.Nframes,Settings.n_background_samples));

% Check if enough memory is free for single batch processing
n_bytes = Settings.Video_width*Settings.Video_heigth*length(frameidx)*8;
mem = memory;
if n_bytes > mem.MaxPossibleArrayBytes
     handles.method = 'sequential';
else
    handles.method = 'single-batch';
end

switch(handles.method)
    case 'single-batch'
        SummedFrames = zeros(Settings.Video_width,Settings.Video_heigth,length(frameidx));

        for FrameNR = 1:length(frameidx)   
            Settings.Current_frame = frameidx(FrameNR);
            frame = LoadFrame(Settings); 
            SummedFrames(:,:,FrameNR) = frame; 
            waitbar(FrameNR/length(frameidx))
        end

        for i = 1:size(SummedFrames,1)
            for j = 1:size(SummedFrames,2)
                if numel(find(SummedFrames(i,j,:) > Settings.object_threshold )) > Settings.object_pixel_ratio*size(SummedFrames,3)
                    Objects(i,j) = 0;
                else
                    Objects(i,j) = 1;
                end
            end
            
            waitbar(i/size(SummedFrames,1),h,'Detecting Objects - Thresholding')    
        end
        
        Kl = conv2(mean(SummedFrames,3), ones(10,10)./100, 'same');
        Ks = conv2(mean(SummedFrames,3), ones(6,6)./36, 'same');
        Edges = zeros(size(Objects));
        Edges(Ks./Kl < 0.95) = 1;
        
        close(h)
        
    case 'sequential'
        Objects = zeros(Settings.Video_width, Settings.Video_heigth);
        
        for FrameNR = 1:length(frameidx)
            disp(frameidx(FrameNR))
            Settings.Current_frame = frameidx(FrameNR);
            frame=  LoadFrame(Settings);
            
            [a,b] = find(frame > Settings.object_threshold);
            Objects(a,b) = 1;
            waitbar(FrameNR/length(frameidx),h,'Detecting Objects - Thr seq')
        end
        close(h)
end

if Settings.costum_background
    Objects = Costumbackground( Settings, Objects);
end


imagesc(handles.axes1,Objects)
colormap('gray')
axis('off')

handles.Edges  = Edges;
handles.Objects = Objects;
handles.SummedFrames = SummedFrames;

set(handles.slider1,'Value',Settings.object_threshold);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ObjectDetection wait for user response (see UIRESUME)
if ~Settings.batch_mode
    uiwait(handles.figure1);
end


% --- Outputs from this function are returned to the command line.
function varargout = ObjectDetection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.Objects;
varargout{2} = handles.mTrackSettings.object_threshold;
varargout{3} = handles.Edges;
close(handles.figure1);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = waitbar(0, 'Updating threshold');
handles.mTrackSettings.object_threshold = get(hObject,'Value');
Settings = handles.mTrackSettings;
SummedFrames = handles.SummedFrames;
for i = 1:size(SummedFrames,1)
    for j = 1:size(SummedFrames,2)
        if numel(find(SummedFrames(i,j,:) > Settings.object_threshold )) > Settings.object_pixel_ratio*size(SummedFrames,3)
            Objects(i,j) = 0;
        else
            Objects(i,j) = 1;
        end
    end
    waitbar(i/size(SummedFrames,1),h,'Updating threshold')    
end

if Settings.costum_background
    Objects = Costumbackground( Settings, Objects);
end

imagesc(handles.axes1, Objects)
axis(handles.axes1,'off')
close(h)
handles.Objects = Objects;
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


function pushbutton2_Callback(hObject, eventdata, handles)
edit Costumbackground.m

function checkbox1_Callback(hObject, eventdata, handles)
handles.mTrackSettings.costum_background = get(hObject, 'Value');
guidata(hObject, handles)

