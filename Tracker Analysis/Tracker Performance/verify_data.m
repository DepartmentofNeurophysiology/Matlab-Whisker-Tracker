function varargout = verify_data(varargin)
%VERIFY_DATA
%
% USAGE
%
%   Opens a GUI to veryify tracking results, on Manual, MWT and
%   Janelia Data.
%
%   GUI FUNCTIONS
%
%   ON STARTUP - check if video has previously been inspected,
%               ask user to continue previous session
%   figure(1) - display current frame with overlayn trace
%   figure(2) - display current frame for reference
%
%   NO - current trace is marked as false trace
%   YES - current trace is marked as true trace
%   checkbox: tip missing - check if tip of current trace is
%                missing
%
%   checkbox: noisy ending - check if tip of current trace
%                 is 'messy'
%   SKIP FRAME - skip current frame, all trace are marked 
%               as false
%   QUIT - exit GUI
%
%
%   OUTPUT
%   Validity is stored in a binary fashion: true trace -1
%   false trace - 0. For all data a cell about:
%   - trace validity
%   - missing tip
%   - noisy tip
%   is available (ie. Manual, Manual_missing_tip, Manual_noisy_tip).
%   Results are saved in 'FILENAME_Evaluation.mat' in the variable
%   'Eval_traces';
%

%% Begin initialization code - DO NOT EDIT
%codegen

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @verify_data_OpeningFcn, ...
    'gui_OutputFcn',  @verify_data_OutputFcn, ...
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


% --- Executes just before verify_data is made visible.
function verify_data_OpeningFcn(hObject, ~, handles, varargin)
handles.output = hObject;

% Select '_Annotations_Tracker' file
[handles.FileName, handles.PathName] = uigetfile('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\*_compiled.mat');
load_data = load(fullfile( handles.PathName, handles.FileName));
handles.Data = load_data.Annotations;
handles.SaveName = fullfile( handles.PathName, [handles.Data.Settings.FileName(1:end-4) '_Evaluation.mat']);
%handles.Data = load(fullfile('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_Annotations_Tracker.mat'));


current_path = handles.PathName(1);
video_path = handles.Data.Settings.Video(1);
if ~strcmp(video_path, current_path)
    fprintf('Directory as saved in data (%s) is not equal to current', video_path);
    fprintf(' working directory (%s)... - updating saved workign directory', current_path);
    handles.Data.Settings.Video = fullfile(handles.PathName, handles.Data.Settings.FileName);
    handles.Data.Settings.outpath = handles.PathName;
    handles.Data.Settings.PathName = handles.PathName;
end
    
    

% Create axes for display
handles.fig = figure(1);
set(gcf,'position',[100 100 round(1*handles.Data.Settings.Video_heigth) ...
    round(1*handles.Data.Settings.Video_width)]);
set(gcf,'Units','pixels')
set(gca,'Units','normalized')
set(gca,'Position',[0 0 1 1])
handles.ax = gca;
colormap(handles.ax,'gray')

handles.fig2 = figure(2);
set(gcf,'position',[100+round(1*handles.Data.Settings.Video_heigth) 100 round(1*handles.Data.Settings.Video_heigth) ...
    round(1*handles.Data.Settings.Video_width)]);
set(gcf,'Units','pixels')
set(gca,'Units','normalized')
set(gca,'Position',[0 0 1 1])
handles.ax2 = gca;
colormap(handles.ax2,'gray')



% Find number of frames to track
%nframes = handles.Data.Tracker.MetaData.NFrames;
nframes = handles.Data.Settings.Nframes;
ntraces_tracker = 0;
%ntraces_manual = 0;
%ntraces_janelia = 0;

for i = 1:nframes
    ntraces_tracker = ntraces_tracker + size(handles.Data.Tracker.Traces_clean{i},2);
    %ntraces_manual = ntraces_manual + size(handles.Data.Manual.Traces{i},2);
    
    %if i < size(handles.Data.Janelia.Traces_clean, 1)
    %    ntraces_janelia = ntraces_janelia + size(handles.Data.Janelia.Traces_clean{i},2);
    %end
end

handles.disp_struct.nframes = nframes;
handles.disp_struct.ntraces_tracker = ntraces_tracker;
%handles.disp_struct.ntraces_manual = ntraces_manual;
%handles.disp_struct.ntraces_janelia = ntraces_janelia;

handles.state.tip_missing = 0;
handles.state.noisy_tip = 0;

% Update handles structure
guidata(hObject, handles);






function verifydata(hObject, ~, handles, varargin)
%% Initialize
Settings = handles.Data.Settings;
nframes = Settings.Nframes;

%Manual = handles.Data.Manual;
Tracker = handles.Data.Tracker;
%Janelia = handles.Data.Janelia;

frame_flag = 1;

% Load data from earlier inspection
if ~exist(handles.SaveName, 'file')
    current_frame = 1;
    
    disp_struct = handles.disp_struct;
    disp_struct.current_frame = 0;
    
    %tracked.manual = 0;
    %tracked.janelia = 0;
    tracked.tracker = 0;
    
    %Eval_traces.Manual = cell(1, nframes);
    %Eval_traces.Manual_missing_tip = cell(1,nframes);
    %Eval_traces.Manual_noisy_tip = cell(1,nframes);
    Eval_traces.Tracker = cell(1, nframes);
    Eval_traces.Tracker_missing_tip = cell(1,nframes);
    Eval_traces.Tracker_noisy_tip = cell(1, nframes);
    %Eval_traces.Janelia = cell(1, nframes);
    %Eval_traces.Janelia_missing_tip = cell(1, nframes);
    %Eval_traces.Janelia_noisy_tip = cell(1, nframes);
    
else
    inpt = inputdlg('Load previous checked data (Y/N)?');
    
    if strcmp(inpt,'Y')
        load(handles.SaveName)
        notated = zeros(1, Tracker.MetaData.NFrames);
        for i = 1:length(notated)
            if ~isempty(Eval_traces.Tracker{i})  %#ok<*OR2,NODEF>
                notated(i) = 1;
            end
        end
        last_frame = find(notated, 1, 'last');
        if ~isempty(last_frame)
            current_frame = last_frame;
        else
            current_frame = 1;
        end
    else
    current_frame = 1;
    
    disp_struct = handles.disp_struct;
    disp_struct.current_frame = 0;
    
    tracked.manual = 0;
    tracked.janelia = 0;
    tracked.tracker = 0;
    
    %Eval_traces.Manual = cell(1, nframes);
    %Eval_traces.Manual_missing_tip = cell(1,nframes);
    %Eval_traces.Manual_noisy_tip = cell(1,nframes);
    Eval_traces.Tracker = cell(1, nframes);
    Eval_traces.Tracker_missing_tip = cell(1,nframes);
    Eval_traces.Tracker_noisy_tip = cell(1, nframes);
    %Eval_traces.Janelia = cell(1, nframes);
    %Eval_traces.Janelia_missing_tip = cell(1, nframes);
    %Eval_traces.Janelia_noisy_tip = cell(1, nframes);
    end
    
end




update_disp_struct(handles,tracked ,disp_struct)
handles.state.skip_frame = 0;
handles.state.quit = 0;
handles.disp_struct = disp_struct;
handles.tracked = tracked;
handles.Eval_traces = Eval_traces;
guidata(hObject, handles)


%% Do inspection

% only verify janelia if there is manual or MWT data

while frame_flag == 1
    
    
    update_disp_struct(handles,tracked ,disp_struct)    
    disp_struct.current_frame = current_frame;
    Settings.Current_frame = current_frame;
    frame = LoadFrame(Settings);
    imagesc(handles.ax, frame);
    imagesc(handles.ax2, frame);    
    %Eval_traces.Manual{current_frame} = zeros(1, size(Manual.Traces{current_frame},2));
    %Eval_traces.Manual_missing_tip{current_frame} = zeros(1, size(Manual.Traces{current_frame},2));
    %Eval_traces.Manual_noisy_tip{current_frame} = zeros(1, size(Manual.Traces{current_frame},2));
    
    %if isempty(Manual.Traces{current_frame}) && isempty(Tracker.Traces_clean{current_frame})
    %    current_frame = current_frame + 1;
    %    continue    
    %end
    
    if current_frame > size(Tracker.Traces_clean, 1)
        break
    end
    
    if isempty(Tracker.Traces_clean{current_frame})
        current_frame = current_frame + 1;
        continue
    end

    %% Manual traces
    
    %{
    for i = 1: size(Manual.Traces{current_frame},2)
        set(handles.noisy_ending,'Value',0);
        set(handles.tip_missing,'Value',0);
        
        guidata(hObject, handles);      
        
        
        trace = Manual.Traces{current_frame}{i};
        if isempty(trace)
            Eval_traces.Manual{current_frame}(i) = NaN;
            continue
        end
        cla(handles.ax);
        imagesc(handles.ax, frame);
        hold(handles.ax, 'on');
        plot(handles.ax, trace(:,2), trace(:,1), 'r')
        uiwait(handles.figure1);
        state = guidata(hObject);
        
        if state.state.skip_frame
            tracked.manual = tracked.manual + size(Manual.Traces{current_frame},2)-i+1;
            break
        end
        
        switch(state.state.true_whisker)
            case 'yes'
                Eval_traces.Manual{current_frame}(i) = 2;
            case 'no'
                Eval_traces.Manual{current_frame}(i) = 1;
        end
        
        if isfield(state.state,'tip_missing') && state.state.tip_missing == 1
            Eval_traces.Manual_missing_tip{current_frame}(i) = 1;
        else
            Eval_traces.Manual_missing_tip{current_frame}(i) = 0;
        end
        
        if isfield(state.state,'noisy_tip') && state.state.noisy_tip == 1
            Eval_traces.Manual_noisy_tip{current_frame}(i) = 1;
        else
            Eval_traces.Manual_noisy_tip{current_frame}(i) = 0;
        end
        
        %disp(state.state.noisy_tip)
        
        
        tracked.manual = tracked.manual  + 1;
         update_disp_struct(handles, tracked, disp_struct)
    end
    
    %}
    
    %% Tracker Traces
    Eval_traces.Tracker{current_frame} = zeros(1, size(Tracker.Traces_clean{current_frame},2));
    Eval_traces.Tracker_missing_tip{current_frame} = zeros(1, size(Tracker.Traces_clean{current_frame},2));
    Eval_traces.Tracker_noisy_tip{current_frame} = zeros(1, size(Tracker.Traces_clean{current_frame},2));
    
    
    
    for i = 1: size(Tracker.Traces_clean{current_frame},2)
        set(handles.noisy_ending,'Value',0);
        set(handles.tip_missing,'Value',0);
        
        guidata(hObject, handles)
        
        trace = Tracker.Traces_clean{current_frame}{i};
        if isempty(trace)
            Eval_traces.Tracker{current_frame}(i) = NaN;
           
            continue
        end
        cla(handles.ax);
        imagesc(handles.ax, frame);
        hold(handles.ax, 'on');
        plot(handles.ax, trace(:,2), trace(:,1), 'r')
        
        if get(handles.autoscale,'Value')
            range(1,1:2) = min(trace);
            range(2,1:2) = max(trace);
            drange = diff(range,1);
            [~, minid] = max(drange);
            switch(minid)
                case 1
                    w_size(1) = drange(1)+50;
                    w_size(2) = w_size(1)*640/512;
                case 2
                    w_size(2) = drange(2)+50;
                    w_size(1) = w_size(2)*512/640;
            end
            dw = (w_size-drange)./2;
            X1 = range(1,1) - dw(1);
            X2 = range(2,1) + dw(1);
            Y1 = range(1,2) -dw(2);
            Y2 = range(2,2) + dw(2);


            xlim(handles.ax, [Y1 Y2])
            ylim(handles.ax, [X1 X2])   
            
        end
        cla(handles.ax2);
        imagesc(handles.ax2, frame);
        
        if get(handles.autoscale,'Value')
            hold(handles.ax2, 'on');
            plot(handles.ax2, trace(:,2), trace(:,1), 'r')
        end
        
        uiwait(handles.figure1);
        state = guidata(hObject);
        
        if state.state.skip_frame
             tracked.tracker = tracked.tracker + size(Tracker.Traces_clean{current_frame},2)-i+1;
            break
        end
        
        if state.state.quit
            return
        end
        
        switch(state.state.true_whisker)
            case 'yes'
                Eval_traces.Tracker{current_frame}(i) = 2;
            case 'no'
                Eval_traces.Tracker{current_frame}(i) = 1;
        end
        
        if isfield(state.state,'tip_missing') && state.state.tip_missing == 1
            Eval_traces.Tracker_missing_tip{current_frame}(i) = 1;
        else
            Eval_traces.Tracker_missing_tip{current_frame}(i) = 0;
        end
        
        if isfield(state.state,'noisy_tip') && state.state.noisy_tip == 1
            Eval_traces.Tracker_noisy_tip{current_frame}(i) = 1;
        else
            Eval_traces.Tracker_noisy_tip{current_frame}(i) = 0;
        end
        
         %disp(state.state.noisy_tip)
         
        tracked.tracker = tracked.tracker+1;
         update_disp_struct(handles, tracked, disp_struct)
    end
    
    
    %% Janelia traces
    
    %{
    Eval_traces.Janelia{current_frame} = zeros(1, size(Janelia.Traces_clean{current_frame},2));
    Eval_traces.Janelia_missing_tip{current_frame} = zeros(1, size(Janelia.Traces_clean{current_frame},2));
    Eval_traces.Janelia_noisy_tip{current_frame} = zeros(1, size(Janelia.Traces_clean{current_frame},2));    
    
    for i = 1: size(Janelia.Traces_clean{current_frame},2)
        set(handles.noisy_ending,'Value',0);
        set(handles.tip_missing,'Value',0);
        
        guidata(hObject, handles)
        
        trace = Janelia.Traces_clean{current_frame}{i};
        if isempty(trace)
            Eval_traces.Janelia{current_frame}(i) = NaN;
            
            continue
        end
        cla(handles.ax);
        imagesc(handles.ax, frame);
        hold(handles.ax, 'on');
        plot(handles.ax, trace(:,2), trace(:,1), 'r')
        uiwait(handles.figure1);
        state = guidata(hObject);
        
        if state.state.skip_frame
            tracked.janelia = tracked.janelia + size(Janelia.Traces_clean{current_frame}, 2)-i+1;
            break
        end
        if state.state.quit
            return
        end
        
        switch(state.state.true_whisker)
            case 'yes'
                Eval_traces.Janelia{current_frame}(i) = 2;
            case 'no'
                Eval_traces.Janelia{current_frame}(i) = 1;
        end
        
        if isfield(state.state,'tip_missing') && state.state.tip_missing == 1
            Eval_traces.Janelia_missing_tip{current_frame}(i) = 1;
        else
            Eval_traces.Janelia_missing_tip{current_frame}(i) = 0;
        end
        
        if isfield(state.state,'noisy_tip') && state.state.noisy_tip == 1
            Eval_traces.Janelia_noisy_tip{current_frame}(i) = 1;
        else
            Eval_traces.Janelia_noisy_tip{current_frame}(i) = 0;
        end
        
         %disp(state.state.noisy_tip)
        
        tracked.janelia = tracked.janelia+1;
         update_disp_struct(handles, tracked, disp_struct)
    end
    %}
    
    current_frame = current_frame + 1;
    
 
    
    if current_frame > Tracker.MetaData.NFrames
        break
    end
    
    if state.state.quit
        return
    end
    
    if state.state.skip_frame
        handles.state.skip_frame = 0;
        guidata(hObject, handles)
    end
    
    
    handles.Eval_traces = Eval_traces;
    handles.tracked = tracked;
    handles.disp_struct = disp_struct;
    guidata(hObject, handles)
    save(handles.SaveName, 'Eval_traces','tracked','disp_struct')
    
end
  save(handles.SaveName, 'Eval_traces','tracked','disp_struct')


guidata(hObject, handles)






function update_disp_struct(handles,tracked, disp_struct)
    %disp_structstr = sprintf('Frame %d/%d\nman: %d/%d\ntrack: %d/%d\njan: %d/%d',...
    %    disp_struct.current_frame, disp_struct.nframes,...
    %    tracked.manual, disp_struct.ntraces_manual, ...
    %    tracked.tracker, disp_struct.ntraces_tracker, ...
    %    tracked.janelia, disp_struct.ntraces_janelia);
    disp_structstr = sprintf('Frame %d/%d\nMWT: %d/%d', ...
        disp_struct.current_frame, disp_struct.nframes,...
        tracked.tracker, disp_struct.ntraces_tracker);
    set(handles.display,'String',disp_structstr)





% --- Outputs from this function are returned to the command line.
function varargout = verify_data_OutputFcn(hObject, ~, handles)
verifydata(hObject, [], handles, [])
close(handles.fig)
close(handles.fig2)
close(handles.figure1)
varargout{1} = handles.output;




function no_button_Callback(hObject, ~, handles) %#ok<*DEFNU>
uiresume
handles.state.true_whisker = 'no';
guidata(hObject, handles)


function yes_button_Callback(hObject, ~, handles)
handles.state.true_whisker = 'yes';
uiresume
guidata(hObject, handles);



function tip_missing_Callback(hObject, ~, handles)
if get(hObject,'Value')
    handles.state.tip_missing = 1;
else
    handles.state.tip_missing = 0;
    
end
guidata(hObject,handles)




function noisy_ending_Callback(hObject, ~, handles)
if get(hObject,'Value')
    handles.state.noisy_tip = 1;
else
    handles.state_noisy_tip = 0;
end
guidata(hObject,handles)



function skip_frame_Callback(hObject, ~ , handles)
handles.state.skip_frame = 1;
uiresume
guidata(hObject, handles);


function quit_Callback(hObject, ~, handles)
handles.state.quit = 1;
uiresume
guidata(hObject, handles);


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


% --- Executes on button press in autoscale.
function autoscale_Callback(hObject, eventdata, handles)
% hObject    handle to autoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoscale
