function PRINT_VIDEO(varargin)
%% 
% Name - Values
% dPath - datapath with videos
% FileIndex - file number as in directory
% FileName  - file name
%
% dNose     - show nose
% dDTT      - distance to target
% 
% dMraw     - show raw manual 
% dMclean   - show clean manual
% dMtouch   - show Manual touch
%
% dTraw     - show raw tracker
% dTclean   - show clean tracker
% dTtouch   - show Tracker touch
% 
% d2A       - use two axes for display
% Max       - 'ax1' or 'ax2', axis for display
% Tax       - 'ax1' or 'ax2', axis for display
%
% FrameSelect - 'full', show all frames,
%               'cut', shw frames in 'FrameRange'
%               'annotated', show frames with tracker annotations
% FrameRange  - [a,b] , show frames in range a,b
%
% dExp      - Export video
% dExpT     - 'avi' or 'gif', export video type
%
% Press 'q' during playback to quit

%%
p = inputParser;

addParameter(p,'dPath', 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance');
addParameter(p,'FileIndex',[])
addParameter(p,'FileName',[])

addParameter(p,'dNose', 0);
addParameter(p,'dDTT',0);

addParameter(p,'dMraw',0);
addParameter(p,'dMclean',0);
addParameter(p,'dMtouch',0);

addParameter(p,'dTraw',0);
addParameter(p,'dTclean',0);
addParameter(p,'dTtouch',0);

addParameter(p,'d2A',0);
addParameter(p,'Max','ax1');
addParameter(p,'Tax','ax1');


addParameter(p,'FrameSelect','full')
addParameter(p,'FrameRange',[1,2])
addParameter(p,'dExp', 0);
addParameter(p,'dExpT', 'avi');
parse(p, varargin{:});

Files = dir(fullfile(p.Results.dPath,'*_compiled.mat'));




%%


if isempty(p.Results.FileIndex) & isempty(p.Results.FileName)
    disp('Select video file... (''FileIndex'')')
    return
end

if ~isempty(p.Results.FileIndex)
    load(fullfile( Files(p.Results.FileIndex).folder, Files(p.Results.FileIndex).name) )
    video_file = fullfile( Files(p.Results.FileIndex).folder, [Files(p.Results.FileIndex).name(1:end-13) '.dat']);
elseif ~isempty(p.Results.FileName)
    fname = [p.Results.FileName '_compiled.mat'];
    load(fullfile( p.Results.dPath, fname))
    video_file = fullfile(p.Results.dPath, [p.Results.FileName '.dat']);
end
disp(video_file)
Settings = Annotations.Settings;
Settings.Video(1) = 'E';



if p.Results.dDTT
    Dist_to_target = Annotations.Tracker.dist_nose_target;
    Nose = Annotations.Tracker.Nose;
    switch(Annotations.Tracker.Direction)
        case 'Up'
            target_heigth = Annotations.Tracker.gapinfo.edge_1;
        case 'Down'
            target_heigth = Annotations.Tracker.gapinfo.edge_2;
    end
end

if p.Results.dNose;    Nose = Annotations.Tracker.Nose; end
if p.Results.dMraw;    Manual_raw = Annotations.Manual.RawNotations; end
if p.Results.dMclean;  Manual_clean = Annotations.Manual.Traces; end
if p.Results.dMtouch;  Manual_touch = Annotations.Manual.Touch; end
if p.Results.dTraw;    Tracker_raw = Annotations.Tracker.Traces; end
if p.Results.dTclean;  Tracker_clean = Annotations.Tracker.Traces_clean; end
if p.Results.dTtouch;  Tracker_touch = Annotations.Tracker.Touch; end


if isfield(Annotations,'Tracker')
    gapwidth = abs(Annotations.Tracker.gapinfo.edge_1 - Annotations.Tracker.gapinfo.edge_2);
    nframes = size(Annotations.Output.Nose, 1);
else
    gapwidth = 200;
    keyboard % implement 'nframes'
end

Colors = makeColor(gapwidth);


display.fig = figure(1);

if ~p.Results.d2A
    display.fig_width = round(1*Annotations.Settings.Video_heigth);
    display.two_ax = 0;
elseif p.Results.d2A
    display.fig_width = 2*round(1*Annotations.Settings.Video_heigth);
    display.two_ax = 1;
end
display.fig_heigth = round(1*Annotations.Settings.Video_width);

set(display.fig, 'position', [150 150 display.fig_width display.fig_heigth]);
set(display.fig,'Units','pixels')

display.ax1 = axes();
set(display.ax1,'Units','normalized')
if ~p.Results.d2A
    set(display.ax1,'Position',[0 0 1 1])
elseif p.Results.d2A
    set(display.ax1,'Position',[0 0 0.5 1])
end
colormap(display.ax1, 'gray')

if p.Results.d2A
    display.ax2 = axes();
    set(display.ax2,'Units','normalized')
    set(display.ax2,'Position',[0.5 0 0.5 1])
    colormap(display.ax2, 'gray');
end


switch(p.Results.FrameSelect)
    case 'full'
        display.show_frames = 1:Annotations.Settings.Nframes;
        
    case 'cut'
        if length(p.Results.FrameRange) == 2
            display.show_frames = p.Results.FrameRange(1):p.Results.FrameRange(2);
        elseif length(p.Results.FrameRange) == 1
            display.show_frames = p.Results.FrameRange;
        end
        
    case 'annotated'
        idx = zeros(1, size(Annotations.Tracker.Traces_clean, 1));
        for ii = 1:length(idx)
            if ~isempty(Annotations.Tracker.Traces_clean{ii})
                idx(ii) = 1;
            end
        end
        id = [find(idx == 1, 1, 'first') find(idx==1, 1, 'last')];  
        display.show_frames = id(1):id(2);
end



if p.Results.dExp
    switch(p.Results.dExpT)
        case 'avi'
            vidname = [Files(p.Results.FileIndex).name(1:end-13) '_Annotated'];
            vidout = VideoWriter(fullfile(Files(p.Results.FileIndex).folder, vidname), 'Motion JPEG AVI');
            open(vidout)
            
        case 'gif'
            vidname = fullfile(Files(p.Results.FileIndex).folder, ...
                [Files(p.Results.FileIndex).name(1:end-13) '_Annotated.gif']);
            
    end
end






%%




for id = 1:length(display.show_frames)
    frame_index = display.show_frames(id);
    
    Settings.Current_frame = frame_index;
    frame = LoadFrame(Settings);
    
    
    cla(display.ax1);
    imagesc(display.ax1, frame);
    hold(display.ax1, 'on')
    axis(display.ax1,'off')
    
    if p.Results.d2A
        cla(display.ax2);
        imagesc(display.ax2, frame);
        hold(display.ax2, 'on');
        axis(display.ax2, 'off')
    end
    
    
    
    
    if p.Results.dNose
        scatter(display.ax1, Nose(frame_index,2), Nose(frame_index,1), ...
            'MarkerFaceColor', Colors.nose, 'MarkerEdgeColor', Colors.nose)
        
        vec = Angle(frame_index,:);
        quiver(display.ax1, Nose(frame_index,2), Nose(frame_index,1), ...
            vec(2)*40, vec(1)*40, 'color',Colors.nose)
        
        if p.Results.d2A
            scatter(display.ax2, Nose(frame_index,2), Nose(frame_index,1), ...
                'MarkerFaceColor', Colors.nose, 'MarkerEdgeColor', Colors.nose)
            
            vec = Angle(frame_index,:);
            quiver(display.ax2, Nose(frame_index,2), Nose(frame_index,1), ...
                vec(2)*40, vec(1)*40, 'color',Colors.nose)
            
        end
        
    end
    
    
    
    
    
    
    
    if p.Results.dDTT
        if Annotations.Tracker.gapinfo.main_ax == 1
            keyboard
        elseif Annotations.Tracker.gapinfo.main_ax == 2
            if ~isnan(Nose(frame_index,2))
                X = Nose(frame_index,2);
                dist = abs(round( Dist_to_target(frame_index)));
                if dist == 0
                    dist = 1;
                end
                
                switch(Annotations.Tracker.Direction)
                    case 'Up'
                        midpoint = target_heigth + abs((Nose(frame_index,1)-target_heigth)/2);
                    case 'Down'
                        midpoint = target_heigth - abs((Nose(frame_index,1)-target_heigth)/2);
                end
                
                
                
                line(display.ax1,[X X], [Nose(frame_index,1) target_heigth],'color',Colors.dist_nose_target(dist,:))
                line(display.ax1,[X-5 X+5], [target_heigth target_heigth], 'color', Colors.dist_nose_target(dist,:))
                line(display.ax1,[X-5 X+5], [Nose(frame_index,1) Nose(frame_index,1)],'color',...
                    Colors.dist_nose_target(dist,:))
                
                text(display.ax1, X+5,midpoint,sprintf('%4d px', dist),'color',Colors.dist_nose_target(dist,:))
                
                if display.two_ax
                    line(display.ax2,[X X], [Nose(frame_index,1) target_heigth],'color',Colors.dist_nose_target(dist,:))
                    line(display.ax2,[X-5 X+5], [target_heigth target_heigth], 'color', Colors.dist_nose_target(dist,:))
                    line(display.ax2,[X-5 X+5], [Nose(frame_index,1) Nose(frame_index,1)],'color',...
                        Colors.dist_nose_target(dist,:))
                    
                    text(display.ax2, X+5,midpoint,sprintf('%4d px', dist),'color',Colors.dist_nose_target(dist,:))
                    
                end
                
            end
        end
    end
    
    
    
    
    
    
    
    
    if p.Results.dMraw
        for j = 1:size(Manual_raw{frame_index},2)
            if strcmp(Manual_raw{frame_index}{j}{4},'track')
                scatter(display.(p.Results.Max),Manual_raw{frame_index}{j}{1},Manual_raw{frame_index}{j}{2},...
                    'MarkerFaceColor',Colors.manual_dark,'MarkerEdgeColor',Colors.manual_dark)
            end
        end
    end
    if p.Results.dMclean
        for j = 1:size(Manual_clean{frame_index}, 2)
            trace = Manual_clean{frame_index}{j};
            if ~isempty(trace)
                plot(display.(p.Results.Max),trace(:,2), trace(:,1),'color',Colors.manual_light)
            end
        end
    end
    if p.Results.dMtouch
        if ~isempty(Manual_touch.pt{frame_index})
            pts = Manual_touch.pt{frame_index};
            scatter(display.(p.Results.Max),pts(:,2),pts(:,1),'MarkerFaceColor',Colors.manual_touch,...
                'MarkerEdgeColor',Colors.manual_touch,'Marker',Colors.manual_touch_style)
        end
    end
    
    
    
    
    
    if p.Results.dTraw
        for j = 1:size(Tracker_raw{frame_index}, 2)
            trace = Tracker_raw{frame_index}{j};
            plot(display.(p.Results.Tax),trace(:,2), trace(:,1),'color',Colors.tracker_dark)
        end
    end
    if p.Results.dTclean
        for j = 1:size(Tracker_clean{frame_index}, 2)
            trace = Tracker_clean{frame_index}{j};
            plot(display.(p.Results.Tax),trace(:,2), trace(:,1),'color',Colors.tracker_light)
        end
    end
    if p.Results.dTtouch
        idx = find(Tracker_touch{frame_index});
        for j = 1:length(idx)
            pt = Tracker_clean{frame_index}{idx(j)}(end,:);
            scatter(display.(p.Results.Tax),pt(2), pt(1),'MarkerFaceColor',Colors.tracker_touch,...
                'MarkerEdgeColor',Colors.tracker_touch,'Marker',Colors.tracker_touch_style)
        end
    end
    
    text(display.ax1,10,10,num2str(frame_index),'color','r','BackgroundColor','k')
    
    drawnow
    
    
    hold(display.ax1,'off')
    
    if p.Results.d2A
        hold(display.ax2,'off')
    end
    
    if p.Results.dExp
        switch(p.Results.dExpT) %#ok<*UNRCH>
            case 'avi'
                frameout = getframe(display.fig);
                writeVideo(vidout, frameout.cdata);
            case 'gif'
                frameout = getframe(display.fig);
                im = frame2im(frameout);
                [imind, cm] = rgb2ind(im, 256);
                if display.show_frames(1)-frame_index == 0
                    imwrite(imind, cm, vidname, 'gif', 'Loopcount', inf, 'DelayTime',0.05)
                    
                else
                    imwrite(imind, cm, vidname, 'gif', 'WriteMode', 'append','DelayTime',0.05)
                end
        end
    end
    
    
    if strcmp(display.fig.CurrentCharacter,'q')
        close(display.fig)
        break
    end
    
    
    
    
end


if p.Results.dExp
    switch(p.Results.dExpT)
        case 'avi'
            close(vidout)
    end
end

if id == length(display.show_frames)
    close(display.fig)
end



















