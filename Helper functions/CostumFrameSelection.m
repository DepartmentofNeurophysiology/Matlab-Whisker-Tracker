function frame_idx = CostumFrameSelection(Settings, Output)
%%
% Return a binary array of size 1 x nframes, 1 indicates the frame should
% be tracked

%%


switch( Settings.frame_select)
    case 'DEFAULT' 
        frame_idx = ones(1,Settings.Nframes);
        
        
    case 'NOSE_REQUIRED'
        % Only track frames where a nose position is found at least 5px
        % from border
        frame_idx = ones(1,Settings.Nframes);
        
        Nose  = Output.Nose;

        for i = 1:length(frame_idx)
            if Nose(i,1) <= 5 | Nose(i,2) <= 5 | ...
                    Nose(i,1) >= Settings.Video_width-5 | ...
                    Nose(i,2) >= Settings.Video_heigth-5
                
                frame_idx(i) = 0;
            end
        end
        
        
    case 'NOSE_INGAP'
        frame_idx = ones(1,Settings.Nframes);
        Nose = Output.Nose;
        gapinfo = Output.gapinfo;
        
        switch(Output.Direction)
            case 'Down'
                x1 = gapinfo.edge_1;
                x2 = gapinfo.edge_2 + 10;
                
                
            case 'Up'
                x1 = gapinfo.edge_1 - 10;
                x2 = gapinfo.edge_2 ;
                
                
        end
        
        for i = 1: length( frame_idx)
            if isnan(Nose(i,1))
                frame_idx(i) = 0;
                continue
            end
            
            if Nose(i,1) <= x1 | Nose(i,1) >= x2
                frame_idx(i) = 0;
            end
        end
        
        n_cons = 0;
        for i = 1: size(Nose,1)
            if frame_idx(i) == 1
                n_cons = n_cons+1;
            elseif frame_idx(i) == 0
                if n_cons <= Settings.min_consc_frames
                    frame_idx(i-n_cons:i-1) = 0;
                end
            end
        end
        
        
    case 'ManualDataRequired'        
        
        man_file = fullfile(Settings.PathName,[Settings.FileName(1:end-4) '_Annotations.mat']);
        
        if exist(man_file, 'file')
            manual_data = load(man_file);
            frame_idx = zeros(1, Settings.Nframes);
            
            for i = 1:length(frame_idx)
                if i<= size(manual_data.CurvesByFrame,1) && ~isempty(manual_data.CurvesByFrame{i})
                    frame_idx(i) = 1;
                end
            end
            
        else
            disp('manual file does not exist')
            keyboard
        end
        keyboard
        
    case 'use_file'
        load(fullfile(Settings.PathName, 'Selected_frames.mat'));
        for i = 1:size(Output, 2)
            if strcmp(Output(i).Video, Settings.FileName)
                break
            end
        end
        Pairs = Output(i).Pairs;
        frame_idx = zeros(1, Settings.Nframes);
        for i = 1:size(Pairs, 2)
            if any(Pairs{i} > Settings.Nframes)
                disp('wrong pairs for video:')
                disp(Settings.FileName)
                frame_idx(:) = 0;
                return
            end
            frame_idx(Pairs{i}(1):Pairs{i}(2)) = 1;
        end
        
        
end



