function frame = LoadFrame(Settings)
% frame = LoadFrame(Settings)
% Load a video frame using specifications in a Settings struct:
% Settings.Video        - video-filename
% Settings.Video_width  - video resolution (width, optional if resolution is not inherent to video type)
% Settings.Video_heigth - video resolution (heith, optional if resolution is not inherent to video type)
% Settings.Current_frame- frame index to load
%
% Returns double matrix
%%
warning off MATLAB:subscripting:noSubscriptsSpecified

% Get file name
Video = Settings.Video;
% Get frame idx to extract
framenr = Settings.Current_frame;
% Find data type
dotposition = find(Video == '.',1,'last');
extension = Video(dotposition:end);



switch(extension)
    
    % Add costum read functions here
    case '.dat'
        Width = Settings.Video_width;
        Heigth = Settings.Video_heigth;
        f = fopen(Video,'r');
        fdim = [Width, Heigth];
        fseek(f,(framenr-1)*fdim(1)*fdim(2),'bof');
        frame = im2double(fread(f,fdim,'*uint8'));
        fclose(f);
        
        
        
    case '.mat'
        m = matfile(Video);
        fr = m.movf(1, framenr);
        frame = rgb2gray(fr.cdata);
        frame = im2double(frame);
        
        
    otherwise
        Settings.Video_object = VideoReader(Video);
        if isfield(Settings,'Video_object') % Specified in 1st section in ParameterSetup
            frame = read(Settings.Video_object, framenr);
            disp(size(frame))
            %frame = rgb2gray(frame(:,:,:));
            frame = im2double(frame);
            
            
        else
            fprintf('Extension not supported: %s#n',extension)
        
        end
end


