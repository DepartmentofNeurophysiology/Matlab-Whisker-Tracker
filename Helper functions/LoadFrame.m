function frame = LoadFrame(Settings)
% Load a frame from videofile 'Video' (fullpath+extension) at frameposition
% 'framenr'.
Video = Settings.Video;
Width = Settings.Video_width;
Heigth = Settings.Video_heigth;
framenr = Settings.Current_frame;


%%
dotposition = find(Video == '.',1,'last');
extension = Video(dotposition:end);

switch(extension)
    
    % Add costum read functions here
    case '.dat'
        f = fopen(Video,'r');
        fdim = [Width, Heigth];
        fseek(f,(framenr-1)*fdim(1)*fdim(2),'bof');
        frame = im2double(fread(f,fdim,'*uint8'));
        fclose(f);
        %frame = im2double(frame);
        
        
    case '.mat'
        m = matfile(Video);
        fr = m.movf(1, framenr);
        frame = rgb2gray(fr.cdata);
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


