function frame = LoadFrame(Settings)
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

