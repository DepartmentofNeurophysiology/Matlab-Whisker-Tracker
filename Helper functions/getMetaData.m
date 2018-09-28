function Settings = getMetaData(Settings)
%% Settings = getMetaData(Settings)
% Update Settings struct by adding metadata:
% - PathName
% - FileName
% - batch_mode
% - outpath
% - Video_width
% - Video_heigth
% - Nframes
% - Video_object (optional)
%%

Video = Settings.Video;

slash_idx = find(Video == '\',1,'last');
Settings.PathName = Video(1:slash_idx-1);
Settings.FileName = Video(slash_idx+1:end);
Settings.batch_mode = 1;
Settings.outpath = Settings.PathName;

pidx = find(Video == '.');
extension = Video(pidx+1:end);


switch extension
    case 'dat'
        m_file = Video;
        m_file(end-2) = 'm';
        load(m_file)
        
        Settings.Video_width = Data.Resolution(1);
        Settings.Video_heigth = Data.Resolution(2);
        Settings.Nframes = Data.NFrames;
                
    case 'mat'
        m = matfile(Video);
        data = m.movf(1,1);
    
        Settings.Video_width = size(data.cdata,1);
        Settings.Video_heigth = size(data.cdata,2);
        Settings.Nframes = size(m.movf,2);
        Settings.Video_object = m;
        
    otherwise
        vid = VideoReader(Video);        
    
        Settings.Video_width = vid.Height;
        Settings.Video_heigth = vid.Height;
        Settings.Nframes = floor(vid.Duration*vid.FrameRate);
        
        
        
end