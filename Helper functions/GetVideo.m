function GetVideo()
[FileName,PathName] = uigetfile('*.dat','MultiSelect','on');
path = pwd;
%%
for j= 1:length(FileName)
    j
    Name = FileName{j}(1:length(FileName{j})-4);
    
    cd(PathName)
    Namemat = [Name '.mat'];
    Namedat = [Name '.dat'];
    load(Namemat) ; 
    fid = fopen(Namedat,'r');
    BytesPerFrame = Data.NFrames*Data.Resolution(1)*Data.Resolution(2);
    Video = zeros([Data.Resolution,1,Data.NFrames],'uint8');
    D = fread(fid,BytesPerFrame,'*uint8');
    cFrames = 1:Data.NFrames;
    Video(:,:,:,cFrames) = reshape(D,Data.Dims); 
    Name = [Name '_Vid.mat'];
    save (Name, 'Video')
    
end

cd(path)
