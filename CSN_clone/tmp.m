vidpath = 'H:\Data\Julien\Highspeed video\12_Oct'; %clone 2
%vidpath = 'G:\HighspeedVideo\13_Oct';



name = 'Acq_A_006';

PRINT_VIDEO('dPath', fullfile(vidpath, name), 'FileName', 'ImgA', ...
            'dTtouch', 1, 'dTclean', 1, 'dNose', 1, 'ROI', 1)
        
        
        %%
        
        
        data = load(fullfile(vidpath,name,'MWTforPython.mat'))
        
        figure(1)
        clf
        imshow(data.Output.Objects)
        
        
        
 %% Tracng notes
 
 %{
 12
 27 - tactile trial, maar hij touched het object niet
 28 - veel noise
 30 - veel noise, nose tracker zit op kabel
 70 - tactile trial, maar hij raakt het object niet
 105 - tt, maar raakt het object niet
 
 13
 
 
 
 
 %}