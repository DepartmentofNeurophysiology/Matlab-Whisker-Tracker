path = 'H:\Data\Julien\Highspeed video\12_Oct';

vid = 'Acq_A_015';

PRINT_VIDEO('dPath', fullfile(path, vid), 'FileName', 'ImgA', ...
    'dTraw', 1, 'dTtouch', 1, 'dNose', 1, 'ROI', 1)


%%
clear
load(fullfile(path, vid, 'MWTforPython.mat'))
