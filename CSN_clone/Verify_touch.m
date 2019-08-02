%% Rat 16

clear
clc

sensitivity = 2;

path = 'G:\HighspeedVideo\16_Oct';
vids = dir(fullfile(path,'*','*.avi'));

for idx = 4:size(vids,1)

    idx
    load(fullfile(vids(idx).folder,'ImgA_compiled.mat'))

    T = Annotations.Tracker.Touch;
    touch = zeros(1, size(T,2));

    for k = 1:size(T,2)
        touch(k) = sum(T{sensitivity,k});
    end

    save(fullfile(vids(idx).folder,'touches'),'touch')
end