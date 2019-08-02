path = 'G:\HighspeedVideo\16_Oct\';

files = dir(path);
for i = 1:length(files)
    i
    if isdir(fullfile(path,files(i).name))
        fs = dir(fullfile(path, files(i).name));
        for j = 1:length(fs)
            if strcmp(fs(j).name,'ImgA_compiled.mat')
                load(fullfile(fs(j).folder, fs(j).name))
                
                T = Annotations.Tracker.Touch;

                touch = zeros(1, size(T,2));
                sensidx = 4;


                for k = 1:size(T,2)
                    touch(k) = sum(T{sensidx,k});
                end
                
                save(fullfile(fs(j).folder,'touches'), 'touch')
            end
        end
    end
end

%%
file = 'ImgA_compiled.mat';




