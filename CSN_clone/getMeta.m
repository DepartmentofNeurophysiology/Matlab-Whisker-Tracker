path = 'H:\Data\Julien\Highspeed video\22_Oct';

files = dir(path);
for i = 1:length(files)
    disp([i length(files)])
    if isdir(fullfile(path,files(i).name))
        fs = dir(fullfile(path, files(i).name));
        for j = 1:length(fs)
            if strcmp(fs(j).name,'_config.xsv')
                fname = fullfile(fs(j).folder, fs(j).name);
             
                
                fid = fopen(fname);
                line = fgets(fid);

                while ischar(line)    
                    line = fgetl(fid);
                    if strfind(line,'=')
                        res = strsplit(line,'=');
                        val = strsplit(res{2},'\n');
                        try 
                        Meta.(res{1}) = val{1};
                        end
                    end
                end


                save(fullfile(fs(j).folder,'Meta'),'Meta')
              
            end
        end
    end
end
fclose('all')
%%



