function cutRoiDataJulien(session)
folders = dir(session);
tracked = zeros(size(folders, 1));
for i = 1:size(folders, 1)
    fprintf('%d, %d\n', i , size(folders,1))
    if folders(i).isdir
        if exist(fullfile(session, folders(i).name, 'MWTforPython.mat'), 'file')
            load(fullfile(session, folders(i).name,'ImgA_Annotations_Tracker.mat'))
            Settings.Video = fullfile(session, folders(i).name, 'ImgA.avi');
            Settings.ExportName = fullfile(session, folders(i).name, 'ImgA_Annotations_Tracker.mat');
            y = Output.ROI(:,2); 
            ObjectsTouch = Output.ObjectsTouch;
            ymean1 = round(mean(y(1:2)));
            ymean2 = round(mean(y(3:4)));
            nrows = ymean2-ymean1;
            quadrantsize = 0.25*nrows;
            ymax = ymean1 + round(quadrantsize);   
            ObjectsTouch(1:ymean1,:) = 0;
            ObjectsTouch(ymax:end,:) = 0;
            Output.ObjectsTouch = ObjectsTouch;
            save( Settings.ExportName,'Output','Settings')        
            compiledata('file',Settings.Video,'data',{'Tracker'},'overwrite',1)
            makeFileForPython( [Settings.Video(1:end-4) '_compiled.mat'])
        end 
    end
end
