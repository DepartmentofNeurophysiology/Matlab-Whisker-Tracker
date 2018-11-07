   
clear

load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\Data_Figure_Par_Eval.mat')

Data.Protraction_table.Tspeed = (Data.Protraction_table.Tamplitude*1000)./Data.Protraction_table.Tduration;
Data.Protraction_table.Mspeed = (Data.Protraction_table.Mamplitude*1000)./Data.Protraction_table.Mduration;
clc
printEval('verbose',1);

fprintf('\nProtraction\nAmplitude (n=%d):\n', size(Data.Protraction_table,1))
fprintf('Manual : %2.2f (std - %2.2f)\n', mean(Data.Protraction_table.Mamplitude), std(Data.Protraction_table.Mamplitude))
fprintf('Tracker: %2.2f (std - %2.2f)\n\n', mean(Data.Protraction_table.Tamplitude), std(Data.Protraction_table.Tamplitude))
fprintf('Duration:\n')
fprintf('Manual : %2.2f (std - %2.2f)\n', mean(Data.Protraction_table.Mduration), std(Data.Protraction_table.Mduration))
fprintf('Tracker: %2.2f (std - %2.2f)\n\n', mean(Data.Protraction_table.Tduration), std(Data.Protraction_table.Tduration))
fprintf('Speed:\n')
fprintf('Manual : %2.2f (std - %2.2f)\n', mean(Data.Protraction_table.Mspeed), std(Data.Protraction_table.Mspeed))
fprintf('Tracker: %2.2f (std - %2.2f)\n\n\n', mean(Data.Protraction_table.Tspeed), std(Data.Protraction_table.Tspeed))

Data.Retraction_table.Tspeed = (Data.Retraction_table.Tamplitude*1000)./Data.Retraction_table.Tduration;
Data.Retraction_table.Mspeed = (Data.Retraction_table.Mamplitude*1000)./Data.Retraction_table.Mduration;

fprintf('Retraction\nAmplitude (n=%d):\n', size(Data.Retraction_table,1))
fprintf('Manual : %2.2f (std - %2.2f)\n', mean(Data.Retraction_table.Mamplitude), std(Data.Retraction_table.Mamplitude))
fprintf('Tracker: %2.2f (std - %2.2f)\n\n', mean(Data.Retraction_table.Tamplitude), std(Data.Retraction_table.Tamplitude))
fprintf('Duration:\n')
fprintf('Manual : %2.2f (std - %2.2f)\n', mean(Data.Retraction_table.Mduration), std(Data.Retraction_table.Mduration))
fprintf('Tracker: %2.2f (std - %2.2f)\n\n', mean(Data.Retraction_table.Tduration), std(Data.Retraction_table.Tduration))
fprintf('Speed:\n')
fprintf('Manual : %2.2f (std - %2.2f)\n', mean(Data.Retraction_table.Mspeed), std(Data.Retraction_table.Mspeed))
fprintf('Tracker: %2.2f (std - %2.2f)\n\n', mean(Data.Retraction_table.Tspeed), std(Data.Retraction_table.Tspeed))


f=  'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';

files = dir(fullfile(f, '*compiled.mat'));
clc
for i = 1:size(files,1)
    load(fullfile(files(i).folder, files(i).name))
    n  =0 ;
    for j = 1:size(Annotations.Output.Traces , 1)
        if ~isempty(Annotations.Output.Traces{j})
            n = n+1;
        end
    end
    
    res(i,1) = n;
    res(i,2) = Annotations.Output.ProcessingTime;
    res(i,3) = res(i,1)/res(i,2);
    
    fprintf('count: %4d, dur: %4f, fps: %2.2f\n', res(i,1), res(i,2) , res(i,3))
end


few_frames_idx = find(res(:,1) > 30);
fprintf('\n\n Total number of videos: %d\n', length(few_frames_idx))
fprintf('Total number of frames: %d\n', sum(res(few_frames_idx,1)))
fprintf('Mean processing speed: %2.2f FPS (%2.3f) Mpx/s\n', mean(res(few_frames_idx, 3)), mean(res(few_frames_idx,3))*512*640/10^6)


%%
% TRACKER PERFORMANCE
% Running 100 bootstraps on 35486 samples...
% % Correct traces (from all traces): 0.9735, std: 0.0008
% % Correct no artifact (from all traces): 0.7885, std: 0.0022
% 
% Protraction
% Amplitude (n=437):
% Manual : -17.22 (std - 10.45)
% Tracker: -14.73 (std - 11.28)
% 
% Duration:
% Manual : 28.12 (std - 10.18)
% Tracker: 24.54 (std - 10.78)
% 
% Speed:
% Manual : -593.38 (std - 279.44)
% Tracker: -557.69 (std - 368.51)
% 
% 
% Retraction
% Amplitude (n=471):
% Manual : 18.94 (std - 11.09)
% Tracker: 15.94 (std - 11.37)
% 
% Duration:
% Manual : 24.47 (std - 7.30)
% Tracker: 23.10 (std - 9.81)
% 
% Speed:
% Manual : 738.30 (std - 365.95)
% Tracker: 640.11 (std - 385.72)
