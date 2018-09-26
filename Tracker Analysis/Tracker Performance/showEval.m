function showEval(varargin)
%showEval
% name-value input:
% 'file' - absolute path to file
% 'type' - inspect error type ('missing_tip','noisy-end','wrong')
%%

p = inputParser;
p.CaseSensitive = 0;

addParameter(p,'File', '');
addParameter(p,'Type','');
parse(p, varargin{:});

if isempty(p.Results.File)
    return
end


file = p.Results.File;
type = p.Results.Type;


impath = 'C:\Users\Thijs\Desktop\temp';
%%

%clear
%clc
%file = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M47_R05_08_Evaluation.mat';

%type = 'wrong';
load(file)

lidx = find(file == '\', 1, 'last');
file2 = [file(1:lidx)  file(lidx+1:end-15)  '_compiled.mat'];
load(file2)

switch(type)
    case 'missing_tip'
        marks = Eval_traces.Tracker_missing_tip;
        
    case 'noisy-end'
        marks = Eval_traces.Tracker_noisy_tip;
        
    case 'wrong'
        marks = Eval_traces.Tracker;
end

S = Annotations.Settings;
S.Video(1) = file(1);
T = Annotations.Tracker.Traces_clean;

  f1 =  figure(1);
for i = 1:size(marks,2)
    if isempty(marks{i})
        continue
    end
    
    S.Current_frame = i;
    switch(type)
        case 'wrong'
            idx = find(marks{i} == 1);
        otherwise
            idx = find(marks{i});
    end
    f = LoadFrame(S);
    for j = 1:length(idx)
        t = T{i}{idx(j)};
        
      
        clf(f1);
        imagesc(f);
        colormap gray
        hold on
        
        plot(t(:,2), t(:,1), 'r')
        
        xlim([min(t(:,2))-15 max(t(:,2))+30])
        ylim([min(t(:,1))-15 max(t(:,1))+30])
        
        
        w = waitforbuttonpress;
        disp(f1.CurrentCharacter)
        if f1.CurrentCharacter == 's'
            files = dir(fullfile(impath, '*.png'));
            newname = sprintf('%d.png', size(files,1)+1);
            frame = getframe;
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            imwrite(imind, cm, fullfile(impath,newname))
            fprintf('%s - saved\n', newname)
        end
        
        
        
    
    end
end



