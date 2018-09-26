function Envelope = getEnvelope(Params)
%%
%clear
%load('D:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_compiled.mat');

%Params = Annotations.Tracker.Parameters_clean;
%%

for i = 1:size(Params,2)
    

    
    if ~isempty(Params{i})
        loop_angles = Params{i}(:,17);
        
        if ~isempty(loop_angles(loop_angles <= 0))
            
            l_idx = find(loop_angles <= 0);
            
            if length(l_idx) > 5
                [~,sort_idx] = sort(loop_angles(l_idx));
                n_left = length(l_idx);
                mid = floor(n_left/2);
                Envelope.left_min(i) = median(loop_angles(l_idx(sort_idx(1:mid))));
                Envelope.left_max(i) = median(loop_angles(l_idx(sort_idx(mid+1:end))));
            else
                Envelope.left_min(i) = min(loop_angles(l_idx));
                Envelope.left_max(i) = max(loop_angles(l_idx));
            end
                
        end
        
        if ~isempty(loop_angles(loop_angles>0))
            r_idx = find(loop_angles > 0);
            
            
            if length(r_idx) > 5
            
            [~, sort_idx] = sort(loop_angles(r_idx));
            n_right = length(r_idx);
            mid = floor(n_right/2);  
            Envelope.right_min(i) = median(loop_angles(r_idx(sort_idx(1:mid))));            
            Envelope.right_max(i) = median(loop_angles(r_idx(sort_idx(mid+1:end))));
            else
                Envelope.right_min(i) = min(loop_angles(r_idx));
                Envelope.right_max(i) = max(loop_angles(r_idx));
            end
        end
    end
    
end

f_size = 5;

Envelope.left_min = medfilt1(Envelope.left_min, f_size);
Envelope.left_max = medfilt1(Envelope.left_max, f_size);
Envelope.right_min = medfilt1(Envelope.right_min, f_size);
Envelope.right_max = medfilt1(Envelope.right_max, f_size);
