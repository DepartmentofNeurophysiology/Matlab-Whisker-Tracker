function Parameters = getParams(Tracker,type)
% In a frame with n traces, find m parameters and store in [n x m] matrix
% per frame
% m1 - x-position origin
% m2 - y-position origin
% m3 - x-position tip global
% m4 - y-position tip global
% m5 - theta global
% m6 - theta corrected (angle base whiser)
% m7 - length
% m8 - curvature
% m9 - x-position origin-corrected
% m10- y-position origin-corrected
% m11- x-position tip -corrected
% m12- y-position tip - corrected
% m13- headangle
% m14- nosex
% m15- nosey
% m16- theta tip
% m17- theta tip corrected
% m18 - theta tip-base

%%

switch(type)
    case 'raw'
        Traces = Tracker.Traces;
        if isfield(Tracker, 'Headvec')
            nframes = min([size(Traces, 1), size(Tracker.Headvec,1)]);
        else
            nframes = size(Traces, 1);
        end
    case 'clean'
        Traces = Tracker.Traces_clean;
        if isfield(Tracker, 'Headvec')
            nframes = min([size(Traces, 2), size(Tracker.Headvec,1)]);
        else
            nframes = size(Traces,2);
        end
end

nparams = 8;
Parameters = cell(1,nframes);

dist_px = 20;

for i = 1:nframes
    
    ntraces = size(Traces{i},2);
    if isfield(Tracker,'Nose')
        headangle  = atan2d(Tracker.Headvec(i,2), Tracker.Headvec(i,1));
        R = [cosd(headangle), -sind(headangle); sind(headangle), cosd(headangle)];
        nose = Tracker.Nose(i,:);
    end
    
    m = [];
   
    if ntraces > 0
        m(1:ntraces,1:nparams) = NaN;
        for j = 1:ntraces
            trace = Traces{i}{j};
            
            if isempty(trace)
                continue
            end
            
            
            % (x,y) position base
            m(j,1) = Traces{i}{j}(3,1); %trace(1,1);
            m(j,2) = Traces{i}{j}(3,2); %trace(1,2);
            
            % (x,y) position tip
            m(j,3) = trace(end,1);
            m(j,4) = trace(end,2);
            
            % Theta
            dl = 0;
            i1 = 3;
            i2 = size(trace,1);
            for k = i1:size(trace,1)
                dl = dl + sqrt( sum( (trace(k,:)-trace(k-1,:)).^2));
                if dl >= dist_px
                    i2 = k;
                    break
                end
            end                        
            
            v_base = trace(i2,:) - trace(i1,:);
            
            m(j,5) = atan2d(v_base(2), v_base(1));
            if m(j,5) > 180
                m(j,5) = m(j,5) - 360;
            elseif m(j,5) < -180
                m(j,5) = m(j,5) + 360;
            end
            
            % Theta corrected
            if isfield(Tracker,'Nose')
                m(j,6) = (headangle+180) - m(j,5);
                if m(j,6) > 180
                    m(j,6) = m(j,6) - 360;
                elseif m(j,6) < -180
                    m(j,6) = m(j,6) + 360;
                end
            else
                m(j,6) = NaN;
            end
            
            % Length
            m(j,7) = size(trace,1);
            
            % Curvature
            v2 = trace(end,:) - trace(end-4,:);
            m(j,8) = abs(m(j,6)-atan2d(v2(2), v2(1)));
            
            % corrected origin position
            if isfield(Tracker,'Nose')
                pt = m(j,1:2) - nose;
                pt = R*pt';
                m(j,9) = pt(1);
                m(j,10) = pt(2);
            
                % corrected tip position
                pt = m(j,3:4) - nose;
                pt = R*pt';
                m(j,11) = pt(1);
                m(j,12) = pt(2);
            else
                m(j,9:12) = NaN;
            end
            
            % headangle
            if isfield(Tracker,'Nose')
                m(j,13) = headangle+180;
                if m(j,13) > 180
                    m(j,13) =  m(j,13) - 360;
                elseif m(j,13) < -180
                    m(j,13) =  m(j,13) + 360;
                end
            else
                m(j,13)= NaN;
            end
            
            % nose
            
            if isfield(Tracker,'Nose')
                m(j,14) = nose(1);
                m(j,15) = nose(2);
            else
                m(j,14:15) = NaN;
            end
            
            if size(trace, 1) > 5
                % theta tip
                 dl = 0;
                i1 = 3;
                i2 = 1;
                
                for k = i1:size(trace, 1)-1                    
                    dl = dl + sqrt( sum( (trace( size(trace,1)+1-k ,:)-trace( size(trace,1)-k ,:)).^2));
                    if dl >= dist_px
                        i2 = k;
                        break
                    end
                end        



                v_tip = trace( size(trace,1) - i1 ,:) - trace(i2,:);
                m(j,16) = atan2d(v_tip(2), v_tip(1));
                if m(j,16) > 180
                    m(j,16) = m(j,16) - 360;
                elseif m(j,16) < -180
                    m(j,16) = m(j,16) + 360;
                end

                % theta tip corrected
                if isfield(Tracker,'Nose')
                    m(j,17) = m(j,16) + (headangle+180);
                    if m(j,17) > 180
                        m(j,17) = m(j,17) - 360;
                    elseif m(j,17) < -180
                        m(j,17) = m(j,17) + 360;
                    end
                else
                    m(j,17) = NaN;
                end


                % theta tip-base
                m(j,18) = m(j,5) - m(j,16);
            else
                m(j,14:18) = NaN;
            end
            
        end
    end
  
    Parameters{i} = m;
end