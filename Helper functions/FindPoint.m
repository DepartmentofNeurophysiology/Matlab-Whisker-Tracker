function [Trace,TraceState, Settings] = FindPoint(Settings,Shapes,Trace,TraceState)
% [Trace,TraceState] = FindPoint(Data,Settings,Trace,TraceState) finds a 
% next point for intput ‘Trace’. Its direction of tracing is controlled by 
% TraceState (‘ToSnout’ or ‘ToTip’) by prediciting an ROI in the defined 
% direction. Finds a next point using peakdetection and validates the 
% point using ‘ValidatePoint’. If it is a valid point, it is added to Trace, 
% else TraceState is updated to ‘ToTip’ or ‘Finished’.
%%


% Extract processed frame
frame = Shapes.Frame;
Silhouette = Shapes.small;
Objects = Shapes.Objects;
SL = Shapes.large;
Tracked = Shapes.Tracked;

LastPointIDX = size(Trace,1);
NewPointIDX = LastPointIDX + 1;
HalfCircle = Settings.circle_start:1:Settings.circle_end;
Ref = [0 1];


switch(TraceState)
    case 'ToSnout'    
        
        if LastPointIDX == 1            
            tosnoutsetpsize = 3;
            CircleX = ceil(Trace(1,1) +tosnoutsetpsize*sind(1:10:360));
            CircleY = ceil(Trace(1,2) + tosnoutsetpsize*cosd(1:10:360));            
            CircleIDX = sub2ind(size(frame),CircleX,CircleY);           
            
            keep_idx = find(Shapes.large(CircleIDX));
          
            CircleX = CircleX(keep_idx);
            CircleY = CircleY(keep_idx);
            CircleIDX = CircleIDX(keep_idx);
            id = find(SL(CircleIDX) == 1);
            Profile = frame(CircleIDX(id));
            PointIDX = find(Profile == min(Profile),1,'first');           
           if ~isempty(PointIDX)               
            Point(1) = CircleX(id(PointIDX));
            Point(2) = CircleY(id(PointIDX));
            else
                Point = [0,0];
           end
             [flag,~] = ValidatePoint(Settings, Shapes, Point);
            if flag
                Trace(NewPointIDX,1:2) = Point;
            else
                TraceState = 'ToTip';    
                Trace = flip(Trace,1);
            end   
            
            
            
            
            
        elseif LastPointIDX > 1
            Direction = Trace(LastPointIDX,:) - Trace(LastPointIDX-1,:);
            Angle = atan2d(Direction(2),Direction(1)) - ...
                atan2d(Ref(1),Ref(2));
            Theta = HalfCircle - Angle;
            CircleX = ceil(Trace(LastPointIDX,1) + Settings.stepsize...
                *sind(Theta));
            CircleY = ceil(Trace(LastPointIDX,2) + Settings.stepsize...
                *cosd(Theta));            
            CircleIDX = sub2ind(size(frame),CircleX,CircleY);
            Profile = frame(CircleIDX);
            PointIDX = find(Profile == min(Profile),1,'first');
            if ~isempty(PointIDX)               
            Point(1) = CircleX(PointIDX);
            Point(2) = CircleY(PointIDX);
            else
                Point = [0,0];
            end
             [flag,~] = ValidatePoint(Settings, Shapes, Point);
            if flag
                Trace(NewPointIDX,1:2) = Point;
            else
                TraceState = 'ToTip';                
                Trace = flip(Trace,1);
            end           
            
        end
        
    case 'ToTip'
        
      
        if LastPointIDX ==  1
            CircleX = ceil(Trace(1,1) + Settings.stepsize*sind(1:10:360));
            CircleY = ceil(Trace(1,2) + Settings.stepsize*cosd(1:10:360));
            CircleIDX = sub2ind(size(frame),CircleX,CircleY);
            id = find(SL(CircleIDX) == 0);
            Profile = frame(id); %#ok<FNDSB>
            PointIDX = find(Profile == min(Profile),1,'first');
            if isempty(PointIDX)
                Point = [0,0];
            else
                Point(1) = CircleX(id(PointIDX));
                Point(2) = CircleY(id(PointIDX));
            end
             [flag,Settings] = ValidatePoint(Settings, Shapes, Point);
            if flag
                Trace(NewPointIDX,1:2) = Point;
            else                
                TraceState = 'Finished';               
                Trace = flip(Trace,1);
            end
            
            
            
            
            
        else
            Direction = Trace(LastPointIDX,:) - Trace(LastPointIDX - 1,:);
            Angle = atan2d(Direction(2),Direction(1)) - ...
                atan2(Ref(2),Ref(1));
            Theta = HalfCircle - Angle;
            CircleX = ceil(Trace(LastPointIDX,1) + Settings.stepsize...
                *sind(Theta));
            CircleY = ceil(Trace(LastPointIDX,2) + Settings.stepsize...
                *cosd(Theta));
            CircleIDX = sub2ind(size(frame),CircleX,CircleY);
            CircleIDX = unique(CircleIDX);
            Profile = frame(CircleIDX);
            
            [vals, idx] = findpeaks(-Profile);
            
            if ~isempty(idx)
            pointidx = CircleIDX(idx(find(vals == max(vals),1,'first')));
            elseif isempty(idx) & ~isempty(Profile)
                pointidx = CircleIDX(find(Profile == min(Profile),1,'first'));
            else
                pointidx = [];
            end
            
            if ~isempty(pointidx)
            [Point(1),Point(2)] = ind2sub(size(frame), pointidx);
            else
                Point = [1,1];
            end
           
           
            %for ii = 1:length(idx)
            %    [loop_point(1),loop_point(2)] = ind2sub(size(frame), CircleIDX(idx(ii)));
            %    loop_dir = loop_point - Trace(LastPointIDX,:);
            %    loop_angle(ii) = abs( (atan2d(loop_dir(2), loop_dir(1)) - ...
            %        atan2d(Ref(2), Ref(1))) - Angle);
            %end
            
            
            
            
            %PointIDX = find(Profile == min(Profile),1,'first');
            %if ~isempty(PointIDX)               
            %Point(1) = CircleX(PointIDX);
            %Point(2) = CircleY(PointIDX);
            %else
            %    Point = [0,0];
            %end
                       
             [flag,Settings] = ValidatePoint(Settings, Shapes, Point);
            if  flag
                Trace(NewPointIDX,1:2) = Point;

            else
                
                
                
                % try quadratic extrapolation
                if size(Trace, 1) <= 5
                    xdata = Trace(:,1);
                    ydata = Trace(:,2);
                else
                    xdata = Trace(end-5:end,1);
                    ydata = Trace(end-5:end,2);
                end
                
                px = polyfit(1:length(xdata),xdata',2);
                py = polyfit(1:length(ydata),ydata',2);
                fitax = 1:length(xdata)+10; % extrapolate 10 points
                loop_trace(:,1) = round(polyval(px,fitax));
                loop_trace(:,2) = round(polyval(py,fitax));
                
                check = [];
                for i = 1:size(loop_trace,1)
                    [flag, Settings] = ValidatePoint(Settings, Shapes, loop_trace(i,:));
                    if flag
                        check(i) = 1;
                    else
                        check(i) = 0;
                    end
                end
                
                
                newidx = length(xdata) + find(check(length(xdata)+1:end),1,'first');
                
                if ~isempty(newidx)
                    Trace(NewPointIDX,1:2) = loop_trace(newidx,:);
                else
                    TraceState = 'Finished';
                end
                
                
                
                
                
                
                       
            end
        end
        
        
        

      
end


