function [flag, Settings] = ValidatePoint(Settings,Shapes,Point)
%Flag = ValidatePoint(Data,Settings,Point) validates returns 1 if Point is
% a valid point or 0 If not, based on the requirements that the point is not:
%  -	within a small distance from border
%  -	on an object
%  -	on the mouse
%  -	indistinguishable from noise
%  -	tracked before 
%%
Objects = Shapes.Objects;
SilhouetteSmall = Shapes.normal;
Frame = Shapes.Frame;
Tracked = Shapes.Tracked;

 max_height = Settings.Video_width;
max_width = Settings.Video_heigth;

% Controls if Point is a valid point on any whisker.
flag = 1; % Valid point by default



% Stop conditions for tracking
Error.border = 0;
Error.object = 0;
Error.mouse = 0;
Error.noise = 0;
Error.tracked = 0; 




if Settings.dist_from_edge < Settings.stepsize
    Settings.dist_from_edge = Settings.stepsize;
end



%% The point is outside a specific range from border
if Point(1) <= Settings.dist_from_edge+1|| ...
        Point(1) >= size(Objects,1) - Settings.dist_from_edge-1 || ...
        Point(2) <= Settings.dist_from_edge+1 || ...
        Point(2) >= size(Objects,2) - Settings.dist_from_edge-1
    flag = 0;
    Error.border = 1;
end

%% The point lies on an object
if flag == 1 && Objects(Point(1),Point(2)) 
    if ~isfield(Settings, 'object_tick') || Settings.object_tick == 0
        Settings.object_tick = 1;
    end
    
    
    Settings.object_tick = 1;
    flag = 0;
    Error.object = 1;
        
        
     
    
   
end

%% The point lies on the mouse
if flag == 1 && SilhouetteSmall(Point(1),Point(2))
    flag = 0;
    Error.mouse = 1;    
end

%% The point is indistinguishable from noise


if flag == 1 
    
    %%
    if isfield(Settings,'current_length') && Settings.current_length < 15
        gridsize = 4;
    else
        
    gridsize = 4;
    end
    
    X1 = Point(1) - gridsize; if X1 < 1; X1 = 1;end
    X2 = Point(1) + gridsize; if X2 > size(Objects,1); X2 = max_height;end
    
    Y1 = Point(2) - gridsize; if Y1 < 1; Y1 = 1;end
    Y2 = Point(2) + gridsize; if Y2 > size(Objects,2); Y2 = max_width;end
    
    GridX = X1:X2;
    GridY = Y1:Y2;
    
    
    [mGridX,mGridY] = meshgrid(GridX,GridY);
    mGridX = mGridX(:);
    mGridY = mGridY(:);
    
    GridIDX = sub2ind(size(Objects),mGridX,mGridY);
    GridIDX = GridIDX(Objects(GridIDX) == 0);
    GridIDX = GridIDX(SilhouetteSmall(GridIDX) == 0);
    
    toggle = 1;
    switch(toggle)
        case 1  
          
            X1 = Point(1) - 1; if X1 < 1; X1 = 1;end
            X2 = Point(1) + 1; if X2 > max_height; X2 = max_height;end
            
            Y1 = Point(2) - 1; if Y1 < 1; Y1 = 1;end
            Y2 = Point(2) + 1; if Y2 > max_width; Y2 = max_width;end
            X = X1:X2;
            Y = Y1:Y2;
            
            
            
            [D1,D2] = meshgrid(X,Y);
            DiaID = [D1(:),D2(:)];
            dind = sub2ind(size(Frame),DiaID(:,1),DiaID(:,2));
            dind = dind(Objects(dind) == 0);
            
            
        case 2
            % alternate diamond
            X1 = Point(1) - 1; if X1 < 1; X1 = 1;end
            X2 = Point(1) + 1; if X2 > max_height; X2 = max_height;end
            
            Y1 = Point(2) - 1; if Y1 < 1; Y1 = 1;end
            Y2 = Point(2) + 1; if Y2 > max_width; Y2 = max_width;end
            [mD1,mD2] = meshgrid(X1:X2,Y1:Y2);
            mD1 = mD1(:);
            mD2 = mD2(:);
            Diamond(:,1) = mD1;
            Diamond(:,2) = mD2;
        
        case 3
            Diamond = Point;
    end
    
    
    
  
    
    MeanIntensity = sum(Frame(GridIDX))/size(GridIDX,1);
    DiamondIntensity = sum(Frame(dind))/size(dind,1);
    
    IntensityRatio = DiamondIntensity/MeanIntensity;
    
    if IntensityRatio > Settings.trace_threshold
        flag = 0;
        Error.noise = 1;
    end
  
end


%% The point has allready been tracked
if flag == 1 && Tracked(Point(1),Point(2)) == 1
 
  flag = 0;
  Error.tracked = 1;
  Settings.prev_tracked = 1;

end










%% Print stop reason

print_reason = 0;

if print_reason == 1 && flag == 0
    fprintf('Stopped due to -- ')
    if Error.border
        fprintf('border ')
    end
    if Error.object
        fprintf('object ')
    end
    
    if Error.mouse
        fprintf('mouse ')
    end
    
    if Error.noise
        fprintf('noise ')
    end
    
    if Error.tracked
        fprintf('tracked ')
    end
    
    fprintf('\n')
end


Settings.stop = Error;





























