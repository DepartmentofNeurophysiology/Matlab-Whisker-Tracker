function [Touch, Cross] = DetectTouch(Tracker, vidSettings)
%%
Settings = makeAnalyseSettings();

%fprintf(['ASSUMING: \n\t-no objects within gap, gap is determined by intensity' ...
%    ' profile along y-axis at border of frame\n'])
%fprintf('\t-Only touching within gap (touching above platform is neglected)\n')




Objects = Tracker.Objects;
if isfield(Tracker,'Traces_clean')
Traces = Tracker.Traces_clean;
else
Traces = Tracker.Traces';
end
ncols = 5;
nframes = max( size( Tracker.Traces ));


% Set only target area as detectable
if Settings.target_area_only
    switch(Tracker.Direction)
        case 'Up'
            edge = Tracker.gapinfo.edge_1;
            
        case 'Down'
            edge = Tracker.gapinfo.edge_2;
    end
    
    mask = zeros(size(Objects));
    mask(edge-Settings.objects_nrows:edge+Settings.objects_nrows,:) = 1;
    Objects = Objects.*mask;
    
end


opts = [];
[opts(:,1), opts(:,2) ] = find(Objects);


Touch = cell(1,nframes);
Cross = cell(1,nframes);


h= waitbar(0,'Detecting touch');

for i = 1:nframes
    looptouch = zeros(1, size(Traces{i},2));
    crossflag = zeros(1, size(Traces{i},2));
    
    for j = 1:size(Traces{i},2)
        t = Traces{i}{j};
        if isempty(t)
            continue
        end
        dist = [];
        dist(:,:,1) = opts(:,1) - t(end,1)';
        dist(:,:,2) = opts(:,2) - t(end,2)';
        dist = dist.^2;
        dist = sqrt( sum( dist,3));
        [~, tpt] = find(dist <= Settings.dist_object);
        
        if ~isempty(tpt)
            
            vidSettings.Current_frame = i;
            f = LoadFrame(vidSettings);
            
            looptouch(j) = 1;  
            
            touchpt = round(t(end,:));
            snapsize = 10;

            x1 = touchpt(1)-snapsize;
            if x1 < 1; x1 =1; end
            x2 = touchpt(1)+snapsize;
            if x2 > size(f,1); x2 = size(f,1); end
            y1 = touchpt(2)-snapsize;
            if y1 < 1; y1 =1; end
            y2 = touchpt(2)+snapsize;
            if y2 > size(f,2); y2 = size(f,2); end
            
            snap = f(x1:x2, y1:y2);
            
            osnap = Tracker.Objects(x1:x2, y1:y2);
            oshadow = sum(osnap,2)./size(osnap,2);
            switch(Tracker.Direction)
                case 'Up'
                    idx = find(oshadow > 0.9, 1,'first')-3;
                case 'Down'
                    idx = find(oshadow > 0.9, 1,'last')+3;
            end
            
            
            if idx > size(snap,1) | idx < 1
                crossflag = 0;
                continue
            end
                 
            d = snap(idx,:);
            d = d - median(d);
            
            nhigh = numel(find(abs(d) > 0.04));
            if sum(d)/length(d) > 0.01
                crossflag(j)= 1;                
             
                
            else
                crossflag(j) = 0;
            end
            
            
        end
   
        
      
        
    end
    
    Touch{i} = looptouch;
    Cross{i} = crossflag;

    waitbar( i/nframes);

end
close(h);


