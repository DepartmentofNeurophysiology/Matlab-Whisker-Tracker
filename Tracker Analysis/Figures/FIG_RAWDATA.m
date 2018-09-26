
clc

if exist('Tracker','var')
    printtracker = 1;
else
    printtracker = 0;
end

if exist('Manual','var')
    printmanual = 1;
else
    printmanual = 0;
end

%x1_tracker = find(General.keep, 1, 'first');
%x2_tracker = find(General.keep, 1, 'last');

%x1_manual = find( General.manual_keep == 1,1,'first');
%x2_manual = find( General.manual_keep == 1,1,'last');



warning('off')
cmap1 = cbrewer('div','RdYlBu',11);
cmap2 = cbrewer('div','RdYlGn',11);
cmap3 = cbrewer('div','PiYG',11);
cmap4 = cbrewer('div','BrBG',11);
cmap5 = cbrewer('seq','YlOrBr',11);

hx = size(Tracker.Objects,1);
hy = size(Tracker.Objects,2);

rawcolor = [0 0 0];
markersize = 1;

%{
% Generate colormap
if length(General.tracker_labels) <= 6
    cc(1,:) = cmap1(2,:);
    cc(2,:) = cmap1(10,:);
    cc(3,:) = cmap2(10,:);
    cc(4,:) = cmap3(2,:);
    cc(5,:) = cmap5(7,:);
    cc(6,:) = cmap4(2,:);
else
    keyboard
end
%}

x1_tracker = find(General.Tracker_ax, 1, 'first');
x2_tracker = find(General.Tracker_ax, 1, 'last');

% Setup figure
f1 = figure('Units','pixels','Position',[100 100 1400 900]);

ny = 3;
nx = 4;

margy = 0.05;
margx = 0.05;

spacey = 0.05;
spacex = 0.05;

widthx = (1-2*margx-(nx-1)*spacex)/nx;
heighty = (1-2*margy-(ny-1)*spacey)/ny;


xtick = 1;
ytick = 1;
for i = 1:ny*nx
    
    xpos = margx + (xtick-1)*(widthx+spacex);
    ypos = 1-margy-heighty -(ytick-1)*(heighty+spacey);
    
    ax{i} = axes('Units','normalized','Position',[xpos, ypos, widthx, heighty]);
    
    xtick = xtick+1;
    
    if xtick > nx
        xtick = 1;
        ytick = ytick+1;
    end
    
    hold(ax{i},'on')
    xlim(ax{i},[x1_tracker x2_tracker])
    
    
end


% Plot raw parameters

if printtracker
    pos_xo = [];
    pos_yo = [];
    pos_xoc = [];
    pos_yoc = [];
    pos_xt = [];
    pos_yt = [];
    pos_xtc = [];
    pos_ytc = [];
    theta = [];
    thetac = [];
    nose_x = [];
    nose_y = [];
    theta_h = [];
    
    for i = x1_tracker:x2_tracker
        % Generate data for scatterplot
        ntraces = size(Tracker.Parameters{i},1);
        
        if ntraces > 0
            pos_xo(end+1:end+ntraces,1:2) = [ones(ntraces,1)*i, Tracker.Parameters{i}(:,1)];
            pos_yo(end+1:end+ntraces,1:2) = [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,2)];
            pos_xoc(end+1:end+ntraces,1:2) = [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,9)];
            pos_yoc(end+1:end+ntraces,1:2) = [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,10)];
            pos_xt(end+1:end+ntraces,1:2) = [ones(ntraces,1)*i, Tracker.Parameters{i}(:,3)];
            pos_yt(end+1:end+ntraces,1:2) = [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,4)];
            pos_xtc(end+1:end+ntraces,1:2) = [ones(ntraces,1)*i, Tracker.Parameters{i}(:,11)];
            pos_ytc(end+1:end+ntraces,1:2) = [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,12)];
            theta(end+1:end+ntraces,1:2) = [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,5)];
            thetac(end+1:end+ntraces,1:2)= [ones(ntraces, 1)*i, Tracker.Parameters{i}(:,6)];
            nose_x(end+1,1:2) = [i Tracker.Parameters{i}(1,14)];
            nose_y(end+1,1:2) = [i Tracker.Parameters{i}(1,15)];
            theta_h(end+1,1:2) = [i Tracker.Parameters{i}(1,13)];
        end
        
        
        
        
        
    end
    
    scatter(ax{1}, pos_xo(:,1), pos_xo(:,2),markersize, 'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
   scatter(ax{5}, pos_yo(:,1), pos_yo(:,2),markersize, 'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
        scatter(ax{3}, pos_xt(:,1), pos_xt(:,2),markersize, 'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
   scatter(ax{7}, pos_yt(:,1), pos_yt(:,2),markersize, 'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
   scatter(ax{9}, theta(:,1), theta(:,2),markersize, 'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
       scatter(ax{10}, thetac(:,1), thetac(:,2),markersize, 'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
    scatter(ax{2}, pos_xoc(:,1), pos_xoc(:,2),markersize,'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
    scatter(ax{6}, pos_yoc(:,1), pos_yoc(:,2),markersize,'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
     scatter(ax{4}, pos_xtc(:,1), pos_xtc(:,2),markersize,'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
    scatter(ax{8}, pos_ytc(:,1), pos_ytc(:,2),markersize,'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor',rawcolor);
    scatter(ax{11}, nose_x(:,1), nose_x(:,2),markersize,'MarkerFaceColor','r',...
        'MarkerEdgeColor','r');
    scatter(ax{11}, nose_y(:,1), nose_y(:,2),markersize,'MarkerFaceColor','b',...
        'MarkerEdgeColor','b');
    scatter(ax{12}, theta_h(:,1), theta_h(:,2), markersize,'MarkerFaceColor',rawcolor,...
        'MarkerEdgeColor','b');
end



title(ax{1},'Origin-X');                ylim(ax{1},[0,hx]);
title(ax{2},'Origin-X-corrected');      ylim(ax{2},[-200 200])      
title(ax{3},'Tip-X');                   ylim(ax{3},[0,hx]);
title(ax{4},'Tip-X-corrected');         ylim(ax{4},[-200 200])
title(ax{5},'Origin-Y');                ylim(ax{5},[0,hy]);  
title(ax{6},'Origin-Y-corrected');      ylim(ax{6},[-200 200])
title(ax{7},'Tip-Y');                   ylim(ax{7},[0,hy]);
title(ax{8},'Tip-Y-corrected');         ylim(ax{8},[-200 200])
title(ax{9},'theta');                   ylim(ax{9},[-180 180]);
title(ax{10},'theta-corrected');        ylim(ax{10},[-180 180]);
title(ax{11},'Nose position');          ylim(ax{11},[0 hy]);
title(ax{12},'Head angle');             ylim(ax{12},[-180 180]);












