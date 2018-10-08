warning('off')
cmap = cbrewer('div','BrBG',30);
c.Manual = cmap(5,:);
c.Tracker = cmap(25,:);

cmap = cbrewer('seq','Blues',20);
c.Blue = cmap(19,:);

cmap = cbrewer('div','RdBu',200);
cmap = flip( cmap(floor(size(cmap,1)/2):end,:), 1);
c.BlueMap = cmap;

cmap = cbrewer('seq','YlOrBr',200);
c.OrangeMap = flip(cmap,1);

%data_fig
%% Make figure

close all

figure_width = 750;
figure_heigth = 500;

marg_x = 50;
marg_y = 40;


im_x_range = [120 620];
im_y_range = [220 480];
xwidth = abs(diff(im_x_range));
ywidth = abs(diff(im_y_range));
im_ax_heigth = 110;
im_ax_width = im_ax_heigth*(xwidth/ywidth);

theta_ax_heigth = 130;
theta_ax_width = 300;

corr_ax_heigth = 80;
corr_ax_width = 80;

inset_width = 180;
inset_heigth = 90;

par_heigth = 80;
par_width = 80;

touch_width = 100;
touch_heigth =100;

f = figure('Units','points','Position',[500 50 figure_width figure_heigth]);

x = marg_x;
y = figure_heigth-marg_y-im_ax_heigth;

% Axis for frame
im_x_range = [120 620];
im_y_range = [220 480];
xwidth = abs(diff(im_x_range));
ywidth = abs(diff(im_y_range));
ax1 = axes(f,'Units','points',...
    'Position',[x y+10 im_ax_width im_ax_heigth],...
    'Visible','off');

% Axis for thetas
nleft = 70;
x = marg_x+im_ax_width+nleft;
ax2 = axes(f,'Units','points',...
    'Position',[x y theta_ax_width theta_ax_heigth],...
    'Visible','on');



% Axis for paramter set-in
ntop = 0;
y = figure_heigth-marg_y-im_ax_heigth-ntop-theta_ax_heigth;
nleft = 30;
x = marg_x;
ax4 = axes(f,'Units','points',...
    'Position',[x y inset_width inset_heigth],...
    'Visible','on');


% Axis for amplitude measure
ntop = 0;
y = figure_heigth-marg_y-im_ax_heigth-ntop-theta_ax_heigth-10;
x = marg_x + im_ax_width+80;

ax5 = axes(f,'Units','points',...
    'Position',[x y par_width par_heigth],...
    'Visible','on');


% Axis for duration measure
nleft = 30;
x = x+par_width+nleft;
ax6 = axes(f,'Units','points',...
    'Position',[x y par_width par_heigth],...
    'Visible','on');

% Axis for speed measure
x = x+nleft+par_width;
ax7 = axes(f,'Units','points',...
    'Position',[x y par_width par_heigth],...
    'Visible','on');
cbar7 = colorbar(ax7, 'Location','Manual','Units','points','Position',[0 0 5 10]);


% Axis for amplitude
ntop = 30;
y = y-ntop-par_heigth;
x = marg_x + im_ax_width + 80;
ax8 = axes(f,'Units','points',...
    'Position',[x y par_width par_heigth],...
    'Visible','on');

% Axis for duration
x = x+par_width+nleft;
ax9 = axes(f,'Units','points',...
    'Position',[x y par_width par_heigth],...
    'Visible','on');


% Axis for speed
x = x+par_width+nleft;
ax10 = axes(f,'Units','points',...
    'Position',[x y par_width par_heigth],...
    'Visible','on');

% Axis for touch
x = marg_x;
ax11 = axes(f,'Units','points',...
    'Position',[x,y-20, touch_width, touch_heigth],...
    'Visible','on');

% Axis for correlation

x = marg_x+ 130;
ax3 = axes(f,'Units','points',...
    'Position',[x y corr_ax_width corr_ax_heigth],...
    'Visible','on');

%% AX1 - EXAMPLE FRAME

imagesc(ax1, Data.Frame)
colormap gray
caxis(ax1, [0 0.6])
xlim(ax1,im_x_range)
ylim(ax1,im_y_range)
set(ax1,'Visible','off')

hold(ax1, 'on')
for i = 1:size(Data.TTraces,2)
    if Data.TCRflag(i)
        plot(ax1, Data.TTraces{i}(:,2), Data.TTraces{i}(:,1), 'color',c.Tracker,...
            'LineStyle','--','LineWidth',1)
    else
        plot(ax1, Data.TTraces{i}(:,2), Data.TTraces{i}(:,1), 'color',c.Tracker,...
            'LineStyle','-','LineWidth',1)
    end
end

for i = 1:length(Data.MCRflag)
    if Data.MCRflag(i) == 1 || Data.MCRflag(i) == 0
        plot(ax1,Data.MTraces{i}(:,2), Data.MTraces{i}(:,1), 'color',c.Manual,...
            'LineStyle','-','LineWidth',1)
    elseif Data.MCRflag(i) == 2
        plot(ax1,Data.MTraces{i}(:,2), Data.MTraces{i}(:,1), 'color',c.Manual,...
            'LineStyle','--','LineWidth',1)
    end
end



%% AX2 - EXAMPLE THETAS

hold(ax2,'on')
plot(ax2, Data.xax, Data.Angles.Tracker.r_max_filtered,'color',c.Tracker,'LineStyle','-')
plot(ax2, Data.xax, Data.Angles.Tracker.r_min_filtered,'color',c.Tracker,'LineStyle','--')
plot(ax2, Data.xax, Data.Angles.Tracker.l_min_filtered,'color',c.Tracker,'LineStyle','--')
plot(ax2, Data.xax, Data.Angles.Tracker.l_max_filtered,'color',c.Tracker,'LineStyle','-')

plot(ax2, Data.xax, Data.Angles.Manual.r_max_filtered,'color',c.Manual,'LineStyle','-')
plot(ax2, Data.xax, Data.Angles.Manual.r_min_filtered,'color',c.Manual,'LineStyle','--')
plot(ax2, Data.xax, Data.Angles.Manual.l_min_filtered,'color',c.Manual,'LineStyle','--')
plot(ax2, Data.xax, Data.Angles.Manual.l_max_filtered,'color',c.Manual,'LineStyle','-')

xlim(ax2,[4300 4550])
ylim(ax2,[-120 120])
xlabel(ax2, 'Time [ms]')
ylabel(ax2, '\theta')
ax2.XTick = [4300:100:4500];

apos = ax2.Position;
% annotation(f, 'line',[apos(1)+theta_ax_width, apos(1)+theta_ax_width+7]/figure_width,...
%     [ones(1,2)*(apos(2)+160)]/figure_heigth)
% annotation(f, 'line',[apos(1)+theta_ax_width, apos(1)+theta_ax_width+7]/figure_width,...
%     [ones(1,2)*(apos(2)+135)]/figure_heigth)
% annotation(f, 'line',[apos(1)+theta_ax_width, apos(1)+theta_ax_width+7]/figure_width,...
%     [ones(1,2)*(apos(2)+75)]/figure_heigth)
% annotation(f, 'line',[apos(1)+theta_ax_width, apos(1)+theta_ax_width+7]/figure_width,...
%     [ones(1,2)*(apos(2)+50)]/figure_heigth)
% 
% text(ax2, theta_ax_width+8, 120, sprintf('R: %1.2f',Data.corr.r_max.R),...
%     'Units','points','FontSize',8)
% text(ax2, theta_ax_width+8, 95, sprintf('R: %1.2f',Data.corr.r_min.R),...
%     'Units','points','FontSize',8)
% text(ax2, theta_ax_width+8, 40, sprintf('R: %1.2f',Data.corr.l_min.R),...
%     'Units','points','FontSize',8)
% text(ax2, theta_ax_width+8, 15, sprintf('R: %1.2f',Data.corr.l_max.R),...
%     'Units','points','FontSize',8)
%% AX3 - CORRELATION

edges = [0:0.05:1];
h = hist(abs(Data.Correlation.Total),edges);
b = bar(ax3,edges,h/sum(h), 'EdgeColor',c.Blue,'FaceColor',c.Blue,'BarWidth',1);
xlim(ax3,[-0.05 1])
xlabel(ax3, 'R')


%% AX4 - Theta inset
cla(ax4)
hold(ax4,'on')
plot(ax4, Data.xax, Data.Angles.Tracker.r_max_filtered,'color',c.Tracker,'LineStyle','-')
plot(ax4, Data.xax, Data.Angles.Manual.r_max_filtered,'color',c.Manual,'LineStyle','-')

scatter(ax4, Data.xax(Data.Ttrogh2), Data.Angles.Tracker.r_max_filtered(Data.Ttrogh2), ...
      'MarkerFaceColor',c.Tracker,'MarkerEdgeColor','k')
scatter(ax4, Data.xax(Data.TPeak), Data.Angles.Tracker.r_max_filtered(Data.TPeak),...
    'MarkerFaceColor',c.Tracker,'MarkerEdgeColor','k')
scatter(ax4, Data.xax(Data.Ttrogh), Data.Angles.Tracker.r_max_filtered(Data.Ttrogh),...
    'MarkerFaceColor',c.Tracker,'MarkerEdgeColor','k')

scatter(ax4, Data.xax(Data.MPeak), Data.Angles.Manual.r_max_filtered(Data.MPeak),...
    'MarkerFaceColor',c.Manual,'MarkerEdgeColor','k')
scatter(ax4, Data.xax(Data.Mtrogh),Data.Angles.Manual.r_max_filtered(Data.Mtrogh),...
    'MarkerFaceColor',c.Manual,'MarkerEdgeColor','k')
scatter(ax4, Data.xax(Data.Mtrogh2), Data.Angles.Manual.r_max_filtered(Data.Mtrogh2), ...
     'MarkerFaceColor',c.Manual,'MarkerEdgeColor','k')


line(ax4, [Data.xax(Data.MPeak)-2, Data.xax(Data.Mtrogh)], ones(1,2)*Data.Angles.Manual.r_max_filtered(Data.Mtrogh),...
    'color',c.Manual,'LineStyle','--')
line(ax4, ones(1,2)*Data.xax(Data.MPeak)-1, [Data.Angles.Manual.r_max_filtered(Data.MPeak),...
    Data.Angles.Manual.r_max_filtered(Data.Mtrogh)],  'color',c.Manual,'LineStyle','--')
line(ax4, [Data.xax(Data.MPeak)+2, Data.xax(Data.Mtrogh2)], ones(1,2)*Data.Angles.Manual.r_max_filtered(Data.Mtrogh2),...
    'color',c.Manual,'LineStyle','--')
line(ax4, ones(1,2)*Data.xax(Data.MPeak)+1, [Data.Angles.Manual.r_max_filtered(Data.MPeak),...
    Data.Angles.Manual.r_max_filtered(Data.Mtrogh2)],  'color',c.Manual,'LineStyle','--')


line(ax4, ones(1,2)*Data.xax(Data.Ttrogh), [Data.Angles.Tracker.r_max_filtered(Data.TPeak),...
    Data.Angles.Tracker.r_max_filtered(Data.Ttrogh)],  'color',c.Tracker,'LineStyle','--')
line(ax4, Data.xax([Data.TPeak, Data.Ttrogh]), ones(1,2)*Data.Angles.Tracker.r_max_filtered(Data.MPeak),...
    'color',c.Tracker,'LineStyle','--')
line(ax4, ones(1,2)*Data.xax(Data.Ttrogh2), [Data.Angles.Tracker.r_max_filtered(Data.TPeak),...
    Data.Angles.Tracker.r_max_filtered(Data.Ttrogh2)],  'color',c.Tracker,'LineStyle','--')
line(ax4, Data.xax([Data.TPeak, Data.Ttrogh2]), ones(1,2)*Data.Angles.Tracker.r_max_filtered(Data.MPeak),...
    'color',c.Tracker,'LineStyle','--')


text(ax4, Data.xax(Data.Ttrogh)+3, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 5,...
    '\it Dur.')
text(ax4, Data.xax(Data.Ttrogh)-3, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 15,...
    '\it Amp.','Rotation',90)
text(ax4, Data.xax(Data.Ttrogh)+3, Data.Angles.Manual.r_max_filtered(Data.Mtrogh) + 3,...
    '\it Dur.')
text(ax4, Data.xax(Data.MPeak)-4, Data.Angles.Manual.r_max_filtered(Data.MPeak) + 3,...
    '\it Amp.','Rotation',90)

text(ax4, Data.xax(Data.TPeak)+3, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 5,...
    '\it Dur.')
text(ax4, Data.xax(Data.Ttrogh2)-3, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 7,...
    '\it Amp.','Rotation',90)
text(ax4, Data.xax(Data.TPeak)+3, Data.Angles.Manual.r_max_filtered(Data.Mtrogh2) + 3,...
    '\it Dur.')
text(ax4, Data.xax(Data.MPeak)+3, Data.Angles.Manual.r_max_filtered(Data.MPeak) + 2,...
    '\it Amp.','Rotation',90)


text(ax4, (Data.xax(Data.MPeak)+Data.xax(Data.Mtrogh))/2 - 10 , ...
    Data.Angles.Manual.r_max_filtered(Data.Mtrogh)+10, 'Protraction')
text(ax4, (Data.xax(Data.MPeak)+Data.xax(Data.Mtrogh2))/2 - 8 , ...
    Data.Angles.Manual.r_max_filtered(Data.Mtrogh2)+10, 'Retraction')

xlim(ax4,[4400 4500])
ylim(ax4,[50 120])
xlabel(ax4, 'Time [ms]')
ylabel(ax4, '\theta')
ax4.XTick = [4400 4500];
ax4.YTick = [50 120];




%% AX5 - Measured amplitudes

imagesc(ax5, Data.BIN_AMPLITUDE_PRO)
set(ax5, 'YDir','normal')
colormap(ax5, c.BlueMap)
ax5.XTick = Data.BIN_AMPLITUDExtickax;
ax5.XTickLabel = Data.BIN_AMPLITUDExticklabels;
ax5.YTick = Data.BIN_AMPLITUDExtickax;
ax5.YTickLabel = Data.BIN_AMPLITUDExticklabels;
title(ax5, 'Amplitude [deg]')
ylabel(ax5, 'Tracker')
text(ax5, -0.5, 0.1, 'PROTRACTION','Units','normalized', 'Rotation',90);


%% AX6 - Measured durations
imagesc(ax6, Data.BIN_DURATION_PRO)
set(ax6, 'YDir', 'normal')
colormap(ax6, c.BlueMap)
ax6.XTick = Data.BIN_DURATIONxtickax;
ax6.XTickLabel = Data.BIN_DURATIONxticklabels;
ax6.YTick = Data.BIN_DURATIONxtickax;
ax6.YTickLabel = Data.BIN_DURATIONxticklabels;




title(ax6,'Duration [ms]')


%% AX7 - Measured speeds
imagesc(ax7, Data.BIN_SPEED_PRO)
set(ax7, 'YDir', 'normal')
colormap(ax7, c.BlueMap)
ax7.XTick = Data.BIN_SPEEDxtickax;
ax7.XTickLabel = Data.BIN_SPEEDxticklabels;
ax7.YTick = Data.BIN_SPEEDxtickax;
ax7.YTickLabel = Data.BIN_SPEEDxticklabels;

title(ax7 ,'Speed [deg/s]')



cpos(1:2) = ax7.Position(1:2) + [par_width+10, par_heigth/2-50];
cpos(3:4) = [10 100];
cbar = colorbar(ax7, 'Units','points','Location','Manual',...
    'Position',cpos);
cbar.Label.String = 'normalized count';




%% AX8 - Measured amplitudes

imagesc(ax8, Data.BIN_AMPLITUDE_CON)
set(ax8, 'YDir','normal')
colormap(ax8, c.BlueMap)
ax8.XTick = Data.BIN_AMPLITUDExtickax;
ax8.XTickLabel = Data.BIN_AMPLITUDExticklabels;
ax8.YTick = Data.BIN_AMPLITUDExtickax;
ax8.YTickLabel = Data.BIN_AMPLITUDExticklabels;




xlabel(ax8, 'Manual')
ylabel(ax8, 'Tracker')
text(ax8, -0.5, 0.1, 'RETRACTION','Units','normalized', 'Rotation',90);


%% AX9 - duration contraction

imagesc(ax9, Data.BIN_DURATION_CON)
set(ax9, 'YDir', 'normal')
colormap(ax9, c.BlueMap)
ax9.XTick = Data.BIN_DURATIONxtickax;
ax9.XTickLabel = Data.BIN_DURATIONxticklabels;
ax9.YTick = Data.BIN_DURATIONxtickax;
ax9.YTickLabel = Data.BIN_DURATIONxticklabels;



xlabel(ax9, 'Manual')


%% AX10 - duration contraction

imagesc(ax10, Data.BIN_SPEED_CON)
set(ax10, 'YDir', 'normal')
colormap(ax10, c.BlueMap)
ax10.XTick = Data.BIN_SPEEDxtickax;
ax10.XTickLabel = Data.BIN_SPEEDxticklabels;
ax10.YTick = Data.BIN_SPEEDxtickax;
ax10.YTickLabel = Data.BIN_SPEEDxticklabels;



xlabel(ax10, 'Manual')


%%

imagesc(ax11, Data.TOUCHCOUNT)
set(ax11, 'YDir','normal')
colormap(ax11, c.OrangeMap)
ax11.XTick = Data.TOUCHtickax;
ax11.XTickLabel = Data.TOUCHticklabel;
ax11.YTick = Data.TOUCHtickax;
ax11.YTickLabel = Data.TOUCHticklabel;

title(ax11, 'Touch')
xlabel(ax11, 'Manual')
ylabel(ax11, 'Tracker')









