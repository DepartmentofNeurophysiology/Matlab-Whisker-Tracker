warning('off')
c = makeColor();

VideoDatabase =  'F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\';

if exist(fullfile(VideoDatabase,'Tracker Performance\Data_Figure_Par_Eval.mat'),'file')
    load(fullfile(VideoDatabase,'Tracker Performance\Data_Figure_Par_Eval.mat'))
else    
    Data_Figure_Parameter_Validation;
end


Rs.FontWeight = 'bold';
Rs.Units =  'normalized';
Rs.FontSize  = 9;
Rs.color = 'r';
Fonts.Rs = Rs;

%% Make figure
close all

figure_width = 780;
figure_heigth = 500;

marg_x = 50;
marg_y = 10;


im_x_range = [120 620];
im_y_range = [220 480];
xwidth = abs(diff(im_x_range));
ywidth = abs(diff(im_y_range));
im_ax_heigth = 130;
im_ax_width = im_ax_heigth*(xwidth/ywidth);

theta_ax_heigth = 130;
theta_ax_width = 300;

corr_ax_heigth = 95;
corr_ax_width = 95;

inset_width = 200;
inset_heigth = 90;

par_heigth =100;
par_width = 100;

touch_width = 95;
touch_heigth =95;

f = figure('Units','points','Position',[10 50 figure_width figure_heigth]);

x = marg_x+10;
y = figure_heigth-marg_y-im_ax_heigth;

% Axis for frame
im_x_range = [120 620];
im_y_range = [220 480];
xwidth = abs(diff(im_x_range));
ywidth = abs(diff(im_y_range));
ax1 = axes(f,'Units','points',...
    'Position',[x y im_ax_width im_ax_heigth],...
    'Visible','on');

% Axis for thetas
nleft = 100;
x = marg_x+im_ax_width+nleft;
ax2 = axes(f,'Units','points',...
    'Position',[x y theta_ax_width theta_ax_heigth],...
    'Visible','on');



% Axis for paramter set-in
ntop = 10;
y = figure_heigth-marg_y-im_ax_heigth-ntop-theta_ax_heigth;
nleft = 30;
x = marg_x+40;
ax4 = axes(f,'Units','points',...
    'Position',[x-5 y inset_width inset_heigth],...
    'Visible','on');


% Axis for amplitude measure
ntop = 30;
y = figure_heigth-marg_y-im_ax_heigth-ntop-theta_ax_heigth-20;
x = marg_x + inset_width +120;

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
%cbar7 = colorbar(ax7, 'Location','Manual','Units','points','Position',[0 0 0 10]);


% Axis for amplitude
ntop = 20;
y = y-ntop-par_heigth;
x = marg_x + inset_width +120;
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

ntop = 25;
y = figure_heigth-marg_y-im_ax_heigth-ntop-theta_ax_heigth-10;
y = y-ntop-par_heigth;

% Axis for touch
x = marg_x+150;
ax11 = axes(f,'Units','points',...
    'Position',[x,y, touch_width, touch_heigth],...
    'Visible','on');

% Axis for correlation
x = marg_x;
ax3 = axes(f,'Units','points',...
    'Position',[x+10 y corr_ax_width corr_ax_heigth],...
    'Visible','on');

%% AX1 - EXAMPLE FRAME

imagesc(ax1, Data.Frame)
colormap gray
caxis(ax1, [0 1])
xlim(ax1,im_x_range)
ylim(ax1,im_y_range)
set(ax1,'Visible','off')

hold(ax1, 'on')
for i = 1:size(Data.TTraces,2)
    if Data.TCRflag(i) == 0
        plot(ax1, Data.TTraces{i}(:,2), Data.TTraces{i}(:,1), 'color',c.Tracker,...
            'LineStyle','--','LineWidth',1)
    else
        plot(ax1, Data.TTraces{i}(:,2), Data.TTraces{i}(:,1), 'color',c.Tracker,...
            'LineStyle','-','LineWidth',1)
    end
end

for i = 1:length(Data.MCRflag)
    if Data.MCRflag(i) == 2
        plot(ax1,Data.MTraces{i}(:,2), Data.MTraces{i}(:,1), 'color',c.Manual,...
            'LineStyle','-','LineWidth',1)
    elseif Data.MCRflag(i) == 1 || Data.MCRflag(i) == 0
        plot(ax1,Data.MTraces{i}(:,2), Data.MTraces{i}(:,1), 'color',c.Manual,...
            'LineStyle','--','LineWidth',1)
    end
end

%rectangle(ax1, 'Position',[im_x_range(1) im_y_range(1) 40 40], 'FaceColor',c.Gray,'EdgeColor','none')
text(ax1, im_x_range(1) - 30, im_y_range(1)+20, 'A','FontSize',16)

%% AX2 - EXAMPLE THETAS
cla(ax2)
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

%rectangle(ax2, 'Position',[4300 75 20 60], 'FaceColor',c.Gray,'EdgeColor','none')
text(ax2, 4250, 110, 'B','FontSize',16)


xlabel(ax2, 'Time [ms]','FontWeight','bold')
ax2.XLabel.Position(2) = ax2.XLabel.Position(2)+5;

ylabel(ax2, 'Whisker angle (\theta)','FontWeight','bold')
ax2.YLabel.Position(1) = ax2.YLabel.Position(1)+5;
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
ylabel(ax3, 'count')
title('Correlation')
ax3.YTick = [0 0.15];
ax3.YLabel.Position(1) = ax3.YLabel.Position(1)+0.2;

text(ax3, -0.45, 0.17, 'D','FontSize',16)

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


text(ax4, Data.xax(Data.Ttrogh)+3, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 8,...
    '\it Dur.')
text(ax4, Data.xax(Data.Ttrogh)-5, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 5,...
    '\it Amp.','Rotation',90)
text(ax4, Data.xax(Data.Ttrogh)+3, Data.Angles.Manual.r_max_filtered(Data.Mtrogh) + 8,...
    '\it Dur.')
text(ax4, Data.xax(Data.MPeak)-6, Data.Angles.Manual.r_max_filtered(Data.MPeak) + 5,...
    '\it Amp.','Rotation',90)

text(ax4, Data.xax(Data.TPeak)+5, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 7,...
    '\it Dur.')
text(ax4, Data.xax(Data.Ttrogh2)-5, Data.Angles.Tracker.r_max_filtered(Data.TPeak) + 7,...
    '\it Amp.','Rotation',90)
text(ax4, Data.xax(Data.TPeak)+3, Data.Angles.Manual.r_max_filtered(Data.Mtrogh2) + 7,...
    '\it Dur.')
text(ax4, Data.xax(Data.MPeak)+3, Data.Angles.Manual.r_max_filtered(Data.MPeak) + 5,...
    '\it Amp.','Rotation',90)


text(ax4, (Data.xax(Data.MPeak)+Data.xax(Data.Mtrogh))/2 - 10 , ...
    Data.Angles.Manual.r_max_filtered(Data.Mtrogh)+14, 'Protraction')
text(ax4, (Data.xax(Data.MPeak)+Data.xax(Data.Mtrogh2))/2 - 5 , ...
    Data.Angles.Manual.r_max_filtered(Data.Mtrogh2)+14, 'Retraction')



xlim(ax4,[4400 4510])
ylim(ax4,[40 120])
xlabel(ax4, 'Time [ms]')
ax4.XLabel.Position(2) = ax4.XLabel.Position(2)+10;
ylabel(ax4, '\theta')
ax4.YLabel.Position(1) = ax4.YLabel.Position(1) + 5;
ax4.XTick = [4400 4500];
ax4.YTick = [50 120];


text(ax4, 4375, 115, 'C','FontSize',16)


%% AX5 - Measured amplitudes

imagesc(ax5, Data.histImProA)
set(ax5, 'YDir','normal')
colormap(ax5, c.BlueMap)
bins = Data.Amplitude_bins;
stepsize = mean(diff(bins));
ax5.XTick = [1 length(bins)/2 12];
ax5.XTickLabel = [0 max(bins)/2 max(bins)];
ax5.YTick = [1 length(bins)/2 12];
ax5.YTickLabel = [0 max(bins)/2 max(bins)];
title(ax5, 'Amplitude [deg]')
ylabel(ax5, 'Tracker','Position',[-1.2 6.5 1])
text(ax5, -0.35, 0.14, 'PROTRACTION','Units','normalized', 'Rotation',90,'FontWeight','bold');

text(ax5, 0.65, 0.9, sprintf('r^2: %1.2f', Data.histImProA_Rs), Fonts.Rs)

text(ax5,0-5, 16, 'F','FontSize',16)




%% AX6 - Measured durations
imagesc(ax6, Data.histImProD)
set(ax6, 'YDir', 'normal')
colormap(ax6, c.BlueMap)

bins = Data.Duration_bins;

ax6.XTick = [1 length(bins)/2 length(bins)-1];
ax6.XTickLabel = [0 max(bins)/2 max(bins)];
ax6.YTick = [1 length(bins)/2 length(bins)-1];
ax6.YTickLabel = [0 max(bins)/2 max(bins)];

title(ax6,'Duration [ms]')
text(ax6, 0.65, 0.9, sprintf('r^2: %1.2f', Data.histImProD_Rs), Fonts.Rs)

%% AX7 - Measured speeds
imagesc(ax7, Data.histImProS)
set(ax7, 'YDir', 'normal')
colormap(ax7, c.BlueMap)

bins = Data.Speed_bins;

ax7.XTick = [1 length(bins)/2 length(bins)-1];
ax7.XTickLabel = [0 max(bins)/2 max(bins)];
ax7.YTick = [1 length(bins)/2 length(bins)-1];
ax7.YTickLabel = [0 max(bins)/2 max(bins)];

title(ax7 ,'Speed [deg/s]')

text(ax7, 0.65, 0.9, sprintf('r^2: %1.2f', Data.histImProS_Rs), Fonts.Rs)

cpos(1:2) = ax7.Position(1:2) + [par_width+10, par_heigth/2-110];
cpos(3:4) = [10 100];
cbar = colorbar(ax7, 'Units','points','Location','Manual',...
    'Position',cpos,'Ticks',[0 1]);
cbar.Label.String = 'normalized count';




%% AX8 - Measured amplitudes

imagesc(ax8, Data.histImRetA)
set(ax8, 'YDir','normal')
colormap(ax8, c.BlueMap)
bins = Data.Amplitude_bins;
stepsize = mean(diff(bins));
ax8.XTick = [1 length(bins)/2 12];
ax8.XTickLabel = [0 max(bins)/2 max(bins)];
ax8.YTick = [1 length(bins)/2 12];
ax8.YTickLabel = [0 max(bins)/2 max(bins)];

text(ax8, 0.65, 0.9, sprintf('r^2: %1.2f', Data.histImRetA_Rs), Fonts.Rs)

xlabel(ax8, 'Manual')
ylabel(ax8, 'Tracker','Position',[-1.2 6.5 1])
text(ax8, -0.35, 0.12, 'RETRACTION','Units','normalized', 'Rotation',90,'FontWeight','bold');


%% AX9 - duration contraction

imagesc(ax9, Data.histImRetD)
set(ax9, 'YDir', 'normal')
colormap(ax9, c.BlueMap)

bins = Data.Duration_bins;

ax9.XTick = [1 length(bins)/2 length(bins)-1];
ax9.XTickLabel = [0 max(bins)/2 max(bins)];
ax9.YTick = [1 length(bins)/2 length(bins)-1];
ax9.YTickLabel = [0 max(bins)/2 max(bins)];

xlabel(ax9, 'Manual')
text(ax9, 0.65, 0.9, sprintf('r^2: %1.2f', Data.histImRetD_Rs), Fonts.Rs)


%% AX10 - speed contraction

imagesc(ax10, Data.histImRetS)
set(ax10, 'YDir', 'normal')
colormap(ax10, c.BlueMap)

bins = Data.Speed_bins;

ax10.XTick = [1 length(bins)/2 length(bins)-1];
ax10.XTickLabel = [0 max(bins)/2 max(bins)];
ax10.YTick = [1 length(bins)/2 length(bins)-1];
ax10.YTickLabel =  [0 max(bins)/2 max(bins)];
text(ax10, 0.65, 0.9, sprintf('r^2: %1.2f', Data.histImRetS_Rs), Fonts.Rs)



xlabel(ax10, 'Manual')


%% AZ 11 - touch
cla(ax11)
% imagesc(ax11, Data.TOUCHCOUNT)
% set(ax11, 'YDir','normal')
% colormap(ax11, c.OrangeMap)
% ax11.XTick = Data.TOUCHtickax;
% ax11.XTickLabel = Data.TOUCHticklabel;
% ax11.YTick = Data.TOUCHtickax;
% ax11.YTickLabel = Data.TOUCHticklabel;
% 

% 


edges = [-5:1:5];
h = hist(Data.touch_diff,edges);
b = bar(ax11,edges,h/sum(h), 'EdgeColor',c.Blue,'FaceColor',c.Blue,'BarWidth',1);
xlim(ax11,[-5 5])
ylim(ax11,[0 0.3])
xlabel(ax11, 'R')
ylabel(ax11, 'count')
title('Correlation')
ax11.YTick = [0 0.3];
ax11.YLabel.Position(1) = ax11.YLabel.Position(1)+1;

title(ax11, 'Touchcount')
xlabel(ax11, 'Tracker-Manual')
%ylabel(ax11, 'Tracker')
text(ax11, -2, 11, 'E','FontSize',16)

