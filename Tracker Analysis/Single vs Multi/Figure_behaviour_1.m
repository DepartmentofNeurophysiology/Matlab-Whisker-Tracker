close all
ax_width = 200;
ax_heigth = 200;

marg_x = 50;
x_space = 70;
% fig_width = 3*ax_width + 2*marg_x + 2*x_space+50;
fig_width = 2*ax_width + 2*marg_x + 1*x_space+50;

marg_y = 40;
y_space = 60;
fig_heigth = 2*ax_heigth+y_space + 2*marg_y;

f = figure('Units','points','Position',[10 100 fig_width fig_heigth]);

x = marg_x+10;
y = fig_heigth - ax_heigth-  marg_y+10;
ax1 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);
x = x+ax_width+x_space;
ax2 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);
% x = x+ax_width+x_space;
% ax3 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);

x = marg_x + 10;
y = y - y_space - ax_heigth;
ax4 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);

x = x+ax_width+x_space;
ax5 = axes(f,'Units', 'points', 'Position', [x y ax_width ax_heigth]);





%

load('F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi\DataFigBeh.mat')

c = makeColor;



%% Ax 1


%Data.table_sorted_1.GroupCount = Data.table_sorted_1.GroupCount./max(...
%   Data.table_sorted_1.GroupCount);

hold(ax1, 'on')
s1 = scatter(ax1, Data.ax1.single_width, Data.ax1.single_touch, ...
    'MarkerFaceColor', c.Single, 'MarkerEdgeColor', c.Single);
s2 = scatter(ax1, Data.ax1.multi_width , Data.ax1.multi_touch, ...
    'MarkerFaceColor',c.Multi,'MarkerEdgeColor',c.Multi);

plot(ax1, Data.ax1.single_fit_ax, Data.ax1.single_fit_plot, 'color', c.Single)
plot(ax1, Data.ax1.multi_fit_ax, Data.ax1.multi_fit_plot, 'color', c.Multi)

text(ax1, Data.ax1.single_fit_ax(end-7), Data.ax1.single_fit_plot(end-7)+45, sprintf('R: %1.2f',Data.ax1.single_RS),...
    'color',c.Single)
text(ax1, Data.ax1.multi_fit_ax(end-7), Data.ax1.multi_fit_plot(end-7) + 45, sprintf('R: %1.2f', Data.ax1.multi_RS),...
    'color', c.Multi)

ax1.XTick = [20 40 60];
%ax1.YTick = [0 0.5 1];
xlabel(ax1, 'Gap width (mm)')
ylabel(ax1, '#Touch')
xlim(ax1, [20 60])


text(ax1, 21, 950, 'A', 'FontSize',16)
%% Ax 2

cla(ax2)
hold(ax2, 'on')
scatter(ax2, Data.ax2.single_width, Data.ax2.single_duration*1000, ...
    'MarkerFaceColor', c.Single, 'MarkerEdgeColor', c.Single);
scatter(ax2, Data.ax2.multi_width, Data.ax2.multi_duration*1000, ...
    'MarkerFaceColor', c.Multi, 'MarkerEdgeColor', c.Multi);

plot(ax2, Data.ax2.single_fit_ax, Data.ax2.single_fit_plot*1000, 'color', c.Single)
plot(ax2, Data.ax2.multi_fit_ax, Data.ax2.multi_fit_plot*1000, 'color', c.Multi)

text(ax2, Data.ax2.single_fit_ax(end-7), Data.ax2.single_fit_plot(end-7)*1000+200, sprintf('R: %1.2f', Data.ax2.single_RS),...
    'color', c.Single)
text(ax2, Data.ax2.multi_fit_ax(end-7), Data.ax2.multi_fit_plot(end-7)*1000-50, sprintf('R: %1.2f', Data.ax2.multi_RS), ...
    'color', c.Multi)


xlabel(ax2, 'Gap widt (mm)')
ylabel(ax2, 'Duration (ms)')

xlim(ax2, [20 60])
ax2.XTick = [20 40 60];

text(ax2, 21, 1450, 'B', 'FontSize',16)
%% Ax 3

cla(ax4)
hold(ax4, 'on')
% plot(ax4, Data.ax3.single_dist, Data.ax3.single_touch, 'color', c.Single)
% plot(ax4, Data.ax3.multi_dist,  Data.ax3.multi_touch,  'color', c.Multi)

area(ax4, Data.ax3.multi_dist, Data.ax3.multi_touch, ...
    'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'FaceAlpha', 0.5, ...
    'LineWidth',2);
area(ax4, Data.ax3.single_dist, Data.ax3.single_touch, ...
    'FaceColor', [c.Single], 'EdgeColor', c.Single,'FaceAlpha',0.5,...
    'LineWidth',2);

set(ax4, 'XDir', 'reverse')

xlabel(ax4, 'Distance from target (mm)')
ylabel(ax4, '#Touch')


xlim(ax4, [-20 40])
ax4.YTick = [0:50:100];
ax4.XTick = [-20 0 40];


text(ax4, 38, 95, 'C', 'FontSize', 16)
%% Ax4 


cla(ax5)
hold(ax5, 'on')
% plot(ax5, Data.ax4.single_dist, Data.ax4.single_duration*1000, 'color', c.Single)
% plot(ax5, Data.ax4.multi_dist,  Data.ax4.multi_duration*1000,  'color', c.Multi)

area(ax5, Data.ax4.multi_dist, Data.ax4.multi_duration*1000, ...
    'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'FaceAlpha', 0.5, ...
    'LineWidth',2);
area(ax5, Data.ax4.single_dist, Data.ax4.single_duration*1000, ...
    'FaceColor', c.Single, 'EdgeColor', c.Single, 'FaceAlpha', 0.5, ...
    'LineWidth',2);


set(ax5, 'XDir', 'reverse')
xlabel(ax5, 'Distance from target (mm)')
ylabel(ax5, 'Duration (ms)')

xlim(ax5, [-20 40])



ax5.YTick = [0:50:100];
ax5.XTick = [-20 0 40];





text(ax5, 38, 115, 'D', 'FontSize', 16)



%%


l = legend(ax2, [s2 s1], {'Full','Deprived'},'Position',[0.88 0.5 0.1 0.1]);