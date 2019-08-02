close all
ax_width = 200;
ax_heigth = 200;

marg_x = 50;
x_space = 50;
fig_width = 3*ax_width + 2*marg_x + 2*x_space+50;

marg_y = 40;
y_space = 10;
fig_heigth = 2*ax_heigth+y_space + 2*marg_y;

f = figure('Units','points','Position',[10 100 fig_width fig_heigth]);

x = marg_x+40;
y = fig_heigth - ax_heigth-  marg_y+10;
ax1 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);
x = x+ax_width+x_space;
ax2 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);
x = x+ax_width+x_space;
ax3 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);

x = marg_x + 40;
y = y - y_space - ax_heigth;
ax4 = axes(f,'Units','points','Position',[x y ax_width ax_heigth]);

x = x+ax_width+x_space;
ax5 = axes(f,'Units', 'points', 'Position', [x y ax_width ax_heigth]);
x = x+ax_width+x_space;
ax6 = axes(f,'Units', 'points', 'Position', [x y ax_width ax_heigth]);




%

load('F:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi\DataFigPar.mat')

c = makeColor;



%% Axes 1

hold(ax1, 'on')
% plot(ax1, Data.Abin_ax, Data.Abin.deprived_pro, 'color', c.Single)
% plot(ax1, Data.Abin_ax, Data.Abin.full_pro, 'color', c.Multi)
area(ax1, Data.Abin_ax, Data.Abin.full_pro, 'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'LineWidth', 2,'FaceAlpha',0.5)

area(ax1, Data.Abin_ax, Data.Abin.deprived_pro, 'FaceColor', c.Single, 'EdgeColor', c.Single, 'LineWidth', 2,'FaceAlpha',0.5)

ax1.XTick = [];
ax1.YTick = [0 0.3];

ylabel(ax1, 'p')


text(ax1, -0.25, 0.3, 'PROTRACTION','Units','normalized', 'Rotation',90,'FontWeight','bold');

%%

hold(ax2, 'on')
% plot(ax2, Data.Tbin_ax, Data.Tbin.deprived_pro, 'color', c.Single)
% plot(ax2, Data.Tbin_ax, Data.Tbin.full_pro, 'color', c.Multi)
area(ax2, Data.Tbin_ax, Data.Tbin.full_pro, 'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'LineWidth', 2,'FaceAlpha',0.5)

area(ax2, Data.Tbin_ax, Data.Tbin.deprived_pro, 'FaceColor', c.Single, 'EdgeColor', c.Single, 'LineWidth', 2,'FaceAlpha',0.5)

ax2.XTick = [];
ax2.YTick = [0 0.2];

%%

hold(ax3, 'on')
% plot(ax3, Data.Sbin_ax, Data.Sbin.deprived_pro, 'color', c.Single)
% plot(ax3, Data.Sbin_ax, Data.Sbin.full_pro, 'color', c.Multi)
area(ax3, Data.Sbin_ax, Data.Sbin.full_pro, 'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'LineWidth', 2,'FaceAlpha',0.5)

area(ax3, Data.Sbin_ax, Data.Sbin.deprived_pro, 'FaceColor', c.Single, 'EdgeColor', c.Single, 'LineWidth', 2,'FaceAlpha',0.5)

ax3.XTick = [];
ax3.YTick = [0 0.1];
ylim(ax3, [0 0.1])


%% Axes 4

hold(ax4, 'on')
% plot(ax4, Data.Abin_ax, Data.Abin.deprived_ret, 'color', c.Single)
% plot(ax4, Data.Abin_ax, Data.Abin.full_ret, 'color', c.Multi)
area(ax4, Data.Abin_ax, Data.Abin.full_ret, 'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'LineWidth', 2,'FaceAlpha',0.5)

area(ax4, Data.Abin_ax, Data.Abin.deprived_ret, 'FaceColor', c.Single, 'EdgeColor', c.Single, 'LineWidth', 2,'FaceAlpha',0.5)

xlabel(ax4, 'Amplitude (deg)')
ylabel(ax4, 'p')

ax4.YTick = [0 0.3];
text(ax4, -0.25, 0.3, 'RETRACTION','Units','normalized', 'Rotation',90,'FontWeight','bold');

%%

hold(ax5, 'on')
% plot(ax5, Data.Tbin_ax, Data.Tbin.deprived_ret, 'color', c.Single)
% plot(ax5, Data.Tbin_ax, Data.Tbin.full_ret, 'color', c.Multi)
area(ax5, Data.Tbin_ax, Data.Tbin.full_ret, 'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'LineWidth', 2,'FaceAlpha',0.5)
area(ax5, Data.Tbin_ax, Data.Tbin.deprived_ret, 'FaceColor', c.Single, 'EdgeColor', c.Single, 'LineWidth', 2,'FaceAlpha',0.5)

ax5.YTick = [0 0.15];
ylim(ax5, [0 0.15])
xlabel(ax5, 'Duration (ms)')

%%

hold(ax6, 'on')
% plot(ax6, Data.Sbin_ax, Data.Sbin.deprived_ret, 'color', c.Single)
% plot(ax6, Data.Sbin_ax, Data.Sbin.full_ret, 'color', c.Multi)
area(ax6, Data.Sbin_ax, Data.Sbin.full_ret, 'FaceColor', c.Multi, 'EdgeColor', c.Multi, 'LineWidth', 2,'FaceAlpha',0.5)
area(ax6, Data.Sbin_ax, Data.Sbin.deprived_ret, 'FaceColor', c.Single, 'EdgeColor', c.Single, 'LineWidth', 2,'FaceAlpha',0.5)

ax6.YTick = [0 0.1];
xlabel(ax6, 'Speed (deg/s)')