close all
ax_width = 200;
ax_heigth = 200;

marg_x = 50;
x_space = 50;
fig_width = 3*ax_width + 2*marg_x + 2*x_space+50;

marg_y = 40;
fig_heigth = ax_heigth + 2*marg_y;

f = figure('Units','points','Position',[10 100 fig_width fig_heigth]);

x = marg_x+10;
y = fig_heigth - ax_heigth-  marg_y+10;
ax1 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);
x = x+ax_width+x_space;
ax2 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);
x = x+ax_width+x_space;
ax3 = axes(f, 'Units','points','Position',[x y ax_width ax_heigth]);





%%


load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi\DataFigBeh.mat')

c = makeColor;



%% Ax 1


Data.table_sorted_1.GroupCount = Data.table_sorted_1.GroupCount./max(...
    [Data.table_sorted_1.GroupCount(sidx); Data.table_sorted_1.GroupCount(midx)]);

hold(ax1, 'on')
s1 = scatter(ax1, Data.show.single.xax, Data.show.single.data, ...
    'MarkerFaceColor', c.Single, 'MarkerEdgeColor', c.Single);
s2 = scatter(ax1, Data.show.multi.xax , Data.show.multi.data, ...
    'MarkerFaceColor',c.Multi,'MarkerEdgeColor',c.Multi);




plot(ax1, Data.fit.single.showax, Data.fit.single.showfit, 'color', c.Single)
plot(ax1, Data.fit.multi.showax, Data.fit.multi.showfit, 'color', c.Multi)

ax1.XTick = [25 55];
ax1.YTick = [0 0.5 1];
xlabel(ax1, 'Gap width (mm)')
ylabel(ax1, 'Nomalized touch count')
xlim(ax1, [x0 x1])

%% Ax 2

hold(ax2, 'on')
plot(ax2,Data.dist_v_ntouch(:,1), Data.dist_v_ntouch(:,2), 'color', c.Single);
plot(ax2,Data.dist_v_ntouch(:,1), Data.dist_v_ntouch(:,3), 'color',c.Multi);
set(ax2, 'XDir','reverse')
xlim(ax2, [-20 30])
xlabel(ax2, 'Distance from target')
ylabel(ax2, '%Touch')

%% Ax 3

hold(ax3, 'on')
plot(ax3, Data.dist_v_count(:,1), Data.dist_v_count(:,2), 'color', c.Single);
 plot(ax3, Data.dist_v_count(:,1), Data.dist_v_count(:,3), 'color', c.Multi);
set(ax3, 'XDir', 'reverse')
xlabel(ax3, 'Distance from target')
xlim(ax3, [-20 60])
ylabel(ax3, 'p(Time)')


l = legend(ax3, [s1 s2],{'Single','Multi'},'Location','east');
