close all
ax_width = 200;
ax_heigth = 200;

marg_x = 50;
x_space = 50;
% fig_width = 3*ax_width + 2*marg_x + 2*x_space+50;
fig_width = 2*ax_width + 2*marg_x + 1*x_space+50;

marg_y = 40;
y_space = 40;
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

load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Singe vs Multi\DataFigBeh.mat')

c = makeColor;



%% Ax 1


%Data.table_sorted_1.GroupCount = Data.table_sorted_1.GroupCount./max(...
%   Data.table_sorted_1.GroupCount);

hold(ax1, 'on')
s1 = scatter(ax1, Data.ax1.single_width, Data.ax1.single_touch, ...
    'MarkerFaceColor', c.Single, 'MarkerEdgeColor', c.Single);
s2 = scatter(ax1, Data.ax1.multi_width , Data.ax1.multi_touch, ...
    'MarkerFaceColor',c.Multi,'MarkerEdgeColor',c.Multi);




%plot(ax1, Data.fit.single.showax, Data.fit.single.showfit, 'color', c.Single)
%plot(ax1, Data.fit.multi.showax, Data.fit.multi.showfit, 'color', c.Multi)

%ax1.XTick = [25 55];
%ax1.YTick = [0 0.5 1];
xlabel(ax1, 'Gap width (mm)')
ylabel(ax1, 'touch count')
xlim(ax1, [24 56])

%% Ax 2

% hold(ax2, 'on')
% plot(ax2,Data.dist_v_ntouch(:,1), Data.dist_v_ntouch(:,2), 'color', c.Single);
% plot(ax2,Data.dist_v_ntouch(:,1), Data.dist_v_ntouch(:,3), 'color',c.Multi);
% set(ax2, 'XDir','reverse')
% xlim(ax2, [-20 30])
% xlabel(ax2, 'Distance from target')
% ylabel(ax2, '%Touch')

hold(ax2, 'on')
scatter(ax2, Data.ax2.single_width, Data.ax2.single_duration, ...
    'MarkerFaceColor', c.Single, 'MarkerEdgeColor', c.Single);
scatter(ax2, Data.ax2.multi_width, Data.ax2.multi_duration, ...
    'MarkerFaceColor', c.Multi, 'MarkerEdgeColor', c.Multi);
xlabel(ax2, 'Gap widt (mm)')
ylabel(ax2, 'Duration (s)')
%% Ax 3

% hold(ax3, 'on')
% plot(ax3, Data.dist_v_count(:,1), Data.dist_v_count(:,2), 'color', c.Single);
%  plot(ax3, Data.dist_v_count(:,1), Data.dist_v_count(:,3), 'color', c.Multi);
% 
% 
% xlim(ax3, [-20 60])
% ylabel(ax3, 'p(Time)')
% 
% 
% l = legend(ax3, [s1 s2],{'Single','Multi'},'Location','east');
cla(ax4)
hold(ax4, 'on')
plot(ax4, Data.ax3.single_dist, Data.ax3.single_touch, 'color', c.Single)
plot(ax4, Data.ax3.multi_dist,  Data.ax3.multi_touch,  'color', c.Multi)
set(ax4, 'XDir', 'reverse')
xlabel(ax4, 'Distance from target (mm)')
ylabel(ax4, 'average Touchcount per trial')
%% Ax4 
% 
% hold(ax4, 'on')
% d = varfun(@sum, Data.table, 'InputVariables', {'nTouch'}, 'GroupingVariables',{'Video','Type'});
% 
% idx = find(strcmp(d.Type, 'Single') & d.GroupCount <= 600);
% scatter(ax4, d.GroupCount(idx), d.sum_nTouch(idx), 'MarkerFaceColor', c.Single, ...
%     'MarkerEdgeColor', c.Single);
% 
% idx = find(strcmp(d.Type, 'Multi') & d.GroupCount <= 600);
% scatter(ax4, d.GroupCount(idx), d.sum_nTouch(idx), 'MarkerFaceColor', c.Multi, ...
%     'MarkerEdgeColor', c.Multi)
% 
% xlabel(ax4, 'Trial duration')
% ylabel(ax4, 'Touch count')

cla(ax5)
hold(ax5, 'on')
plot(ax5, Data.ax4.single_dist, Data.ax4.single_duration*1000, 'color', c.Single)
plot(ax5, Data.ax4.multi_dist,  Data.ax4.multi_duration*1000,  'color', c.Multi)
set(ax5, 'XDir', 'reverse')
xlabel(ax5, 'Distance from target (mm)')
ylabel(ax5, 'average Time per trial (ms)')