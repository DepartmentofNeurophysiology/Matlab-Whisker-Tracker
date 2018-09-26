clear
close all

datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
Files = dir(fullfile(datapath,'*_Annotations_Tracker.mat'));

if ~exist(fullfile(datapath,'Touch'), 'file')
    mkdir(fullfile(datapath, 'Touch'))
end


for fileidx = 1:size(Files,1)

close all
file = fullfile(Files(fileidx).folder, Files(fileidx).name);
load(file)
width = size(Tracker.Objects,1);
heigth = size(Tracker.Objects,2);
r = width/heigth;
aw = 300;
ah = 150;
axh = 100;

f3 = figure('Units', 'points', 'Position',[100 axh 700 150+aw]);


ax1 = axes(f3,'Units','points','Position',[50 axh r*aw aw]);
ax2 = axes(f3, 'Units', 'points', 'Position',[50+r*aw+50 axh+aw-r*ah 300 r*ah]) ;  
ax3 = axes(f3, 'Units', 'points', 'Position',[50+r*aw+50 axh 300 ah]);
                                                                                                                                         

hold(ax1, 'on')
hold(ax2,'on')
hold(ax3,'on');

xlim(ax1,[0 size(Tracker.Objects,1)]);
ylim(ax1,[0 size(Tracker.Objects,2)]);

Settings.Current_frame = round(Settings.Nframes/2);
f = LoadFrame(Settings);
imagesc(ax1, Tracker.Objects')
colormap gray
hold(ax1, 'on')

idx = zeros(1, size(Manual.Traces,1));
for i = 1:length(idx)
    if ~isempty(Manual.Traces{i})
        idx(i) = 1;
    end
end
x1 = find(idx,1,'first');
x2 = find(idx,1,'last');
xlim(ax2,[x1 x2])
xlim(ax3,[x1 x2])
ylim(ax2,[1 size(Output.Objects,1)])
ylim(ax3,[1 size(Output.Objects,2)])


xlabel(ax1,'Position X (px)')
xlabel(ax3,'Frame')

ylabel(ax1,'Position Y (px)')
ylabel(ax2,'Position X (px)')
ylabel(ax3,'Position Y (px)')


for i = x1:x2
   
    
    
    
 if  i <= size(Manual.Touch.pt,2) && ~isempty(Manual.Touch.pt{i})
        for j = 1:size(Manual.Touch.pt{i},1)
            pt = Manual.Touch.pt{i}(j,:);
            scatter(ax2,i,pt(1),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha',0.5)
            scatter(ax3,i,pt(2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha',0.5)
            scatter(ax1,pt(1),pt(2),'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0], 'MarkerFaceAlpha',0.5)
        end
    end
end
    
    saveas(gcf, fullfile(datapath,'Touch',[Files(fileidx).name '.png']))

end