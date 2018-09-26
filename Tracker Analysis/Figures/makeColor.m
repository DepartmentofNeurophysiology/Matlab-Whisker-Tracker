function colors = makeColor(gapwidth)
%%
warning('off')
% Color pallet for visual output
cc = cbrewer('seq','YlGn',20);
colors.tracker_light = cc(9,:);
colors.tracker_dark = cc(12,:);
colors.tracker_touch = cc(9,:);
colors.tracker_touch_style = 'o';
cc = cbrewer('div','RdYlGn',12);
colors.manual_light = cc(4,:);
colors.manual_dark = cc(2,:);
colors.manual_touch = cc(4,:);
colors.manual_touch_style = '^';

cc = cbrewer('seq','Oranges',20);
colors.janelia_light = cc(9,:);
colors.janelia_dark = cc(12,:);
colors.janelia_touch = cc(9,:);


cc =cbrewer('seq','Oranges',12);
colors.nose = cc(6,:);

cc= cbrewer('seq','OrRd',12);
colors.raw = cc(10,:);


if exist('gapwidth','var')
cc = cbrewer('div','RdYlGn',gapwidth);
colors.dist_nose_target = [cc; ones(200,3).*cc(end,:)];
end

cc = cbrewer('div','RdGy',100);
colors.distmap = cc;



warning('on')




