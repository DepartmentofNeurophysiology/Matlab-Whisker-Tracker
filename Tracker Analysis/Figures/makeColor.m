function colors = makeColor(varargin)
%%
p = inputParser;
p.CaseSensitive = 0;
defaultGap = [];

addParameter(p, 'gapwidth', defaultGap);

parse(p, varargin{:});



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

cmap = cbrewer('seq','Blues',20);
colors.Blue = cmap(19,:);


cmap = cbrewer('div','RdBu',200);
cmap = flip( cmap(floor(size(cmap,1)/2):end,:), 1);
colors.BlueMap = cmap;

cmap = cbrewer('div','BrBG',30);
colors.Manual = cmap(5,:);
colors.Tracker = cmap(25,:);

cmap = cbrewer('seq','YlOrBr',200);
colors.OrangeMap = flip(cmap,1);

if ~isempty(p.Results.gapwidth)
cc = cbrewer('div','RdYlGn',p.Results.gapwidth);
colors.dist_nose_target = [cc; ones(200,3).*cc(end,:)];
end

cc = cbrewer('div','RdGy',100);
colors.distmap = cc;


cmap = cbrewer('seq','Reds',20);
colors.Objects =[0 0 0; cmap(17,:)];

cmap = cbrewer('div','RdBu',20);

colors.Nose = cmap(17,:);
colors.Roi = colors.Nose;
colors.Roi2 = cmap(19,:);
colors.Red = cmap(5,:);
cmap = cbrewer('div','PuOr',20);
colors.Snout = cmap(5,:);
colors.Snout2 = cmap(2,:);
colors.Interest = cmap(17,:);
cmap = cbrewer('seq','Greens',20);
colors.Trace = cmap(10,:);
colors.Trace2 = cmap(15,:);

cmap =cbrewer('seq', 'Greys', 20);
colors.Gray = cmap(10,:);

cmap = cbrewer('div','PRGn',20);
colors.Single = cmap(3,:);
cmap = cbrewer('div','BrBG',20);
colors.Multi = cmap(17,:);

warning('on')




