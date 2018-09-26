function Trace = getTrace(OB,Settings)

%% TrackFrame
thr = Settings.Silhouettethreshold*255;

Shape = zeros(Settings.Video_width, Settings.Video_heigth);

KL = ones(9,9)./81;
KS = ones(3,3)./9;

f = fopen(Settings.Video, 'r');

fseek(f, (Settings.Current_frame-1)*Settings.Video_width*Settings.Video_heigth,'bof');


frame = fread(f, [Settings.Video_width, Settings.Video_heigth],'*uint8');

LPl = conv2(frame, KL,'same');
LPs = conv2(frame, KS,'same');


Shape(:) = 0;
Shape(frame <= thr) = 1;
Shape(OB>0) = 0;


ShapeSmall = imerode(Shape, strel('diamond', 2));
ShapeLarge = imdilate(ShapeSmall, strel('disk', Settings.Dilationsize));




SeedEdge = find(edge(ShapeLarge));

[Sx, Sy] = ind2sub(size(OB), SeedEdge);

unmarked = 2:length(SeedEdge);
marked = 1;
temp = zeros(length(SeedEdge), 2);
temp(1,:) = [Sx(1),Sx(2)];
idx = 2;
while ~isempty(unmarked)
    d = sqrt(sum( ([Sx(unmarked),Sy(unmarked)] - [temp(end,1),temp(end,2)]).^2, 2));
    [~, id] = min(d);
    temp(idx,:) = [Sx(unmarked(id)), Sy(unmarked(id))];
    marked = [marked, unmarked(id)];
    unmarked = [unmarked(1:id-1), unmarked(id+1:end)];
    idx = idx+1;
end
SeedEdge = sub2ind(size(frame), temp(:,1), temp(:,2));
SeedProfile = double(frame(SeedEdge));

%LP = medfilt1(SeedProfile, 500);
LP = conv(SeedProfile, ones(1,500)./500, 'same');
SeedProfile = abs(SeedProfile - LP);

[~, SeedIdx] = findpeaks(SeedProfile, 'Threshold', Settings.Origin_threshold);

SeedsX = temp(SeedIdx, 1);
SeedsY = temp(SeedIdx, 2);    


Trace = [];
angles = repmat(1:20:360,[length(SeedsX),1]);
theta = [];
theta(1,:,:) = 3*cosd(angles);
theta(2,:,:) = 3*sind(angles);
Trace(1,:,1) = SeedsX;
Trace(2,:,1) = SeedsY;
T = repmat(Trace(:,:,1), [1, 1, size(angles,2)]);
res = round(theta+T);
res(res<1) = 1;
res(res(1,:,:) > 512, :, :) = 1;
res(res(2,:,:) > 640,:,:) = 1;
res(res>512) = 1;
Cidx = squeeze(sub2ind(size(frame),  res(1,:,:), res(2,:,:)));
Proi = frame(Cidx);
Proi(find(ShapeLarge(Cidx))) = 1;
[~, minidx] = min(Proi, [], 2);
f = sub2ind(size(Cidx),1:length(minidx),minidx');
[Trace(1,:,2), Trace(2,:,2)] = ind2sub(size(frame), Cidx(f));

H = -30:30;

frame(1,1) = 999;

for i = 1:100
    loopidx = size(Trace,3);
    if ~any(any(Trace(:,:,end) ~= 1))
        break
    end
    
    vt = Trace(:,:,loopidx) - Trace(:,:,loopidx-1);
    Angles = atan2d(vt(2,:), vt(1,:));
    angles = Angles'+H;
    theta = [];
    theta(1,:,:) = 5*cosd(angles);
    theta(2,:,:) = 5*sind(angles);
    T = repmat(Trace(:,:,loopidx), [1, 1, size(angles,2)]);
    res = round(theta+T);
    res(res<1) = 1;
    res(res(1,:,:) > 640, :, :) = 1;
    res(res(2,:,:) > 512,:,:) = 1;
    res(res>512) = 1;
    
    Cidx = squeeze(sub2ind(size(frame), res(1,:,:), res(2,:,:)));
    Proi = frame(Cidx);
    [~, minidx] = min(Proi, [], 2);
    f=  sub2ind(size(Cidx), 1:length(minidx), minidx');
    nanidx = [find(OB(Cidx(f))), find(LPs(Cidx(f))./LPl(Cidx(f)) > Settings.trace_threshold)];

    [Trace(1,:,loopidx+1), Trace(2,:,loopidx+1)] = ind2sub(size(frame), Cidx(f));
    Trace(1:2,nanidx, loopidx+1) = 1;
end

Trace(Trace == 1 ) = NaN;

