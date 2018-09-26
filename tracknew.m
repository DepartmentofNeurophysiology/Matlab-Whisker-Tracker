%%
clc

tic
SL = Shapes.large;
OB = Objects;
OB(OB > 0) = 1;

KL = ones(9,9)./81;
KS = ones(3,3)./9;
LPl = conv2(Shapes.Frame, KL,'same');
LPs = conv2(Shapes.Frame, KS,'same');

Trace = [];
angles = repmat(1:20:360,[33,1]);
theta = [];
theta(1,:,:) = 3*cosd(angles);
theta(2,:,:) = 3*sind(angles);
Trace(:,:,1) = Origins(:,1:2)';
T = repmat(Trace(:,:,1), [1, 1, size(angles,2)]);
res = round(theta+T);
Cidx = squeeze(sub2ind(size(Shapes.Frame),  res(1,:,:), res(2,:,:)));
Proi = Shapes.Frame(Cidx);
Proi(find(SL(Cidx))) = 1;
[~, minidx] = min(Proi, [], 2);
f = sub2ind(size(Cidx),1:length(minidx),minidx');
[Trace(1,:,2), Trace(2,:,2)] = ind2sub(size(Shapes.Frame), Cidx(f));

H = -30:30;

Shapes.Frame(1,1) = 999;

for i = 1:100
    loopidx = size(Trace,3);
    if ~any(any(Trace(:,:,end) ~= 1))
        break
    end
    
    vt = Trace(:,:,loopidx) - Trace(:,:,loopidx-1);
    Angles = atan2d(vt(2,:), vt(1,:));
    angles = Angles'+H;
    theta = [];
    theta(1,:,:) = Settings.stepsize*cosd(angles);
    theta(2,:,:) = Settings.stepsize*sind(angles);
    T = repmat(Trace(:,:,loopidx), [1, 1, size(angles,2)]);
    res = round(theta+T);
    res(res<1) = 1;
    res(res(1,:,:) > f_heigth, :, :) = 1;
    res(res(2,:,:) > f_width,:,:) = 1;
   
    
    Cidx = squeeze(sub2ind(size(Shapes.Frame), res(1,:,:), res(2,:,:)));
    Proi = Shapes.Frame(Cidx);
    [~, minidx] = min(Proi, [], 2);
    f=  sub2ind(size(Cidx), 1:length(minidx), minidx');
    nanidx = [find(OB(Cidx(f))), find(LPs(Cidx(f))./LPl(Cidx(f)) > Settings.trace_threshold)];
    
    [Trace(1,:,loopidx+1), Trace(2,:,loopidx+1)] = ind2sub(size(Shapes.Frame), Cidx(f));
    Trace(1:2,nanidx, loopidx+1) = 1;
end

Trace(Trace == 1 ) = NaN;




