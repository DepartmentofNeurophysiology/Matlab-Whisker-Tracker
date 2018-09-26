%%
clc
SL = Shapes.large;


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



%%
H = -30:30;
vt = Trace(:,:,2) - Trace(:,:,1);
Angles = atan2d(vt(2,:), vt(1,:));
angles = Angles'+H;
theta = [];
theta(1,:,:) = 3*cosd(angles);
theta(2,:,:) = 3*sind(angles);
T = repmat(Trace(:,:,2), [1, 1, size(angles,2)]);
res = round(theta+T);
Cidx = squeeze(sub2ind(size(Shapes.Frame), res(1,:,:), res(2,:,:)));
Proi = Shapes.Frame(Cidx);
[~, minidx] = min(Proi, [], 2);
f=  sub2ind(size(Cidx), 1:length(minidx), minidx');
[Trace(1,:,3), Trace(2,:,3)] = ind2sub(size(Shapes.Frame), Cidx(f));




%%

figure(1)
clf
imshow(Shapes.Frame)
hold on

for i = 1:size(Trace, 2)
    plot(squeeze(Trace(2,i,:)), squeeze(Trace(1,i,:)), 'b')
end


