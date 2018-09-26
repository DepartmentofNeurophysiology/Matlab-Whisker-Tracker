clear


load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_compiled.mat')

%%
fnr = 1325;

Envelope_man = getEnvelope(Annotations.Manual.Parameters);
Envelope_mwt = getEnvelope(Annotations.Tracker.Parameters_clean);


S = Annotations.Settings;
N = round(Annotations.Tracker.Nose(fnr,:));
S.Current_frame = fnr;
S.Video(1) = 'E';
f = LoadFrame(S);

fb = ~imbinarize(f, 0.3);
msk = zeros(size(fb));
msk(N(1)-80:N(1)+80,N(2)-80:N(2)+80) = 1;
fb =fb.*msk;
fb = bwareafilt(logical(fb), 1);

s_x = sum(fb,1);
s_x = s_x./sum(s_x);
ix = find(s_x);
C(2) = sum(s_x(ix).*ix);

sy = sum(fb,2);
sy = sy./sum(sy);
ix = find(sy);
C(1) = sum(sy(ix).*ix);


%%

cc = cbrewer('div','RdBu',20);
c_man = cc(4,:);
c_mwt = cc(18,:);
cc =cbrewer('div','PRGn',20);
c_jan = cc(18,:);

xrange = [80 600];
yrange = [230 460];

f1 = figure(1);
clf(f1);

f1.Position = [100 50 diff(xrange)+10 2*(10+diff(yrange))];
f1.Units = 'points';
ax1 = axes(f1);
ax1.Position = [0 0*(1/2) 1 (1/2)];
imagesc(ax1, f)
colormap(ax1, 'gray')
hold(ax1, 'on')
ax2 = axes(f1);
ax2.Position = [0 1*(1/2) 1 (1/2)];
imagesc(ax2, f)
colormap(ax2, 'gray')
hold(ax2, 'on')

scatter(ax1, C(2), C(1), 'g','filled')
M = Annotations.Manual.Traces;
for i = 1:size(M{fnr}, 2)
    t = M{fnr}{i};
    plot(ax1, t(:,2), t(:,1), 'r')
end
%drawPie(ax1, C, Envelope_man.left_min(fnr), Envelope_man.left_max(fnr), 200, c_man)
%drawPie(ax1, C, Envelope_man.right_max(fnr), Envelope_man.right_min(fnr), 200, c_man)



T = Annotations.Tracker.Traces_clean;
for i = 1:size(T{fnr}, 2)
    t = T{fnr}{i};
    plot(ax2, t(:,2), t(:,1), 'b')
    text(ax2, t(end,2), t(end,1), num2str(Annotations.Tracker.Parameters_clean{fnr}(i,5)),'color','k')
end
%drawPie(ax2, C, Envelope_mwt.left_min(fnr), Envelope_mwt.left_max(fnr), 200, c_mwt)
%drawPie(ax2, C, Envelope_mwt.right_max(fnr), Envelope_mwt.right_min(fnr), 200, c_mwt)


xlim(ax1, [79 603])
ylim(ax1, [233 457])
xlim(ax2, [79 603])
ylim(ax2, [233 457])





 




