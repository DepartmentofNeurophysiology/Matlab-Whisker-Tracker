function output = compareThetas(T,Tpeak,M,Mpeak, framerange)
%%

if length(T) > length(M)
    T = T(1:length(M));
end

if length(framerange) > length(T)
    framerange = framerange(1:length(T));
end

IDX = find(~isnan(T) & ~isnan(M) & framerange);

T(isnan(M)) = NaN;
if ~any(find(~isnan(M)))  
    output.data = 'manual_missing';
    output.R = NaN;
    output.P = NaN;
    output.res = NaN;
    output.miss = NaN;
    return
    
elseif ~any(find(~isnan(T)))   
    output.data = 'tracker_missing';
    output.R = NaN;
    output.P = NaN;
    output.res = NaN;
    output.miss = NaN;
    return
    
elseif length(IDX) < 5  
        output.R = NaN;
    output.P = NaN;
    output.res = NaN;
    output.miss = NaN;
    return
end

[R,P] = corrcoef(T(IDX), M(IDX));

Tpeak = Tpeak(ismember(Tpeak, IDX));
Mpeak = Mpeak(ismember(Mpeak, IDX));


dpeaks = abs(Mpeak' - Tpeak);
res = mean(min(dpeaks, [] ,2));
miss = length(Mpeak) - length(Tpeak);

% if miss > 4
%     keyboard
% end
% 
% p1 = plot(ax,  T, 'r');
% p2 = plot(ax,  M, 'b');
% 
% y1 = floor(min([T, M])/10)*10;
% y2 = ceil(max([T, M])/10)*10;
% 
% yrange = 180;
% minY = (y1+y2)/2 - yrange/2;
% maxY = (y1+y2)/2 + yrange/2;

% for i = 1:length(Tpeak)
%     line(ax, [Tpeak(i) Tpeak(i)], [minY, maxY], 'color', 'r','LineStyle','--');
%     scatter(ax, Tpeak(i), T(Tpeak(i)), 'g', 'filled')
% end
% 
% for i = 1:length(Mpeak)
%     line(ax, [Mpeak(i) Mpeak(i)], [minY, maxY], 'color', 'b','LineStyle','--');
%     scatter(ax, Mpeak(i), M(Mpeak(i)), 'g', 'filled')
% end
% 
% title(ax, sprintf('R: %1.3f P: %1.3f', abs(R(1,2)), P(1,2)))
% xlabel(ax, sprintf('Mean distance peaks (MWT-Manual): %2.2f, difference peak count: %d', res, miss))
% 
% ylim(ax, [minY, maxY])
% xlim(ax, [IDX(1), IDX(end)])
% legend(ax, [p1, p2],{'MWT', 'Manual'})


output.R = R(1,2);
output.P = P(1,2);
output.res = res;
output.miss = miss;
output.data = 'available';
