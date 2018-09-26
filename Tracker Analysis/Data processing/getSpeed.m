function [speed_pro, speed_retr] = getSpeed(Theta, Times)
dT = diff(Theta);
negidx = find(dT < 0);

if ~isempty(negidx)
    negT = dT(negidx);
    t = Times([negidx negidx(end)+1]);
    dt = diff(t)';
    speed_pro = 1000*mean(negT./dt);
else
    speed_pro = NaN;
end

posidx = find(dT > 0);

if ~isempty(posidx)
posT = dT(posidx);
t = Times([posidx posidx(end)+1]);
dt = diff(t)';
speed_retr = 1000*mean(posT./dt);
else
    speed_retr = NaN;
end

