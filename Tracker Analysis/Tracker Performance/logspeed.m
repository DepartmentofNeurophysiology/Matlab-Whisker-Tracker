function timer = logspeed(timer, buffersize)
if isempty(timer)
    timer.buffer(1:buffersize) = NaN;
    return
end


time = clock;
timer.buffer = circshift(timer.buffer, -1);
timer.buffer(end) = time(4)*3600 + time(5)*60 + time(6);

timer.total = timer.buffer(end) - timer.buffer(1);

if any(isnan(timer.buffer))
    timer.speed = NaN;
else
    timer.speed =length(timer.buffer)/timer.total;
end
