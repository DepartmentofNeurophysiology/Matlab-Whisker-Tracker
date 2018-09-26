clear



load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M46_R01_02_compiled.mat')

Angles = getAngles(Annotations);

i = 1;
figure(1)
clf
ax = subplot(2,2,1);
Out(i).l_min = compareThetas(Angles.Tracker.l_min_filtered, Angles.Tracker.l_min_peaks, ...
    Angles.Manual.r_min_filtered, Angles.Manual.r_min_peaks, ax);

ax = subplot(2,2,2);
Out(i).l_max = compareThetas(Angles.Tracker.l_max_filtered, Angles.Tracker.l_max_peaks, ...
    Angles.Manual.r_max_filtered, Angles.Manual.r_max_peaks, ax);

ax = subplot(2,2,3);
Out(i).r_min = compareThetas(Angles.Tracker.r_min_filtered, Angles.Tracker.r_min_peaks, ...
    Angles.Manual.l_min_filtered, Angles.Manual.l_min_peaks, ax);

ax = subplot(2,2,4);
Out(i).r_max = compareThetas(Angles.Tracker.r_max_filtered, Angles.Tracker.r_max_peaks, ...
    Angles.Manual.l_max_filtered, Angles.Manual.l_max_peaks, ax);













