function[ disttarget, exploring] = getDistTarget(Tracker)
%%

%load('E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance\M47_R04_07_Annotations_Tracker.mat')
%Tracker.gapinfo = detectGap(Tracker.Objects);

%%
switch(Tracker.Direction)
    case 'Up'
        target = Tracker.gapinfo.edge_1;
        sign = 1;
    case 'Down'
        target = Tracker.gapinfo.edge_2;
        sign = -1;

    otherwise
        target = NaN;
end


if  Tracker.gapinfo.main_ax == 2
    idx = 1;
elseif Tracker.gapinfo.main_ax == 1
    idx = 2;
end

Nose = Tracker.Nose(:,idx);
disttarget = sign* (Nose-target);



gapsize = abs(Tracker.gapinfo.edge_1 - Tracker.gapinfo.edge_2);
idx = find(~isnan(disttarget) & disttarget <= gapsize & disttarget >= -40);
exploring = zeros(1, length(disttarget));
exploring(idx) = 1;



