function FIG_ERROR_DIST(varargin)
% Number of artifacts as function of distance from target platform
p = inputParser;
p.CaseSensitive = 0;

defaultPath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
addParameter(p,'Path', defaultPath);

parse(p, varargin{:});

datapath = p.Results.Path;

%%

eval_files = dir(fullfile(datapath, '*_Evaluation.mat'));
if exist(fullfile(datapath,'gapwidths.mat'), 'file')
    load(fullfile(datapath, 'gapwidths.mat'))
end


bin_edges = -10:5:50;
noisy_dist = zeros(size(eval_files,1), length(bin_edges)-1);
tip_dist = zeros(size(eval_files,1), length(bin_edges)-1);


f = figure(1);
clf(f);
f.Units = 'points';
f.Position = [200 200 500 320];

ax1 = axes();
ax1.Position = [0.11 0.13 0.85 0.83];
hold(ax1, 'on');

for i = 1:size(eval_files, 1)
   InputData = load(fullfile(eval_files(i).folder, eval_files(i).name)); 
   vid_file = fullfile(eval_files(i).folder, [eval_files(i).name(1:end-15) '_compiled.mat']);
   TrackData = load(vid_file);
   
   idx = find(strcmp(gapwidth, [eval_files(i).name(1:end-15) '.dat']));
   
   gwidth = gapwidth{idx,2};
   gapinfo = TrackData.Annotations.Tracker.gapinfo;
   npix = abs(gapinfo.edge_1 - gapinfo.edge_2);
   mm_per_pix = gwidth/npix; % mm px^-1
   
   dist_to_target = TrackData.Annotations.Tracker.dist_nose_target.*mm_per_pix; % dist to target in mm
   
   
   for j = 1:length(bin_edges)-1
       idx = find(dist_to_target >= bin_edges(j) & dist_to_target < bin_edges(j+1));
       
       counttip = 0;
       countnoisy = 0;
       total = 0;
       for k = 1:length(idx)
           if ~isempty(TrackData.Annotations.Tracker.Traces{idx(k)})
               counttip = counttip + numel(find( InputData.Eval_traces.Tracker_missing_tip{idx(k)} ));
               countnoisy = countnoisy + numel(find( InputData.Eval_traces.Tracker_noisy_tip{idx(k)}));
               total = total+ numel(find( InputData.Eval_traces.Tracker{idx(k)} == 2));
           end
       end
       
       if total > 0
           noisy_dist(i,j) = countnoisy/total;
           tip_dist(i,j) = counttip/total;
       else
           noisy_dist(i,j) = 0;
           tip_dist(i,j) =0 ;
       end
   end
   
   plot(ax1, bin_edges(1:end-1)+0.5*mean(diff(bin_edges)), noisy_dist(i,:),'color',[0 0 1 0.5])
   plot(ax1, bin_edges(1:end-1)+0.5*mean(diff(bin_edges)), tip_dist(i,:),'color',[1 0 0 0.5])

   
end

p1 = plot(ax1, bin_edges(1:end-1)+0.5*mean(diff(bin_edges)), mean(noisy_dist,1), 'color', 'b','LineWidth',2);
p2 = plot(ax1, bin_edges(1:end-1)+0.5*mean(diff(bin_edges)), mean(tip_dist,1), 'color', 'r','LineWidth', 2);

xlabel(ax1, 'Distance to target (mm)');
ylabel(ax1, 'Artefacts (%)');


legend(ax1, [p1, p2], {'Noisy Endings','Missing Tip'}, 'Location','northeast')








