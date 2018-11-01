function Evaluation = printEval(varargin)
%% printEval(varargin)
%(optional) name -values:
%   'path' - absolute path to directory with '*Evaluation.mat' files
%   'verbose' - toggle print function (0 - none, 1 - total results, 2 - per
%               video)
p = inputParser;
p.CaseSensitive = 0;
defaultPath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
defaultVerbose = 0;


addParameter(p,'Path', defaultPath);
addParameter(p,'Verbose',defaultVerbose);

parse(p, varargin{:});


datapath = p.Results.Path;
Verbose = p.Results.Verbose;

eval_files = dir(fullfile(datapath, '*_Evaluation.mat'));

vididx = [];
Correct = [];
Tip_missing = [];
Noisy_ending = [];

% append all evaluation results
for i = 1:size(eval_files, 1)
   vididx(end+1) = length(Correct) + 1;
    
   load(fullfile(eval_files(i).folder, eval_files(i).name)) 
   
   nframes = size(Eval_traces.Tracker, 2);
   for j = 1:nframes
      if isempty(Eval_traces.Tracker{j})
          continue
      end
      
      Correct = [Correct, Eval_traces.Tracker{j}-1]; %#ok<*AGROW>
      Tip_missing = [Tip_missing, Eval_traces.Tracker_missing_tip{j}];
      Noisy_ending = [Noisy_ending, Eval_traces.Tracker_noisy_tip{j}];
   end
    
end

idx = find(Correct & Tip_missing == 0 & Noisy_ending == 0);
Correct_no_artifact = zeros(length(Correct), 1);
Correct_no_artifact(idx) = 1; %#ok<FNDSB>

n_traps = 100;


if Verbose == 1 || Verbose == 2    
    fprintf('TRACKER PERFORMANCE\n')
    fprintf('Running %d bootstraps on %d samples...\n', n_traps, length(Correct));
end


Correct_btstrp = bootstrp(n_traps, @mean, Correct);
Correct_no_artifact_btstrp =  bootstrp(n_traps, @mean, Correct_no_artifact);

if Verbose == 1 || Verbose == 2
    fprintf('%% Correct traces (from all traces): %1.4f, std: %1.4f\n', mean(Correct_btstrp), std(Correct_btstrp))
    fprintf('%% Correct no artifact (from all traces): %1.4f, std: %1.4f\n', mean(Correct_no_artifact_btstrp), std(Correct_no_artifact_btstrp))
end

Evaluation.Correct = Correct;
Evaluation.Correct_btstrp = Correct_btstrp;
Evaluation.Correct_no_artifact = Correct_no_artifact;
Evaluation.Correct_no_artifact_btstrp = Correct_no_artifact_btstrp;





















