% The functions in this path are used to evaluate trackerperformance
%
Tracker_speed = 1;
Tracker_evaluation = 1;
Parameter_evaluation = 1;

%% Tracker speed
% To evaluate tracker processing speed, MWT was run on 10 videos

if Tracker_speed
    Datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
    Videos = dir(fullfile(Datapath, '*.dat'));
    n_to_test = 10;
    idx = randperm(size(Videos,1), n_to_test);
    for i = 1:length(idx)
        Test_Videos{1,i} = fullfile(Videos(idx(i)).folder, Videos(idx(i)).name);
    end
    Timing = timeMWT('files', Test_Videos);
end

%% Trace Evaluation
% Using the function 'verify_data' MWT output was manually evaluated
% The function showEval can be used to inspect the evaluation results
if Tracker_evaluation
    Evaluation = printEval('Verbose',1);
end

% Show dependency of # missing tips & distance to target
%FIG_ERROR_DIST;


%% Parameter estimation



if Parameter_evaluation
    overwrite = 0;
    Datapath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';
    
    load(fullfile(Datapath, 'gapwidths.mat'))
    
    comp_files = dir(fullfile(Datapath, '*compiled.mat'));
    outpath = fullfile(Datapath, 'Parameter estimation');
    if ~exist(outpath, 'file')
        mkdir(outpath)
    end
    
    if exist(fullfile(Datapath, 'Frames_with_nose.mat'), 'file')
        load(fullfile(Datapath, 'Frames_with_nose.mat'))
    else
        keyboard % Make frames_with_nose using 'frame_select'
    end
    
    
    
    % Extract cross correlation
    h = waitbar(0, 'doing things');
    
    if exist(fullfile(outpath, 'correlation.mat'), 'file') & overwrite == 0
        load(fullfile(outpath, 'correlation.mat'))
    else
    
        for i = 1:size(comp_files, 1)

            load(fullfile(comp_files(i).folder, comp_files(i).name))

            vidname = [comp_files(i).name(1:end-13) '.dat'];
            if ~strcmp(vidname, Output(i).Video)
                keyboard % check framerange file
            else
                nframes = size(Annotations.Output.Traces,1);
                frange = zeros(1, nframes);
                pairs = Output(i).Pairs;
                for j = 1:size(pairs, 2)
                    frange(pairs{j}(1):pairs{j}(2)) = 1;
                end

                if numel(find(frange)) < 25
                    Out(i).l_min = NaN;
                    Out(i).l_max = NaN;
                    Out(i).r_min = NaN;
                    Out(i).r_max = NaN;
                    continue
                end

            end




            f = figure(1);
            clf(f);
            f.Units = 'points';
            f.Position = [10 10 1200 800];

            ax3 = axes('Units','points','Position',[100 75 475 300]);
            ax1 = axes('Units','points','Position',[100 425 475 300]);
            ax2 = axes('Units','points','Position',[625 425 475 300]);
            ax4 = axes('Units','points','Position',[625 75 475 300]);

            Angles = getAngles(Annotations);

            if strcmp(Annotations.Tracker.Direction, 'Down')
                Out(i).l_min = compareThetas(Angles.Tracker.l_min_filtered, Angles.Tracker.l_min_peaks, ...
                    Angles.Manual.r_min_filtered, Angles.Manual.r_min_peaks, frange, ax1);
                Out(i).l_max = compareThetas(Angles.Tracker.l_max_filtered, Angles.Tracker.l_max_peaks, ...
                    Angles.Manual.r_max_filtered, Angles.Manual.r_max_peaks, frange,ax2);
                Out(i).r_min = compareThetas(Angles.Tracker.r_min_filtered, Angles.Tracker.r_min_peaks, ...
                    Angles.Manual.l_min_filtered, Angles.Manual.l_min_peaks,frange, ax3);
                Out(i).r_max = compareThetas(Angles.Tracker.r_max_filtered, Angles.Tracker.r_max_peaks, ...
                    Angles.Manual.l_max_filtered, Angles.Manual.l_max_peaks, frange,ax4);
            elseif strcmp(Annotations.Tracker.Direction, 'Up')
                Out(i).l_min = compareThetas(Angles.Tracker.l_min_filtered, Angles.Tracker.l_min_peaks, ...
                    Angles.Manual.l_min_filtered, Angles.Manual.l_min_peaks, frange,ax1);
                Out(i).l_max = compareThetas(Angles.Tracker.l_max_filtered, Angles.Tracker.l_max_peaks, ...
                    Angles.Manual.l_max_filtered, Angles.Manual.l_max_peaks, frange,ax2);
                Out(i).r_min = compareThetas(Angles.Tracker.r_min_filtered, Angles.Tracker.r_min_peaks, ...
                    Angles.Manual.r_min_filtered, Angles.Manual.r_min_peaks,frange, ax3);
                Out(i).r_max = compareThetas(Angles.Tracker.r_max_filtered, Angles.Tracker.r_max_peaks, ...
                    Angles.Manual.r_max_filtered, Angles.Manual.r_max_peaks,frange, ax4);
            end

            saveas(gcf, fullfile(outpath,[comp_files(i).name(1:end-13) '.png']))

            waitbar(i/size(comp_files, 1))

        end
    end
    
    close(h)
    
    Rlmin(1:size(Out, 2)) = NaN;
    Rlmax(1:size(Out, 2)) = NaN;
    Rrmin(1:size(Out, 2)) = NaN;
    Rrmax(1:size(Out, 2)) = NaN;
    
    for i =1:size(Out, 2)
        if isstruct(Out(i).l_min)
            Rlmin(i) = Out(i).l_min.R;
            Rlmax(i) = Out(i).l_max.R;
            Rrmin(i) = Out(i).r_min.R;
            Rrmax(i) = Out(i).r_max.R;
        end        
    end
    
    edges = [0:0.05:1];
    
    fprintf('Correlation between measured thetas (Manual, Tracker):\n')
    fprintf('Left rostral: %1.2f +- %1.2f\n',mean(Rlmin,'omitnan'), std(Rlmin,'omitnan'))
    fprintf('Left caudal : %1.2f +- %1.2f\n',mean(Rlmax,'omitnan'), std(Rlmax,'omitnan'))
    fprintf('Right rostral: %1.2f +- %1.2f\n', mean(Rrmin,'omitnan'), std(Rrmin, 'omitnan'))
    fprintf('Right caudal: %1.2f +- %1.2f\n', mean(Rrmax,'omitnan'), std(Rrmax, 'omitnan'))
    
    figure(1)
    clf
    subplot(2,2,1)
    histogram(Rlmin,edges)
    title('Left rostral')
    xlabel('Correlation coefficient (R)')
    ylabel('Count')
    
    subplot(2,2,2)
    histogram(Rlmax, edges)
    title('Left caudal')
    xlabel('Correlation coefficient (R)')
    ylabel('Count')
    
    subplot(2,2,3)
    histogram(Rrmin, edges)
    title('Right rostral')
    xlabel('Correlation coefficient (R)')
    ylabel('Count')
    
    subplot(2,2,4)
    histogram(Rrmax, edges)
    title('Right caudal')
    xlabel('Correlation coefficient (R)')
    ylabel('Count')
    
    saveas(gcf, fullfile(outpath,['histograms.png']))
    
    save(fullfile(outpath, 'correlation.mat'),'Out')
    
    %}
    
    h = waitbar(0, 'angles');
    for i = 1:size(comp_files, 1)
        load(fullfile(comp_files(i).folder, comp_files(i).name))
        Angles = getAngles(Annotations);
        
        % Tracker
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Tracker.r_min_peaks);   
        Frmin(i) = 1000/median(diff(peaktimes)); %#ok<*SAGROW>
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Tracker.l_min_peaks);
        Flmin(i) = 1000/median(diff(peaktimes));
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Tracker.r_max_peaks);
        Frmax(i) = 1000/median(diff(peaktimes));
        Flmax(i) = 1000/median(diff(peaktimes));
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Tracker.l_max_peaks);
        
        [Slminpro(i), Slminretr(i)] = getSpeed(Angles.Tracker.l_min_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        [Srminpro(i), Srminretr(i)] = getSpeed(Angles.Tracker.r_min_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        [Slmaxpro(i), Slmaxretr(i)] = getSpeed(Angles.Tracker.l_max_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        [Srmaxpro(i), Srmaxretr(i)] = getSpeed(Angles.Tracker.r_max_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);

        
        
        % Manual
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Manual.r_min_peaks);
        ManualFrmin(i) = 1000/median(diff(peaktimes));
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Manual.l_min_peaks);
        ManualFlmin(i) = 1000/median(diff(peaktimes));
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Manual.r_max_peaks);
        ManualFrmax(i) = 1000/median(diff(peaktimes));
        peaktimes = Annotations.Tracker.MetaData.TimeMS( Angles.Manual.l_max_peaks);
        ManualFlmax(i) = 1000/median(diff(peaktimes));
        
        [ManualSlminpro(i), ManualSlminretr(i)] = getSpeed(Angles.Manual.l_min_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        [ManualSrminpro(i), ManualSrminretr(i)] = getSpeed(Angles.Manual.r_min_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        [ManualSlmaxpro(i), ManualSlmaxretr(i)] = getSpeed(Angles.Manual.l_max_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        [ManualSrmaxpro(i), ManualSrmaxretr(i)] = getSpeed(Angles.Manual.r_max_filtered, ...
            Annotations.Tracker.MetaData.TimeMS);
        
        waitbar(i/size( comp_files, 1));
        
        
    end
    close(h)
    
    
    %%
    ymin = min([min(Frmin), min(Flmin), min(Frmax), min(Flmax), min(ManualFrmin), ...
        min(ManualFlmin), min(ManualFrmax), min(ManualFlmax)]);
    ymax =  max([max(Frmin), max(Flmin), max(Frmax), max(Flmax), max(ManualFrmin), ...
        max(ManualFlmin), max(ManualFrmax), max(ManualFlmax)]);
    ymin = floor(ymin*10)/10;
    ymax = ceil(ymax*10)/10;
    
    clc
   fprintf('Measured whisking frequencies:\n')
    fprintf('\tTracker:\n')
    fprintf('\tLeft  rostral: %3.0fHz +- %1.2fHz, Left  caudal: %3.0fHz +- %1.2fHz\n',...
        mean(Flmin,'omitnan'), std(Flmin,'omitnan'), mean(Flmax,'omitnan'), std(Flmax, 'omitnan'))
    fprintf('\tRight rostral: %3.0fHz +- %1.2fHz, Right caudal: %3.0fHz +- %1.2fHz\n',...
        mean(Frmin,'omitnan'), std(Frmin,'omitnan'), mean(Frmax,'omitnan'), std(Frmax, 'omitnan'))
    fprintf('\n\tManual:\n')
    fprintf('\tLeft  rostral: %3.0fHz +- %1.2fHz, Left  caudal: %3.0fHz +- %1.2fHz\n',...
        mean(ManualFlmin,'omitnan'), std(ManualFlmin,'omitnan'), mean(ManualFlmax,'omitnan'), std(ManualFlmax, 'omitnan'))
    fprintf('\tRight rostral: %3.0fHz +- %1.2fHz, Right caudal: %3.0fHz +- %1.2fHz\n',...
        mean(ManualFrmin,'omitnan'), std(ManualFrmin,'omitnan'), mean(ManualFrmax,'omitnan'), std(ManualFrmax, 'omitnan'))   
    
    
    fprintf('\nMeasured whisking speeds (deg/s):\n')
    fprintf('\tTracker protraction:\n')
    fprintf('\tLeft rostral : %4.0f +- %3.0f, Left caudal : %4.0f +- %3.0f\n',...
        median(Slminpro,'omitnan'), std(Slminpro, 'omitnan'), median(Slmaxpro, 'omitnan'), std(Slmaxpro,'omitnan'))
    fprintf('\tRight rostral: %4.0f +- %3.0f, Right caudal: %4.0f +- %3.0f\n',...
        median(Srminpro,'omitnan'), std(Srminpro,'omitnan') , median(Srmaxpro, 'omitnan'), std(Srmaxpro, 'omitnan'))
    fprintf('\n\tManual retraction:\n')
    fprintf('\tLeft rostral : %4.0f +- %3.0f, Left caudal : %4.0f +- %3.0f\n',...
        median(ManualSlminpro,'omitnan'), std(ManualSlminpro, 'omitnan'), median(ManualSlmaxpro, 'omitnan'), std(ManualSlmaxpro,'omitnan'))
    fprintf('\tRight rostral: %4.0f +- %3.0f, Right caudal: %4.0f +- %3.0f\n',...
        median(ManualSrminpro,'omitnan'), std(ManualSrminpro,'omitnan') , median(ManualSrmaxpro, 'omitnan'), std(ManualSrmaxpro, 'omitnan'))
    
    fprintf('\n\tTracker retraction:\n')
     fprintf('\tLeft rostral : %4.0f +- %3.0f, Left caudal : %4.0f +- %3.0f\n',...
        median(Slminretr,'omitnan'), std(Slminretr, 'omitnan'), median(Slmaxretr, 'omitnan'), std(Slmaxretr,'omitnan'))
    fprintf('\tRight rostral: %4.0f +- %3.0f, Right caudal: %4.0f +- %3.0f\n',...
        median(Srminretr,'omitnan'), std(Srminretr,'omitnan') , median(Srmaxretr, 'omitnan'), std(Srmaxretr, 'omitnan'))
    fprintf('\n\tManual protraction:\n')
    fprintf('\tLeft rostral : %4.0f +- %3.0f, Left caudal : %4.0f +- %3.0f\n',...
        median(ManualSlminretr,'omitnan'), std(ManualSlminretr, 'omitnan'), median(ManualSlmaxretr, 'omitnan'), std(ManualSlmaxretr,'omitnan'))
    fprintf('\tRight rostral: %4.0f +- %3.0f, Right caudal: %4.0f +- %3.0f\n',...
        median(ManualSrminretr,'omitnan'), std(ManualSrminretr,'omitnan') , median(ManualSrmaxretr, 'omitnan'), std(ManualSrmaxretr, 'omitnan'))
    
    
    %%
    figure(2)
    clf
    subplot(1, 2, 1)
    hold on
    scatter(1:size(comp_files, 1), Frmin, 'b', 'filled')
    scatter(1:size(comp_files, 1), Flmin, 'r', 'filled')
    scatter(1:size(comp_files, 1), Frmax, 'b')
    scatter(1:size(comp_files, 1), Flmax, 'r')
    ylim([ymin ymax])
    title('Tracker')
    ylabel('\theta')
    xlabel('video')
    
    
        subplot(1, 2, 2)
    hold on
    s1 = scatter(1:size(comp_files, 1), ManualFrmin, 'b', 'filled');
    s2 = scatter(1:size(comp_files, 1), ManualFlmin, 'r', 'filled');
    s3 = scatter(1:size(comp_files, 1), ManualFrmax, 'b');
    s4 = scatter(1:size(comp_files, 1), ManualFlmax, 'r');
    ylim([ymin ymax])
    legend([s1, s2, s3, s4], {'Rmin', 'Lmin', 'Rmax', 'Lmax'},'location','northwest')
    title('Manual')
    ylabel('\theta')
    xlabel('video')
    
    
    yminretr = min([Slminretr, Slmaxretr, Srminretr, Srmaxretr, ...
        ManualSlminretr, ManualSlmaxretr, ManualSrminretr, ManualSrmaxretr]);
    ymaxretr = max([Slminretr, Slmaxretr, Srminretr, Srmaxretr, ...
        ManualSlminretr, ManualSlmaxretr, ManualSrminretr, ManualSrmaxretr]);
    yminpro = min([Slminpro, Slmaxpro, Srminpro, Srmaxpro, ...
        ManualSlminpro, ManualSlmaxpro, ManualSrminpro, ManualSrmaxpro])
   ymaxpro = max([Slminpro, Slmaxpro, Srminpro, Srmaxpro, ...
        ManualSlminpro, ManualSlmaxpro, ManualSrminpro, ManualSrmaxpro])
    
    
    figure(3)
    clf
    subplot(2,2,1)
    hold on
    scatter(1:size(comp_files, 1), Slminpro, 'r', 'filled')
    scatter(1:size(comp_files, 1), Srminpro, 'b', 'filled')
    scatter(1:size(comp_files, 1), Slmaxpro, 'r')
    scatter(1:size(comp_files, 1), Srmaxpro, 'b')
    ylim([yminpro ymaxpro])
    title('Tracker protraction')
    
   subplot(2,2,3)
    hold on
    scatter(1:size(comp_files, 1), Slminretr, 'r', 'filled')
    scatter(1:size(comp_files, 1), Srminretr, 'b', 'filled')
    scatter(1:size(comp_files, 1), Slmaxretr, 'r')
    scatter(1:size(comp_files, 1), Srmaxretr, 'b')
    title('Tracker retraction')
    ylim([yminretr ymaxretr])
    
        subplot(2,2,2)
    hold on
    scatter(1:size(comp_files, 1), ManualSlminpro, 'r', 'filled')
    scatter(1:size(comp_files, 1), ManualSrminpro, 'b', 'filled')
    scatter(1:size(comp_files, 1), ManualSlmaxpro, 'r')
    scatter(1:size(comp_files, 1), ManualSrmaxpro, 'b')
        ylim([yminpro ymaxpro])

    title('Manual protraction')
    
   subplot(2,2,4)
    hold on
    scatter(1:size(comp_files, 1), ManualSlminretr, 'r', 'filled')
    scatter(1:size(comp_files, 1), ManualSrminretr, 'b', 'filled')
    scatter(1:size(comp_files, 1), ManualSlmaxretr, 'r')
    scatter(1:size(comp_files, 1), ManualSrmaxretr, 'b')
    title('Manual retraction')
    ylim([yminretr ymaxretr])

end






