function Settings = makeAnalyseSettings()
% Filter parameters for frame rejection
Settings.min_nose_dist = 40; % minimum required distance of nose from border
Settings.frame_density = 0.75; % percentage of neighbouring frames with traces
Settings.max_accel = 5; % Maximum acceleration of nose (px/frameinterval)
Settings.max_distance = 120; % Maximum distance from estimated midpoint, to be included

Settings.max_gap_fill = 20; % gaps of up to 20 frames are marked as valid




% Filter parameters for trace rejection
Settings.max_dist_trace_nose = 150; % Only inlude traces close to nose
Settings.min_trace_length = 5;


% Fitting settings
Settings.fit_degree_small = 1;
Settings.fit_degree_medium = 2;
Settings.fit_degree_large = 3;


Settings.n_spline_samples = 30;



Settings.cutoff = 3; % cutoff length on both endings of trace before fitting

% Touch detection settings
Settings.dist_object = 3;

% Only detect touches around target platform (not home platform)
Settings.target_area_only = 1;

% Set number of rows around target area that are allowed for touch
% detection
Settings.objects_nrows = 5;



