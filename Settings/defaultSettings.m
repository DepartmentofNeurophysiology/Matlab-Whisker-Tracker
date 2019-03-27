% Generate default settings
%% Default settings
% these are not updated by parametersetup and regard data and processing
% type

% Extension (for automated video detection in directories)
Settings.video_extension = '.avi';

% Use parallel processing
Settings.use_parfor = 1;

% Format (regexp string of filename format, optional, can be left empty
Settings.format = '';
%Settings.format = '(?<date>\d+)_M(?<mouse>\d+)_R(?<trial>\d+)_(?<session>\d+)_(?<t2>\d+).mat';
% Track nose
Settings.track_nose = 1;

% Frame select (see 'CostumFrameSelect')
Settings.frame_select = 'DEFAULT';

% ROI arc size
Settings.cirlce_start = -25;
Settings.circle_end = 25;
Settings.stepsize = 5;
Settings.extrapollationsize = 7;
Settings.dist_from_edge = 5;




%% Tracker settings
% these are updated using 'ParameterSetup'
Settings.Gaussian_kernel_size = 3;
Settings.doGaussian = 1;
Settings.Gamma =  1;
Settings.Background_threshold = 0.5;
Settings.Edges_kernel_large = 10;
Settings.Edges_kernel_small = 6;
Settings.Edges_threshold = 0.35;
Settings.Shape_threshold = 0.21;
Settings.Dilation = 5;
Settings.Seed_threshold = 0.05;
Settings.Trace_threshold = 0.49;
Settings.Trace_kernel_large = 5;
Settings.Trace_kernel_small = 1;
Settings.FullRoi = 1;

