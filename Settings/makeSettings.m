%% Introduction
% This file contains settings for raw whisker tracking. Tweaking its
% parameters might improve the tracker results. The performance of the
% current settings can be checked with the GUI controlled by
% ParameterSetup.m. Any settings updated by the GUI will be saved in
% 'Settings.mat', which is loaded into this script as 'costum'.


%#### Interface Settings ##################################################

% Load Manual updated settings, these override settings native in this
% script. Delete 'Settings.mat' or prevent file loading here to reset the
% Settings.
if exist('Settings.mat')
    costum = load('Settings.mat');
end



%#### Path Details ########################################################
% Set the file extension, the GUI's will only recognize files with this
% extension:
Settings.video_extension = '.dat';

% regexp expression for name format (can be left empty)
Settings.format = 'M(?<MOUSE>\d+)_R(?<SESSION>\d+)_(?<TRIAL>\d+).dat';

% Set a path to which the GUI's will automatically direct:
Settings.default_video_path = 'E:\Studie\Stage Neurobiologie\Videos\Mouse 47';

% As explained in ParameterSetup, some dataformats (those not readable by
% VideoReader) will need external metadata to load the video. Make sure the
% LoadFrame function will handle any costum format properly
Settings.use_external_specfile = 1; 

% The background detection calls costumBackground.m, a plugin that adds
% additional operations to extract a background. 
Settings.costum_background = 1;

% The Object detection function can be called in single or batch mode, the
% first will display a GUI the second only return the extracted background.
Settings.batch_mode = 0; % don't touch this one                                    
                                    
       
%#### Export details ##################################################
% In the current build all data is stored in a single path; the one
% specified below:
Settings.outpath = 'E:\Studie\Stage Neurobiologie\Videos\VideoDatabase\Tracker Performance';


%#### Parameters Tracking #################################################
% Objects are detected by estimating the % of frames wherein a pixel value
% falls below a threshold:
Settings.object_pixel_ratio = 0.05;
Settings.object_threshold = 0.7; 

% The video is sampled at n points:
Settings.n_background_samples = 30; % number of sample frames to extract background

% Tracking the nose is optional, thracking whiskers does not depend on it
Settings.TrackNose = 1;

% Allow nose only in the gap (turn off in non-gap based task)
Settings.nose_only_in_gap = 1;

% When the gap door opens, it is sometimes tracked as nose. It can be
% filtered by limiting the maximum speed allowed for a nose
Settings.max_nose_speed = 5; % px/frame

% Nose position is tracked every nth frame, as specified in this variable:
Settings.nose_interval = 5; 

% Frames to track can be filtered using nose tracking. If nose trackng is
% requried, set to 'NOSE_REQUIRED', if the nose has to be within the gap
% set to, 'NOSE_INGAP'. If all frames should be tracked, set to 'DEFAULT'
Settings.frame_select = 'ManualDataRequired';

% Minimum number of consecutive frames with nose
Settings.min_consc_frames = 50;

% In whisker tracking a silhouette is extracted based on a threshold, to
% generate a ROI for tracking seeds.If ParameterSetup doesnt return proper
% tracking seeds, adjusting this trheshold will improve tracked seeds;
Settings.Silhouettethreshold = 0.3; % Also used in whisker tracking

% Seeds are found on an expansion of the rodent silhouette in a frame,
% increasing dilation size will track seeds further away from the fur
Settings.Dilationsize = 20;

% Sensitivity towards detecting Trace Seeds, increasing this threshold will
% descrease the number of seeds found, slightly improving tracking speed.
Settings.Origin_threshold = 0.05;


% Trace propagation iteratevely creates a ROI to detect a new point on a
% trace. The ROI is characterized as an arc of size circle_end-cirlce_start
% degrees, under an angle which is derived from previous trace steps:
Settings.circle_start = -25;
Settings.circle_end = 25;

% The distance of an ROI w.r.t. previous tracked point is set as:
Settings.stepsize = 5;

% When encountering an object during tracing, the algortihm will try to
% continue tracking behind the object in an area at distance
% stepsize*extrapolationsize:
Settings.extrapolationsize = 7;


% During tracking each now point is validated against a few criteria,
% including a minimum distance from the edge, a minimum trace length and a
% signal to noise ratio (trace_threshold).
Settings.dist_from_edge = 5;
Settings.minimum_traclength = 8;
Settings.trace_threshold = 0.99; % stop criterium for single trace tracking



%#### Update manual adjusted settings #####################################

% If there are manual updates on the settings, update them
if exist('costum','var')
    Settings.object_threshold = costum.Settings.object_threshold;
    Settings.Dilationsize = costum.Settings.Dilationsize;
    Settings.Origin_threshold = costum.Settings.Origin_threshold;
    Settings.trace_threshold = costum.Settings.trace_threshold;
end

