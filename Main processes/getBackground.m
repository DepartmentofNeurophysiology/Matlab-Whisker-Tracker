function Results = getBackground(Settings)
%Results = getBackground(Settings)
% Extract background and edges using background averaging and
% a laplacian of gaussian approach.
%
% Input:
% Settings struct with fields:
% -state (optional, used in parameter setup to prevent video sampling on
%  each run)
% -Nframes  (number of frames in video)
% -Edges_kernel_large (size of box kernel for edge detection)
% -Edges_kernel_small (size of delta kernel for edge detection)
% -Background_threshold (threshold for binarization)
%
% Output:
% Results struct with fields:
% -Edges  (binary matrix marking edges)
% -Objects (binary matrix marking objects)
% -gapinfo (information on gapsize
%%

if isfield(Settings, 'state') && strcmp(Settings.state, 'setup')
    sumFrames = Settings.sumFrames;
    meanFrames = Settings.meanFrames;

else    
    idx = round(linspace(2,Settings.Nframes,100));
    Frames = zeros(512, 640, length(idx));
    for i = 1:length(idx)
        Settings.Current_frame = idx(i)-1;
        Frames(:,:,i)= LoadFrame(Settings);
    end
    sumFrames = sum(Frames,3);
    meanFrames = mean(Frames,3);
end

KL = conv2(meanFrames, ...
    ones(Settings.Edges_kernel_large,Settings.Edges_kernel_large)./...
    Settings.Edges_kernel_large,'same');
KS = conv2(meanFrames, ...
    ones(Settings.Edges_kernel_small, Settings.Edges_kernel_small)./...
    Settings.Edges_kernel_small^2, 'same');

H = KS./KL;
H = H-min(min(H));
H = H./max(max(H));

tax = 0:0.01:1;
for i = 1:length(tax)
    count(i) = numel(find(H < tax(i))); %#ok<AGROW>
end
[~, idx] = max(diff(count));
edge_threshold = tax(idx)-0.02;


Edges = zeros(512, 640);
%Edges(H < Settings.Edges_threshold) = 1;
Edges(H < edge_threshold) = 1;


Background = sumFrames./100;
Background = Background - min(min(Background));
Background = Background ./ max(max(Background));
Background(Background > Settings.Background_threshold) = 0;
Objects = zeros(size(Background));

Objects(Background > 0) = 1;
Results = Costumbackground(Objects, Edges);



