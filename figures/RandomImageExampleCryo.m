%% Docs
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2022
% ***********************************************************
%% RNG seed
rng(0, 'twister');

%% File setup
dt = datetime;
dt = datestr(dt, 'yyyy-mm-dd-HHMM');
fnNOEXT = ['RandomImageExampleCryo-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

imageSize = 101;
imageNo = 8;

save(fn, 'imageSize', 'imageNo', 'fnNOEXT');

printBegEndMsg('Setup parameters', false);

%% Generate images
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol_true = cryo_downsample(f.volref, imageSize*ones(1, 3));

rotation = rand_rots(imageNo);
im = cryo_project(vol_true, rotation, imageSize, 'double');

save(fn, 'im', 'rotation', '-append');

%% Figure of images
perRow = 4;

fig = figure;
t = tiledlayout(ceil(size(im, 3)/perRow), perRow, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

caxislims = [min(im, [], 'all'), max(im, [], 'all')];

for J=1:size(im, 3)
    nexttile;
    imagesc(im(:, :, J));
    colormap('hot');
    colorbar;
    caxis(caxislims);
    xticks([]);
    yticks([]);
end

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
