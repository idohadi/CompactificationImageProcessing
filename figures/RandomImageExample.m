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
fnNOEXT = ['RandomImageExample-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

imageSize = 101;
imageNo = 8;

genFunc = @() imageByUpsampling(14, 16, imageSize, round(0.2*imageSize));

save(fn, 'imageSize', 'imageNo', 'fnNOEXT');

printBegEndMsg('Setup parameters', false);

%% Generate images
im = zeros(imageSize, imageSize, imageNo);
for J=1:imageNo
    im(:, :, J) = genFunc();
end

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
