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
fnNOEXT = ['NoiseStatistics-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

sigma = 1;
trialNo = 10^4;

imageSize = 101;

bandlimit = 50;
tDesign = loadtd(2*bandlimit);
interval = cos(pi/4)*[-1, 1];
scalingParam = 1;

save(fn, 'sigma', 'imageSize', 'bandlimit', 'interval', 'scalingParam');

[~, sh] = image2shc(randn(imageSize, imageSize), bandlimit, ...
    tDesign, interval, scalingParam);

printBegEndMsg('Setup parameters', false);

%% Run test
image = zeros(imageSize, imageSize, trialNo);
imagePowSpec = zeros(imageSize, imageSize, trialNo);
shc = zeros((bandlimit+1)^2, trialNo);
shcPowSpec = zeros(bandlimit+1, trialNo);

% Calculate the images
printBegEndMsg('Generate white noise images', true);
image(:, :, :) = sigma^2 * randn(imageSize, imageSize, trialNo);
printBegEndMsg('Generate white noise images', false);

% Calculate power spectrum of images
printBegEndMsg('Calculate power spectrum of images', true);
imagePowSpec = fft2d(image, [1, 2]);
imagePowSpec = abs(imagePowSpec).^2;
printBegEndMsg('Calculate power spectrum of images', false);

% Project the images onto the sphere and compute their power spectrum
printBegEndMsg('Project images to sphere and compute power spectrum', true);
parfor J=1:size(image, 3)
    shc = image2shc(image(:, :, J), bandlimit, tDesign, ...
        interval, scalingParam, sh);
    shcPowSpec(:, J) = powerSpectrum(shc, bandlimit);
end
printBegEndMsg('Project images to sphere and compute power spectrum', false);

%% Generate the mean power spectrum
printBegEndMsg('Calculating the mean power spectrums', true);
shcPowSpecMean = squeeze(mean(shcPowSpec, 2))./(2*(0:bandlimit)'+1);
imagePowSpecMean = squeeze(mean(imagePowSpec, 3));

% Save result
save(fn, 'shcPowSpecMean', 'imagePowSpecMean', '-append');

printBegEndMsg('Calculating the mean power spectrums', false);

%% Produce figure
fig = figure;
t = tiledlayout(1, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Power Spectrum of White Noise');

nexttile;
imagesc(imagePowSpecMean);
colormap('hot');
colorbar;

title('Image');

nexttile;
stem(shcPowSpecMean);

title('Projection');
xlabel('Frequency');
ylabel('Power Spectrum');

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;

%% Utility functions
function Y = fft2d(X, dims)
for J=dims
    Y = fft(X, [], J)/size(X, J);
end
end
