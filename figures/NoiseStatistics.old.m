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
image = zeros(imageSize, imageSize, trialNo, length(sigma));
imagePowSpec = zeros(imageSize, imageSize, trialNo, length(sigma));
shc = zeros((bandlimit+1)^2, trialNo, length(sigma));
shcPowSpec = zeros(bandlimit+1, trialNo, length(sigma));

for s=1:length(sigma)
    printBegEndMsg(num2str([s, length(sigma), sigma(s)^2], ...
        'Variance %d of %d (Sigma^2 = %.3f)'), true);

    % Calculate the images
    printBegEndMsg('Generate white noise images', true);
    image(:, :, :, s) = sigma(s)^2 * randn(imageSize, imageSize, trialNo);
    printBegEndMsg('Generate white noise images', false);
    
    % Calculate power spectrum of images
    printBegEndMsg('Calculate power spectrum of images', true);
    imagePowSpec = fft2d(image, [1, 2]);
    imagePowSpec = abs(imagePowSpec).^2;
    printBegEndMsg('Calculate power spectrum of images', false);
    
    % Project the images onto the sphere and compute their power spectrum
    printBegEndMsg('Project images to sphere and compute power spectrum', true);
    parfor J=1:size(image, 3)
        shc(:, J, s) = image2shc(image(:, :, J, s), bandlimit, tDesign, ...
            interval, scalingParam, sh);
        shcPowSpec(:, J, s) = powerSpectrum(shc(:, J, s), bandlimit);
    end
    printBegEndMsg('Project images to sphere and compute power spectrum', false);

    printBegEndMsg(num2str([s, length(sigma), sigma(s)^2], ...
        'Variance %d of %d (Sigma^2 = %.3f)'), false);
end

% Save result
%save(fn, 'image', 'imagePowSpec', 'shc', 'shcPowSpec', '-append');

%% Generate the mean power spectrum
printBegEndMsg('Calculating the mean power spectrums', true);
shcAbsMean = squeeze(mean(abs(shc).^2, 2));
shcPowSpecMean = squeeze(mean(shcPowSpec, 2));
imagePowSpecMean = squeeze(mean(imagePowSpec, 3));

% Save result
save(fn, 'shcAbsMean', 'shcPowSpecMean', 'imagePowSpecMean', '-append');

printBegEndMsg('Calculating the mean power spectrums', false);

%% Produce figure
fig = figure;
t = tiledlayout(1, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

shcPowSpecMeanNormalized = shcPowSpecMean./(2*(0:bandlimit)'+1);
shcNorms = zeros(size(shcAbsMean, 1), 1);
for l=0:bandlimit
    for m=-l:l
        shcNorms(l^2 + l + m +1) = l;
    end
end
shcAbsMeanNormalized = shcAbsMean./shcNorms;
rescell = {shcAbsMean, shcPowSpecMeanNormalized, imagePowSpecMean};
caxislims = [minCell(rescell), maxCell(rescell)];

nexttile;



savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;

%% Utility functions
function Y = fft2d(X, dims)
for J=dims
    Y = fft(X, [], J)/size(X, J);
end
end

function maxX = maxCell(C)
maxX = 0;
for J=1:numel(C)
    tmp2 = max(C{J}, [], 'all');
    maxX = max(maxX, tmp2);
end
end

function minX = minCell(C)
minX = 0;
for J=1:numel(C)
    tmp2 = min(C{J}, [], 'all');
    minX = min(minX, tmp2);
end
end