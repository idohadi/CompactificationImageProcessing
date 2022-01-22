 %% Docs
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2022
% ***********************************************************

%% Pool setting
poolSize = 32;
try
    parpool(poolSize);
catch
    
end

%% RNG seed
rng(0, 'twister');

%% File setup
dt = datetime;
dt = datestr(dt, 'yyyy-mm-dd-HHMM');
fnNOEXT = ['MRACryo-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

snr = 1;
maxTranslation = 5;
sampleSize = 10^4;
imageSize = 101;

bandlimit = 70;
loadCGTable(bandlimit);
global CGs;
tDesign = loadtd(2*bandlimit);
interval = cos(pi/4)*[-1, 1];
scalingParam = 1;

clear bispectrum;
clear bispectrum_mex;
clear buildU;
clear buildU;
clear buildK;
clear buildK_mex;

save(fn, 'snr', 'maxTranslation', 'sampleSize', 'imageSize', ...
    'bandlimit', 'interval', 'scalingParam', 'fnNOEXT', '-v7.3');

% Computing denoising matrix
K = buildBispectrumDebiasingMatrix(imageSize, bandlimit, tDesign, ...
    interval, scalingParam);

save(fn, 'K', '-append');

printBegEndMsg('Setup parameters', false);

%% Generate the dataset
printBegEndMsg('Dataset generation', true);

% Image by cryo-EM simulation
printBegEndMsg('Generate representatives', true);
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol_true = cryo_downsample(f.volref, imageSize*ones(1, 3));

rotation = rand_rots(1);
im = cryo_project(vol_true, rotation, imageSize, 'double');

save(fn, 'rotation', 'im', '-append');
printBegEndMsg('Generate representatives', false);

% Generate images
printBegEndMsg('Generate image dataset', true);
rotations = 360*rand(1, sampleSize);
translationAngle = 2*pi*rand(1, sampleSize);
translationSize = rand(1, sampleSize);

dataset = zeros(imageSize, imageSize, sampleSize);
for J=1:sampleSize
    dataset(:, :, J) = imrotate(im, rotations(J), 'bicubic', 'crop');
    translation = maxTranslation*translationSize(J)...
        *[cos(translationAngle(J)), sin(translationAngle(J))];
    dataset(:, :, J) = imtranslate(dataset(:, :, J), translation, ...
        'cubic', 'OutputView', 'same');
end
save(fn, 'translationAngle', 'translationSize', '-append');
printBegEndMsg('Generate image dataset', false);


% Generate noise
printBegEndMsg('Add noise to dataset', true);
signal = norm(im, 'fro')^2/imageSize^2;
sigma = sqrt(signal/snr); % variance of noise
noise = sigma*randn(size(dataset));
noisyDataset = dataset + noise;

save(fn, 'noise', '-append');
clear noise;
clear dataset;
printBegEndMsg('Add noise to dataset', false);

printBegEndMsg('Dataset generation', false);

%% Run test
printBegEndMsg('Running test', true);

% Compute the length of the bispectrum vector
[shc, sh] = image2shc(noisyDataset(:, :, 1), bandlimit, tDesign, interval, ...
    scalingParam);
b = bispectrum(shc, bandlimit, CGs);
blen = size(b, 1);

bispectra = zeros(blen, sampleSize);

printBegEndMsg('Calculating bispectrum', true);
for J=1:size(noisyDataset, 3)
    shc = image2shc(noisyDataset(:, :, J), bandlimit, tDesign, interval, ...
        scalingParam, sh);
    bispectra(:, J) = bispectrum(shc, bandlimit, CGs) - sigma^2*K*cSHC2rSHC(shc);
end
save(fn, 'bispectra', '-append');
printBegEndMsg('Calculating bispectrum', false);


printBegEndMsg('Estimating bispectrum', true);
bispEstimator = mean(bispectra, 2);
save(fn, 'bispEstimator', '-append');
printBegEndMsg('Estimating bispectrum', false);


printBegEndMsg('Inversion', true);
printBegEndMsg('Computing initial guess.', true);
[FBsPCA_data, z,zcost,info, Timing, imager, image0, ZA] ...
    = hMRA_uniform(noisyDataset, im, sigma);
x0 = image2shc(im, bandlimit, tDesign, interval, scalingParam, sh);
initialRelError = norm(im - imager, 'fro')/norm(im, 'fro');
initialAbsError = norm(im - imager, 'fro');
save(fn, 'imager', 'x0', 'initialRelError', 'initialAbsError', '-append');

printBegEndMsg(num2str(initialRelError, ...
    'Computing initial guess.\n\tInit guess rel err = %.3e.'), ...
    false);

try
    parpool(poolSize);
catch
    
end

printBegEndMsg('Invert the bispectrum.', true);
[shcEstimator, rootedResidual, inversionOutput] = ibispectrum(bispEstimator, bandlimit, x0);
save(fn, 'shcEstimator', 'rootedResidual', 'inversionOutput', '-append');
printBegEndMsg('Invert the bispectrum.', false);

printBegEndMsg('Align inverted bispectrum and estimate SHC.', true);
[shcRelError, alignedSHCEst, optimalRotation, alignmentOutput] ...
    = alignSHC(shc, shcEstimator, bandlimit, tDesign);
save(fn, 'shcRelError', 'alignedSHCEst', 'optimalRotation', ...
    'alignmentOutput', '-append');
printBegEndMsg(num2str(shcRelError, ...
    'Align inverted bispectrum and estimate SHC.\n\tSHC rel err = %.3e.'), false);

printBegEndMsg('Inversion', false);


printBegEndMsg('Image estimation.', true);
imageEst ...
    = shc2image(alignedSHCEst, bandlimit, imageSize, ...
                interval, scalingParam);
imageEst = real(imageEst);
imageRelError = norm(im - imageEst, 'fro')/norm(im, 'fro');
imageBackProjection = shc2image(shc, bandlimit, imageSize, ...
                                interval, scalingParam);
imageRelErrorCleaned = imageRelError ...
    - norm(im - imageBackProjection)/norm(im, 'fro');
printBegEndMsg(num2str([imageRelError, imageRelErrorCleaned], ...
    'Image estimation.\n\tIm rel err = %.3e.\n\tIm rel err cln = %.3e.'), false);
save(fn, 'imageEst', 'imageRelError', 'imageRelErrorCleaned', '-append');
printBegEndMsg('Running test', true);

%% Produce figure
fig = figure;



savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
