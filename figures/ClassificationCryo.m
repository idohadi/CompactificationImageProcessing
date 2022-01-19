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
fnNOEXT = ['ClassificationCryo-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

sigma = 0.5;
maxTranslation = 10;
sampleSize = 10^4;
imageSize = 101;
classNo = 100;
classProb = ones(classNo, 1)/classNo;
k = 100;

bandlimit = 50;
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

save(fn, 'sigma', 'maxTranslation', 'sampleSize', 'imageSize', ...
    'classNo', 'k', 'bandlimit', 'interval', 'scalingParam', 'fnNOEXT', ...
    'classProb', '-v7.3');

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

rotationsRepresentatives = rand_rots(classNo);
classRepresenatives = cryo_project(vol_true, rotationsRepresentatives, imageSize, 'double');

save(fn, 'rotationsRepresentatives', 'classRepresenatives', '-append');
printBegEndMsg('Generate representatives', false);

% Generate images
printBegEndMsg('Generate image dataset', true);
classMembership = randsample(1:classNo, sampleSize, true, classProb);
rotations = 360*rand(1, sampleSize);
translationAngle = 2*pi*rand(1, sampleSize);
translationSize = rand(1, sampleSize);

dataset = zeros(imageSize, imageSize, sampleSize);
parfor J=1:sampleSize
    dataset(:, :, J) = imrotate( ...
        classRepresenatives(:, :, classMembership(J)), ...
        rotations(J), ...
        'bicubic', 'crop');
    translation = maxTranslation*translationSize(J)...
        *[cos(translationAngle(J)), sin(translationAngle(J))];
    dataset(:, :, J) = imtranslate(dataset(:, :, J), translation, ...
        'cubic', 'OutputView', 'same');
end
printBegEndMsg('Generate image dataset', false);

save(fn, 'translationAngle', 'translationSize', '-append');


% Generate noise
printBegEndMsg('Add noise to dataset', true);
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
parfor J=1:size(noisyDataset, 3)
    shc = image2shc(noisyDataset(:, :, J), bandlimit, tDesign, interval, ...
        scalingParam, sh);
    bispectra(:, J) = bispectrum(shc, bandlimit, CGs) - sigma^2*K*cSHC2rSHC(shc);
end
save(fn, 'bispectra', '-append');
printBegEndMsg('Calculating bispectrum', false);


printBegEndMsg('Calculating distance matrix', true);
distanceMatrix = pdist(bispectra');
distanceMatrix = squareform(distanceMatrix);
save(fn, 'distanceMatrix', '-append');
printBegEndMsg('Calculating distance matrix', false);

printBegEndMsg('Finding k nearest neighrbors', true);
[B, I] = mink(distanceMatrix, k+1);
save(fn, 'B', 'I', '-append');
printBegEndMsg('Finding k nearest neighrbors', false);

printBegEndMsg('Running test', true);

%% Produce figure
fig = figure;


% savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
