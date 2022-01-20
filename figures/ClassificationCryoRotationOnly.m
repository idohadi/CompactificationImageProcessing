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
fnNOEXT = ['ClassificationCryoRotationOnly-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

snr = 1;
maxTranslation = 5;
sampleSize = 5*10^3;
imageSize = 101;
classNo = 7;
classProb = ones(classNo, 1)/classNo;
k = 50;

save(fn, 'snr', 'maxTranslation', 'sampleSize', 'imageSize', ...
    'classNo', 'k', 'fnNOEXT', ...
    'classProb', '-v7.3');

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
classRepresentatives = cryo_project(vol_true, rotationsRepresentatives, imageSize, 'double');
classRepresentativesMeans = zeros(classNo, 1);
for J=1:classNo
    classRepresentativesMeans(J) = norm(classRepresentatives(:, :, J), 'fro')^2/imageSize^2;
end
signal = mean(classRepresentativesMeans);
sigma = sqrt(signal/snr);

save(fn, 'rotationsRepresentatives', 'classRepresentatives', ...
    'classRepresentativesMeans', 'signal', 'sigma', '-append');
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
        classRepresentatives(:, :, classMembership(J)), ...
        rotations(J), ...
        'bicubic', 'crop');
    translation = maxTranslation*translationSize(J)...
        *[cos(translationAngle(J)), sin(translationAngle(J))];
    dataset(:, :, J) = imtranslate(dataset(:, :, J), translation, ...
        'cubic', 'OutputView', 'same');
end
printBegEndMsg('Generate image dataset', false);

save(fn, 'rotations', 'classMembership', 'translationAngle', ...
    'translationSize', '-append');


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

% nv is "noise variance"
printBegEndMsg('Perform sPCA', true);
[sPCA_data, denoised_images] = sPCA_PSWF(noisyDataset, sigma^2);
save(fn, 'sPCA_data', '-append');
printBegEndMsg('Perform sPCA', false);

printBegEndMsg('Initial classification', true);
[class, class_refl, rot, corr] ...
    = Initial_classification_FD_update(sPCA_data, k);
save(fn, 'class', 'class_refl', 'rot', 'corr', '-append');
printBegEndMsg('Initial classification', false);


printBegEndMsg('VDM classification', true);
flag = 0;
[class_VDM, class_VDM_refl, angle] = ...
    VDM(class, corr, rot, class_refl, k, flag, k);
save(fn, 'flag', 'class_VDM', 'class_VDM_refl', 'angle', '-append');
printBegEndMsg('VDM classification', false);

% %% Calcualting node specificity
% printBegEndMsg('Calculating node specificity', true);
% classSpecificity = zeros(sampleSize, 1);
% for J=1:sampleSize
%     sameClassNo = classMembership(I(:, J))==classMembership(J);
%     classSpecificity(J) = (sum(sameClassNo)-1)/k;
% end
% printBegEndMsg('Calculating node specificity', false);
% 
% %% Produce figure
% fig = figure;
% 
% bins = 0:0.025:1;
% fa = 0.5;
% ea = 0.6;
% 
% histogram(classSpecificity, bins, 'FaceAlpha', fa, 'EdgeAlpha', ea);
% set(gca, 'yscale', 'log');
% 
% ylim([10^0, 10^4]);
% xlim([0, 1]);
% xticks(0:0.2:1);
% xlabel('Node score');
% 
% savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
