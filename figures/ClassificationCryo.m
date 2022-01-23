%% Docs
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2022
% ***********************************************************

%% Pool setting
poolSize = 48;
try
    parpool(poolSize);
catch
    
end

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

snr = 1;
maxTranslation = [0, 2.5, 5, 7.5, 10];
sampleSize = 10^4;
imageSize = 101;
classNo = 100;
classProb = ones(classNo, 1)/classNo;
k = 50;

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
for J=1:sampleSize
    dataset(:, :, J) = imrotate( ...
        classRepresentatives(:, :, classMembership(J)), ...
        rotations(J), ...
        'bicubic', 'crop');
end
printBegEndMsg('Generate image dataset', false);

save(fn, 'rotations', 'classMembership', 'translationAngle', ...
    'translationSize', '-append');

printBegEndMsg('Dataset generation', false);

%% Run test
printBegEndMsg('Running test', true);

results = struct();

% Compute the length of the bispectrum vector
[shc, sh] = image2shc(dataset(:, :, 1), bandlimit, tDesign, interval, ...
    scalingParam);
b = bispectrum(shc, bandlimit, CGs);
blen = size(b, 1);

bispectra = zeros(blen, sampleSize);

for mt=1:length(maxTranslation)
    maxT = maxTranslation(mt);
    printBegEndMsg(num2str([mt, length(maxTranslation), maxT], ...
        'Translation %d of %d (max trans = %.3f)'), true);
    
    printBegEndMsg('Translating image', true);
    noisyDataset = zeros(size(dataset));
    for J=1:size(dataset, 3)
        translation = maxT*translationSize(J)...
            *[cos(translationAngle(J)), sin(translationAngle(J))];
        noisyDataset(:, :, J) = imtranslate(dataset(:, :, J), translation, ...
            'cubic', 'OutputView', 'same');
    end
    printBegEndMsg('Translating image', false);
    
    
    printBegEndMsg('Adding noise', true);
    noisyDataset = noisyDataset + sigma*randn(size(noisyDataset));
    printBegEndMsg('Adding noise', false);
    
    
    
    printBegEndMsg('Rotation and translation invariant algorithm', true);
    t = tic;
    
    printBegEndMsg('Calculating bispectrum', true);
    for J=1:size(noisyDataset, 3)
        shc = image2shc(noisyDataset(:, :, J), bandlimit, tDesign, interval, ...
            scalingParam, sh);
        bispectra(:, J) = bispectrum(shc, bandlimit, CGs) - sigma^2*K*cSHC2rSHC(shc);
    end
    printBegEndMsg('Calculating bispectrum', false);


    printBegEndMsg('Calculating distance matrix', true);
    distanceMatrix = pdist(bispectra');
    distanceMatrix = squareform(distanceMatrix);
    printBegEndMsg('Calculating distance matrix', false);
    
    printBegEndMsg('Finding k nearest neighrbors', true);
    [B, I] = mink(distanceMatrix, k+1);
    printBegEndMsg('Finding k nearest neighrbors', false);
    
    printBegEndMsg('Calculating node specificity', true);
    classSpecificityRotAndTrans = zeros(sampleSize, 1);
    for J=1:sampleSize
        sameClassNo = classMembership(I(:, J))==classMembership(J);
        classSpecificityRotAndTrans(J) = (sum(sameClassNo)-1)/k;
    end
    printBegEndMsg('Calculating node specificity', false);

    
    bispRuntime = toc(t);
    results(mt).maxT = maxT;
    results(mt).distanceMatrix = distanceMatrix;
    results(mt).B = B;
    results(mt).I = I;
    results(mt).classSpecificityRotAndTrans = classSpecificityRotAndTrans;
    results(mt).bispRuntime = bispRuntime;
    
    printBegEndMsg('Rotation and translation invariant algorithm', false);
    
    
    
    printBegEndMsg('Rotation invariant algorithm', true);
    t = tic;
    
    printBegEndMsg('Perform sPCA', true);
    %sPCA_data =  data_sPCA(noisyDataset, sigma^2);
    sPCA_data = sPCA_PSWF(noisyDataset, sigma^2);
    printBegEndMsg('Perform sPCA', false);

    printBegEndMsg('Initial classification', true);
    [class, class_refl, rot, corr] ...
        = Initial_classification_FD_update(sPCA_data, k);
    printBegEndMsg('Initial classification', false);

    printBegEndMsg('VDM classification', true);
    flag = 0;
    num_eig = 24;
    [class_VDM, class_VDM_refl, angle] = ...
        VDM(class, corr, rot, class_refl, k, flag, num_eig, k);

    rotRuntime = toc(t);
    
    I = class_VDM.';
    
    printBegEndMsg('Calculating node specificity', true);
    classSpecificityRotOnly = zeros(sampleSize, 1);
    for J=1:sampleSize
        sameClassNo = classMembership(I(:, J))==classMembership(J);
        classSpecificityRotOnly(J) = (sum(sameClassNo))/k;
    end
    printBegEndMsg('Calculating node specificity', false);

    
    results(mt).sPCA_data = sPCA_data;
    results(mt).class = class;
    results(mt).class_refl = class_refl;
    results(mt).rot = rot;
    results(mt).corr = corr;
    results(mt).flag = flag;
    results(mt).num_eig = num_eig;
    results(mt).class_VDM = class_VDM;
    results(mt).class_VDM_refl = class_VDM_refl;
    results(mt).angle = angle;
    results(mt).classSpecificityRotOnly = classSpecificityRotOnly;

    printBegEndMsg('Rotation invariant algorithm', false);
    
    
    printBegEndMsg(num2str([mt, length(maxTranslation), maxT], ...
        'Translation %d of %d (max trans = %.3f)'), false);
end

save(fn, 'results', '-append');

printBegEndMsg('Running test', false);


%% Produce figure
fig = figure;

bins = 0:0.025:1;

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Distribution of Node Score');

fa = 0.5;
ea = 0.6;
for mt=1:length(maxTranslation)
    nexttile;
    hold;
    histogram(results(mt).classSpecificityRotOnly, bins, ...
        'FaceAlpha', fa, 'EdgeAlpha', ea);
    histogram(results(mt).classSpecificityRotAndTrans, bins, ...
        'FaceAlpha', fa, 'EdgeAlpha', ea);
    set(gca, 'yscale', 'log');
    
    title(num2str(maxTranslation(mt), 'Max translation size = %.1f'));
    xlabel('Mean node score');
    
    ylim([10^0, sampleSize]);
    xlim([0, 1]);
    xticks(0:0.2:1);
    
    if mt==length(maxTranslation)
    legend({'Rotation only', 'Rotation and translation'}, ...
            'Orientation', 'vertical', ...
            'Location', 'northwest');
    end
    hold off;
    
end

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
