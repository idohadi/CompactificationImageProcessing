%% Docs
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2021
% ***********************************************************

%% RNG seed
rng(0, 'twister');

%% File setup
dt = datetime;
dt = datestr(dt, 'yyyy-mm-dd-HHMM');
fnNOEXT = ['Classification-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
snr = 1;
maxTranslation = [0, 2.5, 5, 7.5, 10];
trialNo = 10;
sampleSize = 5*10^3;
classesNo = 7;
classProb = ones(classesNo, 1)/classesNo;
imageSize = 101;

%% Rotation and translation algorithm setup
printBegEndMsg('Rotation and translation algorithm setup', true);
bandlimit = 30;
loadCGTable(bandlimit);
tDesign = loadtd(2*bandlimit);
interval = cos(pi/4)*[-1, 1];
scalingParam = 1;
wpass = 0.05;

clear bispectrum;
clear bispectrum_mex;
clear buildU;
clear buildU;
clear buildK;
clear buildK_mex;

save(fn, 'bandlimit', 'interval', 'scalingParam', 'wpass');

% Computing denoising matrix
K = buildBispectrumDebiasingMatrix(imageSize, bandlimit, tDesign, ...
    interval, scalingParam);

printBegEndMsg('Rotation and translation algorithm setup', false);

%% Rotation only algorithm setup
printBegEndMsg('Rotation only algorithm setup', true);

truncation = 2; 

save(fn, 'truncation', '-append');

printBegEndMsg('Rotation only algorithm setup', false);

%% Generate the dataset
printBegEndMsg('Dataset generation', true);
genFunc = @() imageByUpsampling(14, 16, imageSize, round(0.2*imageSize));

% Generate class representatives and compute the noise variance
printBegEndMsg('Generating class representatives.', true);
classRepresentatives = zeros([imageSize, imageSize, classesNo]);
classRepresentativesMeans = zeros(1, classesNo);
for J=1:classesNo
    classRepresentatives(:, :, J) = genFunc();
    classRepresentativesMeans(J) = norm(classRepresentatives(:, :, J), 'fro')^2/imageSize^2;
end
signal = mean(classRepresentativesMeans);
sigma = sqrt(signal/snr);

save(fn, 'classRepresentatives', 'sigma', 'genFunc', 'signal', ...
    'classRepresentativesMeans', '-append');
printBegEndMsg('Generating class representatives.', false);

% Generate image paramters
printBegEndMsg('Genearting image paramters.', true);
classMembership = randsample(1:classesNo, sampleSize, true, classProb);
rotations = 360*rand(1, sampleSize);
translationAngle = 2*pi*rand(1, sampleSize);
translationSize = rand(1, sampleSize);

save(fn, 'classMembership', 'rotations', 'translationAngle', ...
    'translationSize', '-append');
printBegEndMsg('Genearting image paramters.', false);

% Generate images (no noise, no translation, only rotations applied)
printBegEndMsg('Genearting images.', true);
% dataset = zeros([imageSize, imageSize, sampleSize]);
dataset = classRepresentatives(:, :, classMembership);
parfor J=1:sampleSize
    dataset(:, :, J) = imrotate(dataset(:, :, J), rotations(J), ...
        'bicubic', 'crop');
end

save(fn, 'dataset', '-append');
printBegEndMsg('Genearting images.', false);

printBegEndMsg('Dataset generation', false);

%% Run test
% Result structure
results = struct();

for mt=1:length(maxTranslation)
    printBegEndMsg(num2str([mt, length(maxTranslation), maxTranslation(mt)], ...
        'Max translation size %d of %d. max translation size = %.2f'), true);
    
    % Applying translations to images, if needed
    printBegEndMsg('Applying translations to images.', true);
    if maxTranslation(mt)>2*eps
        maxT = maxTranslation(mt);
        parfor J=1:sampleSize
            translation = maxT*translationSize(J)...
                *[cos(translationAngle(J)), sin(translationAngle(J))];
            dataset(:, :, J) = imtranslate(dataset(:, :, J), translation, ...
                'cubic', 'OutputView', 'same');
        end
    end
    printBegEndMsg('Applying translations to images.', false);
    
    for J=1:trialNo
        % Print output
        printBegEndMsg(num2str([J, trialNo], 'Trial %d of %d.'), true);
        
        % Add noise
        printBegEndMsg('Adding noise.', true);
        noisyDataset = dataset + sigma*randn(size(dataset));
        printBegEndMsg('Adding noise.', false);
            
        % Run rotation and translation algorithm
        printBegEndMsg('Rotation and translation algorithm.', true);
        t = tic;
        [GRotAndTrans, DRotAndTrans] ...
            = classifyImages(noisyDataset, bandlimit, 'K', K, ...
            'wpass', wpass, 'sigma2', sigma^2, 'lowRank', 400, ...
            'scalingParam', scalingParam);
        classificationRotAndTransRuntime = toc(t);
        printBegEndMsg('Rotation and translation algorithm.', false);
        
        % Run rotation only algorithm
        printBegEndMsg('Rotation only algorithm.', true);
        t = tic;
        [GRotOnly, DRotOnly] ...
            = classifyImagesRotInv(noisyDataset, truncation, 'wpass', wpass);
        classificationRotOnlyRuntime = toc(t);
        printBegEndMsg('Rotation only algorithm.', false);
        
        % Compute node specificity
        printBegEndMsg('Node specificity computation.', true);
        t = tic;
        classSpecificityRotAndTrans = zeros(sampleSize, 1);
        classSpecificityRotOnly = zeros(sampleSize, 1);
        for I=1:sampleSize
            sameClassRotAndTrans ...
                = classMembership(logical(GRotAndTrans(I, :)))==classMembership(I);
            classSpecificityRotAndTrans(I) ...
                = sum(sameClassRotAndTrans)/sum(GRotAndTrans(I, :));
            
            sameClassRotOnly ...
                = classMembership(logical(GRotOnly(I, :)))==classMembership(I);
            classSpecificityRotOnly(I) ...
                = sum(sameClassRotOnly)/sum(GRotOnly(I, :));
        end
        measureRuntime = toc(t);
        printBegEndMsg('Node specificity computation.', false);

        % Saving data
        results(mt, J).GRotAndTrans = GRotAndTrans;
        results(mt, J).DRotAndTrans = DRotAndTrans;
        results(mt, J).GRotOnly = GRotOnly;
        results(mt, J).DRotOnly = DRotOnly;
        results(mt, J).classSpecificityRotAndTrans = classSpecificityRotAndTrans;
        results(mt, J).classSpecificityRotAndTrans95th = prctile(classSpecificityRotAndTrans, 95);
        results(mt, J).classSpecificityRotOnly = classSpecificityRotOnly;
        results(mt, J).classSpecificityRotOnly95th = prctile(classSpecificityRotOnly, 95);
        results(mt, J).classificationRotAndTransRuntime = classificationRotAndTransRuntime;
        results(mt, J).classificationRotOnlyRuntime = classificationRotOnlyRuntime;
        results(mt, J).measureRuntime = measureRuntime;
        
        % Save intermediate results to file
        save(fn, 'results', '-append');
        % Print output
        printBegEndMsg(num2str([J, trialNo], 'Trial %d of %d.'), true);
        
        strOutput = num2str([results(mt, J).classSpecificityRotAndTrans95th, ...
            results(mt, J).classificationRotAndTransRuntime/60, ...
            results(mt, J).classSpecificityRotOnly95th, ...
            results(mt, J).classificationRotOnlyRuntime/60], ...
            ['Results saved to disk. Summary statistics: \\n', ...
            '\t\t\t95th prtle\t\tRuntime (min)\\n', ...
            '\tRot & trans \t\t %.3f \t\t\t%.3f\\n', ...
            '\tRot only \t\t %.3f \t\t\t%.3f\\n']);
        fprintf(strOutput);
    end
    
    printBegEndMsg(num2str([mt, length(maxTranslation), maxTranslation(mt)], ...
        'Max translation size %d of %d. max translation size = %.2f'), false);
end

% Save result
save(fn, 'results', '-append');

%% Produce figure
fig = figure;

bins = 0:0.025:1;

pos = [1000, 523, 694, 455];

classSpecificityRotAndTrans = zeros([size(results), sampleSize]);
classSpecificityRotOnly = zeros([size(results), sampleSize]);
for mt=1:length(maxTranslation)
    for J=1:trialNo
        classSpecificityRotAndTrans(mt, J, :) = results(mt, J).classSpecificityRotAndTrans;
        classSpecificityRotOnly(mt, J, :) = results(mt, J).classSpecificityRotOnly;
    end
end

classSpecificityRotAndTrans = squeeze(mean(classSpecificityRotAndTrans, 2));
classSpecificityRotOnly = squeeze(mean(classSpecificityRotOnly, 2));

t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Distribution of Mean Node Score');

fa = 0.5;
ea = 0.6;
for mt=1:length(maxTranslation)
    nexttile;
    hold;
    rotOnlyH{mt} = histogram(classSpecificityRotOnly(mt, :), bins, ...
        'FaceAlpha', fa, 'EdgeAlpha', ea);
    rotAndTransH{mt} = histogram(classSpecificityRotAndTrans(mt, :), bins, ...
        'FaceAlpha', fa, 'EdgeAlpha', ea);
    set(gca, 'yscale', 'log');
    title(num2str([maxTranslation(mt)], 'Max translation size = %.1f'));
    if mt==length(maxTranslation)
    legend({'Rotation only', 'Rotation and translation'}, ...
            'Orientation', 'vertical', ...
            'Location', 'northwest');
    end
    hold off;
    
    ylim([10^0, 10^4]);
    xlim([0, 1]);
    xticks(0:0.2:1);
    xlabel('Mean node score');
end

% savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
