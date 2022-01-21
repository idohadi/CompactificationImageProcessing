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
fnNOEXT = ['BispectrumInvarianceCryo-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

save(fn, 'fnNOEXT');

%% Image parameters
bandlimit = 50;
tDesign = loadtd(2*bandlimit);
interval = cos(pi/4)*[-1, 1];
imageSize = 101;

%% Scaling paramter choice
scalingParam = 1;
save(fn, 'scalingParam', '-append');

%% Image generation
printBegEndMsg('Image generation.', true);
rotation = rand_rots(imageNo);
im = cryo_project(vol_true, rotation, imageSize, 'double');
[shc, sh] = image2shc(im, bandlimit, tDesign, interval, scalingParam);
loadCGTable(bandlimit);
global CGs;
bisp = bispectrum(shc, bandlimit, CGs);
save(fn, 'im', 'shc', 'bisp', '-append');
printBegEndMsg('Image generation.', false);

%% Rotation only
printBegEndMsg('Rotations only test.', true);

% One angle rotation grid
rotationGrid = 1:359;

% Measure outcome
rotOnlyBispDist = zeros(length(rotationGrid), 1);
for R=rotationGrid
    % Rotate original image
    rotatedIm = imrotate(im, R, 'bicubic', 'crop');
    rotatedSHC = image2shc(rotatedIm, bandlimit, tDesign, interval, ...
        scalingParam, sh);
    % Calculate rotated image bispectrum and compare to original
    transformedBisp = bispectrum(rotatedSHC, bandlimit, CGs);
    rotOnlyBispDist(R) = norm(transformedBisp - bisp);
end
save(fn, 'rotOnlyBispDist', '-append');

printBegEndMsg('Rotations only test.', false);

%% Translation only
printBegEndMsg('Translation only test.', true);

% Translation sampling paramters
maxTranslationSize = 15;
translationSizeGrid = 0:0.5:maxTranslationSize;
transAnglesNo = 10;

% Measure outcome
transOnlyBispDist = zeros(length(translationSizeGrid), transAnglesNo);
for S=1:length(translationSizeGrid)
    % Fix translation size
    ts = translationSizeGrid(S);
    % Sample rotation angles
    rotAnglesSamples = 360*rand(transAnglesNo, 1);
    
    printBegEndMsg(num2str([S, length(translationSizeGrid), translationSizeGrid(S)], 'Translation size %d of %d (size = %.3f).'), true);
    for R=1:length(rotAnglesSamples)
        % Build translation vector
        trans = ts * [cos(rotAnglesSamples(R)), sin(rotAnglesSamples(R))];
        % Translate image by trans
        translatedIm = imtranslate(im, trans, 'cubic', 'OutputView', 'same');
        translatedSHC = image2shc(translatedIm, bandlimit, tDesign, ...
            interval, scalingParam, sh);
        % Calculate translated image bispectrum and compare to original
        transformedBisp = bispectrum(translatedSHC, bandlimit, CGs);
        transOnlyBispDist(S, R) = norm(transformedBisp - bisp);
    end
    printBegEndMsg(num2str([S, length(translationSizeGrid), translationSizeGrid(S)], 'Translation size %d of %d (size = %.3f).'), false);
end
save(fn, 'transOnlyBispDist', '-append');

printBegEndMsg('Translation only test.', false);


%% Rotation and translation
printBegEndMsg('Rotation and translation test.', true);
% Translation and rotation sampling paramters
maxTranslationSize = 15;
translationSizeGrid = 0:0.5:maxTranslationSize;
sampleNo = 15; % Number of samples from [0, 2*pi]^2

% Measure outcome
rotAndTransBispDist = zeros(length(translationSizeGrid), sampleNo);
for S=1:length(translationSizeGrid)
    printBegEndMsg(num2str([S, length(translationSizeGrid), translationSizeGrid(S)], 'Translation size %d of %d (size = %.3f).'), true);

    % Fix translation size
    ts = translationSizeGrid(S);
    % Sample angles
    rotAnglesSamples = 360*rand(sampleNo, 1); % In degrees
    translationAnglesSamples = 2*pi*rand(sampleNo, 1); % In radians
    for R=1:length(rotAnglesSamples)
        % Rotate the original image
        transformedIm = imrotate(im, rotAnglesSamples(R), 'bicubic', 'crop');
        % Build translation vector
        trans = ts * [cos(translationAnglesSamples(R)), sin(translationAnglesSamples(R))];
        % Translate image by trans
        transformedIm = imtranslate(transformedIm, trans, 'cubic', 'OutputView', 'same');
        transformedSHC = image2shc(transformedIm, bandlimit, tDesign, ...
            interval, scalingParam, sh);
        % Calculate translated image bispectrum and compare to original
        transformedBisp = bispectrum(transformedSHC, bandlimit, CGs);
        rotAndTransBispDist(S, R) = norm(transformedBisp - bisp);
    end
    
    printBegEndMsg(num2str([S, length(translationSizeGrid), translationSizeGrid(S)], 'Translation size %d of %d (size = %.3f).'), false);
end
save(fn, 'rotAndTransBispDist', '-append');

printBegEndMsg('Rotation and translation test.', false);

%% Plot the results
rotationGrid = 1:359;
maxTranslationSize = 15;
translationSizeGrid = 0:0.5:maxTranslationSize;

fig = figure;
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

yLim = [0, 0.13];

% Rotation only 
nexttile;
normFactor = norm(bisp);

plot(rotationGrid, rotOnlyBispDist/normFactor);

xlabel('Rotation (angle)');
ylabel('Relative error');
title('Rotation only');

yticks(yLim(1):0.01:yLim(end));
xlim([0, 360]);
ylim(yLim);

% Translation only
nexttile;
hold;
normFactor = norm(bisp);

polX = [translationSizeGrid, flip(translationSizeGrid)];
prctl95 = [prctile(transOnlyBispDist, 97.5, 2)/normFactor, ...
    flip(prctile(transOnlyBispDist, 2.5, 2))/normFactor];
c = [17, 17, 17]/255;
fill(polX(:), prctl95(:), c, 'FaceAlpha', 0.3, 'LineStyle', 'none');

y = mean(transOnlyBispDist, 2);
plot(translationSizeGrid, y/normFactor, 'Color', [0 0.4470 0.7410]);
hold off;

xlabel('Translation size (pixels)');
ylabel('Relative distance');
title('Translation only');

yticks(yLim(1):0.01:yLim(end));
xlim([0, 15]);
ylim(yLim);

% Rotation and translation
nexttile;
hold;
normFactor = norm(bisp);

polX = [translationSizeGrid, flip(translationSizeGrid)];
prctl95 = [prctile(rotAndTransBispDist, 97.5, 2)/normFactor, ...
    flip(prctile(rotAndTransBispDist, 2.5, 2))/normFactor];
c = [17, 17, 17]/255;
fill(polX(:), prctl95(:), c, 'FaceAlpha', 0.3, 'LineStyle', 'none');

y = mean(rotAndTransBispDist, 2);
plot(translationSizeGrid, y/normFactor, 'Color', [0 0.4470 0.7410]);
hold off;

xlabel('Translation size (pixels)');
ylabel('Relative distance');
title('Rotation and translation');

yticks(yLim(1):0.01:yLim(end));
xlim([0, 15]);
ylim(yLim);

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
