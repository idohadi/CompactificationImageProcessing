%% Figure
% Dependencies
%   SphericalInvariantAnalysis  commit bf9324d6969949c37b3b0e6d9d2412b92a3494a3 (HEAD -> ibispWeight)
% 
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% Created   19/10/2020, 20:11
% ***********************************************************

%% DOCS: Experimental procedure
% 
% Load data saved from code2020_09_24_0500.mat:
% * results
% * N
% * groundTruth
% 
% Choose a sample size (using an index of N). Of the 50 repeats of this 
% sample size, take the one with maximum relative error.
% 
% For the chosen repeat, repeat 15 times:
% 1. Generate a direction, a random complex vector with norm 1.
% 2. Use said direction to generate an initial bispetrum inversion guess by 
%       x0 = shc + initRelError*norm(shc)*direction
%    such that x0 has initRelError relative error (in l2-norm).
% 3. Invert the estimator of the bispectrum of the chosen repeat and sample
%    size.
% 5. Align the resulting SHC to the ground truth SHC.
% 6. Estimate the image and measure the relative error.
% 

%% CODE: Load data
% Data from code2020_09_24_0500.mat
data = load('code2020_09_24_0500.mat');
N = data.N;
previousResults = data.results;
groundTruth = data.groundTruth;

im = groundTruth.im;
shc = groundTruth.shc;
ps = groundTruth.ps;
bisp = groundTruth.bisp;
% t-design
td = loadtd(2*bandlimit+2);

% Constants
n = 19; 
sampleSize = N(n);
w = [1; 2];
initRelError = 0.6;

% Chosing a repeat
[maxSHCRelErr, trial] = max([previousResults(19, :).shcRelError]);
maxImageRelErr = previousResults(19, trial).imageRelError;
bispEst = previousResults(n, trial).bispEst;
psEst = previousResults(n, trial).psEst;

%% CODE: Generate data
tAll = tic;
trials = 15;
results = struct(); % results of this experiment
for J=1:trials
    % Initial guess
    direction = randSHC(bandlimit);
    direction = direction/norm(direction);
    x0 = shc + norm(shc)*initRelError*direction;

    [invertedSHC, rootedResidual, invOutput] = ibispectrum(bispEst, psEst, bandlimit, w, x0);
    t = tic;
    [relativeDistance, alignedSHC2, optimalRotation, alOutput] = alignSHC(groundTruth.shc, invertedSHC, bandlimit, td);
    t = toc(t);
    
    % Save data
    results(J).shcRelDist = relativeDistance;
    results(J).alignedSHC = alignedSHC2;
    results(J).alOutput = alOutput;
    results(J).optimalRotation = optimalRotation;
    results(J).invertedSHC = invertedSHC;
    results(J).rootedResidual = rootedResidual;
    results(J).invFirstOrdOpt = invOutput.firstorderopt;
    results(J).direction = direction;
    results(J).x0 = x0;
    
    % To image
    shc2 = results(J).alignedSHC;
    imageEst ...
        = shc2image(shc2, bandlimit, imageSize, ...
                    interval, scalingParam);
    imageEst = real(imageEst);

    imageRelError = norm(im - imageEst, 'fro')/norm(im, 'fro');
    imageBackProjection = shc2image(shc, bandlimit, imageSize, ...
                                    interval, scalingParam);
    
    % Save data
    results(J).imageRelDistClean = imageRelError ...
        - norm(im - imageBackProjection)/norm(im, 'fro');
    results(J).imageRelDist = imageRelError;
    results(J).imageEst = imageEst;
    
    disp([num2str(J), ' of ', num2str(trials), ' completed. shc rel d = ', ...
        num2str(relativeDistance), ', im rel d = ', ...
        num2str(results(J).imageRelDistClean)]);
end

[minSHCRelErr, minSHCRelErrInd] = min([results.shcRelDist]);
[minImageRelErr, minImageRelErrInd] = min([results.imageRelDistClean]);
tAll = toc(tAll);

save('code2020_10_19_2030.mat', 'results', 'im', 'shc', 'ps', 'bisp', ...
    'sampleSize', 'w', 'maxSHCRelErr', 'trial', 'psEst', 'bispEst', ...
    'initRelError', 'minSHCRelErr', 'minSHCRelErrInd', ...
    'minImageRelErr', 'minImageRelErrInd', 'tAll');

%% CODE: Plotting
fig = figure; 
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact'); 

% Final results of the new test
nexttile; 

y = [maxSHCRelErr, minSHCRelErr, maxImageRelErr, minImageRelErr];
x = 1:numel(y);
labels = {'Original SHC', 'Best new SHC', 'Original image', 'Best new image'};

bar(x, y); 
title(['Original relative error vs. new (sample size ', ...
    num2str(sampleSize, '%.0e'), ', trial ', num2str(trial), ')']); 
xticks(x);
xticklabels(labels);
ylabel('Relative error');

% SHC/image relative error for all repeats of current tests vs. inversion 
% residual
nexttile;

xSHC = [results.shcRelDist];
xImage = [results.imageRelDistClean];
y = [results.rootedResidual];

hold;
scatter(xSHC, y);
scatter(xImage, y);
hold off;

title('Relative error per repeat vs. inversion residual');
xlabel('Relative error');
ylabel('Rooted residual');
legend({'SHC', 'Image'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'vertical');

% SHC/image relative error for all repeats of current tests vs. inversion
% first order optimality
nexttile;
y = [results.invFirstOrdOpt];

hold;
scatter(xSHC, y);
scatter(xImage, y);
hold off;

title('Relative error per repeat vs. first order optimality');
xlabel('Relative error');
ylabel('Rooted residual');
legend({'SHC', 'Image'}, ...
    'Location', 'southoutside', ...
    'Orientation', 'vertical');
