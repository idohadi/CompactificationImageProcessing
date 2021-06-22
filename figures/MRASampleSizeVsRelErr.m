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
fnNOEXT = ['MRASampleSizeVsRelErr-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file
save(fn, 'dt', 'fnNOEXT');

%% Setup parameters
clear bispectrum;
clear bispectrum_mex;
clear buildU;
clear buildU;
clear buildK;
clear buildK_mex;

bandlimit = 16;
scalingParam = 1;
interval = cos(pi/4)*[-1, 1];
snr = 0.1;
imageSize = 101;

N = [10^2*(1:9), 10^3*(1:9), 10^4*(1:9)]; % Sample sizes
repNo = 15; % Repeatition per sample size

tDesign = loadtd(2*bandlimit);
loadCGTable(bandlimit);
global CGs;

% Computing denoising matrix
K = buildBispectrumDebiasingMatrix(imageSize, bandlimit, ...
    tDesign, interval, scalingParam);

save(fn, 'bandlimit', 'scalingParam', 'interval', 'snr', 'N', 'repNo', ...
    'imageSize', 'K', '-append');

%% Generate the ground truth
printBegEndMsg('Ground truth image generation.', true);
% Random image
im = imageByUpsampling(14, 16, imageSize, round(0.2*imageSize));

% Calculate bispectrum
[shc, sh] = image2shc(im, bandlimit, tDesign, interval, scalingParam);
bisp = bispectrum(shc, bandlimit, CGs);

% Calculate variance of noise
signal = norm(im, 'fro')^2/imageSize^2;
sigma2 = signal/snr; % variance of noise

save(fn, 'im', 'shc', 'sh', 'bisp', 'signal', 'sigma2', '-append');

printBegEndMsg('Ground truth image generation.', false);

printBegEndMsg('Generate denoised dataset.', true);
denoisedData = zeros([imageSize, imageSize, max(N)]);
parfor J=1:max(N)
    denoisedData(:, :, J) = imrotate(im, 360*rand(), 'bicubic', 'crop');
    denoisedData(:, :, J) = imtranslate(denoisedData(:, :, J) , ...
        randTranslation(5), ...
        'bicubic', 'OutputView', 'same');
end
printBegEndMsg('Generate denoised dataset.', false);

%% Tester main body
printBegEndMsg('Tester main body.', true);

results = struct();

% Initialize estimators
shcEst = zeros(length(shc), length(N));
bispEst = zeros(length(bisp), length(N));

for r=1:repNo
    printBegEndMsg(num2str([r, repNo], 'Rep %d of %d.'), true);
    
    % Generate dataset with noise
    printBegEndMsg('Generate noisy dataset.', true);
    data = denoisedData + sqrt(sigma2)*randn(size(denoisedData));
    printBegEndMsg('Generate noisy dataset.', false);
    
    % Loop over sample sizes
    for n=1:length(N)
        printBegEndMsg(num2str([n, length(N), N(n)], ...
            'Sample size %d of %d (N = %.3e).'), true);
        
        % Bispectrum estimation
        printBegEndMsg('Bispectrum estimation.', true);
        if n==1
            bispEst(:, n) = estimateBispectrum(false, ...
                data(:, :, 1:N(n)), bandlimit, tDesign, ...
                interval, scalingParam, K, sigma2, sh);
        else
            bispEst(:, n) = estimateBispectrum(false, ...
                data(:, :, N(n-1)+1:N(n)), bandlimit, tDesign, ...
                interval, scalingParam, K, sigma2, sh);
            bispEst(:, n) = ((N(n) - N(n-1))/N(n))*bispEst(:, n) ...
                + (N(n-1)/N(n))*bispEst(:, n-1);
        end
        
        % Save results
        results(n, r).bispEst = bispEst(:, n);
        results(n, r).bispEstAbsError = norm(bispEst(:, n) - bisp);
        results(n, r).bispEstRelError = norm(bispEst(:, n) - bisp)/norm(bisp);
        
        save(fn, 'results', 'r', 'n', '-append');
        printBegEndMsg(num2str(results(n, r).bispEstRelError, ...
            'Bispectrum estimation.\n\tBisp rel err = %.3e.'), ...
            false);
        
        
        % Bispectrum inversion initialization
        printBegEndMsg('Computing initial guess.', true);
        [FBsPCA_data, z,zcost,info, Timing, imager, image0, ZA] ...
            = hMRA_uniform(data(:, :, 1:N(n)), im, sqrt(sigma2));
        x0 = image2shc(imager, bandlimit, tDesign, interval, scalingParam, sh);
%         x0 = shc;

        % Save results
        results(n, r).imager = imager;
        results(n, r).x0 = x0;
        results(n, r).initialRelError = norm(im - imager, 'fro')/norm(im, 'fro');
        results(n, r).initialAbsError = norm(im - imager, 'fro');
        
        save(fn, 'results', 'r', 'n', '-append');
        printBegEndMsg(num2str(results(n, r).initialRelError, ...
            'Computing initial guess.\n\tInit guess rel err = %.3e.'), ...
            false);
        
        
        % Bispectrum inversion
        printBegEndMsg('Invert the bispectrum.', true);
        [invertedSHC, rootedResidual, output] ...
            = ibispectrum(bispEst(:, n), bandlimit, x0);
        
        % Save results
        results(n, r).invertedSHC = invertedSHC;
        results(n, r).rootedResidual = rootedResidual;
        results(n, r).inversionOutput = output;
        
        save(fn, 'results', 'r', 'n', '-append');
        printBegEndMsg(num2str([rootedResidual, output.firstorderopt], ...
            'Invert the bispectrum.\n\tRooted res = %.3e. \n\t1st ord opt = %.3e.'), false);
        
        
        % SHC alignment and estimation
        printBegEndMsg('Align inverted bispectrum and estimate SHC.', true);
        [shcRelError, alignedSHCEst, optimalRotation, output] ...
            = alignSHC(shc, invertedSHC, bandlimit, tDesign);
        
        % Save results
        results(n, r).shcRelError = shcRelError;
        results(n, r).alignedSHCEst = alignedSHCEst;
        results(n, r).optimalRotation = optimalRotation;
        results(n, r).alignmentOutput = output;
        
        printBegEndMsg(num2str(shcRelError, ...
            'Align inverted bispectrum and estimate SHC.\n\tSHC rel err = %.3e.'), false);
        
        % Image estimation
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
        
        % Save results
        results(n, r).imageEst = imageEst;
        results(n, r).imageRelError = imageRelError;
        results(n, r).imageBackProjection = imageBackProjection;
        results(n, r).imageRelErrorCleaned = imageRelErrorCleaned;
        
        save(fn, 'results', 'r', 'n', '-append');
        printBegEndMsg(num2str([imageRelError, imageRelErrorCleaned], ...
            'Image estimation.\n\tIm rel err = %.3e.\n\tIm rel err cln = %.3e.'), false);
        
        printBegEndMsg(num2str([n, length(N), N(n)], ...
            'Sample size %d of %d (N = %.3e).'), false);
    end
    
    printBegEndMsg(num2str([r, repNo], 'Rep %d of %d.'), false);
end

printBegEndMsg('Tester main body.', false);

%% Plot the results
% Data collation
estimators = zeros([size(results), 3]);
% Bispectrum estimator error
estimators(:, :, 1) = reshape([results.bispEstRelError], size(results));
% Initial guess error
estimators(:, :, 2) = reshape([results.initialRelError], size(results));
% Image estimator error
estimators(:, :, 3) = reshape([results.imageRelErrorCleaned], size(results));

% Average estimators over tries
meanEstimators = squeeze(mean(estimators, 2));

% Find slope of log-log sequence
slope = zeros(size(meanEstimators, 2), 2);
for J=1:size(meanEstimators, 2)
    slope(J, :) = polyfit(log(N), log(meanEstimators(:, J)), 1);
end
slope = slope(:, 1);

% Setup figure
fig = figure;

mrSz = 5; % marker size
basicLabels = {'Bispectrum', 'Initial guess', 'Image'};
legendLabels = cell(size(basicLabels));

hold;
for J=1:size(meanEstimators, 2)
    plot(N, meanEstimators(:, J), ...
        'Marker', 'o', 'MarkerSize', mrSz);
    
    slopeStr = num2str(slope(J), '%.3f');
    legendLabels{J} = [basicLabels{J}, ' (', slopeStr, ')'];
end

imageBackProjection = shc2image(shc, bandlimit, imageSize, ...
                                        interval, scalingParam);
imageRelErrClean = norm(im - imageBackProjection, 'fro')/norm(im, 'fro');
yline(imageRelErrClean, '--');
legendLabels{end+1} = 'Back-projection bound';
hold off;

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

xlabel('Sample size');
ylabel('Relative error');

legend(legendLabels, 'Location', 'southwest');

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
