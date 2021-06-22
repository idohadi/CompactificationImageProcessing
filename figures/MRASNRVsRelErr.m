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
fnNOEXT = ['MRASNRVsRelErr-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file
save(fn);

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
snr = 1./[1, 2, 3, 4, 5, 10, 25, 50, 100];
imageSize = 101;

N = 10^4; % Sample sizes
repNo = 15; % Repeatition per sample size

tDesign = loadtd(2*bandlimit+2);
loadCGTable(bandlimit);
global CGs;

% Computing denoising matrix
K = buildBispectrumDebiasingMatrix(imageSize, bandlimit, tDesign, ...
    interval, scalingParam);

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
sigma2 = signal./snr; % variance of noise

save(fn, 'im', 'shc', 'sh', 'bisp', 'signal', 'sigma2', '-append');

printBegEndMsg('Ground truth image generation.', false);

%% Tester main body
printBegEndMsg('Tester main body.', true);

results = struct();

% Initialize estimators
shcEst = zeros(length(shc), length(sigma2));
bispEst = zeros(length(bisp), length(sigma2));

for r=1:repNo
    printBegEndMsg(num2str([r, repNo], 'Rep %d of %d.'), true);
    
    % Loop over SNRs
    for s=1:length(sigma2)
        printBegEndMsg(num2str([s, length(snr), snr(s)], ...
            'SNR %d of %d (SNR = %.3f).'), true);
        
        
        % Generate dataset, without noise
        printBegEndMsg('Generate dataset without noise.', true);
        data = zeros([imageSize, imageSize, N]);
        parfor J=1:N
            data(:, :, J) = imrotate(im, 360*rand(), 'bicubic', 'crop');
            data(:, :, J) = imtranslate(data(:, :, J) , ...
                randTranslation(5), ...
                'bicubic', 'OutputView', 'same');
        end
        printBegEndMsg('Generate dataset without noise.', false);
        
        % Add noise to the data
        printBegEndMsg('Add noise to dataset.', true);
        data = data + sqrt(sigma2(s))*randn(size(data));
        printBegEndMsg('Add noise to dataset.', false);
        
        % Bispectrum estimation
        printBegEndMsg('Bispectrum estimation.', true);
        bispEst(:, s) = estimateBispectrum(false, data, bandlimit, ...
            tDesign, interval, scalingParam, K, sigma2(s), sh);
        
        % Save results
        results(s, r).bispEst = bispEst(:, s);
        results(s, r).bispEstAbsError = norm(bispEst(:, s) - bisp);
        results(s, r).bispEstRelError = norm(bispEst(:, s) - bisp)/norm(bisp);
        
        save(fn, 'results', 'r', 's', '-append');
        printBegEndMsg(num2str(results(s, r).bispEstRelError, ...
            'Bispectrum estimation.\n\tBisp rel err = %.3e.'), ...
            false);
        
        
        % Bispectrum inversion initialization
        printBegEndMsg('Computing initial guess.', true);
        [FBsPCA_data, z,zcost,info, Timing, imager, image0, ZA] ...
            = hMRA_uniform(data, im, sqrt(sigma2(s)));
        x0 = image2shc(imager, bandlimit, tDesign, interval, scalingParam, sh);
        
        % Save results
        results(s, r).imager = imager;
        results(s, r).x0 = x0;
        results(s, r).initialRelError = norm(im - imager, 'fro')/norm(im, 'fro');
        results(s, r).initialAbsError = norm(im - imager, 'fro');
        
        save(fn, 'results', 'r', 's', '-append');
        printBegEndMsg(num2str(results(s, r).initialRelError, ...
            'Computing initial guess.\n\tInit guess rel err = %.3e.'), ...
            false);
        
        
        % Bispectrum inversion
        printBegEndMsg('Invert the bispectrum.', true);
        [invertedSHC, rootedResidual, output] ...
            = ibispectrum(bispEst(:, s), bandlimit, x0);
        
        % Save results
        results(s, r).invertedSHC = invertedSHC;
        results(s, r).rootedResidual = rootedResidual;
        results(s, r).inversionOutput = output;
        
        save(fn, 'results', 'r', 's', '-append');
        printBegEndMsg(num2str([rootedResidual, output.firstorderopt], ...
            'Invert the bispectrum.\n\tRooted res = %.3e. \n\t1st ord opt = %.3e.'), false);
        
        
        % SHC alignment and estimation
        printBegEndMsg('Align inverted bispectrum and estimate SHC.', true);
        [shcRelError, alignedSHCEst, optimalRotation, output] ...
            = alignSHC(shc, invertedSHC, bandlimit, tDesign);
        
        % Save results
        results(s, r).shcRelError = shcRelError;
        results(s, r).alignedSHCEst = alignedSHCEst;
        results(s, r).optimalRotation = optimalRotation;
        results(s, r).alignmentOutput = output;
        
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
        results(s, r).imageEst = imageEst;
        results(s, r).imageRelError = imageRelError;
        results(s, r).imageBackProjection = imageBackProjection;
        results(s, r).imageRelErrorCleaned = imageRelErrorCleaned;
        
        save(fn, 'results', 'r', 's', '-append');
        printBegEndMsg(num2str([imageRelError, imageRelErrorCleaned], ...
            'Image estimation.\n\tIm rel err = %.3e.\n\tIm rel err cln = %.3e.'), false);
        
        printBegEndMsg(num2str([s, length(snr), snr(s)], ...
            'SNR %d of %d (SNR = %.3f).'), false);
        
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
estimators(:, :, 3) = reshape([results.imageRelError], size(results));

% Average estimators over tries
meanEstimators = squeeze(mean(estimators, 2));

% Find slope of log-log sequence
slope = zeros(size(meanEstimators, 2), 2);
for J=1:size(meanEstimators, 2)
    slope(J, :) = polyfit(log(snr), log(meanEstimators(:, J)), 1);
end
slope = slope(:, 1);

% Setup figure
fig = figure;

mrSz = 5; % marker size
basicLabels = {'Bispectrum', 'Initial guess', 'Image'};
legendLabels = cell(size(basicLabels));

hold;
for J=1:size(meanEstimators, 2)
    plot(snr, meanEstimators(:, J), ...
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

xlabel('SNR');
ylabel('Relative error');

legend(legendLabels, 'Location', 'southwest');

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
