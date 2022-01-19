%% RNG seed
rng(0, 'twister');

%% File setup
dt = datetime;
dt = datestr(dt, 'yyyy-mm-dd-HHMM');
fnNOEXT = ['LambdaChoiceCryoFiltered-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Load images
printBegEndMsg('Load saved parameters and images', true);

load('LambdaChoice-2022-01-17-0003.mat', 'imCryo', ...
    'bandlimit', 'imageNo', 'imageSize', 'interval')

lambda = 1;

sigma = 0.5:0.5:5;

save(fn, 'imCryo', 'lambda', 'sigma', 'bandlimit', 'imageNo', ...
    'imageSize', 'interval', 'fnNOEXT');

printBegEndMsg('Load saved parameters and images', false);

%% Run test
printBegEndMsg('Running test', true);

cryoBackProjErr = zeros(length(sigma), length(bandlimit), imageNo);
cryoBackProjRelErr = zeros(length(sigma), length(bandlimit), imageNo);

for b=1:length(bandlimit)
    printBegEndMsg(num2str([b, length(bandlimit), bandlimit(b)], ...
        'Bandlimit %d of %d (banlimit = %d)'), true);

    bl = bandlimit(b);
    tDesign = loadtd(2*bl);
    
    for s=1:length(sigma)
        printBegEndMsg(num2str([s, length(sigma), sigma(s)], ...
            'Sigma %d of %d (sigma = %.3f)'), true);
    
        scalingParam = lambda;
        [shc, sh] = image2shc(randn(imageSize, imageSize), bl, ...
            tDesign, interval, scalingParam);
        [~, sh2] = shc2image(shc, bl, imageSize, interval, scalingParam);
        
        ims = imgaussfilt(imCryo, sigma(s));
        parfor t=1:imageNo
            % Cryo
            cryoSHC = image2shc(ims(:, :, t), bl, ...
                tDesign, interval, scalingParam, sh);
            cryoIm = shc2image(cryoSHC, bl, imageSize, interval, scalingParam, sh2);
            cryoBackProjErr(s, b, t) = norm(cryoIm - ims(:, :, t), 'fro');
            cryoBackProjRelErr(s, b, t) = cryoBackProjErr(s, b, t)/norm(ims(:, :, t), 'fro');
        end
        
        printBegEndMsg(num2str([s, length(sigma), sigma(s)], ...
            'Sigma %d of %d (sigma = %.3f)'), false);
    end
    
    printBegEndMsg(num2str([b, length(bandlimit), bandlimit(b)], ...
        'Bandlimit %d of %d (banlimit = %d)'), false);
end

% Save result
save(fn, 'cryoBackProjErr', ...
    'cryoBackProjRelErr', ...
    '-append');

printBegEndMsg('Running test', false);


%% Generate the mean for every image
printBegEndMsg('Calculating the mean over all images', true);

cryoBackProjErrMean = squeeze(mean(cryoBackProjErr, 3));
cryoBackProjRelErrMean = squeeze(mean(cryoBackProjRelErr, 3));

% Save result
save(fn, 'cryoBackProjErrMean', ...
    'cryoBackProjRelErrMean', ...
    '-append');

printBegEndMsg('Calculating the mean over all images', false);

%% Produce figure
fig = figure;
t = tiledlayout(1, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Relationship of Scaling Parameter and Bandlimit');

nexttile;
imagesc(cryoBackProjRelErrMean);
colormap('hot');
colorbar;

xt = 5:5:length(bandlimit);
xticks(xt);
xticklabels(bandlimit(xt));

yt = 2:2:length(sigma);
yticks(yt);
yticklabels(sigma(yt));

title('Cryo-EM images');
xlabel('Bandlimit');
ylabel('Sigma');

set(gca,'ColorScale','log')

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
