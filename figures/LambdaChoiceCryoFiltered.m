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
t = tiledlayout(3, 5, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'Relationship of Smoothing and Bandlimit');

nexttile([1, 4]);
imagesc(cryoBackProjRelErrMean);
colormap('hot');
c = colorbar;
c.Ticks = [10^-2, 10^-1, 10^-0];

xt = 5:5:length(bandlimit);
xticks(xt);
xticklabels(bandlimit(xt));

yt = 2:2:length(sigma);
yticks(yt);
yticklabels(sigma(yt));

title('Cryo-EM images');
xlabel('Bandlimit');
ylabel('\sigma');

set(gca,'ColorScale','log')


N = 1;
imCryoFiltered = zeros(size(imCryo, 1), size(imCryo, 2), length(sigma)+1);
imCryoFiltered(:, :, 1) = imCryo(:, :, N);

for s=1:length(sigma)
    imCryoFiltered(:, :, s+1) = imgaussfilt(imCryoFiltered(:, :, 1), ...
        sigma(s));
end

caxismin = min(imCryoFiltered, [], 'all');
caxismax = max(imCryoFiltered, [], 'all');
caxislims = [caxismin, caxismax];

nexttile;
imagesc(imCryoFiltered(:, :, 1));
colormap('hot');
colorbar;
caxis(caxislims);

title('Original');
xticks([]);
yticks([]);

for J=2:size(imCryoFiltered, 3)
    nexttile;
    
    imagesc(imCryoFiltered(:, :, J));
    colormap('hot');
    colorbar;
    caxis(caxislims);

    title(num2str(sigma(J-1), '\\sigma = %.1f'));
    xticks([]);
    yticks([]);
end

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
