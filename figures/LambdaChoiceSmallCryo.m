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
fnNOEXT = ['LambdaChoiceSmallCryo-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

lambda = 1:0.1:5;
bandlimit = 2:2:70;
imageNo = 50;

imageSize = 43;

interval = cos(pi/4)*[-1, 1];
scalingParam = 1;


save(fn, 'lambda', 'bandlimit', 'imageNo', 'imageSize', 'interval', ...
    'scalingParam', 'fnNOEXT');

printBegEndMsg('Setup parameters', false);

%% Image generation functions
printBegEndMsg('Image generation functions', true);

% Image by upsampling
genFunc1 = @(tn) imagebyUps(tn, imageSize);

% Image by cryo-EM simulation
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol_true = cryo_downsample(f.volref, imageSize*ones(1, 3));

genFunc2 = @(tn) cryo_project(vol_true, rand_rots(tn), imageSize, 'double');


genFuncs = {genFunc1, genFunc2};

printBegEndMsg('Image generation functions', false);

%% Generate ground truth images
printBegEndMsg('Generate ground truth images', true);

imUps = genFunc1(imageNo);
imCryo = genFunc2(imageNo);

save(fn, 'imUps', 'imCryo', '-append');

printBegEndMsg('Generate ground truth images', false);

%% Run test
printBegEndMsg('Running test', true);

cryoBackProjErr = zeros(length(lambda), length(bandlimit), imageNo);
cryoBackProjRelErr = zeros(length(lambda), length(bandlimit), imageNo);

for b=1:length(bandlimit)
    printBegEndMsg(num2str([b, length(bandlimit), bandlimit(b)], ...
        'Bandlimit %d of %d (banlimit = %d)'), true);

    bl = bandlimit(b);
    tDesign = loadtd(2*bl);
    
    for l=1:length(lambda)
        printBegEndMsg(num2str([l, length(lambda), lambda(l)], ...
            'Lambda %d of %d (lambda = %.3f)'), true);
    
        scalingParam = lambda(l);
        [shc, sh] = image2shc(randn(imageSize, imageSize), bl, ...
            tDesign, interval, scalingParam);
        [~, sh2] = shc2image(shc, bl, imageSize, interval, scalingParam);
        
        parfor t=1:imageNo
            % Cryo
            cryoSHC = image2shc(imCryo(:, :, t), bl, ...
                tDesign, interval, scalingParam, sh);
            cryoIm = shc2image(cryoSHC, bl, imageSize, interval, scalingParam, sh2);
            cryoBackProjErr(l, b, t) = norm(cryoIm - imCryo(:, :, t), 'fro');
            cryoBackProjRelErr(l, b, t) = cryoBackProjErr(l, b, t)/norm(imCryo(:, :, t), 'fro');
        end
        
        printBegEndMsg(num2str([l, length(lambda), lambda(l)], ...
            'Lambda %d of %d (lambda = %.3f)'), false);
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

yt = 1:5:length(lambda);
yticks(yt);
yticklabels(lambda(yt));

title('Cryo-EM images');
xlabel('Bandlimit');
ylabel('Scaling parameter');

set(gca,'ColorScale','log')

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;

%% Utility functions
function im = imagebyUps(tn, imSz)
im = zeros(imSz, imSz, tn);
for J=1:tn
    im(:, :, J) = imageByUpsampling(14, 16, imSz, round(0.2*imSz));
end
end
