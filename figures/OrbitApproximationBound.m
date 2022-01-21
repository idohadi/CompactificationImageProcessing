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
fnNOEXT = ['OrbitApproximationBound-', dt]; 
diary([fnNOEXT, '.log']); % Log file
fn = [fnNOEXT, '.mat']; % Output file

%% Setup parameters
printBegEndMsg('Setup parameters', true);

imageNo = 1;
imageSize = 101;

bandlimit = 50;
tDesign = loadtd(2*bandlimit);
interval = cos(pi/4)*[-1, 1];
scalingParam = 1:0.1:5;

sampleSize = 100;

save(fn, 'imageNo', 'imageSize', 'bandlimit', 'interval', ...
    'scalingParam', 'sampleSize');

printBegEndMsg('Setup parameters', false);

%% Generate images
printBegEndMsg('Dataset generation', true);

% Image by cryo-EM simulation
printBegEndMsg('Generate representatives', true);
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol_true = cryo_downsample(f.volref, imageSize*ones(1, 3));

rotation = rand_rots(imageNo);
im = cryo_project(vol_true, rotation, imageSize, 'double');

save(fn, 'rotation', 'im', '-append');

printBegEndMsg('Dataset generation', false);

%% Rotation and translation
printBegEndMsg('Generate rotation and translation', true);

rotation = 2*pi*rand(sampleSize, 1);
transRot = 2*pi*rand();
maxTrans = 10;
translation = maxTrans*rand(sampleSize, 1)*[cos(transRot), sin(transRot)];

save(fn, 'rotation', 'transRot', 'maxTrans', 'translation', '-append');

printBegEndMsg('Generate rotation and translation', false);


%% Run test
printBegEndMsg('Run test', true);

LHS = zeros(length(scalingParam), sampleSize);
RHS = zeros(length(scalingParam), sampleSize);

for n=1:sampleSize
    for J=1:length(scalingParam)
        shc = image2shc(im, bandlimit, tDesign, interval, scalingParam(J));
        P = contractionMap(translation(n), rotation(n), scalingParam(J));

        rotatedSHC = rotateSHC(shc, bandlimit, P.', tDesign);

        R = [cos(rotation(n)), sin(rotation(n)); ...
            -sin(rotation(n)), cos(rotation(n))];
        bb = - (R.' * translation(n).').';
        rr = -360*rotation(n)/(2*pi);
        imTmp = imrotate(im, rr, 'bicubic', 'crop');
        imTmp = imtranslate(imTmp, bb, 'cubic', 'OutputView', 'same');

        rotatedImSHC = image2shc(imTmp, bandlimit, tDesign, interval, ...
            scalingParam(J));

        LHS(J, n) = norm(rotatedImSHC - rotatedSHC, 2);
    end
end

save(fn, 'LHS', 'RHS', '-append');

printBegEndMsg('Run test', false);


%% Plot the results
fig = figure;

hold;

plot(scalingParam, mean(LHS, 2), 'MarkerFaceColor', 'blue');

polX = [scalingParam, flip(scalingParam)];
prctl95 = [prctile(LHS, 97.5, 2), ...
    flip(prctile(LHS, 2.5, 2))];
c = [17, 17, 17]/255;
fill(polX(:), prctl95(:), c, 'FaceAlpha', 0.3, 'LineStyle', 'none');

hold off;

xlabel('$\widetilde{\lambda}$', 'Interpreter', 'latex');
ylabel('2-norm approx.');
title('Eq. (2.11) test');

set(gca, 'yscale', 'log');

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;

%% Utility function
function out = logbound(b, lambda, lambdatilde)
out = 0.5*log(1 - cos(lambda*pi/lambdatilde));
out = out + log(4*pi);
out = out + norm(b) + lambda*pi;
out = out - 2*log(lambdatilde);
end