%% Docs
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2022
% ***********************************************************

%% Pool setting
poolSize = 32;
try
    parpool(poolSize);
catch
    
end

%% RNG seed
rng(0, 'twister');

%% File setup
dt = datetime;
dt = datestr(dt, 'yyyy-mm-dd-HHMM');
fnNOEXT = ['OrbitApproximationBoundUps-', dt]; 
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
    'scalingParam', 'sampleSize', 'fnNOEXT');

printBegEndMsg('Setup parameters', false);

%% Generate images
printBegEndMsg('Dataset generation', true);

im = imageByUpsampling(14, 16, imageSize, round(0.2*imageSize));

save(fn, 'rotation', 'im', '-append');

printBegEndMsg('Dataset generation', false);

%% Rotation and translation
printBegEndMsg('Generate rotation and translation', true);

rotation = 2*pi*zeros(sampleSize, 1);
transRot = 2*pi*rand();
maxTrans = 10;
translation = maxTrans*rand(sampleSize, 1)*[cos(transRot), sin(transRot)];

save(fn, 'rotation', 'transRot', 'maxTrans', 'translation', '-append');

printBegEndMsg('Generate rotation and translation', false);


%% Run test
printBegEndMsg('Run test', true);

LHS = zeros(length(scalingParam), sampleSize);

for n=1:sampleSize
    printBegEndMsg(num2str([n, sampleSize], 'Iteration %d of %d'), true);
    
    tt1 = translation(n, :);
    rr1 = rotation(n);
    parfor J=1:length(scalingParam)
        shc = image2shc(im, bandlimit, tDesign, interval, scalingParam(J));
        P = contractionMap(cos(pi/4)*tt1/(0.5*imageSize), ...
            rr1, scalingParam(J));
        
        rotatedSHC = rotateSHC(shc, bandlimit, P, tDesign);

        R = [cos(rr1), sin(rr1); ...
            -sin(rr1), cos(rr1)];
        bb = -tt1;
        rr = - 360*rr1/(2*pi);
        imTmp = im;
        imTmp = imtranslate(imTmp, bb, 'cubic', 'OutputView', 'same');
        imTmp = imrotate(imTmp, rr, 'bicubic', 'crop');

        rotatedImSHC = image2shc(imTmp, bandlimit, tDesign, interval, ...
            scalingParam(J));

        LHS(J, n) = norm(rotatedImSHC - rotatedSHC, 2);
    end
    
    printBegEndMsg(num2str([n, sampleSize], 'Iteration %d of %d'), false);
end

save(fn, 'LHS', '-append');

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

% Draw polynomial
P = polyfit(log(scalingParam), log(mean(LHS, 2)), 1);
x = log(scalingParam);
y = polyval(P, x);
y = exp(y);

plot(scalingParam, y, '--', 'Color', 'red');

hold off;

xlabel('$\widetilde{\lambda}$', 'Interpreter', 'latex');
ylabel('2-norm approx.');
title('Eq. (2.11) test');

legend({'2-norm error (2.11 LHS)', ...
    '95% conf. intvl.', ...
    num2str(P(1), 'best fit line (slope = %.2f)')});

set(gca, 'xscale', 'log');
set(gca, 'yscale', 'log');

savefig(fig, [fnNOEXT, '.fig']);

%% Shut down the diary
diary off;
