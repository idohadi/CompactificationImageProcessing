%% im2shc_demo
% TODO: eliminate the use in ASPIRE

imSz = 101;  % Size of image

% Generate a volume
root = aspire_root();
file_name = fullfile(root, 'projections', 'simulation', 'maps', 'cleanrib.mat');
f = load(file_name);
vol_true = cryo_downsample(f.volref, imSz*ones(1, 3));

% Project it to produce the image
rots_true = rand_rots(1);
im = cryo_project(vol_true, rots_true, imSz, 'double');

%% Prepare to use im2shc
scalingParam = 1;
N = 50; 

[assocLegMat, omega, theta, phi] = assocLegendreMatrices(N);
evalPoints = calcEvalPoints(theta, phi, scalingParam);

%% Estimate SHC using im2shc
tic;
shc = im2shc(im, evalPoints, N, assocLegMat, omega, phi);
toc;

fig = plotSphericalFunction(shc, N);

figure; imagesc(im); 
colormap('hot');
