%%
sampleSize = 10^4;

rotationAngles = 2*pi*rand(sampleSize, 1);
rotationAngles = sort(rotationAngles);
[translations, translationAngle, translationSize] ...
    = randTranslation(sampleSize);

path = 'groel.mrc';

imageSize = 101;
SNR = 1/50;

%%
[sampledImages, sampleVariances, basicViewOut] = generateImages(path, imageSize, rotationAngles, translations);
noisyImages = addNoise(sampledImages, sampleVariances, SNR);

fig1 = figure;
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

A = randperm(sampleSize, 9);
for J=1:9
    nexttile;
    imagesc(sampledImages(:, :, A(J))); 
    colormap('hot'); 
    colorbar;
end

fig2 = figure;
tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for J=1:sampleSize
    nexttile;
    imagesc(noisyImages(:, :, A(J))); 
    colormap('hot'); 
    colorbar;
end

save('RawData.mat', 'SNR', 'rotationAngles', 'translations', ...
    'sampledImages', 'noisyImages', 'basicViewOut');
