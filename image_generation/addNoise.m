function noisyImages = addNoise(sampledImages, sampleVariances, SNR)
%%
% Call format
%   generateImages(filename, imageSize, rotationAngles, translations)
%   generateImages(filename, imageSize, rotationAngles, translations, basicView)
%   sampledImages = generateImages(__)
%   [sampledImages, sampleVariances] = generateImages(__)
%   [sampledImages, sampleVariances, basicViewOut] = generateImages(__)
% 
% Generate a set of images, each is a Radon transform along the direction
% defined by basicView. The i-th image is rotated by rotationAngles(i) and
% translated by translation(:, i).
% 
% Input arguments
%   sampledImages       double      k x k x n array, n k x k images.
%   sampleVariances     double      n x 1 array, the sample variance of the
%                                   kth image.
%   SNR                 double      scalar, a number in (0,1), the desired
%                                   signal-to-noise ratio.
%   
% Output arguments
%   noisyImages         double      k x k x n array, sampled images with 
%                                   noise.
% 
% Notes
%   This function depends on the packages ASPIRE and SmallRotationToolbox.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Handle input
assert(size(sampledImages, 1)==size(sampledImages, 2), ...
    'Images must have equal width and height.');
assert(size(sampledImages, 3)==size(sampleVariances, 1), ...
    'Number of variances and number of images must be equal.');
assert(SNR>0 & SNR<1 & isscalar(SNR) & isnumeric(SNR), ...
    'SNR must be a double scalar in (0, 1).');

%% Add noise to images
% Calculate the standard deviation required to achieve the desired SNR
sampleVariances = squeeze(var(sampledImages, 0, [1, 2]));
sigma = sqrt(sampleVariances/SNR);
sigma = reshape(sigma, [1, 1, length(sigma)]);

% Compute the noise for each pixel of the image and add it to the images
noisyImages = sampledImages + sigma.*randn(size(sampledImages));
