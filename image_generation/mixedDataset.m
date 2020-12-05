function [dataset, classRepresentatives, classMembership, ...
    denoisedDataset, translations, rotations] ...
        = mixedDataset(sampleSize, classesNo, filename, varargin)
%%
% Call format
%   mixedDataset(sampleSize, classesNo, filename)
%   mixedDataset(sampleSize, classesNo, filename, __)
%   [dataset, classRepresentatives, classMembership] = mixedDataset(__)
%   [dataset, classRepresentatives, classMembership, denoisedDataset] = mixedDataset(__)
%   [dataset, classRepresentatives, classMembership, denoisedDataset, translations] = mixedDataset(__)
%   [dataset, classRepresentatives, classMembership, denoisedDataset, translations, rotations] = mixedDataset(__)
% 
% Generate a dataset of noisy images. Each image in the dataset is one of
% classesNo images, rotated and translated and then noised.
% 
% The dataset is generated in two steps:
%   1. Class generation. Generate classesNo of images using the image 
%      generation function in imageGenFunc.
%   2. Dataset generation. Generate the dataset using the following steps:
%       a. Sample k in {1, 2, ..., classesNo} from the probability 
%          distribution defined by classProb. Repeat this sampleSize times.
%          Save the result in classMembership.
%       b. For every n in {1, 2, ..., sampleSize}, take class
%          representative classMembership(k), rotate it and then translate 
%          it by some random rotation and translation. Save the 
%          corresponding rotation and translation in rotations and 
%          translations, respectively. Save the result in denoisedDataset. 
%          Generate Gaussian noise with i.i.d. coordinates and standard 
%          deviation sigma and add it to the denoised image. Save the 
%          result in dataset.
% 
% 
% Input arguments
%   sampleSize              double      positive integer, sample size. 
%   classesNo               double      positive integer, number of classes.
%   filename                char        charecter vector, path to save the
%                                       dataset.
% 
% Output arguments
%   dataset                 double      imageSize x imageSize x sampleSize
%                                       array, dataset of noisy images.
%   classRepresentatives    double      imageSize x imageSize x classesNo
%                                       array, class representatives.
%   classMembership         double      sampleSize x 1 array,
%                                       classMembersihp(j) in {1, ..., classesNo}
%                                       means that the image dataset(:, :, j) 
%                                       was generated from classRepresentatives(:, :, j).
%   denoisedDataset         double      imageSize x imageSize x sampleSize
%                                       array, data of dataset images
%                                       without the noise.
%   translations            double      sampleSize x 2 array,
%                                       translations(j, :) is the
%                                       translation applied to the class
%                                       representative in order to generate
%                                       the j-th dataset image.
%   rotations               double      sampleSize x 1 array, 
%                                       rotations(j, :) is the
%                                       rotation applied to the class
%                                       representative in order to generate
%                                       the j-th dataset image.
% 
% Optional arguments
%   classProb       classNo x 1 double vector of positive numbers summing
%                   up to 1. A probability distribution for class images.
%                   classProb(j) is the probability of j.
%   imageGenFunc    a structure containing two fields:
%                   func    a function handle for a function returning an
%                           image, given constant arguments.
%                   args    the constant arguments.
%                   This is a function generating an image of size
%                   imageSize x imageSize.
%   sigma           the standard deviation of the added noise components.
%   lambda     parameter for the function randTranslation.
%   maxPixels  parameter for the function randTranslation.
% 
% Default optional arguments
%   classProb               classesNo x 1 double array, 
%                               classProb(k) = 1/classesNo
%   imageGenFunc            imageGenFunc.func = imageByUpsampling, 
%                           imageGenFunc.args = {14, 16, 101, round(0.2*101)}
%   sigma                   1
%   lambda             1.5
%   maxPixels          5
% 
% Notes
%   None
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input handling
% Input checks
assert(sampleSize>=1 & round(sampleSize)==sampleSize, ...
    'Sample size must be a positive integer.');
assert(classesNo>=1 & round(classesNo)==classesNo, ...
    'Number of classes must be a positive integer.');
assert(ischar(filename) & exist(filename, 'file')==0, ...
    'File name must be a char array and must point to a non-existent folder.');

% Setting up optional input handling (name, value pairs)
p = inputParser;
addParameter(p, 'classProb', ones(classesNo, 1)/classesNo, ...
    @(x) sum(x)==1 & all(x>=0));
addParameter(p, 'imageGenFunc', ...
    struct('func', @imageByUpsampling, 'args', {14, 16, 101, round(0.2*101)}), ...
    @(x) isfield(x, 'func') & isfield(x, 'args'));
addParameter(p, 'sigma', 1, @(x) isscalar(x) & x>=0);
addParameter(p, 'lambda', 1.5);
addParameter(p, 'maxPixels', 5);

% Process the optional input
parse(p, varargin);
classProb = p.Results.classProb;
imageGenFunc = p.Results.imageGenFunc;
sigma = p.Results.sigma;
lambda = p.Results.lambda;
maxPixels = p.Results.maxPixels;

%% Generate dataset
% Generate one class representative and extract image size
genFunc = @() imageGenFunc.func(imageGenFunc.args{:});
im = genFunc();
imageSize = size(im, 1);
assert(size(im, 1)==size(im, 2), ...
    'Image generation function does not create rectangular image.');

% Generate class representatives
classRepresentatives = zeros([imageSize, imageSize, classesNo]);
classRepresentatives(:, :, 1) = im;
clear im;
for J=2:classesNo
    classRepresentatives(:, :, J) = genFunc();
end

% Generate class membership vector
classMembership = randsample(1:classesNo, sampleSize, true, classProb);
classMembership = classMembership(:);

% Generate dataset
dataset = zeros(imageSize, imageSize, sampleSize);
denoisedDataset = zeros(imageSize, imageSize, sampleSize);
translations = zeros(sampleSize, 2);
rotations = 2*pi*rand(sampleSize, 1);

for J=1:sampleSize
    m = classMembership(J);
    t = randTranslation(maxPixels, lambda);
    
    denoisedIm = imrotate(classRepresentatives(:, :, m), ...
        360*rotations(J)/(2*pi), 'bicubic', 'crop');
    denoisedIm = imtranslate(denoisedIm, t, 'cubic', 'OutputView', 'same');
    im = denoisedIm + sigma*randn(imageSize, imageSize);
    
    % Save results
    translations(J, :) = t;
    denoisedDataset(:, :, J) = denoisedIm;
    dataset(:, :, J) = im;
end

% Save the result to file
save(filename, 'dataset', 'classRepresentatives', 'classMembership', ...
    'denoisedDataset', 'translations', 'rotations');
