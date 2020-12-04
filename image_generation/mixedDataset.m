function [dataset, classRepresentatives, classMembership, ...
    denoisedDataset, translations, rotations] ...
        = mixedDataset(sampleSize, classesNo, filename, imageGenFunc, ...
            classProp, transMaxPixels, transLambda, sigma)
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
%          distribution defined by classProp. Repeat this sampleSize times.
%          Save the result in classMembership.
%       b. For every n in {1, 2, ..., sampleSize}, take class
%          representative classMembership(k), rotate it and translate it by
%          some random rotation and translation. Save the corresponding
%          rotation and translation in rotations and translations,
%          respectively. Save the result in denoisedDataset. Generate
%          Gaussian noise with i.i.d. coordinates and standard deviation
%          sigma and add it to the denoised image. Save the result in dataset.
% 
% Input arguments
%   sampleSize              double      positive integer, sample size. 
%   classesNo               double      positive integer, number of classes.
%   filename                char        charecter vector, path to save the
%                                       dataset.
% 
% Output arguments
%   dataset                 double      
%   classRepresentatives    double      
%   classMembership         double      
%   denoisedDataset         double      
%   ranslations             double      
%   rotations               double      
% 
% Optional arguments
%   classProp       classNo x 1 double vector of positive numbers summing
%                   up to 1. A probability distribution for class images.
%                   classProp(j) is the probability of j.
%   imageGenFunc    a structure containing two fields:
%                   func    a function handle for a function returning an
%                           image, given constant arguments.
%                   args    the constant arguments.
%   sigma           the standard deviation of the added noise components.
%   transLambda     parameter for the function randTranslation.
%   transMaxPixels  parameter for the function randTranslation.
% 
% Default optional arguments
%   classProp               classesNo x 1 double array, 
%                               classProp(k) = 1/classesNo
%   imageGenFunc            imageGenFunc.func = imageByUpsampling, 
%                           imageGenFunc.args = {14, 16, 101, round(0.2*101)}
%   sigma                   1
%   transLambda             1.5
%   transMaxPixels          5
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
% TODO

% Setting up optional input handling (name, value pairs)
% TODO

%% Generate dataset
% TODO
