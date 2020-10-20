function [translation, translationAngle, translationSize] ...
    = randTranslation(sampleSize, maxPixels, lambda)
%%
% Call format
%   randTranslation()
%   randTranslation(sampleSize)
%   randTranslation(sampleSize, maxPixels)
%   randTranslation(sampleSize, maxPixels, lambda)
%   translation = randTranslation(__)
%   [translation, translationAngle] = randTranslation(__)
%   [translation, translationAngle, translationSize] = randTranslation(__)
% 
% Compute random translations.
% 
% A translation is produced using the following procedure:
%   1. Uniformly sample a translate angle from [0, 2pi].
%   2. Samplea  translation size by repeatedly sampling the distribution
%       Exp[lambda] until a value below maxPixels is obtained.
%   3. Compute the translation by 
%       translation = translationSize * [cos(translateionAngle); sin(translateionAngle)]
% 
% Input arguments
%   sampleSize          double      positive integer, number of
%                                   translations to sample.
%   maxPixels           double      positive scalar, the maximal number of 
%                                   pixels.
%   lambda              double      positive scalar, defining parameter of
%                                   the exponential distribution.
% Defaults
%   sampleSize          1
%   maxPixels           5
%   lambda              1.5
% 
% Output arguments
%   translation         double      2 x sampleSize array, described above.
%   translationAngle    double      scalar, a number in [0,2pi].
%   translationSize     double      scalar, described above.
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
narginchk(0, 3);

if nargin==0
    sampleSize = 1;
    maxPixels = 5;
    lambda = 1.5;
elseif nargin==1
    maxPixels = 5;
    lambda = 1.5;
elseif nargin==2
    lambda = 1.5;
end

%% Generate a random translation
translationAngle = 2*pi*rand(1, sampleSize);

translationSize = exprnd(lambda, 1, sampleSize);
indicator = translationSize>maxPixels;
while any(indicator)
    translationSize(indicator) = exprnd(lambda);
    indicator = translationSize>maxPixels;
end

translation ...
    = translationSize.*[cos(translationAngle); sin(translationAngle)];
