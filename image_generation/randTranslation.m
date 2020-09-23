function [translation, translationAngle, translationSize] ...
    = randTranslation(maxPixels, lambda)
%%
% Call format
%   randTranslation()
%   randTranslation(maxPixels)
%   randTranslation(maxPixels, lambda)
%   translation = randTranslation(__)
%   [translation, translationAngle] = randTranslation(__)
%   [translation, translationAngle, translationSize] = randTranslation(__)
% 
% Compute a random translation. The translation is later used to transform
% translate an image with MATLAB's imtranslate function.
% 
% The translation is produced using the following procedure:
%   1. Uniformly sample a translate angle from [0, 2pi].
%   2. Samplea  translation size by repeatedly sampling the distribution
%       Exp[lambda] until a value below maxPixels is obtained.
%   3. Compute the translation by 
%       translation = translationSize * [cos(translateionAngle), sin(translateionAngle)]
% 
% Input arguments
%   maxPixels           double      positive scalar, the maximal number of 
%                                   pixels.
%   lambda              double      positive scalar, defining parameter of
%                                   the exponential distribution.
% Defaults
%   maxPixels           5
%   lambda              1.5
% 
% Output arguments
%   translation         double      1 x 2 array, described above.
%   translationAngle    double      scalar, a number in [0,2pi].
%   translationSize     double      scalar, described above.
% 
% Notes
%   This function performs no input checks.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input handling
narginchk(0, 2);

if nargin==0
    maxPixels = 5;
    lambda = 1.5;
elseif nargin==1
    lambda = 1.5;
end

%% Generate a random translation
translationAngle = 2*pi*rand();
translationSize = Inf;
while translationSize>maxPixels
    translationSize = exprnd(lambda);
end
translation ...
    = translationSize*[cos(translationAngle), sin(translationAngle)];
