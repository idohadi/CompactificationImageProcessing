function im = imageByUpsampling(lowBandlimit, highBandlimit, imageSize, paddingSize)
%%
% Call format
%   im = imageByUpsampling(lowBandlimit, highBandlimit)
%   im = imageByUpsampling(lowBandlimit, highBandlimit, imageSize)
%   im = imageByUpsampling(lowBandlimit, highBandlimit, imageSize, paddingSize)
% 
% Generate an image of size 
% (imageSize+2*paddingSize) x (imageSize+2*paddingSize) with black padding 
% paddingSize by the following procedure:
%   1.  Randomly genearting spherical harmonics coefficients up to bandlimit
%       lowBandlimit. 
%   2.  Back-project a region around the north pole to R^2, forming an 
%       image of size imageSize x imageSize. 
%   3.  Pad the image with black pixels in all four directions.
%   4.  Smooth this image using a Gaussian filter.
%   5.  Project this image onto the sphere and estimate its spehrical
%       harmonics coefficients up to bandlimit highBandlimit.
%   6.  Back-project the resulting image.
% 
% Input arguments
%   lowBandlimit    double      positive integer
%   highBandlimit   double      positive integer
%   imageSize       double      positive integer
%   paddingSize     double      positive integer, the number of black
%                               pixels to add to each dimension
% 
% Defaults
%   imageSize       150
%   paddingSize     50
%   
% Output arguments
%   im              double      the final image
% 
% Notes
%   This function relies heavily on my old code base.
% 
% Reference
%   None
% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input validation
narginchk(2, 4);
assert(round(lowBandlimit)==lowBandlimit ...
    && round(highBandlimit)==highBandlimit ...
    && lowBandlimit>0 && highBandlimit>0, ...
    'Bandlimit must be a positive integer.');

if nargin==2
    imageSize = 150;
    paddingSize = 50;
elseif nargin==3
    paddingSize = 50;
end
assert(round(imageSize)==imageSize ...
    && imageSize>0 ...
    && numel(imageSize)==1, ...
    'Image size must be a positive integer.');
assert(round(paddingSize)==paddingSize ...
    && paddingSize>0 ...
    && numel(paddingSize)==1, ...
    'Padding size must be a positive integer scalar.');

%% Image generation
% Step 1
shc = randomNormalizedSHC(lowBandlimit, 0);
shc = r2c(shc, lowBandlimit);

% Step 2
interval = [-0.5, 0.5];
a = 1.5;
im = shc2image(shc, lowBandlimit, imageSize, KondorProj(a), interval, interval);

% Step 3
im(1:paddingSize, :) = 0;
im(end-paddingSize:end, :) = 0;
im(:, 1:paddingSize) = 0;
im(:, end-paddingSize:end) = 0;

% Step 4
sigma = 5;
im = imgaussfilt(real(im), sigma);

% Step 5
td = loadtd('sf050.01302');
shc = image2shc(im, highBandlimit, td, KondorBackProj(a), interval, interval);

% Step 6
im = shc2image(shc, highBandlimit, imageSize, KondorProj(a), interval, interval);
im = real(im);
