function [avgedData, nearestNeighbors] = calssifyImages(data, bandlimit, varargin)
%%
% Call format
%   [avgedData, nearestNeighbors] = calssifyImages(data, bandlimit)
%   [avgedData, nearestNeighbors] = calssifyImages(data, bandlimit, Nneighbors)
% 
% Identify close neighbors of images in data and denoise each image, using
% its neighbors.
% 
% 
% Input arguments
%   data                double      imageSize x imageSize x sampleSize
%                                   array, 
%   bandlimit           double      scalar, positive integer, bandlimit for
%                                   projection.
% 
% Output arguments
%   avgedData           double      
%   nearestNeighbors    double      
% 
% Optional arguments
%   interval        Parameter for image2shc.
%   Nneighbors      Number of nearest neighbors to find.
%   scalingParam    Projection scaling parameter to use in image2shc.
% 
% Default optional arguments
%   interval        [-0.5, 0.5]
%   Nneighbors      20
%   scalingParam    1.5
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
assert(size(data, 1)==size(data, 2), 'Images must be rectangular.');
imageSize = size(data, 1);
sampleSize = size(data, 3);
assert(isscalar(bandlimit) & bandlimit>=1 & round(bandlimit)==bandlimit, ...
    'Bandlimit must be a positive integer scalar.');

% Setting up optional input handling (name, value pairs)
p = inputParser;
addParameter(p, 'interval', [-0.5, 0.5], @(x) numel(x)==2 & x(1)<x(2));
addParameter(p, 'Nneighbors', 30, @(x) isscalar(x) & x>=1);
addParameter(p, 'scalingParam', 1.5, @(x) isscalar(x) & x>0);

% Process the optional input
parse(p, varargin{:});
interval = p.Results.interval;
Nneighbors = p.Results.Nneighbors;
scalingParam = p.Results.scalingParam;

%% Classify images
% Compute bispectrum
td = loadtd(2*bandlimit + 2);
shc = image2shc(data(:, :, 1), bandlimit, td, 
for n=1:sampleSize
    
end

% TODO

% Compute nearest neighbors
% TODO

% Denoise data using the averaging
% TODO
