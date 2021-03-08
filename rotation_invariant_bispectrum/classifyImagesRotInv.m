function [G, D] = classifyImagesRotInv(data, truncation, varargin)
%% 
% Call format
%   [G, D] = classifyImagesRotInv(data, truncation, angularLimits)
%   [G, D] = classifyImagesRotInv(data, truncation, angularLimits, __)
% 
% Identify close neighbors of images in data and denoise each image, using
% its neighbors.
% 
% 
% Input arguments
%   data                double      imageSize x imageSize x sampleSize
%                                   array, 
%   truncation          double      scalar, positive integer, bandlimit for
%                                   projection.
%   angularLimits       double 
% 
% Output arguments
%   G                   double      the nearest neighbors graph in a
%                                   matrix.
%   D                   double      The distances from MATLAB's nearest
%                                   neighbors function.
%                                   
% 
% Optional arguments
%   lowRank             The rank to use when applying the low rank
%                       approximation to the rotation-invariant 
%                       representation matrix.
%   Nneighbors          Number of nearest neighbors to find.
%   wpass               MATLAB's lowpass function parameter.
%                       If wpass=0, no low-pass filter is used.
% 
% Default optional arguments
%   lowRank             400
%   Nneighbors          50
%   wpass               0.05
% 
% Notes
%   None
% 
% Reference
%   [1] Halko, N., Martinsson, P.-G., Shkolnisky, Y., & Tygert, M. (2011). 
%       An Algorithm for the Principal Component Analysis of Large Data 
%       Sets. SIAM Journal on Scientific Computing, 33(5), 2580–2594. 
%       https://doi.org/10.1137/100804139

% ***********************************************************
% Author    Ido Hadi
% Email     idohadi@mail.tau.ac.il
% Year      2020
% ***********************************************************

%% Input handling
assert(size(data, 1)==size(data, 2), 'Images must be rectangular.');
imageSize = size(data, 1);
sampleSize = size(data, 3);
assert(isscalar(truncation) & truncation>=1 & round(truncation)==truncation, ...
    'Bandlimit must be a positive integer scalar.');

% Setting up optional input handling (name, value pairs)
p = inputParser;
addParameter(p, 'lowRank', 400, @(x) isscalar(x) & x>=1 & round(x)==x);
addParameter(p, 'Nneighbors', 50, @(x) isscalar(x) & x>=1);
addParameter(p, 'wpass', 0.05, @(x) isscalar(x) & x>=0 & x<1);

% Process the optional input
parse(p, varargin{:});
k = p.Results.lowRank;
Nneighbors = p.Results.Nneighbors;
wpass = p.Results.wpass;

%% Classify images
basis = ffb_basis(imageSize*ones(1, 2), truncation);
angularLimits = basis.k_max;

% Compute the length of the bispectrum vector
bLen = rotationInvariantBispectrumLength(truncation, angularLimits);

% Compute bispectra of data
clear bispectrum;
clear bispectrum_mex;

if wpass==0     % Don't use a low-pass filter
    rowFunc = struct('rowFunc', @rowFuncWODenoising, ...
            'args', {{data, basis, truncation, angularLimits, bLen}}, ...
            'm', sampleSize, ...
            'n', 2*bLen+angularLimits(1));
    [U, S, ~] = outOfCoreRandomizedSVD(rowFunc, k);
    b = U*S;
    clear U;
    clear S;
else            % Use a low-pass filter
    rowFunc = struct('rowFunc', @rowFuncWDenoising, ...
            'args', {{data, basis, truncation, angularLimits, bLen, wpass}}, ...
            'm', sampleSize, ...
            'n', 2*bLen+angularLimits(1));
    [U, S, ~] = outOfCoreRandomizedSVD(rowFunc, k);
    b = U*S;
    clear U;
    clear S;
end

% Compute nearest neighbors
[idx, D] = knnsearch(b, b, 'K', Nneighbors);

% Construct similarity matrix
G = sparse(repelem((1:sampleSize)', Nneighbors), ...
    reshape(idx', [numel(idx), 1]), 1, ...
    sampleSize, sampleSize, numel(idx));
G = G - diag(diag(G));

end

%% Utility functions
function row = rowFuncWDenoising(N, data, basis, truncation, angularLimits, bLen, wpass)
% With denoising
coeff = basis.expand(lowpass(data(:, :, N), wpass));
row = rotationInvariantBispectrum(coeff, truncation, angularLimits, bLen);
end

function row = rowFuncWODenoising(N, data, basis, truncation, angularLimits, bLen)
% Without denoising
coeff = basis.expand(data(:, :, N));
row = rotationInvariantBispectrum(coeff, truncation, angularLimits, bLen);
end
